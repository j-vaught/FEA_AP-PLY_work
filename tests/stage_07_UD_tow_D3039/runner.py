"""Stage 07 - ASTM D3039 UD tow tension verification.

Author: J.C. Vaught

Post-matrix update: composite solids use the verified LAW12 + TYPE6/SOL_ORTH
row from references/openradioss_law_compatibility_matrix.md. Uniform off-axis
orientation is bound to references/openradioss_orientation_convention.md:
TYPE6 Ip=3, Iorth=0, Phi=theta.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


THIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parents[1]
RUNS_DIR = THIS_DIR / "runs"
RESULTS_DIR = THIS_DIR / "results"
FIGURES_DIR = THIS_DIR / "figures"
RUN_LOG = THIS_DIR / "run.log"
CARD_PATH = ROOT_DIR / "references" / "material_cards" / "im7_8552.json"

OR_ROOT = Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"
ENGINE = OR_ROOT / "exec" / "engine_linux64_gf"
ANIM_TO_VTK = OR_ROOT / "exec" / "anim_to_vtk_linux64_gf"
N_THREADS = int(os.environ.get("RAD_NT", "16"))


@dataclass(frozen=True)
class Material:
    rho: float
    e1: float
    e2: float
    e3: float
    nu12: float
    nu13: float
    nu23: float
    g12: float
    g13: float
    g23: float
    xt: float
    xc: float
    yt: float
    yc: float
    s12: float


@dataclass(frozen=True)
class Coupon:
    name: str
    theta_deg: float
    length: float
    width: float
    thickness: float
    nx: int
    ny: int
    nz: int


COUPONS = (
    Coupon("7A", 0.0, 0.250, 0.015, 0.001, 16, 4, 1),
    Coupon("7B", 45.0, 0.250, 0.025, 0.001, 16, 6, 1),
)


def fmt_f(*values: float) -> str:
    return "".join(f"{value:20.12g}" for value in values)


def fmt_i(*values: int) -> str:
    return "".join(f"{value:10d}" for value in values)


def load_material() -> Material:
    data = json.loads(CARD_PATH.read_text(encoding="utf-8"))
    elastic = data["elastic"]
    strength = data["strength"]
    return Material(
        rho=float(elastic["density_kg_m3"]),
        e1=float(elastic["E1_Pa"]),
        e2=float(elastic["E2_Pa"]),
        e3=float(elastic["E3_Pa"]),
        nu12=float(elastic["nu12"]),
        nu13=float(elastic["nu13"]),
        nu23=float(elastic["nu23"]),
        g12=float(elastic["G12_Pa"]),
        g13=float(elastic["G13_Pa"]),
        g23=float(elastic["G23_Pa"]),
        xt=float(strength["Xt_Pa"]),
        xc=float(strength["Xc_Pa"]),
        yt=float(strength["Yt_Pa"]),
        yc=float(strength["Yc_Pa"]),
        s12=float(strength["S12_Pa"]),
    )


def ex_offaxis(theta_deg: float, mat: Material) -> float:
    theta = math.radians(theta_deg)
    c2 = math.cos(theta) ** 2
    s2 = math.sin(theta) ** 2
    inv = (
        c2 * c2 / mat.e1
        + s2 * s2 / mat.e2
        + (1.0 / mat.g12 - 2.0 * mat.nu12 / mat.e1) * c2 * s2
    )
    return 1.0 / inv


def radioss_env() -> dict[str, str]:
    env = os.environ.copy()
    env["OR"] = str(OR_ROOT)
    env["RAD_CFG_PATH"] = str(OR_ROOT / "hm_cfg_files")
    env["RAD_H3D_PATH"] = str(OR_ROOT / "extlib" / "h3d" / "lib" / "linux64")
    reader = str(OR_ROOT / "extlib" / "hm_reader" / "linux64")
    env["LD_LIBRARY_PATH"] = reader + ":" + env.get("LD_LIBRARY_PATH", "")
    return env


def run_cmd(cmd: list[str], cwd: Path, log_lines: list[str]) -> subprocess.CompletedProcess[str]:
    log_lines.append("$ " + " ".join(cmd))
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        env=radioss_env(),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_lines.append(proc.stdout)
    return proc


def law12_block(mat: Material) -> list[str]:
    """LAW12 + TYPE6 is the verified solid-composite substitution."""
    nu31 = mat.nu13 * mat.e3 / mat.e1
    return [
        "/MAT/LAW12/1",
        "IM7_8552_canonical_Soden_WWFEII",
        "#              RHO_I",
        fmt_f(mat.rho),
        "#             MAT_EA              MAT_EB              MAT_EC",
        fmt_f(mat.e1, mat.e2, mat.e3),
        "#           MAT_PRAB            MAT_PRBC            MAT_PRCA",
        fmt_f(mat.nu12, mat.nu23, nu31),
        "#            MAT_GAB             MAT_GBC             MAT_GCA",
        fmt_f(mat.g12, mat.g23, mat.g13),
        "#           sigma_t1            sigma_t2            sigma_t3               delta",
        fmt_f(mat.xt, mat.yt, mat.yt, 0.05),
        "#           MAT_BETA                   n                fmax               Wpref",
        fmt_f(1.0, 1.0, 1.0, 1.0),
        "#          sigma_1yt           sigma_2yt           sigma_1yc           sigma_2yc",
        fmt_f(mat.xt, mat.yt, mat.xc, mat.yc),
        "#         sigma_12yt          sigma_12yc          sigma_23yt          sigma_23yc",
        fmt_f(mat.s12, mat.s12, mat.s12, mat.s12),
        "#          sigma_3yt           sigma_3yc          sigma_13yt          sigma_13yc",
        fmt_f(mat.yt, mat.yc, mat.s12, mat.s12),
        "#              alpha                  Ef                   c          EPS_RATE_0   STRFLAG",
        fmt_f(0.0, 0.0, 0.0, 0.0) + fmt_i(1),
    ]


def mesh_blocks(coupon: Coupon) -> tuple[list[str], list[int], list[int]]:
    nodes: list[str] = ["/NODE"]
    node_id: dict[tuple[int, int, int], int] = {}
    nid = 1
    for k in range(coupon.nz + 1):
        for j in range(coupon.ny + 1):
            for i in range(coupon.nx + 1):
                node_id[(i, j, k)] = nid
                x = coupon.length * i / coupon.nx
                y = coupon.width * (j / coupon.ny - 0.5)
                z = coupon.thickness * (k / coupon.nz - 0.5)
                nodes.append(f"{nid:10d}{x:20.12g}{y:20.12g}{z:20.12g}")
                nid += 1

    bricks = ["/BRICK/1"]
    eid = 1
    for k in range(coupon.nz):
        for j in range(coupon.ny):
            for i in range(coupon.nx):
                conn = [
                    node_id[(i, j, k)],
                    node_id[(i + 1, j, k)],
                    node_id[(i + 1, j + 1, k)],
                    node_id[(i, j + 1, k)],
                    node_id[(i, j, k + 1)],
                    node_id[(i + 1, j, k + 1)],
                    node_id[(i + 1, j + 1, k + 1)],
                    node_id[(i, j + 1, k + 1)],
                ]
                bricks.append(fmt_i(eid, *conn))
                eid += 1

    left = [node_id[(0, j, k)] for k in range(coupon.nz + 1) for j in range(coupon.ny + 1)]
    right = [node_id[(coupon.nx, j, k)] for k in range(coupon.nz + 1) for j in range(coupon.ny + 1)]
    return nodes + ["/PART/1", "ud_coupon", fmt_i(1, 1, 0)] + bricks, left, right


def group_block(group_id: int, name: str, node_ids: list[int]) -> list[str]:
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for i in range(0, len(node_ids), 10):
        lines.append(fmt_i(*node_ids[i : i + 10]))
    return lines


def type6_property(theta_deg: float) -> list[str]:
    return [
        "/PROP/TYPE6/1",
        "law12_type6_phi_ip3_orientation",
        "#   Isolid    Ismstr               Icpre  Itetra10     Inpts   Itetra4    Iframe                  Dn",
        fmt_i(24, 4) + f"{1:20d}{0:10d}{0:10d}{0:10d}{2:10d}{0.0:20.12g}",
        "#                 qa                  qb                   h",
        fmt_f(0.0, 0.0, 0.0),
        "#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth",
        fmt_f(1.0, 0.0, 0.0) + fmt_i(0, 3, 0),
        "#                Phi                 Px                  Py                  Pz",
        fmt_f(theta_deg, 0.0, 0.0, 0.0),
        "#             dt_min   istrain      IHKT",
        fmt_f(0.0) + fmt_i(0, 0),
    ]


def write_starter(coupon: Coupon, mat: Material) -> Path:
    mesh, left, right = mesh_blocks(coupon)
    anchor_yz = [left[0]]
    anchor_z = [left[-1]]
    job = f"stage07_{coupon.name}_law12_type6"
    lines = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        "      2023         0",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 07 {coupon.name} theta {coupon.theta_deg:.1f} deg LAW12 TYPE6",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 4) + f"{0:20d}{2:40d}",
    ]
    lines.extend(law12_block(mat))
    lines.extend(mesh)
    lines.extend(type6_property(coupon.theta_deg))
    lines.extend(
        [
            "/BCS/1",
            "left_grip_x_fixed",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   100 000{0:10d}{100:10d}",
            "/BCS/2",
            "left_anchor_yz",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   011 000{0:10d}{102:10d}",
            "/BCS/3",
            "left_anchor_z",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   001 000{0:10d}{103:10d}",
            "/FUNCT/1",
            "unit_ramp",
            fmt_f(0.0, 0.0),
            fmt_f(1.0e-3, 1.0),
            "/IMPDISP/1",
            "right_grip_x",
            "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
            f"{1:10d}{'X':>10}{0:10d}{0:10d}{101:10d}{0:10d}{0:10d}",
            "#            Scale_x             Scale_y              Tstart               Tstop",
            fmt_f(1.0, 0.002 * coupon.length, 0.0, 0.0),
        ]
    )
    lines.extend(group_block(100, "left_grip", left))
    lines.extend(group_block(101, "right_grip", right))
    lines.extend(group_block(102, "left_anchor_yz", anchor_yz))
    lines.extend(group_block(103, "left_anchor_z", anchor_z))
    lines.extend(["/END", ""])
    path = RUNS_DIR / f"{job}_0000.rad"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_engine(job: str) -> Path:
    path = RUNS_DIR / f"{job}_0001.rad"
    run_time = 1.0e-3
    lines = [
        "#RADIOSS ENGINE",
        "/ANIM/DT",
        fmt_f(run_time, run_time),
        "/ANIM/VECT/DISP",
        "/ANIM/VECT/FREAC",
        "/ANIM/BRICK/TENS/STRESS/ALL",
        "/ANIM/BRICK/TENS/STRAIN/ALL",
        "/ANIM/GZIP",
        "/TFILE/4",
        fmt_f(run_time / 20.0),
        "/RFILE",
        fmt_i(1000),
        "/PRINT/-100/55",
        f"/RUN/{job}/1",
        fmt_f(run_time),
        "/VERS/2023",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def convert_anim_to_vtk(job: str, log_lines: list[str]) -> Path:
    anim = RUNS_DIR / f"{job}A001"
    gz = RUNS_DIR / f"{job}A001.gz"
    if gz.exists():
        with gzip.open(gz, "rb") as src, anim.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    if not anim.exists():
        raise FileNotFoundError(f"animation frame not found for {job}")
    cmd = [str(ANIM_TO_VTK), str(anim)]
    log_lines.append("$ " + " ".join(cmd))
    proc = subprocess.run(
        cmd,
        cwd=str(RUNS_DIR),
        env=radioss_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.stderr:
        log_lines.append(proc.stderr.decode("utf-8", errors="replace"))
    log_lines.append(f"[exit {proc.returncode}]")
    vtk = RUNS_DIR / f"{job}A001.vtk"
    if proc.returncode != 0:
        raise RuntimeError(f"anim_to_vtk failed for {job}")
    if proc.stdout.lstrip().startswith(b"# vtk"):
        vtk.write_bytes(proc.stdout[proc.stdout.find(b"# vtk") :])
    if not vtk.exists():
        raise FileNotFoundError(f"anim_to_vtk produced no VTK for {job}")
    return vtk


def stress_array(grid):
    for name in grid.cell_data.keys():
        if "Strs" in name:
            return grid.cell_data[name]
    raise KeyError(f"no stress array in VTK cell data: {list(grid.cell_data.keys())}")


def strain_array(grid):
    for name in grid.cell_data.keys():
        if "Stra" in name:
            return grid.cell_data[name]
    raise KeyError(f"no strain array in VTK cell data: {list(grid.cell_data.keys())}")


def _point_array(grid, contains: str):
    for name in grid.point_data.keys():
        if contains in name:
            return grid.point_data[name]
    raise KeyError(f"no point data containing {contains}; available={list(grid.point_data.keys())}")


def extract_modulus(coupon: Coupon, vtk_path: Path, e_ref: float) -> dict[str, object]:
    import numpy as np
    import pyvista as pv  # type: ignore[import-not-found]

    grid = pv.read(str(vtk_path))
    stress = np.asarray(stress_array(grid), dtype=float)
    centers = grid.cell_centers().points
    x0 = centers[:, 0].min()
    x1 = centers[:, 0].max()
    gauge = (centers[:, 0] >= x0 + 0.35 * (x1 - x0)) & (centers[:, 0] <= x0 + 0.65 * (x1 - x0))
    sigma_x = float(np.mean(stress[gauge, 0]))
    disp = np.asarray(_point_array(grid, "Displacement"), dtype=float)
    points = np.asarray(grid.points, dtype=float)
    original_x = points[:, 0] - disp[:, 0]
    left = original_x < original_x.min() + 1.0e-8
    right = original_x > original_x.max() - 1.0e-8
    eps_x = abs(float(np.mean(disp[right, 0]) - np.mean(disp[left, 0]))) / coupon.length
    reaction = np.asarray(_point_array(grid, "Reaction"), dtype=float)
    force_x = float(np.sum(reaction[right, 0]))
    reaction_sigma_x = force_x / (coupon.width * coupon.thickness)
    modulus = abs(sigma_x) / eps_x
    rel_err = abs(modulus - e_ref) / e_ref
    return {
        "solver_modulus_pa": modulus,
        "relative_error_pct": 100.0 * rel_err,
        "pass": rel_err <= 0.02,
        "sigma_x_mean_pa": sigma_x,
        "stress_measure": "mean_gauge_sigma_x",
        "reaction_force_x_n": force_x,
        "reaction_sigma_x_pa": reaction_sigma_x,
        "strain_x": eps_x,
        "strain_measure": "engineering_displacement_difference",
        "vtk_path": str(vtk_path),
    }


def run_stage() -> tuple[dict[str, object], list[dict[str, object]], list[str]]:
    mat = load_material()
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    log_lines = [f"material_card={CARD_PATH}", f"starter={STARTER}", f"engine={ENGINE}"]
    rows: list[dict[str, object]] = []
    all_started = True
    all_engine = True
    all_pass = True

    for coupon in COUPONS:
        e_ref = mat.e1 if coupon.theta_deg == 0.0 else ex_offaxis(coupon.theta_deg, mat)
        starter = write_starter(coupon, mat)
        job = starter.name.removesuffix("_0000.rad")
        engine = write_engine(job)
        starter_proc = run_cmd([str(STARTER), "-i", starter.name, "-nt", str(N_THREADS)], RUNS_DIR, log_lines)
        all_started = all_started and starter_proc.returncode == 0
        engine_rc: int | None = None
        extracted: dict[str, object] = {
            "solver_modulus_pa": "",
            "relative_error_pct": "",
            "pass": False,
            "sigma_x_mean_pa": "",
            "stress_measure": "",
            "reaction_force_x_n": "",
            "reaction_sigma_x_pa": "",
            "strain_x": "",
            "strain_measure": "",
            "vtk_path": "",
        }
        if starter_proc.returncode == 0:
            engine_proc = run_cmd([str(ENGINE), "-i", engine.name, "-nt", str(N_THREADS)], RUNS_DIR, log_lines)
            engine_rc = engine_proc.returncode
            all_engine = all_engine and engine_proc.returncode == 0
            if engine_proc.returncode == 0:
                vtk = convert_anim_to_vtk(job, log_lines)
                extracted = extract_modulus(coupon, vtk, e_ref)
        all_pass = all_pass and bool(extracted["pass"])

        rows.append(
            {
                "run": coupon.name,
                "theta_deg": coupon.theta_deg,
                "reference_modulus_pa": e_ref,
                "material_property": "LAW12_TYPE6_SOL_ORTH",
                "starter_rc": starter_proc.returncode,
                "engine_rc": "" if engine_rc is None else engine_rc,
                "solver_modulus_pa": extracted["solver_modulus_pa"],
                "relative_error_pct": extracted["relative_error_pct"],
                "sigma_x_mean_pa": extracted["sigma_x_mean_pa"],
                "stress_measure": extracted["stress_measure"],
                "reaction_force_x_n": extracted["reaction_force_x_n"],
                "reaction_sigma_x_pa": extracted["reaction_sigma_x_pa"],
                "strain_x": extracted["strain_x"],
                "strain_measure": extracted["strain_measure"],
                "vtk_path": extracted["vtk_path"],
                "verdict": "PASS" if extracted["pass"] else "FAIL",
            }
        )

    metrics = {
        "canonical_material_property": "LAW12 + TYPE6/SOL_ORTH",
        "matrix_evidence": "LAW12 row: TYPE6/SOL_ORTH solid OK; TYPE14 solid = B3047",
        "orientation_evidence": "references/openradioss_orientation_convention.md uniform property angle row: Ip=3, Iorth=0, Phi=theta",
        "strain_measurement": "engineering strain from displacement difference; VTK cell Stra[0] is not used for modulus",
        "starter_all_ok": all_started,
        "engine_all_ok": all_engine,
        "d3039_modulus_gate_evaluated": True,
        "d3039_modulus_gate_pass": all_pass,
    }
    return metrics, rows, log_lines


def write_outputs(metrics: dict[str, object], rows: list[dict[str, object]], log_lines: list[str], wall_s: float) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    timeseries = RESULTS_DIR / "timeseries.csv"
    with timeseries.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    git_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=str(ROOT_DIR),
        text=True,
        stdout=subprocess.PIPE,
        check=False,
    ).stdout.strip()
    verdict = "PASS" if bool(metrics["d3039_modulus_gate_pass"]) else "FAIL"
    results = {
        "stage": 7,
        "verdict": verdict,
        "metrics": {**metrics, "wall_clock_s": wall_s},
        "reference": {
            "material_card": str(CARD_PATH),
            "on_axis_E1_Pa": load_material().e1,
            "off_axis_45deg_Ex_Pa": ex_offaxis(45.0, load_material()),
        },
        "tolerance": {"modulus_relative_error": 0.02},
        "git_sha": git_sha,
    }
    (RESULTS_DIR / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    (FIGURES_DIR / "stage07_d3039_probe.typ").write_text(
        "\n".join(
            [
                '#set page(width: 170mm, height: auto, margin: 10mm)',
                '#let rows = csv("../results/timeseries.csv")',
                '#text(size: 12pt, weight: "bold")[Stage 07 D3039 starter probe]',
                '#v(5pt)',
                '#table(',
                '  columns: (22mm, 22mm, 42mm, 32mm, 30mm),',
                '  stroke: rgb("#5C5C5C"),',
                '  [Run], [Theta], [Reference Pa], [Error %], [Verdict],',
                '  ..rows.map(r => ([#r.at(0)], [#r.at(1)], [#r.at(2)], [#r.at(7)], [#r.at(15)])).flatten(),',
                ')',
                "",
            ]
        ),
        encoding="utf-8",
    )

    blocker = THIS_DIR / "blocker.md"
    if verdict == "PASS" and blocker.exists():
        blocker.unlink()
    elif verdict != "PASS":
        blocker.write_text(
            "\n".join(
                [
                    "# Stage 07 Blocker - TYPE6 off-axis D3039 mismatch",
                    "",
                    "Author: J.C. Vaught",
                    "",
                    "The post-matrix `LAW12 + TYPE6/SOL_ORTH` D3039 decks use the verified property-level orientation recipe from `references/openradioss_orientation_convention.md`: `Ip=3`, `Iorth=0`, `Phi=theta`. They start, run, write animation output, convert through `anim_to_vtk`, and post-process with PyVista. The remaining blocker is numerical, not starter parsing.",
                    "",
                    "Observed results from `results/timeseries.csv`:",
                    "",
                    *[
                        f"- Run {row['run']}, {row['theta_deg']} deg: `E_FEM = {float(row['solver_modulus_pa']) / 1.0e9:.3f} GPa` vs reference `{float(row['reference_modulus_pa']) / 1.0e9:.3f} GPa`, error `{float(row['relative_error_pct']):.3f}%`, {row['verdict']}."
                        for row in rows
                    ],
                    "",
                    "The modulus calculation uses mean gauge `sigma_x` divided by engineering strain from the displacement field. It does not use VTK cell `Stra[0]` for the gate.",
                    "",
                    f"Verdict: `{verdict}` after the post-rotation convention update.",
                    "",
                ]
            ),
            encoding="utf-8",
        )
    RUN_LOG.write_text("\n".join(log_lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if not args.all:
        parser.error("use --all")
    start = time.perf_counter()
    metrics, rows, log_lines = run_stage()
    wall_s = time.perf_counter() - start
    write_outputs(metrics, rows, log_lines, wall_s)
    verdict = "PASS" if bool(metrics["d3039_modulus_gate_pass"]) else "FAIL"
    print(json.dumps({"stage": 7, "verdict": verdict, "wall_clock_s": wall_s}, indent=2))
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
