"""Stage 10 - UD tow-wise direct mesoscale verification.

Author: J.C. Vaught

This stage is independent of the composite orientation convention because the
constituents are isotropic fiber and matrix solids. The runner builds a
scripted HEXA8 mesoscale cell with a circular fiber bundle, runs three
OpenRadioss cases, converts the final animation frame to VTK, and compares
effective E1, E2, and G12 against Halpin-Tsai / rule-of-mixtures references.
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

import numpy as np


THIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parents[1]
RUNS_DIR = THIS_DIR / "runs"
RESULTS_DIR = THIS_DIR / "results"
FIGURES_DIR = THIS_DIR / "figures"
RUN_LOG = THIS_DIR / "run.log"

OR_ROOT = Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"
ENGINE = OR_ROOT / "exec" / "engine_linux64_gf"
ANIM_TO_VTK = OR_ROOT / "exec" / "anim_to_vtk_linux64_gf"
N_THREADS = int(os.environ.get("RAD_NT", "16"))

LX = 100.0e-6
LY = 100.0e-6
LZ = 50.0e-6
NX = 20
NY = 20
NZ = 4
FIBER_RADIUS = 0.432 * LX
STRAIN = 1.0e-3
RUN_TIME = 2.0e-7

E_F = 230.0e9
NU_F = 0.20
RHO_F = 1780.0
E_M = 4.08e9
NU_M = 0.39
RHO_M = 1300.0

ZETA_E2 = 2.0
ZETA_G12 = 1.0
HT_REL_TOL = 0.05


@dataclass(frozen=True)
class MeshData:
    nodes: list[str]
    fiber_bricks: list[tuple[int, list[int]]]
    matrix_bricks: list[tuple[int, list[int]]]
    groups: dict[str, list[int]]
    vf_actual: float


@dataclass(frozen=True)
class LoadCase:
    name: str
    target_key: str
    component: str
    denominator: float


LOAD_CASES = (
    LoadCase("axial_z", "E1_ROM_Pa", "sigma_zz", STRAIN),
    LoadCase("transverse_x", "E2_HT_Pa", "sigma_xx", STRAIN),
    LoadCase("shear_xz", "G12_HT_Pa", "sigma_xz", STRAIN),
)


def fmt_f(*values: float) -> str:
    return "".join(f"{value:20.12g}" for value in values)


def fmt_i(*values: int) -> str:
    return "".join(f"{value:10d}" for value in values)


def shear_modulus(e: float, nu: float) -> float:
    return e / (2.0 * (1.0 + nu))


def halpin_tsai(phase: float, matrix: float, vf: float, zeta: float) -> float:
    ratio = phase / matrix
    eta = (ratio - 1.0) / (ratio + zeta)
    return matrix * (1.0 + zeta * eta * vf) / (1.0 - eta * vf)


def reference_targets(vf: float) -> dict[str, float]:
    gf = shear_modulus(E_F, NU_F)
    gm = shear_modulus(E_M, NU_M)
    return {
        "E1_ROM_Pa": vf * E_F + (1.0 - vf) * E_M,
        "E2_HT_Pa": halpin_tsai(E_F, E_M, vf, ZETA_E2),
        "G12_HT_Pa": halpin_tsai(gf, gm, vf, ZETA_G12),
    }


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


def build_mesh() -> MeshData:
    node_id: dict[tuple[int, int, int], int] = {}
    nodes = ["/NODE"]
    nid = 1
    for k in range(NZ + 1):
        for j in range(NY + 1):
            for i in range(NX + 1):
                node_id[(i, j, k)] = nid
                x = LX * i / NX
                y = LY * j / NY
                z = LZ * k / NZ
                nodes.append(f"{nid:10d}{x:20.12g}{y:20.12g}{z:20.12g}")
                nid += 1

    fiber: list[tuple[int, list[int]]] = []
    matrix: list[tuple[int, list[int]]] = []
    cx = 0.5 * LX
    cy = 0.5 * LY
    eid = 1
    for k in range(NZ):
        for j in range(NY):
            for i in range(NX):
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
                xc = LX * (i + 0.5) / NX
                yc = LY * (j + 0.5) / NY
                if (xc - cx) ** 2 + (yc - cy) ** 2 <= FIBER_RADIUS**2:
                    fiber.append((eid, conn))
                else:
                    matrix.append((eid, conn))
                eid += 1

    groups = {
        "x0": sorted(node_id[(0, j, k)] for k in range(NZ + 1) for j in range(NY + 1)),
        "x1": sorted(node_id[(NX, j, k)] for k in range(NZ + 1) for j in range(NY + 1)),
        "z0": sorted(node_id[(i, j, 0)] for j in range(NY + 1) for i in range(NX + 1)),
        "z1": sorted(node_id[(i, j, NZ)] for j in range(NY + 1) for i in range(NX + 1)),
        "all_y": sorted(node_id[(i, j, k)] for k in range(NZ + 1) for j in range(NY + 1) for i in range(NX + 1)),
        "anchor_xy": [node_id[(0, 0, 0)]],
        "anchor_yz": [node_id[(NX, 0, 0)]],
    }
    for k in range(NZ + 1):
        groups[f"zlayer_{k}"] = sorted(node_id[(i, j, k)] for j in range(NY + 1) for i in range(NX + 1))

    vf_actual = len(fiber) / (NX * NY * NZ)
    return MeshData(nodes=nodes, fiber_bricks=fiber, matrix_bricks=matrix, groups=groups, vf_actual=vf_actual)


def group_block(group_id: int, name: str, node_ids: list[int]) -> list[str]:
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for i in range(0, len(node_ids), 10):
        lines.append(fmt_i(*node_ids[i : i + 10]))
    return lines


def material_and_property_blocks() -> list[str]:
    return [
        "/MAT/ELAST/1",
        "IM7_fiber_isotropic",
        fmt_f(RHO_F, 0.0),
        fmt_f(E_F, NU_F),
        "/MAT/ELAST/2",
        "8552_matrix_isotropic",
        fmt_f(RHO_M, 0.0),
        fmt_f(E_M, NU_M),
        "/PROP/SOLID/1",
        "fiber_solid",
        "#   Isolid    Ismstr               Icpre               Inpts    Itetra    Iframe                  dn",
        fmt_i(24, 0) + f"{1:20d}{0:20d}{0:10d}{0:10d}{0:20d}",
        "#                q_a                 q_b                   h            LAMBDA_V                MU_V",
        fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
        "#             dt_min   istrain      IHKT",
        fmt_f(0.0) + fmt_i(0, 0),
        "/PROP/SOLID/2",
        "matrix_solid",
        "#   Isolid    Ismstr               Icpre               Inpts    Itetra    Iframe                  dn",
        fmt_i(24, 0) + f"{1:20d}{0:20d}{0:10d}{0:10d}{0:20d}",
        "#                q_a                 q_b                   h            LAMBDA_V                MU_V",
        fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
        "#             dt_min   istrain      IHKT",
        fmt_f(0.0) + fmt_i(0, 0),
    ]


def bcs_and_loads(case: LoadCase, mesh: MeshData) -> list[str]:
    base_groups = {
        "x0": 100,
        "x1": 101,
        "z0": 102,
        "z1": 103,
        "all_y": 104,
        "anchor_xy": 105,
        "anchor_yz": 106,
    }
    lines: list[str] = [
        "/FUNCT/1",
        "unit_ramp",
        fmt_f(0.0, 0.0),
        fmt_f(RUN_TIME, 1.0),
    ]
    if case.name == "axial_z":
        lines.extend(
            [
                "/BCS/1",
                "z0_fixed_z",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   001 000{0:10d}{base_groups['z0']:10d}",
                "/BCS/2",
                "anchor_xy",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   110 000{0:10d}{base_groups['anchor_xy']:10d}",
                "/BCS/3",
                "anchor_y",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   010 000{0:10d}{base_groups['anchor_yz']:10d}",
                "/IMPDISP/1",
                "z1_axial_z",
                "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
                f"{1:10d}{'Z':>10}{0:10d}{0:10d}{base_groups['z1']:10d}{0:10d}{0:10d}",
                "#            Scale_x             Scale_y              Tstart               Tstop",
                fmt_f(1.0, STRAIN * LZ, 0.0, 0.0),
            ]
        )
    elif case.name == "transverse_x":
        lines.extend(
            [
                "/BCS/1",
                "x0_fixed_x",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   100 000{0:10d}{base_groups['x0']:10d}",
                "/BCS/2",
                "anchor_yz",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   011 000{0:10d}{base_groups['anchor_xy']:10d}",
                "/BCS/3",
                "anchor_z",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   001 000{0:10d}{base_groups['anchor_yz']:10d}",
                "/IMPDISP/1",
                "x1_transverse_x",
                "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
                f"{1:10d}{'X':>10}{0:10d}{0:10d}{base_groups['x1']:10d}{0:10d}{0:10d}",
                "#            Scale_x             Scale_y              Tstart               Tstop",
                fmt_f(1.0, STRAIN * LX, 0.0, 0.0),
            ]
        )
    elif case.name == "shear_xz":
        lines.extend(
            [
                "/BCS/1",
                "z0_fixed_x",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   100 000{0:10d}{base_groups['z0']:10d}",
                "/BCS/2",
                "anchor_yz",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   011 000{0:10d}{base_groups['anchor_xy']:10d}",
                "/BCS/3",
                "anchor_y",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   010 000{0:10d}{base_groups['anchor_yz']:10d}",
                "/IMPDISP/1",
                "z1_shear_x",
                "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
                f"{1:10d}{'X':>10}{0:10d}{0:10d}{base_groups['z1']:10d}{0:10d}{0:10d}",
                "#            Scale_x             Scale_y              Tstart               Tstop",
                fmt_f(1.0, STRAIN * LZ, 0.0, 0.0),
            ]
        )
    else:
        raise ValueError(case.name)

    for name, gid in base_groups.items():
        lines.extend(group_block(gid, name, mesh.groups[name]))
    return lines


def write_starter(case: LoadCase, mesh: MeshData) -> Path:
    job = f"stage10_{case.name}"
    lines = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        "      2023         0",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 10 UD mesoscale {case.name}",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 0) + f"{0:20d}{2:40d}",
    ]
    lines.extend(material_and_property_blocks())
    lines.extend(mesh.nodes)
    lines.extend(["/PART/1", "fiber_part", fmt_i(1, 1, 0), "/BRICK/1"])
    for eid, conn in mesh.fiber_bricks:
        lines.append(fmt_i(eid, *conn))
    lines.extend(["/PART/2", "matrix_part", fmt_i(2, 2, 0), "/BRICK/2"])
    for eid, conn in mesh.matrix_bricks:
        lines.append(fmt_i(eid, *conn))
    lines.extend(bcs_and_loads(case, mesh))
    lines.extend(["/END", ""])
    path = RUNS_DIR / f"{job}_0000.rad"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_engine(job: str) -> Path:
    path = RUNS_DIR / f"{job}_0001.rad"
    lines = [
        "#RADIOSS ENGINE",
        "/ANIM/DT",
        fmt_f(RUN_TIME, RUN_TIME),
        "/ANIM/VECT/DISP",
        "/ANIM/BRICK/TENS/STRESS/ALL",
        "/ANIM/BRICK/TENS/STRAIN/ALL",
        "/ANIM/GZIP",
        "/TFILE/4",
        fmt_f(RUN_TIME / 20.0),
        "/RFILE",
        fmt_i(1000),
        "/PRINT/-100/55",
        f"/RUN/{job}/1",
        fmt_f(RUN_TIME),
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


def _cell_array(grid, contains: str):
    for name in grid.cell_data.keys():
        if contains in name:
            return grid.cell_data[name]
    raise KeyError(f"no cell data containing {contains}; available={list(grid.cell_data.keys())}")


def extract_modulus(case: LoadCase, vtk_path: Path) -> float:
    import pyvista as pv  # type: ignore[import-not-found]

    grid = pv.read(str(vtk_path))
    stress = np.asarray(_cell_array(grid, "Strs"), dtype=float)
    if stress.shape[1] >= 9:
        sigma = {
            "sigma_xx": stress[:, 0],
            "sigma_zz": stress[:, 8],
            "sigma_xz": 0.5 * (stress[:, 2] + stress[:, 6]),
        }[case.component]
    else:
        sigma = {
            "sigma_xx": stress[:, 0],
            "sigma_zz": stress[:, 2],
            "sigma_xz": stress[:, 5],
        }[case.component]
    return abs(float(np.mean(sigma))) / case.denominator


def run_stage() -> tuple[dict[str, object], list[dict[str, object]], list[str]]:
    mesh = build_mesh()
    refs = reference_targets(mesh.vf_actual)
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    log_lines = [
        f"starter={STARTER}",
        f"engine={ENGINE}",
        f"mesh={NX}x{NY}x{NZ}",
        f"fiber_bricks={len(mesh.fiber_bricks)}",
        f"matrix_bricks={len(mesh.matrix_bricks)}",
        f"vf_actual={mesh.vf_actual:.6f}",
    ]
    rows: list[dict[str, object]] = []
    all_started = True
    all_engine = True
    all_pass = True

    for case in LOAD_CASES:
        starter = write_starter(case, mesh)
        job = starter.name.removesuffix("_0000.rad")
        engine = write_engine(job)
        starter_proc = run_cmd([str(STARTER), "-i", starter.name, "-nt", str(N_THREADS)], RUNS_DIR, log_lines)
        all_started = all_started and starter_proc.returncode == 0
        engine_rc: int | str = ""
        modulus: float | str = ""
        rel_err: float | str = ""
        vtk_path = ""
        passed = False
        if starter_proc.returncode == 0:
            engine_proc = run_cmd([str(ENGINE), "-i", engine.name, "-nt", str(N_THREADS)], RUNS_DIR, log_lines)
            engine_rc = engine_proc.returncode
            all_engine = all_engine and engine_proc.returncode == 0
            if engine_proc.returncode == 0:
                vtk = convert_anim_to_vtk(job, log_lines)
                vtk_path = str(vtk)
                modulus = extract_modulus(case, vtk)
                target = refs[case.target_key]
                rel_err = abs(modulus - target) / target
                passed = rel_err <= HT_REL_TOL
        all_pass = all_pass and passed
        rows.append(
            {
                "case": case.name,
                "quantity": case.target_key,
                "solver_Pa": modulus,
                "target_Pa": refs[case.target_key],
                "relative_error_pct": "" if rel_err == "" else 100.0 * float(rel_err),
                "starter_rc": starter_proc.returncode,
                "engine_rc": engine_rc,
                "vtk_path": vtk_path,
                "verdict": "PASS" if passed else "FAIL",
            }
        )

    metrics = {
        "mesh": f"{NX}x{NY}x{NZ} HEXA8 voxel cell",
        "fiber_volume_fraction": mesh.vf_actual,
        "fiber_bricks": len(mesh.fiber_bricks),
        "matrix_bricks": len(mesh.matrix_bricks),
        "starter_all_ok": all_started,
        "engine_all_ok": all_engine,
        "modulus_gate_evaluated": True,
        "modulus_gate_pass": all_pass,
        "max_relative_error_pct": max(
            float(row["relative_error_pct"]) for row in rows if row["relative_error_pct"] != ""
        ),
    }
    return metrics, rows, log_lines


def _format_fail_rows(rows: list[dict[str, object]]) -> list[str]:
    out: list[str] = []
    for row in rows:
        if row["verdict"] == "FAIL":
            if row["solver_Pa"] == "":
                out.append(f"- {row['case']}: solver did not complete; starter_rc={row['starter_rc']}, engine_rc={row['engine_rc']}.")
            else:
                out.append(
                    f"- {row['case']}: FEM `{float(row['solver_Pa']) / 1.0e9:.3f} GPa`, "
                    f"target `{float(row['target_Pa']) / 1.0e9:.3f} GPa`, "
                    f"error `{float(row['relative_error_pct']):.3f}%`."
                )
    return out


def write_outputs(metrics: dict[str, object], rows: list[dict[str, object]], log_lines: list[str], wall_s: float) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    with (RESULTS_DIR / "effective_moduli.csv").open("w", newline="", encoding="utf-8") as f:
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
    verdict = "PASS" if bool(metrics["modulus_gate_pass"]) else "FAIL"
    results = {
        "stage": 10,
        "verdict": verdict,
        "metrics": {**metrics, "wall_clock_s": wall_s},
        "reference": {
            "model": "rule of mixtures for E1; Halpin-Tsai zeta=2 for E2 and zeta=1 for G12",
            "fiber": {"E_Pa": E_F, "nu": NU_F, "rho": RHO_F},
            "matrix": {"E_Pa": E_M, "nu": NU_M, "rho": RHO_M},
        },
        "tolerance": {"modulus_relative_error": HT_REL_TOL},
        "git_sha": git_sha,
    }
    (RESULTS_DIR / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    (FIGURES_DIR / "stage10_ud_mesoscale.typ").write_text(
        "\n".join(
            [
                '#set page(width: 170mm, height: auto, margin: 10mm)',
                '#let rows = csv("../results/effective_moduli.csv")',
                '#text(size: 12pt, weight: "bold")[Stage 10 UD mesoscale moduli]',
                '#v(5pt)',
                '#table(',
                '  columns: (34mm, 34mm, 34mm, 28mm, 24mm),',
                '  stroke: rgb("#5C5C5C"),',
                '  [Case], [FEM Pa], [Target Pa], [Error %], [Verdict],',
                '  ..rows.map(r => ([#r.at(0)], [#r.at(2)], [#r.at(3)], [#r.at(4)], [#r.at(8)])).flatten(),',
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
                    "# Stage 10 Blocker - UD mesoscale modulus mismatch",
                    "",
                    "Author: J.C. Vaught",
                    "",
                    "The isotropic fiber/matrix HEXA8 mesoscale cell starts, runs, converts to VTK, and post-processes, but one or more effective moduli missed the 5 percent gate.",
                    "",
                    *(_format_fail_rows(rows) or ["- No completed solver rows were available."]),
                    "",
                    f"Verdict: `{verdict}`.",
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
    verdict = "PASS" if bool(metrics["modulus_gate_pass"]) else "FAIL"
    print(json.dumps({"stage": 10, "verdict": verdict, "wall_clock_s": wall_s}, indent=2))
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
