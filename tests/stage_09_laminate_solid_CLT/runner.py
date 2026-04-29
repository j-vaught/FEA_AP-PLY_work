"""Stage 09 - solid laminate CLT verification.

Author: J.C. Vaught

Composite plies use the verified LAW12 + TYPE6/SOL_ORTH row from
references/openradioss_law_compatibility_matrix.md. Per-ply orientation is
bound to references/openradioss_orientation_convention.md using TYPE6
Ip=3, Iorth=0, Phi=theta_ply. Three grouped affine strain cases recover the
laminate A-matrix from VTK stress resultants.
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
CARD_PATH = ROOT_DIR / "references" / "material_cards" / "im7_8552.json"

OR_ROOT = Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"
ENGINE = OR_ROOT / "exec" / "engine_linux64_gf"
ANIM_TO_VTK = OR_ROOT / "exec" / "anim_to_vtk_linux64_gf"
N_THREADS = int(os.environ.get("RAD_NT", "16"))

PLY_T = 0.18e-3
LX = 0.020
LY = 0.020
NX_INPLANE = 10
NY_INPLANE = 10
STRAIN = 5.0e-4
RUN_TIME = 2.0e-4
LAYUPS = {
    "A_crossply": (0.0, 90.0, 90.0, 0.0),
    "B_quasiiso": (0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0),
}
COMPONENTS = ("xx", "yy", "xy")


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
    zt: float
    zc: float
    s12: float
    s13: float
    s23: float


@dataclass(frozen=True)
class MeshData:
    nodes: list[str]
    bricks: list[tuple[int, int, list[int]]]
    groups: dict[str, list[int]]


@dataclass(frozen=True)
class LoadCase:
    name: str
    column: int
    denominator: float


LOAD_CASES = (
    LoadCase("exx", 0, STRAIN),
    LoadCase("eyy", 1, STRAIN),
    LoadCase("gxy", 2, STRAIN),
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
        zt=float(strength["Zt_Pa"]),
        zc=float(strength["Zc_Pa"]),
        s12=float(strength["S12_Pa"]),
        s13=float(strength["S13_Pa"]),
        s23=float(strength["S23_Pa"]),
    )


def q_matrix(mat: Material) -> np.ndarray:
    nu21 = mat.nu12 * mat.e2 / mat.e1
    denom = 1.0 - mat.nu12 * nu21
    q = np.zeros((3, 3), dtype=float)
    q[0, 0] = mat.e1 / denom
    q[1, 1] = mat.e2 / denom
    q[0, 1] = q[1, 0] = mat.nu12 * mat.e2 / denom
    q[2, 2] = mat.g12
    return q


def qbar(q: np.ndarray, theta_deg: float) -> np.ndarray:
    th = math.radians(theta_deg)
    c = math.cos(th)
    s = math.sin(th)
    c2 = c * c
    s2 = s * s
    c4 = c2 * c2
    s4 = s2 * s2
    q11, q22, q12, q66 = q[0, 0], q[1, 1], q[0, 1], q[2, 2]
    out = np.zeros((3, 3), dtype=float)
    out[0, 0] = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * s4
    out[1, 1] = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * c4
    out[0, 1] = out[1, 0] = (q11 + q22 - 4.0 * q66) * s2 * c2 + q12 * (s4 + c4)
    out[2, 2] = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * s2 * c2 + q66 * (s4 + c4)
    out[0, 2] = out[2, 0] = (q11 - q12 - 2.0 * q66) * s * c**3 + (q12 - q22 + 2.0 * q66) * s**3 * c
    out[1, 2] = out[2, 1] = (q11 - q12 - 2.0 * q66) * s**3 * c + (q12 - q22 + 2.0 * q66) * s * c**3
    return out


def a_matrix(mat: Material, thetas: tuple[float, ...]) -> np.ndarray:
    q = q_matrix(mat)
    a = np.zeros((3, 3), dtype=float)
    for theta in thetas:
        a += qbar(q, theta) * PLY_T
    return a


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


def law12_block(mat: Material, mat_id: int) -> list[str]:
    nu31 = mat.nu13 * mat.e3 / mat.e1
    return [
        f"/MAT/LAW12/{mat_id}",
        f"IM7_8552_ply_{mat_id}",
        "#              RHO_I",
        fmt_f(mat.rho),
        "#                E11                 E22                 E33",
        fmt_f(mat.e1, mat.e2, mat.e3),
        "#               NU12                NU23                NU31",
        fmt_f(mat.nu12, mat.nu23, nu31),
        "#                G12                 G23                 G31",
        fmt_f(mat.g12, mat.g23, mat.g13),
        "#           SIGMA_T1            SIGMA_T2            SIGMA_T3               DELTA",
        fmt_f(mat.xt, mat.yt, mat.zt, 0.05),
        "#                  B                   n                fmax               Wpref",
        fmt_f(1.0, 1.0, 1.0, 1.0),
        "#          sigma_1yt           sigma_2yt           sigma_1yc           sigma_2yc",
        fmt_f(mat.xt, mat.yt, mat.xc, mat.yc),
        "#         sigma_12yt          sigma_12yc          sigma_23yt          sigma_23yc",
        fmt_f(mat.s12, mat.s12, mat.s23, mat.s23),
        "#          sigma_3yt           sigma_3yc          sigma_13yt          sigma_13yc",
        fmt_f(mat.zt, mat.zc, mat.s13, mat.s13),
        "#              alpha                  Ef                   c          EPS_RATE_0   STRFLAG",
        fmt_f(0.0, 0.0, 0.0, 0.0) + fmt_i(1),
    ]


def type6_property(prop_id: int, theta_deg: float) -> list[str]:
    return [
        f"/PROP/TYPE6/{prop_id}",
        f"ply_{prop_id}_type6_phi_ip3_theta_{theta_deg:+.0f}",
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


def stack_mesh(nplies: int) -> MeshData:
    node_id: dict[tuple[int, int, int], int] = {}
    nodes = ["/NODE"]
    nid = 1
    for k in range(nplies + 1):
        for j in range(NY_INPLANE + 1):
            for i in range(NX_INPLANE + 1):
                node_id[(i, j, k)] = nid
                x = LX * i / NX_INPLANE
                y = LY * j / NY_INPLANE
                nodes.append(f"{nid:10d}{x:20.12g}{y:20.12g}{PLY_T * k:20.12g}")
                nid += 1

    bricks: list[tuple[int, int, list[int]]] = []
    eid = 1
    for k in range(nplies):
        for j in range(NY_INPLANE):
            for i in range(NX_INPLANE):
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
                bricks.append((k + 1, eid, conn))
                eid += 1

    all_nodes = sorted(node_id.values())
    groups = {
        "all": all_nodes,
        "x0": sorted(node_id[(0, j, k)] for k in range(nplies + 1) for j in range(NY_INPLANE + 1)),
        "x1": sorted(node_id[(NX_INPLANE, j, k)] for k in range(nplies + 1) for j in range(NY_INPLANE + 1)),
        "y0": sorted(node_id[(i, 0, k)] for k in range(nplies + 1) for i in range(NX_INPLANE + 1)),
        "y1": sorted(node_id[(i, NY_INPLANE, k)] for k in range(nplies + 1) for i in range(NX_INPLANE + 1)),
        "z_anchor": [node_id[(0, 0, 0)]],
    }
    return MeshData(nodes=nodes, bricks=bricks, groups=groups)


def group_block(group_id: int, name: str, node_ids: list[int]) -> list[str]:
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for i in range(0, len(node_ids), 10):
        lines.append(fmt_i(*node_ids[i : i + 10]))
    return lines


def bcs_and_loads(case: LoadCase, mesh: MeshData) -> list[str]:
    group_ids = {"x0": 100, "x1": 101, "y0": 102, "y1": 103, "all": 104, "z_anchor": 105}
    lines: list[str] = []
    if case.name == "exx":
        lines.extend(
            [
                "/BCS/1",
                "x0_fixed_x",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   100 000{0:10d}{group_ids['x0']:10d}",
                "/BCS/2",
                "all_fixed_y",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   010 000{0:10d}{group_ids['all']:10d}",
                "/IMPDISP/1",
                "x1_exx",
                "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
                f"{1:10d}{'X':>10}{0:10d}{0:10d}{group_ids['x1']:10d}{0:10d}{0:10d}",
                "#            Scale_x             Scale_y              Tstart               Tstop",
                fmt_f(1.0, STRAIN * LX, 0.0, 0.0),
            ]
        )
    elif case.name == "eyy":
        lines.extend(
            [
                "/BCS/1",
                "y0_fixed_y",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   010 000{0:10d}{group_ids['y0']:10d}",
                "/BCS/2",
                "all_fixed_x",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   100 000{0:10d}{group_ids['all']:10d}",
                "/IMPDISP/1",
                "y1_eyy",
                "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
                f"{1:10d}{'Y':>10}{0:10d}{0:10d}{group_ids['y1']:10d}{0:10d}{0:10d}",
                "#            Scale_x             Scale_y              Tstart               Tstop",
                fmt_f(1.0, STRAIN * LY, 0.0, 0.0),
            ]
        )
    elif case.name == "gxy":
        lines.extend(
            [
                "/BCS/1",
                "y0_fixed_x",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   100 000{0:10d}{group_ids['y0']:10d}",
                "/BCS/2",
                "all_fixed_y",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   010 000{0:10d}{group_ids['all']:10d}",
                "/IMPDISP/1",
                "y1_gxy_x",
                "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
                f"{1:10d}{'X':>10}{0:10d}{0:10d}{group_ids['y1']:10d}{0:10d}{0:10d}",
                "#            Scale_x             Scale_y              Tstart               Tstop",
                fmt_f(1.0, STRAIN * LY, 0.0, 0.0),
            ]
        )
    else:
        raise ValueError(case.name)

    lines.extend(
        [
            "/BCS/3",
            "z_anchor",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   001 000{0:10d}{group_ids['z_anchor']:10d}",
            "/FUNCT/1",
            "unit_ramp",
            fmt_f(0.0, 0.0),
            fmt_f(RUN_TIME, 1.0),
        ]
    )
    for name, gid in group_ids.items():
        lines.extend(group_block(gid, name, mesh.groups[name]))
    return lines


def write_starter(layup_name: str, thetas: tuple[float, ...], case: LoadCase, mat: Material) -> Path:
    job = f"stage09_{layup_name}_{case.name}_law12_type6"
    mesh = stack_mesh(len(thetas))
    lines = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        "      2023         0",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 09 {layup_name} {case.name} LAW12 TYPE6",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 4) + f"{0:20d}{2:40d}",
    ]
    for ply, theta in enumerate(thetas, start=1):
        lines.extend(law12_block(mat, ply))
        lines.extend(type6_property(ply, theta))
    lines.extend(mesh.nodes)
    for ply, _theta in enumerate(thetas, start=1):
        lines.extend([f"/PART/{ply}", f"ply_{ply:02d}", fmt_i(ply, ply, 0), f"/BRICK/{ply}"])
        for brick_ply, eid, conn in mesh.bricks:
            if brick_ply == ply:
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


def extract_resultant(vtk_path: Path) -> np.ndarray:
    import pyvista as pv  # type: ignore[import-not-found]

    grid = pv.read(str(vtk_path))
    stress = np.asarray(_cell_array(grid, "Strs"), dtype=float)
    if stress.shape[1] >= 9:
        sigma_xx = stress[:, 0]
        sigma_yy = stress[:, 4]
        sigma_xy = 0.5 * (stress[:, 1] + stress[:, 3])
    else:
        sigma_xx = stress[:, 0]
        sigma_yy = stress[:, 1]
        sigma_xy = stress[:, 3]
    scale = PLY_T / (NX_INPLANE * NY_INPLANE)
    return np.array(
        [
            float(np.sum(sigma_xx)) * scale,
            float(np.sum(sigma_yy)) * scale,
            float(np.sum(sigma_xy)) * scale,
        ],
        dtype=float,
    )


def run_case(layup_name: str, thetas: tuple[float, ...], case: LoadCase, mat: Material, log_lines: list[str]) -> tuple[np.ndarray | None, dict[str, object]]:
    starter = write_starter(layup_name, thetas, case, mat)
    job = starter.name.removesuffix("_0000.rad")
    engine = write_engine(job)
    starter_proc = run_cmd([str(STARTER), "-i", starter.name, "-nt", str(N_THREADS)], RUNS_DIR, log_lines)
    row = {"starter_rc": starter_proc.returncode, "engine_rc": "", "vtk_path": ""}
    if starter_proc.returncode != 0:
        return None, row
    engine_proc = run_cmd([str(ENGINE), "-i", engine.name, "-nt", str(N_THREADS)], RUNS_DIR, log_lines)
    row["engine_rc"] = engine_proc.returncode
    if engine_proc.returncode != 0:
        return None, row
    vtk = convert_anim_to_vtk(job, log_lines)
    row["vtk_path"] = str(vtk)
    return extract_resultant(vtk), row


def component_pass(reference: float, solver: float, scale: float) -> tuple[bool, float]:
    if abs(reference) > 1.0e-9 * scale:
        err = abs(solver - reference) / abs(reference)
        return err <= 0.02, 100.0 * err
    err = abs(solver) / scale
    return err <= 0.02, 100.0 * err


def run_stage() -> tuple[dict[str, object], list[dict[str, object]], list[str]]:
    mat = load_material()
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    log_lines = [f"material_card={CARD_PATH}", f"starter={STARTER}", f"engine={ENGINE}"]
    rows: list[dict[str, object]] = []
    all_started = True
    all_engine = True
    all_pass = True
    max_error = 0.0

    for layup_name, thetas in LAYUPS.items():
        ref_a = a_matrix(mat, thetas)
        fem_a = np.full((3, 3), np.nan)
        case_info: dict[str, dict[str, object]] = {}
        for case in LOAD_CASES:
            resultant, info = run_case(layup_name, thetas, case, mat, log_lines)
            case_info[case.name] = info
            all_started = all_started and info["starter_rc"] == 0
            all_engine = all_engine and info["engine_rc"] == 0
            if resultant is not None:
                fem_a[:, case.column] = resultant / case.denominator

        scale = float(np.max(np.abs(ref_a)))
        for i, row_name in enumerate(COMPONENTS):
            for j, col_name in enumerate(COMPONENTS):
                solver = fem_a[i, j]
                passed = False
                rel_err_pct: float | str = ""
                if not math.isnan(float(solver)):
                    passed, rel_err_pct = component_pass(float(ref_a[i, j]), float(solver), scale)
                    max_error = max(max_error, float(rel_err_pct))
                all_pass = all_pass and passed
                rows.append(
                    {
                        "layup": layup_name,
                        "component": f"A_{row_name}{col_name}",
                        "reference_N_per_m": ref_a[i, j],
                        "solver_N_per_m": "" if math.isnan(float(solver)) else solver,
                        "relative_error_pct": rel_err_pct,
                        "verdict": "PASS" if passed else "FAIL",
                    }
                )

        for case_name, info in case_info.items():
            rows.append(
                {
                    "layup": layup_name,
                    "component": f"case_{case_name}",
                    "reference_N_per_m": "",
                    "solver_N_per_m": "",
                    "relative_error_pct": "",
                    "verdict": f"starter_rc={info['starter_rc']}; engine_rc={info['engine_rc']}; vtk={info['vtk_path']}",
                }
            )

    metrics = {
        "canonical_material_property": "LAW12 + TYPE6/SOL_ORTH",
        "matrix_evidence": "LAW12 row: TYPE6/SOL_ORTH solid OK; TYPE14 solid = B3047",
        "orientation_evidence": "references/openradioss_orientation_convention.md per-ply property angle: Ip=3, Iorth=0, Phi=theta_ply",
        "layup_count": len(LAYUPS),
        "load_case_count": len(LOAD_CASES),
        "mesh_evidence": f"{NX_INPLANE}x{NY_INPLANE} in-plane HEXA8 cells with one element through each ply",
        "starter_all_ok": all_started,
        "engine_all_ok": all_engine,
        "clt_A_gate_evaluated": True,
        "clt_A_gate_pass": all_pass,
        "max_component_error_pct": max_error,
    }
    return metrics, rows, log_lines


def _format_fail_rows(rows: list[dict[str, object]]) -> list[str]:
    out: list[str] = []
    for row in rows:
        if row["verdict"] == "FAIL":
            out.append(
                f"- {row['layup']} {row['component']}: FEM `{float(row['solver_N_per_m']):.6e}` N/m, "
                f"CLT `{float(row['reference_N_per_m']):.6e}` N/m, error `{float(row['relative_error_pct']):.3f}%`."
            )
    return out


def write_outputs(metrics: dict[str, object], rows: list[dict[str, object]], log_lines: list[str], wall_s: float) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    with (RESULTS_DIR / "timeseries.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    mat = load_material()
    for name, thetas in LAYUPS.items():
        np.savetxt(
            RESULTS_DIR / f"cltA_{name}_reference.csv",
            a_matrix(mat, thetas),
            delimiter=",",
            header=f"CLT A-matrix [N/m] for layup {name}, plies {thetas} deg, t_ply = {PLY_T} m",
        )

    git_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=str(ROOT_DIR),
        text=True,
        stdout=subprocess.PIPE,
        check=False,
    ).stdout.strip()
    verdict = "PASS" if bool(metrics["clt_A_gate_pass"]) else "FAIL"
    results = {
        "stage": 9,
        "verdict": verdict,
        "metrics": {**metrics, "wall_clock_s": wall_s},
        "reference": {
            "material_card": str(CARD_PATH),
            "ply_thickness_m": PLY_T,
            "layups": {name: list(thetas) for name, thetas in LAYUPS.items()},
        },
        "tolerance": {"A_matrix_component_relative_error": 0.02},
        "git_sha": git_sha,
    }
    (RESULTS_DIR / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    (FIGURES_DIR / "stage09_clt_probe.typ").write_text(
        "\n".join(
            [
                '#set page(width: 175mm, height: auto, margin: 10mm)',
                '#let rows = csv("../results/timeseries.csv")',
                '#text(size: 12pt, weight: "bold")[Stage 09 laminate CLT sweep]',
                '#v(5pt)',
                '#table(',
                '  columns: (34mm, 24mm, 34mm, 34mm, 24mm),',
                '  stroke: rgb("#5C5C5C"),',
                '  [Layup], [Component], [CLT N/m], [FEM N/m], [Verdict],',
                '  ..rows.map(r => ([#r.at(0)], [#r.at(1)], [#r.at(2)], [#r.at(3)], [#r.at(5)])).flatten(),',
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
        fail_rows = _format_fail_rows(rows)
        blocker.write_text(
            "\n".join(
                [
                    "# Stage 09 Blocker - laminate CLT mismatch",
                    "",
                    "Author: J.C. Vaught",
                    "",
                    "The laminate stack uses `LAW12 + TYPE6/SOL_ORTH` with per-ply `Ip=3`, `Iorth=0`, `Phi=theta_ply` on `/PROP/TYPE6`.",
                    "",
                    "Observed failing A-matrix components from `results/timeseries.csv`:",
                    "",
                    *(fail_rows or ["- Solver did not complete one or more load cases."]),
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
    verdict = "PASS" if bool(metrics["clt_A_gate_pass"]) else "FAIL"
    print(json.dumps({"stage": 9, "verdict": verdict, "wall_clock_s": wall_s}, indent=2))
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
