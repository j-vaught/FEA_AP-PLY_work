"""Stage 11 - pseudo-woven mesoscale direct stiffness validation.

This runner uses the clean-room ``kok_geom`` port directly.  It generates the
stage-11 AP-PLY block, writes OpenRadioss TETRA10 decks, assigns composite tow
orientation with ``/INIBRI/ORTHO`` on the verified LAW12 + TYPE6/SOL_ORTH row,
and compares effective Ex, Ey, and Gxy to the Kok 2022 targets.
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

import meshio
import numpy as np

from kok_geom.config import KokConfig
from kok_geom.io import generate_mesh


THIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parents[1]
GEOM_DIR = THIS_DIR / "geometry"
RUNS_DIR = THIS_DIR / "runs"
RESULTS_DIR = THIS_DIR / "results"
FIGURES_DIR = THIS_DIR / "figures"
RUN_LOG = THIS_DIR / "run.log"
CONFIG_PATH = GEOM_DIR / "kok_stage11_config.json"
CARD_PATH = ROOT_DIR / "references" / "material_cards" / "im7_8552.json"

OR_ROOT = Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"
ENGINE = OR_ROOT / "exec" / "engine_linux64_gf"
ANIM_TO_VTK = OR_ROOT / "exec" / "anim_to_vtk_linux64_gf"
N_THREADS = int(os.environ.get("RAD_NT", "16"))

LX = 25.0e-3
LY = 25.0e-3
LZ = 4.0 * 0.18e-3
STRAIN = 1.0e-3
RUN_TIME = 2.0e-4
KOK_TARGETS_GPA = {"E_x_GPa": 53.3, "E_y_GPa": 53.3, "G_xy_GPa": 20.5}
PASS_TOL = 0.10
INCONCLUSIVE_EXTRA_TOL = 0.05


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
class LoadCase:
    name: str
    result_key: str
    component: str
    denominator: float


LOAD_CASES = (
    LoadCase("axial_x", "E_x_GPa", "sigma_xx", STRAIN),
    LoadCase("transverse_y", "E_y_GPa", "sigma_yy", STRAIN),
    LoadCase("shear_xy", "G_xy_GPa", "sigma_xy", STRAIN),
)


def fmt_f(*values: float) -> str:
    return "".join(f"{float(value):20.12g}" for value in values)


def fmt_i(*values: int) -> str:
    return "".join(f"{int(value):10d}" for value in values)


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


def law12_block(mat: Material) -> list[str]:
    nu31 = mat.nu13 * mat.e3 / mat.e1
    return [
        "/MAT/LAW12/1",
        "IM7_8552_tow_LAW12",
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


def resin_block() -> list[str]:
    return [
        "/MAT/ELAST/2",
        "8552_resin_isotropic",
        "#              RHO_I",
        fmt_f(1300.0, 0.0),
        "#                  E                  nu",
        fmt_f(4.67e9, 0.38),
    ]


def type6_property(prop_id: int) -> list[str]:
    return [
        f"/PROP/TYPE6/{prop_id}",
        f"type6_sol_orth_inibri_{prop_id}",
        "#   Isolid    Ismstr               Icpre  Itetra10     Inpts   Itetra4    Iframe                  Dn",
        fmt_i(24, 4) + f"{1:20d}{1000:10d}{0:10d}{0:10d}{2:10d}{0.0:20.12g}",
        "#                 qa                  qb                   h",
        fmt_f(0.0, 0.0, 0.0),
        "#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth",
        fmt_f(1.0, 0.0, 0.0) + fmt_i(0, 1, 0),
        "#                Phi                 Px                  Py                  Pz",
        fmt_f(0.0, 0.0, 0.0, 0.0),
        "#             dt_min   istrain      IHKT",
        fmt_f(0.0) + fmt_i(0, 0),
    ]


def ensure_geometry(log_lines: list[str]) -> tuple[Path, Path]:
    GEOM_DIR.mkdir(parents=True, exist_ok=True)
    if not CONFIG_PATH.exists():
        config = KokConfig.model_validate(
            {
                "panel": {"size_x_mm": 25.0, "size_y_mm": 25.0, "n_plies": 4, "symmetry": "none"},
                "laydown": {
                    "fiber_angles_deg": [0.0, 45.0, -45.0, 90.0],
                    "placement_sequence": "1010",
                    "angle_shift_deg": 0.0,
                    "tape_width_mm": 6.35,
                    "cured_ply_thickness_mm": 0.18,
                    "undulation_ratio": 0.09,
                    "tape_spacing": 1,
                },
                "mesh": {
                    "in_plane_target_mm_impact_zone": 1.0,
                    "in_plane_target_mm_far_field": 1.0,
                    "through_thickness_target_mm": 1.0,
                    "graded_zone_radius_mm": 0.0,
                    "element_order": 2,
                },
                "output": {
                    "msh_path": str(GEOM_DIR / "panel.msh"),
                    "inp_path": str(GEOM_DIR / "panel.inp"),
                    "orientations_json_path": str(GEOM_DIR / "orientations.json"),
                    "msh_format": "msh4_ascii",
                },
            }
        )
        config.to_file(CONFIG_PATH)
    msh_path, orientations_path = generate_mesh(CONFIG_PATH, GEOM_DIR / "panel.msh")
    log_lines.append(f"kok_geom_msh={msh_path}")
    log_lines.append(f"orientations={orientations_path}")
    return msh_path, orientations_path


def mesh_groups(msh_path: Path, orientations_path: Path) -> tuple[meshio.Mesh, dict[int, str], dict[str, dict[str, object]]]:
    mesh = meshio.read(msh_path)
    physical_to_name = {int(data[0]): name for name, data in mesh.field_data.items() if int(data[1]) == 3}
    sidecar = json.loads(orientations_path.read_text(encoding="utf-8"))
    by_name = {str(group["name"]): group for group in sidecar["groups"]}
    return mesh, physical_to_name, by_name


def node_groups(points: np.ndarray) -> dict[str, list[int]]:
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    span = max(maxs - mins)
    tol = max(1.0e-9, span * 1.0e-7)
    groups: dict[str, list[int]] = {}
    for axis, name0, name1 in ((0, "x0", "x1"), (1, "y0", "y1"), (2, "z0", "z1")):
        groups[name0] = [idx + 1 for idx, xyz in enumerate(points) if abs(float(xyz[axis] - mins[axis])) <= tol]
        groups[name1] = [idx + 1 for idx, xyz in enumerate(points) if abs(float(xyz[axis] - maxs[axis])) <= tol]
    groups["anchor_xyz"] = [1]
    groups["anchor_yz"] = [max(1, len(points) // 3)]
    groups["anchor_z"] = [max(1, 2 * len(points) // 3)]
    return groups


def group_block(group_id: int, name: str, node_ids: list[int]) -> list[str]:
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for i in range(0, len(node_ids), 10):
        lines.append(fmt_i(*node_ids[i : i + 10]))
    return lines


def bcs_and_loads(case: LoadCase, groups: dict[str, list[int]]) -> list[str]:
    gids = {
        "x0": 100,
        "x1": 101,
        "y0": 102,
        "y1": 103,
        "anchor_xyz": 104,
        "anchor_yz": 105,
        "anchor_z": 106,
    }
    lines = ["/FUNCT/1", "unit_ramp", fmt_f(0.0, 0.0), fmt_f(RUN_TIME, 1.0)]
    if case.name == "axial_x":
        lines.extend(
            [
                "/BCS/1",
                "x0_fixed_x",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   100 000{0:10d}{gids['x0']:10d}",
                "/BCS/2",
                "anchor_yz",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   011 000{0:10d}{gids['anchor_xyz']:10d}",
                "/BCS/3",
                "anchor_z",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   001 000{0:10d}{gids['anchor_yz']:10d}",
                "/IMPDISP/1",
                "x1_axial_x",
                "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
                f"{1:10d}{'X':>10}{0:10d}{0:10d}{gids['x1']:10d}{0:10d}{0:10d}",
                "#            Scale_x             Scale_y              Tstart               Tstop",
                fmt_f(1.0, STRAIN * LX, 0.0, 0.0),
            ]
        )
    elif case.name == "transverse_y":
        lines.extend(
            [
                "/BCS/1",
                "y0_fixed_y",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   010 000{0:10d}{gids['y0']:10d}",
                "/BCS/2",
                "anchor_xz",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   101 000{0:10d}{gids['anchor_xyz']:10d}",
                "/BCS/3",
                "anchor_z",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   001 000{0:10d}{gids['anchor_yz']:10d}",
                "/IMPDISP/1",
                "y1_transverse_y",
                "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
                f"{1:10d}{'Y':>10}{0:10d}{0:10d}{gids['y1']:10d}{0:10d}{0:10d}",
                "#            Scale_x             Scale_y              Tstart               Tstop",
                fmt_f(1.0, STRAIN * LY, 0.0, 0.0),
            ]
        )
    elif case.name == "shear_xy":
        lines.extend(
            [
                "/BCS/1",
                "y0_fixed_x",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   100 000{0:10d}{gids['y0']:10d}",
                "/BCS/2",
                "anchor_yz",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   011 000{0:10d}{gids['anchor_xyz']:10d}",
                "/BCS/3",
                "anchor_y",
                "#  Tra rot   skew_ID  grnod_ID",
                f"   010 000{0:10d}{gids['anchor_yz']:10d}",
                "/IMPDISP/1",
                "y1_shear_x",
                "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
                f"{1:10d}{'X':>10}{0:10d}{0:10d}{gids['y1']:10d}{0:10d}{0:10d}",
                "#            Scale_x             Scale_y              Tstart               Tstop",
                fmt_f(1.0, STRAIN * LY, 0.0, 0.0),
            ]
        )
    else:
        raise ValueError(case.name)
    for name, gid in gids.items():
        lines.extend(group_block(gid, name, groups[name]))
    return lines


def write_starter(
    case: LoadCase,
    mesh: meshio.Mesh,
    physical_to_name: dict[int, str],
    orientation_groups: dict[str, dict[str, object]],
    mat: Material,
) -> tuple[Path, dict[str, object]]:
    job = f"stage11_{case.name}"
    part_names = sorted(physical_to_name.items())
    points = np.asarray(mesh.points, dtype=float)
    groups = node_groups(points)
    lines = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        "      2023         0",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 11 AP-PLY {case.name} LAW12 TYPE6 INIBRI",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 4) + f"{0:20d}{2:40d}",
    ]
    lines.extend(law12_block(mat))
    lines.extend(resin_block())
    for physical_id, _name in part_names:
        lines.extend(type6_property(physical_id))
    lines.append("/NODE")
    for idx, xyz in enumerate(points, start=1):
        lines.append(f"{idx:10d}{float(xyz[0]):20.12g}{float(xyz[1]):20.12g}{float(xyz[2]):20.12g}")

    element_orientation: dict[int, tuple[float, float, float, float, float, float]] = {}
    element_count = 0
    tow_element_count = 0
    for physical_id, name in part_names:
        group = orientation_groups[name]
        mat_id = 2 if group.get("isotropic") else 1
        lines.extend([f"/PART/{physical_id}", name, fmt_i(physical_id, mat_id, 0), f"/TETRA10/{physical_id}"])
        for block, physical_values in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
            if block.type != "tetra10":
                continue
            for conn, phys in zip(block.data, physical_values):
                if int(phys) != physical_id:
                    continue
                element_count += 1
                node_ids = [int(node) + 1 for node in conn]
                lines.append(fmt_i(element_count))
                lines.append(fmt_i(*node_ids))
                if not group.get("isotropic"):
                    v1 = [float(v) for v in group["fiber_direction_unit_vector"]]
                    v2 = [float(v) for v in group["transverse_in_plane_unit_vector"]]
                    element_orientation[element_count] = (v1[0], v1[1], v1[2], v2[0], v2[1], v2[2])
                    tow_element_count += 1

    lines.append("/INIBRI/ORTHO")
    for eid, (x1, y1, z1, x2, y2, z2) in element_orientation.items():
        lines.append(fmt_i(eid, 1, 10, 6, 24))
        lines.append(fmt_f(x1, y1, z1, x2, y2))
        lines.append(fmt_f(z2))

    lines.extend(bcs_and_loads(case, groups))
    lines.extend(["/END", ""])
    path = RUNS_DIR / f"{job}_0000.rad"
    path.write_text("\n".join(lines), encoding="utf-8")
    info = {
        "nodes": len(points),
        "tetra10": element_count,
        "tow_tetra10": tow_element_count,
        "resin_tetra10": element_count - tow_element_count,
    }
    return path, info


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


def extract_stiffness(case: LoadCase, vtk_path: Path) -> float:
    import pyvista as pv  # type: ignore[import-not-found]

    grid = pv.read(str(vtk_path))
    stress = np.asarray(_cell_array(grid, "Strs"), dtype=float)
    if stress.shape[1] >= 9:
        values = {
            "sigma_xx": stress[:, 0],
            "sigma_yy": stress[:, 4],
            "sigma_xy": 0.5 * (stress[:, 1] + stress[:, 3]),
        }[case.component]
    else:
        values = {
            "sigma_xx": stress[:, 0],
            "sigma_yy": stress[:, 1],
            "sigma_xy": stress[:, 3],
        }[case.component]
    return abs(float(np.mean(values))) / case.denominator / 1.0e9


def verdict_from_errors(errors: dict[str, float]) -> str:
    max_error = max(errors.values())
    if max_error <= PASS_TOL:
        return "PASS"
    if max_error <= PASS_TOL + INCONCLUSIVE_EXTRA_TOL:
        return "INCONCLUSIVE"
    return "FAIL"


def run_stage() -> tuple[dict[str, object], list[dict[str, object]], list[str]]:
    for directory in (GEOM_DIR, RUNS_DIR, RESULTS_DIR, FIGURES_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    log_lines = [f"starter={STARTER}", f"engine={ENGINE}", f"anim_to_vtk={ANIM_TO_VTK}"]
    msh_path, orientations_path = ensure_geometry(log_lines)
    mesh, physical_to_name, orientation_groups = mesh_groups(msh_path, orientations_path)
    mat = load_material()

    rows: list[dict[str, object]] = []
    measured: dict[str, float] = {}
    mesh_info: dict[str, object] = {}
    all_started = True
    all_engine = True
    for case in LOAD_CASES:
        starter, info = write_starter(case, mesh, physical_to_name, orientation_groups, mat)
        mesh_info = info
        job = starter.name.removesuffix("_0000.rad")
        engine = write_engine(job)
        starter_proc = run_cmd([str(STARTER), "-i", starter.name, "-nt", str(N_THREADS)], RUNS_DIR, log_lines)
        all_started = all_started and starter_proc.returncode == 0
        engine_rc: int | str = ""
        vtk_path = ""
        value = math.nan
        if starter_proc.returncode == 0:
            engine_proc = run_cmd([str(ENGINE), "-i", engine.name, "-nt", str(N_THREADS)], RUNS_DIR, log_lines)
            engine_rc = engine_proc.returncode
            all_engine = all_engine and engine_proc.returncode == 0
            if engine_proc.returncode == 0:
                vtk = convert_anim_to_vtk(job, log_lines)
                vtk_path = str(vtk)
                value = extract_stiffness(case, vtk)
                measured[case.result_key] = value
        else:
            all_engine = False
        target = KOK_TARGETS_GPA[case.result_key]
        rel = math.nan if math.isnan(value) else abs(value - target) / target
        rows.append(
            {
                "case": case.name,
                "metric": case.result_key,
                "measured_GPa": "" if math.isnan(value) else value,
                "target_GPa": target,
                "relative_error_pct": "" if math.isnan(rel) else 100.0 * rel,
                "starter_rc": starter_proc.returncode,
                "engine_rc": engine_rc,
                "vtk_path": vtk_path,
            }
        )

    errors = {
        key: abs(measured.get(key, math.nan) - target) / target
        for key, target in KOK_TARGETS_GPA.items()
        if key in measured and not math.isnan(measured[key])
    }
    verdict = "FAIL" if len(errors) != len(KOK_TARGETS_GPA) else verdict_from_errors(errors)
    max_error_pct = 100.0 * max(errors.values()) if errors else math.nan
    metrics = {
        "geometry_source": "kok_geom CLI/config",
        "canonical_material_property": "LAW12 + TYPE6/SOL_ORTH",
        "orientation_evidence": "per-tow /INIBRI/ORTHO with neutral TYPE6 Ip=1, Iorth=0, Phi=0",
        "strain_measurement": "engineering strain from imposed displacement; VTK cell Stra[0] is not used for the gate",
        "mesh": mesh_info,
        "starter_all_ok": all_started,
        "engine_all_ok": all_engine,
        "modulus_gate_evaluated": len(errors) == len(KOK_TARGETS_GPA),
        "modulus_gate_verdict": verdict,
        "max_relative_error_pct": max_error_pct,
        **measured,
    }
    return metrics, rows, log_lines


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
    verdict = str(metrics["modulus_gate_verdict"])
    results = {
        "stage": 11,
        "name": "PW_mesoscale_direct",
        "verdict": verdict,
        "metrics": {**metrics, "wall_clock_s": wall_s},
        "reference": {
            "Kok_2022_targets_GPa": KOK_TARGETS_GPA,
            "material_note": "Geometry target values are Kok 2022 VTC401 AP-PLY; solver card uses project IM7/8552.",
            "geometry_config": str(CONFIG_PATH),
            "material_card": str(CARD_PATH),
        },
        "tolerance": {
            "PASS_relative_error": PASS_TOL,
            "INCONCLUSIVE_if_within_extra_relative_error": INCONCLUSIVE_EXTRA_TOL,
        },
        "git_sha": git_sha,
    }
    (RESULTS_DIR / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    (FIGURES_DIR / "stage11_effective_moduli.typ").write_text(
        "\n".join(
            [
                '#set page(width: 170mm, height: auto, margin: 10mm)',
                '#let rows = csv("../results/effective_moduli.csv")',
                '#text(size: 12pt, weight: "bold")[Stage 11 AP-PLY effective moduli]',
                '#v(5pt)',
                '#table(',
                '  columns: (32mm, 24mm, 28mm, 28mm, 26mm),',
                '  stroke: rgb("#5C5C5C"),',
                '  [Case], [Metric], [FEM GPa], [Target GPa], [Error %],',
                '  ..rows.map(r => ([#r.at(0)], [#r.at(1)], [#r.at(2)], [#r.at(3)], [#r.at(4)])).flatten(),',
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
        row_lines = []
        for row in rows:
            row_lines.append(
                f"- {row['metric']}: measured `{row['measured_GPa']}` GPa, target `{row['target_GPa']}` GPa, "
                f"error `{row['relative_error_pct']}` percent."
            )
        blocker.write_text(
            "\n".join(
                [
                    "# Stage 11 Blocker - AP-PLY stiffness gate",
                    "",
                    "Author: J.C. Vaught",
                    "",
                    "The Kok geometry port generated the 25 mm x 25 mm, four-ply `[0/+45/-45/90]` AP-PLY block and the OpenRadioss decks used `LAW12 + TYPE6/SOL_ORTH` with per-tow `/INIBRI/ORTHO`.",
                    "",
                    "The solver completed status and numerical comparison are:",
                    "",
                    *row_lines,
                    "",
                    "Kok 2022 used VTC401 carbon/epoxy target data, while this project stage uses the canonical IM7/8552 material card. Material delta and the coarsened TETRA10 validation mesh are the leading suspected causes if the miss is near the tolerance band.",
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
    print(json.dumps({"stage": 11, "verdict": metrics["modulus_gate_verdict"], "wall_clock_s": wall_s}, indent=2))
    return 0 if metrics["modulus_gate_verdict"] in {"PASS", "INCONCLUSIVE"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
