#!/usr/bin/env python3
"""Empirical LAW12/LAW14 + TYPE6 solid orientation probe.

This is intentionally separate from the stage runners. It generates small
HEXA8 decks for candidate /PROP/TYPE6 orientation conventions, runs
OpenRadioss, converts the animation output to VTK, and compares recovered
E_x against the unidirectional lamina transformation.
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
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[2]
CARD_PATH = ROOT / "references" / "material_cards" / "im7_8552.json"
THIS_DIR = ROOT / "references" / "orientation_probe_decks"
RUNS_DIR = THIS_DIR / "runs"
RESULTS_CSV = THIS_DIR / "orientation_probe_results.csv"
RESULTS_JSON = THIS_DIR / "orientation_probe_results.json"

OR_ROOT = Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss"))
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"
ENGINE = OR_ROOT / "exec" / "engine_linux64_gf"
ANIM_TO_VTK = OR_ROOT / "exec" / "anim_to_vtk_linux64_gf"
N_THREADS = int(os.environ.get("RAD_NT", "4"))

ANGLES = (0.0, 30.0, 45.0, 60.0, 90.0)
LAWS = ("LAW12", "LAW14")
IORTHS = (0, 1)


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
class MeshSpec:
    length: float = 0.040
    width: float = 0.010
    thickness: float = 0.010
    nx: int = 4
    ny: int = 4
    nz: int = 1
    strain: float = 5.0e-4
    run_time: float = 2.0e-4


@dataclass(frozen=True)
class PropState:
    vx: float = 1.0
    vy: float = 0.0
    vz: float = 0.0
    skew_id: int = 0
    ip: int = 1
    phi: float = 0.0
    px: float = 0.0
    py: float = 0.0
    pz: float = 0.0
    inibri: bool = False


@dataclass(frozen=True)
class Mechanism:
    name: str
    family: str
    ip_values: tuple[int, ...]
    description: str
    state: Callable[[float, int], PropState]


@dataclass
class MeshData:
    lines: list[str]
    left_nodes: list[int]
    right_nodes: list[int]
    element_ids: list[int]


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


def analytic_ex(theta_deg: float, mat: Material) -> float:
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


def run_cmd(cmd: list[str], cwd: Path, timeout: int = 180) -> tuple[int, str]:
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            env=radioss_env(),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=timeout,
        )
        return proc.returncode, proc.stdout
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout if isinstance(exc.stdout, str) else ""
        return 124, output + f"\nTIMEOUT after {timeout}s\n"


def law_block(law: str, mat: Material) -> list[str]:
    nu31 = mat.nu13 * mat.e3 / mat.e1
    header = f"/MAT/{law}/1"
    common = [
        header,
        f"IM7_8552_{law}_orientation_probe",
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
    ]
    if law == "LAW12":
        return common + [
            "#          sigma_3yt           sigma_3yc          sigma_13yt          sigma_13yc",
            fmt_f(mat.zt, mat.zc, mat.s13, mat.s13),
            "#              alpha                  Ef                   c          EPS_RATE_0   STRFLAG",
            fmt_f(0.0, 0.0, 0.0, 0.0) + fmt_i(1),
        ]
    if law == "LAW14":
        return common + [
            "#              ALPHA                 E_f                   c          EPS_RATE_0       ICC",
            fmt_f(0.0, 0.0, 0.0, 0.0) + fmt_i(1),
        ]
    raise ValueError(f"unsupported law {law}")


def mesh_blocks(mesh: MeshSpec) -> MeshData:
    nodes: list[str] = ["/NODE"]
    node_id: dict[tuple[int, int, int], int] = {}
    nid = 1
    for k in range(mesh.nz + 1):
        for j in range(mesh.ny + 1):
            for i in range(mesh.nx + 1):
                node_id[(i, j, k)] = nid
                x = mesh.length * i / mesh.nx
                y = mesh.width * (j / mesh.ny - 0.5)
                z = mesh.thickness * (k / mesh.nz - 0.5)
                nodes.append(f"{nid:10d}{x:20.12g}{y:20.12g}{z:20.12g}")
                nid += 1

    bricks = ["/BRICK/1"]
    element_ids: list[int] = []
    eid = 1
    for k in range(mesh.nz):
        for j in range(mesh.ny):
            for i in range(mesh.nx):
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
                element_ids.append(eid)
                eid += 1

    left = [node_id[(0, j, k)] for k in range(mesh.nz + 1) for j in range(mesh.ny + 1)]
    right = [node_id[(mesh.nx, j, k)] for k in range(mesh.nz + 1) for j in range(mesh.ny + 1)]
    lines = nodes + ["/PART/1", "orientation_probe_block", fmt_i(1, 1, 0)] + bricks
    return MeshData(lines=lines, left_nodes=left, right_nodes=right, element_ids=element_ids)


def group_block(group_id: int, name: str, node_ids: list[int]) -> list[str]:
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for i in range(0, len(node_ids), 10):
        lines.append(fmt_i(*node_ids[i : i + 10]))
    return lines


def skew_block(theta_deg: float) -> list[str]:
    theta = math.radians(theta_deg)
    c = math.cos(theta)
    s = math.sin(theta)
    return [
        "/SKEW/FIX/1",
        f"skew_x_theta_{theta_deg:.1f}",
        fmt_f(0.0, 0.0, 0.0),
        fmt_f(c, s, 0.0),
        fmt_f(-s, c, 0.0),
    ]


def type6_property(state: PropState, iorth: int) -> list[str]:
    return [
        "/PROP/TYPE6/1",
        "type6_sol_orth_orientation_probe",
        "#   Isolid    Ismstr               Icpre  Itetra10     Inpts   Itetra4    Iframe                  Dn",
        fmt_i(24, 4) + f"{1:20d}{0:10d}{0:10d}{0:10d}{2:10d}{0.0:20.12g}",
        "#                 qa                  qb                   h",
        fmt_f(0.0, 0.0, 0.0),
        "#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth",
        fmt_f(state.vx, state.vy, state.vz) + fmt_i(state.skew_id, state.ip, iorth),
        "#                Phi                 Px                  Py                  Pz",
        fmt_f(state.phi, state.px, state.py, state.pz),
        "#             dt_min   istrain      IHKT",
        fmt_f(0.0) + fmt_i(0, 0),
    ]


def inibri_ortho_block(theta_deg: float, element_ids: list[int]) -> list[str]:
    theta = math.radians(theta_deg)
    c = math.cos(theta)
    s = math.sin(theta)
    g1 = (c, s, 0.0)
    g2 = (-s, c, 0.0)
    lines = [
        "/INIBRI/ORTHO",
        "#  brick_ID  Nb_layer  Isolnod Prop_type    Isolid",
        "#                  X1                  Y1                  Z1                  X2                  Y2",
        "#                  Z2",
    ]
    for eid in element_ids:
        lines.append(fmt_i(eid, 1, 8, 6, 24))
        lines.append(fmt_f(g1[0], g1[1], g1[2], g2[0], g2[1]))
        lines.append(fmt_f(g2[2]))
    return lines


def control_blocks(mesh: MeshSpec, mesh_data: MeshData) -> list[str]:
    left = mesh_data.left_nodes
    right = mesh_data.right_nodes
    lines = [
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
        fmt_f(mesh.run_time, 1.0),
        "/IMPDISP/1",
        "right_grip_x",
        "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
        f"{1:10d}{'X':>10}{0:10d}{0:10d}{101:10d}{0:10d}{0:10d}",
        "#            Scale_x             Scale_y              Tstart               Tstop",
        fmt_f(1.0, mesh.strain * mesh.length, 0.0, 0.0),
    ]
    lines.extend(group_block(100, "left_grip", left))
    lines.extend(group_block(101, "right_grip", right))
    lines.extend(group_block(102, "left_anchor_yz", [left[0]]))
    lines.extend(group_block(103, "left_anchor_z", [left[-1]]))
    return lines


def write_starter(
    case_dir: Path,
    job: str,
    law: str,
    mechanism: Mechanism,
    theta_deg: float,
    iorth: int,
    mesh: MeshSpec,
    mat: Material,
) -> Path:
    mesh_data = mesh_blocks(mesh)
    state = mechanism.state(theta_deg, mechanism.ip_values[0])
    lines = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        "      2023         0",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"{law} TYPE6 {mechanism.name} theta {theta_deg:.1f} Iorth {iorth}",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 4) + f"{0:20d}{2:40d}",
    ]
    lines.extend(law_block(law, mat))
    lines.extend(mesh_data.lines)
    lines.extend(type6_property(state, iorth))
    if state.skew_id:
        lines.extend(skew_block(theta_deg))
    if state.inibri:
        lines.extend(inibri_ortho_block(theta_deg, mesh_data.element_ids))
    lines.extend(control_blocks(mesh, mesh_data))
    lines.extend(["/END", ""])
    path = case_dir / f"{job}_0000.rad"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_engine(case_dir: Path, job: str, mesh: MeshSpec) -> Path:
    path = case_dir / f"{job}_0001.rad"
    lines = [
        "#RADIOSS ENGINE",
        "/ANIM/DT",
        fmt_f(mesh.run_time, mesh.run_time),
        "/ANIM/VECT/DISP",
        "/ANIM/VECT/FREAC",
        "/ANIM/BRICK/TENS/STRESS/ALL",
        "/ANIM/BRICK/TENS/STRAIN/ALL",
        "/ANIM/GZIP",
        "/TFILE/4",
        fmt_f(mesh.run_time / 20.0),
        "/RFILE",
        fmt_i(1000),
        "/PRINT/-100/55",
        f"/RUN/{job}/1",
        fmt_f(mesh.run_time),
        "/VERS/2023",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def convert_anim_to_vtk(case_dir: Path, job: str) -> Path:
    anim = case_dir / f"{job}A001"
    gz = case_dir / f"{job}A001.gz"
    if gz.exists():
        with gzip.open(gz, "rb") as src, anim.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    if not anim.exists():
        raise FileNotFoundError(f"animation frame not found: {anim.name}")

    proc = subprocess.run(
        [str(ANIM_TO_VTK), str(anim)],
        cwd=str(case_dir),
        env=radioss_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=60,
    )
    if proc.returncode != 0:
        stderr = proc.stderr.decode("utf-8", errors="replace")
        raise RuntimeError(f"anim_to_vtk failed rc={proc.returncode}: {stderr[:800]}")
    vtk = case_dir / f"{job}A001.vtk"
    if proc.stdout.lstrip().startswith(b"# vtk"):
        vtk.write_bytes(proc.stdout[proc.stdout.find(b"# vtk") :])
    if not vtk.exists():
        raise FileNotFoundError("anim_to_vtk produced no VTK stdout")
    return vtk


def _array_by_name(mapping, contains: str):
    for name in mapping.keys():
        if contains in name:
            return mapping[name]
    raise KeyError(f"no array containing {contains!r}; available={list(mapping.keys())}")


def extract_modulus(vtk_path: Path, mesh: MeshSpec) -> dict[str, float | str]:
    import numpy as np
    import pyvista as pv  # type: ignore[import-not-found]

    grid = pv.read(str(vtk_path))
    stress = np.asarray(_array_by_name(grid.cell_data, "Strs"), dtype=float)
    strain = np.asarray(_array_by_name(grid.cell_data, "Stra"), dtype=float)
    centers = grid.cell_centers().points
    x0 = float(centers[:, 0].min())
    x1 = float(centers[:, 0].max())
    gauge = (centers[:, 0] >= x0 + 0.25 * (x1 - x0)) & (centers[:, 0] <= x0 + 0.75 * (x1 - x0))
    sigma_x = float(np.mean(stress[gauge, 0]))
    eps_x_cell = abs(float(np.mean(strain[gauge, 0])))
    e_cell_stress = abs(sigma_x) / eps_x_cell if eps_x_cell > 0.0 else math.nan

    e_global_stress = math.nan
    e_global_reaction = math.nan
    force_x = math.nan
    reaction_sigma_x = math.nan
    eps_x_eng = math.nan
    try:
        disp = np.asarray(_array_by_name(grid.point_data, "Displacement"), dtype=float)
        reaction = np.asarray(_array_by_name(grid.point_data, "Reaction"), dtype=float)
        points = np.asarray(grid.points, dtype=float)
        original_x = points[:, 0] - disp[:, 0]
        left = original_x < original_x.min() + 1.0e-9
        right = original_x > original_x.max() - 1.0e-9
        eps_x_eng = abs(float(np.mean(disp[right, 0]) - np.mean(disp[left, 0]))) / mesh.length
        force_x = float(np.sum(reaction[right, 0]))
        reaction_sigma_x = force_x / (mesh.width * mesh.thickness)
        if eps_x_eng > 0.0:
            e_global_stress = abs(sigma_x) / eps_x_eng
            e_global_reaction = abs(reaction_sigma_x) / eps_x_eng
    except Exception:
        pass

    return {
        "e_global_stress_pa": e_global_stress,
        "e_global_reaction_pa": e_global_reaction,
        "e_cell_stress_pa": e_cell_stress,
        "sigma_x_mean_pa": sigma_x,
        "strain_x_cell": eps_x_cell,
        "strain_x_engineering": eps_x_eng,
        "reaction_force_x_n": force_x,
        "reaction_sigma_x_pa": reaction_sigma_x,
        "vtk_path": str(vtk_path.relative_to(ROOT)),
    }


def mechanisms() -> list[Mechanism]:
    def unit(_: float) -> tuple[float, float, float]:
        return (1.0, 0.0, 0.0)

    def fiber(theta: float) -> tuple[float, float, float]:
        r = math.radians(theta)
        return (math.cos(r), math.sin(r), 0.0)

    def transverse(theta: float) -> tuple[float, float, float]:
        r = math.radians(theta)
        return (-math.sin(r), math.cos(r), 0.0)

    return [
        Mechanism(
            "M1_phi_ip1",
            "M1_phi",
            (1,),
            "PROP Phi on Ip=1 (r,s plane)",
            lambda theta, ip: PropState(vx=1.0, vy=0.0, vz=0.0, ip=ip, phi=theta),
        ),
        Mechanism(
            "M1_phi_ip3",
            "M1_phi",
            (3,),
            "PROP Phi on Ip=3 (t,r plane)",
            lambda theta, ip: PropState(vx=1.0, vy=0.0, vz=0.0, ip=ip, phi=theta),
        ),
        Mechanism(
            "M2_skew_requested_ip1",
            "M2_skew_requested",
            (1,),
            "Prior-style skew_ID set but Ip=1, so CFG/source indicate skew is ignored",
            lambda theta, ip: PropState(vx=1.0, vy=0.0, vz=0.0, skew_id=1, ip=ip, phi=0.0),
        ),
        Mechanism(
            "M2_skew_requested_ip3",
            "M2_skew_requested",
            (3,),
            "skew_ID set but Ip=3, retained for requested cross product",
            lambda theta, ip: PropState(vx=1.0, vy=0.0, vz=0.0, skew_id=1, ip=ip, phi=0.0),
        ),
        Mechanism(
            "M2_skew_ip0",
            "M2_skew",
            (0,),
            "documented skew_ID use: Ip=0 references /SKEW/FIX",
            lambda theta, ip: PropState(vx=1.0, vy=0.0, vz=0.0, skew_id=1, ip=ip, phi=0.0),
        ),
        Mechanism(
            "M3_vector_requested_ip1",
            "M3_vector_requested",
            (1,),
            "fiber vector populated but Ip=1, retained for requested cross product",
            lambda theta, ip: PropState(vx=fiber(theta)[0], vy=fiber(theta)[1], vz=fiber(theta)[2], ip=ip),
        ),
        Mechanism(
            "M3_vector_requested_ip3",
            "M3_vector_requested",
            (3,),
            "fiber vector populated but Ip=3, retained for requested cross product",
            lambda theta, ip: PropState(vx=fiber(theta)[0], vy=fiber(theta)[1], vz=fiber(theta)[2], ip=ip),
        ),
        Mechanism(
            "M3_vector_ip11",
            "M3_vector_project",
            (11,),
            "project V on plane normal to local 3-axis",
            lambda theta, ip: PropState(vx=fiber(theta)[0], vy=fiber(theta)[1], vz=fiber(theta)[2], ip=ip),
        ),
        Mechanism(
            "M3_vector_ip13",
            "M3_vector_project",
            (13,),
            "project V on plane normal to local 2-axis",
            lambda theta, ip: PropState(vx=fiber(theta)[0], vy=fiber(theta)[1], vz=fiber(theta)[2], ip=ip),
        ),
        Mechanism(
            "M4_vector_point_ip24",
            "M4_vector_point",
            (24,),
            "source-driven Vj+Pt mode; V is transverse axis, P shifts the normal",
            lambda theta, ip: PropState(
                vx=transverse(theta)[0],
                vy=transverse(theta)[1],
                vz=transverse(theta)[2],
                ip=ip,
                px=-math.sin(math.radians(theta)),
                py=math.cos(math.radians(theta)),
                pz=1.0,
            ),
        ),
        Mechanism(
            "M5_inibri_ortho_ip1",
            "M5_inibri_ortho",
            (1,),
            "/INIBRI/ORTHO global g1/g2 axes with neutral Ip=1 property",
            lambda theta, ip: PropState(vx=1.0, vy=0.0, vz=0.0, ip=ip, inibri=True),
        ),
        Mechanism(
            "M5_inibri_ortho_ip20",
            "M5_inibri_ortho",
            (20,),
            "/INIBRI/ORTHO global g1/g2 axes with source-documented Ip=20",
            lambda theta, ip: PropState(vx=1.0, vy=0.0, vz=0.0, ip=ip, inibri=True),
        ),
        Mechanism(
            "M6_point_ip21",
            "M6_source_point",
            (21,),
            "source-driven Pt mode: fiber from point P toward element centroid",
            lambda theta, ip: PropState(
                vx=1.0,
                vy=0.0,
                vz=0.0,
                ip=ip,
                px=-math.cos(math.radians(theta)),
                py=-math.sin(math.radians(theta)),
                pz=0.0,
            ),
        ),
        Mechanism(
            "M6_vphi_ip23",
            "M6_source_vector_phi",
            (23,),
            "source-driven Vj+Phi mode",
            lambda theta, ip: PropState(vx=0.0, vy=0.0, vz=1.0, ip=ip, phi=theta),
        ),
    ]


def job_name(law: str, mech: str, iorth: int, angle: float) -> str:
    angle_tag = f"{int(round(angle)):03d}"
    return f"ori_{law.lower()}_{mech}_io{iorth}_a{angle_tag}"[:110]


def tail(text: str, limit: int = 3000) -> str:
    return text[-limit:] if len(text) > limit else text


def error_summary(text: str) -> str:
    return tail(text).replace("\r", "").replace("\n", "\\n")


def run_case(
    law: str,
    mechanism: Mechanism,
    angle: float,
    iorth: int,
    mat: Material,
    mesh: MeshSpec,
    keep_logs: bool,
) -> dict[str, object]:
    job = job_name(law, mechanism.name, iorth, angle)
    case_dir = RUNS_DIR / job
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    write_starter(case_dir, job, law, mechanism, angle, iorth, mesh, mat)
    write_engine(case_dir, job, mesh)

    e_ref = analytic_ex(angle, mat)
    row: dict[str, object] = {
        "law": law,
        "mechanism": mechanism.name,
        "family": mechanism.family,
        "description": mechanism.description,
        "iorth": iorth,
        "ip": mechanism.ip_values[0],
        "angle_deg": angle,
        "analytic_gpa": e_ref / 1.0e9,
        "job": job,
        "starter_check_rc": "",
        "starter_rc": "",
        "engine_rc": "",
        "e_global_stress_gpa": "",
        "e_global_reaction_gpa": "",
        "e_cell_stress_gpa": "",
        "sigma_x_mean_pa": "",
        "strain_x_cell": "",
        "strain_x_engineering": "",
        "reaction_force_x_n": "",
        "reaction_sigma_x_pa": "",
        "rel_err_pct": "",
        "angle_pass": False,
        "combo_pass_count": 0,
        "verdict": "BLOCKED",
        "error": "",
        "vtk_path": "",
    }

    check_rc, check_out = run_cmd(
        [str(STARTER), "-i", f"{job}_0000.rad", "-nt", str(N_THREADS), "-check"], case_dir
    )
    row["starter_check_rc"] = check_rc
    (case_dir / "starter_check.log").write_text(check_out, encoding="utf-8", errors="replace")
    if check_rc != 0:
        row["error"] = error_summary(check_out)
        return row

    starter_rc, starter_out = run_cmd([str(STARTER), "-i", f"{job}_0000.rad", "-nt", str(N_THREADS)], case_dir)
    row["starter_rc"] = starter_rc
    (case_dir / "starter.log").write_text(starter_out, encoding="utf-8", errors="replace")
    if starter_rc != 0:
        row["error"] = error_summary(starter_out)
        return row

    engine_rc, engine_out = run_cmd([str(ENGINE), "-i", f"{job}_0001.rad", "-nt", str(N_THREADS)], case_dir)
    row["engine_rc"] = engine_rc
    (case_dir / "engine.log").write_text(engine_out, encoding="utf-8", errors="replace")
    if engine_rc != 0:
        row["error"] = error_summary(engine_out)
        return row

    try:
        vtk_path = convert_anim_to_vtk(case_dir, job)
        values = extract_modulus(vtk_path, mesh)
        e_global_stress = float(values["e_global_stress_pa"])
        rel_err = abs(e_global_stress - e_ref) / e_ref
        row.update(
            {
                "e_global_stress_gpa": e_global_stress / 1.0e9,
                "e_global_reaction_gpa": float(values["e_global_reaction_pa"]) / 1.0e9,
                "e_cell_stress_gpa": float(values["e_cell_stress_pa"]) / 1.0e9,
                "sigma_x_mean_pa": values["sigma_x_mean_pa"],
                "strain_x_cell": values["strain_x_cell"],
                "strain_x_engineering": values["strain_x_engineering"],
                "reaction_force_x_n": values["reaction_force_x_n"],
                "reaction_sigma_x_pa": values["reaction_sigma_x_pa"],
                "rel_err_pct": 100.0 * rel_err,
                "angle_pass": rel_err <= 0.02,
                "verdict": "FAIL",
                "vtk_path": values["vtk_path"],
            }
        )
    except Exception as exc:
        row["error"] = repr(exc)
        row["verdict"] = "BLOCKED"

    if not keep_logs:
        for pattern in ("*.out", "*.sta", "*.T??", "*.rst", "*.ctl", "*.log"):
            for path in case_dir.glob(pattern):
                path.unlink(missing_ok=True)
    return row


def finalize_verdicts(rows: list[dict[str, object]]) -> None:
    keys = {(r["law"], r["mechanism"], r["iorth"], r["ip"]) for r in rows}
    pass_counts: dict[tuple[object, ...], int] = {}
    for key in keys:
        pass_counts[key] = sum(
            1
            for r in rows
            if (r["law"], r["mechanism"], r["iorth"], r["ip"]) == key and bool(r["angle_pass"])
        )
    for row in rows:
        key = (row["law"], row["mechanism"], row["iorth"], row["ip"])
        row["combo_pass_count"] = pass_counts[key]
        if row["verdict"] == "BLOCKED":
            continue
        row["verdict"] = "PASS" if pass_counts[key] >= 3 else "FAIL"


def write_results(rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    finalize_verdicts(rows)
    fieldnames = list(rows[0].keys())
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    with RESULTS_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    RESULTS_JSON.write_text(json.dumps(rows, indent=2), encoding="utf-8")


def parse_csv_set(value: str | None, default: tuple[str, ...]) -> tuple[str, ...]:
    if not value:
        return default
    return tuple(part.strip() for part in value.split(",") if part.strip())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="store_true", help="run a short source-guided pilot")
    parser.add_argument("--laws", help="comma-separated laws, default LAW12,LAW14")
    parser.add_argument("--mechanisms", help="comma-separated mechanism names")
    parser.add_argument("--angles", help="comma-separated angles")
    parser.add_argument("--iorths", help="comma-separated Iorth values")
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--keep-logs", action="store_true")
    args = parser.parse_args()

    if not STARTER.exists() or not ENGINE.exists() or not ANIM_TO_VTK.exists():
        raise SystemExit(f"missing OpenRadioss executable under {OR_ROOT}")

    mat = load_material()
    mesh = MeshSpec()
    all_mechs = mechanisms()
    by_name = {mech.name: mech for mech in all_mechs}

    laws = parse_csv_set(args.laws, LAWS)
    angles = tuple(float(x) for x in parse_csv_set(args.angles, tuple(str(a) for a in ANGLES)))
    iorths = tuple(int(x) for x in parse_csv_set(args.iorths, tuple(str(i) for i in IORTHS)))
    if args.mechanisms:
        selected = [by_name[name] for name in parse_csv_set(args.mechanisms, ())]
    elif args.pilot:
        selected = [
            by_name["M1_phi_ip1"],
            by_name["M2_skew_ip0"],
            by_name["M3_vector_ip11"],
            by_name["M5_inibri_ortho_ip1"],
            by_name["M6_vphi_ip23"],
        ]
        laws = ("LAW12",)
        angles = (0.0, 45.0, 90.0)
        iorths = (0, 1)
    else:
        selected = all_mechs

    cases: list[tuple[str, Mechanism, float, int]] = [
        (law, mechanism, angle, iorth)
        for law in laws
        for mechanism in selected
        for iorth in iorths
        for angle in angles
    ]
    if args.max_cases:
        cases = cases[: args.max_cases]

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    start = time.time()
    total = len(cases)
    for idx, (law, mechanism, angle, iorth) in enumerate(cases, start=1):
        print(
            f"[{idx:04d}/{total:04d}] {law} {mechanism.name} Iorth={iorth} angle={angle:g}",
            flush=True,
        )
        row = run_case(law, mechanism, angle, iorth, mat, mesh, args.keep_logs)
        rows.append(row)
        write_results(rows)
        if row["verdict"] == "BLOCKED":
            print(f"  BLOCKED check/start/engine/extract: {str(row['error'])[:160]}", flush=True)
        else:
            print(
                "  E_global={:.4g} GPa analytic={:.4g} GPa err={:.2f}%".format(
                    float(row["e_global_stress_gpa"]),
                    float(row["analytic_gpa"]),
                    float(row["rel_err_pct"]),
                ),
                flush=True,
            )

    write_results(rows)
    elapsed = time.time() - start
    print(f"wrote {RESULTS_CSV.relative_to(ROOT)} ({len(rows)} rows) in {elapsed:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
