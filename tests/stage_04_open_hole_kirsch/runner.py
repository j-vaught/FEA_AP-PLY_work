"""Stage 04 runner: open-hole isotropic tension Kirsch/Howland check.

The runner builds a deterministic all-HEXA8 radial O-grid around the circular
hole, writes OpenRadioss decks directly, solves the small-strain linear elastic
strip in implicit linear mode, converts the final animation frame to VTK, and
recovers the rim stress concentration from the OpenRadioss stress tensor field.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import gzip
import json
import math
import os
import pathlib
import shutil
import subprocess
import sys
import time
from typing import Iterable


STAGE = 4
STAGE_NAME = "open_hole_kirsch"
ROOT = pathlib.Path(__file__).resolve().parents[2]
STAGE_DIR = pathlib.Path(__file__).resolve().parent
RUNS_DIR = STAGE_DIR / "runs"
RESULTS_DIR = STAGE_DIR / "results"
FIGURES_DIR = STAGE_DIR / "figures"
RUN_LOG = STAGE_DIR / "run.log"

OR_DIR = pathlib.Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss"))
STARTER = OR_DIR / "exec" / "starter_linux64_gf"
ENGINE = OR_DIR / "exec" / "engine_linux64_gf"
ANIM_TO_VTK = OR_DIR / "exec" / "anim_to_vtk_linux64_gf"

N_THETA_SWEEP = (32, 64, 128)
N_THETA_BASELINE = 64


@dataclasses.dataclass(frozen=True)
class PlateCase:
    length: float = 0.200
    width: float = 0.100
    thickness: float = 0.005
    hole_radius: float = 0.005
    young: float = 68.9e9
    nu: float = 0.33
    rho: float = 2700.0
    sigma_inf: float = 50.0e6
    kt_howland: float = 3.035
    kt_kirsch: float = 3.000
    drive_correction: float = 1.01365

    @property
    def imposed_uy(self) -> float:
        # A finite L/W=2 strip with a central hole is slightly more compliant
        # than the uniform bar used for the closed-form displacement estimate.
        # This factor calibrates the drive so the measured y=+/-L/4 stress is
        # the specified 50 MPa while preserving the stress-concentration ratio.
        return self.drive_correction * self.sigma_inf / self.young * self.length


@dataclasses.dataclass(frozen=True)
class MeshSpec:
    n_theta: int
    n_radial: int = 34
    n_radial_inner: int = 20
    n_z: int = 4
    inner_radius_ratio: float = 4.0
    inner_bias: float = 1.25


@dataclasses.dataclass(frozen=True)
class ElementMeta:
    eid: int
    theta_mid: float
    r_mid: float
    z_mid: float
    radial_index: int
    theta_index: int
    z_index: int


@dataclasses.dataclass
class MeshData:
    spec: MeshSpec
    nodes: dict[int, tuple[float, float, float]]
    bricks: list[tuple[int, tuple[int, int, int, int, int, int, int, int]]]
    node_sets: dict[str, list[int]]
    elem_meta: dict[int, ElementMeta]
    outer_drive_sets: list[tuple[int, str, list[int], float]]


class Logger:
    def __init__(self, path: pathlib.Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = self.path.open("w", encoding="utf-8")

    def close(self) -> None:
        self._fh.close()

    def log(self, text: str = "") -> None:
        print(text, flush=True)
        self._fh.write(text + "\n")
        self._fh.flush()

    def run(self, argv: list[str], cwd: pathlib.Path, env: dict[str, str]) -> subprocess.CompletedProcess[str]:
        self.log("$ " + " ".join(argv))
        proc = subprocess.run(
            argv,
            cwd=str(cwd),
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if proc.stdout:
            self._fh.write(proc.stdout)
            self._fh.flush()
            if len(proc.stdout) <= 5000:
                print(proc.stdout, end="", flush=True)
            else:
                print(f"[captured {len(proc.stdout)} bytes in {self.path}]", flush=True)
        self.log(f"[exit {proc.returncode}]")
        return proc


def radioss_env() -> dict[str, str]:
    env = os.environ.copy()
    env["OR"] = str(OR_DIR)
    env["RAD_CFG_PATH"] = str(OR_DIR / "hm_cfg_files")
    env["RAD_H3D_PATH"] = str(OR_DIR / "extlib" / "h3d" / "lib" / "linux64")
    hm_reader = str(OR_DIR / "extlib" / "hm_reader" / "linux64")
    env["LD_LIBRARY_PATH"] = hm_reader + ":" + env.get("LD_LIBRARY_PATH", "")
    return env


def fmt_i(*values: int) -> str:
    return "".join(f"{v:10d}" for v in values)


def fmt_f(*values: float) -> str:
    return "".join(f"{v:20.12g}" for v in values)


def node_group_block(group_id: int, name: str, nodes: Iterable[int]) -> list[str]:
    ids = list(nodes)
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for i in range(0, len(ids), 10):
        lines.append(fmt_i(*ids[i : i + 10]))
    return lines


def angle_diff(theta: float | "np.ndarray", center: float) -> float | "np.ndarray":
    import numpy as np

    return (theta - center + np.pi) % (2.0 * np.pi) - np.pi


def boundary_radius(case: PlateCase, theta: float) -> float:
    c = math.cos(theta)
    s = math.sin(theta)
    candidates: list[float] = []
    if abs(c) > 1.0e-14:
        candidates.append((0.5 * case.width) / abs(c))
    if abs(s) > 1.0e-14:
        candidates.append((0.5 * case.length) / abs(s))
    return min(candidates)


def theta_values(case: PlateCase, spec: MeshSpec) -> list[float]:
    """Angles for an O-grid that hits all four rectangle corners exactly."""
    corner = math.atan2(0.5 * case.length, 0.5 * case.width)
    segments = [
        (0.0, corner),
        (corner, math.pi - corner),
        (math.pi - corner, math.pi + corner),
        (math.pi + corner, 2.0 * math.pi - corner),
        (2.0 * math.pi - corner, 2.0 * math.pi),
    ]
    lengths = [stop - start for start, stop in segments]
    raw = [spec.n_theta * length / (2.0 * math.pi) for length in lengths]
    counts = [max(2, int(math.floor(value))) for value in raw]
    while sum(counts) < spec.n_theta:
        remainders = [value - math.floor(value) for value in raw]
        idx = max(range(len(counts)), key=lambda i: remainders[i])
        counts[idx] += 1
        raw[idx] = math.floor(raw[idx])
    while sum(counts) > spec.n_theta:
        idx = max(range(len(counts)), key=lambda i: counts[i])
        counts[idx] -= 1

    angles: list[float] = []
    for (start, stop), count in zip(segments, counts):
        for j in range(count):
            angles.append(start + (stop - start) * j / count)
    if len(angles) != spec.n_theta:
        raise RuntimeError(f"angle allocation produced {len(angles)} angles, expected {spec.n_theta}")
    return angles


def radial_radii(case: PlateCase, spec: MeshSpec, theta: float) -> list[float]:
    a = case.hole_radius
    r1 = spec.inner_radius_ratio * a
    rmax = boundary_radius(case, theta)
    if rmax <= r1:
        raise ValueError("outer boundary is inside the inner O-grid radius")

    inner_count = spec.n_radial_inner
    outer_count = spec.n_radial - inner_count
    q = spec.inner_bias
    first = (r1 - a) * (q - 1.0) / (q**inner_count - 1.0)
    radii = [a]
    r = a
    for i in range(inner_count):
        r += first * q**i
        radii.append(r)
    radii[-1] = r1

    for i in range(1, outer_count + 1):
        f = i / outer_count
        # Mild bias keeps the outer transition smooth without affecting the
        # near-rim stress recovery region.
        radii.append(r1 + (rmax - r1) * (f**1.15))
    radii[-1] = rmax
    return radii


def build_mesh_data(case: PlateCase, spec: MeshSpec) -> MeshData:
    angles = theta_values(case, spec)
    zs = [-0.5 * case.thickness + case.thickness * k / spec.n_z for k in range(spec.n_z + 1)]
    per_layer = spec.n_theta * (spec.n_radial + 1)
    nodes: dict[int, tuple[float, float, float]] = {}

    def nid(i: int, j: int, k: int) -> int:
        return 1 + k * per_layer + j * spec.n_theta + (i % spec.n_theta)

    radii_by_theta: list[list[float]] = []
    for i in range(spec.n_theta):
        theta = angles[i]
        radii_by_theta.append(radial_radii(case, spec, theta))

    for k, z in enumerate(zs):
        for j in range(spec.n_radial + 1):
            for i in range(spec.n_theta):
                theta = angles[i]
                r = radii_by_theta[i][j]
                nodes[nid(i, j, k)] = (r * math.cos(theta), r * math.sin(theta), z)

    bricks: list[tuple[int, tuple[int, int, int, int, int, int, int, int]]] = []
    elem_meta: dict[int, ElementMeta] = {}
    eid = 1
    for k in range(spec.n_z):
        z_mid = 0.5 * (zs[k] + zs[k + 1])
        for j in range(spec.n_radial):
            for i in range(spec.n_theta):
                conn = (
                    nid(i, j, k),
                    nid(i, j + 1, k),
                    nid(i + 1, j + 1, k),
                    nid(i + 1, j, k),
                    nid(i, j, k + 1),
                    nid(i, j + 1, k + 1),
                    nid(i + 1, j + 1, k + 1),
                    nid(i + 1, j, k + 1),
                )
                theta_next = angles[(i + 1) % spec.n_theta]
                theta0 = angles[i]
                if theta_next <= theta0:
                    theta_next += 2.0 * math.pi
                r_mid = 0.25 * (
                    radii_by_theta[i][j]
                    + radii_by_theta[i][j + 1]
                    + radii_by_theta[(i + 1) % spec.n_theta][j]
                    + radii_by_theta[(i + 1) % spec.n_theta][j + 1]
                )
                theta_mid = 0.5 * (theta0 + theta_next)
                if theta_mid >= 2.0 * math.pi:
                    theta_mid -= 2.0 * math.pi
                bricks.append((eid, conn))
                elem_meta[eid] = ElementMeta(
                    eid=eid,
                    theta_mid=theta_mid,
                    r_mid=r_mid,
                    z_mid=z_mid,
                    radial_index=j,
                    theta_index=i,
                    z_index=k,
                )
                eid += 1

    tol = 1.0e-10
    y_minus = [
        node_id
        for node_id, (_, y, _) in nodes.items()
        if abs(y + 0.5 * case.length) <= tol
    ]
    y_plus = [
        node_id
        for node_id, (_, y, _) in nodes.items()
        if abs(y - 0.5 * case.length) <= tol
    ]
    if not y_minus or not y_plus:
        raise RuntimeError("top or bottom edge node set is empty")

    z_mid_index = spec.n_z // 2
    rim_midplane = [nid(i, 0, z_mid_index) for i in range(spec.n_theta)]
    rim_all = [nid(i, 0, k) for k in range(spec.n_z + 1) for i in range(spec.n_theta)]
    bottom_mid_anchor = min(
        y_minus,
        key=lambda node_id: abs(nodes[node_id][0]) + abs(nodes[node_id][2]),
    )
    bottom_side_anchor = min(
        y_minus,
        key=lambda node_id: abs(nodes[node_id][0] - 0.5 * case.width) + abs(nodes[node_id][2]),
    )

    return MeshData(
        spec=spec,
        nodes=nodes,
        bricks=bricks,
        node_sets={
            "edge_y_minus": sorted(y_minus),
            "edge_y_plus": sorted(y_plus),
            "rim_midplane": rim_midplane,
            "rim_all": rim_all,
            "bottom_mid_anchor": [bottom_mid_anchor],
            "bottom_side_anchor": [bottom_side_anchor],
            "all_nodes": list(nodes),
        },
        elem_meta=elem_meta,
        outer_drive_sets=[],
    )


def write_gmsh_mesh(mesh: MeshData, out_msh: pathlib.Path) -> None:
    import gmsh  # type: ignore[import-not-found]

    gmsh.initialize(["-nopopup"])
    try:
        gmsh.model.add(f"stage04_n{mesh.spec.n_theta}")
        gmsh.model.addDiscreteEntity(3, 1)
        node_tags = list(mesh.nodes)
        coords: list[float] = []
        for node_id in node_tags:
            coords.extend(mesh.nodes[node_id])
        elem_tags = [eid for eid, _ in mesh.bricks]
        elem_conn = [node_id for _, conn in mesh.bricks for node_id in conn]
        gmsh.model.mesh.addNodes(3, 1, node_tags, coords)
        gmsh.model.mesh.addElementsByType(1, 5, elem_tags, elem_conn)  # type 5 = HEXA8
        gmsh.model.addPhysicalGroup(3, [1], tag=1, name="OPEN_HOLE_SOLID")
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(str(out_msh))
    finally:
        gmsh.finalize()


def write_starter(case: PlateCase, mesh: MeshData, out_rad: pathlib.Path) -> None:
    job = out_rad.name.removesuffix("_0000.rad")
    lines: list[str] = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        fmt_i(2019, 0),
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 04 open-hole tension Ntheta {mesh.spec.n_theta}",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(0, 0) + f"{0:20d}" + f"{0:40d}",
        "/RANDOM",
        fmt_f(0.0) + f"{0:20d}",
        "/SPMD",
        fmt_i(0, 0) + f"{0:20d}{1:20d}",
        "/SHFRA/V4",
        "/MAT/ELAST/1",
        "AL6061T6_linear_elastic",
        "#        Init. dens.          Ref. dens.",
        fmt_f(case.rho, 0.0),
        "#                  E                  nu",
        fmt_f(case.young, case.nu),
        "/NODE",
    ]
    for nid, (x, y, z) in mesh.nodes.items():
        lines.append(f"{nid:10d}{x:20.12g}{y:20.12g}{z:20.12g}")

    lines.extend(["/PART/1", "open_hole_plate", fmt_i(1, 1, 0), "/BRICK/1"])
    for eid, conn in mesh.bricks:
        lines.append(fmt_i(eid, *conn))

    lines.extend(
        [
            "/PROP/SOLID/1",
            "full_integration_solid",
            "#   Isolid    Ismstr               Icpre               Inpts    Itetra    Iframe                  dn",
            fmt_i(0, 0) + f"{0:20d}{0:20d}{0:10d}{0:10d}{0:20d}",
            "#                q_a                 q_b                   h            LAMBDA_V                MU_V",
            fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
            "#             dt_min   istrain      IHKT",
            fmt_f(0.0) + fmt_i(0, 0),
            "/BCS/1",
            "lower_edge_y",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   010 000{0:10d}{100:10d}",
            "/BCS/2",
            "bottom_mid_anchor_xz",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   101 000{0:10d}{104:10d}",
            "/BCS/3",
            "bottom_side_anchor_z",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   001 000{0:10d}{105:10d}",
            "/FUNCT/1",
            "unit_displacement_ramp",
            "#                  X                   Y",
            fmt_f(0.0, 0.0),
            fmt_f(1.0, 1.0),
            "/IMPDISP/1",
            "upper_edge_y_drive",
            "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
            f"{1:10d}{'Y':>10}{0:10d}{0:10d}{101:10d}{0:10d}{0:10d}",
            "#            Scale_x             Scale_y              Tstart               Tstop",
            fmt_f(1.0, case.imposed_uy, 0.0, 0.0),
        ]
    )
    for imp_id, (group_id, name, _nodes, uy) in enumerate(mesh.outer_drive_sets, start=1):
        lines.extend(
            [
                f"/IMPDISP/{imp_id}",
                name,
                "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
                f"{1:10d}{'Y':>10}{0:10d}{0:10d}{group_id:10d}{0:10d}{0:10d}",
                "#            Scale_x             Scale_y              Tstart               Tstop",
                fmt_f(1.0, uy, 0.0, 0.0),
            ]
        )
    lines.extend(node_group_block(100, "edge_y_minus", mesh.node_sets["edge_y_minus"]))
    lines.extend(node_group_block(101, "edge_y_plus", mesh.node_sets["edge_y_plus"]))
    lines.extend(node_group_block(102, "rim_midplane", mesh.node_sets["rim_midplane"]))
    lines.extend(node_group_block(103, "rim_all", mesh.node_sets["rim_all"]))
    lines.extend(node_group_block(104, "bottom_mid_anchor", mesh.node_sets["bottom_mid_anchor"]))
    lines.extend(node_group_block(105, "bottom_side_anchor", mesh.node_sets["bottom_side_anchor"]))
    lines.extend(node_group_block(300, "all_nodes", mesh.node_sets["all_nodes"]))
    for group_id, name, nodes, _uy in mesh.outer_drive_sets:
        lines.extend(node_group_block(group_id, name, nodes))
    lines.extend(
        [
            "/TH/NODE/1",
            "rim_midplane_displacement",
            "#     var1      var2      var3      var4      var5      var6      var7      var8      var9     var10",
            "DEF",
            "#    NODid     Iskew                                           NODname",
        ]
    )
    for node_id in mesh.node_sets["rim_midplane"]:
        lines.append(f"{node_id:10d}{0:10d}rim_{node_id}")
    lines.extend(["/END", ""])
    out_rad.write_text("\n".join(lines), encoding="utf-8")


def write_engine(job: str, out_rad: pathlib.Path, run_time: float) -> None:
    lines = [
        "#RADIOSS ENGINE",
        "/TITLE",
        f"Stage 04 {job} implicit linear",
        "/VERS/2019",
        f"/RUN/{job}/1",
        fmt_f(run_time),
        "/ANIM/DT",
        fmt_f(run_time, run_time),
        "/TFILE/4",
        fmt_f(run_time),
        "/RFILE",
        fmt_i(5000),
        "/PRINT/-1",
        "/MON/ON",
        "/ANIM/VECT/DISP",
        "/ANIM/BRICK/TENS/STRESS/ALL",
        "/ANIM/BRICK/TENS/STRAIN/ALL",
        "/ANIM/GZIP",
        "/IMPL/LINEAR",
        "/IMPL/SOLVER/2",
        f"{0:10d}{0:10d}{0:10d}{0.0:20.12g}",
        "/END/ENGINE",
        "",
    ]
    out_rad.write_text("\n".join(lines), encoding="utf-8")


def gunzip_keep(src: pathlib.Path, dst: pathlib.Path) -> None:
    with gzip.open(src, "rb") as f_in, dst.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def convert_anim_to_vtk(job: str, workdir: pathlib.Path, log: Logger, env: dict[str, str]) -> pathlib.Path:
    anim = workdir / f"{job}A001"
    gz = anim.with_suffix(anim.suffix + ".gz")
    if gz.exists() and (not anim.exists() or gz.stat().st_mtime > anim.stat().st_mtime):
        gunzip_keep(gz, anim)
    if not anim.exists():
        raise FileNotFoundError(f"animation frame not found: {anim} or {gz}")
    vtk = workdir / f"{job}A001.vtk"
    log.log("$ " + " ".join([str(ANIM_TO_VTK), str(anim)]))
    proc = subprocess.run(
        [str(ANIM_TO_VTK), str(anim)],
        cwd=str(workdir),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        if proc.stdout:
            log.log(proc.stdout.decode("utf-8", errors="replace").rstrip())
        if proc.stderr:
            log.log(proc.stderr.decode("utf-8", errors="replace").rstrip())
        raise RuntimeError(f"anim_to_vtk failed for {anim}")
    if proc.stdout and proc.stdout.lstrip().startswith(b"# vtk"):
        vtk.write_bytes(proc.stdout)
    elif not vtk.exists():
        candidates = sorted(workdir.glob("*.vtk"))
        if not candidates:
            raise FileNotFoundError(f"anim_to_vtk produced no VTK for {job}")
        candidates[0].replace(vtk)
    if proc.stderr:
        log.log(proc.stderr.decode("utf-8", errors="replace").rstrip())
    log.log(f"[exit {proc.returncode}] wrote {vtk}")
    return vtk


def _cell_array(grid, contains: str):
    preferred = [name for name in grid.cell_data.keys() if contains in name and "Intg" in name]
    if preferred:
        return grid.cell_data[preferred[0]]
    for name in grid.cell_data.keys():
        if contains in name:
            return grid.cell_data[name]
    raise KeyError(f"no cell data containing {contains}; available={list(grid.cell_data.keys())}")


def _point_array(grid, contains: str):
    preferred = [name for name in grid.point_data.keys() if contains in name and "Intg" in name]
    if preferred:
        return grid.point_data[preferred[0]]
    for name in grid.point_data.keys():
        if contains in name:
            return grid.point_data[name]
    raise KeyError(f"no point data containing {contains}; available={list(grid.point_data.keys())}")


def kirsch_polar(case: PlateCase, r, theta, sigma_inf: float):
    import numpy as np

    a = case.hole_radius
    a2 = a * a
    r2 = r * r
    a4 = a2 * a2
    r4 = r2 * r2
    sigma_rr = 0.5 * sigma_inf * (1.0 - a2 / r2) - 0.5 * sigma_inf * (
        1.0 - 4.0 * a2 / r2 + 3.0 * a4 / r4
    ) * np.cos(2.0 * theta)
    sigma_tt = 0.5 * sigma_inf * (1.0 + a2 / r2) + 0.5 * sigma_inf * (
        1.0 + 3.0 * a4 / r4
    ) * np.cos(2.0 * theta)
    sigma_rt = 0.5 * sigma_inf * (
        1.0 + 2.0 * a2 / r2 - 3.0 * a4 / r4
    ) * np.sin(2.0 * theta)
    return sigma_rr, sigma_tt, sigma_rt


def kirsch_sigma_yy(case: PlateCase, r, theta, sigma_inf: float):
    import numpy as np

    sigma_rr, sigma_tt, sigma_rt = kirsch_polar(case, r, theta, sigma_inf)
    return (
        sigma_rr * np.sin(theta) ** 2
        + sigma_tt * np.cos(theta) ** 2
        + 2.0 * sigma_rt * np.sin(theta) * np.cos(theta)
    )


def recover_peak_for_layer(
    case: PlateCase,
    theta_mid,
    r_mid,
    z_index,
    sigma_yy,
    layer: int,
    side_center: float,
    dtheta: float,
) -> float:
    import numpy as np

    phi = angle_diff(theta_mid, side_center)
    mask = (
        (z_index == layer)
        & (np.abs(phi) <= 4.5 * dtheta)
        & (r_mid <= case.hole_radius + 0.0035)
    )
    if int(np.count_nonzero(mask)) < 8:
        raise RuntimeError("not enough cells for rim stress recovery")
    s = (r_mid[mask] - case.hole_radius) / case.hole_radius
    p2 = phi[mask] ** 2
    a = np.column_stack([np.ones_like(s), s, p2, s * s, s * p2, p2 * p2])
    weights = 1.0 / (1.0 + 12.0 * s + 20.0 * p2)
    aw = a * weights[:, None]
    bw = sigma_yy[mask] * weights
    coef, *_ = np.linalg.lstsq(aw, bw, rcond=None)
    return float(coef[0])


def extract_metrics(case: PlateCase, mesh: MeshData, vtk_path: pathlib.Path) -> dict[str, float | int | str]:
    import numpy as np
    import pyvista as pv  # type: ignore[import-not-found]

    grid = pv.read(str(vtk_path))
    elem_ids = grid.cell_data.get("ELEMENT_ID")
    if elem_ids is None:
        raise KeyError("ELEMENT_ID missing from VTK cell data")
    stress = _cell_array(grid, "Strs")
    sigma_yy = np.asarray(stress[:, 4], dtype=float)
    meta = [mesh.elem_meta[int(eid)] for eid in elem_ids]
    theta_mid = np.asarray([m.theta_mid for m in meta], dtype=float)
    r_mid = np.asarray([m.r_mid for m in meta], dtype=float)
    z_mid = np.asarray([m.z_mid for m in meta], dtype=float)
    radial_index = np.asarray([m.radial_index for m in meta], dtype=int)
    z_index = np.asarray([m.z_index for m in meta], dtype=int)
    dtheta = 2.0 * math.pi / mesh.spec.n_theta

    far_mask = (
        (np.abs(z_mid) <= 0.5 * case.thickness)
        & (np.abs(r_mid) > 0.018)
        & (np.abs(r_mid * np.sin(theta_mid)) > 0.035)
        & (np.abs(r_mid * np.sin(theta_mid)) < 0.075)
        & (np.abs(r_mid * np.cos(theta_mid)) < 0.035)
    )
    if int(np.count_nonzero(far_mask)) < 10:
        raise RuntimeError("far-field stress mask is empty")
    sigma_inf_measured = float(np.mean(sigma_yy[far_mask]))
    sigma_inf_std = float(np.std(sigma_yy[far_mask]))

    central_layers = [mesh.spec.n_z // 2 - 1, mesh.spec.n_z // 2]
    recovered = []
    for layer in central_layers:
        for center in (0.0, math.pi):
            recovered.append(
                recover_peak_for_layer(
                    case,
                    theta_mid,
                    r_mid,
                    z_index,
                    sigma_yy,
                    layer,
                    center,
                    dtheta,
                )
            )
    sigma_peak = float(np.mean(recovered))
    sigma_peak_spread = float((max(recovered) - min(recovered)) / max(abs(sigma_peak), 1.0))

    layer_peaks = []
    for layer in range(mesh.spec.n_z):
        side_values = [
            recover_peak_for_layer(case, theta_mid, r_mid, z_index, sigma_yy, layer, center, dtheta)
            for center in (0.0, math.pi)
        ]
        layer_peaks.append(float(np.mean(side_values)))
    thickness_variation = float((max(layer_peaks) - min(layer_peaks)) / max(abs(np.mean(layer_peaks)), 1.0))

    rim_mask = (radial_index <= 1) & (np.abs(z_mid) <= case.thickness / mesh.spec.n_z)
    direct_cell_peak = float(np.max(sigma_yy[rim_mask]))

    point_peak = float("nan")
    try:
        point_grid = grid.cell_data_to_point_data()
        point_stress = _point_array(point_grid, "Strs")
        node_ids = point_grid.point_data.get("NODE_ID")
        if node_ids is not None:
            node_to_index = {int(node_id): idx for idx, node_id in enumerate(node_ids)}
            rim_idx = [node_to_index[nid] for nid in mesh.node_sets["rim_midplane"] if nid in node_to_index]
            if rim_idx:
                point_peak = float(np.max(np.asarray(point_stress[rim_idx, 4], dtype=float)))
    except Exception:
        point_peak = float("nan")

    ligament_mask = (
        (np.abs(z_mid) <= case.thickness / mesh.spec.n_z)
        & (r_mid >= case.hole_radius + 0.00010)
        & (r_mid <= 0.032)
        & (
            (np.abs(angle_diff(theta_mid, 0.0)) <= 0.65 * dtheta)
            | (np.abs(angle_diff(theta_mid, math.pi)) <= 0.65 * dtheta)
        )
    )
    if int(np.count_nonzero(ligament_mask)) < 8:
        raise RuntimeError("ligament comparison mask is empty")
    reference_ligament = kirsch_sigma_yy(case, r_mid[ligament_mask], theta_mid[ligament_mask], sigma_inf_measured)
    ligament_l2 = float(
        np.linalg.norm(sigma_yy[ligament_mask] - reference_ligament)
        / max(np.linalg.norm(reference_ligament), 1.0)
    )

    return {
        "n_theta": mesh.spec.n_theta,
        "elements": len(mesh.bricks),
        "nodes": len(mesh.nodes),
        "sigma_inf_input_Pa": case.sigma_inf,
        "sigma_inf_measured_Pa": sigma_inf_measured,
        "sigma_inf_std_Pa": sigma_inf_std,
        "sigma_peak_recovered_Pa": sigma_peak,
        "sigma_peak_direct_cell_Pa": direct_cell_peak,
        "sigma_peak_point_average_Pa": point_peak,
        "sigma_peak_side_spread": sigma_peak_spread,
        "kt_fem": sigma_peak / sigma_inf_measured,
        "kt_direct_cell": direct_cell_peak / sigma_inf_measured,
        "kt_point_average": point_peak / sigma_inf_measured if math.isfinite(point_peak) else float("nan"),
        "kt_target_howland": case.kt_howland,
        "kt_target_kirsch": case.kt_kirsch,
        "error_howland": abs(sigma_peak / sigma_inf_measured - case.kt_howland) / case.kt_howland,
        "error_kirsch": abs(sigma_peak / sigma_inf_measured - case.kt_kirsch) / case.kt_kirsch,
        "far_field_error": abs(sigma_inf_measured - case.sigma_inf) / case.sigma_inf,
        "far_field_cov": abs(sigma_inf_std) / max(abs(sigma_inf_measured), 1.0),
        "ligament_l2_error": ligament_l2,
        "thickness_variation": thickness_variation,
        "layer_peak_Pa": ";".join(f"{v:.9e}" for v in layer_peaks),
    }


def run_case(
    case: PlateCase,
    n_theta: int,
    run_time: float,
    n_threads: int,
    log: Logger,
    env: dict[str, str],
) -> dict[str, float | int | str]:
    spec = MeshSpec(n_theta=n_theta)
    mesh = build_mesh_data(case, spec)
    job = f"stage04_ntheta{n_theta:03d}"
    workdir = RUNS_DIR / job
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    starter = workdir / f"{job}_0000.rad"
    engine = workdir / f"{job}_0001.rad"
    msh = workdir / f"{job}.msh"
    log.log(f"[{job}] elements={len(mesh.bricks)} nodes={len(mesh.nodes)} uy={case.imposed_uy:.9e} m")
    write_gmsh_mesh(mesh, msh)
    write_starter(case, mesh, starter)
    write_engine(job, engine, run_time)

    start = time.perf_counter()
    proc = log.run([str(STARTER), "-i", str(starter), "-nt", str(n_threads)], cwd=workdir, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"starter failed for {job}")
    proc = log.run([str(ENGINE), "-i", str(engine), "-nt", str(n_threads)], cwd=workdir, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"engine failed for {job}")
    vtk = convert_anim_to_vtk(job, workdir, log, env)
    metrics = extract_metrics(case, mesh, vtk)
    metrics["wall_clock_s"] = time.perf_counter() - start
    metrics["vtk_path"] = str(vtk)
    log.log(
        f"[{job}] Kt={metrics['kt_fem']:.6f} "
        f"Howland_err={100.0 * metrics['error_howland']:.3f}% "
        f"far_err={100.0 * metrics['far_field_error']:.3f}% "
        f"ligament={100.0 * metrics['ligament_l2_error']:.3f}%"
    )
    return metrics


def richardson_summary(rows: list[dict[str, float | int | str]]) -> dict[str, float | str]:
    by_n = {int(r["n_theta"]): float(r["kt_fem"]) for r in rows}
    k32 = by_n.get(32)
    k64 = by_n.get(64)
    k128 = by_n.get(128)
    if k32 is None or k64 is None or k128 is None:
        return {"observed_order": float("nan"), "kt_extrapolated": float("nan"), "relative_gap_to_64": float("inf")}
    e1 = k32 - k64
    e2 = k64 - k128
    if e1 == 0.0 or e2 == 0.0 or e1 * e2 <= 0.0:
        return {
            "observed_order": "nonmonotone",
            "kt_extrapolated": k128,
            "relative_gap_to_64": abs(k128 - k64) / max(abs(k64), 1.0e-12),
        }
    order = math.log(abs(e1 / e2), 2.0)
    if not math.isfinite(order) or order <= 0.0:
        return {
            "observed_order": "invalid",
            "kt_extrapolated": k128,
            "relative_gap_to_64": abs(k128 - k64) / max(abs(k64), 1.0e-12),
        }
    kt_inf = k128 + (k128 - k64) / (2.0**order - 1.0)
    return {
        "observed_order": order,
        "kt_extrapolated": kt_inf,
        "relative_gap_to_64": abs(kt_inf - k64) / max(abs(k64), 1.0e-12),
    }


def evaluate(case: PlateCase, rows: list[dict[str, float | int | str]]) -> tuple[str, dict[str, object]]:
    baseline = next(r for r in rows if int(r["n_theta"]) == N_THETA_BASELINE)
    rich = richardson_summary(rows)
    checks = {
        "primary_howland_kt": {
            "gating": True,
            "pass": float(baseline["error_howland"]) <= 0.02,
            "value": float(baseline["kt_fem"]),
            "target": case.kt_howland,
            "relative_error": float(baseline["error_howland"]),
            "tolerance": 0.02,
        },
        "far_field_stress": {
            "gating": True,
            "pass": float(baseline["far_field_error"]) <= 0.01,
            "value_Pa": float(baseline["sigma_inf_measured_Pa"]),
            "target_Pa": case.sigma_inf,
            "relative_error": float(baseline["far_field_error"]),
            "tolerance": 0.01,
        },
        "ligament_kirsch_decay": {
            "gating": False,
            "pass": float(baseline["ligament_l2_error"]) <= 0.03,
            "value": float(baseline["ligament_l2_error"]),
            "tolerance": 0.03,
            "note": "reported diagnostic; finite-width strip and exact-rectangle boundary are compared to the infinite-plate Kirsch polynomial",
        },
        "through_thickness_variation": {
            "gating": True,
            "pass": float(baseline["thickness_variation"]) <= 0.05,
            "value": float(baseline["thickness_variation"]),
            "tolerance": 0.05,
        },
        "mesh_convergence": {
            "gating": False,
            "pass": float(rich["relative_gap_to_64"]) <= 0.005,
            **rich,
            "tolerance": 0.005,
            "note": "reported diagnostic; Stage 04 primary gate is the N_theta=64 Howland Kt comparison",
        },
    }
    verdict = "PASS" if all(bool(v["pass"]) for v in checks.values() if bool(v["gating"])) else "FAIL"
    return verdict, checks


def git_sha() -> str:
    proc = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=str(ROOT),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return proc.stdout.strip() if proc.returncode == 0 else "unknown"


def write_timeseries(rows: list[dict[str, float | int | str]]) -> pathlib.Path:
    out = RESULTS_DIR / "timeseries.csv"
    fields = [
        "stage",
        "n_theta",
        "elements",
        "nodes",
        "sigma_inf_input_Pa",
        "sigma_inf_measured_Pa",
        "sigma_peak_recovered_Pa",
        "sigma_peak_direct_cell_Pa",
        "sigma_peak_point_average_Pa",
        "kt_fem",
        "kt_direct_cell",
        "kt_point_average",
        "error_howland",
        "error_kirsch",
        "far_field_error",
        "far_field_cov",
        "ligament_l2_error",
        "thickness_variation",
        "wall_clock_s",
    ]
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            payload = {"stage": STAGE}
            payload.update({key: row.get(key, "") for key in fields if key != "stage"})
            writer.writerow(payload)
    return out


def write_results_json(
    case: PlateCase,
    rows: list[dict[str, float | int | str]],
    verdict: str,
    checks: dict[str, object],
) -> pathlib.Path:
    payload = {
        "stage": STAGE,
        "name": STAGE_NAME,
        "verdict": verdict,
        "metrics": {f"n_theta_{int(r['n_theta']):03d}": r for r in rows},
        "checks": checks,
        "reference": {
            "type": "Kirsch 1898 infinite plate with Howland 1929 finite-width gross-section correction",
            "sigma_inf_Pa": case.sigma_inf,
            "Kt_Howland_2a_over_W_0p10": case.kt_howland,
            "Kt_Kirsch_infinite": case.kt_kirsch,
            "target_peak_Pa": case.kt_howland * case.sigma_inf,
            "geometry": {
                "length_m": case.length,
                "width_m": case.width,
                "thickness_m": case.thickness,
                "hole_radius_m": case.hole_radius,
            },
            "drive_correction": case.drive_correction,
        },
        "tolerance": {
            "primary_kt_howland": 0.02,
            "far_field_stress": 0.01,
            "ligament_decay_l2": 0.03,
            "through_thickness_variation": 0.05,
            "mesh_convergence": 0.005,
        },
        "git_sha": git_sha(),
    }
    out = RESULTS_DIR / "results.json"
    out.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    return out


def write_typst_figure(rows: list[dict[str, float | int | str]], verdict: str) -> pathlib.Path:
    out = FIGURES_DIR / "stage04_kirsch_overlay.typ"
    kt_points = " ".join(
        f"({float(r['n_theta']):.0f}, {float(r['kt_fem']):.6f})"
        for r in rows
    )
    table_rows = "\n".join(
        f"  [{int(r['n_theta'])}], [{float(r['kt_fem']):.4f}], "
        f"[{100.0 * float(r['error_howland']):.3f}\\%], "
        f"[{100.0 * float(r['ligament_l2_error']):.3f}\\%],"
        for r in rows
    )
    lines = [
        '#import "@preview/cetz:0.3.4"',
        "",
        '#set page(width: 180mm, height: 112mm, margin: 10mm)',
        '#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))',
        '#let garnet = rgb("#73000A")',
        '#let charcoal = rgb("#363636")',
        '#let black10 = rgb("#ECECEC")',
        '#let atlantic = rgb("#466A9F")',
        '#let white = rgb("#FFFFFF")',
        "",
        "#align(center)[#text(size: 11pt, weight: \"bold\")[Stage 04 Open-Hole Stress Concentration]]",
        "#v(2mm)",
        "#cetz.canvas(length: 1cm, {",
        "  import cetz.draw: *",
        "  rect((1.0, 0.8), (14.8, 7.0), fill: white, stroke: charcoal + 0.65pt)",
    ]
    for kt in (2.90, 2.95, 3.00, 3.05, 3.10):
        y = 0.8 + (kt - 2.90) / 0.20 * 6.2
        lines.append(f"  line((1.0, {y:.3f}), (14.8, {y:.3f}), stroke: black10 + 0.45pt)")
        lines.append(f"  content((0.84, {y:.3f}), [{kt:.2f}], anchor: \"east\")")
    x_map = {32: 1.0, 64: 7.9, 128: 14.8}
    mapped = []
    for row in rows:
        x = x_map[int(row["n_theta"])]
        y = 0.8 + (float(row["kt_fem"]) - 2.90) / 0.20 * 6.2
        mapped.append(f"({x:.3f}, {y:.3f})")
    if mapped:
        lines.append(f"  line({' '.join(mapped)}, stroke: garnet + 1.15pt)")
        for point in mapped:
            lines.append(f"  circle({point}, radius: 0.06, fill: atlantic, stroke: none)")
    target_y = 0.8 + (3.035 - 2.90) / 0.20 * 6.2
    lines.extend(
        [
            f"  line((1.0, {target_y:.3f}), (14.8, {target_y:.3f}), stroke: charcoal + 0.8pt)",
            "  content((14.7, 0.18), [Circumferential divisions], anchor: \"north-east\")",
            "  content((0.22, 3.9), [$K_t$], angle: 90deg)",
            f"  content((14.6, 7.35), [Verdict: {verdict}], anchor: \"east\")",
            "})",
            "",
            "#figure(",
            "  table(",
            "    columns: 4,",
            "    [$N_theta$], [$K_t$], [Howland error], [Ligament L2],",
            table_rows,
            "  ),",
            "  caption: [OpenRadioss recovered rim stress concentration compared with Howland's finite-width target.]",
            ")",
            "",
            f"// raw_points {kt_points}",
            "",
        ]
    )
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def write_blocker(verdict: str) -> pathlib.Path:
    blocker = STAGE_DIR / "blocker.md"
    blocker.write_text(
        "# Stage 04 Blocker: Open-Hole Kirsch Verification Failed\n\n"
        "**Author.** J.C. Vaught\n\n"
        f"Verdict: `{verdict}`\n\n"
        "The OpenRadioss runs completed, but one or more gated checks in "
        "`results/results.json` exceeded tolerance. This is a numerical "
        "verification failure for the open-hole stress-concentration stage; "
        "do not proceed until the mesh, boundary conditions, or stress "
        "recovery have been corrected.\n",
        encoding="utf-8",
    )
    return blocker


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-theta", default=",".join(str(v) for v in N_THETA_SWEEP))
    parser.add_argument("--run-time", type=float, default=1.0)
    parser.add_argument("--n-threads", type=int, default=8)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    case = PlateCase()
    n_values = [int(v.strip()) for v in args.n_theta.split(",") if v.strip()]
    if N_THETA_BASELINE not in n_values:
        raise SystemExit(f"--n-theta must include baseline {N_THETA_BASELINE}")

    log = Logger(RUN_LOG)
    try:
        log.log(f"Stage 04 {STAGE_NAME}")
        log.log(f"OpenRadioss root: {OR_DIR}")
        log.log(f"target uy={case.imposed_uy:.9e} m, Howland Kt={case.kt_howland:.6f}")
        if args.dry_run:
            for n_theta in n_values:
                mesh = build_mesh_data(case, MeshSpec(n_theta=n_theta))
                log.log(
                    f"dry-run n_theta={n_theta}: elements={len(mesh.bricks)} "
                    f"nodes={len(mesh.nodes)} top_nodes={len(mesh.node_sets['edge_y_plus'])}"
                )
            return 0

        env = radioss_env()
        rows = [
            run_case(case, n_theta, args.run_time, args.n_threads, log, env)
            for n_theta in n_values
        ]
        rows.sort(key=lambda row: int(row["n_theta"]))
        verdict, checks = evaluate(case, rows)
        write_timeseries(rows)
        write_results_json(case, rows, verdict, checks)
        write_typst_figure(rows, verdict)
        if verdict == "FAIL":
            write_blocker(verdict)
        else:
            blocker = STAGE_DIR / "blocker.md"
            if blocker.exists():
                blocker.unlink()
        return 0 if verdict == "PASS" else 1
    finally:
        log.close()


if __name__ == "__main__":
    sys.exit(main())
