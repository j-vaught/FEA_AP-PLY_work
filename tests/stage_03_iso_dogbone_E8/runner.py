"""Stage 03 runner: isotropic ASTM E8-style dogbone tension.

This runner builds a structured HEXA8 dogbone, writes OpenRadioss decks
directly, runs a small load sweep with /MAT/PLAS_JOHNS, converts the final
animation frames to VTK, extracts gauge stress/strain with PyVista, and writes
the required stage artifacts.
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


STAGE = 3
STAGE_NAME = "iso_dogbone_E8"
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

LOAD_FACTORS = (0.50, 0.80, 1.00, 1.20)


@dataclasses.dataclass(frozen=True)
class DogboneCase:
    gauge_length: float = 0.050
    gauge_width: float = 0.0125
    thickness: float = 0.006
    reduced_length: float = 0.060
    fillet_length: float = 0.0125
    grip_width: float = 0.020
    grip_length: float = 0.050
    young: float = 2.0e11
    nu: float = 0.30
    rho: float = 7850.0
    sigma_y: float = 2.50e8
    sigma_u: float = 4.00e8
    jc_b: float = 2.75e8
    jc_n: float = 0.36

    @property
    def total_length(self) -> float:
        return self.reduced_length + 2.0 * (self.grip_length + self.fillet_length)

    @property
    def area_gauge(self) -> float:
        return self.gauge_width * self.thickness

    @property
    def yield_load(self) -> float:
        return self.sigma_y * self.area_gauge

    def width_at(self, x: float) -> float:
        """Smooth E8-style dogbone width as a function of axial coordinate."""
        half_reduced = 0.5 * self.reduced_length
        ax = abs(x)
        if ax <= half_reduced:
            return self.gauge_width
        if ax >= half_reduced + self.fillet_length:
            return self.grip_width
        s = (ax - half_reduced) / self.fillet_length
        # Cosine blend gives zero slope at the gauge and grip tangencies.
        return self.gauge_width + 0.5 * (1.0 - math.cos(math.pi * s)) * (
            self.grip_width - self.gauge_width
        )


@dataclasses.dataclass(frozen=True)
class MeshSpec:
    label: str
    n_grip: int
    n_fillet: int
    n_reduced: int
    ny: int
    nz: int


MESHES = {
    "coarse": MeshSpec("coarse", n_grip=10, n_fillet=5, n_reduced=60, ny=6, nz=3),
    "medium": MeshSpec("medium", n_grip=16, n_fillet=6, n_reduced=96, ny=8, nz=4),
}


@dataclasses.dataclass
class MeshData:
    spec: MeshSpec
    nodes: dict[int, tuple[float, float, float]]
    bricks: list[tuple[int, tuple[int, int, int, int, int, int, int, int]]]
    node_sets: dict[str, list[int]]
    gauge_elements: list[int]
    gauge_x_neg: float
    gauge_x_pos: float


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
            if len(proc.stdout) <= 4000:
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


def axis_points(case: DogboneCase, spec: MeshSpec) -> list[float]:
    half = 0.5 * case.total_length
    half_reduced = 0.5 * case.reduced_length
    left_grip_end = -half_reduced - case.fillet_length
    right_grip_start = half_reduced + case.fillet_length

    segments = [
        (-half, left_grip_end, spec.n_grip),
        (left_grip_end, -half_reduced, spec.n_fillet),
        (-half_reduced, half_reduced, spec.n_reduced),
        (half_reduced, right_grip_start, spec.n_fillet),
        (right_grip_start, half, spec.n_grip),
    ]
    xs: list[float] = []
    for start, stop, count in segments:
        pts = [start + (stop - start) * i / count for i in range(count)]
        xs.extend(pts)
    xs.append(half)
    return xs


def build_mesh_data(case: DogboneCase, spec: MeshSpec) -> MeshData:
    xs = axis_points(case, spec)
    z0 = -0.5 * case.thickness
    dz = case.thickness / spec.nz
    nodes: dict[int, tuple[float, float, float]] = {}

    def nid(i: int, j: int, k: int) -> int:
        return 1 + k * (spec.ny + 1) * len(xs) + j * len(xs) + i

    for k in range(spec.nz + 1):
        z = z0 + dz * k
        for j in range(spec.ny + 1):
            eta = -1.0 + 2.0 * j / spec.ny
            for i, x in enumerate(xs):
                y = 0.5 * case.width_at(x) * eta
                nodes[nid(i, j, k)] = (x, y, z)

    bricks: list[tuple[int, tuple[int, int, int, int, int, int, int, int]]] = []
    gauge_elements: list[int] = []
    eid = 1
    half_gauge = 0.5 * case.gauge_length
    half_reduced = 0.5 * case.reduced_length
    for k in range(spec.nz):
        for j in range(spec.ny):
            for i in range(len(xs) - 1):
                conn = (
                    nid(i, j, k),
                    nid(i + 1, j, k),
                    nid(i + 1, j + 1, k),
                    nid(i, j + 1, k),
                    nid(i, j, k + 1),
                    nid(i + 1, j, k + 1),
                    nid(i + 1, j + 1, k + 1),
                    nid(i, j + 1, k + 1),
                )
                bricks.append((eid, conn))
                cx = 0.5 * (xs[i] + xs[i + 1])
                if -half_gauge <= cx <= half_gauge:
                    gauge_elements.append(eid)
                eid += 1

    i_left = 0
    i_right = len(xs) - 1
    i_gauge_neg = min(range(len(xs)), key=lambda idx: abs(xs[idx] + half_gauge))
    i_gauge_pos = min(range(len(xs)), key=lambda idx: abs(xs[idx] - half_gauge))
    if abs(xs[i_gauge_neg] + half_gauge) > 1.0e-10 or abs(xs[i_gauge_pos] - half_gauge) > 1.0e-10:
        raise RuntimeError("mesh does not contain exact gauge endpoint planes")

    fixed = [nid(i_left, j, k) for k in range(spec.nz + 1) for j in range(spec.ny + 1)]
    loaded = [nid(i_right, j, k) for k in range(spec.nz + 1) for j in range(spec.ny + 1)]
    gauge_neg = [nid(i_gauge_neg, j, k) for k in range(spec.nz + 1) for j in range(spec.ny + 1)]
    gauge_pos = [nid(i_gauge_pos, j, k) for k in range(spec.nz + 1) for j in range(spec.ny + 1)]

    return MeshData(
        spec=spec,
        nodes=nodes,
        bricks=bricks,
        node_sets={
            "fixed": fixed,
            "loaded": loaded,
            "gauge_neg": gauge_neg,
            "gauge_pos": gauge_pos,
            "all_nodes": list(nodes),
        },
        gauge_elements=gauge_elements,
        gauge_x_neg=xs[i_gauge_neg],
        gauge_x_pos=xs[i_gauge_pos],
    )


def write_gmsh_mesh(mesh: MeshData, out_msh: pathlib.Path) -> None:
    import gmsh  # type: ignore[import-not-found]

    gmsh.initialize(["-nopopup"])
    try:
        gmsh.model.add(f"stage03_{mesh.spec.label}")
        gmsh.model.addDiscreteEntity(3, 1)
        node_tags = list(mesh.nodes)
        coords: list[float] = []
        for node_id in node_tags:
            coords.extend(mesh.nodes[node_id])
        elem_tags = [eid for eid, _ in mesh.bricks]
        elem_conn = [node for _, conn in mesh.bricks for node in conn]
        gmsh.model.mesh.addNodes(3, 1, node_tags, coords)
        gmsh.model.mesh.addElementsByType(1, 5, elem_tags, elem_conn)
        gmsh.model.addPhysicalGroup(3, [1], tag=1, name="DOGBONE_SOLID")
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(str(out_msh))
    finally:
        gmsh.finalize()


def node_group_block(group_id: int, name: str, nodes: Iterable[int]) -> list[str]:
    ids = list(nodes)
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for i in range(0, len(ids), 10):
        lines.append(fmt_i(*ids[i : i + 10]))
    return lines


def brick_group_block(group_id: int, name: str, elems: Iterable[int]) -> list[str]:
    ids = list(elems)
    lines = [f"/GRBRIC/BRIC/{group_id}", name]
    for i in range(0, len(ids), 10):
        lines.append(fmt_i(*ids[i : i + 10]))
    return lines


def write_starter(
    case: DogboneCase,
    mesh: MeshData,
    load_factor: float,
    out_rad: pathlib.Path,
    run_time: float,
    damping_alpha: float,
) -> None:
    job = out_rad.name.removesuffix("_0000.rad")
    total_load = load_factor * case.yield_load
    force_per_node = total_load / len(mesh.node_sets["loaded"])
    ramp_time = 0.8 * run_time
    lines: list[str] = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        fmt_i(2019, 0),
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 03 A36 dogbone load factor {load_factor:g}",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 0) + f"{0:20d}" + f"{2:40d}",
        "/RANDOM",
        fmt_f(0.0) + f"{0:20d}",
        "/SPMD",
        fmt_i(0, 0) + f"{0:20d}{1:20d}",
        "/SHFRA/V4",
        "/DAMP/1",
        "whole_dogbone_mass_damping",
        "#              alpha                beta  grnod_ID   skew_ID              Tstart               Tstop",
        fmt_f(damping_alpha, 0.0) + fmt_i(300, 0) + fmt_f(0.0, run_time),
        "/MAT/PLAS_JOHNS/1",
        "A36_steel_JC_rate_off",
        "#        Init. dens.          Ref. dens.",
        fmt_f(case.rho, case.rho),
        "#                  E                  Nu",
        fmt_f(case.young, case.nu),
        "#                  a                   b                   n             Eps_max              sigmax",
        fmt_f(case.sigma_y, case.jc_b, case.jc_n, 0.0, 0.0),
        "#                  c                EPS0       Icc   Fsmooth               F_CUT",
        fmt_f(0.0, 1.0) + fmt_i(0, 0) + fmt_f(0.0),
        "#                  m              T_melt               rhoCp                 T_i",
        fmt_f(0.0, 0.0, 0.0, 0.0),
        "/NODE",
    ]
    for nid, (x, y, z) in mesh.nodes.items():
        lines.append(f"{nid:10d}{x:20.12g}{y:20.12g}{z:20.12g}")
    lines.extend(["/PART/1", "dogbone_part", fmt_i(1, 1, 0), "/BRICK/1"])
    for eid, conn in mesh.bricks:
        lines.append(fmt_i(eid, *conn))
    lines.extend(
        [
            "/PROP/SOLID/1",
            "heph_constant_pressure",
            "#   Isolid    Ismstr               Icpre               Inpts    Itetra    Iframe                  dn",
            fmt_i(24, 0) + f"{1:20d}{0:20d}{0:10d}{0:10d}{0:20d}",
            "#                q_a                 q_b                   h            LAMBDA_V                MU_V",
            fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
            "#             dt_min   istrain      IHKT",
            fmt_f(0.0) + fmt_i(0, 0),
            "/BCS/1",
            "fixed_grip",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   111 000{0:10d}{100:10d}",
            "/BCS/2",
            "loaded_grip_lateral_restraint",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   011 000{0:10d}{101:10d}",
            "/FUNCT/1",
            "smooth_force_ramp",
            "#                  X                   Y",
            fmt_f(0.0, 0.0),
            fmt_f(ramp_time, 1.0),
            fmt_f(run_time, 1.0),
            "/CLOAD/1",
            "loaded_grip_tension_x",
            "#funct_IDT       Dir   skew_ID sensor_ID  grnod_ID                       Ascalex             Fscaley",
            f"{1:10d}{'X':>10}{0:10d}{0:10d}{101:10d}{1.0:30.12g}{force_per_node:20.12g}",
        ]
    )
    lines.extend(node_group_block(100, "fixed_grip", mesh.node_sets["fixed"]))
    lines.extend(node_group_block(101, "loaded_grip", mesh.node_sets["loaded"]))
    lines.extend(node_group_block(102, "gauge_neg", mesh.node_sets["gauge_neg"]))
    lines.extend(node_group_block(103, "gauge_pos", mesh.node_sets["gauge_pos"]))
    lines.extend(node_group_block(300, "all_nodes", mesh.node_sets["all_nodes"]))
    lines.extend(brick_group_block(200, "gauge_bricks", mesh.gauge_elements))
    lines.extend(
        [
            "/TH/NODE/1",
            "gauge_negative_plane",
            "DEF",
        ]
    )
    for node_id in mesh.node_sets["gauge_neg"]:
        lines.append(f"{node_id:10d}{0:10d}gneg_{node_id}")
    lines.extend(["/TH/NODE/2", "gauge_positive_plane", "DEF"])
    for node_id in mesh.node_sets["gauge_pos"]:
        lines.append(f"{node_id:10d}{0:10d}gpos_{node_id}")
    lines.extend(["/END", ""])
    out_rad.write_text("\n".join(lines), encoding="utf-8")


def write_engine(job: str, out_rad: pathlib.Path, run_time: float) -> None:
    lines = [
        "#RADIOSS ENGINE",
        "/ANIM/DT",
        fmt_f(run_time, run_time),
        "/ANIM/VECT/DISP",
        "/ANIM/VECT/VEL",
        "/ANIM/BRICK/TENS/STRESS/ALL",
        "/ANIM/BRICK/TENS/STRAIN/ALL",
        "/ANIM/ELEM/EPSP",
        "/ANIM/GZIP",
        "/TFILE/4",
        fmt_f(run_time / 100.0),
        "/RFILE",
        fmt_i(5000),
        "/PRINT/-100/55",
        f"/RUN/{job}/1",
        fmt_f(run_time),
        "/VERS/2019",
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


def extract_metrics(case: DogboneCase, mesh: MeshData, vtk_path: pathlib.Path, load_factor: float) -> dict[str, float]:
    import numpy as np
    import pyvista as pv  # type: ignore[import-not-found]

    grid = pv.read(str(vtk_path))
    node_ids = grid.point_data.get("NODE_ID")
    disp = grid.point_data.get("Displacement")
    if node_ids is None or disp is None:
        raise KeyError(f"NODE_ID or Displacement missing from {vtk_path}")
    node_to_idx = {int(node_id): i for i, node_id in enumerate(node_ids)}
    gneg_idx = [node_to_idx[n] for n in mesh.node_sets["gauge_neg"] if n in node_to_idx]
    gpos_idx = [node_to_idx[n] for n in mesh.node_sets["gauge_pos"] if n in node_to_idx]
    if not gneg_idx or not gpos_idx:
        raise KeyError("gauge endpoint nodes missing from VTK")
    u_neg = float(np.mean(disp[gneg_idx, 0]))
    u_pos = float(np.mean(disp[gpos_idx, 0]))
    eps_disp = (u_pos - u_neg) / (mesh.gauge_x_pos - mesh.gauge_x_neg)

    elem_ids = grid.cell_data.get("ELEMENT_ID")
    if elem_ids is None:
        raise KeyError("ELEMENT_ID missing from VTK")
    gauge_mask = np.isin(elem_ids, np.asarray(mesh.gauge_elements, dtype=elem_ids.dtype))
    stress = _cell_array(grid, "Strs")
    strain = _cell_array(grid, "Stra")
    sigma_x = np.asarray(stress[gauge_mask, 0], dtype=float)
    eps_x = np.asarray(strain[gauge_mask, 0], dtype=float)
    sigma_mean = float(np.mean(sigma_x))
    sigma_std = float(np.std(sigma_x))
    eps_mean = float(np.mean(eps_x))
    sigma_ref = load_factor * case.sigma_y
    eps_ref_elastic = sigma_ref / case.young
    epsp_est = max(0.0, abs(eps_disp) - abs(sigma_mean) / case.young)
    sigma_jc = case.sigma_y + case.jc_b * (epsp_est ** case.jc_n) if epsp_est > 0.0 else case.sigma_y
    return {
        "load_factor": load_factor,
        "applied_load_N": load_factor * case.yield_load,
        "sigma_ref_Pa": sigma_ref,
        "sigma_mean_Pa": sigma_mean,
        "sigma_std_Pa": sigma_std,
        "eps_disp": eps_disp,
        "eps_mean": eps_mean,
        "eps_ref_elastic": eps_ref_elastic,
        "epsp_est": epsp_est,
        "sigma_jc_from_epsp_Pa": sigma_jc,
        "stress_error": abs(abs(sigma_mean) - sigma_ref) / sigma_ref,
        "strain_error": abs(abs(eps_mean) - abs(eps_disp)) / max(abs(eps_disp), 1.0e-12),
        "uniformity": abs(sigma_std) / max(abs(sigma_mean), 1.0),
    }


def run_case(
    case: DogboneCase,
    mesh_label: str,
    load_factor: float,
    run_time: float,
    damping_alpha: float,
    n_threads: int,
    log: Logger,
    env: dict[str, str],
) -> dict[str, float]:
    mesh = build_mesh_data(case, MESHES[mesh_label])
    factor_label = f"{load_factor:.2f}".replace(".", "p")
    job = f"stage03_{mesh_label}_lf{factor_label}"
    workdir = RUNS_DIR / job
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    starter = workdir / f"{job}_0000.rad"
    engine = workdir / f"{job}_0001.rad"
    msh = workdir / f"{job}.msh"
    log.log(
        f"[{job}] elements={len(mesh.bricks)} nodes={len(mesh.nodes)} "
        f"load={load_factor * case.yield_load:.6g} N"
    )
    write_gmsh_mesh(mesh, msh)
    write_starter(case, mesh, load_factor, starter, run_time, damping_alpha)
    write_engine(job, engine, run_time)
    start = time.perf_counter()
    proc = log.run([str(STARTER), "-i", str(starter), "-nt", str(n_threads)], cwd=workdir, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"starter failed for {job}")
    proc = log.run([str(ENGINE), "-i", str(engine), "-nt", str(n_threads)], cwd=workdir, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"engine failed for {job}")
    vtk = convert_anim_to_vtk(job, workdir, log, env)
    metrics = extract_metrics(case, mesh, vtk, load_factor)
    metrics["wall_clock_s"] = time.perf_counter() - start
    metrics["elements"] = len(mesh.bricks)
    metrics["nodes"] = len(mesh.nodes)
    log.log(
        f"[{job}] sigma={metrics['sigma_mean_Pa']:.6e} Pa "
        f"ref={metrics['sigma_ref_Pa']:.6e} Pa "
        f"stress_err={100.0 * metrics['stress_error']:.3f}% "
        f"strain_err={100.0 * metrics['strain_error']:.3f}%"
    )
    return metrics


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


def evaluate(rows: list[dict[str, float]], case: DogboneCase) -> tuple[str, dict[str, object]]:
    by_factor = {round(r["load_factor"], 2): r for r in rows}
    elastic = by_factor[0.80]
    yield_row = by_factor[1.00]
    post = by_factor[1.20]
    e_rows = [by_factor[0.50], by_factor[0.80]]
    ds = e_rows[1]["sigma_mean_Pa"] - e_rows[0]["sigma_mean_Pa"]
    de = e_rows[1]["eps_disp"] - e_rows[0]["eps_disp"]
    e_app = abs(ds / de) if abs(de) > 0.0 else float("nan")
    c1 = elastic["uniformity"] <= 0.01
    c2 = elastic["stress_error"] <= 0.01
    c3 = elastic["strain_error"] <= 0.01
    c4 = abs(abs(yield_row["sigma_mean_Pa"]) - case.sigma_y) / case.sigma_y <= 0.01
    s1 = math.isfinite(e_app) and abs(e_app - case.young) / case.young <= 0.01
    post_plastic = post["epsp_est"] > 1.0e-5 and abs(post["sigma_mean_Pa"]) > case.sigma_y
    verdict = "PASS" if all((c1, c2, c3, c4, s1, post_plastic)) else "FAIL"
    checks = {
        "C1_uniformity": {"pass": c1, "value": elastic["uniformity"], "tolerance": 0.01},
        "C2_stress_vs_load_over_area": {"pass": c2, "value": elastic["stress_error"], "tolerance": 0.01},
        "C3_strain_vtk_vs_extensometer": {"pass": c3, "value": elastic["strain_error"], "tolerance": 0.01},
        "C4_yield_onset_stress": {
            "pass": c4,
            "observed_Pa": yield_row["sigma_mean_Pa"],
            "target_Pa": case.sigma_y,
            "relative_error": abs(abs(yield_row["sigma_mean_Pa"]) - case.sigma_y) / case.sigma_y,
            "tolerance": 0.01,
        },
        "S1_apparent_E": {
            "pass": s1,
            "observed_Pa": e_app,
            "target_Pa": case.young,
            "relative_error": abs(e_app - case.young) / case.young if math.isfinite(e_app) else float("inf"),
            "tolerance": 0.01,
        },
        "S3_post_yield_plasticity_observed": {
            "pass": post_plastic,
            "epsp_est": post["epsp_est"],
            "sigma_mean_Pa": post["sigma_mean_Pa"],
        },
    }
    return verdict, checks


def write_timeseries(rows: list[dict[str, float]]) -> pathlib.Path:
    out = RESULTS_DIR / "timeseries.csv"
    fields = [
        "stage",
        "load_factor",
        "applied_load_N",
        "sigma_ref_Pa",
        "sigma_mean_Pa",
        "sigma_std_Pa",
        "eps_disp",
        "eps_mean",
        "eps_ref_elastic",
        "epsp_est",
        "stress_error",
        "strain_error",
        "uniformity",
        "wall_clock_s",
    ]
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow({"stage": STAGE, **{k: r.get(k, "") for k in fields if k != "stage"}})
    return out


def write_results_json(rows: list[dict[str, float]], verdict: str, checks: dict[str, object], case: DogboneCase) -> pathlib.Path:
    payload = {
        "stage": STAGE,
        "name": STAGE_NAME,
        "verdict": verdict,
        "metrics": {
            f"load_factor_{r['load_factor']:.2f}": r for r in rows
        },
        "checks": checks,
        "reference": {
            "type": "closed-form uniaxial stress/strain plus input LAW2 yield stress",
            "E_Pa": case.young,
            "sigma_y_Pa": case.sigma_y,
            "area_gauge_m2": case.area_gauge,
            "yield_load_N": case.yield_load,
        },
        "tolerance": {
            "stress_uniformity": 0.01,
            "stress_vs_load_over_area": 0.01,
            "strain_vtk_vs_extensometer": 0.01,
            "yield_stress": 0.01,
            "apparent_E": 0.01,
        },
        "git_sha": git_sha(),
    }
    out = RESULTS_DIR / "results.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return out


def write_typst_figure(rows: list[dict[str, float]], verdict: str) -> pathlib.Path:
    out = FIGURES_DIR / "stage03_stress_strain.typ"
    pts = " ".join(
        f"({r['eps_disp']:.6g}, {r['sigma_mean_Pa'] * 1.0e-6:.6g})"
        for r in rows
    )
    table_rows = "\n".join(
        f"  [{r['load_factor']:.2f}], [{r['eps_disp']:.6f}], [{r['sigma_mean_Pa'] * 1.0e-6:.2f}], [{100*r['stress_error']:.3f}\\%],"
        for r in rows
    )
    lines = [
        '#import "@preview/cetz:0.3.4"',
        "",
        '#set page(width: 180mm, height: 118mm, margin: 10mm)',
        '#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))',
        '#let garnet = rgb("#73000A")',
        '#let charcoal = rgb("#363636")',
        '#let black10 = rgb("#ECECEC")',
        '#let atlantic = rgb("#466A9F")',
        '#let white = rgb("#FFFFFF")',
        "",
        "#align(center)[#text(size: 11pt, weight: \"bold\")[Stage 03 A36 Dogbone Stress-Strain]]",
        "#v(2mm)",
        "#cetz.canvas(length: 1cm, {",
        "  import cetz.draw: *",
        "  rect((1.0, 0.8), (14.6, 7.0), fill: white, stroke: charcoal + 0.65pt)",
    ]
    for stress in (0, 100, 200, 300, 400):
        y = 0.8 + stress / 400 * 6.2
        lines.append(f"  line((1.0, {y:.3f}), (14.6, {y:.3f}), stroke: black10 + 0.45pt)")
        lines.append(f"  content((0.84, {y:.3f}), [{stress}], anchor: \"east\")")
    if rows:
        max_eps = max(abs(r["eps_disp"]) for r in rows) * 1.1
        mapped = []
        for r in rows:
            x = 1.0 + abs(r["eps_disp"]) / max_eps * 13.6
            y = 0.8 + abs(r["sigma_mean_Pa"]) * 1.0e-6 / 400 * 6.2
            mapped.append(f"({x:.3f}, {y:.3f})")
        lines.append(f"  line({' '.join(mapped)}, stroke: garnet + 1.15pt)")
        for r, p in zip(rows, mapped):
            lines.append(f"  circle({p}, radius: 0.06, fill: atlantic, stroke: none)")
    lines.extend(
        [
            "  content((7.8, 0.18), [Gauge strain], anchor: \"north\")",
            "  content((0.22, 3.9), [Gauge stress (MPa)], angle: 90deg)",
            f"  content((14.5, 7.35), [Verdict: {verdict}], anchor: \"east\")",
            "})",
            "",
            "#figure(",
            "  table(",
            "    columns: 4,",
            "    [Load factor], [Gauge strain], [Stress (MPa)], [Stress error],",
            table_rows,
            "  ),",
            "  caption: [Stage 03 extracted gauge response from OpenRadioss final frames.]",
            ")",
            "",
        ]
    )
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", choices=sorted(MESHES), default="medium")
    parser.add_argument("--load-factors", default=",".join(f"{v:g}" for v in LOAD_FACTORS))
    parser.add_argument("--run-time", type=float, default=0.002)
    parser.add_argument("--damping-alpha", type=float, default=10000.0)
    parser.add_argument("--n-threads", type=int, default=4)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    case = DogboneCase()
    factors = [float(v.strip()) for v in args.load_factors.split(",") if v.strip()]
    required = {0.50, 0.80, 1.00, 1.20}
    if not required.issubset({round(v, 2) for v in factors}):
        raise SystemExit(f"load factors must include {sorted(required)}")

    log = Logger(RUN_LOG)
    try:
        log.log(f"Stage 03 {STAGE_NAME}")
        log.log(f"OpenRadioss root: {OR_DIR}")
        log.log(f"Gauge area={case.area_gauge:.8e} m^2, yield load={case.yield_load:.6f} N")
        if args.dry_run:
            mesh = build_mesh_data(case, MESHES[args.mesh])
            log.log(f"dry-run mesh={args.mesh}: {len(mesh.bricks)} elements, {len(mesh.nodes)} nodes")
            return 0

        env = radioss_env()
        rows = [
            run_case(case, args.mesh, lf, args.run_time, args.damping_alpha, args.n_threads, log, env)
            for lf in factors
        ]
        rows.sort(key=lambda r: r["load_factor"])
        verdict, checks = evaluate(rows, case)
        write_timeseries(rows)
        write_results_json(rows, verdict, checks, case)
        write_typst_figure(rows, verdict)
        if verdict == "FAIL":
            blocker = STAGE_DIR / "blocker.md"
            blocker.write_text(
                "# Stage 03 Blocker: Dogbone Tensile Verification Failed\n\n"
                "**Author.** J.C. Vaught\n\n"
                f"Verdict: `{verdict}`\n\n"
                "The OpenRadioss run completed, but at least one gated stress, strain, "
                "yield-onset, or apparent-modulus check exceeded tolerance. See "
                "`results/results.json` and `run.log` for the measured values.\n",
                encoding="utf-8",
            )
        return 0 if verdict == "PASS" else 1
    finally:
        log.close()


if __name__ == "__main__":
    sys.exit(main())
