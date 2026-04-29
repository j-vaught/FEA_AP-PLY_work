"""Stage 02 runner: cantilever large-deflection geometric nonlinearity.

Builds a structured solid HEXA8 cantilever, writes OpenRadioss decks, runs a
nonlinear implicit solve when available, falls back to explicit dynamic
relaxation on request/auto-failure, converts the final animation frame to VTK,
and compares tip motion with the Bisshopp-Drucker elastica reference.
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

from scipy import optimize, special  # type: ignore[import-not-found]


STAGE = 2
STAGE_NAME = "cantilever_large"
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

PASS_ALPHAS = (1.0, 3.0, 5.0)
PLOT_ALPHAS = (0.05, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0)
PASS_TOL = 0.02


@dataclasses.dataclass(frozen=True)
class CantileverCase:
    length: float = 1.000
    width_z: float = 0.0250
    height_y: float = 0.0030
    young: float = 2.0e11
    nu: float = 0.30
    rho: float = 7850.0

    @property
    def inertia_z(self) -> float:
        return self.width_z * self.height_y**3 / 12.0

    @property
    def ei(self) -> float:
        return self.young * self.inertia_z

    def load_for_alpha(self, alpha: float) -> float:
        return alpha * self.ei / self.length**2


@dataclasses.dataclass(frozen=True)
class ElasticaRef:
    alpha: float
    theta_L: float
    delta_y_over_L: float
    delta_x_over_L: float


def _alpha_of_theta(theta_L: float) -> float:
    k = math.sqrt((1.0 + math.sin(theta_L)) / 2.0)
    sin_phi0 = max(-1.0, min(1.0, 1.0 / (k * math.sqrt(2.0))))
    phi0 = math.asin(sin_phi0)
    big_K = float(special.ellipk(k * k))
    F_phi0 = float(special.ellipkinc(phi0, k * k))
    return (big_K - F_phi0) ** 2


def elastica_reference(alpha: float) -> ElasticaRef:
    if alpha <= 0.0:
        return ElasticaRef(alpha, 0.0, 0.0, 0.0)
    theta_L = optimize.brentq(
        lambda th: _alpha_of_theta(th) - alpha,
        1.0e-6,
        math.pi / 2.0 - 1.0e-6,
        xtol=1.0e-12,
    )
    k = math.sqrt((1.0 + math.sin(theta_L)) / 2.0)
    phi0 = math.asin(max(-1.0, min(1.0, 1.0 / (k * math.sqrt(2.0)))))
    big_K = float(special.ellipk(k * k))
    big_E = float(special.ellipe(k * k))
    F_phi0 = float(special.ellipkinc(phi0, k * k))
    E_phi0 = float(special.ellipeinc(phi0, k * k))
    sqrt_alpha = big_K - F_phi0
    delta_x_over_L = 1.0 - math.sqrt(2.0 * math.sin(theta_L)) / sqrt_alpha
    delta_y_over_L = 1.0 - 2.0 * (big_E - E_phi0) / sqrt_alpha
    return ElasticaRef(alpha, theta_L, delta_y_over_L, delta_x_over_L)


@dataclasses.dataclass(frozen=True)
class MeshSpec:
    label: str
    nx: int
    ny_height: int
    nz_width: int


MESHES = {
    "coarse": MeshSpec("coarse", nx=40, ny_height=3, nz_width=5),
    "baseline": MeshSpec("baseline", nx=80, ny_height=4, nz_width=6),
    "fine": MeshSpec("fine", nx=160, ny_height=4, nz_width=8),
}


@dataclasses.dataclass
class MeshData:
    spec: MeshSpec
    nodes: dict[int, tuple[float, float, float]]
    bricks: list[tuple[int, tuple[int, int, int, int, int, int, int, int]]]
    node_sets: dict[str, list[int]]
    probe_nodes: list[int]


@dataclasses.dataclass
class SweepRow:
    alpha: float
    mesh: str
    solver_mode: str
    elements: int
    nodes: int
    dx_fem_over_L: float
    dy_fem_over_L: float
    dx_ref_over_L: float
    dy_ref_over_L: float
    err_x: float
    err_y: float
    passed: bool
    gated: bool
    wall_clock_s: float


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


def parse_csv_list(text: str, cast: type = str) -> tuple:
    return tuple(cast(item.strip()) for item in text.split(",") if item.strip())


def build_mesh_data(case: CantileverCase, spec: MeshSpec) -> MeshData:
    xs = [case.length * i / spec.nx for i in range(spec.nx + 1)]
    ys = [-case.height_y / 2.0 + case.height_y * j / spec.ny_height for j in range(spec.ny_height + 1)]
    zs = [-case.width_z / 2.0 + case.width_z * k / spec.nz_width for k in range(spec.nz_width + 1)]
    nodes: dict[int, tuple[float, float, float]] = {}

    def nid(i: int, j: int, k: int) -> int:
        return 1 + k * (spec.ny_height + 1) * (spec.nx + 1) + j * (spec.nx + 1) + i

    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                nodes[nid(i, j, k)] = (x, y, z)

    bricks: list[tuple[int, tuple[int, int, int, int, int, int, int, int]]] = []
    eid = 1
    for k in range(spec.nz_width):
        for j in range(spec.ny_height):
            for i in range(spec.nx):
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
                eid += 1

    clamp = [nid(0, j, k) for k in range(spec.nz_width + 1) for j in range(spec.ny_height + 1)]
    tip = [nid(spec.nx, j, k) for k in range(spec.nz_width + 1) for j in range(spec.ny_height + 1)]
    tip_by_centroid_distance = sorted(tip, key=lambda node_id: nodes[node_id][1] ** 2 + nodes[node_id][2] ** 2)
    if math.hypot(nodes[tip_by_centroid_distance[0]][1], nodes[tip_by_centroid_distance[0]][2]) <= 1.0e-12:
        probe_nodes = [tip_by_centroid_distance[0]]
    else:
        probe_nodes = tip_by_centroid_distance[:4]
    return MeshData(
        spec=spec,
        nodes=nodes,
        bricks=bricks,
        node_sets={
            "clamp": clamp,
            "tip": tip,
            "probe": probe_nodes,
            "all_nodes": list(nodes),
        },
        probe_nodes=probe_nodes,
    )


def write_gmsh_mesh(mesh: MeshData, out_msh: pathlib.Path) -> None:
    import gmsh  # type: ignore[import-not-found]

    gmsh.initialize(["-nopopup"])
    try:
        gmsh.model.add(mesh.spec.label)
        gmsh.model.addDiscreteEntity(3, 1)
        node_tags = list(mesh.nodes)
        coords: list[float] = []
        for nid in node_tags:
            coords.extend(mesh.nodes[nid])
        gmsh.model.mesh.addNodes(3, 1, node_tags, coords)
        elem_tags = [eid for eid, _ in mesh.bricks]
        elem_conn = [node for _, conn in mesh.bricks for node in conn]
        gmsh.model.mesh.addElementsByType(1, 5, elem_tags, elem_conn)
        gmsh.write(str(out_msh))
    finally:
        gmsh.finalize()


def node_group_block(group_id: int, name: str, ids: Iterable[int]) -> list[str]:
    ids = list(ids)
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for i in range(0, len(ids), 10):
        lines.append(fmt_i(*ids[i : i + 10]))
    return lines


def write_starter(
    case: CantileverCase,
    mesh: MeshData,
    alpha: float,
    out_rad: pathlib.Path,
    solver_mode: str,
    run_time: float,
    damping_alpha: float,
) -> None:
    job = out_rad.name.removesuffix("_0000.rad")
    total_load = case.load_for_alpha(alpha)
    force_per_tip_node = -total_load / len(mesh.node_sets["tip"])
    ramp_end = run_time if solver_mode == "implicit" else min(run_time, 0.20 * run_time)
    ramp_points = [(0.0, 0.0), (ramp_end, 1.0)]
    if not math.isclose(ramp_end, run_time):
        ramp_points.append((run_time, 1.0))
    lines: list[str] = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        fmt_i(2019, 0),
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 02 alpha {alpha:g} {solver_mode} large-deflection cantilever",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(14, 11) + f"{0:20d}" + f"{2:40d}",
        "/RANDOM",
        fmt_f(0.0) + f"{0:20d}",
        "/SPMD",
        fmt_i(0, 0) + f"{0:20d}{1:20d}",
        "/SHFRA/V4",
    ]
    if solver_mode == "explicit":
        lines.extend(
            [
                "/DAMP/1",
                "whole_cantilever_mass_damping",
                "#              alpha                beta  grnod_ID   skew_ID              Tstart               Tstop",
                fmt_f(damping_alpha, 0.0) + fmt_i(300, 0) + fmt_f(0.0, run_time),
            ]
        )
    lines.extend(
        [
            "/MAT/ELAST/1",
            "linear_spring_steel",
            fmt_f(case.rho, 0.0),
            fmt_f(case.young, case.nu),
            "/NODE",
        ]
    )
    for nid, (x, y, z) in mesh.nodes.items():
        lines.append(f"{nid:10d}{x:20.12g}{y:20.12g}{z:20.12g}")
    lines.extend(["/PART/1", "cantilever_part", fmt_i(1, 1, 0), "/BRICK/1"])
    for eid, conn in mesh.bricks:
        lines.append(fmt_i(eid, *conn))
    lines.extend(
        [
            "/PROP/SOLID/1",
            "solid14_corotational",
            "#   Isolid    Ismstr               Icpre               Inpts    Itetra    Iframe                  dn",
            fmt_i(14, 11) + f"{1:20d}{0:20d}{0:10d}{0:10d}{0:20d}",
            "#                q_a                 q_b                   h            LAMBDA_V                MU_V",
            fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
            "#             dt_min   istrain      IHKT",
            fmt_f(0.0) + fmt_i(0, 0),
            "/BCS/1",
            "clamp_x0",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   111 000{0:10d}{100:10d}",
        ]
    )
    lines.extend(
        [
            "/FUNCT/1",
            "dead_load_ramp",
            "#                  X                   Y",
        ]
    )
    for t, scale in ramp_points:
        lines.append(fmt_f(t, scale))
    lines.extend(
        [
            "/CLOAD/1",
            "tip_dead_load_y",
            "#funct_IDT       Dir   skew_ID sensor_ID  grnod_ID                       Ascalex             Fscaley",
            f"{1:10d}{'Y':>10}{0:10d}{0:10d}{200:10d}{1.0:30.12g}{force_per_tip_node:20.12g}",
        ]
    )
    lines.extend(node_group_block(100, "clamp_x0", mesh.node_sets["clamp"]))
    lines.extend(node_group_block(200, "tip_face", mesh.node_sets["tip"]))
    lines.extend(node_group_block(201, "tip_probe", mesh.node_sets["probe"]))
    lines.extend(node_group_block(300, "all_nodes", mesh.node_sets["all_nodes"]))
    lines.extend(
        [
            "/TH/NODE/1",
            "tip_probe_displacement",
            "#     var1      var2      var3      var4      var5      var6      var7      var8      var9     var10",
            "DEF",
            "#    NODid     Iskew                                           NODname",
        ]
    )
    th_probe_nodes = mesh.probe_nodes
    for node_id in th_probe_nodes:
        lines.append(f"{node_id:10d}{0:10d}probe_{node_id}")
    if solver_mode == "implicit":
        lines.extend(
            [
                "/TH/NODE/2",
                "clamp_reaction_y",
                "#     var1      var2      var3      var4      var5      var6      var7      var8      var9     var10",
                "REACY",
                "#    NODid     Iskew                                           NODname",
            ]
        )
        for node_id in mesh.node_sets["clamp"]:
            lines.append(f"{node_id:10d}{0:10d}clamp_{node_id}")
    lines.extend(["/END", ""])
    out_rad.write_text("\n".join(lines), encoding="utf-8")


def write_engine(job: str, out_rad: pathlib.Path, solver_mode: str, run_time: float, dt_noda: float) -> None:
    if solver_mode == "implicit":
        lines = [
            "#RADIOSS ENGINE",
            "/TITLE",
            f"Stage 02 {job} nonlinear implicit",
            "/VERS/2019",
            f"/RUN/{job}/1",
            fmt_f(run_time),
            "/ANIM/DT",
            fmt_f(run_time, run_time),
            "/TFILE/4",
            fmt_f(run_time / 200.0),
            "/RFILE",
            fmt_i(5000),
            "/PRINT/-1",
            "/MON/ON",
            "/ANIM/VECT/DISP",
            "/ANIM/VECT/VEL",
            "/ANIM/BRICK/TENS/STRESS/ALL",
            "/ANIM/BRICK/TENS/STRAIN/ALL",
            "/ANIM/GZIP",
            "/IMPL/PRINT/NONL/-1",
            "/IMPL/NONLIN/KTCON",
            "/IMPL/NONLIN/1",
            f"{25:10d}{2:10d}{0.02:20.12g}",
            "/IMPL/SOLVER/2",
            f"{0:10d}{0:10d}{0:10d}{0.0:20.12g}",
            "/IMPL/DTINI",
            fmt_f(0.02 * run_time),
            "/IMPL/DT/STOP",
            fmt_f(1.0e-8 * run_time, 0.20 * run_time),
            "/IMPL/DT/2",
            f"{18:10d}{0.0:10.1f}{30:10d}{0.67:10.2f}{1.35:10.2f}",
            "/END/ENGINE",
            "",
        ]
    else:
        lines = [
            "#RADIOSS ENGINE",
            "/ANIM/DT",
            fmt_f(run_time, run_time),
            "/ANIM/VECT/DISP",
            "/ANIM/VECT/VEL",
            "/ANIM/BRICK/TENS/STRESS/ALL",
            "/ANIM/BRICK/TENS/STRAIN/ALL",
            "/ANIM/GZIP",
            "/TFILE/4",
            fmt_f(run_time / 200.0),
            "/RFILE",
            fmt_i(5000),
            "/PRINT/-100/55",
            "/KEREL",
            "/DT/NODA/CST",
            fmt_f(0.0, dt_noda),
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
        log.log(proc.stdout.decode("utf-8", errors="replace"))
        log.log(proc.stderr.decode("utf-8", errors="replace"))
        raise RuntimeError(f"anim_to_vtk failed for {anim}")
    vtk.write_bytes(proc.stdout)
    log.log(f"[exit 0] wrote {vtk}")
    return vtk


def extract_tip_displacement(vtk_path: pathlib.Path, probe_nodes: list[int]) -> tuple[float, float, float]:
    import numpy as np
    import pyvista as pv  # type: ignore[import-not-found]

    grid = pv.read(str(vtk_path))
    disp = grid.point_data.get("Displacement", grid.point_data.get("DISP"))
    if disp is None:
        raise KeyError(f"displacement field missing; available point data: {list(grid.point_data.keys())}")
    node_ids = grid.point_data.get("NODE_ID")
    samples = []
    if node_ids is not None:
        id_to_index = {int(node_id): idx for idx, node_id in enumerate(node_ids)}
        for node_id in probe_nodes:
            if node_id in id_to_index:
                samples.append(disp[id_to_index[node_id]])
    if not samples:
        raise KeyError("none of the probe nodes were found in VTK NODE_ID data")
    avg = np.asarray(samples, dtype=float).mean(axis=0)
    return float(avg[0]), float(avg[1]), float(avg[2])


def run_case(
    case: CantileverCase,
    mesh_label: str,
    alpha: float,
    solver_mode: str,
    run_time: float,
    dt_noda: float,
    damping_alpha: float,
    n_threads: int,
    log: Logger,
    env: dict[str, str],
) -> SweepRow:
    spec = MESHES[mesh_label]
    mesh = build_mesh_data(case, spec)
    alpha_label = f"{alpha:.2f}".replace(".", "p")
    job = f"stage02_{solver_mode}_{mesh_label}_a{alpha_label}"
    workdir = RUNS_DIR / job
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    msh = workdir / f"{job}.msh"
    starter = workdir / f"{job}_0000.rad"
    engine = workdir / f"{job}_0001.rad"

    total_load = case.load_for_alpha(alpha)
    log.log(
        f"[{job}] mesh={mesh_label} elements={len(mesh.bricks)} nodes={len(mesh.nodes)} "
        f"alpha={alpha:g} P={total_load:.8g} N mode={solver_mode}"
    )
    write_gmsh_mesh(mesh, msh)
    write_starter(case, mesh, alpha, starter, solver_mode, run_time, damping_alpha)
    write_engine(job, engine, solver_mode, run_time, dt_noda)

    start = time.perf_counter()
    proc = log.run([str(STARTER), "-i", str(starter), "-nt", str(n_threads)], cwd=workdir, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"starter failed for {job}")
    proc = log.run([str(ENGINE), "-i", str(engine), "-nt", str(n_threads)], cwd=workdir, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"engine failed for {job}")
    vtk = convert_anim_to_vtk(job, workdir, log, env)
    wall = time.perf_counter() - start

    ux, uy, _uz = extract_tip_displacement(vtk, mesh.probe_nodes)
    ref = elastica_reference(alpha)
    dx = max(0.0, -ux) / case.length
    dy = abs(uy) / case.length
    err_x = abs(dx - ref.delta_x_over_L) / max(ref.delta_x_over_L, 1.0e-12)
    err_y = abs(dy - ref.delta_y_over_L) / max(ref.delta_y_over_L, 1.0e-12)
    gated = mesh_label == "baseline" and any(math.isclose(alpha, gate) for gate in PASS_ALPHAS)
    passed = err_x <= PASS_TOL and err_y <= PASS_TOL
    log.log(
        f"[{job}] dx/L={dx:.8e} ref={ref.delta_x_over_L:.8e} err={100*err_x:.3f}% | "
        f"dy/L={dy:.8e} ref={ref.delta_y_over_L:.8e} err={100*err_y:.3f}% "
        f"verdict={'PASS' if passed else 'FAIL'}"
    )
    return SweepRow(
        alpha=alpha,
        mesh=mesh_label,
        solver_mode=solver_mode,
        elements=len(mesh.bricks),
        nodes=len(mesh.nodes),
        dx_fem_over_L=dx,
        dy_fem_over_L=dy,
        dx_ref_over_L=ref.delta_x_over_L,
        dy_ref_over_L=ref.delta_y_over_L,
        err_x=err_x,
        err_y=err_y,
        passed=passed,
        gated=gated,
        wall_clock_s=wall,
    )


def run_sweep(
    case: CantileverCase,
    mesh_labels: tuple[str, ...],
    alphas: tuple[float, ...],
    solver_mode: str,
    run_time: float,
    dt_noda: float,
    damping_alpha: float,
    n_threads: int,
    log: Logger,
    env: dict[str, str],
) -> list[SweepRow]:
    rows: list[SweepRow] = []
    for mesh_label in mesh_labels:
        for alpha in alphas:
            rows.append(
                run_case(
                    case,
                    mesh_label,
                    alpha,
                    solver_mode,
                    run_time,
                    dt_noda,
                    damping_alpha,
                    n_threads,
                    log,
                    env,
                )
            )
    return rows


def write_timeseries(rows: list[SweepRow]) -> pathlib.Path:
    out = RESULTS_DIR / "timeseries.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "stage",
                "alpha",
                "mesh",
                "solver_mode",
                "elements",
                "nodes",
                "dx_fem_over_L",
                "dy_fem_over_L",
                "dx_ref_over_L",
                "dy_ref_over_L",
                "err_x",
                "err_y",
                "gated",
                "verdict",
                "wall_clock_s",
            ]
        )
        for r in rows:
            writer.writerow(
                [
                    STAGE,
                    r.alpha,
                    r.mesh,
                    r.solver_mode,
                    r.elements,
                    r.nodes,
                    r.dx_fem_over_L,
                    r.dy_fem_over_L,
                    r.dx_ref_over_L,
                    r.dy_ref_over_L,
                    r.err_x,
                    r.err_y,
                    int(r.gated),
                    "PASS" if r.passed else "FAIL",
                    r.wall_clock_s,
                ]
            )
    return out


def write_results_json(rows: list[SweepRow], solver_mode: str) -> pathlib.Path:
    gated_rows = [r for r in rows if r.gated]
    have_all_pass_alphas = all(any(math.isclose(r.alpha, alpha) for r in gated_rows) for alpha in PASS_ALPHAS)
    if not have_all_pass_alphas:
        verdict = "INCONCLUSIVE"
    else:
        verdict = "PASS" if all(r.passed for r in gated_rows) else "FAIL"
    metrics = {
        f"{r.mesh}_alpha_{r.alpha:g}": {
            "solver_mode": r.solver_mode,
            "dx_fem_over_L": r.dx_fem_over_L,
            "dy_fem_over_L": r.dy_fem_over_L,
            "dx_ref_over_L": r.dx_ref_over_L,
            "dy_ref_over_L": r.dy_ref_over_L,
            "err_x": r.err_x,
            "err_y": r.err_y,
            "elements": r.elements,
            "nodes": r.nodes,
            "wall_clock_s": r.wall_clock_s,
            "gated": r.gated,
        }
        for r in rows
    }
    payload = {
        "stage": STAGE,
        "name": STAGE_NAME,
        "verdict": verdict,
        "solver_mode": solver_mode,
        "toolchain_note": (
            "Explicit dynamic relaxation is a user-selected exploratory fallback; gated "
            "Stage 02 verification expects nonlinear implicit completion."
            if solver_mode == "explicit"
            else "MUMPS-linked nonlinear implicit path completed the recorded run."
        ),
        "metrics": metrics,
        "reference": {
            "type": "Bisshopp-Drucker elastica via elliptic integrals",
            "pass_alphas": list(PASS_ALPHAS),
        },
        "tolerance": {"relative_error_dx": PASS_TOL, "relative_error_dy": PASS_TOL, "gated_mesh": "baseline"},
        "git_sha": git_sha(),
    }
    out = RESULTS_DIR / "results.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return out


def write_typst_figure(rows: list[SweepRow]) -> pathlib.Path:
    out = FIGURES_DIR / "stage02_elastica.typ"
    baseline = [r for r in rows if r.mesh == "baseline"]
    if not baseline:
        baseline = rows
    max_alpha = max(r.alpha for r in baseline)
    ymax = max(max(r.dy_ref_over_L, r.dy_fem_over_L) for r in baseline) * 1.10
    x0, y0, width, height = 1.15, 0.85, 13.4, 6.0

    def xmap(alpha: float) -> float:
        return x0 + (alpha / max_alpha) * width

    def ymap(value: float) -> float:
        return y0 + (value / ymax) * height

    ref_points = " ".join(f"({xmap(r.alpha):.3f}, {ymap(r.dy_ref_over_L):.3f})" for r in baseline)
    fem_points = " ".join(f"({xmap(r.alpha):.3f}, {ymap(r.dy_fem_over_L):.3f})" for r in baseline)
    markers = []
    for r in baseline:
        markers.extend(
            [
                f"  circle(({xmap(r.alpha):.3f}, {ymap(r.dy_fem_over_L):.3f}), radius: 0.055, fill: garnet, stroke: none)",
                f"  content(({xmap(r.alpha):.3f}, {y0 - 0.32:.3f}), [{r.alpha:g}], anchor: \"north\")",
            ]
        )
    table_rows = "\n".join(
        f"  [{r.alpha:g}], [{r.mesh}], [{100*r.err_y:.2f}\\%], [{100*r.err_x:.2f}\\%], [{'PASS' if r.passed else 'FAIL'}],"
        for r in rows
    )
    lines = [
        '#import "@preview/cetz:0.3.4"',
        "",
        '#set page(width: 180mm, height: 118mm, margin: 10mm)',
        '#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))',
        '#let garnet = rgb("#73000A")',
        '#let charcoal = rgb("#363636")',
        '#let black70 = rgb("#5C5C5C")',
        '#let black10 = rgb("#ECECEC")',
        '#let atlantic = rgb("#466A9F")',
        '#let horseshoe = rgb("#65780B")',
        '#let white = rgb("#FFFFFF")',
        "",
        "#align(center)[#text(size: 11pt, weight: \"bold\")[Stage 02 Cantilever Elastica Check]]",
        "#v(2mm)",
        "#cetz.canvas(length: 1cm, {",
        "  import cetz.draw: *",
        f"  rect(({x0:.3f}, {y0:.3f}), ({x0 + width:.3f}, {y0 + height:.3f}), fill: white, stroke: charcoal + 0.65pt)",
    ]
    for tick in (0.0, 0.2, 0.4, 0.6, 0.8):
        if tick <= ymax:
            y = ymap(tick)
            lines.append(f"  line(({x0:.3f}, {y:.3f}), ({x0 + width:.3f}, {y:.3f}), stroke: black10 + 0.45pt)")
            lines.append(f"  content(({x0 - 0.16:.3f}, {y:.3f}), [{tick:.1f}], anchor: \"east\")")
    if len(baseline) >= 2:
        lines.extend(
            [
                f"  line({ref_points}, stroke: atlantic + 0.95pt)",
                f"  line({fem_points}, stroke: garnet + 1.15pt)",
            ]
        )
    else:
        only = baseline[0]
        lines.extend(
            [
                f"  circle(({xmap(only.alpha):.3f}, {ymap(only.dy_ref_over_L):.3f}), radius: 0.055, fill: atlantic, stroke: none)",
            ]
        )
    lines.extend(
        [
            *markers,
            f"  content(({x0 + width / 2:.3f}, {y0 - 0.78:.3f}), [Dimensionless load alpha], anchor: \"north\")",
            f"  content(({x0 - 0.88:.3f}, {y0 + height / 2:.3f}), [Tip sag, dy/L], angle: 90deg)",
            f"  line(({x0 + 8.9:.3f}, {y0 + height + 0.40:.3f}), ({x0 + 9.65:.3f}, {y0 + height + 0.40:.3f}), stroke: atlantic + 0.95pt)",
            f"  content(({x0 + 9.75:.3f}, {y0 + height + 0.40:.3f}), [Elastica], anchor: \"west\")",
            f"  line(({x0 + 11.2:.3f}, {y0 + height + 0.40:.3f}), ({x0 + 11.95:.3f}, {y0 + height + 0.40:.3f}), stroke: garnet + 1.15pt)",
            f"  content(({x0 + 12.05:.3f}, {y0 + height + 0.40:.3f}), [OpenRadioss], anchor: \"west\")",
            "})",
            "",
            '#text(size: 8pt, fill: black70)[Source CSV: #raw("../results/timeseries.csv")]',
            "",
            "#figure(",
            "  table(",
            "    columns: 5,",
            "    [alpha], [mesh], [dy error], [dx error], [verdict],",
            table_rows,
            "  ),",
            "  caption: [Stage 02 tip displacement error by load and mesh.]",
            ")",
            "",
        ]
    )
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["auto", "implicit", "explicit"], default="implicit")
    parser.add_argument("--meshes", default="baseline", help="comma-separated subset of coarse,baseline,fine")
    parser.add_argument("--alphas", default="1,3,5", help="comma-separated alpha values")
    parser.add_argument("--run-time", type=float, default=1.0)
    parser.add_argument("--dt-noda", type=float, default=1.0e-4)
    parser.add_argument("--damping-alpha", type=float, default=40.0)
    parser.add_argument("--n-threads", type=int, default=1)
    parser.add_argument("--ref-only", action="store_true")
    args = parser.parse_args(argv)

    if args.ref_only:
        print(f"{'alpha':>6} {'theta_L':>10} {'dy/L':>10} {'dx/L':>10}")
        for alpha in PLOT_ALPHAS:
            ref = elastica_reference(alpha)
            print(f"{ref.alpha:>6.2f} {ref.theta_L:>10.4f} {ref.delta_y_over_L:>10.4f} {ref.delta_x_over_L:>10.4f}")
        return 0

    mesh_labels = parse_csv_list(args.meshes, str)
    alphas = parse_csv_list(args.alphas, float)
    for mesh_label in mesh_labels:
        if mesh_label not in MESHES:
            raise SystemExit(f"unknown mesh {mesh_label}")

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    log = Logger(RUN_LOG)
    case = CantileverCase()
    env = radioss_env()
    solver_mode = "implicit" if args.mode == "auto" else args.mode
    try:
        log.log(f"Stage 02 {STAGE_NAME}")
        log.log(f"OpenRadioss root: {OR_DIR}")
        log.log(f"EI={case.ei:.8g} N m^2, inertia_z={case.inertia_z:.8e} m^4")
        rows: list[SweepRow]
        rows = run_sweep(
            case,
            mesh_labels,
            alphas,
            solver_mode,
            args.run_time,
            args.dt_noda,
            args.damping_alpha,
            args.n_threads,
            log,
            env,
        )

        write_timeseries(rows)
        results_path = write_results_json(rows, solver_mode)
        write_typst_figure(rows)
        payload = json.loads(results_path.read_text(encoding="utf-8"))
        return 0 if payload["verdict"] == "PASS" else 1
    finally:
        log.close()


if __name__ == "__main__":
    sys.exit(main())
