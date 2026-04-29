"""Stage 05 runner: notched dogbone with isotropic ductile damage.

This runner builds structured all-HEXA8 notched dogbones, writes native
OpenRadioss decks using /MAT/LAW22, solves three mesh densities with explicit
dynamic relaxation, extracts load-displacement histories, and evaluates the
pre-onset mesh-objectivity RMSE. Post-peak divergence is reported, not gated,
because LAW22 is a local damage model without nonlocal regularization.
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
import re
import shutil
import subprocess
import sys
import time
from typing import Iterable


STAGE = 5
STAGE_NAME = "dogbone_damage"
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
TH_TO_CSV = OR_DIR / "exec" / "th_to_csv_linux64_gf"


@dataclasses.dataclass(frozen=True)
class Material:
    name: str = "A36_LAW22"
    rho: float = 7850.0
    young: float = 2.0e11
    nu: float = 0.30
    sigma_y: float = 2.50e8
    hard_b: float = 2.75e8
    hard_n: float = 0.36
    eps_max: float = 0.50
    sigma_max: float = 4.50e8
    eps_damage: float = 0.05
    damage_softening_slope: float = -5.0e9


@dataclasses.dataclass(frozen=True)
class Geometry:
    length: float = 0.200
    gauge_length: float = 0.060
    gauge_width: float = 0.0125
    net_width: float = 0.0080
    notch_radius: float = 0.00225
    grip_width: float = 0.025
    fillet_length: float = 0.0125
    thickness: float = 0.006
    u_end: float = 0.006

    @property
    def grip_length(self) -> float:
        return 0.5 * (self.length - self.gauge_length - 2.0 * self.fillet_length)

    @property
    def net_area(self) -> float:
        return self.net_width * self.thickness

    def width_at(self, x: float) -> float:
        """Dogbone width at centered axial coordinate x."""
        half_gauge = 0.5 * self.gauge_length
        ax = abs(x)
        if ax <= self.notch_radius:
            cut = math.sqrt(max(self.notch_radius**2 - ax**2, 0.0))
            return self.gauge_width - 2.0 * cut
        if ax <= half_gauge:
            return self.gauge_width
        if ax >= half_gauge + self.fillet_length:
            return self.grip_width
        s = (ax - half_gauge) / self.fillet_length
        return self.gauge_width + 0.5 * (1.0 - math.cos(math.pi * s)) * (
            self.grip_width - self.gauge_width
        )


@dataclasses.dataclass(frozen=True)
class MeshSpec:
    label: str
    n_grip: int
    n_fillet: int
    n_gauge_side: int
    n_notch: int
    ny: int
    nz: int
    h_ligament_mm: float


MESHES = {
    "coarse": MeshSpec("coarse", 16, 6, 36, 16, 12, 6, 1.00),
    "medium": MeshSpec("medium", 24, 8, 60, 28, 16, 8, 0.50),
    "fine": MeshSpec("fine", 32, 10, 80, 36, 20, 10, 0.25),
}


@dataclasses.dataclass
class MeshData:
    spec: MeshSpec
    nodes: dict[int, tuple[float, float, float]]
    bricks: list[tuple[int, tuple[int, int, int, int, int, int, int, int]]]
    node_sets: dict[str, list[int]]
    ligament_elements: list[int]
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


def brick_group_block(group_id: int, name: str, elems: Iterable[int]) -> list[str]:
    ids = list(elems)
    lines = [f"/GRBRIC/BRIC/{group_id}", name]
    for i in range(0, len(ids), 10):
        lines.append(fmt_i(*ids[i : i + 10]))
    return lines


def axis_points(geom: Geometry, spec: MeshSpec) -> list[float]:
    half = 0.5 * geom.length
    half_gauge = 0.5 * geom.gauge_length
    rn = geom.notch_radius
    segments = [
        (-half, -half_gauge - geom.fillet_length, spec.n_grip),
        (-half_gauge - geom.fillet_length, -half_gauge, spec.n_fillet),
        (-half_gauge, -rn, spec.n_gauge_side),
        (-rn, rn, spec.n_notch),
        (rn, half_gauge, spec.n_gauge_side),
        (half_gauge, half_gauge + geom.fillet_length, spec.n_fillet),
        (half_gauge + geom.fillet_length, half, spec.n_grip),
    ]
    xs: list[float] = []
    for start, stop, count in segments:
        xs.extend(start + (stop - start) * i / count for i in range(count))
    xs.append(half)
    return xs


def build_mesh_data(geom: Geometry, spec: MeshSpec) -> MeshData:
    xs = axis_points(geom, spec)
    zs = [-0.5 * geom.thickness + geom.thickness * k / spec.nz for k in range(spec.nz + 1)]
    nodes: dict[int, tuple[float, float, float]] = {}

    def nid(i: int, j: int, k: int) -> int:
        return 1 + k * (spec.ny + 1) * len(xs) + j * len(xs) + i

    for k, z in enumerate(zs):
        for j in range(spec.ny + 1):
            eta = -1.0 + 2.0 * j / spec.ny
            for i, x in enumerate(xs):
                y = 0.5 * geom.width_at(x) * eta
                nodes[nid(i, j, k)] = (x, y, z)

    bricks: list[tuple[int, tuple[int, int, int, int, int, int, int, int]]] = []
    ligament_elements: list[int] = []
    eid = 1
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
                if abs(cx) <= 1.5 * geom.notch_radius:
                    ligament_elements.append(eid)
                eid += 1

    i_left = 0
    i_right = len(xs) - 1
    i_gauge_neg = min(range(len(xs)), key=lambda idx: abs(xs[idx] + 0.5 * geom.gauge_length))
    i_gauge_pos = min(range(len(xs)), key=lambda idx: abs(xs[idx] - 0.5 * geom.gauge_length))
    j_mid = spec.ny // 2
    k_mid = spec.nz // 2

    left = [nid(i_left, j, k) for k in range(spec.nz + 1) for j in range(spec.ny + 1)]
    right = [nid(i_right, j, k) for k in range(spec.nz + 1) for j in range(spec.ny + 1)]
    right_anchor = [nid(i_right, j_mid, k_mid)]
    right_probe = [nid(i_right, j_mid, k_mid)]
    all_nodes = list(nodes)
    return MeshData(
        spec=spec,
        nodes=nodes,
        bricks=bricks,
        node_sets={
            "left": left,
            "right": right,
            "right_anchor": right_anchor,
            "right_probe": right_probe,
            "all_nodes": all_nodes,
        },
        ligament_elements=ligament_elements,
        gauge_x_neg=xs[i_gauge_neg],
        gauge_x_pos=xs[i_gauge_pos],
    )


def write_gmsh_mesh(mesh: MeshData, out_msh: pathlib.Path) -> None:
    import gmsh  # type: ignore[import-not-found]

    gmsh.initialize(["-nopopup"])
    try:
        gmsh.model.add(f"stage05_{mesh.spec.label}")
        gmsh.model.addDiscreteEntity(3, 1)
        node_tags = list(mesh.nodes)
        coords: list[float] = []
        for node_id in node_tags:
            coords.extend(mesh.nodes[node_id])
        elem_tags = [eid for eid, _ in mesh.bricks]
        elem_conn = [node for _, conn in mesh.bricks for node in conn]
        gmsh.model.mesh.addNodes(3, 1, node_tags, coords)
        gmsh.model.mesh.addElementsByType(1, 5, elem_tags, elem_conn)
        gmsh.model.addPhysicalGroup(3, [1], tag=1, name="DOGBONE_DAMAGE_SOLID")
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(str(out_msh))
    finally:
        gmsh.finalize()


def write_starter(
    geom: Geometry,
    mat: Material,
    mesh: MeshData,
    out_rad: pathlib.Path,
    run_time: float,
    damping_alpha: float,
) -> None:
    job = out_rad.name.removesuffix("_0000.rad")
    lines: list[str] = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        fmt_i(2019, 0),
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 05 LAW22 notched dogbone {mesh.spec.label}",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 4) + f"{0:20d}" + f"{2:40d}",
        "/RANDOM",
        fmt_f(0.0) + f"{0:20d}",
        "/SPMD",
        fmt_i(0, 0) + f"{0:20d}{1:20d}",
        "/SHFRA/V4",
        "/DAMP/1",
        "global_mass_damping",
        "#              alpha                beta  grnod_ID   skew_ID              Tstart               Tstop",
        fmt_f(damping_alpha, 0.0) + fmt_i(300, 0) + fmt_f(0.0, run_time),
        "/MAT/LAW22/1",
        mat.name,
        "#        Init. dens.          Ref. dens.",
        fmt_f(mat.rho, 0.0),
        "#                  E                  Nu",
        fmt_f(mat.young, mat.nu),
        "#                  a                   b                   n             Eps_max           SIGMA_max",
        fmt_f(mat.sigma_y, mat.hard_b, mat.hard_n, mat.eps_max, mat.sigma_max),
        "#                  c           Eps_dot_0       ICC",
        fmt_f(0.0, 0.0) + fmt_i(0),
        "#            Eps_dam                 E_t",
        fmt_f(mat.eps_damage, mat.damage_softening_slope),
        "/NODE",
    ]
    for nid, (x, y, z) in mesh.nodes.items():
        lines.append(f"{nid:10d}{x:20.12g}{y:20.12g}{z:20.12g}")

    lines.extend(["/PART/1", "dogbone_damage_part", fmt_i(1, 1, 0), "/BRICK/1"])
    for eid, conn in mesh.bricks:
        lines.append(fmt_i(eid, *conn))

    lines.extend(
        [
            "/PROP/SOLID/1",
            "heph_large_strain",
            "#   Isolid    Ismstr               Icpre               Inpts    Itetra    Iframe                  dn",
            fmt_i(24, 4) + f"{1:20d}{0:20d}{0:10d}{0:10d}{0:20d}",
            "#                q_a                 q_b                   h            LAMBDA_V                MU_V",
            fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
            "#             dt_min   istrain      IHKT",
            fmt_f(0.0) + fmt_i(0, 0),
            "/BCS/1",
            "left_encastre",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   111 000{0:10d}{100:10d}",
            "/BCS/2",
            "right_anchor_yz",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   011 000{0:10d}{102:10d}",
            "/FUNCT/1",
            "right_displacement_ramp",
            "#                  X                   Y",
            fmt_f(0.0, 0.0),
            fmt_f(0.25, 0.15625),
            fmt_f(0.50, 0.50),
            fmt_f(0.75, 0.84375),
            fmt_f(1.0, 1.0),
            "/IMPDISP/1",
            "right_grip_pull_x",
            "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
            f"{1:10d}{'X':>10}{0:10d}{0:10d}{101:10d}{0:10d}{0:10d}",
            "#            Scale_x             Scale_y              Tstart               Tstop",
            fmt_f(run_time, geom.u_end, 0.0, 0.0),
        ]
    )
    lines.extend(node_group_block(100, "left_grip", mesh.node_sets["left"]))
    lines.extend(node_group_block(101, "right_grip", mesh.node_sets["right"]))
    lines.extend(node_group_block(102, "right_anchor", mesh.node_sets["right_anchor"]))
    lines.extend(node_group_block(103, "right_probe", mesh.node_sets["right_probe"]))
    lines.extend(node_group_block(300, "all_nodes", mesh.node_sets["all_nodes"]))
    lines.extend(brick_group_block(200, "notch_ligament", mesh.ligament_elements))
    lines.extend(
        [
            "/TH/NODE/1",
            "right_probe",
            "DEF",
        ]
    )
    for node_id in mesh.node_sets["right_probe"]:
        lines.append(f"{node_id:10d}{0:10d}right_probe_{node_id}")
    lines.extend(["/TH/NODE/2", "right_reaction_x", "REACX"])
    for node_id in mesh.node_sets["right"]:
        lines.append(f"{node_id:10d}{0:10d}right_{node_id}")
    lines.extend(["/TH/NODE/3", "left_reaction_x", "REACX"])
    for node_id in mesh.node_sets["left"]:
        lines.append(f"{node_id:10d}{0:10d}left_{node_id}")
    lines.extend(["/END", ""])
    out_rad.write_text("\n".join(lines), encoding="utf-8")


def write_engine(job: str, out_rad: pathlib.Path, run_time: float, dt_noda: float) -> None:
    anim_dt = run_time / 60.0
    lines = [
        "#RADIOSS ENGINE",
        "/ANIM/DT",
        fmt_f(anim_dt, anim_dt),
        "/ANIM/VECT/DISP",
        "/ANIM/VECT/VEL",
        "/ANIM/VECT/FREAC",
        "/ANIM/VECT/FINT",
        "/ANIM/BRICK/TENS/STRESS/ALL",
        "/ANIM/BRICK/TENS/STRAIN/ALL",
        "/ANIM/BRICK/TENS/DAMA",
        "/ANIM/ELEM/EPSP",
        "/ANIM/GZIP",
        "/TFILE/4",
        fmt_f(run_time / 240.0),
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


def _anim_frame_indices(job: str, workdir: pathlib.Path) -> list[int]:
    pattern = re.compile(re.escape(job) + r"A(\d{3})(?:\.gz)?$")
    indices: set[int] = set()
    for path in workdir.glob(f"{job}A*"):
        match = pattern.match(path.name)
        if match:
            indices.add(int(match.group(1)))
    return sorted(indices)


def convert_anim_frame_to_vtk(
    job: str,
    workdir: pathlib.Path,
    frame_index: int,
    log: Logger,
    env: dict[str, str],
) -> pathlib.Path:
    anim = workdir / f"{job}A{frame_index:03d}"
    gz = anim.with_suffix(anim.suffix + ".gz")
    if gz.exists() and (not anim.exists() or gz.stat().st_mtime > anim.stat().st_mtime):
        gunzip_keep(gz, anim)
    if not anim.exists():
        raise FileNotFoundError(f"animation frame not found: {anim} or {gz}")
    vtk = workdir / f"{job}A{frame_index:03d}.vtk"
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
    log.log(f"[exit {proc.returncode}] wrote {vtk}")
    return vtk


def convert_all_anim_to_vtk(job: str, workdir: pathlib.Path, log: Logger, env: dict[str, str]) -> list[pathlib.Path]:
    indices = _anim_frame_indices(job, workdir)
    if not indices:
        indices = [1]
    log.log(f"[{job}] converting {len(indices)} animation frame(s) to VTK")
    return [convert_anim_frame_to_vtk(job, workdir, idx, log, env) for idx in indices]


def convert_t01_to_csv(job: str, workdir: pathlib.Path, log: Logger, env: dict[str, str]) -> pathlib.Path:
    t01 = workdir / f"{job}T01"
    gz = t01.with_suffix(t01.suffix + ".gz")
    if gz.exists() and (not t01.exists() or gz.stat().st_mtime > t01.stat().st_mtime):
        gunzip_keep(gz, t01)
    if not t01.exists():
        candidates = sorted(workdir.glob("*T01"))
        if not candidates:
            raise FileNotFoundError(f"T01 file not found for {job}")
        t01 = candidates[0]
    log.log("$ " + " ".join([str(TH_TO_CSV), str(t01)]))
    proc = subprocess.run(
        [str(TH_TO_CSV), str(t01)],
        cwd=str(workdir),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if proc.stdout:
        log.log(proc.stdout.rstrip())
    if proc.returncode != 0:
        raise RuntimeError(f"th_to_csv failed for {t01}")
    csv_path = pathlib.Path(str(t01) + ".csv")
    if not csv_path.exists():
        raise FileNotFoundError(f"th_to_csv did not write {csv_path}")
    return csv_path


def read_force_curve(csv_path: pathlib.Path, geom: Geometry, run_time: float) -> dict[str, object]:
    with csv_path.open(newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        columns = reader.fieldnames or []
        reaction_cols = [c for c in columns if "right_reaction_x" in c]
        if reaction_cols:
            force_source = "right_reaction_x"
        else:
            reaction_cols = [c for c in columns if "left_reaction_x" in c]
            force_source = "left_reaction_x"
        if not reaction_cols:
            raise KeyError("reaction columns not found in T01 CSV")
        time_s: list[float] = []
        disp_m: list[float] = []
        force_n: list[float] = []
        for row in reader:
            t = float(row["time"])
            reaction = sum(float(row[c]) for c in reaction_cols if row.get(c, ""))
            time_s.append(t)
            disp_m.append(geom.u_end * min(max(t / run_time, 0.0), 1.0))
            force_n.append(abs(reaction))
    if len(time_s) < 3:
        raise RuntimeError(f"not enough T01 rows in {csv_path}")
    return {"time_s": time_s, "displacement_m": disp_m, "force_N": force_n, "force_source": force_source}


def read_anim_force_curve(
    vtk_paths: list[pathlib.Path],
    mesh: MeshData,
    geom: Geometry,
    run_time: float,
) -> dict[str, object]:
    import numpy as np
    import pyvista as pv  # type: ignore[import-not-found]

    time_s = [0.0]
    disp_m = [0.0]
    force_n = [0.0]
    force_source = ""
    right_nodes = set(mesh.node_sets["right"])
    for frame_i, vtk_path in enumerate(vtk_paths, start=1):
        grid = pv.read(str(vtk_path))
        node_ids = grid.point_data.get("NODE_ID")
        disp = grid.point_data.get("Displacement")
        if node_ids is None or disp is None:
            raise KeyError(f"NODE_ID or Displacement missing from {vtk_path}")
        right_idx = [i for i, node_id in enumerate(node_ids) if int(node_id) in right_nodes]
        if not right_idx:
            raise KeyError(f"right-grip nodes missing from {vtk_path}")

        vector_name = None
        for token in ("reaction", "freac"):
            for name, array in grid.point_data.items():
                if token in name.lower() and getattr(array, "ndim", 0) == 2 and array.shape[1] >= 3:
                    vector_name = name
                    break
            if vector_name:
                break
        if vector_name is None:
            for token in ("internal", "fint", "force"):
                for name, array in grid.point_data.items():
                    lname = name.lower()
                    if token in lname and getattr(array, "ndim", 0) == 2 and array.shape[1] >= 3:
                        vector_name = name
                        break
                if vector_name:
                    break
        if vector_name is None:
            raise KeyError(f"no nodal force vector in {vtk_path}; point data={list(grid.point_data.keys())}")

        force_vec = np.asarray(grid.point_data[vector_name], dtype=float)
        u_vec = np.asarray(disp, dtype=float)
        force_source = vector_name
        time_s.append(run_time * frame_i / max(len(vtk_paths), 1))
        disp_m.append(float(np.mean(u_vec[right_idx, 0])))
        force_n.append(abs(float(np.sum(force_vec[right_idx, 0]))))

    order = np.argsort(np.asarray(disp_m))
    return {
        "time_s": [float(time_s[i]) for i in order],
        "displacement_m": [float(disp_m[i]) for i in order],
        "force_N": [float(force_n[i]) for i in order],
        "force_source": force_source,
    }


def _cell_array(grid, contains: str):
    preferred = [name for name in grid.cell_data.keys() if contains in name and "Intg" in name]
    if preferred:
        return grid.cell_data[preferred[0]]
    for name in grid.cell_data.keys():
        if contains in name:
            return grid.cell_data[name]
    raise KeyError(f"no cell data containing {contains}; available={list(grid.cell_data.keys())}")


def extract_final_metrics(vtk_path: pathlib.Path, mesh: MeshData) -> dict[str, float]:
    import numpy as np
    import pyvista as pv  # type: ignore[import-not-found]

    grid = pv.read(str(vtk_path))
    elem_ids = grid.cell_data.get("ELEMENT_ID")
    if elem_ids is None:
        raise KeyError("ELEMENT_ID missing from VTK")
    ligament_mask = np.isin(elem_ids, np.asarray(mesh.ligament_elements, dtype=elem_ids.dtype))
    stress = _cell_array(grid, "Strs")
    epsp = None
    for name in grid.cell_data.keys():
        if "Plastic" in name or "EPSP" in name or "Epsp" in name:
            epsp = grid.cell_data[name]
            break
    sigma_x = np.asarray(stress[ligament_mask, 0], dtype=float)
    payload = {
        "ligament_sigma_x_mean_Pa": float(np.mean(sigma_x)),
        "ligament_sigma_x_max_Pa": float(np.max(sigma_x)),
    }
    if epsp is not None:
        epsp_lig = np.asarray(epsp[ligament_mask], dtype=float)
        payload["ligament_epsp_max"] = float(np.max(epsp_lig))
        payload["ligament_epsp_mean"] = float(np.mean(epsp_lig))
    else:
        payload["ligament_epsp_max"] = float("nan")
        payload["ligament_epsp_mean"] = float("nan")
    return payload


def run_case(
    geom: Geometry,
    mat: Material,
    spec: MeshSpec,
    run_time: float,
    damping_alpha: float,
    dt_noda: float,
    n_threads: int,
    log: Logger,
    env: dict[str, str],
) -> dict[str, object]:
    mesh = build_mesh_data(geom, spec)
    job = f"stage05_{spec.label}"
    workdir = RUNS_DIR / job
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    starter = workdir / f"{job}_0000.rad"
    engine = workdir / f"{job}_0001.rad"
    msh = workdir / f"{job}.msh"
    log.log(f"[{job}] elements={len(mesh.bricks)} nodes={len(mesh.nodes)} h={spec.h_ligament_mm:.2f} mm")
    write_gmsh_mesh(mesh, msh)
    write_starter(geom, mat, mesh, starter, run_time, damping_alpha)
    write_engine(job, engine, run_time, dt_noda)

    start = time.perf_counter()
    proc = log.run([str(STARTER), "-i", str(starter), "-nt", str(n_threads)], cwd=workdir, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"starter failed for {job}")
    proc = log.run([str(ENGINE), "-i", str(engine), "-nt", str(n_threads)], cwd=workdir, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"engine failed for {job}")
    vtk_paths = convert_all_anim_to_vtk(job, workdir, log, env)
    th_csv = convert_t01_to_csv(job, workdir, log, env)
    try:
        curve = read_anim_force_curve(vtk_paths, mesh, geom, run_time)
    except Exception as exc:
        log.log(f"[{job}] animation force extraction failed ({exc}); falling back to T01 reactions")
        curve = read_force_curve(th_csv, geom, run_time)
    final = extract_final_metrics(vtk_paths[-1], mesh)
    force = curve["force_N"]
    disp = curve["displacement_m"]
    peak_i = max(range(len(force)), key=lambda i: force[i])
    result: dict[str, object] = {
        "label": spec.label,
        "h_ligament_mm": spec.h_ligament_mm,
        "elements": len(mesh.bricks),
        "nodes": len(mesh.nodes),
        "wall_clock_s": time.perf_counter() - start,
        "curve": curve,
        "force_source": curve["force_source"],
        "peak_force_N": force[peak_i],
        "peak_displacement_m": disp[peak_i],
        "final_force_N": force[-1],
        "final_displacement_m": disp[-1],
        "vtk_path": str(vtk_paths[-1]),
        "vtk_first_frame_path": str(vtk_paths[0]),
        "vtk_frame_count": len(vtk_paths),
        **final,
    }
    log.log(
        f"[{job}] peak={force[peak_i]:.6g} N at u={disp[peak_i]:.6g} m "
        f"final_epsp={result['ligament_epsp_max']:.6g}"
    )
    return result


def interp_force(row: dict[str, object], u_grid):
    import numpy as np

    curve = row["curve"]
    assert isinstance(curve, dict)
    return np.interp(u_grid, curve["displacement_m"], curve["force_N"])


def windowed_rmse(a: dict[str, object], b: dict[str, object], u_hi: float, f_norm: float) -> float:
    import numpy as np

    if u_hi <= 0.0:
        return float("nan")
    u_grid = np.linspace(0.0, u_hi, 256)
    fa = interp_force(a, u_grid)
    fb = interp_force(b, u_grid)
    return float(np.sqrt(np.mean(((fa - fb) / max(abs(f_norm), 1.0)) ** 2)))


def post_rmse(a: dict[str, object], b: dict[str, object], u_lo: float, u_hi: float, f_norm: float) -> float:
    import numpy as np

    if u_hi <= u_lo:
        return float("nan")
    u_grid = np.linspace(u_lo, u_hi, 256)
    fa = interp_force(a, u_grid)
    fb = interp_force(b, u_grid)
    return float(np.sqrt(np.mean(((fa - fb) / max(abs(f_norm), 1.0)) ** 2)))


def evaluate(rows: list[dict[str, object]], mat: Material, geom: Geometry) -> tuple[str, dict[str, object]]:
    by_label = {str(row["label"]): row for row in rows}
    medium = by_label["medium"]
    # The damage tensor is not consistently exposed in converted VTK on this
    # build. Use the LAW22 plastic-strain damage threshold as the pre-onset
    # window proxy, mapped through the 60 mm gauge length.
    u_d = min(geom.gauge_length * (mat.sigma_y / mat.young + mat.eps_damage), geom.u_end)
    f_peak = float(medium["peak_force_N"])
    rmse_cm = windowed_rmse(by_label["coarse"], by_label["medium"], u_d, f_peak)
    rmse_mf = windowed_rmse(by_label["medium"], by_label["fine"], u_d, f_peak)
    rmse_post_cm = post_rmse(by_label["coarse"], by_label["medium"], u_d, geom.u_end, f_peak)
    rmse_post_mf = post_rmse(by_label["medium"], by_label["fine"], u_d, geom.u_end, f_peak)
    final_epsp = float(medium.get("ligament_epsp_max", float("nan")))
    damage_reached = math.isfinite(final_epsp) and final_epsp >= mat.eps_damage
    pre_pass = rmse_cm <= 0.05 and rmse_mf <= 0.05
    verdict = "PASS" if pre_pass and damage_reached else "FAIL"
    checks = {
        "pre_onset_rmse_coarse_medium": {
            "pass": rmse_cm <= 0.05,
            "value": rmse_cm,
            "tolerance": 0.05,
            "gating": True,
        },
        "pre_onset_rmse_medium_fine": {
            "pass": rmse_mf <= 0.05,
            "value": rmse_mf,
            "tolerance": 0.05,
            "gating": True,
        },
        "damage_region_reached": {
            "pass": damage_reached,
            "value_epsp_max_medium": final_epsp,
            "eps_damage": mat.eps_damage,
            "gating": True,
        },
        "post_onset_rmse_coarse_medium": {
            "pass": rmse_post_cm <= 0.05,
            "value": rmse_post_cm,
            "tolerance": 0.05,
            "gating": False,
            "note": "reported only; LAW22 local damage softening is expected to be mesh dependent after localization",
        },
        "post_onset_rmse_medium_fine": {
            "pass": rmse_post_mf <= 0.05,
            "value": rmse_post_mf,
            "tolerance": 0.05,
            "gating": False,
            "note": "reported only; LAW22 local damage softening is expected to be mesh dependent after localization",
        },
    }
    return verdict, {
        "u_damage_onset_proxy_m": u_d,
        "F_peak_medium_N": f_peak,
        "checks": checks,
    }


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


def write_timeseries(rows: list[dict[str, object]]) -> pathlib.Path:
    import numpy as np

    out = RESULTS_DIR / "timeseries.csv"
    u_max = min(float(row["final_displacement_m"]) for row in rows)
    u_grid = np.linspace(0.0, u_max, 400)
    by_label = {str(row["label"]): row for row in rows}
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["stage", "displacement_m", "F_coarse_N", "F_medium_N", "F_fine_N"])
        for u, fc, fm, ff in zip(
            u_grid,
            interp_force(by_label["coarse"], u_grid),
            interp_force(by_label["medium"], u_grid),
            interp_force(by_label["fine"], u_grid),
        ):
            writer.writerow([STAGE, u, fc, fm, ff])
    mesh_csv = RESULTS_DIR / "mesh_objectivity.csv"
    shutil.copy(out, mesh_csv)
    return out


def write_results_json(
    rows: list[dict[str, object]],
    verdict: str,
    evaluation: dict[str, object],
    mat: Material,
    geom: Geometry,
) -> pathlib.Path:
    compact_rows = []
    for row in rows:
        compact = {k: v for k, v in row.items() if k != "curve"}
        compact_rows.append(compact)
    payload = {
        "stage": STAGE,
        "name": STAGE_NAME,
        "verdict": verdict,
        "metrics": {str(row["label"]): row for row in compact_rows},
        "checks": evaluation["checks"],
        "reference": {
            "type": "LAW22 mesh self-consistency through pre-onset/pre-localization window",
            "material": dataclasses.asdict(mat),
            "geometry": dataclasses.asdict(geom),
            "u_damage_onset_proxy_m": evaluation["u_damage_onset_proxy_m"],
            "F_peak_medium_N": evaluation["F_peak_medium_N"],
            "damage_onset_note": (
                "The LAW22 damage tensor is not consistently exposed by this "
                "build's anim_to_vtk output; the material plastic-strain damage "
                "threshold mapped to gauge displacement is used as the "
                "pre-localization window proxy."
            ),
        },
        "tolerance": {"pre_onset_rmse": 0.05},
        "git_sha": git_sha(),
    }
    out = RESULTS_DIR / "results.json"
    out.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    summary = RESULTS_DIR / "mesh_objectivity_summary.json"
    summary.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    return out


def write_typst_figure(rows: list[dict[str, object]], verdict: str, evaluation: dict[str, object]) -> pathlib.Path:
    out = FIGURES_DIR / "stage05_mesh_objectivity.typ"
    checks = evaluation["checks"]
    lines = [
        '#import "@preview/cetz:0.3.4"',
        "",
        '#set page(width: 180mm, height: 116mm, margin: 10mm)',
        '#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))',
        '#let garnet = rgb("#73000A")',
        '#let charcoal = rgb("#363636")',
        '#let black10 = rgb("#ECECEC")',
        '#let atlantic = rgb("#466A9F")',
        '#let white = rgb("#FFFFFF")',
        "",
        "#align(center)[#text(size: 11pt, weight: \"bold\")[Stage 05 LAW22 Mesh Objectivity]]",
        "#v(2mm)",
        "#figure(",
        "  table(",
        "    columns: 5,",
        "    [Mesh], [Elements], [Peak force (N)], [Final EPSP max], [Wall-clock (s)],",
    ]
    for row in rows:
        lines.append(
            f"    [{row['label']}], [{int(row['elements'])}], "
            f"[{float(row['peak_force_N']):.2f}], "
            f"[{float(row['ligament_epsp_max']):.4f}], "
            f"[{float(row['wall_clock_s']):.2f}],"
        )
    lines.extend(
        [
            "  ),",
            "  caption: [OpenRadioss LAW22 notched-dogbone mesh sweep.]",
            ")",
            "",
            "#figure(",
            "  table(",
            "    columns: 4,",
            "    [Check], [Value], [Tolerance], [Gating],",
            f"    [Coarse-medium pre RMSE], [{100*checks['pre_onset_rmse_coarse_medium']['value']:.3f}\\%], [5.000\\%], [yes],",
            f"    [Medium-fine pre RMSE], [{100*checks['pre_onset_rmse_medium_fine']['value']:.3f}\\%], [5.000\\%], [yes],",
            f"    [Coarse-medium post RMSE], [{100*checks['post_onset_rmse_coarse_medium']['value']:.3f}\\%], [reported], [no],",
            f"    [Medium-fine post RMSE], [{100*checks['post_onset_rmse_medium_fine']['value']:.3f}\\%], [reported], [no],",
            "  ),",
            f"  caption: [Verdict: {verdict}. Post-onset divergence is reported because LAW22 is local CDM.]",
            ")",
            "",
        ]
    )
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def write_blocker(verdict: str) -> None:
    (STAGE_DIR / "blocker.md").write_text(
        "# Stage 05 Blocker: LAW22 Mesh-Objectivity Gate Failed\n\n"
        "**Author.** J.C. Vaught\n\n"
        f"Verdict: `{verdict}`\n\n"
        "The OpenRadioss LAW22 runs completed, but the pre-onset RMSE gate or "
        "damage-region reachability check failed. See `results/results.json` "
        "and `run.log` for the measured values.\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meshes", default="coarse,medium,fine")
    parser.add_argument("--run-time", type=float, default=0.03)
    parser.add_argument("--dt-noda", type=float, default=1.0e-7)
    parser.add_argument("--damping-alpha", type=float, default=8000.0)
    parser.add_argument("--n-threads", type=int, default=16)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    geom = Geometry()
    mat = Material()
    labels = [v.strip() for v in args.meshes.split(",") if v.strip()]
    if set(labels) != set(MESHES):
        raise SystemExit("--meshes must include coarse,medium,fine for the stage verdict")

    log = Logger(RUN_LOG)
    try:
        log.log(f"Stage 05 {STAGE_NAME}")
        log.log(f"OpenRadioss root: {OR_DIR}")
        log.log(f"LAW22 material: {mat.name}, eps_damage={mat.eps_damage:g}")
        if args.dry_run:
            for label in labels:
                mesh = build_mesh_data(geom, MESHES[label])
                log.log(f"dry-run {label}: elements={len(mesh.bricks)} nodes={len(mesh.nodes)}")
            return 0

        env = radioss_env()
        rows = [
            run_case(
                geom,
                mat,
                MESHES[label],
                args.run_time,
                args.damping_alpha,
                args.dt_noda,
                args.n_threads,
                log,
                env,
            )
            for label in labels
        ]
        rows.sort(key=lambda row: ("coarse", "medium", "fine").index(str(row["label"])))
        verdict, evaluation = evaluate(rows, mat, geom)
        write_timeseries(rows)
        write_results_json(rows, verdict, evaluation, mat, geom)
        write_typst_figure(rows, verdict, evaluation)
        if verdict == "FAIL":
            write_blocker(verdict)
        else:
            blocker = STAGE_DIR / "blocker.md"
            if blocker.exists():
                blocker.unlink()
        log.log(f"Stage 05 verdict: {verdict}")
        return 0 if verdict == "PASS" else 1
    finally:
        log.close()


if __name__ == "__main__":
    sys.exit(main())
