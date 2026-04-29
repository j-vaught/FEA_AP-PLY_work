"""Stage 01 runner: linear-elastic 3-point and 4-point solid-beam bending.

The runner builds a structured HEXA8 mesh, writes OpenRadioss block decks
directly from the same mesh, runs starter/engine, converts animation frames to
legacy VTK, extracts midspan displacement, and writes the stage artifacts.
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


STAGE = 1
STAGE_NAME = "beam_bending"
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


@dataclasses.dataclass(frozen=True)
class BeamCase:
    length: float = 0.200
    width: float = 0.020
    height: float = 0.010
    young: float = 6.89e10
    nu: float = 0.33
    rho: float = 2700.0
    total_load: float = 1000.0

    @property
    def inertia_y(self) -> float:
        return self.width * self.height**3 / 12.0

    @property
    def shear_modulus(self) -> float:
        return self.young / (2.0 * (1.0 + self.nu))

    def eb_deflection(self, load_case: str) -> float:
        if load_case == "3pt":
            return self.total_load * self.length**3 / (48.0 * self.young * self.inertia_y)
        if load_case == "4pt":
            a = self.length / 3.0
            return (
                (self.total_load / 2.0)
                * a
                * (3.0 * self.length**2 - 4.0 * a**2)
                / (24.0 * self.young * self.inertia_y)
            )
        raise ValueError(load_case)

    def timoshenko_deflection(self, load_case: str) -> float:
        # Shear terms are reported, not gated. 3-point expression is exact for
        # the point-load case; 4-point uses the constant-shear support spans.
        kappa = 5.0 / 6.0
        area = self.width * self.height
        if load_case == "3pt":
            shear = self.total_load * self.length / (4.0 * kappa * self.shear_modulus * area)
        elif load_case == "4pt":
            a = self.length / 3.0
            shear = (self.total_load / 2.0) * a / (kappa * self.shear_modulus * area)
        else:
            raise ValueError(load_case)
        return self.eb_deflection(load_case) + shear


@dataclasses.dataclass(frozen=True)
class MeshSpec:
    label: str
    nx_each_sixth: int
    ny: int
    nz: int

    @property
    def nx(self) -> int:
        return 6 * self.nx_each_sixth


MESHES = {
    "M0": MeshSpec("M0", nx_each_sixth=7, ny=4, nz=5),
    "M1": MeshSpec("M1", nx_each_sixth=18, ny=12, nz=14),
    "M2": MeshSpec("M2", nx_each_sixth=28, ny=16, nz=20),
}


@dataclasses.dataclass
class MeshData:
    spec: MeshSpec
    nodes: dict[int, tuple[float, float, float]]
    bricks: list[tuple[int, tuple[int, int, int, int, int, int, int, int]]]
    node_sets: dict[str, list[int]]
    probe_node: int


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
            print(proc.stdout, end="", flush=True)
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


def grid_coordinates(case: BeamCase, spec: MeshSpec) -> tuple[list[float], list[float], list[float]]:
    # Six equal x segments make x=L/3, L/2, and 2L/3 exact node planes.
    xs = [case.length * i / spec.nx for i in range(spec.nx + 1)]
    ys = [-case.width / 2.0 + case.width * j / spec.ny for j in range(spec.ny + 1)]
    zs = [-case.height / 2.0 + case.height * k / spec.nz for k in range(spec.nz + 1)]
    return xs, ys, zs


def build_mesh_data(case: BeamCase, spec: MeshSpec) -> MeshData:
    xs, ys, zs = grid_coordinates(case, spec)
    nodes: dict[int, tuple[float, float, float]] = {}

    def nid(i: int, j: int, k: int) -> int:
        return 1 + k * (spec.ny + 1) * (spec.nx + 1) + j * (spec.nx + 1) + i

    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                nodes[nid(i, j, k)] = (x, y, z)

    bricks: list[tuple[int, tuple[int, int, int, int, int, int, int, int]]] = []
    eid = 1
    for k in range(spec.nz):
        for j in range(spec.ny):
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

    i_l3 = spec.nx // 3
    i_mid = spec.nx // 2
    i_2l3 = 2 * spec.nx // 3
    j_mid = spec.ny // 2
    k_bot = 0
    k_top = spec.nz

    support_left = [nid(0, j, k_bot) for j in range(spec.ny + 1)]
    support_right = [nid(spec.nx, j, k_bot) for j in range(spec.ny + 1)]
    anchor_left = [nid(0, j_mid, k_bot)]
    load_3pt = [nid(i_mid, j, k_top) for j in range(spec.ny + 1)]
    load_4pt_a = [nid(i_l3, j, k_top) for j in range(spec.ny + 1)]
    load_4pt_b = [nid(i_2l3, j, k_top) for j in range(spec.ny + 1)]
    all_nodes = list(nodes)
    probe = nid(i_mid, j_mid, k_bot)
    return MeshData(
        spec=spec,
        nodes=nodes,
        bricks=bricks,
        node_sets={
            "support_left": support_left,
            "support_right": support_right,
            "anchor_left": anchor_left,
            "load_3pt": load_3pt,
            "load_4pt_a": load_4pt_a,
            "load_4pt_b": load_4pt_b,
            "all_nodes": all_nodes,
        },
        probe_node=probe,
    )


def write_gmsh_mesh(case: BeamCase, mesh: MeshData, out_msh: pathlib.Path) -> None:
    import gmsh

    gmsh.initialize()
    try:
        gmsh.model.add(f"stage01_{mesh.spec.label}")
        xs, ys, zs = grid_coordinates(case, mesh.spec)
        point_tags: dict[tuple[int, int, int], int] = {}
        for k, z in enumerate(zs):
            for j, y in enumerate(ys):
                for i, x in enumerate(xs):
                    point_tags[(i, j, k)] = gmsh.model.geo.addPoint(x, y, z)
        # A real transfinite OCC/grid construction is overkill for the deck,
        # but this .msh records the same generated HEXA8 topology for audit.
        node_tags = list(mesh.nodes)
        coords: list[float] = []
        for nid in node_tags:
            coords.extend(mesh.nodes[nid])
        elem_tags = [eid for eid, _ in mesh.bricks]
        elem_conn: list[int] = []
        for _, conn in mesh.bricks:
            elem_conn.extend(conn)
        gmsh.model.addDiscreteEntity(3, 1)
        gmsh.model.mesh.addNodes(3, 1, node_tags, coords)
        gmsh.model.mesh.addElementsByType(1, 5, elem_tags, elem_conn)  # type 5 = HEXA8
        gmsh.model.addPhysicalGroup(3, [1], tag=1, name="BEAM_SOLID")
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


def write_starter(
    case: BeamCase,
    mesh: MeshData,
    load_case: str,
    out_rad: pathlib.Path,
    run_time: float,
) -> None:
    job = out_rad.name.removesuffix("_0000.rad")
    if load_case == "3pt":
        load_sets = [(200, "load_3pt", -case.total_load / len(mesh.node_sets["load_3pt"]))]
    elif load_case == "4pt":
        load_sets = [
            (210, "load_4pt_a", -0.5 * case.total_load / len(mesh.node_sets["load_4pt_a"])),
            (211, "load_4pt_b", -0.5 * case.total_load / len(mesh.node_sets["load_4pt_b"])),
        ]
    else:
        raise ValueError(load_case)

    lines: list[str] = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        fmt_i(2019, 0),
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 01 {load_case} linear elastic beam bending",
        "/ANALY",
        fmt_i(0, 0, 0),
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 0) + f"{0:20d}" + f"{2:40d}",
        "/IOFLAG",
        fmt_i(0, 0, 0, 0, 0, 0),
        "/RANDOM",
        fmt_f(0.0) + f"{0:20d}",
        "/SPMD",
        fmt_i(0, 0) + f"{0:20d}{1:20d}",
        "/SHFRA/V4",
        "/DAMP/1",
        "whole_beam_mass_damping",
        "#              alpha                beta  grnod_ID   skew_ID              Tstart               Tstop",
        fmt_f(6000.0, 0.0) + fmt_i(300, 0) + fmt_f(0.0, run_time),
        "/MAT/ELAST/1",
        "AL6061T6",
        fmt_f(case.rho, 0.0),
        fmt_f(case.young, case.nu),
        "/NODE",
    ]
    for nid, (x, y, z) in mesh.nodes.items():
        lines.append(f"{nid:10d}{x:20.12g}{y:20.12g}{z:20.12g}")

    ramp_time = 0.75 * run_time
    lines.extend(
        [
            "/PART/1",
            "beam_part",
            fmt_i(1, 1, 0),
            "/BRICK/1",
        ]
    )
    for eid, conn in mesh.bricks:
        lines.append(fmt_i(eid, *conn))

    lines.extend(
        [
            "/PROP/SOLID/1",
            "heph_constant_pressure",
            "#   Isolid    Ismstr               Icpre               Inpts    Itetra    Iframe                  dn",
            fmt_i(0, 0) + f"{1:20d}{0:20d}{0:10d}{0:10d}{0:20d}",
            "#                q_a                 q_b                   h            LAMBDA_V                MU_V",
            fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
            "#             dt_min   istrain      IHKT",
            fmt_f(0.0) + fmt_i(0, 0),
            "/BCS/1",
            "support_left_vertical",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   001 000{0:10d}{100:10d}",
            "/BCS/2",
            "support_right_vertical",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   001 000{0:10d}{101:10d}",
            "/BCS/3",
            "left_midwidth_anchor_xy",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   110 000{0:10d}{102:10d}",
            "/FUNCT/1",
            "smooth_load_ramp",
            "#                  X                   Y",
            fmt_f(0.0, 0.0),
            fmt_f(ramp_time, 1.0),
            fmt_f(run_time, 1.0),
        ]
    )
    for cid, (gid, set_name, force_per_node) in enumerate(load_sets, start=1):
        lines.extend(
            [
                f"/CLOAD/{cid}",
                f"{load_case}_{set_name}",
                "#funct_IDT       Dir   skew_ID sensor_ID  grnod_ID                       Ascalex             Fscaley",
                f"{1:10d}{'Z':>10}{0:10d}{0:10d}{gid:10d}{1.0:30.12g}{force_per_node:20.12g}",
            ]
        )

    lines.extend(node_group_block(100, "support_left", mesh.node_sets["support_left"]))
    lines.extend(node_group_block(101, "support_right", mesh.node_sets["support_right"]))
    lines.extend(node_group_block(102, "anchor_left", mesh.node_sets["anchor_left"]))
    lines.extend(node_group_block(200, "load_3pt", mesh.node_sets["load_3pt"]))
    lines.extend(node_group_block(210, "load_4pt_a", mesh.node_sets["load_4pt_a"]))
    lines.extend(node_group_block(211, "load_4pt_b", mesh.node_sets["load_4pt_b"]))
    lines.extend(node_group_block(300, "all_nodes", mesh.node_sets["all_nodes"]))
    lines.extend(
        [
            "/TH/NODE/1",
            "midspan_probe",
            "#     var1      var2      var3      var4      var5      var6      var7      var8      var9     var10",
            "DEF",
            "#    NODid     Iskew                                           NODname",
            f"{mesh.probe_node:10d}{0:10d}midspan",
            "/END",
            "",
        ]
    )
    out_rad.write_text("\n".join(lines), encoding="utf-8")


def write_engine(job: str, out_rad: pathlib.Path, run_time: float) -> None:
    lines = [
        "#RADIOSS ENGINE",
        "/ANIM/DT",
        "#   TSTART     TFREQ",
        fmt_f(run_time, run_time),
        "/ANIM/VECT/DISP",
        "/ANIM/VECT/VEL",
        "/ANIM/BRICK/TENS/STRESS/ALL",
        "/ANIM/BRICK/TENS/STRAIN/ALL",
        "/ANIM/GZIP",
        "/TFILE/4",
        "#            dT_HIS",
        fmt_f(run_time / 200.0),
        "/RFILE",
        "#   NCYCLE",
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
    anim_plain = workdir / f"{job}A001"
    anim_gz = workdir / f"{job}A001.gz"
    if anim_gz.exists():
        gunzip_keep(anim_gz, anim_plain)
    if not anim_plain.exists():
        raise FileNotFoundError(f"OpenRadioss animation frame not found: {anim_plain} or {anim_gz}")

    log.log("$ " + " ".join([str(ANIM_TO_VTK), str(anim_plain)]))
    proc = subprocess.run(
        [str(ANIM_TO_VTK), str(anim_plain)],
        cwd=str(workdir),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        if proc.stderr:
            log.log(proc.stderr.rstrip())
        raise RuntimeError(f"anim_to_vtk failed for {anim_plain}")
    final_vtk = workdir / f"{job}A001.vtk"
    if proc.stdout.lstrip().startswith("# vtk"):
        final_vtk.write_text(proc.stdout, encoding="utf-8")
        if proc.stderr:
            log.log(proc.stderr.rstrip())
        log.log(f"[exit {proc.returncode}] wrote {final_vtk}")
        return final_vtk
    candidates = sorted(workdir.glob(f"{job}A001*.vtk")) + sorted(workdir.glob("*.vtk"))
    if not candidates:
        raise FileNotFoundError(f"anim_to_vtk produced no .vtk in {workdir}")
    vtk = candidates[0]
    if vtk != final_vtk:
        vtk.replace(final_vtk)
    log.log(f"[exit {proc.returncode}] wrote {final_vtk}")
    return final_vtk


def extract_midspan_uz(vtk_path: pathlib.Path, probe_node: int, probe_xyz: tuple[float, float, float]) -> float:
    import pyvista as pv
    import numpy as np

    grid = pv.read(str(vtk_path))
    node_ids = grid.point_data.get("NODE_ID")
    if node_ids is not None:
        matches = np.flatnonzero(node_ids == probe_node)
        if len(matches):
            idx = int(matches[0])
        else:
            idx = grid.find_closest_point(probe_xyz)
    else:
        idx = grid.find_closest_point(probe_xyz)
    for name in ("DISPLACEMENT", "Displacement", "disp", "DISP", "U"):
        arr = grid.point_data.get(name)
        if arr is not None:
            return float(arr[idx][2])
    # The built-in converter may write deformed coordinates only. In that case,
    # use the difference between the final point and the original probe location.
    point = grid.points[idx]
    return float(point[2] - probe_xyz[2])


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
    out = RESULTS_DIR / "timeseries.csv"
    fields = [
        "stage",
        "load_case",
        "mesh",
        "elements",
        "nodes",
        "uz_fem_m",
        "delta_eb_m",
        "delta_timoshenko_m",
        "rel_error_eb",
        "rel_error_timoshenko",
        "verdict",
        "wall_clock_s",
    ]
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return out


def write_results_json(rows: list[dict[str, object]], case: BeamCase) -> pathlib.Path:
    gated = [r for r in rows if r["mesh"] == "M1"]
    verdict = "PASS" if gated and all(r["verdict"] == "PASS" for r in gated) else "FAIL"
    metrics = {
        f"{r['load_case']}_{r['mesh']}": {
            "uz_fem_m": r["uz_fem_m"],
            "delta_eb_m": r["delta_eb_m"],
            "delta_timoshenko_m": r["delta_timoshenko_m"],
            "rel_error_eb": r["rel_error_eb"],
            "rel_error_timoshenko": r["rel_error_timoshenko"],
            "elements": r["elements"],
            "nodes": r["nodes"],
            "wall_clock_s": r["wall_clock_s"],
        }
        for r in rows
    }
    payload = {
        "stage": STAGE,
        "name": STAGE_NAME,
        "verdict": verdict,
        "metrics": metrics,
        "reference": {
            "type": "Euler-Bernoulli closed form",
            "three_point_delta_m": case.eb_deflection("3pt"),
            "four_point_delta_m": case.eb_deflection("4pt"),
        },
        "tolerance": {"relative_error_eb": 0.01, "gated_mesh": "M1", "load_cases": ["3pt", "4pt"]},
        "git_sha": git_sha(),
    }
    out = RESULTS_DIR / "results.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return out


def write_typst_figure(rows: list[dict[str, object]]) -> pathlib.Path:
    out = FIGURES_DIR / "stage01_convergence.typ"
    csv_rel = pathlib.Path("..") / "results" / "timeseries.csv"
    ymax = max(0.025, 1.15 * max(float(r["rel_error_eb"]) for r in rows))
    x0 = 1.25
    y0 = 0.85
    width = 13.25
    height = 5.8
    spacing = width / len(rows)
    bar_width = 0.68

    def ymap(value: float) -> float:
        return y0 + (value / ymax) * height

    tick_lines: list[str] = []
    for tick in (0.0, 0.005, 0.010, 0.015, 0.020, 0.025):
        if tick > ymax * 1.001:
            continue
        y = ymap(tick)
        stroke = "horseshoe + 0.9pt" if math.isclose(tick, 0.010) else "black10 + 0.45pt"
        tick_lines.extend(
            [
                f"  line(({x0:.3f}, {y:.3f}), ({x0 + width:.3f}, {y:.3f}), stroke: {stroke})",
                f"  content(({x0 - 0.16:.3f}, {y:.3f}), [{tick * 100:.1f}\\%], anchor: \"east\")",
            ]
        )

    bar_lines: list[str] = []
    table_rows = "\n".join(
        f"  [{r['load_case']}], [{r['mesh']}], [{100.0 * float(r['rel_error_eb']):.3f}\\%], [{r['verdict']}],"
        for r in rows
    )
    for idx, r in enumerate(rows):
        err = float(r["rel_error_eb"])
        x = x0 + spacing * (idx + 0.5)
        yerr = ymap(err)
        color = "horseshoe" if r["verdict"] == "PASS" else "garnet"
        label = f"{r['load_case']} {r['mesh']}"
        bar_lines.extend(
            [
                f"  rect(({x - bar_width / 2:.3f}, {y0:.3f}), ({x + bar_width / 2:.3f}, {yerr:.3f}), fill: {color}, stroke: none)",
                f"  line(({x - bar_width / 2:.3f}, {yerr:.3f}), ({x + bar_width / 2:.3f}, {yerr:.3f}), stroke: charcoal + 0.6pt)",
                f"  content(({x:.3f}, {yerr + 0.24:.3f}), [{100.0 * err:.3f}\\%], anchor: \"south\")",
                f"  content(({x:.3f}, {y0 - 0.32:.3f}), [{label}], anchor: \"north\")",
            ]
        )

    lines = [
        '#import "@preview/cetz:0.3.4"',
        "",
        '#set page(width: 180mm, height: 118mm, margin: 10mm)',
        '#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))',
        "",
        '#let garnet = rgb("#73000A")',
        '#let charcoal = rgb("#363636")',
        '#let black70 = rgb("#5C5C5C")',
        '#let black50 = rgb("#A2A2A2")',
        '#let black10 = rgb("#ECECEC")',
        '#let horseshoe = rgb("#65780B")',
        '#let atlantic = rgb("#466A9F")',
        '#let white = rgb("#FFFFFF")',
        "",
        "#align(center)[#text(size: 11pt, weight: \"bold\")[Stage 01 Beam Bending Error Gate]]",
        "#v(2mm)",
        "",
        "#cetz.canvas(length: 1cm, {",
        "  import cetz.draw: *",
        "",
        f"  rect(({x0:.3f}, {y0:.3f}), ({x0 + width:.3f}, {y0 + height:.3f}), fill: white, stroke: charcoal + 0.65pt)",
        *tick_lines,
        *bar_lines,
        f"  content(({x0 + width / 2:.3f}, {y0 - 0.78:.3f}), [Mesh / load case], anchor: \"north\")",
        f"  content(({x0 - 0.95:.3f}, {y0 + height / 2:.3f}), [Euler-Bernoulli relative error], angle: 90deg)",
        f"  content(({x0 + 0.18:.3f}, {ymap(0.010) + 0.24:.3f}), [1\\% pass gate], anchor: \"west\")",
        f"  rect(({x0 + 8.50:.3f}, {y0 + height + 0.28:.3f}), ({x0 + 8.95:.3f}, {y0 + height + 0.58:.3f}), fill: garnet, stroke: none)",
        f"  content(({x0 + 9.07:.3f}, {y0 + height + 0.43:.3f}), [Coarse mesh], anchor: \"west\")",
        f"  rect(({x0 + 11.05:.3f}, {y0 + height + 0.28:.3f}), ({x0 + 11.50:.3f}, {y0 + height + 0.58:.3f}), fill: horseshoe, stroke: none)",
        f"  content(({x0 + 11.62:.3f}, {y0 + height + 0.43:.3f}), [Gated mesh], anchor: \"west\")",
        "})",
        "",
        f"#text(size: 8pt, fill: black70)[Source CSV: #raw(\"{csv_rel.as_posix()}\")]",
        "",
        "#figure(",
        "  table(",
        "    columns: 4,",
        "    [Load case], [Mesh], [Euler-Bernoulli error], [Verdict],",
        table_rows,
        "  ),",
        "  caption: [Stage 01 midspan displacement error by load case and mesh.]",
        ")",
        "",
    ]
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def run_case(
    case: BeamCase,
    mesh_label: str,
    load_case: str,
    run_time: float,
    n_threads: int,
    log: Logger,
    env: dict[str, str],
) -> dict[str, object]:
    spec = MESHES[mesh_label]
    mesh = build_mesh_data(case, spec)
    job = f"stage01_{load_case}_{mesh_label}"
    workdir = RUNS_DIR / job
    workdir.mkdir(parents=True, exist_ok=True)
    msh = workdir / f"{job}.msh"
    starter = workdir / f"{job}_0000.rad"
    engine = workdir / f"{job}_0001.rad"

    log.log(f"[{job}] generating structured HEXA8 mesh: {len(mesh.bricks)} elements, {len(mesh.nodes)} nodes")
    write_gmsh_mesh(case, mesh, msh)
    write_starter(case, mesh, load_case, starter, run_time)
    write_engine(job, engine, run_time)

    start = time.perf_counter()
    proc = log.run([str(STARTER), "-i", str(starter), "-nt", str(n_threads)], cwd=workdir, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"starter failed for {job}")
    proc = log.run([str(ENGINE), "-i", str(engine), "-nt", str(n_threads)], cwd=workdir, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"engine failed for {job}")
    vtk = convert_anim_to_vtk(job, workdir, log, env)
    wall = time.perf_counter() - start

    probe_xyz = mesh.nodes[mesh.probe_node]
    uz = extract_midspan_uz(vtk, mesh.probe_node, probe_xyz)
    delta_eb = -case.eb_deflection(load_case)
    delta_tim = -case.timoshenko_deflection(load_case)
    rel_eb = abs(uz - delta_eb) / abs(delta_eb)
    rel_tim = abs(uz - delta_tim) / abs(delta_tim)
    verdict = "PASS" if rel_eb <= 0.01 else "FAIL"
    log.log(
        f"[{job}] uz={uz:.8e} m, EB={delta_eb:.8e} m, "
        f"error={100.0 * rel_eb:.3f}%, verdict={verdict}"
    )
    return {
        "stage": STAGE,
        "load_case": load_case,
        "mesh": mesh_label,
        "elements": len(mesh.bricks),
        "nodes": len(mesh.nodes),
        "uz_fem_m": uz,
        "delta_eb_m": delta_eb,
        "delta_timoshenko_m": delta_tim,
        "rel_error_eb": rel_eb,
        "rel_error_timoshenko": rel_tim,
        "verdict": verdict,
        "wall_clock_s": wall,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meshes", default="M0,M1", help="comma-separated subset of M0,M1,M2")
    parser.add_argument("--load-cases", default="3pt,4pt", help="comma-separated subset of 3pt,4pt")
    parser.add_argument("--run-time", type=float, default=0.006, help="explicit simulation time in seconds")
    parser.add_argument("--n-threads", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    case = BeamCase()
    mesh_labels = [m.strip() for m in args.meshes.split(",") if m.strip()]
    load_cases = [c.strip() for c in args.load_cases.split(",") if c.strip()]
    for label in mesh_labels:
        if label not in MESHES:
            raise SystemExit(f"unknown mesh {label}")
    for load_case in load_cases:
        if load_case not in {"3pt", "4pt"}:
            raise SystemExit(f"unknown load case {load_case}")

    log = Logger(RUN_LOG)
    try:
        log.log(f"Stage 01 {STAGE_NAME}")
        log.log(f"OpenRadioss root: {OR_DIR}")
        log.log(f"Euler-Bernoulli 3pt: {-case.eb_deflection('3pt'):.8e} m")
        log.log(f"Euler-Bernoulli 4pt: {-case.eb_deflection('4pt'):.8e} m")
        if args.dry_run:
            return 0

        env = radioss_env()
        rows: list[dict[str, object]] = []
        for mesh_label in mesh_labels:
            for load_case in load_cases:
                rows.append(run_case(case, mesh_label, load_case, args.run_time, args.n_threads, log, env))

        write_timeseries(rows)
        results_path = write_results_json(rows, case)
        write_typst_figure(rows)
        payload = json.loads(results_path.read_text(encoding="utf-8"))
        return 0 if payload["verdict"] == "PASS" else 1
    finally:
        log.close()


if __name__ == "__main__":
    sys.exit(main())
