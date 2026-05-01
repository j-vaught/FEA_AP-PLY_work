#!/usr/bin/env python3
"""Stage 16 post-M5 ballistic runner.

This runner replaces the stale pre-M5 Stage 16 path. It reads the current
`kok_geom` mesh and orientation sidecar directly, splits the 24-ply panel on
the 23 interior ply planes so `/INTER/TYPE2` cohesive interfaces are real
surface pairs, writes a native OpenRadioss deck on the verified
`LAW12 + TYPE6 + /INIBRI/ORTHO + /FAIL/HASHIN` solid-composite route, and can
execute the Phase B single-shot wall-clock probe at `-nt 32`.

The current CLI is intentionally narrow:

    python tests/stage_16_PW_panel_ballistic/runner.py

writes the 200 mm Phase B deck, runs starter + engine unless disabled, converts
`T01` to CSV when available, and writes a JSON summary into the run directory.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import signal
import shutil
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import meshio
import numpy as np


THIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parents[1]
GEOM_DIR = THIS_DIR / "geometry"
RUNS_DIR = THIS_DIR / "runs" / "stage_16_post_m5"
RESULTS_DIR = THIS_DIR / "results"

CARD_PATH = ROOT_DIR / "references" / "material_cards" / "im7_8552.json"
RESULTS_JSON = RESULTS_DIR / "results.json"

OR_ROOT = Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"
ENGINE = OR_ROOT / "exec" / "engine_linux64_gf"
TH_TO_CSV = OR_ROOT / "exec" / "th_to_csv_linux64_gf"

SIM_DURATION_S = 0.5e-3
ANIM_DT_S = 5.0e-6
TH_DT_S = 1.0e-6
ROUND_DIGITS = 12
GLOBAL_TH_COLUMNS = 16
NC_RE = re.compile(
    r"NC=\s*(?P<cycle>\d+)\s+T=\s*(?P<time>[0-9.E+-]+)\s+DT=\s*(?P<dt>[0-9.E+-]+)"
)
ELAPSED_RE = re.compile(
    r"ELAPSED TIME=\s*(?P<elapsed>[0-9.E+-]+)\s*s\s+REMAINING TIME=\s*(?P<remaining>[0-9.E+-]+)\s*s"
)


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
class Cohesive:
    gic: float
    giic: float
    tn: float
    ts: float
    eta_bk: float
    penalty: float


@dataclass(frozen=True)
class ProjectileSpec:
    radius_m: float = 6.35e-3
    mass_kg: float = 13.4e-3
    density_kg_m3: float = 7850.0
    elastic_modulus_pa: float = 210.0e9
    poisson: float = 0.30
    gap_to_panel_m: float = 0.50e-3
    cross_cells: int = 8
    axial_layers: int = 8


@dataclass
class StructuredPanel:
    mesh_path: Path
    orientations_path: Path
    points: np.ndarray
    tetra10: np.ndarray
    physicals: np.ndarray
    physical_to_name: dict[int, str]
    group_by_name: dict[str, dict[str, Any]]
    x_values: np.ndarray
    y_values: np.ndarray
    z_values: np.ndarray
    node_i: np.ndarray
    node_j: np.ndarray
    node_k: np.ndarray
    node_grid: np.ndarray
    interface_plane_ks: tuple[int, ...]
    strike_plane_k: int
    back_plane_k: int
    centre_i: int
    centre_j: int
    base_node_count: int


@dataclass
class ProjectileMesh:
    points: np.ndarray
    tetra4: np.ndarray
    slave_node_ids: list[int]
    probe_local_node_id: int
    length_m: float
    voxel_area_m2: float


@dataclass
class DeckMetadata:
    starter_path: Path
    engine_path: Path
    base_nodes: int
    duplicated_interface_nodes: int
    total_nodes: int
    panel_tetra10: int
    panel_composite_tetra10: int
    panel_resin_tetra10: int
    projectile_tetra4: int
    interface_count: int
    interface_surface_segments: int
    strike_nodes: int
    clamp_nodes: int
    projectile_length_m: float


@dataclass
class PhaseBRunSummary:
    starter_path: str
    engine_path: str
    threads: int
    impact_velocity_m_s: float
    starter_rc: int
    engine_rc: int | None
    starter_wall_clock_s: float
    engine_wall_clock_s: float
    total_wall_clock_s: float
    stop_at_cycle: int | None
    checkpoint_cycle: int | None
    checkpoint_time_s: float | None
    checkpoint_dt_s: float | None
    checkpoint_elapsed_s: float | None
    checkpoint_remaining_s: float | None
    t01_csv_written: bool
    parse_status: str
    residual_velocity_m_s: float | None
    peak_contact_force_n: float | None
    contact_duration_us: float | None
    back_face_deflection_mm: float | None
    hourglass_ratio_max: float | None
    time_end_s: float | None
    deck: dict[str, Any]


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


def load_materials() -> tuple[Material, Cohesive]:
    data = json.loads(CARD_PATH.read_text(encoding="utf-8"))
    elastic = data["elastic"]
    strength = data["strength"]
    cohesive = data["cohesive"]
    material = Material(
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
    cohesive_card = Cohesive(
        gic=float(cohesive["GIc_J_m2"]),
        giic=float(cohesive["GIIc_J_m2"]),
        tn=float(cohesive["Tn_Pa"]),
        ts=float(cohesive["Ts_Pa"]),
        eta_bk=float(cohesive["eta_BK"]),
        penalty=float(cohesive["penalty_stiffness_N_m3"]),
    )
    return material, cohesive_card


def load_panel(mesh_path: Path, orientations_path: Path) -> StructuredPanel:
    mesh = meshio.read(mesh_path)
    tetra_block_index = next(index for index, block in enumerate(mesh.cells) if block.type == "tetra10")
    tetra10 = np.asarray(mesh.cells[tetra_block_index].data, dtype=np.int32)
    physicals = np.asarray(mesh.cell_data["gmsh:physical"][tetra_block_index], dtype=np.int32)
    points = np.asarray(mesh.points, dtype=np.float64)

    physical_to_name = {
        int(meta[0]): name
        for name, meta in mesh.field_data.items()
        if int(meta[1]) == 3
    }
    sidecar = json.loads(orientations_path.read_text(encoding="utf-8"))
    group_by_name = {str(group["name"]): group for group in sidecar["groups"]}

    x_values = np.unique(np.round(points[:, 0], ROUND_DIGITS))
    y_values = np.unique(np.round(points[:, 1], ROUND_DIGITS))
    z_values = np.unique(np.round(points[:, 2], ROUND_DIGITS))
    x_index = {float(value): idx for idx, value in enumerate(x_values)}
    y_index = {float(value): idx for idx, value in enumerate(y_values)}
    z_index = {float(value): idx for idx, value in enumerate(z_values)}

    base_node_count = points.shape[0]
    node_i = np.zeros(base_node_count + 1, dtype=np.int16)
    node_j = np.zeros(base_node_count + 1, dtype=np.int16)
    node_k = np.zeros(base_node_count + 1, dtype=np.int16)
    node_grid = np.zeros((len(z_values), len(y_values), len(x_values)), dtype=np.int32)

    for node_id, xyz in enumerate(points, start=1):
        ix = x_index[round(float(xyz[0]), ROUND_DIGITS)]
        iy = y_index[round(float(xyz[1]), ROUND_DIGITS)]
        iz = z_index[round(float(xyz[2]), ROUND_DIGITS)]
        node_i[node_id] = ix
        node_j[node_id] = iy
        node_k[node_id] = iz
        node_grid[iz, iy, ix] = node_id

    interface_plane_ks = tuple(range(2, len(z_values) - 1, 2))
    strike_plane_k = len(z_values) - 1
    back_plane_k = 0
    centre_i = len(x_values) // 2
    centre_j = len(y_values) // 2

    return StructuredPanel(
        mesh_path=mesh_path,
        orientations_path=orientations_path,
        points=points,
        tetra10=tetra10,
        physicals=physicals,
        physical_to_name=physical_to_name,
        group_by_name=group_by_name,
        x_values=x_values,
        y_values=y_values,
        z_values=z_values,
        node_i=node_i,
        node_j=node_j,
        node_k=node_k,
        node_grid=node_grid,
        interface_plane_ks=interface_plane_ks,
        strike_plane_k=strike_plane_k,
        back_plane_k=back_plane_k,
        centre_i=centre_i,
        centre_j=centre_j,
        base_node_count=base_node_count,
    )


def duplicate_interface_planes(panel: StructuredPanel) -> tuple[dict[int, np.ndarray], int]:
    duplicate_planes: dict[int, np.ndarray] = {}
    next_node_id = panel.base_node_count + 1
    for plane_k in panel.interface_plane_ks:
        template = panel.node_grid[plane_k]
        count = template.size
        ids = np.arange(next_node_id, next_node_id + count, dtype=np.int32).reshape(template.shape)
        duplicate_planes[plane_k] = ids
        next_node_id += count
    return duplicate_planes, next_node_id - 1


def build_projectile_mesh(strike_z: float, spec: ProjectileSpec) -> ProjectileMesh:
    cross_cells = spec.cross_cells
    axial_layers = spec.axial_layers
    dx = (2.0 * spec.radius_m) / cross_cells

    selected_cells: list[tuple[int, int]] = []
    for iy in range(cross_cells):
        y0 = -spec.radius_m + iy * dx
        yc = y0 + 0.5 * dx
        for ix in range(cross_cells):
            x0 = -spec.radius_m + ix * dx
            xc = x0 + 0.5 * dx
            if xc * xc + yc * yc <= spec.radius_m * spec.radius_m:
                selected_cells.append((ix, iy))

    if not selected_cells:
        raise RuntimeError("projectile cross-section discretization selected no cells")

    voxel_area = len(selected_cells) * dx * dx
    dz = spec.mass_kg / (spec.density_kg_m3 * voxel_area * axial_layers)
    length = dz * axial_layers
    z0 = strike_z + spec.gap_to_panel_m

    node_lookup: dict[tuple[int, int, int], int] = {}
    points: list[tuple[float, float, float]] = []
    tetra4: list[tuple[int, int, int, int]] = []

    def node_id(ix: int, iy: int, iz: int) -> int:
        key = (ix, iy, iz)
        existing = node_lookup.get(key)
        if existing is not None:
            return existing
        x = -spec.radius_m + ix * dx
        y = -spec.radius_m + iy * dx
        z = z0 + iz * dz
        new_id = len(points) + 1
        node_lookup[key] = new_id
        points.append((x, y, z))
        return new_id

    for iz in range(axial_layers):
        for ix, iy in selected_cells:
            n000 = node_id(ix, iy, iz)
            n100 = node_id(ix + 1, iy, iz)
            n110 = node_id(ix + 1, iy + 1, iz)
            n010 = node_id(ix, iy + 1, iz)
            n001 = node_id(ix, iy, iz + 1)
            n101 = node_id(ix + 1, iy, iz + 1)
            n111 = node_id(ix + 1, iy + 1, iz + 1)
            n011 = node_id(ix, iy + 1, iz + 1)
            tetra4.extend(
                (
                    (n000, n100, n110, n111),
                    (n000, n110, n010, n111),
                    (n000, n010, n011, n111),
                    (n000, n011, n001, n111),
                    (n000, n001, n101, n111),
                    (n000, n101, n100, n111),
                )
            )

    slave_node_ids = list(range(1, len(points) + 1))
    probe_local_node_id = node_lookup[(cross_cells // 2, cross_cells // 2, axial_layers // 2)]
    return ProjectileMesh(
        points=np.asarray(points, dtype=np.float64),
        tetra4=np.asarray(tetra4, dtype=np.int32),
        slave_node_ids=slave_node_ids,
        probe_local_node_id=probe_local_node_id,
        length_m=length,
        voxel_area_m2=voxel_area,
    )


def law12_block(material: Material) -> list[str]:
    nu31 = material.nu13 * material.e3 / material.e1
    return [
        "/MAT/LAW12/1",
        "IM7_8552_tow_LAW12",
        fmt_f(material.rho),
        fmt_f(material.e1, material.e2, material.e3),
        fmt_f(material.nu12, material.nu23, nu31),
        fmt_f(material.g12, material.g23, material.g13),
        fmt_f(material.xt, material.yt, material.zt, 0.05),
        fmt_f(1.0, 1.0, 1.0, 1.0),
        fmt_f(material.xt, material.yt, material.xc, material.yc),
        fmt_f(material.s12, material.s12, material.s23, material.s23),
        fmt_f(material.zt, material.zc, material.s13, material.s13),
        fmt_f(0.0, 0.0, 0.0, 0.0) + fmt_i(1),
    ]


def hashin_block(material: Material, pthickfail: float = 1.0) -> list[str]:
    return [
        "/FAIL/HASHIN/1",
        fmt_i(1, 0, 1) + f"{pthickfail:20.12g}",
        fmt_f(material.xt, material.yt, material.zt, material.xc, material.yc),
        fmt_f(material.zc, material.s12, material.s12, material.s23, material.s13),
        fmt_f(0.0, 1.0, 0.0, 0.0),
    ]


def resin_block() -> list[str]:
    return [
        "/MAT/ELAST/2",
        "8552_resin_isotropic",
        fmt_f(1300.0, 0.0),
        fmt_f(4.67e9, 0.38),
    ]


def projectile_block(spec: ProjectileSpec) -> list[str]:
    return [
        "/MAT/ELAST/3",
        "projectile_steel_elastic",
        fmt_f(spec.density_kg_m3, 0.0),
        fmt_f(spec.elastic_modulus_pa, spec.poisson),
    ]


def type6_property(prop_id: int) -> list[str]:
    return [
        f"/PROP/TYPE6/{prop_id}",
        f"type6_sol_orth_inibri_{prop_id}",
        fmt_i(24, 4) + f"{1:20d}{1000:10d}{0:10d}{0:10d}{2:10d}{0.0:20.12g}",
        fmt_f(0.0, 0.0, 0.0),
        fmt_f(1.0, 0.0, 0.0) + fmt_i(0, 1, 0),
        fmt_f(0.0, 0.0, 0.0, 0.0),
        fmt_f(0.0) + fmt_i(0, 0),
    ]


def solid_property(prop_id: int) -> list[str]:
    return [
        f"/PROP/SOLID/{prop_id}",
        f"solid_prop_{prop_id}",
        fmt_i(0, 0) + f"{0:20d}{0:10d}{0:10d}{3:10d}{0:10d}{0.0:20.12g}",
        fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
        fmt_f(0.0),
    ]


def write_lines(handle, lines: list[str]) -> None:
    for line in lines:
        handle.write(line)
        handle.write("\n")


def group_block(group_id: int, name: str, node_ids: list[int]) -> list[str]:
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for start in range(0, len(node_ids), 10):
        lines.append(fmt_i(*node_ids[start : start + 10]))
    return lines


def remap_panel_element_nodes(
    panel: StructuredPanel,
    duplicate_planes: dict[int, np.ndarray],
    conn: np.ndarray,
    bottom_plane_k: int,
) -> list[int]:
    node_ids = [int(node) + 1 for node in conn]
    if bottom_plane_k <= 0:
        return node_ids
    remapped: list[int] | None = None
    for local_idx, node_id in enumerate(node_ids):
        if int(panel.node_k[node_id]) != bottom_plane_k:
            continue
        if remapped is None:
            remapped = list(node_ids)
        remapped[local_idx] = int(
            duplicate_planes[bottom_plane_k][
                int(panel.node_j[node_id]),
                int(panel.node_i[node_id]),
            ]
        )
    return remapped if remapped is not None else node_ids


def boundary_node_ids(panel: StructuredPanel, duplicate_planes: dict[int, np.ndarray]) -> list[int]:
    base_ids = np.unique(
        np.concatenate(
            [
                panel.node_grid[:, :, 0].ravel(),
                panel.node_grid[:, :, -1].ravel(),
                panel.node_grid[:, 0, :].ravel(),
                panel.node_grid[:, -1, :].ravel(),
            ]
        )
    )
    duplicate_boundary_ids: list[np.ndarray] = []
    for plane_ids in duplicate_planes.values():
        duplicate_boundary_ids.extend(
            [
                plane_ids[:, 0].ravel(),
                plane_ids[:, -1].ravel(),
                plane_ids[0, :].ravel(),
                plane_ids[-1, :].ravel(),
            ]
        )
    if duplicate_boundary_ids:
        all_ids = np.unique(np.concatenate([base_ids, *duplicate_boundary_ids]))
    else:
        all_ids = base_ids
    return [int(node_id) for node_id in all_ids]


def write_surface_seg(handle, surf_id: int, name: str, plane_ids: np.ndarray, reverse: bool) -> int:
    handle.write(f"/SURF/SEG/{surf_id}\n")
    handle.write(f"{name}\n")
    seg_id = 1
    ny, nx = plane_ids.shape
    for j in range(0, ny - 1, 2):
        for i in range(0, nx - 1, 2):
            n1 = int(plane_ids[j, i])
            n2 = int(plane_ids[j, i + 2])
            n3 = int(plane_ids[j + 2, i + 2])
            n4 = int(plane_ids[j + 2, i])
            if reverse:
                handle.write(fmt_i(seg_id, n1, n4, n3, n2) + "\n")
            else:
                handle.write(fmt_i(seg_id, n1, n2, n3, n4) + "\n")
            seg_id += 1
    return seg_id - 1


def write_starter(
    panel: StructuredPanel,
    material: Material,
    cohesive: Cohesive,
    projectile_spec: ProjectileSpec,
    starter_path: Path,
    impact_velocity_m_s: float,
    pthickfail: float = 1.0,
) -> DeckMetadata:
    duplicate_planes, last_duplicate_node_id = duplicate_interface_planes(panel)
    projectile = build_projectile_mesh(float(panel.z_values[panel.strike_plane_k]), projectile_spec)
    projectile_node_offset = last_duplicate_node_id
    projectile_slave_global = [projectile_node_offset + node_id for node_id in projectile.slave_node_ids]
    projectile_probe_global = projectile_node_offset + projectile.probe_local_node_id
    total_nodes = projectile_node_offset + len(projectile.points)

    strike_nodes = [int(node_id) for node_id in panel.node_grid[panel.strike_plane_k].ravel()]
    clamp_nodes = boundary_node_ids(panel, duplicate_planes)
    back_face_center_node = int(panel.node_grid[panel.back_plane_k, panel.centre_j, panel.centre_i])

    physical_meta: dict[int, dict[str, Any]] = {}
    for physical_id, name in panel.physical_to_name.items():
        group = panel.group_by_name[name]
        ply_index = int(group["ply_index"])
        physical_meta[physical_id] = {
            "name": name,
            "ply_index": ply_index,
            "bottom_plane_k": 2 * (ply_index - 1),
            "isotropic": bool(group.get("isotropic", False)),
            "v1": tuple(group.get("fiber_direction_unit_vector", (0.0, 0.0, 0.0))),
            "v2": tuple(group.get("transverse_in_plane_unit_vector", (0.0, 0.0, 0.0))),
        }

    starter_path.parent.mkdir(parents=True, exist_ok=True)
    inibri_tmp = Path(tempfile.mkstemp(prefix="stage16_inibri_", suffix=".tmp", dir=str(starter_path.parent))[1])
    composite_count = 0
    resin_count = 0
    panel_element_count = len(panel.tetra10)
    interface_segments = 0

    try:
        with starter_path.open("w", encoding="utf-8") as handle, inibri_tmp.open("w", encoding="utf-8") as inibri:
            write_lines(
                handle,
                [
                    "#RADIOSS STARTER",
                    "/BEGIN",
                    starter_path.stem.removesuffix("_0000"),
                    "      2026         0",
                    f"{'kg':>20}{'m':>20}{'s':>20}",
                    f"{'kg':>20}{'m':>20}{'s':>20}",
                    "/TITLE",
                    "Stage 16 P1-TWT post-M5 Phase B single-shot",
                    "/DEF_SOLID",
                    fmt_i(24, 4) + f"{0:20d}{2:40d}",
                ],
            )
            write_lines(handle, law12_block(material))
            write_lines(handle, hashin_block(material, pthickfail))
            write_lines(handle, resin_block())
            write_lines(handle, projectile_block(projectile_spec))
            write_lines(handle, type6_property(1))
            write_lines(handle, type6_property(2))
            write_lines(handle, solid_property(3))

            handle.write("/NODE\n")
            for node_id, xyz in enumerate(panel.points, start=1):
                handle.write(
                    f"{node_id:10d}{float(xyz[0]):20.12g}{float(xyz[1]):20.12g}{float(xyz[2]):20.12g}\n"
                )
            for plane_k in panel.interface_plane_ks:
                original_plane = panel.node_grid[plane_k]
                duplicate_plane = duplicate_planes[plane_k]
                for j in range(original_plane.shape[0]):
                    for i in range(original_plane.shape[1]):
                        original_id = int(original_plane[j, i])
                        duplicate_id = int(duplicate_plane[j, i])
                        xyz = panel.points[original_id - 1]
                        handle.write(
                            f"{duplicate_id:10d}{float(xyz[0]):20.12g}{float(xyz[1]):20.12g}{float(xyz[2]):20.12g}\n"
                        )
            for local_id, xyz in enumerate(projectile.points, start=1):
                node_id = projectile_node_offset + local_id
                handle.write(
                    f"{node_id:10d}{float(xyz[0]):20.12g}{float(xyz[1]):20.12g}{float(xyz[2]):20.12g}\n"
                )

            element_id = 1

            write_lines(handle, ["/PART/1", "laminate_composite", fmt_i(1, 1, 0), "/TETRA10/1"])
            for conn, physical_id in zip(panel.tetra10, panel.physicals):
                meta = physical_meta[int(physical_id)]
                if meta["isotropic"]:
                    continue
                node_ids = remap_panel_element_nodes(
                    panel=panel,
                    duplicate_planes=duplicate_planes,
                    conn=conn,
                    bottom_plane_k=int(meta["bottom_plane_k"]),
                )
                handle.write(fmt_i(element_id) + "\n")
                handle.write(fmt_i(*node_ids) + "\n")
                v1 = meta["v1"]
                v2 = meta["v2"]
                inibri.write(fmt_i(element_id, 1, 10, 6, 24) + "\n")
                inibri.write(fmt_f(float(v1[0]), float(v1[1]), float(v1[2]), float(v2[0]), float(v2[1])) + "\n")
                inibri.write(fmt_f(float(v2[2])) + "\n")
                composite_count += 1
                element_id += 1

            write_lines(handle, ["/PART/2", "laminate_resin", fmt_i(2, 2, 0), "/TETRA10/2"])
            for conn, physical_id in zip(panel.tetra10, panel.physicals):
                meta = physical_meta[int(physical_id)]
                if not meta["isotropic"]:
                    continue
                node_ids = remap_panel_element_nodes(
                    panel=panel,
                    duplicate_planes=duplicate_planes,
                    conn=conn,
                    bottom_plane_k=int(meta["bottom_plane_k"]),
                )
                handle.write(fmt_i(element_id) + "\n")
                handle.write(fmt_i(*node_ids) + "\n")
                resin_count += 1
                element_id += 1

            write_lines(handle, ["/PART/3", "projectile_rigid_part", fmt_i(3, 3, 0), "/TETRA4/3"])
            for tetra in projectile.tetra4:
                global_nodes = [projectile_node_offset + int(node_id) for node_id in tetra]
                handle.write(fmt_i(element_id, *global_nodes) + "\n")
                element_id += 1

            handle.write("/INIBRI/ORTHO\n")
            with inibri_tmp.open("r", encoding="utf-8") as inibri_readback:
                shutil.copyfileobj(inibri_readback, handle)

            lower_surf_ids: list[int] = []
            upper_surf_ids: list[int] = []
            for interface_idx, plane_k in enumerate(panel.interface_plane_ks, start=1):
                lower_surf_id = 1000 + interface_idx
                upper_surf_id = 2000 + interface_idx
                lower_surf_ids.append(lower_surf_id)
                upper_surf_ids.append(upper_surf_id)
                interface_segments = write_surface_seg(
                    handle,
                    lower_surf_id,
                    f"iface_{interface_idx:02d}_lower",
                    panel.node_grid[plane_k],
                    reverse=False,
                )
                write_surface_seg(
                    handle,
                    upper_surf_id,
                    f"iface_{interface_idx:02d}_upper",
                    duplicate_planes[plane_k],
                    reverse=True,
                )

            write_lines(handle, ["/SURF/PART/ALL/3001", "projectile_all_surface", fmt_i(3)])
            write_lines(handle, group_block(4001, "clamped_cut_edges", clamp_nodes))
            write_lines(handle, group_block(4002, "panel_strike_face", strike_nodes))
            for interface_idx, plane_k in enumerate(panel.interface_plane_ks, start=1):
                write_lines(
                    handle,
                    group_block(
                        6000 + interface_idx,
                        f"iface_{interface_idx:02d}_secondary_nodes",
                        [int(node_id) for node_id in duplicate_planes[plane_k].ravel()],
                    ),
                )

            write_lines(
                handle,
                [
                    "/RBODY/1",
                    "projectile_rigid_body",
                    fmt_i(0, 0, 0, 0) + f"{projectile_spec.mass_kg:20.12g}{4003:10d}{0:10d}{0:10d}{0:10d}",
                    fmt_f(0.0, 0.0, 0.0),
                    fmt_f(0.0, 0.0, 0.0),
                    fmt_i(0, 0, 0),
                ],
            )
            write_lines(handle, group_block(4003, "projectile_slave_nodes", projectile_slave_global))

            write_lines(
                handle,
                [
                    "/INIVEL/TRA/1",
                    "projectile_initial_velocity",
                    fmt_f(0.0, 0.0, -impact_velocity_m_s) + fmt_i(4004, 0),
                ],
            )
            write_lines(handle, group_block(4004, "projectile_velocity_nodes", projectile_slave_global))

            write_lines(
                handle,
                [
                    "/BCS/1",
                    "panel_cut_edge_clamp",
                    f"   111 000{0:10d}{4001:10d}",
                ],
            )

            stfac1 = cohesive.penalty / 1.0e10
            for interface_idx, (upper_surf_id, lower_surf_id) in enumerate(zip(upper_surf_ids, lower_surf_ids), start=1):
                write_lines(
                    handle,
                    [
                        f"/INTER/TYPE2/{5000 + interface_idx}",
                        f"COH_PLY_{interface_idx:02d}_{interface_idx + 1:02d}",
                        fmt_i(6000 + interface_idx, lower_surf_id, 1000, 25, 0, 2, 1000, upper_surf_id)
                        + f"{0.0:20.12g}",
                        fmt_f(stfac1, 0.05) + f"{0:30d}",
                        fmt_f(cohesive.tn, cohesive.ts, cohesive.gic, cohesive.giic, cohesive.eta_bk),
                    ],
                )

            write_lines(
                handle,
                [
                    "/INTER/TYPE7/1",
                    "projectile_to_panel_contact",
                    fmt_i(4002, 3001, 0, 0, 2, 0, 0, 0, 0),
                    fmt_f(0.0, 0.0, 0.0),
                    fmt_f(0.0, 0.0, 0.0, 0.0) + fmt_i(0, 0),
                    fmt_f(0.0, 0.3, 0.0, 0.0, 0.0),
                    "       000" + f"{6:30d}{0.05:20.12g}{1.0:20.12g}{0.0:20.12g}",
                    fmt_i(0, 0) + f"{0.0:20.12g}{2:10d}{0:10d}{0:10d}{0.0:20.12g}{0:10d}",
                ],
            )

            write_lines(
                handle,
                [
                    "/TH/NODE/1",
                    "projectile_probe_history",
                    "DEF",
                    f"{projectile_probe_global:10d}{0:10d}projectile_probe",
                    "/TH/NODE/2",
                    "back_face_center_history",
                    "DEF",
                    f"{back_face_center_node:10d}{0:10d}back_face_center",
                    "/END",
                    "",
                ],
            )
    finally:
        if inibri_tmp.exists():
            inibri_tmp.unlink()

    duplicated_interface_nodes = last_duplicate_node_id - panel.base_node_count
    return DeckMetadata(
        starter_path=starter_path,
        engine_path=starter_path.with_name(starter_path.name.replace("_0000.rad", "_0001.rad")),
        base_nodes=panel.base_node_count,
        duplicated_interface_nodes=duplicated_interface_nodes,
        total_nodes=total_nodes,
        panel_tetra10=panel_element_count,
        panel_composite_tetra10=composite_count,
        panel_resin_tetra10=resin_count,
        projectile_tetra4=len(projectile.tetra4),
        interface_count=len(panel.interface_plane_ks),
        interface_surface_segments=interface_segments * len(panel.interface_plane_ks),
        strike_nodes=len(strike_nodes),
        clamp_nodes=len(clamp_nodes),
        projectile_length_m=projectile.length_m,
    )


def write_engine(engine_path: Path) -> Path:
    job = engine_path.stem.removesuffix("_0001")
    lines = [
        "#RADIOSS ENGINE",
        "/ANIM/DT",
        fmt_f(0.0, ANIM_DT_S),
        "/ANIM/GZIP",
        "/MON/ON",
        "/PARITH/ON",
        "/PRINT/-100/55",
        f"/RUN/{job}/1",
        fmt_f(SIM_DURATION_S),
        "/TFILE/4",
        fmt_f(TH_DT_S),
        "/VERS/2026",
        "",
    ]
    engine_path.write_text("\n".join(lines), encoding="utf-8")
    return engine_path


def parse_engine_progress(console_log: Path) -> dict[str, float | int] | None:
    if not console_log.exists():
        return None

    latest: dict[str, float | int] | None = None
    pending: dict[str, float | int] | None = None
    for line in console_log.read_text(encoding="utf-8", errors="replace").splitlines():
        nc_match = NC_RE.search(line)
        if nc_match is not None:
            pending = {
                "cycle": int(nc_match.group("cycle")),
                "time_s": float(nc_match.group("time")),
                "dt_s": float(nc_match.group("dt")),
            }
            continue
        elapsed_match = ELAPSED_RE.search(line)
        if elapsed_match is not None and pending is not None:
            latest = {
                **pending,
                "elapsed_s": float(elapsed_match.group("elapsed")),
                "remaining_s": float(elapsed_match.group("remaining")),
            }
            pending = None
    return latest


def run_solver(
    starter_path: Path,
    engine_path: Path,
    threads: int,
    starter_only: bool,
    stop_at_cycle: int | None = None,
) -> tuple[int, int | None, float, float, dict[str, float | int] | None]:
    env = radioss_env()
    workdir = starter_path.parent

    starter_cmd = [str(STARTER), "-i", starter_path.name, "-nt", str(threads)]
    starter_log = workdir / (starter_path.stem + ".console.log")
    starter_t0 = time.monotonic()
    with starter_log.open("w", encoding="utf-8") as log_handle:
        starter_proc = subprocess.run(
            starter_cmd,
            cwd=str(workdir),
            env=env,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    starter_wall = time.monotonic() - starter_t0
    if starter_proc.returncode != 0 or starter_only:
        return starter_proc.returncode, None, starter_wall, 0.0, None

    engine_cmd = [str(ENGINE), "-i", engine_path.name, "-nt", str(threads)]
    engine_log = workdir / (engine_path.stem + ".console.log")
    engine_t0 = time.monotonic()
    if stop_at_cycle is None:
        with engine_log.open("w", encoding="utf-8") as log_handle:
            engine_proc = subprocess.run(
                engine_cmd,
                cwd=str(workdir),
                env=env,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                check=False,
            )
        engine_wall = time.monotonic() - engine_t0
        return starter_proc.returncode, engine_proc.returncode, starter_wall, engine_wall, None

    checkpoint: dict[str, float | int] | None = None
    stop_sent = False
    with engine_log.open("w", encoding="utf-8") as log_handle:
        engine_proc = subprocess.Popen(
            engine_cmd,
            cwd=str(workdir),
            env=env,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
        )
        while True:
            rc = engine_proc.poll()
            checkpoint = parse_engine_progress(engine_log)
            if not stop_sent and checkpoint is not None and int(checkpoint["cycle"]) >= stop_at_cycle:
                engine_proc.send_signal(signal.SIGINT)
                stop_sent = True
            if rc is not None:
                break
            time.sleep(2.0)
        if stop_sent and engine_proc.returncode is None:
            try:
                engine_proc.wait(timeout=15.0)
            except subprocess.TimeoutExpired:
                engine_proc.terminate()
                try:
                    engine_proc.wait(timeout=10.0)
                except subprocess.TimeoutExpired:
                    engine_proc.kill()
                    engine_proc.wait()
        checkpoint = parse_engine_progress(engine_log)
    engine_wall = time.monotonic() - engine_t0
    return starter_proc.returncode, engine_proc.returncode, starter_wall, engine_wall, checkpoint


def convert_t01_to_csv(t01_path: Path) -> Path | None:
    if not t01_path.exists():
        return None
    proc = subprocess.run(
        [str(TH_TO_CSV), t01_path.name],
        cwd=str(t01_path.parent),
        env=radioss_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return None
    csv_path = t01_path.with_suffix(t01_path.suffix + ".csv")
    return csv_path if csv_path.exists() else None


def read_numeric_csv(csv_path: Path) -> tuple[list[str], np.ndarray]:
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows: list[list[float]] = []
        for row in reader:
            values = [float(value) for value in row if value.strip()]
            if values:
                rows.append(values)
    if not rows:
        raise RuntimeError(f"no numeric rows found in {csv_path}")
    return header, np.asarray(rows, dtype=np.float64)


def contact_duration_us(time_s: np.ndarray, force_n: np.ndarray) -> float:
    peak = float(np.max(np.abs(force_n))) if force_n.size else 0.0
    if peak <= 0.0:
        return 0.0
    mask = np.abs(force_n) >= 0.10 * peak
    if not np.any(mask):
        return 0.0
    active_times = time_s[mask]
    return float((active_times[-1] - active_times[0]) * 1.0e6)


def parse_phase_b_csv(csv_path: Path) -> dict[str, float]:
    _header, data = read_numeric_csv(csv_path)
    if data.shape[1] < GLOBAL_TH_COLUMNS + 12:
        raise RuntimeError(
            f"{csv_path.name} has {data.shape[1]} columns; expected at least {GLOBAL_TH_COLUMNS + 12}"
        )

    time_s = data[:, 0]
    internal = data[:, 1]
    hourglass = data[:, 12]
    offset = GLOBAL_TH_COLUMNS

    projectile_history = data[:, offset : offset + 6]
    back_face_history = data[:, offset + 6 : offset + 12]
    projectile_z = projectile_history[:, 2]
    projectile_vz = projectile_history[:, 5]
    back_face_z = back_face_history[:, 2]

    ratio = np.zeros_like(hourglass)
    positive_internal = internal > 0.0
    ratio[positive_internal] = hourglass[positive_internal] / internal[positive_internal]

    return {
        "time_end_s": float(time_s[-1]),
        "residual_velocity_m_s": float(max(0.0, -projectile_vz[-1])),
        "back_face_deflection_mm": float(np.max(np.abs(back_face_z - back_face_z[0])) * 1000.0),
        "hourglass_ratio_max": float(np.max(ratio)),
        "projectile_final_z_m": float(projectile_z[-1]),
    }


def run_phase_b(
    mesh_path: Path,
    orientations_path: Path,
    run_dir: Path,
    impact_velocity_m_s: float,
    threads: int,
    write_only: bool,
    starter_only: bool,
    stop_at_cycle: int | None,
    pthickfail: float = 1.0,
) -> PhaseBRunSummary:
    run_dir.mkdir(parents=True, exist_ok=True)

    material, cohesive = load_materials()
    panel = load_panel(mesh_path, orientations_path)
    projectile_spec = ProjectileSpec()

    starter_path = run_dir / "stage16_phase_b_single_shot_0000.rad"
    engine_path = run_dir / "stage16_phase_b_single_shot_0001.rad"

    deck = write_starter(
        panel=panel,
        material=material,
        cohesive=cohesive,
        projectile_spec=projectile_spec,
        starter_path=starter_path,
        impact_velocity_m_s=impact_velocity_m_s,
        pthickfail=pthickfail,
    )
    write_engine(engine_path)

    starter_rc = 0
    engine_rc: int | None = None
    starter_wall = 0.0
    engine_wall = 0.0
    checkpoint = None
    t01_csv_written = False
    parse_status = "NOT_RUN"
    residual_velocity = None
    peak_contact_force = None
    contact_duration = None
    back_face = None
    hourglass_ratio = None
    time_end_s = None

    if not write_only:
        starter_rc, engine_rc, starter_wall, engine_wall, checkpoint = run_solver(
            starter_path=starter_path,
            engine_path=engine_path,
            threads=threads,
            starter_only=starter_only,
            stop_at_cycle=stop_at_cycle,
        )
        if starter_rc == 0 and not starter_only and engine_rc == 0:
            t01_path = run_dir / "stage16_phase_b_single_shotT01"
            csv_path = convert_t01_to_csv(t01_path)
            if csv_path is not None:
                t01_csv_written = True
                try:
                    metrics = parse_phase_b_csv(csv_path)
                except Exception as exc:
                    parse_status = f"PARSE_FAILED: {exc}"
                else:
                    parse_status = "OK"
                    residual_velocity = metrics["residual_velocity_m_s"]
                    back_face = metrics["back_face_deflection_mm"]
                    hourglass_ratio = metrics["hourglass_ratio_max"]
                    time_end_s = metrics["time_end_s"]
            else:
                parse_status = "T01_TO_CSV_FAILED"
        elif starter_rc != 0:
            parse_status = "STARTER_FAILED"
        elif starter_only:
            parse_status = "STARTER_ONLY"
        elif checkpoint is not None:
            parse_status = "ENGINE_STOPPED_AT_CYCLE"
        else:
            parse_status = "ENGINE_FAILED"

    summary = PhaseBRunSummary(
        starter_path=str(starter_path),
        engine_path=str(engine_path),
        threads=threads,
        impact_velocity_m_s=impact_velocity_m_s,
        starter_rc=starter_rc,
        engine_rc=engine_rc,
        starter_wall_clock_s=starter_wall,
        engine_wall_clock_s=engine_wall,
        total_wall_clock_s=starter_wall + engine_wall,
        stop_at_cycle=stop_at_cycle,
        checkpoint_cycle=None if checkpoint is None else int(checkpoint["cycle"]),
        checkpoint_time_s=None if checkpoint is None else float(checkpoint["time_s"]),
        checkpoint_dt_s=None if checkpoint is None else float(checkpoint["dt_s"]),
        checkpoint_elapsed_s=None if checkpoint is None else float(checkpoint["elapsed_s"]),
        checkpoint_remaining_s=None if checkpoint is None else float(checkpoint["remaining_s"]),
        t01_csv_written=t01_csv_written,
        parse_status=parse_status,
        residual_velocity_m_s=residual_velocity,
        peak_contact_force_n=peak_contact_force,
        contact_duration_us=contact_duration,
        back_face_deflection_mm=back_face,
        hourglass_ratio_max=hourglass_ratio,
        time_end_s=time_end_s,
        deck=asdict(deck),
    )
    (run_dir / "phase_b_summary.json").write_text(
        json.dumps(asdict(summary), indent=2, default=str),
        encoding="utf-8",
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage 16 post-M5 Phase B runner")
    parser.add_argument("--mesh", type=Path, default=GEOM_DIR / "p1_twt_200.msh")
    parser.add_argument("--orientations", type=Path, default=GEOM_DIR / "orientations.json")
    parser.add_argument("--run-dir", type=Path, default=RUNS_DIR / "phase_b_single_shot")
    parser.add_argument("--velocity", type=float, default=100.0)
    parser.add_argument("--threads", type=int, default=int(os.environ.get("RAD_NT", "32")))
    parser.add_argument("--write-only", action="store_true")
    parser.add_argument("--starter-only", action="store_true")
    parser.add_argument("--stop-at-cycle", type=int, default=None)
    parser.add_argument("--pthickfail", type=float, default=1.0,
                        help="Hashin erosion threshold (1.0 = all IPs fail; 0.5 = half IPs fail).")
    args = parser.parse_args(argv)

    summary = run_phase_b(
        mesh_path=args.mesh,
        orientations_path=args.orientations,
        run_dir=args.run_dir,
        impact_velocity_m_s=args.velocity,
        threads=args.threads,
        write_only=args.write_only,
        starter_only=args.starter_only,
        stop_at_cycle=args.stop_at_cycle,
        pthickfail=args.pthickfail,
    )

    print(json.dumps(asdict(summary), indent=2, default=str))
    return 0 if summary.starter_rc == 0 and (summary.engine_rc in (None, 0)) else 1


if __name__ == "__main__":
    raise SystemExit(main())
