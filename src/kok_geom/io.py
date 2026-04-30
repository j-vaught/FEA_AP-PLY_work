"""M1 mesh and sidecar export."""

from __future__ import annotations

import json
import math
from pathlib import Path

import meshio
import numpy as np

from kok_geom.config import KokConfig, MM
from kok_geom.geometry.laminate import Laminate, LaminateSolid
from kok_geom.geometry.ply import Ply, PlySolid
from kok_geom.geometry.tow import Tow, UndulationProfile
from kok_geom.mesh import configure_msh_format, generate_volume_mesh
from kok_geom.occ_backend import GmshSession, OCCBackend


def write_orientations(
    path: str | Path,
    *,
    config: KokConfig,
    tow: Tow | None = None,
    physical_name: str | None = None,
    ply: PlySolid | None = None,
    laminate: LaminateSolid | None = None,
    groups: list[dict[str, object]] | None = None,
) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    orientation_groups: list[dict[str, object]] = []
    if groups is not None:
        orientation_groups = groups
    elif laminate is not None:
        orientation_groups = laminate.orientation_groups()
    elif ply is not None:
        for ply_tow in ply.tows:
            orientation_groups.append(
                {
                    "name": ply_tow.name,
                    "kind": "straight_tow",
                    "ply_index": ply.ply_index,
                    "nominal_angle_deg": config.laydown.fiber_angles_deg[0],
                    "fiber_direction_unit_vector": list(ply_tow.orientation),
                    "transverse_in_plane_unit_vector": [
                        -ply_tow.orientation[1],
                        ply_tow.orientation[0],
                        0.0,
                    ],
                    "through_thickness_unit_vector": [0.0, 0.0, 1.0],
                }
            )
        orientation_groups.append(
            {
                "name": ply.resin.name,
                "kind": "resin_pocket",
                "ply_index": ply.ply_index,
                "isotropic": True,
            }
        )
    elif tow is not None and physical_name is not None:
        orientation_groups.append(
            {
                "name": physical_name,
                "kind": "straight_tow",
                "ply_index": 1,
                "nominal_angle_deg": tow.angle_deg,
                "fiber_direction_unit_vector": list(tow.orientation),
                "transverse_in_plane_unit_vector": list(tow.transverse_in_plane),
                "through_thickness_unit_vector": [0.0, 0.0, 1.0],
            }
        )
    else:
        raise ValueError("write_orientations requires either ply or tow+physical_name")
    data: dict[str, object] = {
        "schema_version": 1,
        "panel_config_hash": config.config_hash(),
        "groups": orientation_groups,
    }
    if laminate is not None:
        data["undulations"] = laminate.undulation_metadata()
    out.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    return out


def convert_msh_to_inp(msh_path: str | Path, inp_path: str | Path) -> Path:
    mesh = meshio.read(msh_path)
    # meshio exposes signed GMSH bounding-entity ids as a cell set; Abaqus INP
    # has no matching concept and its reader treats the negative ids as bad set
    # references on round-trip. Physical-group ELSETs are preserved without it.
    mesh.cell_sets.pop("gmsh:bounding_entities", None)
    out = Path(inp_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(out, mesh, file_format="abaqus")
    return out


def generate_mesh(config_path: str | Path, out: str | Path | None = None) -> tuple[Path, Path]:
    """Build the M1 single-tow mesh and orientation sidecar."""

    config = KokConfig.from_file(config_path)
    if out is not None:
        out_path = Path(out)
        if out_path.suffix:
            config.output.msh_path = str(out_path)
            config.output.orientations_json_path = str(out_path.with_name("orientations.json"))
            if config.output.inp_path is None and out_path.suffix.lower() == ".msh":
                config.output.inp_path = str(out_path.with_suffix(".inp"))
        else:
            config.output.msh_path = str(out_path / Path(config.output.msh_path).name)
            config.output.orientations_json_path = str(out_path / Path(config.output.orientations_json_path).name)
            if config.output.inp_path is None:
                config.output.inp_path = str(out_path / Path(config.output.msh_path).with_suffix(".inp").name)

    msh_path = Path(config.output.msh_path)
    orientations_path = Path(config.output.orientations_json_path)

    if config.panel.n_plies > 1:
        _generate_structured_laminate_mesh(config, msh_path, orientations_path)
        if config.output.inp_path:
            convert_msh_to_inp(msh_path, config.output.inp_path)
        return msh_path, orientations_path

    laminate = Laminate.from_config(config)

    with GmshSession("kok_geom_m1"):
        backend = OCCBackend()
        laminate_solid = laminate.build_occ(backend)
        configure_msh_format(config.output.msh_format)
        generate_volume_mesh(
            target_size_m=config.mesh.target_size_m,
            element_order=config.mesh.element_order,
        )
        backend.write(msh_path)

    write_orientations(
        orientations_path,
        config=config,
        laminate=laminate_solid,
    )
    if config.output.inp_path:
        convert_msh_to_inp(msh_path, config.output.inp_path)
    return msh_path, orientations_path


def _generate_structured_laminate_mesh(
    config: KokConfig,
    msh_path: Path,
    orientations_path: Path,
) -> None:
    """Write a clean TETRA10 mesh classified by the AP-PLY analytic geometry."""

    laminate = Laminate.from_config(config)
    events_by_tow, records = laminate._undulation_events()
    size_x = config.panel.size_x_m
    size_y = config.panel.size_y_m
    ply_t = config.laydown.cured_ply_thickness_m
    n_plies = config.panel.n_plies
    target = config.mesh.target_size_m
    nx = max(2, int(math.ceil(size_x / target)))
    ny = max(2, int(math.ceil(size_y / target)))
    nz = max(n_plies, int(math.ceil((n_plies * ply_t) / target)))
    if nz % n_plies:
        nz += n_plies - (nz % n_plies)

    xs = np.linspace(-0.5 * size_x, 0.5 * size_x, nx + 1)
    ys = np.linspace(-0.5 * size_y, 0.5 * size_y, ny + 1)
    zs = np.linspace(0.0, n_plies * ply_t, nz + 1)
    points: list[tuple[float, float, float]] = [
        (float(x), float(y), float(z)) for z in zs for y in ys for x in xs
    ]

    def vid(i: int, j: int, k: int) -> int:
        return k * (ny + 1) * (nx + 1) + j * (nx + 1) + i

    midpoint_index: dict[tuple[int, int], int] = {}

    def mid(a: int, b: int) -> int:
        key = (a, b) if a < b else (b, a)
        if key not in midpoint_index:
            pa = points[a]
            pb = points[b]
            midpoint_index[key] = len(points)
            points.append(
                (
                    0.5 * (pa[0] + pb[0]),
                    0.5 * (pa[1] + pb[1]),
                    0.5 * (pa[2] + pb[2]),
                )
            )
        return midpoint_index[key]

    group_ids: dict[str, int] = {}
    orientation_groups: dict[str, dict[str, object]] = {}

    def group_id(name: str, data: dict[str, object]) -> int:
        if name not in group_ids:
            group_ids[name] = len(group_ids) + 1
            orientation_groups[name] = data
        return group_ids[name]

    ply_offsets = [ply.center_offsets_m() for ply in laminate.plies]
    cells: list[list[int]] = []
    physical: list[int] = []
    tet_pattern = (
        (0, 1, 3, 7),
        (0, 3, 2, 7),
        (0, 2, 6, 7),
        (0, 6, 4, 7),
        (0, 4, 5, 7),
        (0, 5, 1, 7),
    )
    tet_edges = ((0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3))

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                hex_vertices = (
                    vid(i, j, k),
                    vid(i + 1, j, k),
                    vid(i, j + 1, k),
                    vid(i + 1, j + 1, k),
                    vid(i, j, k + 1),
                    vid(i + 1, j, k + 1),
                    vid(i, j + 1, k + 1),
                    vid(i + 1, j + 1, k + 1),
                )
                for tet in tet_pattern:
                    corners = [hex_vertices[idx] for idx in tet]
                    centroid = tuple(
                        sum(points[node][axis] for node in corners) / 4.0 for axis in range(3)
                    )
                    name, data = _classify_structured_cell(
                        laminate,
                        ply_offsets,
                        events_by_tow,
                        centroid,
                    )
                    cells.append(corners + [mid(corners[a], corners[b]) for a, b in tet_edges])
                    physical.append(group_id(name, data))

    field_data = {name: np.array([phys_id, 3], dtype=int) for name, phys_id in group_ids.items()}
    mesh = meshio.Mesh(
        points=np.asarray(points, dtype=float),
        cells=[("tetra10", np.asarray(cells, dtype=int))],
        cell_data={"gmsh:physical": [np.asarray(physical, dtype=int)]},
        field_data=field_data,
    )
    msh_path.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(msh_path, mesh, file_format="gmsh22", binary=False)
    write_orientations(
        orientations_path,
        config=config,
        groups=list(orientation_groups.values()),
    )
    data = json.loads(orientations_path.read_text(encoding="utf-8"))
    data["undulations"] = [
        {
            "name": record.name,
            "kind": "undulation",
            "ply_indices": [record.lower_ply_index, record.upper_ply_index],
            "nominal_angles_deg": [record.lower_angle_deg, record.upper_angle_deg],
            "point_m": list(record.point_m),
            "phi_avg_deg": record.phi_avg_deg,
            "fiber_direction_unit_vector": list(record.fiber_direction_unit_vector),
        }
        for record in records
    ]
    orientations_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _classify_structured_cell(
    laminate: Laminate,
    ply_offsets: list[tuple[float, ...]],
    events_by_tow: dict[tuple[int, int], list[object]],
    centroid: tuple[float, float, float],
) -> tuple[str, dict[str, object]]:
    x, y, z = centroid
    ply_t = laminate.plies[0].cured_ply_thickness_m
    ply_idx = min(len(laminate.plies) - 1, max(0, int(z / ply_t)))
    ply = laminate.plies[ply_idx]
    ux, uy = ply.direction
    nx, ny = ply.transverse
    d = x * nx + y * ny
    offsets = ply_offsets[ply_idx]
    tow_idx = min(range(len(offsets)), key=lambda idx: abs(offsets[idx] - d))
    offset = offsets[tow_idx]
    if abs(d - offset) > 0.5 * ply.tape_width_m:
        name = f"RESIN_PLY_{ply.ply_index}"
        return name, {"name": name, "kind": "resin_pocket", "ply_index": ply.ply_index, "isotropic": True}

    s = x * ux + y * uy
    profile = UndulationProfile(
        cured_ply_thickness_m=ply.cured_ply_thickness_m,
        undulation_ratio=ply.undulation_ratio,
    )
    for event in events_by_tow.get((ply.ply_index, tow_idx), ()):
        center = float(event.center_s_m)
        lu = profile.undulation_length_m
        if center - lu <= s < center:
            return _structured_undulation_group(ply, event, sign=1.0, name=event.ramp_up_name)
        if center <= s <= center + lu:
            return _structured_undulation_group(ply, event, sign=-1.0, name=event.ramp_down_name)

    name = f"TOW_PLY_{ply.ply_index}_TAG_{tow_idx}"
    return name, {
        "name": name,
        "kind": "straight_tow",
        "ply_index": ply.ply_index,
        "nominal_angle_deg": ply.angle_deg,
        "fiber_direction_unit_vector": [ux, uy, 0.0],
        "transverse_in_plane_unit_vector": [nx, ny, 0.0],
        "through_thickness_unit_vector": [0.0, 0.0, 1.0],
    }


def _structured_undulation_group(ply: Ply, event: object, *, sign: float, name: str) -> tuple[str, dict[str, object]]:
    profile = UndulationProfile(
        cured_ply_thickness_m=ply.cured_ply_thickness_m,
        undulation_ratio=ply.undulation_ratio,
    )
    phi = profile.phi_avg_rad
    ux, uy = ply.direction
    vx, vy = ply.transverse
    fiber = [ux * math.cos(phi), uy * math.cos(phi), sign * math.sin(phi)]
    transverse = [vx, vy, 0.0]
    through = [
        fiber[1] * transverse[2] - fiber[2] * transverse[1],
        fiber[2] * transverse[0] - fiber[0] * transverse[2],
        fiber[0] * transverse[1] - fiber[1] * transverse[0],
    ]
    norm = math.sqrt(sum(value * value for value in through))
    through = [value / norm for value in through]
    return name, {
        "name": name,
        "kind": "undulation",
        "ply_index": ply.ply_index,
        "ply_indices": [event.record.lower_ply_index, event.record.upper_ply_index],
        "nominal_angle_deg": ply.angle_deg,
        "nominal_angles_deg": [event.record.lower_angle_deg, event.record.upper_angle_deg],
        "point_m": list(event.record.point_m),
        "phi_avg_deg": math.degrees(phi),
        "fiber_direction_unit_vector": fiber,
        "transverse_in_plane_unit_vector": transverse,
        "through_thickness_unit_vector": through,
    }
