"""M1 mesh and sidecar export."""

from __future__ import annotations

import json
from pathlib import Path

import meshio

from kok_geom.config import KokConfig, MM
from kok_geom.geometry.laminate import Laminate, LaminateSolid
from kok_geom.geometry.ply import Ply, PlySolid
from kok_geom.geometry.tow import Tow
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
    data = {
        "schema_version": 1,
        "panel_config_hash": config.config_hash(),
        "groups": orientation_groups,
    }
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
        else:
            config.output.msh_path = str(out_path / Path(config.output.msh_path).name)
            config.output.orientations_json_path = str(out_path / Path(config.output.orientations_json_path).name)

    laminate = Laminate.from_config(config)
    msh_path = Path(config.output.msh_path)
    orientations_path = Path(config.output.orientations_json_path)

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
