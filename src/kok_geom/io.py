"""M1 mesh and sidecar export."""

from __future__ import annotations

import json
from pathlib import Path

import meshio

from kok_geom.config import KokConfig, MM
from kok_geom.geometry.tow import Tow
from kok_geom.mesh import configure_msh_format, generate_volume_mesh
from kok_geom.occ_backend import GmshSession, OCCBackend


def write_orientations(
    path: str | Path,
    *,
    config: KokConfig,
    tow: Tow,
    physical_name: str,
) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "schema_version": 1,
        "panel_config_hash": config.config_hash(),
        "groups": [
            {
                "name": physical_name,
                "kind": "straight_tow",
                "ply_index": 1,
                "nominal_angle_deg": tow.angle_deg,
                "fiber_direction_unit_vector": list(tow.orientation),
                "transverse_in_plane_unit_vector": list(tow.transverse_in_plane),
                "through_thickness_unit_vector": [0.0, 0.0, 1.0],
            }
        ],
    }
    out.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    return out


def convert_msh_to_inp(msh_path: str | Path, inp_path: str | Path) -> Path:
    mesh = meshio.read(msh_path)
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

    angle = config.laydown.fiber_angles_deg[0]
    tow = Tow(
        length_m=config.panel.size_x_m,
        angle_deg=angle,
        tape_width_m=config.laydown.tape_width_for_ply_mm(0) * MM,
        cured_ply_thickness_m=config.laydown.cured_ply_thickness_m,
        undulation_ratio=config.laydown.undulation_ratio,
        center_m=(0.0, 0.0, 0.5 * config.laydown.cured_ply_thickness_m),
        name="TOW_PLY_1_TAG_0",
    )
    msh_path = Path(config.output.msh_path)
    orientations_path = Path(config.output.orientations_json_path)

    with GmshSession("kok_geom_m1"):
        backend = OCCBackend()
        tow.build_occ(backend, add_physical_group=True)
        configure_msh_format(config.output.msh_format)
        generate_volume_mesh(
            target_size_m=config.mesh.target_size_m,
            element_order=config.mesh.element_order,
        )
        backend.write(msh_path)

    write_orientations(
        orientations_path,
        config=config,
        tow=tow,
        physical_name=tow.name,
    )
    if config.output.inp_path:
        convert_msh_to_inp(msh_path, config.output.inp_path)
    return msh_path, orientations_path
