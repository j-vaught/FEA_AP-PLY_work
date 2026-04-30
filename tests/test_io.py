import json
import math

import meshio

from kok_geom.config import KokConfig
from kok_geom.io import convert_msh_to_inp, generate_mesh


def _small_laminate_config(tmp_path):
    return KokConfig.model_validate(
        {
            "panel": {"size_x_mm": 4.0, "size_y_mm": 4.0, "n_plies": 2, "symmetry": "none"},
            "laydown": {
                "fiber_angles_deg": [0.0, 90.0],
                "placement_sequence": "10",
                "tape_width_mm": 1.0,
                "cured_ply_thickness_mm": 0.18,
                "undulation_ratio": 0.09,
                "tape_spacing": 1,
            },
            "mesh": {
                "in_plane_target_mm_impact_zone": 1.0,
                "in_plane_target_mm_far_field": 1.0,
                "through_thickness_target_mm": 0.18,
                "element_order": 2,
            },
            "output": {
                "msh_path": str(tmp_path / "panel.msh"),
                "orientations_json_path": str(tmp_path / "orientations.json"),
                "msh_format": "msh4_ascii",
            },
        }
    )


def test_orientation_sidecar_matches_msh_physical_groups(tmp_path):
    config = _small_laminate_config(tmp_path)
    config_path = tmp_path / "kok_config.json"
    config.to_file(config_path)
    msh_path, orientations_path = generate_mesh(config_path)

    mesh = meshio.read(msh_path)
    sidecar = json.loads(orientations_path.read_text())
    physical_names = set(mesh.field_data)
    sidecar_names = {group["name"] for group in sidecar["groups"]}
    assert sidecar_names == physical_names
    assert "undulations" in sidecar


def test_orientation_vectors_are_unit_norm(tmp_path):
    config = _small_laminate_config(tmp_path)
    config_path = tmp_path / "kok_config.json"
    config.to_file(config_path)
    _, orientations_path = generate_mesh(config_path)

    sidecar = json.loads(orientations_path.read_text())
    for section in ("groups", "undulations"):
        for group in sidecar.get(section, []):
            vector = group.get("fiber_direction_unit_vector")
            if vector is None:
                continue
            norm = math.sqrt(sum(float(v) * float(v) for v in vector))
            assert math.isclose(norm, 1.0, rel_tol=0.0, abs_tol=1.0e-12)


def test_msh_to_inp_roundtrip(tmp_path):
    config = _small_laminate_config(tmp_path)
    config_path = tmp_path / "kok_config.json"
    config.to_file(config_path)
    msh_path, _ = generate_mesh(config_path)

    inp_path = convert_msh_to_inp(msh_path, tmp_path / "panel.inp")
    reread = meshio.read(inp_path)
    assert sum(len(block.data) for block in reread.cells) > 0
