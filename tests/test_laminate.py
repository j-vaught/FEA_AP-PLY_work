import json
import math

from kok_geom.config import KokConfig
from kok_geom.geometry.laminate import Laminate, expand_stacking_sequence
from kok_geom.io import generate_mesh
from kok_geom.occ_backend import GmshSession, OCCBackend


def test_expand_stacking_sequence_midplane():
    assert expand_stacking_sequence([0.0, 45.0], n_plies=5) == (0.0, 45.0, 0.0, 45.0, 0.0)
    assert expand_stacking_sequence([0.0, 45.0], n_plies=4, symmetry="midplane") == (
        0.0,
        45.0,
        45.0,
        0.0,
    )


def test_angle_shift_translates_repeated_period_by_one_tape_width():
    config = KokConfig.model_validate(
        {
            "panel": {"size_x_mm": 40.0, "size_y_mm": 40.0, "n_plies": 8, "symmetry": "none"},
            "laydown": {
                "fiber_angles_deg": [0.0, 45.0, -45.0, 90.0],
                "placement_sequence": "1010",
                "angle_shift_deg": 1.0,
                "tape_width_mm": 6.35,
                "cured_ply_thickness_mm": 0.18,
                "undulation_ratio": 0.09,
                "tape_spacing": 1,
            },
        }
    )
    laminate = Laminate.from_config(config)
    ply1 = laminate.plies[0]
    ply5 = laminate.plies[4]
    assert ply1.angle_deg == ply5.angle_deg == 0.0
    deltas = [b - a for a, b in zip(ply1.center_offsets_m(), ply5.center_offsets_m())]
    assert deltas
    assert all(math.isclose(delta, 6.35e-3, rel_tol=0.0, abs_tol=1.0e-12) for delta in deltas)


def test_laminate_orientation_groups_include_each_physical_tow_and_resin():
    config = KokConfig.model_validate(
        {
            "panel": {"size_x_mm": 8.0, "size_y_mm": 8.0, "n_plies": 2, "symmetry": "none"},
            "laydown": {
                "fiber_angles_deg": [0.0, 90.0],
                "placement_sequence": "10",
                "tape_width_mm": 1.0,
                "cured_ply_thickness_mm": 0.18,
                "undulation_ratio": 0.09,
                "tape_spacing": 1,
            },
        }
    )
    with GmshSession("test_laminate_orientation_groups"):
        solid = Laminate.from_config(config).build_occ(OCCBackend())
        groups = solid.orientation_groups()
        names = {group["name"] for group in groups}
        assert "RESIN_PLY_1" in names
        assert "RESIN_PLY_2" in names
        assert any(name.startswith("TOW_PLY_1_TAG_") for name in names)
        assert any(name.startswith("TOW_PLY_2_TAG_") for name in names)
        assert any(name.startswith("UNDUL_PLY_1_PLY_2_TAG_") for name in names)
        for group in groups:
            vector = group.get("fiber_direction_unit_vector")
            if vector is not None:
                norm = math.sqrt(sum(float(v) * float(v) for v in vector))
                assert math.isclose(norm, 1.0, rel_tol=0.0, abs_tol=1.0e-12)


def test_cli_mult_ply_orientation_sidecar(tmp_path):
    config = KokConfig.model_validate(
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
                "msh_path": str(tmp_path / "laminate.msh"),
                "orientations_json_path": str(tmp_path / "orientations.json"),
                "msh_format": "msh4_ascii",
            },
        }
    )
    config_path = tmp_path / "kok_config.json"
    config.to_file(config_path)
    _, orientations_path = generate_mesh(config_path)
    sidecar = json.loads(orientations_path.read_text())
    kinds = {group["kind"] for group in sidecar["groups"]}
    assert {"straight_tow", "resin_pocket", "undulation"}.issubset(kinds)
