import json
import math

from kok_geom.config import KokConfig, parse_shorthand


def stage11_dict():
    return {
        "panel": {
            "size_x_mm": 25.0,
            "size_y_mm": 25.0,
            "n_plies": 4,
            "symmetry": "none",
        },
        "laydown": {
            "fiber_angles_deg": [0.0, 45.0, -45.0, 90.0],
            "placement_sequence": "1010",
            "angle_shift_deg": 0.0,
            "tape_width_mm": 6.35,
            "cured_ply_thickness_mm": 0.18,
            "undulation_ratio": 0.09,
            "tape_spacing": 1,
        },
        "mesh": {
            "in_plane_target_mm_impact_zone": 0.5,
            "in_plane_target_mm_far_field": 0.5,
            "through_thickness_target_mm": 0.18,
            "graded_zone_radius_mm": 0.0,
            "element_order": 2,
        },
        "output": {
            "msh_path": "build/stage11_block.msh",
            "orientations_json_path": "build/stage11_orientations.json",
            "msh_format": "msh4_ascii",
        },
    }


def test_roundtrip():
    config = KokConfig.model_validate(stage11_dict())
    restored = KokConfig.model_validate(json.loads(config.canonical_json()))
    assert restored.model_dump(mode="json") == config.model_dump(mode="json")


def test_uofsc_shorthand():
    data = parse_shorthand("[0,45,-45,90][1010][0][6.35]")
    config = KokConfig.model_validate(data)
    assert config.panel.n_plies == 4
    assert config.laydown.fiber_angles_deg == [0.0, 45.0, -45.0, 90.0]
    assert config.laydown.placement_sequence == "1010"
    assert config.laydown.tape_spacing == 1
    assert math.isclose(config.laydown.tape_width_mm, 6.35)


def test_nagelsmit_shorthand():
    data = parse_shorthand("[(45/-45)1x0.5/(90/0)1x0.5]_3S")
    config = KokConfig.model_validate(data)
    assert config.panel.symmetry == "midplane"
    assert config.panel.n_plies == 24
    assert config.laydown.fiber_angles_deg == [45.0, -45.0, 90.0, 0.0]
    assert config.laydown.tape_spacing == 1
    assert math.isclose(config.laydown.tape_width_mm, 12.7)
