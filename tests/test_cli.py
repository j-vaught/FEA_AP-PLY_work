import json
import subprocess
import sys

import meshio

from kok_geom.config import KokConfig


def test_cli_writes_msh_and_orientation_sidecar(tmp_path):
    config = KokConfig.model_validate(
        {
            "panel": {"size_x_mm": 4.0, "size_y_mm": 4.0, "n_plies": 1, "symmetry": "none"},
            "laydown": {
                "fiber_angles_deg": [30.0],
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
                "msh_path": str(tmp_path / "ignored.msh"),
                "orientations_json_path": str(tmp_path / "ignored_orientations.json"),
                "msh_format": "msh4_ascii",
            },
        }
    )
    config_path = tmp_path / "kok_config.json"
    config.to_file(config_path)

    out_path = tmp_path / "cli_panel.msh"
    proc = subprocess.run(
        [sys.executable, "-m", "kok_geom", "--config", str(config_path), "--out", str(out_path)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout
    assert out_path.exists()
    sidecar = out_path.with_name("orientations.json")
    assert sidecar.exists()
    assert any(block.type == "tetra10" for block in meshio.read(out_path).cells)
    assert json.loads(sidecar.read_text())["groups"]
