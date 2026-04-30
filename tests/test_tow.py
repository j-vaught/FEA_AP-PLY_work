import json
import math

import gmsh
import meshio
import numpy as np

from kok_geom.config import KokConfig
from kok_geom.geometry.tow import Tow, UndulationProfile
from kok_geom.io import generate_mesh
from kok_geom.mesh import configure_msh_format, generate_volume_mesh
from kok_geom.occ_backend import GmshSession, OCCBackend


def test_single_tow_orientation_unit_vector():
    tow = Tow.from_mm(length_mm=25.0, angle_deg=45.0)
    vx, vy, vz = tow.orientation
    assert math.isclose(math.sqrt(vx * vx + vy * vy + vz * vz), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(vx, math.sqrt(0.5), rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(vy, math.sqrt(0.5), rel_tol=0.0, abs_tol=1.0e-12)
    assert vz == 0.0


def test_single_tow_volume():
    tow = Tow.from_mm(length_mm=25.0, tape_width_mm=1.0, cured_ply_thickness_mm=0.18)
    with GmshSession("test_single_tow_volume"):
        backend = OCCBackend()
        solid = tow.build_occ(backend)
        backend.synchronize()
        assert math.isclose(
            backend.volume(solid.dim_tag),
            tow.expected_volume_m3,
            rel_tol=1.0e-3,
            abs_tol=0.0,
        )


def test_phi_avg_integral():
    profile = UndulationProfile(cured_ply_thickness_m=0.18e-3, undulation_ratio=0.09)
    x = np.linspace(0.0, profile.half_length_m, 200_001)
    expected = float(np.trapezoid(profile.phi_rad(x), x) / profile.half_length_m)
    assert math.isclose(profile.phi_avg_rad, expected, rel_tol=0.0, abs_tol=1.0e-9)
    assert 0.0 < profile.phi_avg_deg < 10.0


def test_m1_meshio_roundtrip(tmp_path):
    config = KokConfig.model_validate(
        {
            "panel": {"size_x_mm": 4.0, "size_y_mm": 2.0, "n_plies": 1, "symmetry": "none"},
            "laydown": {
                "fiber_angles_deg": [0.0],
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
                "msh_path": str(tmp_path / "tow.msh"),
                "orientations_json_path": str(tmp_path / "orientations.json"),
                "msh_format": "msh4_ascii",
            },
        }
    )
    config_path = tmp_path / "kok_config.json"
    config.to_file(config_path)

    msh_path, orientations_path = generate_mesh(config_path)
    mesh = meshio.read(msh_path)
    assert any(block.type == "tetra10" for block in mesh.cells)
    assert not any(block.type == "tetra" for block in mesh.cells)

    inp_path = tmp_path / "tow.inp"
    meshio.write(inp_path, mesh, file_format="abaqus")
    reread = meshio.read(inp_path)
    assert sum(len(block.data) for block in reread.cells) > 0

    orientations = json.loads(orientations_path.read_text())
    assert orientations["groups"][0]["name"] == "TOW_PLY_1_TAG_0"
    assert orientations["groups"][0]["fiber_direction_unit_vector"] == [1.0, 0.0, 0.0]
