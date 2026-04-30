import math

from kok_geom.geometry.ply import Ply
from kok_geom.occ_backend import GmshSession, OCCBackend


def _volume(body):
    return body.dim_tag if hasattr(body, "dim_tag") else body


def test_interlace_pattern_s1_s2_s3():
    for spacing in (1, 2, 3):
        pitch_mm = (spacing + 1) * 2.0
        ply = Ply.from_mm(
            size_x_mm=20.0,
            size_y_mm=4.0 * pitch_mm,
            tape_width_mm=2.0,
            cured_ply_thickness_mm=0.2,
            tape_spacing=spacing,
            tow_coverage_fraction=0.85,
        )
        assert math.isclose(ply.active_pitch_m, pitch_mm * 1.0e-3, rel_tol=0.0, abs_tol=1.0e-12)
        assert math.isclose(ply.lane_pitch_m, 2.0e-3 / 0.85, rel_tol=0.0, abs_tol=1.0e-12)
        assert ply.expected_tow_count() >= 8
        assert len(ply.center_offsets_m()) == ply.expected_tow_count()

        with GmshSession(f"test_ply_spacing_{spacing}"):
            solid = ply.build_occ(OCCBackend())
            assert len(solid.tows) == ply.expected_tow_count()
            assert all(tow.name.startswith("TOW_PLY_1_TAG_") for tow in solid.tows)
            assert solid.resin.name == "RESIN_PLY_1"


def test_single_ply_resin_volume_closure():
    ply = Ply.from_mm(
        size_x_mm=20.0,
        size_y_mm=20.0,
        tape_width_mm=2.0,
        cured_ply_thickness_mm=0.2,
        tape_spacing=1,
        tow_coverage_fraction=0.85,
    )
    with GmshSession("test_single_ply_resin_volume_closure"):
        backend = OCCBackend()
        solid = ply.build_occ(backend)
        tow_volume = sum(backend.volume(tow.dim_tag) for tow in solid.tows)
        resin_volume = sum(backend.volume(tag) for tag in solid.resin.dim_tags)
        assert math.isclose(
            tow_volume + resin_volume,
            ply.nominal_volume_m3,
            rel_tol=1.0e-8,
            abs_tol=1.0e-15,
        )
        assert tow_volume > 0.0
        assert resin_volume > 0.0
