"""Undulation bookkeeping for crossing tow centerlines."""

from __future__ import annotations

import math
from dataclasses import dataclass

from kok_geom.geometry.ply import Ply
from kok_geom.geometry.tow import UndulationProfile


@dataclass(frozen=True)
class UndulationRecord:
    name: str
    lower_ply_index: int
    upper_ply_index: int
    point_m: tuple[float, float]
    lower_angle_deg: float
    upper_angle_deg: float
    phi_avg_deg: float
    fiber_direction_unit_vector: tuple[float, float, float]


def _det2(a: tuple[float, float], b: tuple[float, float]) -> float:
    return a[0] * b[1] - a[1] * b[0]


def _sub2(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
    return (a[0] - b[0], a[1] - b[1])


def crossing_records(lower: Ply, upper: Ply, *, start_index: int = 0) -> tuple[UndulationRecord, ...]:
    """Return centerline crossings between adjacent plies inside the panel box."""

    u1 = lower.direction
    u2 = upper.direction
    denom = _det2(u1, u2)
    if math.isclose(denom, 0.0, abs_tol=1.0e-12):
        return ()

    n1 = lower.transverse
    n2 = upper.transverse
    profile = UndulationProfile(
        cured_ply_thickness_m=upper.cured_ply_thickness_m,
        undulation_ratio=upper.undulation_ratio,
    )
    phi = profile.phi_avg_rad
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)
    ux, uy = upper.direction
    vec = (ux * cos_phi, uy * cos_phi, sin_phi)
    norm = math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])
    fiber = (vec[0] / norm, vec[1] / norm, vec[2] / norm)

    records: list[UndulationRecord] = []
    tol = 1.0e-10
    x_half = upper.size_x_m / 2.0 + tol
    y_half = upper.size_y_m / 2.0 + tol
    for lower_offset in lower.center_offsets_m():
        p1 = (n1[0] * lower_offset, n1[1] * lower_offset)
        for upper_offset in upper.center_offsets_m():
            p2 = (n2[0] * upper_offset, n2[1] * upper_offset)
            rhs = _sub2(p2, p1)
            along_lower = _det2(rhs, u2) / denom
            x = p1[0] + along_lower * u1[0]
            y = p1[1] + along_lower * u1[1]
            if -x_half <= x <= x_half and -y_half <= y <= y_half:
                idx = start_index + len(records)
                records.append(
                    UndulationRecord(
                        name=f"UNDUL_PLY_{lower.ply_index}_PLY_{upper.ply_index}_TAG_{idx}",
                        lower_ply_index=lower.ply_index,
                        upper_ply_index=upper.ply_index,
                        point_m=(x, y),
                        lower_angle_deg=lower.angle_deg,
                        upper_angle_deg=upper.angle_deg,
                        phi_avg_deg=math.degrees(phi),
                        fiber_direction_unit_vector=fiber,
                    )
                )
    return tuple(records)
