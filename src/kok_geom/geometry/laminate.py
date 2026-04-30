"""Multi-ply AP-PLY laminate assembly."""

from __future__ import annotations

from dataclasses import dataclass

from kok_geom.config import KokConfig, MM
from kok_geom.geometry.ply import Ply, PlySolid
from kok_geom.geometry.undulation import UndulationRecord, crossing_records
from kok_geom.occ_backend import OCCBackend


@dataclass(frozen=True)
class LaminateSolid:
    plies: tuple[PlySolid, ...]
    undulations: tuple[UndulationRecord, ...]

    def orientation_groups(self) -> list[dict[str, object]]:
        groups: list[dict[str, object]] = []
        for ply in self.plies:
            for tow in ply.tows:
                groups.append(
                    {
                        "name": tow.name,
                        "kind": "straight_tow",
                        "ply_index": ply.ply_index,
                        "fiber_direction_unit_vector": list(tow.orientation),
                        "transverse_in_plane_unit_vector": [
                            -tow.orientation[1],
                            tow.orientation[0],
                            0.0,
                        ],
                        "through_thickness_unit_vector": [0.0, 0.0, 1.0],
                    }
                )
            groups.append(
                {
                    "name": ply.resin.name,
                    "kind": "resin_pocket",
                    "ply_index": ply.ply_index,
                    "isotropic": True,
                }
            )
        return groups

    def undulation_metadata(self) -> list[dict[str, object]]:
        return [
            {
                "name": record.name,
                "kind": "undulation",
                "ply_indices": [record.lower_ply_index, record.upper_ply_index],
                "nominal_angles_deg": [record.lower_angle_deg, record.upper_angle_deg],
                "point_m": list(record.point_m),
                "phi_avg_deg": record.phi_avg_deg,
                "fiber_direction_unit_vector": list(record.fiber_direction_unit_vector),
            }
            for record in self.undulations
        ]


@dataclass(frozen=True)
class Laminate:
    plies: tuple[Ply, ...]

    @classmethod
    def from_config(cls, config: KokConfig) -> "Laminate":
        angles = expand_stacking_sequence(
            config.laydown.fiber_angles_deg,
            n_plies=config.panel.n_plies,
            symmetry=config.panel.symmetry,
        )
        period = len(config.laydown.fiber_angles_deg)
        plies: list[Ply] = []
        for idx, angle in enumerate(angles):
            width_m = config.laydown.tape_width_for_ply_mm(idx) * MM
            period_index = idx // period
            origin_offset_m = period_index * config.laydown.angle_shift_deg * width_m
            plies.append(
                Ply(
                    size_x_m=config.panel.size_x_m,
                    size_y_m=config.panel.size_y_m,
                    angle_deg=angle,
                    tape_width_m=width_m,
                    cured_ply_thickness_m=config.laydown.cured_ply_thickness_m,
                    undulation_ratio=config.laydown.undulation_ratio,
                    tape_spacing=config.laydown.tape_spacing,
                    ply_index=idx + 1,
                    z_bottom_m=idx * config.laydown.cured_ply_thickness_m,
                    origin_offset_m=origin_offset_m,
                )
            )
        return cls(plies=tuple(plies))

    def build_occ(self, backend: OCCBackend | None = None) -> LaminateSolid:
        occ = backend or OCCBackend()
        solids: list[PlySolid] = []
        for ply in self.plies:
            solids.append(ply.build_occ(occ))

        undulations: list[UndulationRecord] = []
        for lower, upper in zip(self.plies, self.plies[1:]):
            undulations.extend(crossing_records(lower, upper, start_index=len(undulations)))
        return LaminateSolid(plies=tuple(solids), undulations=tuple(undulations))


def expand_stacking_sequence(
    base_angles_deg: list[float],
    *,
    n_plies: int,
    symmetry: str = "none",
) -> tuple[float, ...]:
    if n_plies < 1:
        raise ValueError("n_plies must be positive")
    if not base_angles_deg:
        raise ValueError("base_angles_deg cannot be empty")

    def repeat_to_count(count: int) -> list[float]:
        return [base_angles_deg[i % len(base_angles_deg)] for i in range(count)]

    if symmetry == "midplane":
        if n_plies % 2:
            raise ValueError("midplane symmetry requires an even ply count")
        half = repeat_to_count(n_plies // 2)
        return tuple(half + list(reversed(half)))
    if symmetry != "none":
        raise ValueError("symmetry must be 'none' or 'midplane'")
    return tuple(repeat_to_count(n_plies))
