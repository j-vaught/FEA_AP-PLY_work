"""Multi-ply AP-PLY laminate assembly."""

from __future__ import annotations

import math
from dataclasses import dataclass

import gmsh

from kok_geom.config import KokConfig, MM
from kok_geom.geometry.ply import Ply, PlySolid
from kok_geom.geometry.undulation import UndulationRecord, crossing_records
from kok_geom.geometry.tow import UndulationProfile
from kok_geom.occ_backend import DimTag, OCCBackend


@dataclass(frozen=True)
class PhysicalVolumeGroup:
    name: str
    kind: str
    dim_tags: tuple[DimTag, ...]
    ply_index: int | None = None
    ply_indices: tuple[int, int] | None = None
    nominal_angle_deg: float | None = None
    nominal_angles_deg: tuple[float, float] | None = None
    fiber_direction_unit_vector: tuple[float, float, float] | None = None
    transverse_in_plane_unit_vector: tuple[float, float, float] | None = None
    through_thickness_unit_vector: tuple[float, float, float] | None = None
    phi_avg_deg: float | None = None
    point_m: tuple[float, float] | None = None
    is_isotropic: bool = False


@dataclass(frozen=True)
class UndulationEvent:
    record: UndulationRecord
    ply_index: int
    tow_index: int
    tow_offset_m: float
    center_s_m: float
    ramp_up_name: str
    ramp_down_name: str


@dataclass(frozen=True)
class LaminateSolid:
    plies: tuple[PlySolid, ...]
    undulations: tuple[UndulationRecord, ...]
    groups: tuple[PhysicalVolumeGroup, ...] = ()

    def orientation_groups(self) -> list[dict[str, object]]:
        if self.groups:
            return [_group_to_orientation_dict(group) for group in self.groups]

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
                    tow_coverage_fraction=config.laydown.tow_coverage_fraction,
                    ply_index=idx + 1,
                    z_bottom_m=idx * config.laydown.cured_ply_thickness_m,
                    origin_offset_m=origin_offset_m,
                )
            )
        return cls(plies=tuple(plies))

    def build_occ(self, backend: OCCBackend | None = None) -> LaminateSolid:
        occ = backend or OCCBackend()
        events_by_tow, records = self._undulation_events()
        groups: list[PhysicalVolumeGroup] = []

        for ply in self.plies:
            ply_clip = occ.add_box(
                x0=-ply.size_x_m / 2.0,
                y0=-ply.size_y_m / 2.0,
                z0=ply.z_bottom_m,
                dx=ply.size_x_m,
                dy=ply.size_y_m,
                dz=ply.cured_ply_thickness_m,
            )
            first_group = len(groups)
            self._build_ply_tows(occ, ply, ply_clip, events_by_tow, groups)
            occupied = [tag.tuple for group in groups[first_group:] for tag in group.dim_tags]
            if occupied:
                resin_out, _ = gmsh.model.occ.cut([ply_clip.tuple], occupied, removeObject=True, removeTool=False)
                resin_tags = tuple(DimTag(dim, tag) for dim, tag in resin_out if dim == 3)
                if resin_tags:
                    groups.append(
                        PhysicalVolumeGroup(
                            name=f"RESIN_PLY_{ply.ply_index}",
                            kind="resin_pocket",
                            dim_tags=resin_tags,
                            ply_index=ply.ply_index,
                            is_isotropic=True,
                        )
                    )

        groups = _fragment_and_remap(groups)
        gmsh.model.occ.synchronize()
        for group in groups:
            if not group.dim_tags:
                continue
            phys = gmsh.model.addPhysicalGroup(3, [tag.tag for tag in group.dim_tags])
            gmsh.model.setPhysicalName(3, phys, group.name)

        return LaminateSolid(plies=(), undulations=records, groups=tuple(groups))

    def _undulation_events(self) -> tuple[dict[tuple[int, int], list[UndulationEvent]], tuple[UndulationRecord, ...]]:
        events_by_tow: dict[tuple[int, int], list[UndulationEvent]] = {}
        records: list[UndulationRecord] = []

        for lower, upper in zip(self.plies, self.plies[1:]):
            u1 = lower.direction
            u2 = upper.direction
            denom = _det2(u1, u2)
            if math.isclose(denom, 0.0, abs_tol=1.0e-12):
                continue
            n1 = lower.transverse
            n2 = upper.transverse
            profile = UndulationProfile(
                cured_ply_thickness_m=lower.cured_ply_thickness_m,
                undulation_ratio=lower.undulation_ratio,
            )
            phi = profile.phi_avg_rad
            x_half = lower.size_x_m / 2.0 + 1.0e-10
            y_half = lower.size_y_m / 2.0 + 1.0e-10
            global_crossing_index = 0
            for lower_idx, lower_offset in enumerate(lower.center_offsets_m()):
                p1 = (n1[0] * lower_offset, n1[1] * lower_offset)
                for upper_offset in upper.center_offsets_m():
                    p2 = (n2[0] * upper_offset, n2[1] * upper_offset)
                    rhs = _sub2(p2, p1)
                    along_lower = _det2(rhs, u2) / denom
                    x = p1[0] + along_lower * u1[0]
                    y = p1[1] + along_lower * u1[1]
                    if not (-x_half <= x <= x_half and -y_half <= y <= y_half):
                        continue
                    if global_crossing_index % lower.tape_spacing != 0:
                        global_crossing_index += 1
                        continue
                    tow_key = (lower.ply_index, lower_idx)
                    existing = events_by_tow.get(tow_key, [])
                    if any(
                        abs(event.center_s_m - along_lower) < 2.05 * profile.undulation_length_m
                        for event in existing
                    ):
                        global_crossing_index += 1
                        continue

                    record_index = len(records)
                    up_name = f"UNDUL_PLY_{lower.ply_index}_PLY_{upper.ply_index}_TAG_{record_index}_UP"
                    down_name = f"UNDUL_PLY_{lower.ply_index}_PLY_{upper.ply_index}_TAG_{record_index}_DOWN"
                    ux, uy = lower.direction
                    cos_phi = math.cos(phi)
                    sin_phi = math.sin(phi)
                    fiber = _normalize3((ux * cos_phi, uy * cos_phi, sin_phi))
                    record = UndulationRecord(
                        name=f"UNDUL_PLY_{lower.ply_index}_PLY_{upper.ply_index}_TAG_{record_index}",
                        lower_ply_index=lower.ply_index,
                        upper_ply_index=upper.ply_index,
                        point_m=(x, y),
                        lower_angle_deg=lower.angle_deg,
                        upper_angle_deg=upper.angle_deg,
                        phi_avg_deg=math.degrees(phi),
                        fiber_direction_unit_vector=fiber,
                    )
                    records.append(record)
                    event = UndulationEvent(
                        record=record,
                        ply_index=lower.ply_index,
                        tow_index=lower_idx,
                        tow_offset_m=lower_offset,
                        center_s_m=along_lower,
                        ramp_up_name=up_name,
                        ramp_down_name=down_name,
                    )
                    events_by_tow.setdefault(tow_key, []).append(event)
                    global_crossing_index += 1
        return events_by_tow, tuple(records)

    def _build_ply_tows(
        self,
        occ: OCCBackend,
        ply: Ply,
        clip: DimTag,
        events_by_tow: dict[tuple[int, int], list[UndulationEvent]],
        groups: list[PhysicalVolumeGroup],
    ) -> None:
        u = ply.direction
        n = ply.transverse
        z_center = ply.z_bottom_m + 0.5 * ply.cured_ply_thickness_m
        profile = UndulationProfile(
            cured_ply_thickness_m=ply.cured_ply_thickness_m,
            undulation_ratio=ply.undulation_ratio,
        )
        lu = profile.undulation_length_m
        s_min = -0.5 * ply.tow_length_for_clipping_m
        s_max = 0.5 * ply.tow_length_for_clipping_m

        for tow_idx, offset in enumerate(ply.center_offsets_m()):
            tow_name = f"TOW_PLY_{ply.ply_index}_TAG_{tow_idx}"
            flat_tags: list[DimTag] = []
            cursor = s_min
            events = sorted(events_by_tow.get((ply.ply_index, tow_idx), ()), key=lambda event: event.center_s_m)
            for event in events:
                ramp_start = max(s_min, event.center_s_m - lu)
                ramp_mid = min(s_max, event.center_s_m)
                ramp_end = min(s_max, event.center_s_m + lu)
                if ramp_start > cursor:
                    flat_tags.extend(_clip_tag(self._add_flat_segment(occ, ply, offset, cursor, ramp_start, z_center), clip))
                if ramp_mid > ramp_start:
                    up_tag = self._add_ramp_segment(
                        occ,
                        ply,
                        event,
                        offset,
                        ramp_start,
                        ramp_mid,
                        z_center,
                        rising=True,
                    )
                    clipped = _clip_tag(up_tag, clip)
                    if clipped:
                        groups.append(_undulation_group(ply, event, clipped, sign=1.0, name=event.ramp_up_name))
                if ramp_end > ramp_mid:
                    down_tag = self._add_ramp_segment(
                        occ,
                        ply,
                        event,
                        offset,
                        ramp_mid,
                        ramp_end,
                        z_center,
                        rising=False,
                    )
                    clipped = _clip_tag(down_tag, clip)
                    if clipped:
                        groups.append(_undulation_group(ply, event, clipped, sign=-1.0, name=event.ramp_down_name))
                cursor = max(cursor, ramp_end)
            if s_max > cursor:
                flat_tags.extend(_clip_tag(self._add_flat_segment(occ, ply, offset, cursor, s_max, z_center), clip))

            if flat_tags:
                ux, uy = ply.direction
                vx, vy = ply.transverse
                groups.append(
                    PhysicalVolumeGroup(
                        name=tow_name,
                        kind="straight_tow",
                        dim_tags=tuple(flat_tags),
                        ply_index=ply.ply_index,
                        nominal_angle_deg=ply.angle_deg,
                        fiber_direction_unit_vector=(ux, uy, 0.0),
                        transverse_in_plane_unit_vector=(vx, vy, 0.0),
                        through_thickness_unit_vector=(0.0, 0.0, 1.0),
                    )
                )

    def _add_flat_segment(
        self,
        occ: OCCBackend,
        ply: Ply,
        offset_m: float,
        start_s_m: float,
        end_s_m: float,
        z_center_m: float,
    ) -> DimTag:
        length = end_s_m - start_s_m
        if length <= 1.0e-12:
            raise ValueError("flat segment length must be positive")
        u = ply.direction
        n = ply.transverse
        center_s = 0.5 * (start_s_m + end_s_m)
        center = (
            u[0] * center_s + n[0] * offset_m,
            u[1] * center_s + n[1] * offset_m,
            z_center_m,
        )
        return occ.add_oriented_box(
            length_m=length,
            width_m=ply.tape_width_m,
            thickness_m=ply.cured_ply_thickness_m,
            angle_rad=ply.angle_rad,
            center_m=center,
        )

    def _add_ramp_segment(
        self,
        occ: OCCBackend,
        ply: Ply,
        event: UndulationEvent,
        offset_m: float,
        start_s_m: float,
        end_s_m: float,
        z_center_m: float,
        *,
        rising: bool,
    ) -> DimTag:
        u = ply.direction
        n = ply.transverse
        length = end_s_m - start_s_m
        if length <= 1.0e-12:
            raise ValueError("ramp segment length must be positive")
        center_s = 0.5 * (start_s_m + end_s_m)
        rise = ply.cured_ply_thickness_m if rising else -ply.cured_ply_thickness_m
        center = (
            u[0] * center_s + n[0] * offset_m,
            u[1] * center_s + n[1] * offset_m,
            z_center_m + 0.5 * ply.cured_ply_thickness_m,
        )
        return occ.add_ramped_box(
            projected_length_m=length,
            rise_m=rise,
            width_m=ply.tape_width_m,
            thickness_m=ply.cured_ply_thickness_m,
            angle_rad=ply.angle_rad,
            projected_center_m=center,
        )


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


def _det2(a: tuple[float, float], b: tuple[float, float]) -> float:
    return a[0] * b[1] - a[1] * b[0]


def _sub2(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
    return (a[0] - b[0], a[1] - b[1])


def _normalize3(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    norm = math.sqrt(sum(value * value for value in vector))
    if norm <= 0.0:
        raise ValueError("cannot normalize zero vector")
    return tuple(value / norm for value in vector)


def _cross3(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _clip_tag(body: DimTag, clip: DimTag) -> tuple[DimTag, ...]:
    out, _ = gmsh.model.occ.intersect([body.tuple], [clip.tuple], removeObject=True, removeTool=False)
    return tuple(DimTag(dim, tag) for dim, tag in out if dim == 3)


def _cut_tags(tags: tuple[DimTag, ...], tool_tuples: list[tuple[int, int]]) -> tuple[DimTag, ...]:
    new_tags: list[DimTag] = []
    for tag in tags:
        out, _ = gmsh.model.occ.cut([tag.tuple], tool_tuples, removeObject=True, removeTool=False)
        new_tags.extend(DimTag(dim, out_tag) for dim, out_tag in out if dim == 3)
    return tuple(new_tags)


def _replace_group_tags(group: PhysicalVolumeGroup, tags: tuple[DimTag, ...]) -> PhysicalVolumeGroup:
    return PhysicalVolumeGroup(
        name=group.name,
        kind=group.kind,
        dim_tags=tags,
        ply_index=group.ply_index,
        ply_indices=group.ply_indices,
        nominal_angle_deg=group.nominal_angle_deg,
        nominal_angles_deg=group.nominal_angles_deg,
        fiber_direction_unit_vector=group.fiber_direction_unit_vector,
        transverse_in_plane_unit_vector=group.transverse_in_plane_unit_vector,
        through_thickness_unit_vector=group.through_thickness_unit_vector,
        phi_avg_deg=group.phi_avg_deg,
        point_m=group.point_m,
        is_isotropic=group.is_isotropic,
    )


def _fragment_and_remap(groups: list[PhysicalVolumeGroup]) -> list[PhysicalVolumeGroup]:
    inputs: list[tuple[int, DimTag]] = []
    for group_idx, group in enumerate(groups):
        for tag in group.dim_tags:
            inputs.append((group_idx, tag))
    if not inputs:
        return groups
    out, maps = gmsh.model.occ.fragment([tag.tuple for _, tag in inputs], [])
    _ = out
    remapped: list[list[DimTag]] = [[] for _ in groups]
    for (group_idx, _tag), mapped in zip(inputs, maps):
        remapped[group_idx].extend(DimTag(dim, tag) for dim, tag in mapped if dim == 3)
    return [_replace_group_tags(group, tuple(tags)) for group, tags in zip(groups, remapped) if tags]


def _undulation_group(
    ply: Ply,
    event: UndulationEvent,
    tags: tuple[DimTag, ...],
    *,
    sign: float,
    name: str,
) -> PhysicalVolumeGroup:
    profile = UndulationProfile(
        cured_ply_thickness_m=ply.cured_ply_thickness_m,
        undulation_ratio=ply.undulation_ratio,
    )
    phi = profile.phi_avg_rad
    ux, uy = ply.direction
    vx, vy = ply.transverse
    fiber = _normalize3((ux * math.cos(phi), uy * math.cos(phi), sign * math.sin(phi)))
    transverse = (vx, vy, 0.0)
    through = _normalize3(_cross3(fiber, transverse))
    return PhysicalVolumeGroup(
        name=name,
        kind="undulation",
        dim_tags=tags,
        ply_index=ply.ply_index,
        ply_indices=(event.record.lower_ply_index, event.record.upper_ply_index),
        nominal_angle_deg=ply.angle_deg,
        nominal_angles_deg=(event.record.lower_angle_deg, event.record.upper_angle_deg),
        fiber_direction_unit_vector=fiber,
        transverse_in_plane_unit_vector=transverse,
        through_thickness_unit_vector=through,
        phi_avg_deg=math.degrees(phi),
        point_m=event.record.point_m,
    )


def _group_to_orientation_dict(group: PhysicalVolumeGroup) -> dict[str, object]:
    if group.is_isotropic:
        data: dict[str, object] = {
            "name": group.name,
            "kind": group.kind,
            "isotropic": True,
        }
        if group.ply_index is not None:
            data["ply_index"] = group.ply_index
        return data

    if group.fiber_direction_unit_vector is None or group.transverse_in_plane_unit_vector is None:
        raise ValueError(f"group {group.name} is missing orientation vectors")
    data = {
        "name": group.name,
        "kind": group.kind,
        "fiber_direction_unit_vector": list(group.fiber_direction_unit_vector),
        "transverse_in_plane_unit_vector": list(group.transverse_in_plane_unit_vector),
        "through_thickness_unit_vector": list(
            group.through_thickness_unit_vector
            if group.through_thickness_unit_vector is not None
            else _normalize3(_cross3(group.fiber_direction_unit_vector, group.transverse_in_plane_unit_vector))
        ),
    }
    if group.ply_index is not None:
        data["ply_index"] = group.ply_index
    if group.ply_indices is not None:
        data["ply_indices"] = list(group.ply_indices)
    if group.nominal_angle_deg is not None:
        data["nominal_angle_deg"] = group.nominal_angle_deg
    if group.nominal_angles_deg is not None:
        data["nominal_angles_deg"] = list(group.nominal_angles_deg)
    if group.phi_avg_deg is not None:
        data["phi_avg_deg"] = group.phi_avg_deg
    if group.point_m is not None:
        data["point_m"] = list(group.point_m)
    return data
