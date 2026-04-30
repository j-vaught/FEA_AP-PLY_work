"""Single-ply tow array and resin-pocket construction."""

from __future__ import annotations

import math
from dataclasses import dataclass

import gmsh

from kok_geom.config import MM
from kok_geom.geometry.tow import Tow
from kok_geom.occ_backend import DimTag, OCCBackend


@dataclass(frozen=True)
class PlyTow:
    name: str
    dim_tag: DimTag
    orientation: tuple[float, float, float]
    physical_group: int
    center_offset_m: float


@dataclass(frozen=True)
class ResinPocket:
    name: str
    dim_tags: tuple[DimTag, ...]
    physical_group: int


@dataclass(frozen=True)
class PlySolid:
    ply_index: int
    tows: tuple[PlyTow, ...]
    resin: ResinPocket

    @property
    def volume_tags(self) -> tuple[DimTag, ...]:
        tags: list[DimTag] = [tow.dim_tag for tow in self.tows]
        tags.extend(self.resin.dim_tags)
        return tuple(tags)


@dataclass(frozen=True)
class Ply:
    """A finite single ply made from parallel slit-tape solids plus resin."""

    size_x_m: float
    size_y_m: float
    angle_deg: float = 0.0
    tape_width_m: float = 6.35 * MM
    cured_ply_thickness_m: float = 0.18 * MM
    undulation_ratio: float = 0.09
    tape_spacing: int = 1
    ply_index: int = 1
    z_bottom_m: float = 0.0
    origin_offset_m: float = 0.0

    @classmethod
    def from_mm(
        cls,
        *,
        size_x_mm: float,
        size_y_mm: float,
        angle_deg: float = 0.0,
        tape_width_mm: float = 6.35,
        cured_ply_thickness_mm: float = 0.18,
        undulation_ratio: float = 0.09,
        tape_spacing: int = 1,
        ply_index: int = 1,
        z_bottom_mm: float = 0.0,
        origin_offset_mm: float = 0.0,
    ) -> "Ply":
        return cls(
            size_x_m=size_x_mm * MM,
            size_y_m=size_y_mm * MM,
            angle_deg=angle_deg,
            tape_width_m=tape_width_mm * MM,
            cured_ply_thickness_m=cured_ply_thickness_mm * MM,
            undulation_ratio=undulation_ratio,
            tape_spacing=tape_spacing,
            ply_index=ply_index,
            z_bottom_m=z_bottom_mm * MM,
            origin_offset_m=origin_offset_mm * MM,
        )

    def __post_init__(self) -> None:
        if self.size_x_m <= 0.0 or self.size_y_m <= 0.0:
            raise ValueError("ply in-plane dimensions must be positive")
        if self.tape_width_m <= 0.0 or self.cured_ply_thickness_m <= 0.0:
            raise ValueError("tape width and ply thickness must be positive")
        if self.undulation_ratio <= 0.0:
            raise ValueError("undulation_ratio must be positive")
        if self.tape_spacing not in {1, 2, 3}:
            raise ValueError("tape_spacing must be one of {1, 2, 3}")

    @property
    def angle_rad(self) -> float:
        return math.radians(self.angle_deg)

    @property
    def direction(self) -> tuple[float, float]:
        return (math.cos(self.angle_rad), math.sin(self.angle_rad))

    @property
    def transverse(self) -> tuple[float, float]:
        ux, uy = self.direction
        return (-uy, ux)

    @property
    def active_pitch_m(self) -> float:
        return (self.tape_spacing + 1) * self.tape_width_m

    @property
    def panel_area_m2(self) -> float:
        return self.size_x_m * self.size_y_m

    @property
    def nominal_volume_m3(self) -> float:
        return self.panel_area_m2 * self.cured_ply_thickness_m

    @property
    def tow_length_for_clipping_m(self) -> float:
        return math.hypot(self.size_x_m, self.size_y_m) + 2.0 * self.tape_width_m

    def transverse_span_m(self) -> float:
        nx, ny = self.transverse
        corners = (
            (-self.size_x_m / 2.0, -self.size_y_m / 2.0),
            (self.size_x_m / 2.0, -self.size_y_m / 2.0),
            (self.size_x_m / 2.0, self.size_y_m / 2.0),
            (-self.size_x_m / 2.0, self.size_y_m / 2.0),
        )
        projections = [x * nx + y * ny for x, y in corners]
        return max(projections) - min(projections)

    def expected_tow_count(self) -> int:
        return max(1, int(math.floor((self.transverse_span_m() + 1.0e-12) / self.active_pitch_m)))

    def center_offsets_m(self) -> tuple[float, ...]:
        count = self.expected_tow_count()
        first = -0.5 * (count - 1) * self.active_pitch_m + self.origin_offset_m
        return tuple(first + i * self.active_pitch_m for i in range(count))

    def build_occ(self, backend: OCCBackend | None = None) -> PlySolid:
        occ = backend or OCCBackend()
        panel = occ.add_box(
            x0=-self.size_x_m / 2.0,
            y0=-self.size_y_m / 2.0,
            z0=self.z_bottom_m,
            dx=self.size_x_m,
            dy=self.size_y_m,
            dz=self.cured_ply_thickness_m,
        )

        clipped_tows: list[tuple[str, DimTag, tuple[float, float, float], float]] = []
        nx, ny = self.transverse
        z_center = self.z_bottom_m + 0.5 * self.cured_ply_thickness_m

        for idx, offset in enumerate(self.center_offsets_m()):
            name = f"TOW_PLY_{self.ply_index}_TAG_{idx}"
            tow = Tow(
                length_m=self.tow_length_for_clipping_m,
                angle_deg=self.angle_deg,
                tape_width_m=self.tape_width_m,
                cured_ply_thickness_m=self.cured_ply_thickness_m,
                undulation_ratio=self.undulation_ratio,
                center_m=(nx * offset, ny * offset, z_center),
                name=name,
            )
            raw = tow.build_occ(occ)
            out, _ = gmsh.model.occ.intersect(
                [raw.dim_tag.tuple],
                [panel.tuple],
                removeObject=True,
                removeTool=False,
            )
            volume_tags = [DimTag(dim, tag) for dim, tag in out if dim == 3]
            if not volume_tags:
                continue
            clipped_tows.append((name, volume_tags[0], tow.orientation, offset))

        tow_tools = [body.tuple for _, body, _, _ in clipped_tows]
        resin_out, _ = gmsh.model.occ.cut([panel.tuple], tow_tools, removeObject=True, removeTool=False)
        resin_tags = tuple(DimTag(dim, tag) for dim, tag in resin_out if dim == 3)
        if not resin_tags:
            raise RuntimeError("OCC cut did not produce a resin pocket volume")

        gmsh.model.occ.synchronize()

        ply_tows: list[PlyTow] = []
        for name, body, orientation, offset in clipped_tows:
            group = gmsh.model.addPhysicalGroup(3, [body.tag])
            gmsh.model.setPhysicalName(3, group, name)
            ply_tows.append(
                PlyTow(
                    name=name,
                    dim_tag=body,
                    orientation=orientation,
                    physical_group=group,
                    center_offset_m=offset,
                )
            )

        resin_name = f"RESIN_PLY_{self.ply_index}"
        resin_group = gmsh.model.addPhysicalGroup(3, [tag.tag for tag in resin_tags])
        gmsh.model.setPhysicalName(3, resin_group, resin_name)
        return PlySolid(
            ply_index=self.ply_index,
            tows=tuple(ply_tows),
            resin=ResinPocket(name=resin_name, dim_tags=resin_tags, physical_group=resin_group),
        )
