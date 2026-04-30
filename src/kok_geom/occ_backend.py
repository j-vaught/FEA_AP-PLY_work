"""Small GMSH-OCC adapter used by the geometry primitives."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from collections.abc import Callable

import gmsh


@dataclass(frozen=True)
class DimTag:
    dim: int
    tag: int

    @property
    def tuple(self) -> tuple[int, int]:
        return (self.dim, self.tag)


class GmshSession:
    """Context manager for isolated GMSH model builds."""

    def __init__(self, model_name: str = "kok_geom", finalize: bool = False) -> None:
        self.model_name = model_name
        self.finalize = finalize
        self._started_here = False

    def __enter__(self) -> "GmshSession":
        if not gmsh.isInitialized():
            gmsh.initialize()
            self._started_here = True
        gmsh.clear()
        gmsh.model.add(self.model_name)
        gmsh.option.setNumber("General.Terminal", 0)
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # type: ignore[no-untyped-def]
        if self.finalize and self._started_here and gmsh.isInitialized():
            gmsh.finalize()


class OCCBackend:
    """Thin wrapper over the subset of gmsh.model.occ used in M1."""

    def add_box(
        self,
        *,
        x0: float,
        y0: float,
        z0: float,
        dx: float,
        dy: float,
        dz: float,
    ) -> DimTag:
        if dx <= 0.0 or dy <= 0.0 or dz <= 0.0:
            raise ValueError("box dimensions must be positive")
        return DimTag(3, gmsh.model.occ.addBox(x0, y0, z0, dx, dy, dz))

    def add_oriented_box(
        self,
        *,
        length_m: float,
        width_m: float,
        thickness_m: float,
        angle_rad: float,
        center_m: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> DimTag:
        if length_m <= 0.0 or width_m <= 0.0 or thickness_m <= 0.0:
            raise ValueError("box dimensions must be positive")

        tag = gmsh.model.occ.addBox(
            -length_m / 2.0,
            -width_m / 2.0,
            -thickness_m / 2.0,
            length_m,
            width_m,
            thickness_m,
        )
        if not math.isclose(angle_rad, 0.0, abs_tol=1.0e-15):
            gmsh.model.occ.rotate([(3, tag)], 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, angle_rad)
        if any(not math.isclose(v, 0.0, abs_tol=1.0e-15) for v in center_m):
            gmsh.model.occ.translate([(3, tag)], *center_m)
        return DimTag(3, tag)

    def add_ramped_box(
        self,
        *,
        projected_length_m: float,
        rise_m: float,
        width_m: float,
        thickness_m: float,
        angle_rad: float,
        projected_center_m: tuple[float, float, float],
    ) -> DimTag:
        """Add a rectangular tow ramp as a tilted OCC box.

        The local x-axis projects to ``projected_length_m`` in the laminate
        plane and rises by ``rise_m`` through thickness.
        """

        if projected_length_m <= 0.0 or width_m <= 0.0 or thickness_m <= 0.0:
            raise ValueError("ramped box dimensions must be positive")
        slope_rad = math.atan2(rise_m, projected_length_m)
        true_length_m = math.hypot(projected_length_m, rise_m)
        tag = gmsh.model.occ.addBox(
            -true_length_m / 2.0,
            -width_m / 2.0,
            -thickness_m / 2.0,
            true_length_m,
            width_m,
            thickness_m,
        )
        if not math.isclose(slope_rad, 0.0, abs_tol=1.0e-15):
            gmsh.model.occ.rotate([(3, tag)], 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -slope_rad)
        if not math.isclose(angle_rad, 0.0, abs_tol=1.0e-15):
            gmsh.model.occ.rotate([(3, tag)], 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, angle_rad)
        gmsh.model.occ.translate([(3, tag)], *projected_center_m)
        return DimTag(3, tag)

    def add_swept_box(
        self,
        *,
        direction_xy: tuple[float, float],
        transverse_xy: tuple[float, float],
        width_m: float,
        thickness_m: float,
        center_at: Callable[[float], tuple[float, float, float]],
        tangent_at: Callable[[float], tuple[float, float, float]],
        s_values_m: list[float],
    ) -> DimTag:
        """Loft a rectangular tow section through sampled centerline stations."""

        if width_m <= 0.0 or thickness_m <= 0.0:
            raise ValueError("swept section dimensions must be positive")
        if len(s_values_m) < 2:
            raise ValueError("at least two sweep stations are required")

        vx, vy = transverse_xy
        width_axis = (vx, vy, 0.0)
        width_norm = math.sqrt(sum(value * value for value in width_axis))
        if width_norm <= 0.0:
            raise ValueError("transverse axis must be nonzero")
        width_axis = tuple(value / width_norm for value in width_axis)

        wires: list[int] = []
        for s in s_values_m:
            cx, cy, cz = center_at(s)
            tx, ty, tz = tangent_at(s)
            tangent_norm = math.sqrt(tx * tx + ty * ty + tz * tz)
            if tangent_norm <= 0.0:
                raise ValueError("sweep tangent must be nonzero")
            tangent = (tx / tangent_norm, ty / tangent_norm, tz / tangent_norm)

            # e_thickness = tangent x width_axis. This is mostly global +Z,
            # with a small in-plane component during the Kok Fig. 4 ramp.
            ex, ey, ez = width_axis
            thickness_axis = (
                tangent[1] * ez - tangent[2] * ey,
                tangent[2] * ex - tangent[0] * ez,
                tangent[0] * ey - tangent[1] * ex,
            )
            thickness_norm = math.sqrt(sum(value * value for value in thickness_axis))
            if thickness_norm <= 0.0:
                raise ValueError("degenerate swept section axes")
            thickness_axis = tuple(value / thickness_norm for value in thickness_axis)

            points: list[int] = []
            for wsign, tsign in ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)):
                x = cx + wsign * 0.5 * width_m * width_axis[0] + tsign * 0.5 * thickness_m * thickness_axis[0]
                y = cy + wsign * 0.5 * width_m * width_axis[1] + tsign * 0.5 * thickness_m * thickness_axis[1]
                z = cz + wsign * 0.5 * width_m * width_axis[2] + tsign * 0.5 * thickness_m * thickness_axis[2]
                points.append(gmsh.model.occ.addPoint(x, y, z))
            lines = [gmsh.model.occ.addLine(points[idx], points[(idx + 1) % 4]) for idx in range(4)]
            wires.append(gmsh.model.occ.addWire(lines, checkClosed=True))

        out = gmsh.model.occ.addThruSections(wires, makeSolid=True, makeRuled=True)
        volumes = [DimTag(dim, tag) for dim, tag in out if dim == 3]
        if not volumes:
            raise RuntimeError("OCC addThruSections did not return a volume")
        return volumes[0]

    def synchronize(self) -> None:
        gmsh.model.occ.synchronize()

    def volume(self, body: DimTag | tuple[int, int]) -> float:
        dim, tag = body.tuple if isinstance(body, DimTag) else body
        if dim != 3:
            raise ValueError("volume is only defined for dim=3 entities")
        return float(gmsh.model.occ.getMass(dim, tag))

    def add_physical_volume(self, body: DimTag, name: str) -> int:
        self.synchronize()
        group = gmsh.model.addPhysicalGroup(3, [body.tag])
        gmsh.model.setPhysicalName(3, group, name)
        return group

    def write(self, path: str | Path) -> None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(path))
