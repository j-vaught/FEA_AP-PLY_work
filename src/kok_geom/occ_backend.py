"""Small GMSH-OCC adapter used by the geometry primitives."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

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
