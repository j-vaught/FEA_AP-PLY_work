"""Single-tow path, cross-section, and undulation primitives."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from kok_geom.config import MM
from kok_geom.occ_backend import DimTag, OCCBackend


@dataclass(frozen=True)
class UndulationProfile:
    """Half-period sinusoidal tow climb used for Kok-style undulations.

    The profile is the clean-room commitment documented in
    ``plan/kok_port_plan.md`` OQ-1:

    ``z(x) = t/2 * (1 - cos(pi x / L_u))`` for ``x`` in ``[0, L_u]``.
    """

    cured_ply_thickness_m: float = 0.18 * MM
    undulation_ratio: float = 0.09
    quadrature_points: int = 96

    def __post_init__(self) -> None:
        if self.cured_ply_thickness_m <= 0.0:
            raise ValueError("cured_ply_thickness_m must be positive")
        if self.undulation_ratio <= 0.0:
            raise ValueError("undulation_ratio must be positive")
        if self.quadrature_points < 16:
            raise ValueError("quadrature_points must be at least 16")

    @property
    def undulation_length_m(self) -> float:
        return self.cured_ply_thickness_m / self.undulation_ratio

    @property
    def half_length_m(self) -> float:
        return 0.5 * self.undulation_length_m

    def z_m(self, x_m: float | np.ndarray) -> float | np.ndarray:
        x = np.asarray(x_m)
        z = 0.5 * self.cured_ply_thickness_m * (
            1.0 - np.cos(np.pi * x / self.undulation_length_m)
        )
        if np.isscalar(x_m):
            return float(z)
        return z

    def dzdx(self, x_m: float | np.ndarray) -> float | np.ndarray:
        x = np.asarray(x_m)
        slope = (
            0.5
            * self.cured_ply_thickness_m
            * np.pi
            / self.undulation_length_m
            * np.sin(np.pi * x / self.undulation_length_m)
        )
        if np.isscalar(x_m):
            return float(slope)
        return slope

    def phi_rad(self, x_m: float | np.ndarray) -> float | np.ndarray:
        value = np.arctan(self.dzdx(x_m))
        if np.isscalar(x_m):
            return float(value)
        return value

    @property
    def phi_avg_rad(self) -> float:
        """Average out-of-plane angle over the rising half-undulation."""

        nodes, weights = np.polynomial.legendre.leggauss(self.quadrature_points)
        a = 0.0
        b = self.half_length_m
        x = 0.5 * (b - a) * nodes + 0.5 * (a + b)
        integral = 0.5 * (b - a) * np.sum(weights * self.phi_rad(x))
        return float(integral / self.half_length_m)

    @property
    def phi_avg_deg(self) -> float:
        return math.degrees(self.phi_avg_rad)

    def centerline_points(
        self,
        *,
        length_m: float,
        samples: int = 33,
        centered: bool = True,
    ) -> np.ndarray:
        """Return a straight tow centerline with the sinusoidal climb embedded."""

        if samples < 3:
            raise ValueError("samples must be at least 3")
        x = np.linspace(0.0, self.undulation_length_m, samples)
        z = self.z_m(x)
        if centered:
            x = x - self.undulation_length_m / 2.0
        if length_m < self.undulation_length_m:
            raise ValueError("length_m must be at least the undulation length")
        return np.column_stack([x, np.zeros_like(x), z])


@dataclass(frozen=True)
class TowSolid:
    dim_tag: DimTag
    orientation: tuple[float, float, float]
    physical_group: int | None = None
    physical_name: str | None = None


@dataclass(frozen=True)
class Tow:
    """One slit-tape represented as an oriented rectangular OCC solid."""

    length_m: float
    angle_deg: float = 0.0
    tape_width_m: float = 6.35 * MM
    cured_ply_thickness_m: float = 0.18 * MM
    undulation_ratio: float = 0.09
    center_m: tuple[float, float, float] = (0.0, 0.0, 0.0)
    name: str = "TOW_PLY_1_TAG_0"

    @classmethod
    def from_mm(
        cls,
        *,
        length_mm: float,
        angle_deg: float = 0.0,
        tape_width_mm: float = 6.35,
        cured_ply_thickness_mm: float = 0.18,
        undulation_ratio: float = 0.09,
        center_mm: tuple[float, float, float] = (0.0, 0.0, 0.0),
        name: str = "TOW_PLY_1_TAG_0",
    ) -> "Tow":
        return cls(
            length_m=length_mm * MM,
            angle_deg=angle_deg,
            tape_width_m=tape_width_mm * MM,
            cured_ply_thickness_m=cured_ply_thickness_mm * MM,
            undulation_ratio=undulation_ratio,
            center_m=tuple(v * MM for v in center_mm),
            name=name,
        )

    def __post_init__(self) -> None:
        if self.length_m <= 0.0:
            raise ValueError("length_m must be positive")
        if self.tape_width_m <= 0.0:
            raise ValueError("tape_width_m must be positive")
        if self.cured_ply_thickness_m <= 0.0:
            raise ValueError("cured_ply_thickness_m must be positive")
        if self.undulation_ratio <= 0.0:
            raise ValueError("undulation_ratio must be positive")

    @property
    def angle_rad(self) -> float:
        return math.radians(self.angle_deg)

    @property
    def orientation(self) -> tuple[float, float, float]:
        vx = math.cos(self.angle_rad)
        vy = math.sin(self.angle_rad)
        norm = math.hypot(vx, vy)
        return (vx / norm, vy / norm, 0.0)

    @property
    def transverse_in_plane(self) -> tuple[float, float, float]:
        vx, vy, _ = self.orientation
        return (-vy, vx, 0.0)

    @property
    def undulation_profile(self) -> UndulationProfile:
        return UndulationProfile(
            cured_ply_thickness_m=self.cured_ply_thickness_m,
            undulation_ratio=self.undulation_ratio,
        )

    @property
    def expected_volume_m3(self) -> float:
        return self.length_m * self.tape_width_m * self.cured_ply_thickness_m

    def build_occ(
        self,
        backend: OCCBackend | None = None,
        *,
        add_physical_group: bool = False,
    ) -> TowSolid:
        occ = backend or OCCBackend()
        body = occ.add_oriented_box(
            length_m=self.length_m,
            width_m=self.tape_width_m,
            thickness_m=self.cured_ply_thickness_m,
            angle_rad=self.angle_rad,
            center_m=self.center_m,
        )
        group = None
        if add_physical_group:
            group = occ.add_physical_volume(body, self.name)
        return TowSolid(
            dim_tag=body,
            orientation=self.orientation,
            physical_group=group,
            physical_name=self.name if group is not None else None,
        )
