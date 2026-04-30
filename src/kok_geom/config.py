"""Typed configuration schema for the AP-PLY geometry port.

The external configuration follows the project documents and accepts common
engineering dimensions in millimetres. Geometry builders convert those fields
to SI metres before touching GMSH.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

MM = 1.0e-3
INCH_TO_MM = 25.4


class PanelConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    size_x_mm: float = Field(default=25.0, gt=0.0)
    size_y_mm: float = Field(default=25.0, gt=0.0)
    n_plies: int = Field(default=1, ge=1)
    polygon_mm: list[tuple[float, float]] | None = None
    symmetry: Literal["none", "midplane"] = "none"

    @model_validator(mode="after")
    def default_polygon(self) -> "PanelConfig":
        if self.polygon_mm is None:
            hx = self.size_x_mm / 2.0
            hy = self.size_y_mm / 2.0
            self.polygon_mm = [(-hx, -hy), (hx, -hy), (hx, hy), (-hx, hy)]
        if len(self.polygon_mm) < 3:
            raise ValueError("panel.polygon_mm must contain at least three points")
        return self

    @property
    def size_x_m(self) -> float:
        return self.size_x_mm * MM

    @property
    def size_y_m(self) -> float:
        return self.size_y_mm * MM


class LaydownConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    fiber_angles_deg: list[float] = Field(default_factory=lambda: [0.0])
    placement_sequence: str = "10"
    angle_shift_deg: float = 0.0
    tape_width_mm: float | list[float] = 6.35
    cured_ply_thickness_mm: float = Field(default=0.18, gt=0.0)
    undulation_ratio: float = Field(default=0.09, gt=0.0)
    tape_spacing: int = Field(default=1, ge=1)
    shorthand: str | None = None

    @field_validator("fiber_angles_deg")
    @classmethod
    def validate_angles(cls, value: list[float]) -> list[float]:
        if not value:
            raise ValueError("laydown.fiber_angles_deg cannot be empty")
        return [float(v) for v in value]

    @field_validator("placement_sequence")
    @classmethod
    def validate_placement_sequence(cls, value: str) -> str:
        if not value or any(ch not in "01" for ch in value):
            raise ValueError("laydown.placement_sequence must be a non-empty 0/1 string")
        if "1" not in value:
            raise ValueError("laydown.placement_sequence must contain at least one active channel")
        return value

    @field_validator("tape_spacing")
    @classmethod
    def validate_tape_spacing(cls, value: int) -> int:
        if value not in {1, 2, 3}:
            raise ValueError("laydown.tape_spacing must be one of {1, 2, 3}")
        return value

    @field_validator("tape_width_mm")
    @classmethod
    def validate_tape_width(cls, value: float | list[float]) -> float | list[float]:
        if isinstance(value, list):
            if not value:
                raise ValueError("laydown.tape_width_mm list cannot be empty")
            if any(v <= 0.0 for v in value):
                raise ValueError("laydown.tape_width_mm values must be positive")
            return [float(v) for v in value]
        if value <= 0.0:
            raise ValueError("laydown.tape_width_mm must be positive")
        return float(value)

    def tape_width_for_ply_mm(self, ply_index: int) -> float:
        if isinstance(self.tape_width_mm, list):
            if ply_index < 0 or ply_index >= len(self.tape_width_mm):
                raise IndexError("ply index outside tape_width_mm list")
            return self.tape_width_mm[ply_index]
        return self.tape_width_mm

    @property
    def cured_ply_thickness_m(self) -> float:
        return self.cured_ply_thickness_mm * MM


class MeshConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    in_plane_target_mm_impact_zone: float = Field(default=0.5, gt=0.0)
    in_plane_target_mm_far_field: float = Field(default=0.5, gt=0.0)
    through_thickness_target_mm: float = Field(default=0.06, gt=0.0)
    graded_zone_radius_mm: float = Field(default=0.0, ge=0.0)
    element_order: Literal[1, 2] = 2

    @property
    def target_size_m(self) -> float:
        return min(
            self.in_plane_target_mm_impact_zone,
            self.in_plane_target_mm_far_field,
            self.through_thickness_target_mm,
        ) * MM


class OutputConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    msh_path: str = "build/panel.msh"
    inp_path: str | None = None
    orientations_json_path: str = "build/orientations.json"
    msh_format: Literal["msh4_ascii", "msh4_binary", "msh2_ascii"] = "msh4_ascii"


class KokConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    panel: PanelConfig = Field(default_factory=PanelConfig)
    laydown: LaydownConfig = Field(default_factory=LaydownConfig)
    mesh: MeshConfig = Field(default_factory=MeshConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)

    @model_validator(mode="after")
    def validate_width_list(self) -> "KokConfig":
        width = self.laydown.tape_width_mm
        if isinstance(width, list) and len(width) != self.panel.n_plies:
            raise ValueError("laydown.tape_width_mm list length must equal panel.n_plies")
        return self

    @classmethod
    def from_file(cls, path: str | Path) -> "KokConfig":
        data = json.loads(Path(path).read_text())
        return cls.model_validate(data)

    @classmethod
    def from_shorthand(cls, shorthand: str, **overrides: Any) -> "KokConfig":
        data = parse_shorthand(shorthand)
        data.update(overrides)
        return cls.model_validate(data)

    def to_file(self, path: str | Path) -> None:
        Path(path).write_text(self.canonical_json() + "\n")

    def canonical_json(self) -> str:
        return json.dumps(self.model_dump(mode="json"), indent=2, sort_keys=True)

    def config_hash(self) -> str:
        return hashlib.sha256(self.canonical_json().encode("utf-8")).hexdigest()


def infer_tape_spacing(placement_sequence: str) -> int:
    """Infer a skip factor from a regular active/inactive placement string."""

    if placement_sequence == "1":
        return 1
    positions = [idx for idx, ch in enumerate(placement_sequence) if ch == "1"]
    if len(positions) >= 2:
        period = positions[1] - positions[0]
        return max(1, min(3, period - 1))
    zeros_after_first = len(placement_sequence) - placement_sequence.index("1") - 1
    return max(1, min(3, zeros_after_first))


def parse_shorthand(shorthand: str) -> dict[str, Any]:
    """Parse the UofSC and Nagelsmit shorthand forms described in the plan."""

    text = shorthand.strip()
    if text.startswith("[("):
        return _parse_nagelsmit_shorthand(text)
    return _parse_uofsc_shorthand(text)


def _parse_uofsc_shorthand(text: str) -> dict[str, Any]:
    match = re.fullmatch(
        r"\[([^\]]+)\]\[([01]+)\]\[([+-]?(?:\d+(?:\.\d*)?|\.\d+))\]\[([+-]?(?:\d+(?:\.\d*)?|\.\d+))\]",
        text,
    )
    if not match:
        raise ValueError(f"Unsupported shorthand: {text}")
    angles = [float(part.strip()) for part in match.group(1).split(",")]
    placement = match.group(2)
    angle_shift = float(match.group(3))
    tape_width_mm = float(match.group(4))
    return {
        "panel": {"size_x_mm": 25.0, "size_y_mm": 25.0, "n_plies": len(angles), "symmetry": "none"},
        "laydown": {
            "fiber_angles_deg": angles,
            "placement_sequence": placement,
            "angle_shift_deg": angle_shift,
            "tape_width_mm": tape_width_mm,
            "cured_ply_thickness_mm": 0.18,
            "undulation_ratio": 0.09,
            "tape_spacing": infer_tape_spacing(placement),
            "shorthand": text,
        },
    }


def _parse_nagelsmit_shorthand(text: str) -> dict[str, Any]:
    pattern = re.compile(
        r"\[\(([^)]+)\)(\d+)[x×]([0-9.]+)/\(([^)]+)\)(\d+)[x×]([0-9.]+)\]_(\d+)S"
    )
    match = pattern.fullmatch(text)
    if not match:
        raise ValueError(f"Unsupported Nagelsmit shorthand: {text}")
    first_angles = [float(part) for part in match.group(1).split("/")]
    second_angles = [float(part) for part in match.group(4).split("/")]
    spacing_a = int(match.group(2))
    spacing_b = int(match.group(5))
    if spacing_a != spacing_b:
        raise ValueError("Mixed tape spacing in one shorthand is not supported")
    width_a_in = float(match.group(3))
    width_b_in = float(match.group(6))
    if not math.isclose(width_a_in, width_b_in, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("Mixed tape widths in one shorthand are not supported")
    repeat_count = int(match.group(7))
    angles = first_angles + second_angles
    return {
        "panel": {
            "size_x_mm": 40.0,
            "size_y_mm": 40.0,
            "n_plies": len(angles) * repeat_count * 2,
            "symmetry": "midplane",
        },
        "laydown": {
            "fiber_angles_deg": angles,
            "placement_sequence": "1" + "0" * spacing_a,
            "angle_shift_deg": 0.0,
            "tape_width_mm": width_a_in * INCH_TO_MM,
            "cured_ply_thickness_mm": 0.18,
            "undulation_ratio": 0.09,
            "tape_spacing": spacing_a,
            "shorthand": text,
        },
    }
