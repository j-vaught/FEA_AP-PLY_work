"""Publication-grade visualization pass for FEA_AP-PLY.

Outputs:
- stage 01-09 2400 x 1600 four-panel composites
- stage 01-06 1920 x 1080 H.264 split-screen motion videos
- stage 10-11 wireframes where VTK exists
- root Typst methodology figures

The script uses PyVista/OSMesa for 3D fields, Typst + CeTZ for 2D charts, PIL
for assembly, and ffmpeg for video muxing. It intentionally does not edit
verification results.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
import pathlib
import shutil
import subprocess
import tempfile
from typing import Iterable

import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pyvista as pv


ROOT = pathlib.Path(__file__).resolve().parents[1]

GARNET = "#73000A"
BLACK = "#000000"
WHITE = "#FFFFFF"
CHARCOAL = "#363636"
BLACK70 = "#5C5C5C"
BLACK50 = "#A2A2A2"
BLACK30 = "#C7C7C7"
BLACK10 = "#ECECEC"
WARM_GREY = "#676156"
SANDSTORM = "#FFF2E3"
ROSE = "#CC2E40"
ATLANTIC = "#466A9F"
CONGAREE = "#1F414D"
HORSESHOE = "#65780B"
GRASS = "#CED318"
HONEYCOMB = "#A49137"

BRAND_RAMP = [SANDSTORM, HONEYCOMB, GARNET, CONGAREE]
PASS_COLOR = HORSESHOE
INCONCLUSIVE_COLOR = HONEYCOMB
FAIL_COLOR = GARNET
NA_COLOR = BLACK50

COMPOSITE_SIZE = (2400, 1600)
PANEL_SIZE = (1110, 640)
VIDEO_SIZE = (1920, 1080)
VIDEO_FPS = 30
VIDEO_MOTION_SECONDS = 8
VIDEO_CARD_SECONDS = 1


@dataclasses.dataclass(frozen=True)
class FieldChoice:
    scalar: str
    location: str
    title: str
    units: str


@dataclasses.dataclass(frozen=True)
class StageConfig:
    stage: int
    slug: str
    title: str
    verdict: str
    key_result: str
    vtk_globs: tuple[str, ...]
    third_field: str
    chart_kind: str
    final_index: int = -1
    view_scale: float = 0.58

    @property
    def stage_dir(self) -> pathlib.Path:
        return ROOT / "tests" / self.slug

    @property
    def figures_dir(self) -> pathlib.Path:
        return self.stage_dir / "figures"

    @property
    def prefix(self) -> str:
        return f"stage_{self.stage:02d}"


@dataclasses.dataclass(frozen=True)
class VideoFrame:
    path: pathlib.Path
    scale: float = 1.0


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/dejavu/DejaVuSans.ttf",
    ]
    for candidate in candidates:
        path = pathlib.Path(candidate)
        if path.exists():
            return ImageFont.truetype(str(path), size)
    return ImageFont.load_default()


FONT_TITLE = font(42, True)
FONT_SUBTITLE = font(24, False)
FONT_PANEL = font(25, True)
FONT_BODY = font(22, False)
FONT_SMALL = font(18, False)
FONT_MONO = font(18, False)


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def stage_configs() -> dict[int, StageConfig]:
    return {
        1: StageConfig(
            1,
            "stage_01_beam_bending",
            "Stage 01 Beam Bending",
            "PASS",
            "3pt and 4pt M1 Euler-Bernoulli errors: 0.968%",
            ("runs/stage01_4pt_M1_richviz/stage01_4pt_M1_richvizA*.vtk",),
            "eq_strain_percent",
            "stage01",
        ),
        2: StageConfig(
            2,
            "stage_02_cantilever_large",
            "Stage 02 Large-Deflection Cantilever",
            "INCONCLUSIVE",
            "alpha=1 passes; alpha=3/5 blocked by nonlinear convergence",
            ("runs/stage02_implicit_baseline_a1p00_richviz/stage02_implicit_baseline_a1p00_richvizA*.vtk",),
            "eq_strain_percent",
            "stage02",
        ),
        3: StageConfig(
            3,
            "stage_03_iso_dogbone_E8",
            "Stage 03 Isotropic ASTM E8 Dogbone",
            "PASS",
            "stress/load error 0.021%; apparent E error 0.044%",
            (
                "runs/stage03_medium_lf1p20_pubviz/stage03_medium_lf1p20_pubvizA*.vtk",
                "runs/stage03_medium_lf1p20/stage03_medium_lf1p20A*.vtk",
            ),
            "plastic_strain_percent",
            "stage03",
        ),
        4: StageConfig(
            4,
            "stage_04_open_hole_kirsch",
            "Stage 04 Open-Hole Kirsch/Howland",
            "PASS",
            "Kt=3.0956 vs 3.035; Howland error 1.996%",
            (
                "runs/stage04_ntheta064_pubviz/stage04_ntheta064_pubvizA*.vtk",
                "runs/stage04_ntheta064/stage04_ntheta064A*.vtk",
            ),
            "eq_strain_percent",
            "stage04",
        ),
        5: StageConfig(
            5,
            "stage_05_dogbone_damage",
            "Stage 05 Notched Dogbone Damage",
            "PASS",
            "mesh-ladder pre-onset RMSE: 4.843% / 3.736%",
            ("runs/stage05_medium/stage05_mediumA*.vtk",),
            "plastic_strain_percent",
            "stage05",
        ),
        6: StageConfig(
            6,
            "stage_06_composite_failure_criteria",
            "Stage 06 Composite Failure Cards",
            "PASS",
            "LAW12 + TYPE6 registers TSAIWU, HASHIN, and PUCK",
            (
                "runs/stage06_TSAIWU_law12_type6_pubviz/stage06_TSAIWU_law12_type6_pubvizA*.vtk",
                "runs/stage06_TSAIWU_law12_type6A*.vtk",
            ),
            "axial_stress_mpa",
            "stage06",
            -1,
            1.08,
        ),
        7: StageConfig(
            7,
            "stage_07_UD_tow_D3039",
            "Stage 07 UD Tow D3039",
            "PASS",
            "0 deg error 0.403%; 45 deg error 0.114%",
            ("runs/stage07_7B_law12_type6A*.vtk", "runs/stage07_7A_law12_type6A*.vtk"),
            "axial_stress_mpa",
            "stage07",
            -1,
            0.38,
        ),
        8: StageConfig(
            8,
            "stage_08_ply_rotation",
            "Stage 08 Ply Rotation Sweep",
            "PASS",
            "max Ex(theta) error: 0.815% over seven angles",
            ("runs/stage08_theta_45_law12_type6A*.vtk", "runs/stage08_theta_00_law12_type6A*.vtk"),
            "axial_stress_mpa",
            "stage08",
        ),
        9: StageConfig(
            9,
            "stage_09_laminate_solid_CLT",
            "Stage 09 Solid Laminate CLT",
            "PASS",
            "max A-matrix component error: 1.085%",
            ("runs/stage09_B_quasiiso_exx_law12_type6A*.vtk", "runs/stage09_A_crossply_exx_law12_type6A*.vtk"),
            "axial_stress_mpa",
            "stage09",
        ),
        10: StageConfig(
            10,
            "stage_10_UD_mesoscale_direct",
            "Stage 10 UD Mesoscale Direct",
            "INCONCLUSIVE",
            "KUBC bias: transverse error 11.189%; shear error 99.934%",
            ("runs/stage10_axial_zA*.vtk",),
            "axial_stress_mpa",
            "stage10",
        ),
        11: StageConfig(
            11,
            "stage_11_PW_mesoscale_direct",
            "Stage 11 AP-PLY Kok Block",
            "FAIL",
            "Ex 16.706 GPa vs target 53.3 GPa; Gxy 0.0229 GPa",
            ("runs/stage11_axial_xA*.vtk",),
            "axial_stress_mpa",
            "stage11",
        ),
    }


def vtk_paths(config: StageConfig) -> list[pathlib.Path]:
    for pattern in config.vtk_globs:
        paths = sorted(config.stage_dir.glob(pattern))
        if paths:
            return paths
    raise FileNotFoundError(f"no VTK paths for stage {config.stage:02d}")


def field_time(grid: pv.DataSet) -> float:
    try:
        return float(np.ravel(grid.field_data["TIME"])[0])
    except Exception:
        return 0.0


def tensor_stack(data: pv.DataSet, prefix: str) -> np.ndarray | None:
    tensors: list[np.ndarray] = []
    for name in data.cell_data.keys():
        if not name.startswith(prefix):
            continue
        arr = np.asarray(data.cell_data[name], dtype=float)
        if arr.ndim == 3 and arr.shape[1:] == (3, 3):
            tensors.append(arr)
        elif arr.ndim == 2 and arr.shape[1] == 9:
            tensors.append(arr.reshape((-1, 3, 3)))
        elif arr.ndim == 2 and arr.shape[1] == 6:
            xx, yy, zz, xy, yz, zx = arr.T
            tensors.append(
                np.stack(
                    [
                        np.stack([xx, xy, zx], axis=1),
                        np.stack([xy, yy, yz], axis=1),
                        np.stack([zx, yz, zz], axis=1),
                    ],
                    axis=1,
                )
            )
    if not tensors:
        return None
    return np.mean(np.stack(tensors, axis=0), axis=0)


def von_mises_from_tensor(tensor: np.ndarray) -> np.ndarray:
    trace = np.trace(tensor, axis1=1, axis2=2) / 3.0
    dev = tensor.copy()
    for i in range(3):
        dev[:, i, i] -= trace
    return np.sqrt(1.5 * np.sum(dev * dev, axis=(1, 2)))


def eq_strain_from_tensor(tensor: np.ndarray) -> np.ndarray:
    trace = np.trace(tensor, axis1=1, axis2=2) / 3.0
    dev = tensor.copy()
    for i in range(3):
        dev[:, i, i] -= trace
    return np.sqrt((2.0 / 3.0) * np.sum(dev * dev, axis=(1, 2)))


def scale_output_fields(grid: pv.DataSet, scale: float) -> None:
    if math.isclose(scale, 1.0):
        return
    for name in list(grid.point_data.keys()):
        if name in {"Displacement", "Velocity", "Reaction_Forces"}:
            grid.point_data[name] = np.asarray(grid.point_data[name], dtype=float) * scale
    for name in list(grid.cell_data.keys()):
        if "Strs" in name or "Stra" in name or "Plastic" in name or "EPSP" in name or "Damage" in name:
            grid.cell_data[name] = np.asarray(grid.cell_data[name], dtype=float) * scale


def add_fields(grid: pv.DataSet) -> None:
    disp = np.asarray(grid.point_data.get("Displacement", np.zeros((grid.n_points, 3))), dtype=float)
    if disp.ndim == 1:
        disp = np.zeros((grid.n_points, 3), dtype=float)
    grid.point_data["displacement_magnitude_mm"] = np.linalg.norm(disp, axis=1) * 1000.0
    stress = tensor_stack(grid, "3DELEM_Strs")
    if stress is not None:
        grid.cell_data["von_mises_mpa"] = von_mises_from_tensor(stress) / 1.0e6
        grid.cell_data["axial_stress_mpa"] = stress[:, 0, 0] / 1.0e6
    strain = tensor_stack(grid, "3DELEM_Stra")
    if strain is not None:
        grid.cell_data["eq_strain_percent"] = eq_strain_from_tensor(strain) * 100.0
    for name in list(grid.cell_data.keys()):
        if "Plastic" in name or "EPSP" in name:
            arr = np.asarray(grid.cell_data[name], dtype=float)
            grid.cell_data["plastic_strain_percent"] = (arr if arr.ndim == 1 else np.linalg.norm(arr, axis=1)) * 100.0
            break
    for name in list(grid.cell_data.keys()):
        if "Damage" in name:
            arr = np.asarray(grid.cell_data[name], dtype=float)
            grid.cell_data["damage_magnitude"] = arr if arr.ndim == 1 else np.linalg.norm(arr, axis=1)
            break


def prepare_grid(path: pathlib.Path, warp_factor: float, scale: float = 1.0) -> pv.DataSet:
    grid = pv.read(path).copy(deep=True)
    scale_output_fields(grid, scale)
    disp = np.asarray(grid.point_data.get("Displacement", np.zeros((grid.n_points, 3))), dtype=float)
    if disp.ndim == 2 and disp.shape[1] == 3:
        grid.points = np.asarray(grid.points, dtype=float) - disp
    add_fields(grid)
    if disp.ndim == 2 and disp.shape[1] == 3 and warp_factor:
        warped = grid.warp_by_vector("Displacement", factor=warp_factor)
        for key in grid.point_data:
            warped.point_data[key] = grid.point_data[key]
        for key in grid.cell_data:
            warped.cell_data[key] = grid.cell_data[key]
        return warped
    return grid


def reference_bounds(paths: Iterable[pathlib.Path]) -> tuple[float, float, float, float, float, float]:
    bounds = []
    for path in paths:
        grid = pv.read(path).copy(deep=True)
        disp = np.asarray(grid.point_data.get("Displacement", np.zeros((grid.n_points, 3))), dtype=float)
        if disp.ndim == 2 and disp.shape[1] == 3:
            grid.points = np.asarray(grid.points, dtype=float) - disp
        bounds.append(tuple(float(v) for v in grid.bounds))
    return (
        min(b[0] for b in bounds),
        max(b[1] for b in bounds),
        min(b[2] for b in bounds),
        max(b[3] for b in bounds),
        min(b[4] for b in bounds),
        max(b[5] for b in bounds),
    )


def characteristic_length(bounds: tuple[float, float, float, float, float, float]) -> float:
    return math.sqrt((bounds[1] - bounds[0]) ** 2 + (bounds[3] - bounds[2]) ** 2 + (bounds[5] - bounds[4]) ** 2)


def automatic_warp(paths: list[pathlib.Path], bounds: tuple[float, float, float, float, float, float]) -> float:
    max_disp = 0.0
    for path in paths:
        grid = pv.read(path)
        disp = np.asarray(grid.point_data.get("Displacement", np.zeros((grid.n_points, 3))), dtype=float)
        if disp.ndim == 2 and disp.shape[1] == 3 and disp.size:
            max_disp = max(max_disp, float(np.nanmax(np.linalg.norm(disp, axis=1))))
    if max_disp <= 0:
        return 0.0
    diag = characteristic_length(bounds) or 1.0
    if max_disp > 0.05 * diag:
        return 1.0
    return float(min(max(0.08 * diag / max_disp, 1.0), 60.0))


def camera_for(bounds: tuple[float, float, float, float, float, float]):
    center = np.array([(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2], dtype=float)
    radius = characteristic_length(bounds) or 1.0
    direction = np.array([1.65, -1.55, 1.05], dtype=float)
    direction /= np.linalg.norm(direction)
    pos = center + direction * radius * 1.75
    return tuple(pos.tolist()), tuple(center.tolist()), (0.0, 0.0, 1.0)


def scalar_range(paths: Iterable[pathlib.Path], scalar: str, location: str, warp: float) -> tuple[float, float]:
    chunks: list[np.ndarray] = []
    for path in paths:
        mesh = prepare_grid(path, warp)
        if location == "point" and scalar in mesh.point_data:
            chunks.append(np.asarray(mesh.point_data[scalar], dtype=float).ravel())
        elif location == "cell" and scalar in mesh.cell_data:
            chunks.append(np.asarray(mesh.cell_data[scalar], dtype=float).ravel())
    if not chunks:
        return (0.0, 1.0)
    values = np.concatenate(chunks)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (0.0, 1.0)
    if scalar in {"axial_stress_mpa"}:
        low = float(np.nanpercentile(values, 0.5))
        high = float(np.nanpercentile(values, 99.5))
    else:
        low = min(0.0, float(np.nanpercentile(values, 0.5)))
        high = float(np.nanpercentile(values, 99.5))
    if not math.isfinite(low):
        low = 0.0
    if not math.isfinite(high) or high <= low:
        high = float(np.nanmax(values)) if values.size else 1.0
    if not math.isfinite(high) or high <= low:
        high = low + 1.0
    return (low, high)


def available_fields(path: pathlib.Path) -> dict[str, FieldChoice]:
    mesh = prepare_grid(path, 0.0)
    fields = {
        "displacement_magnitude_mm": FieldChoice("displacement_magnitude_mm", "point", "Displacement magnitude", "mm"),
    }
    if "von_mises_mpa" in mesh.cell_data:
        fields["von_mises_mpa"] = FieldChoice("von_mises_mpa", "cell", "von Mises stress", "MPa")
    if "eq_strain_percent" in mesh.cell_data:
        fields["eq_strain_percent"] = FieldChoice("eq_strain_percent", "cell", "Equivalent strain", "%")
    if "plastic_strain_percent" in mesh.cell_data:
        fields["plastic_strain_percent"] = FieldChoice("plastic_strain_percent", "cell", "Plastic strain", "%")
    if "damage_magnitude" in mesh.cell_data:
        fields["damage_magnitude"] = FieldChoice("damage_magnitude", "cell", "Damage magnitude", "-")
    if "axial_stress_mpa" in mesh.cell_data:
        fields["axial_stress_mpa"] = FieldChoice("axial_stress_mpa", "cell", "Axial/fiber stress", "MPa")
    return fields


def scalar_bar(title: str, vertical: bool = True) -> dict[str, object]:
    return {
        "title": title,
        "title_font_size": 15,
        "label_font_size": 13,
        "color": CHARCOAL,
        "vertical": vertical,
        "fmt": "%.3g",
        "shadow": False,
        "n_labels": 5,
    }


def bar_label(field: FieldChoice) -> str:
    labels = {
        "displacement_magnitude_mm": "disp (mm)",
        "von_mises_mpa": "vM (MPa)",
        "eq_strain_percent": "eq strain (%)",
        "plastic_strain_percent": "plastic (%)",
        "damage_magnitude": "damage (-)",
        "axial_stress_mpa": "axial (MPa)",
    }
    return labels.get(field.scalar, f"{field.title} ({field.units})")


def setup_plotter(plotter: pv.Plotter, camera, bounds: tuple[float, float, float, float, float, float], view_scale: float = 0.58) -> None:
    plotter.set_background(WHITE)
    plotter.camera_position = camera
    try:
        plotter.camera.parallel_projection = True
        plotter.camera.parallel_scale = view_scale * max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4], 1.0e-9)
    except Exception:
        pass
    plotter.add_axes(color=CHARCOAL, line_width=2)
    plotter.camera.zoom(1.0)


def render_field_png(
    path: pathlib.Path,
    field: FieldChoice,
    clim: tuple[float, float],
    camera,
    bounds: tuple[float, float, float, float, float, float],
    warp: float,
    out_png: pathlib.Path,
    size: tuple[int, int] = PANEL_SIZE,
    scale: float = 1.0,
    view_scale: float = 0.58,
) -> None:
    mesh = prepare_grid(path, warp, scale)
    plotter = pv.Plotter(off_screen=True, window_size=size)
    setup_plotter(plotter, camera, bounds, view_scale)
    show_edges = mesh.n_cells <= 2500
    plotter.add_mesh(
        mesh,
        scalars=field.scalar,
        preference=field.location,
        cmap=BRAND_RAMP,
        clim=clim,
        show_edges=show_edges,
        edge_color=CHARCOAL,
        line_width=0.35 if show_edges else 0.0,
        scalar_bar_args=scalar_bar(bar_label(field)),
    )
    plotter.screenshot(str(out_png))
    plotter.close()


def render_wireframe_png(config: StageConfig, out_png: pathlib.Path) -> None:
    paths = vtk_paths(config)
    bounds = reference_bounds(paths)
    camera = camera_for(bounds)
    mesh = prepare_grid(paths[-1], 0.0)
    plotter = pv.Plotter(off_screen=True, window_size=(1800, 1200))
    setup_plotter(plotter, camera, bounds, config.view_scale)
    if mesh.n_cells <= 5000:
        plotter.add_mesh(mesh, color=BLACK30, show_edges=True, edge_color=GARNET, line_width=1.0)
    else:
        plotter.add_mesh(mesh, style="wireframe", color=GARNET, line_width=0.55)
    plotter.reset_camera_clipping_range()
    plotter.add_text(config.title, position="upper_left", font_size=18, color=CHARCOAL, font="courier")
    plotter.add_text(config.key_result, position="lower_left", font_size=12, color=CHARCOAL, font="courier")
    plotter.screenshot(str(out_png))
    plotter.close()


def read_csv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def fnum(value: str | float | int | None, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    return float(value)


def typst_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def fmt_num(value: float) -> str:
    if abs(value) >= 1000 or (0 < abs(value) < 0.01):
        return f"{value:.3g}"
    return f"{value:.4g}"


def color_name(hex_value: str) -> str:
    return {
        GARNET: "garnet",
        ATLANTIC: "atlantic",
        HORSESHOE: "horseshoe",
        CONGAREE: "congaree",
        HONEYCOMB: "honeycomb",
        ROSE: "rose",
        BLACK50: "black50",
    }.get(hex_value, "garnet")


def chart_header(title: str) -> list[str]:
    return [
        '#import "@preview/cetz:0.3.4"',
        '#set page(width: 146mm, height: 96mm, margin: 4mm)',
        '#set text(font: "Libertinus Serif", size: 8.5pt, fill: rgb("#363636"))',
        f'#let garnet = rgb("{GARNET}")',
        f'#let atlantic = rgb("{ATLANTIC}")',
        f'#let horseshoe = rgb("{HORSESHOE}")',
        f'#let congaree = rgb("{CONGAREE}")',
        f'#let honeycomb = rgb("{HONEYCOMB}")',
        f'#let rose = rgb("{ROSE}")',
        f'#let charcoal = rgb("{CHARCOAL}")',
        f'#let black10 = rgb("{BLACK10}")',
        f'#let black50 = rgb("{BLACK50}")',
        f'#let white = rgb("{WHITE}")',
        f'#align(center)[#text(size: 10.5pt, weight: "bold")[{typst_escape(title)}]]',
        '#v(1mm)',
        '#cetz.canvas(length: 1cm, {',
        '  import cetz.draw: *',
    ]


def chart_footer(note: str) -> list[str]:
    return [
        '})',
        '#v(0.5mm)',
        f'#text(size: 7pt, fill: rgb("{BLACK70}"))[{typst_escape(note)}]',
        '',
    ]


def map_point(x: float, y: float, xlim: tuple[float, float], ylim: tuple[float, float], rect: tuple[float, float, float, float]) -> tuple[float, float]:
    x0, y0, w, h = rect
    xmin, xmax = xlim
    ymin, ymax = ylim
    px = x0 + (x - xmin) / (xmax - xmin) * w if xmax != xmin else x0 + w / 2
    py = y0 + (y - ymin) / (ymax - ymin) * h if ymax != ymin else y0 + h / 2
    return px, py


def write_line_chart(
    out_typ: pathlib.Path,
    title: str,
    xlabel: str,
    ylabel: str,
    series: list[tuple[str, list[tuple[float, float]], str]],
    note: str,
    hlines: list[tuple[float, str, str]] | None = None,
) -> pathlib.Path:
    hlines = hlines or []
    all_points = [pt for _label, pts, _color in series for pt in pts]
    if not all_points:
        all_points = [(0.0, 0.0), (1.0, 1.0)]
    xs = [p[0] for p in all_points]
    ys = [p[1] for p in all_points] + [h[0] for h in hlines]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(0.0, min(ys)), max(ys)
    if math.isclose(xmin, xmax):
        xmin -= 0.5
        xmax += 0.5
    if math.isclose(ymin, ymax):
        ymax = ymin + 1.0
    ypad = 0.08 * (ymax - ymin)
    ymax += ypad
    rect = (1.35, 0.95, 11.7, 5.35)
    lines = chart_header(title)
    lines.append(f'  rect(({rect[0]:.3f}, {rect[1]:.3f}), ({rect[0]+rect[2]:.3f}, {rect[1]+rect[3]:.3f}), fill: white, stroke: charcoal + 0.65pt)')
    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        y = rect[1] + rect[3] * frac
        val = ymin + (ymax - ymin) * frac
        lines.append(f'  line(({rect[0]:.3f}, {y:.3f}), ({rect[0]+rect[2]:.3f}, {y:.3f}), stroke: black10 + 0.45pt)')
        lines.append(f'  content(({rect[0]-0.16:.3f}, {y:.3f}), [{fmt_num(val)}], anchor: "east")')
    for value, label, color in hlines:
        _x, y = map_point(xmin, value, (xmin, xmax), (ymin, ymax), rect)
        lines.append(f'  line(({rect[0]:.3f}, {y:.3f}), ({rect[0]+rect[2]:.3f}, {y:.3f}), stroke: {color_name(color)} + 0.75pt)')
        lines.append(f'  content(({rect[0]+rect[2]-0.1:.3f}, {y+0.18:.3f}), [{typst_escape(label)}], anchor: "east")')
    for label, pts, color in series:
        if not pts:
            continue
        mapped = [map_point(x, y, (xmin, xmax), (ymin, ymax), rect) for x, y in pts]
        point_text = ", ".join(f"({x:.3f}, {y:.3f})" for x, y in mapped)
        if len(mapped) >= 2:
            lines.append(f"  line({point_text}, stroke: {color_name(color)} + 1.05pt)")
        for x, y in mapped:
            lines.append(f"  circle(({x:.3f}, {y:.3f}), radius: 0.045, fill: {color_name(color)}, stroke: none)")
    legend_x = rect[0] + 0.1
    legend_y = rect[1] + rect[3] + 0.45
    for idx, (label, _pts, color) in enumerate(series[:4]):
        x = legend_x + idx * 3.1
        lines.append(f"  line(({x:.3f}, {legend_y:.3f}), ({x+0.42:.3f}, {legend_y:.3f}), stroke: {color_name(color)} + 1.05pt)")
        lines.append(f'  content(({x+0.52:.3f}, {legend_y:.3f}), [{typst_escape(label)}], anchor: "west")')
    lines.extend(
        [
            f'  content(({rect[0]:.3f}, {rect[1]-0.35:.3f}), [{fmt_num(xmin)}], anchor: "north")',
            f'  content(({rect[0]+rect[2]:.3f}, {rect[1]-0.35:.3f}), [{fmt_num(xmax)}], anchor: "north")',
            f'  content(({rect[0]+rect[2]/2:.3f}, {rect[1]-0.65:.3f}), [{typst_escape(xlabel)}], anchor: "north")',
            f'  content(({rect[0]-0.95:.3f}, {rect[1]+rect[3]/2:.3f}), [{typst_escape(ylabel)}], angle: 90deg)',
        ]
    )
    lines.extend(chart_footer(note))
    out_typ.write_text("\n".join(lines), encoding="utf-8")
    return out_typ


def write_bar_chart(
    out_typ: pathlib.Path,
    title: str,
    ylabel: str,
    bars: list[tuple[str, float, str]],
    note: str,
    hline: tuple[float, str, str] | None = None,
) -> pathlib.Path:
    ymax = max([v for _label, v, _color in bars] + ([hline[0]] if hline else [1.0]))
    ymax = ymax * 1.18 if ymax else 1.0
    rect = (1.25, 0.95, 11.8, 5.35)
    lines = chart_header(title)
    lines.append(f'  rect(({rect[0]:.3f}, {rect[1]:.3f}), ({rect[0]+rect[2]:.3f}, {rect[1]+rect[3]:.3f}), fill: white, stroke: charcoal + 0.65pt)')
    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        y = rect[1] + rect[3] * frac
        val = ymax * frac
        lines.append(f'  line(({rect[0]:.3f}, {y:.3f}), ({rect[0]+rect[2]:.3f}, {y:.3f}), stroke: black10 + 0.45pt)')
        lines.append(f'  content(({rect[0]-0.14:.3f}, {y:.3f}), [{fmt_num(val)}], anchor: "east")')
    if hline:
        value, label, color = hline
        y = rect[1] + (value / ymax) * rect[3]
        lines.append(f'  line(({rect[0]:.3f}, {y:.3f}), ({rect[0]+rect[2]:.3f}, {y:.3f}), stroke: {color_name(color)} + 0.8pt)')
        lines.append(f'  content(({rect[0]+0.15:.3f}, {y+0.2:.3f}), [{typst_escape(label)}], anchor: "west")')
    gap = rect[2] / max(len(bars), 1)
    barw = gap * 0.58
    for idx, (label, value, color) in enumerate(bars):
        cx = rect[0] + gap * (idx + 0.5)
        y = rect[1] + (value / ymax) * rect[3]
        lines.append(f'  rect(({cx-barw/2:.3f}, {rect[1]:.3f}), ({cx+barw/2:.3f}, {y:.3f}), fill: {color_name(color)}, stroke: none)')
        lines.append(f'  line(({cx-barw/2:.3f}, {y:.3f}), ({cx+barw/2:.3f}, {y:.3f}), stroke: charcoal + 0.5pt)')
        if len(bars) <= 10 or idx % 2 == 0:
            lines.append(f'  content(({cx:.3f}, {rect[1]-0.27:.3f}), [{typst_escape(label)}], anchor: "north")')
    lines.append(f'  content(({rect[0]-0.95:.3f}, {rect[1]+rect[3]/2:.3f}), [{typst_escape(ylabel)}], angle: 90deg)')
    lines.extend(chart_footer(note))
    out_typ.write_text("\n".join(lines), encoding="utf-8")
    return out_typ


def compile_typst_png(typ_path: pathlib.Path) -> pathlib.Path:
    png_path = typ_path.with_suffix(".png")
    subprocess.run(["typst", "compile", typ_path.name, png_path.name], cwd=str(typ_path.parent), check=True)
    return png_path


def write_history_chart(config: StageConfig) -> pathlib.Path:
    out_typ = config.figures_dir / f"{config.prefix}_history_pubviz.typ"
    if config.chart_kind == "stage01":
        rows = read_csv(config.stage_dir / "results" / "timeseries.csv")
        bars = [(f"{r['load_case']} {r['mesh']}", 100 * fnum(r["rel_error_eb"]), HORSESHOE if r["verdict"] == "PASS" else GARNET) for r in rows]
        typ = write_bar_chart(out_typ, "Beam Deflection Error Gate", "relative error (%)", bars, "Euler-Bernoulli reference; 1% gated M1 tolerance.", (1.0, "1% gate", HORSESHOE))
    elif config.chart_kind == "stage02":
        rows = read_csv(config.stage_dir / "results" / "timeseries.csv")
        r = rows[0]
        series = [
            ("OpenRadioss", [(fnum(r["alpha"]), fnum(r["dy_fem_over_L"]))], GARNET),
            ("Elastica", [(fnum(r["alpha"]), fnum(r["dy_ref_over_L"]))], ATLANTIC),
        ]
        typ = write_line_chart(out_typ, "Cantilever Tip Sag", "dimensionless load alpha", "tip sag dy/L", series, "alpha=1 passes; alpha=3 and 5 remain convergence-blocked.")
    elif config.chart_kind == "stage03":
        rows = read_csv(config.stage_dir / "results" / "timeseries.csv")
        fem = [(100 * fnum(r["eps_disp"]), fnum(r["sigma_mean_Pa"]) / 1.0e6) for r in rows]
        ref = [(100 * fnum(r["eps_disp"]), fnum(r["sigma_ref_Pa"]) / 1.0e6) for r in rows]
        typ = write_line_chart(out_typ, "Dogbone Stress-Strain", "gauge strain (%)", "stress (MPa)", [("OpenRadioss", fem, GARNET), ("load/area", ref, ATLANTIC)], "Medium mesh load-factor sweep.")
    elif config.chart_kind == "stage04":
        rows = read_csv(config.stage_dir / "results" / "timeseries.csv")
        fem = [(fnum(r["n_theta"]), fnum(r["kt_fem"])) for r in rows]
        xs = [x for x, _y in fem]
        howland = [(x, 3.035) for x in xs]
        kirsch = [(x, 3.0) for x in xs]
        typ = write_line_chart(out_typ, "Open-Hole Stress Concentration", "circumferential divisions", "Kt", [("OpenRadioss", fem, GARNET), ("Howland", howland, ATLANTIC), ("Kirsch", kirsch, HORSESHOE)], "Ntheta 64 is the primary Howland gate.")
    elif config.chart_kind == "stage05":
        rows = read_csv(config.stage_dir / "results" / "timeseries.csv")
        series = []
        for key, label, color in [("F_coarse_N", "coarse", ATLANTIC), ("F_medium_N", "medium", GARNET), ("F_fine_N", "fine", CONGAREE)]:
            pts = [(1000 * fnum(r["displacement_m"]), fnum(r[key]) / 1000.0) for r in rows if r.get(key)]
            series.append((label, pts, color))
        typ = write_line_chart(out_typ, "Damage Mesh-Ladder Response", "loaded-end displacement (mm)", "reaction force (kN)", series, "Pre-onset RMSE satisfies the 5% mesh-objectivity gate.")
    elif config.chart_kind == "stage06":
        rows = [r for r in read_csv(config.stage_dir / "results" / "timeseries.csv") if r.get("kind") == "analytic_envelope"]
        series = []
        for criterion, color in [("TSAIWU", GARNET), ("HASHIN", ATLANTIC), ("PUCK", CONGAREE)]:
            pts = [(fnum(r["path_id"]), fnum(r["reported_pa"]) / 1.0e6) for r in rows if r.get("criterion") == criterion]
            series.append((criterion, pts, color))
        typ = write_line_chart(out_typ, "Failure Envelope Probe", "analytic path id", "reported strength (MPa)", series, "LAW12 + TYPE6 overdrive decks register the failure cards.")
    elif config.chart_kind == "stage07":
        rows = read_csv(config.stage_dir / "results" / "timeseries.csv")
        solver = [(fnum(r["theta_deg"]), fnum(r["solver_modulus_pa"]) / 1.0e9) for r in rows]
        ref = [(fnum(r["theta_deg"]), fnum(r["reference_modulus_pa"]) / 1.0e9) for r in rows]
        typ = write_line_chart(out_typ, "UD Coupon Modulus", "fiber angle (deg)", "Ex (GPa)", [("OpenRadioss", solver, GARNET), ("analytic", ref, ATLANTIC)], "Uses displacement-BC engineering strain, not VTK Stra[0].")
    elif config.chart_kind == "stage08":
        rows = read_csv(config.stage_dir / "results" / "timeseries.csv")
        solver = [(fnum(r["theta_deg"]), fnum(r["solver_Ex_Pa"]) / 1.0e9) for r in rows]
        ref = [(fnum(r["theta_deg"]), fnum(r["analytic_Ex_Pa"]) / 1.0e9) for r in rows]
        typ = write_line_chart(out_typ, "Ply Rotation Modulus Sweep", "theta (deg)", "Ex (GPa)", [("OpenRadioss", solver, GARNET), ("analytic", ref, ATLANTIC)], "Maximum relative error is 0.815%.")
    elif config.chart_kind == "stage09":
        rows = [r for r in read_csv(config.stage_dir / "results" / "timeseries.csv") if r.get("component", "").startswith("A_")]
        bars = []
        for r in rows:
            label = r["component"].replace("A_", "")
            color = HORSESHOE if fnum(r["relative_error_pct"]) <= 2.0 else GARNET
            bars.append((label, fnum(r["relative_error_pct"]), color))
        typ = write_bar_chart(out_typ, "CLT A-Matrix Component Error", "relative error (%)", bars, "All reported A-matrix components are within the 2% gate.", (2.0, "2% gate", HORSESHOE))
    else:
        typ = write_bar_chart(out_typ, config.title, "value", [("n/a", 1.0, GARNET)], "")
    return compile_typst_png(typ)


def fit_image(src: Image.Image, size: tuple[int, int]) -> Image.Image:
    src = src.convert("RGB")
    src.thumbnail(size, Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", size, WHITE)
    canvas.paste(src, ((size[0] - src.width) // 2, (size[1] - src.height) // 2))
    return canvas


def draw_wrapped(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, font_obj, fill, max_width: int, line_gap: int = 5) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        test = word if not current else current + " " + word
        if draw.textbbox((0, 0), test, font=font_obj)[2] <= max_width:
            current = test
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=font_obj, fill=fill)
        y += draw.textbbox((0, 0), line, font=font_obj)[3] + line_gap
    return y


def verdict_color(verdict: str) -> str:
    if verdict == "PASS":
        return PASS_COLOR
    if verdict == "FAIL":
        return FAIL_COLOR
    return INCONCLUSIVE_COLOR


def compose_composite(config: StageConfig, panel_paths: list[pathlib.Path], chart_png: pathlib.Path, out_png: pathlib.Path) -> None:
    canvas = Image.new("RGB", COMPOSITE_SIZE, WHITE)
    draw = ImageDraw.Draw(canvas)
    draw.rectangle((0, 0, COMPOSITE_SIZE[0] - 1, COMPOSITE_SIZE[1] - 1), outline=hex_to_rgb(CHARCOAL), width=3)
    draw.text((58, 32), config.title, font=FONT_TITLE, fill=hex_to_rgb(CHARCOAL))
    draw.text((58, 84), "publication composite: deformed fields and gated history", font=FONT_SUBTITLE, fill=hex_to_rgb(BLACK70))
    verdict_fill = hex_to_rgb(verdict_color(config.verdict))
    draw.rectangle((2005, 40, 2315, 88), fill=verdict_fill)
    draw.text((2024, 49), config.verdict, font=FONT_PANEL, fill=hex_to_rgb(WHITE if config.verdict != "INCONCLUSIVE" else BLACK))

    positions = [(60, 145), (1230, 145), (60, 855), (1230, 855)]
    labels = ["Displacement magnitude (mm)", "von Mises stress (MPa)", "Stage-specific strain/stress field", "Typst + CeTZ gate history"]
    images = [Image.open(p) for p in panel_paths] + [Image.open(chart_png)]
    for (x, y), label, img in zip(positions, labels, images):
        draw.text((x, y - 36), label, font=FONT_PANEL, fill=hex_to_rgb(CHARCOAL))
        draw.rectangle((x - 1, y - 1, x + PANEL_SIZE[0], y + PANEL_SIZE[1]), outline=hex_to_rgb(BLACK30), width=2)
        canvas.paste(fit_image(img, PANEL_SIZE), (x, y))

    footer_y = 1512
    draw.rectangle((0, footer_y, COMPOSITE_SIZE[0], COMPOSITE_SIZE[1]), fill=hex_to_rgb(BLACK10))
    draw.rectangle((0, footer_y, COMPOSITE_SIZE[0] - 1, COMPOSITE_SIZE[1] - 1), outline=hex_to_rgb(CHARCOAL), width=2)
    draw.text((58, footer_y + 22), f"Stage {config.stage:02d}", font=FONT_PANEL, fill=hex_to_rgb(CHARCOAL))
    draw.rectangle((205, footer_y + 19, 430, footer_y + 58), fill=verdict_fill)
    draw.text((220, footer_y + 25), config.verdict, font=FONT_SMALL, fill=hex_to_rgb(WHITE if config.verdict != "INCONCLUSIVE" else BLACK))
    draw_wrapped(draw, (465, footer_y + 22), config.key_result, FONT_BODY, hex_to_rgb(CHARCOAL), 1830, 4)
    canvas.save(out_png)


def render_composite(config: StageConfig) -> pathlib.Path:
    pv.OFF_SCREEN = True
    config.figures_dir.mkdir(parents=True, exist_ok=True)
    paths = vtk_paths(config)
    final_path = paths[config.final_index]
    bounds = reference_bounds(paths)
    warp = automatic_warp(paths, bounds)
    camera = camera_for(bounds)
    fields = available_fields(final_path)
    disp_field = fields["displacement_magnitude_mm"]
    vm_field = fields.get("von_mises_mpa", disp_field)
    third = fields.get(config.third_field) or fields.get("eq_strain_percent") or fields.get("plastic_strain_percent") or vm_field
    selected = [disp_field, vm_field, third]
    ranges = {field.scalar: scalar_range(paths, field.scalar, field.location, warp) for field in selected}
    panel_paths: list[pathlib.Path] = []
    with tempfile.TemporaryDirectory(prefix=f"{config.prefix}_composite_", dir=str(config.figures_dir)) as tmp:
        tmpdir = pathlib.Path(tmp)
        for idx, field in enumerate(selected):
            out = tmpdir / f"panel_{idx}.png"
            render_field_png(final_path, field, ranges[field.scalar], camera, bounds, warp, out, view_scale=config.view_scale)
            panel_paths.append(out)
        chart_png = write_history_chart(config)
        out_png = config.figures_dir / f"{config.prefix}_composite.png"
        compose_composite(config, panel_paths, chart_png, out_png)
    print(f"rendered {out_png}")
    return out_png


def render_split_scene(
    frame: VideoFrame,
    config: StageConfig,
    disp_field: FieldChoice,
    vm_field: FieldChoice,
    ranges: dict[str, tuple[float, float]],
    camera,
    bounds,
    warp: float,
    out_png: pathlib.Path,
) -> float:
    mesh = prepare_grid(frame.path, warp, frame.scale)
    t = field_time(mesh) * frame.scale
    plotter = pv.Plotter(off_screen=True, window_size=(1920, 880), shape=(1, 2), border=False)
    for idx, field in enumerate([disp_field, vm_field]):
        plotter.subplot(0, idx)
        setup_plotter(plotter, camera, bounds, config.view_scale)
        show_edges = mesh.n_cells <= 2500
        plotter.add_mesh(
            mesh,
            scalars=field.scalar,
            preference=field.location,
            cmap=BRAND_RAMP,
            clim=ranges[field.scalar],
            show_edges=show_edges,
            edge_color=CHARCOAL,
            line_width=0.35 if show_edges else 0.0,
            scalar_bar_args=scalar_bar(bar_label(field)),
        )
        plotter.add_text(field.title, position="upper_left", font_size=17, color=CHARCOAL, font="courier")
    plotter.screenshot(str(out_png))
    plotter.close()
    return t


def title_card(config: StageConfig, text: str, subtext: str) -> Image.Image:
    canvas = Image.new("RGB", VIDEO_SIZE, WHITE)
    draw = ImageDraw.Draw(canvas)
    draw.rectangle((0, 0, VIDEO_SIZE[0] - 1, VIDEO_SIZE[1] - 1), outline=hex_to_rgb(CHARCOAL), width=4)
    draw.text((100, 155), config.title, font=font(58, True), fill=hex_to_rgb(CHARCOAL))
    draw.rectangle((100, 250, 520, 318), fill=hex_to_rgb(verdict_color(config.verdict)))
    draw.text((124, 267), config.verdict, font=font(32, True), fill=hex_to_rgb(WHITE if config.verdict != "INCONCLUSIVE" else BLACK))
    draw_wrapped(draw, (100, 375), text, font(38, False), hex_to_rgb(CHARCOAL), 1550, 10)
    draw_wrapped(draw, (100, 705), subtext, font(26, False), hex_to_rgb(BLACK70), 1550, 8)
    return canvas


def compose_video_frame(scene_png: pathlib.Path, config: StageConfig, t: float, progress: float) -> Image.Image:
    canvas = Image.new("RGB", VIDEO_SIZE, WHITE)
    draw = ImageDraw.Draw(canvas)
    scene = Image.open(scene_png).convert("RGB")
    canvas.paste(scene, (0, 74))
    draw.rectangle((0, 0, VIDEO_SIZE[0], 74), fill=hex_to_rgb(WHITE))
    draw.text((42, 22), config.title, font=FONT_SUBTITLE, fill=hex_to_rgb(CHARCOAL))
    draw.text((1510, 22), f"t = {t:.4g} s", font=FONT_SUBTITLE, fill=hex_to_rgb(BLACK70))
    draw.line((960, 92, 960, 934), fill=hex_to_rgb(BLACK30), width=2)
    y0 = 990
    draw.rectangle((0, 948, VIDEO_SIZE[0], VIDEO_SIZE[1]), fill=hex_to_rgb(BLACK10))
    draw.text((55, 966), config.key_result, font=FONT_SMALL, fill=hex_to_rgb(CHARCOAL))
    x1, x2 = 90, 1830
    draw.line((x1, y0 + 43, x2, y0 + 43), fill=hex_to_rgb(BLACK50), width=3)
    cx = int(x1 + max(0.0, min(1.0, progress)) * (x2 - x1))
    draw.rectangle((x1, y0 + 36, cx, y0 + 50), fill=hex_to_rgb(verdict_color(config.verdict)))
    draw.line((cx, y0 + 24, cx, y0 + 62), fill=hex_to_rgb(GARNET), width=4)
    draw.text((x1, y0 + 68), "start", font=FONT_SMALL, fill=hex_to_rgb(BLACK70))
    draw.text((x2 - 36, y0 + 68), "end", font=FONT_SMALL, fill=hex_to_rgb(BLACK70))
    return canvas


def selected_video_frames(paths: list[pathlib.Path]) -> list[VideoFrame]:
    if len(paths) <= 2:
        final_path = paths[-1]
        return [VideoFrame(final_path, i / 59.0) for i in range(60)]
    return [VideoFrame(path, 1.0) for path in paths]


def encode_video(frame_dir: pathlib.Path, out_mp4: pathlib.Path) -> None:
    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-framerate",
            str(VIDEO_FPS),
            "-i",
            str(frame_dir / "frame_%04d.png"),
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            "-crf",
            "18",
            "-movflags",
            "+faststart",
            str(out_mp4),
        ],
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )


def render_video(config: StageConfig) -> pathlib.Path:
    config.figures_dir.mkdir(parents=True, exist_ok=True)
    paths = vtk_paths(config)
    frames = selected_video_frames(paths)
    bounds = reference_bounds(paths)
    warp = automatic_warp(paths, bounds)
    camera = camera_for(bounds)
    fields = available_fields(paths[-1])
    disp_field = fields["displacement_magnitude_mm"]
    vm_field = fields.get("von_mises_mpa", disp_field)
    ranges = {
        disp_field.scalar: scalar_range(paths, disp_field.scalar, disp_field.location, warp),
        vm_field.scalar: scalar_range(paths, vm_field.scalar, vm_field.location, warp),
    }
    out_mp4 = config.figures_dir / f"{config.prefix}_motion_hd.mp4"
    with tempfile.TemporaryDirectory(prefix=f"{config.prefix}_motion_hd_", dir=str(config.figures_dir)) as tmp:
        tmpdir = pathlib.Path(tmp)
        unique: list[tuple[pathlib.Path, float]] = []
        for idx, frame_spec in enumerate(frames):
            scene = tmpdir / f"scene_{idx:04d}.png"
            t = render_split_scene(frame_spec, config, disp_field, vm_field, ranges, camera, bounds, warp, scene)
            unique.append((scene, t))
        frame_idx = 0
        start_card = title_card(config, config.key_result, "Split-screen displacement and von Mises fields; color ranges are fixed over the motion.")
        end_card = title_card(config, "Final verdict and key result", config.key_result)
        for _ in range(VIDEO_CARD_SECONDS * VIDEO_FPS):
            start_card.save(tmpdir / f"frame_{frame_idx:04d}.png")
            frame_idx += 1
        motion_count = VIDEO_MOTION_SECONDS * VIDEO_FPS
        for i in range(motion_count):
            uidx = round(i * (len(unique) - 1) / max(motion_count - 1, 1))
            scene, t = unique[uidx]
            compose_video_frame(scene, config, t, i / max(motion_count - 1, 1)).save(tmpdir / f"frame_{frame_idx:04d}.png")
            frame_idx += 1
        for _ in range(VIDEO_CARD_SECONDS * VIDEO_FPS):
            end_card.save(tmpdir / f"frame_{frame_idx:04d}.png")
            frame_idx += 1
        encode_video(tmpdir, out_mp4)
    print(f"rendered {out_mp4}")
    return out_mp4


def write_root_typst(path: pathlib.Path, body_lines: list[str], width_mm: int = 190, height_mm: int = 130) -> pathlib.Path:
    lines = [
        '#import "@preview/cetz:0.3.4"',
        f'#set page(width: {width_mm}mm, height: {height_mm}mm, margin: 7mm)',
        '#set text(font: "Libertinus Serif", size: 8pt, fill: rgb("#363636"))',
        f'#let garnet = rgb("{GARNET}")',
        f'#let horseshoe = rgb("{HORSESHOE}")',
        f'#let honeycomb = rgb("{HONEYCOMB}")',
        f'#let atlantic = rgb("{ATLANTIC}")',
        f'#let congaree = rgb("{CONGAREE}")',
        f'#let black10 = rgb("{BLACK10}")',
        f'#let black30 = rgb("{BLACK30}")',
        f'#let black50 = rgb("{BLACK50}")',
        f'#let charcoal = rgb("{CHARCOAL}")',
        f'#let white = rgb("{WHITE}")',
        "",
    ]
    lines.extend(body_lines)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def compile_typst_pdf(path: pathlib.Path) -> pathlib.Path:
    pdf = path.with_suffix(".pdf")
    subprocess.run(["typst", "compile", path.name, pdf.name], cwd=str(path.parent), check=True)
    return pdf


def summary_figures() -> list[pathlib.Path]:
    out_dir = ROOT / "figures"
    out_dir.mkdir(exist_ok=True)
    produced: list[pathlib.Path] = []
    stages = [
        (1, "PASS", "beam defl."),
        (2, "INCONCLUSIVE", "alpha 1 only"),
        (3, "PASS", "dogbone E8"),
        (4, "PASS", "open hole Kt"),
        (5, "PASS", "damage mesh"),
        (6, "PASS", "FAIL cards"),
        (7, "PASS", "UD coupon"),
        (8, "PASS", "ply angle"),
        (9, "PASS", "CLT A"),
        (10, "INCONCLUSIVE", "KUBC bias"),
        (11, "FAIL", "Kok gap"),
        (12, "INCONCLUSIVE", "cohesive start"),
        (13, "INCONCLUSIVE", "mesh pipe"),
        (14, "INCONCLUSIVE", "no state"),
        (15, "INCONCLUSIVE", "deck syntax"),
        (16, "INCONCLUSIVE", "compute cap"),
    ]
    lines = [
        '#align(center)[#text(size: 12pt, weight: "bold")[FEA_AP-PLY Final Stage Tally]]',
        '#v(2mm)',
        '#cetz.canvas(length: 1cm, {',
        '  import cetz.draw: *',
        '  let x0 = 0.55',
        '  let y0 = 6.1',
        '  let w = 0.92',
        '  let h = 0.58',
        '  content((0.2, 6.42), [result], anchor: "east")',
        '  content((0.2, 5.62), [criterion], anchor: "east")',
    ]
    color_map = {"PASS": "horseshoe", "INCONCLUSIVE": "honeycomb", "FAIL": "garnet"}
    for idx, (stage, verdict, criterion) in enumerate(stages):
        x = 0.55 + idx * 0.96
        color = color_map[verdict]
        short = "PASS" if verdict == "PASS" else ("FAIL" if verdict == "FAIL" else "INC")
        lines.extend(
            [
                f'  content(({x + 0.46:.3f}, 7.05), [{stage:02d}], anchor: "center")',
                f'  rect(({x:.3f}, 6.05), ({x+w:.3f}, 6.63), fill: {color}, stroke: charcoal + 0.35pt)',
                f'  content(({x+0.46:.3f}, 6.34), [{short}], anchor: "center")',
                f'  rect(({x:.3f}, 5.26), ({x+w:.3f}, 5.84), fill: white, stroke: black50 + 0.35pt)',
                f'  content(({x+0.46:.3f}, 5.55), [{typst_escape(criterion)}], anchor: "center")',
            ]
        )
    chips = [
        ("MUMPS rebuild", "atlantic"),
        ("LAW matrix", "garnet"),
        ("orientation recipe", "congaree"),
        ("Kok M1", "horseshoe"),
        ("Kok M2", "horseshoe"),
        ("Kok M3", "horseshoe"),
        ("Kok M4", "honeycomb"),
    ]
    lines.append('  content((0.55, 4.45), [methodology discoveries and Kok port], anchor: "west")')
    x = 0.55
    for label, color in chips:
        width = max(1.65, 0.12 * len(label))
        lines.extend(
            [
                f'  rect(({x:.3f}, 3.72), ({x+width:.3f}, 4.22), fill: {color}, stroke: charcoal + 0.35pt)',
                f'  content(({x+0.08:.3f}, 3.97), [{typst_escape(label)}], anchor: "west")',
            ]
        )
        x += width + 0.22
    lines.append('})')
    typ = write_root_typst(out_dir / "cross_stage_summary.typ", lines, 202, 118)
    produced.extend([typ, compile_typst_pdf(typ)])

    law_lines = [
        '#align(center)[#text(size: 11pt, weight: "bold")[Discovery: LAW x PROP Compatibility]]',
        '#v(1mm)',
        '#cetz.canvas(length: 1cm, {',
        '  import cetz.draw: *',
        '  let laws = ("LAW12", "LAW14", "LAW25", "LAW28", "LAW53", "LAW128")',
        '  let props = ("TYPE14", "TYPE6", "HASHIN", "PUCK", "TSAIWU")',
        '  for i in range(6) { content((0.9 + i * 1.45, 5.9), [#laws.at(i)], anchor: "center") }',
        '  for j in range(5) { content((0.05, 5.15 - j * 0.78), [#props.at(j)], anchor: "west") }',
        '  for i in range(6) {',
        '    rect((0.45 + i * 1.45, 4.85), (1.35 + i * 1.45, 5.45), fill: garnet, stroke: charcoal + 0.3pt)',
        '    content((0.90 + i * 1.45, 5.15), [B3047], anchor: "center")',
        '    for j in range(4) {',
        '      rect((0.45 + i * 1.45, 4.07 - j * 0.78), (1.35 + i * 1.45, 4.67 - j * 0.78), fill: horseshoe, stroke: charcoal + 0.3pt)',
        '      content((0.90 + i * 1.45, 4.37 - j * 0.78), [OK], anchor: "center")',
        '    }',
        '  }',
        '  content((0.55, 0.75), [Result: use LAW12 + TYPE6/SOL_ORTH for solid composite decks; TYPE14 remains blocked.], anchor: "west")',
        '})',
    ]
    typ = write_root_typst(out_dir / "discovery_law_matrix.typ", law_lines, 152, 90)
    produced.extend([typ, compile_typst_pdf(typ)])

    orientation_lines = [
        '#align(center)[#text(size: 11pt, weight: "bold")[Discovery: Orientation Convention]]',
        '#v(1mm)',
        '#cetz.canvas(length: 1cm, {',
        '  import cetz.draw: *',
        '  rect((0.7, 4.8), (3.9, 5.8), fill: horseshoe, stroke: charcoal + 0.4pt)',
        '  content((2.3, 5.3), [Uniform coupon: TYPE6 Phi = theta, Ip=3, Iorth=0], anchor: "center")',
        '  rect((5.0, 4.8), (8.2, 5.8), fill: honeycomb, stroke: charcoal + 0.4pt)',
        '  content((6.6, 5.3), [/SKEW with Ip=1: rejected recipe], anchor: "center")',
        '  rect((9.3, 4.8), (12.8, 5.8), fill: horseshoe, stroke: charcoal + 0.4pt)',
        '  content((11.05, 5.3), [Per-element tow: /INIBRI/ORTHO axes], anchor: "center")',
        '  line((3.9, 5.3), (5.0, 5.3), stroke: charcoal + 0.55pt)',
        '  line((8.2, 5.3), (9.3, 5.3), stroke: charcoal + 0.55pt)',
        '  rect((2.4, 2.8), (11.2, 3.65), fill: black10, stroke: charcoal + 0.4pt)',
        '  content((6.8, 3.23), [Gate measurement: sigma_x / displacement-BC engineering strain, not VTK Stra[0]], anchor: "center")',
        '  content((6.8, 1.75), [Verified maximum match: 0.24% across tested LAW12/LAW14 angle rows.], anchor: "center")',
        '})',
    ]
    typ = write_root_typst(out_dir / "discovery_orientation_convention.typ", orientation_lines, 170, 92)
    produced.extend([typ, compile_typst_pdf(typ)])

    kok_lines = [
        '#align(center)[#text(size: 11pt, weight: "bold")[Discovery: Kok Geometry Port]]',
        '#v(1mm)',
        '#cetz.canvas(length: 1cm, {',
        '  import cetz.draw: *',
        '  let labels = ("M1 tow", "M2 ply", "M3 laminate", "M4 CLI/export")',
        '  for i in range(4) {',
        '    let x = 1.2 + i * 3.1',
        '    rect((x, 4.6), (x + 2.25, 5.45), fill: if i < 3 { horseshoe } else { honeycomb }, stroke: charcoal + 0.4pt)',
        '    content((x + 1.125, 5.02), [#labels.at(i)], anchor: "center")',
        '    if i < 3 { line((x + 2.25, 5.02), (x + 3.1, 5.02), stroke: charcoal + 0.55pt) }',
        '  }',
        '  rect((1.2, 2.35), (10.85, 3.25), fill: black10, stroke: charcoal + 0.4pt)',
        '  line((1.6, 2.8), (4.2, 2.8), stroke: garnet + 1.4pt)',
        '  content((4.35, 2.8), [Stage 11 Ex = 16.7 GPa vs target 53.3 GPa], anchor: "west")',
        '  content((6.2, 1.6), [Port is runnable through OpenRadioss; current clean-room geometry lacks the published stiffness.], anchor: "center")',
        '})',
    ]
    typ = write_root_typst(out_dir / "discovery_kok_port.typ", kok_lines, 166, 90)
    produced.extend([typ, compile_typst_pdf(typ)])
    return produced


def render_stage(stage: int, include_video: bool) -> list[pathlib.Path]:
    configs = stage_configs()
    config = configs[stage]
    produced = []
    if stage <= 9:
        produced.append(render_composite(config))
    elif stage in {10, 11}:
        out = config.figures_dir / f"{config.prefix}_wireframe_pubviz.png"
        render_wireframe_png(config, out)
        produced.append(out)
        print(f"rendered {out}")
    if include_video and stage <= 6:
        produced.append(render_video(config))
    return produced


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=int, action="append")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--video", action="store_true", help="render HD video for selected stage(s) when applicable")
    parser.add_argument("--wireframes", action="store_true", help="render stage 10-11 wireframes")
    parser.add_argument("--summaries", action="store_true")
    args = parser.parse_args(argv)
    if shutil.which("typst") is None:
        raise SystemExit("typst not found")
    if shutil.which("ffmpeg") is None:
        raise SystemExit("ffmpeg not found")
    stages = args.stage or []
    if args.all:
        stages = list(range(1, 10))
    for stage in stages:
        render_stage(stage, args.video)
    if args.wireframes:
        for stage in [10, 11]:
            render_stage(stage, False)
    if args.summaries:
        for path in summary_figures():
            print(f"rendered {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
