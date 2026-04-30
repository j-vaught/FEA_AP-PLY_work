"""Render the Stage 16 Phase A geometry figure from the 200 mm M5 mesh."""

from __future__ import annotations

import math
import pathlib
import re
import tempfile

import meshio
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pyvista as pv
from vtk import VTK_TETRA


ROOT = pathlib.Path(__file__).resolve().parents[1]
STAGE_DIR = ROOT / "tests" / "stage_16_PW_panel_ballistic"
GEOM_DIR = STAGE_DIR / "geometry"
FIGURES_DIR = STAGE_DIR / "figures"
OUT_PNG = FIGURES_DIR / "stage_16_phase_a_geometry.png"

GARNET = "#73000A"
WHITE = "#FFFFFF"
CHARCOAL = "#363636"
BLACK70 = "#5C5C5C"
BLACK50 = "#A2A2A2"
BLACK30 = "#C7C7C7"
BLACK10 = "#ECECEC"
ATLANTIC = "#466A9F"
CONGAREE = "#1F414D"
HORSESHOE = "#65780B"
HONEYCOMB = "#A49137"

ANGLE_SEQUENCE = [45.0, 90.0, -45.0, 0.0]
ANGLE_LABELS = ["0 deg", "+45 deg", "-45 deg", "90 deg"]
ANGLE_COLORS = [GARNET, ATLANTIC, HORSESHOE, CONGAREE]
ANGLE_TO_CODE = {0.0: 0, 45.0: 1, -45.0: 2, 90.0: 3}

COMPOSITE_SIZE = (2400, 1600)
PANEL_VIEW_SIZE = (1110, 1180)
PATCH_VIEW_SIZE = (1110, 1180)


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


FONT_TITLE = font(38, True)
FONT_SUBTITLE = font(23, False)
FONT_PANEL = font(24, True)
FONT_BODY = font(18, False)


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def ply_index_from_name(name: str) -> int:
    match = re.search(r"PLY_(\d+)", name)
    if match is None:
        raise ValueError(f"could not extract ply index from {name}")
    return int(match.group(1))


def angle_for_ply(ply_index: int) -> float:
    return ANGLE_SEQUENCE[(ply_index - 1) % len(ANGLE_SEQUENCE)]


def build_grid() -> pv.UnstructuredGrid:
    mesh = meshio.read(GEOM_DIR / "p1_twt_200.msh")
    tetra10 = np.asarray(next(block.data for block in mesh.cells if block.type == "tetra10"), dtype=np.int64)
    physicals = np.asarray(mesh.cell_data["gmsh:physical"][0], dtype=np.int64)
    physical_to_name = {int(meta[0]): name for name, meta in mesh.field_data.items() if int(meta[1]) == 3}

    linear_tets = tetra10[:, :4]
    cell_array = np.hstack([np.full((linear_tets.shape[0], 1), 4, dtype=np.int64), linear_tets]).ravel()
    celltypes = np.full(linear_tets.shape[0], VTK_TETRA, dtype=np.uint8)
    grid = pv.UnstructuredGrid(cell_array, celltypes, np.asarray(mesh.points, dtype=float))

    angle_codes = np.empty(linear_tets.shape[0], dtype=np.int16)
    ply_indices = np.empty(linear_tets.shape[0], dtype=np.int16)
    for idx, physical_id in enumerate(physicals):
        name = physical_to_name[int(physical_id)]
        ply_index = ply_index_from_name(name)
        ply_indices[idx] = ply_index
        angle_codes[idx] = ANGLE_TO_CODE[angle_for_ply(ply_index)]
    grid.cell_data["angle_code"] = angle_codes
    grid.cell_data["ply_index"] = ply_indices
    return grid


def characteristic_length(bounds: tuple[float, float, float, float, float, float]) -> float:
    return math.sqrt((bounds[1] - bounds[0]) ** 2 + (bounds[3] - bounds[2]) ** 2 + (bounds[5] - bounds[4]) ** 2)


def camera_for(bounds: tuple[float, float, float, float, float, float], scale: float = 1.85, direction: tuple[float, float, float] = (1.6, -1.35, 1.05)):
    center = np.array([(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2], dtype=float)
    radius = characteristic_length(bounds) or 1.0
    direction_vec = np.array(direction, dtype=float)
    direction_vec /= np.linalg.norm(direction_vec)
    position = center + direction_vec * radius * scale
    return tuple(position.tolist()), tuple(center.tolist()), (0.0, 0.0, 1.0)


def setup_plotter(plotter: pv.Plotter, camera, bounds, parallel_scale_factor: float) -> None:
    plotter.set_background(WHITE)
    plotter.camera_position = camera
    plotter.camera.parallel_projection = True
    plotter.camera.parallel_scale = parallel_scale_factor * max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4], 1.0e-9)


def ply_surfaces(grid: pv.UnstructuredGrid, patch_mm: float | None = None, explode_gap_m: float = 0.0) -> list[tuple[pv.PolyData, int]]:
    bounds = tuple(float(value) for value in grid.bounds)
    clip_bounds = None
    if patch_mm is not None:
        half = patch_mm * 0.5e-3
        clip_bounds = (-half, half, -half, half, bounds[4] - 1.0e-6, bounds[5] + 1.0e-6)

    surfaces: list[tuple[pv.PolyData, int]] = []
    for ply_index in range(1, 25):
        cell_ids = np.where(np.asarray(grid.cell_data["ply_index"], dtype=np.int16) == ply_index)[0]
        subgrid = grid.extract_cells(cell_ids.tolist())
        surface = subgrid.extract_surface()
        if clip_bounds is not None:
            surface = surface.clip_box(bounds=clip_bounds, invert=False)
        if surface.n_points == 0:
            continue
        if explode_gap_m:
            surface.points = np.asarray(surface.points, dtype=float).copy()
            surface.points[:, 2] += (ply_index - 12.5) * explode_gap_m
        angle_code = ANGLE_TO_CODE[angle_for_ply(ply_index)]
        surfaces.append((surface, angle_code))
    return surfaces


def render_view(grid: pv.UnstructuredGrid, out_png: pathlib.Path, size: tuple[int, int], patch_mm: float | None = None) -> None:
    surfaces = ply_surfaces(grid, patch_mm=patch_mm, explode_gap_m=0.000035 if patch_mm is None else 0.00016)
    all_bounds = []
    for surface, _code in surfaces:
        all_bounds.append(tuple(float(value) for value in surface.bounds))
    bounds = (
        min(bound[0] for bound in all_bounds),
        max(bound[1] for bound in all_bounds),
        min(bound[2] for bound in all_bounds),
        max(bound[3] for bound in all_bounds),
        min(bound[4] for bound in all_bounds),
        max(bound[5] for bound in all_bounds),
    )
    if patch_mm is None:
        camera = camera_for(bounds, 1.85, (1.6, -1.35, 1.05))
        scale = 0.60
    else:
        camera = camera_for(bounds, 2.25, (1.15, -1.0, 0.82))
        scale = 0.82

    plotter = pv.Plotter(off_screen=True, window_size=size)
    setup_plotter(plotter, camera, bounds, scale)
    for surface, angle_code in surfaces:
        plotter.add_mesh(
            surface,
            color=ANGLE_COLORS[angle_code],
            show_edges=True,
            edge_color=BLACK30,
            line_width=0.25 if patch_mm is None else 0.45,
        )
    plotter.screenshot(str(out_png))
    plotter.close()


def fit_image(src: Image.Image, size: tuple[int, int]) -> Image.Image:
    image = src.convert("RGB")
    image.thumbnail(size, Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", size, WHITE)
    canvas.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return canvas


def add_legend(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    draw.text((x, y), "Angle classes", font=FONT_PANEL, fill=hex_to_rgb(CHARCOAL))
    for idx, (label, color) in enumerate(zip(ANGLE_LABELS, ANGLE_COLORS)):
        y0 = y + 42 + idx * 36
        draw.rectangle((x, y0, x + 24, y0 + 24), fill=hex_to_rgb(color), outline=hex_to_rgb(CHARCOAL), width=1)
        draw.text((x + 38, y0 - 1), label, font=FONT_BODY, fill=hex_to_rgb(CHARCOAL))


def compose(full_png: pathlib.Path, patch_png: pathlib.Path) -> pathlib.Path:
    canvas = Image.new("RGB", COMPOSITE_SIZE, WHITE)
    draw = ImageDraw.Draw(canvas)
    draw.rectangle((0, 0, COMPOSITE_SIZE[0] - 1, COMPOSITE_SIZE[1] - 1), outline=hex_to_rgb(CHARCOAL), width=3)
    draw.rectangle((0, 0, COMPOSITE_SIZE[0], 110), fill=hex_to_rgb(BLACK10))
    draw.text((54, 22), "Stage 16 - Phase A geometry", font=FONT_TITLE, fill=hex_to_rgb(CHARCOAL))
    draw.text((54, 66), "Exploded ply-stack surfaces from the 200 mm M5 mesh and center patch", font=FONT_SUBTITLE, fill=hex_to_rgb(BLACK70))

    draw.text((58, 126), "Full 200 mm panel mesh", font=FONT_PANEL, fill=hex_to_rgb(CHARCOAL))
    draw.text((1230, 126), "Impact-center 25 mm x 25 mm patch", font=FONT_PANEL, fill=hex_to_rgb(CHARCOAL))

    left = fit_image(Image.open(full_png), PANEL_VIEW_SIZE)
    right = fit_image(Image.open(patch_png), PATCH_VIEW_SIZE)
    canvas.paste(left, (50, 170))
    canvas.paste(right, (1230, 170))
    draw.rectangle((49, 169, 50 + PANEL_VIEW_SIZE[0], 170 + PANEL_VIEW_SIZE[1]), outline=hex_to_rgb(BLACK30), width=2)
    draw.rectangle((1229, 169, 1230 + PATCH_VIEW_SIZE[0], 170 + PATCH_VIEW_SIZE[1]), outline=hex_to_rgb(BLACK30), width=2)

    add_legend(draw, 58, 1340)
    add_legend(draw, 1230, 1340)

    footer_y = 1520
    draw.rectangle((0, footer_y, COMPOSITE_SIZE[0], COMPOSITE_SIZE[1]), fill=hex_to_rgb(BLACK10))
    draw.rectangle((0, footer_y, COMPOSITE_SIZE[0] - 1, COMPOSITE_SIZE[1] - 1), outline=hex_to_rgb(CHARCOAL), width=2)
    footer = "Stage 16 - Phase A - 200 mm M5 mesh, 1.44 M TETRA10, post-Kok-M5"
    draw.text((54, footer_y + 25), footer, font=FONT_PANEL, fill=hex_to_rgb(CHARCOAL))
    canvas.save(OUT_PNG)
    return OUT_PNG


def main() -> int:
    pv.OFF_SCREEN = True
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    grid = build_grid()
    with tempfile.TemporaryDirectory(prefix="stage16_phase_a_", dir=str(FIGURES_DIR)) as tmp:
        tmpdir = pathlib.Path(tmp)
        full_png = tmpdir / "full.png"
        patch_png = tmpdir / "patch.png"
        render_view(grid, full_png, PANEL_VIEW_SIZE)
        render_view(grid, patch_png, PATCH_VIEW_SIZE, patch_mm=25.0)
        out = compose(full_png, patch_png)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
