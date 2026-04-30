"""Publication-grade Stage 11 composite and motion render."""

from __future__ import annotations

import csv
import json
import math
import pathlib
import subprocess
import tempfile

import meshio
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pyvista as pv


ROOT = pathlib.Path(__file__).resolve().parents[1]
STAGE_DIR = ROOT / "tests" / "stage_11_PW_mesoscale_direct"
RUNS_DIR = STAGE_DIR / "runs"
FIGURES_DIR = STAGE_DIR / "figures"
RESULTS_DIR = STAGE_DIR / "results"
GEOM_DIR = STAGE_DIR / "geometry"

GARNET = "#73000A"
WHITE = "#FFFFFF"
BLACK = "#000000"
CHARCOAL = "#363636"
BLACK70 = "#5C5C5C"
BLACK50 = "#A2A2A2"
BLACK30 = "#C7C7C7"
BLACK10 = "#ECECEC"
ATLANTIC = "#466A9F"
CONGAREE = "#1F414D"
HORSESHOE = "#65780B"
HONEYCOMB = "#A49137"
SANDSTORM = "#FFF2E3"
ROSE = "#CC2E40"

FIELD_RAMP = [SANDSTORM, GARNET, ATLANTIC]
COMPOSITE_SIZE = (2400, 1600)
PANEL_SIZE = (1110, 640)
VIDEO_SIZE = (1920, 1080)
VIDEO_FPS = 30
VIDEO_CARD_SECONDS = 1
VIDEO_MOTION_SECONDS = 8
TITLE = "Stage 11 PW_mesoscale_direct \u2014 PASS \u2014 post-Kok-M5"


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


FONT_TITLE = font(40, True)
FONT_SUBTITLE = font(24, False)
FONT_PANEL = font(24, True)
FONT_BODY = font(20, False)
FONT_SMALL = font(18, False)
FONT_CARD = font(54, True)
FONT_CARD_BODY = font(30, False)


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def draw_wrapped(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, font_obj, fill, max_width: int, line_gap: int = 4) -> int:
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


def typst_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def vtk_series(case: str) -> list[pathlib.Path]:
    return sorted((RUNS_DIR / f"{case}_pubviz").glob(f"{case}_pubvizA*.vtk"))


def moduli_rows() -> list[dict[str, str]]:
    with (RESULTS_DIR / "effective_moduli.csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def measured_moduli() -> dict[str, float]:
    data = {}
    for row in moduli_rows():
        data[row["metric"]] = float(row["measured_GPa"])
    return data


def typst_chart() -> pathlib.Path:
    rows = moduli_rows()
    metrics = []
    for row in rows:
        metric = row["metric"].replace("_GPa", "").replace("_", "_")
        measured = float(row["measured_GPa"])
        target = float(row["target_GPa"])
        error_pct = float(row["relative_error_pct"])
        metrics.append((metric, measured, target, error_pct))

    typ = FIGURES_DIR / "stage11_effective_moduli.typ"
    lines = [
        '#import "@preview/cetz:0.3.4"',
        '#set page(width: 146mm, height: 96mm, margin: 4mm)',
        '#set text(font: "Libertinus Serif", size: 8.5pt, fill: rgb("#363636"))',
        f'#let garnet = rgb("{GARNET}")',
        f'#let atlantic = rgb("{ATLANTIC}")',
        f'#let horseshoe = rgb("{HORSESHOE}")',
        f'#let honeycomb = rgb("{HONEYCOMB}")',
        f'#let black10 = rgb("{BLACK10}")',
        f'#let black50 = rgb("{BLACK50}")',
        f'#let charcoal = rgb("{CHARCOAL}")',
        f'#let white = rgb("{WHITE}")',
        f'#align(center)[#text(size: 10.5pt, weight: "bold")[{typst_escape(TITLE)}]]',
        '#v(1mm)',
        '#cetz.canvas(length: 1cm, {',
        '  import cetz.draw: *',
        '  let x0 = 1.25',
        '  let y0 = 0.95',
        '  let w = 11.8',
        '  let h = 5.35',
        '  rect((x0, y0), (x0 + w, y0 + h), fill: white, stroke: charcoal + 0.65pt)',
    ]
    ymax = max(max(m, t) for _name, m, t, _e in metrics) * 1.18
    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        y = 0.95 + 5.35 * frac
        val = ymax * frac
        lines.append(f'  line((1.25, {y:.3f}), (13.05, {y:.3f}), stroke: black10 + 0.45pt)')
        lines.append(f'  content((1.12, {y:.3f}), [{val:.1f}], anchor: "east")')
    gap = 11.8 / len(metrics)
    for idx, (name, measured, target, error_pct) in enumerate(metrics):
        cx = 1.25 + gap * (idx + 0.5)
        group_left = cx - gap * 0.22
        bar_w = gap * 0.18
        y_measured = 0.95 + 5.35 * (measured / ymax)
        y_target = 0.95 + 5.35 * (target / ymax)
        lines.append(f'  rect(({group_left - bar_w / 2:.3f}, 0.95), ({group_left + bar_w / 2:.3f}, {y_measured:.3f}), fill: garnet, stroke: none)')
        lines.append(f'  rect(({group_left + bar_w * 1.2:.3f}, 0.95), ({group_left + bar_w * 2.2:.3f}, {y_target:.3f}), fill: atlantic, stroke: none)')
        lines.append(f'  content(({cx:.3f}, 0.63), [{typst_escape(name)}], anchor: "north")')
        lines.append(f'  content(({cx:.3f}, {max(y_measured, y_target) + 0.28:.3f}), [#text(size: 7pt, fill: horseshoe)[{error_pct:.2f}\\%]], anchor: "center")')
    lines.extend(
        [
            '  line((8.6, 6.55), (9.05, 6.55), stroke: garnet + 1.1pt)',
            '  content((9.18, 6.55), [measured], anchor: "west")',
            '  line((10.55, 6.55), (11.0, 6.55), stroke: atlantic + 1.1pt)',
            '  content((11.13, 6.55), [target], anchor: "west")',
            '  content((1.25, 6.55), [green labels = relative error percent], anchor: "west")',
            '  content((0.3, 3.6), [modulus (GPa)], angle: 90deg)',
            '})',
            '#v(0.4mm)',
            '#text(size: 7pt, fill: rgb("#5C5C5C"))[PASS gate: max relative error = 9.88% against the 10% Kok 2022 tolerance.]',
            '',
        ]
    )
    typ.write_text("\n".join(lines), encoding="utf-8")
    subprocess.run(["typst", "compile", typ.name, "stage11_effective_moduli.png"], cwd=str(FIGURES_DIR), check=True)
    subprocess.run(["typst", "compile", typ.name, "stage11_effective_moduli.pdf"], cwd=str(FIGURES_DIR), check=True)
    return FIGURES_DIR / "stage11_effective_moduli.png"


def tensor_stack(grid: pv.DataSet, prefix: str) -> np.ndarray | None:
    tensors: list[np.ndarray] = []
    for name in grid.cell_data.keys():
        if not name.startswith(prefix):
            continue
        arr = np.asarray(grid.cell_data[name], dtype=float)
        if arr.ndim == 2 and arr.shape[1] == 9:
            tensors.append(arr.reshape((-1, 3, 3)))
    if not tensors:
        return None
    return np.mean(np.stack(tensors, axis=0), axis=0)


def reset_to_reference(grid: pv.DataSet) -> None:
    disp = np.asarray(grid.point_data.get("Displacement", np.zeros((grid.n_points, 3))), dtype=float)
    if disp.ndim == 2 and disp.shape[1] == 3:
        grid.points = np.asarray(grid.points, dtype=float) - disp


def build_element_meta() -> dict[int, dict[str, object]]:
    mesh = meshio.read(GEOM_DIR / "panel.msh")
    physicals = np.asarray(mesh.cell_data["gmsh:physical"][0], dtype=int)
    tetra10 = mesh.cells[0].data
    _ = tetra10  # keep correspondence explicit
    field_map = {int(meta[0]): name for name, meta in mesh.field_data.items() if int(meta[1]) == 3}
    sidecar = json.loads((GEOM_DIR / "orientations.json").read_text(encoding="utf-8"))
    groups = {str(group["name"]): group for group in sidecar["groups"]}

    element_meta: dict[int, dict[str, object]] = {}
    element_id = 1
    for physical_id in physicals:
        group = groups[field_map[int(physical_id)]]
        if bool(group.get("isotropic", False)):
            continue
        element_meta[element_id] = group
        element_id += 1
    for physical_id in physicals:
        group = groups[field_map[int(physical_id)]]
        if not bool(group.get("isotropic", False)):
            continue
        element_meta[element_id] = group
        element_id += 1
    return element_meta


def enrich_grid(grid: pv.DataSet, element_meta: dict[int, dict[str, object]], warp_factor: float) -> pv.DataSet:
    mesh = grid.copy(deep=True)
    reset_to_reference(mesh)
    disp = np.asarray(mesh.point_data.get("Displacement", np.zeros((mesh.n_points, 3))), dtype=float)
    mesh.point_data["displacement_magnitude_mm"] = np.linalg.norm(disp, axis=1) * 1000.0

    stress = tensor_stack(mesh, "3DELEM_Strs")
    if stress is None:
        raise RuntimeError("Stage 11 VTK has no stress tensors")
    mean_trace = np.trace(stress, axis1=1, axis2=2) / 3.0
    deviator = stress.copy()
    for axis in range(3):
        deviator[:, axis, axis] -= mean_trace
    mesh.cell_data["von_mises_mpa"] = np.sqrt(1.5 * np.sum(deviator * deviator, axis=(1, 2))) / 1.0e6

    sigma1 = np.full(mesh.n_cells, np.nan, dtype=float)
    element_ids = np.asarray(mesh.cell_data["ELEMENT_ID"], dtype=int)
    for idx, element_id in enumerate(element_ids):
        meta = element_meta.get(int(element_id))
        if meta is None or bool(meta.get("isotropic", False)) or float(meta.get("nominal_angle_deg", 999.0)) != 0.0:
            continue
        v1 = np.asarray(meta["fiber_direction_unit_vector"], dtype=float)
        sigma1[idx] = float(v1 @ stress[idx] @ v1) / 1.0e6
    mesh.cell_data["sigma1_0deg_mpa"] = sigma1

    if warp_factor > 0.0:
        warped = mesh.warp_by_vector("Displacement", factor=warp_factor)
        for key in mesh.point_data:
            warped.point_data[key] = mesh.point_data[key]
        for key in mesh.cell_data:
            warped.cell_data[key] = mesh.cell_data[key]
        return warped
    return mesh


def field_time(grid: pv.DataSet) -> float:
    try:
        return float(np.ravel(grid.field_data["TIME"])[0])
    except Exception:
        return 0.0


def reference_bounds(paths: list[pathlib.Path]) -> tuple[float, float, float, float, float, float]:
    bounds = []
    for path in paths:
        grid = pv.read(path)
        reset_to_reference(grid)
        bounds.append(tuple(float(value) for value in grid.bounds))
    return (
        min(bound[0] for bound in bounds),
        max(bound[1] for bound in bounds),
        min(bound[2] for bound in bounds),
        max(bound[3] for bound in bounds),
        min(bound[4] for bound in bounds),
        max(bound[5] for bound in bounds),
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
    if max_disp <= 0.0:
        return 0.0
    diag = characteristic_length(bounds) or 1.0
    if max_disp > 0.05 * diag:
        return 1.0
    return float(min(max(0.08 * diag / max_disp, 1.0), 60.0))


def camera_for(bounds: tuple[float, float, float, float, float, float]):
    center = np.array([(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2], dtype=float)
    radius = characteristic_length(bounds) or 1.0
    direction = np.array([1.65, -1.45, 1.1], dtype=float)
    direction /= np.linalg.norm(direction)
    position = center + direction * radius * 1.9
    return tuple(position.tolist()), tuple(center.tolist()), (0.0, 0.0, 1.0)


def scalar_range(paths: list[pathlib.Path], scalar: str, location: str, warp: float, element_meta: dict[int, dict[str, object]]) -> tuple[float, float]:
    values: list[np.ndarray] = []
    for path in paths:
        grid = enrich_grid(pv.read(path), element_meta, warp)
        if location == "point":
            arr = np.asarray(grid.point_data[scalar], dtype=float).ravel()
        else:
            arr = np.asarray(grid.cell_data[scalar], dtype=float).ravel()
        arr = arr[np.isfinite(arr)]
        if arr.size:
            values.append(arr)
    if not values:
        return (0.0, 1.0)
    data = np.concatenate(values)
    low = float(np.nanpercentile(data, 0.5))
    high = float(np.nanpercentile(data, 99.5))
    if scalar != "sigma1_0deg_mpa":
        low = min(0.0, low)
    if not math.isfinite(high) or high <= low:
        high = float(np.nanmax(data))
    if not math.isfinite(low):
        low = 0.0
    if not math.isfinite(high) or high <= low:
        high = low + 1.0
    return low, high


def setup_plotter(plotter: pv.Plotter, camera, bounds) -> None:
    plotter.set_background(WHITE)
    plotter.camera_position = camera
    plotter.camera.parallel_projection = True
    plotter.camera.parallel_scale = 0.6 * max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4], 1.0e-9)
    plotter.add_axes(color=CHARCOAL, line_width=2)


def scalar_bar(title: str) -> dict[str, object]:
    return {
        "title": title,
        "title_font_size": 15,
        "label_font_size": 13,
        "color": CHARCOAL,
        "fmt": "%.3g",
        "shadow": False,
        "n_labels": 5,
    }


def render_field_panel(
    vtk_path: pathlib.Path,
    element_meta: dict[int, dict[str, object]],
    scalar: str,
    location: str,
    title: str,
    units: str,
    out_png: pathlib.Path,
    warp: float,
    camera,
    bounds,
    clim: tuple[float, float],
    zero_deg_only: bool = False,
    size: tuple[int, int] = PANEL_SIZE,
) -> None:
    mesh = enrich_grid(pv.read(vtk_path), element_meta, warp)
    plotter = pv.Plotter(off_screen=True, window_size=size)
    setup_plotter(plotter, camera, bounds)
    if zero_deg_only:
        plotter.add_mesh(mesh, color=BLACK10, opacity=0.18, show_edges=False)
        cell_ids = np.where(np.isfinite(np.asarray(mesh.cell_data[scalar], dtype=float)))[0]
        subset = mesh.extract_cells(cell_ids.tolist())
        plotter.add_mesh(
            subset,
            scalars=scalar,
            preference=location,
            cmap=FIELD_RAMP,
            clim=clim,
            show_edges=False,
            scalar_bar_args=scalar_bar(f"{title} ({units})"),
        )
    else:
        plotter.add_mesh(
            mesh,
            scalars=scalar,
            preference=location,
            cmap=FIELD_RAMP,
            clim=clim,
            show_edges=False,
            scalar_bar_args=scalar_bar(f"{title} ({units})"),
        )
    plotter.screenshot(str(out_png))
    plotter.close()


def fit_image(src: Image.Image, size: tuple[int, int]) -> Image.Image:
    image = src.convert("RGB")
    image.thumbnail(size, Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", size, WHITE)
    canvas.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return canvas


def compose_composite(panels: list[pathlib.Path], chart_png: pathlib.Path) -> pathlib.Path:
    canvas = Image.new("RGB", COMPOSITE_SIZE, WHITE)
    draw = ImageDraw.Draw(canvas)
    draw.rectangle((0, 0, COMPOSITE_SIZE[0] - 1, COMPOSITE_SIZE[1] - 1), outline=hex_to_rgb(CHARCOAL), width=3)
    draw.rectangle((0, 0, COMPOSITE_SIZE[0], 118), fill=hex_to_rgb(BLACK10))
    draw.text((58, 24), TITLE, font=FONT_TITLE, fill=hex_to_rgb(CHARCOAL))
    draw.text((58, 74), "axial-x deformed field composite with Kok 2022 modulus gate", font=FONT_SUBTITLE, fill=hex_to_rgb(BLACK70))
    draw.rectangle((2060, 30, 2318, 86), fill=hex_to_rgb(HORSESHOE))
    draw.text((2124, 45), "PASS", font=FONT_PANEL, fill=hex_to_rgb(WHITE), anchor="mm")

    positions = [(60, 145), (1230, 145), (60, 855), (1230, 855)]
    labels = [
        "Displacement magnitude (mm)",
        "von Mises stress (MPa)",
        "0 deg tow local sigma1 (MPa)",
        "Measured vs target effective moduli",
    ]
    images = [Image.open(path) for path in panels] + [Image.open(chart_png)]
    for (x, y), label, image in zip(positions, labels, images):
        draw.text((x, y - 36), label, font=FONT_PANEL, fill=hex_to_rgb(CHARCOAL))
        draw.rectangle((x - 1, y - 1, x + PANEL_SIZE[0], y + PANEL_SIZE[1]), outline=hex_to_rgb(BLACK30), width=2)
        canvas.paste(fit_image(image, PANEL_SIZE), (x, y))

    footer_y = 1512
    draw.rectangle((0, footer_y, COMPOSITE_SIZE[0], COMPOSITE_SIZE[1]), fill=hex_to_rgb(BLACK10))
    draw.rectangle((0, footer_y, COMPOSITE_SIZE[0] - 1, COMPOSITE_SIZE[1] - 1), outline=hex_to_rgb(CHARCOAL), width=2)
    draw.text((58, footer_y + 23), "Measured moduli:", font=FONT_PANEL, fill=hex_to_rgb(CHARCOAL))
    moduli = measured_moduli()
    summary = f"Ex = {moduli['E_x_GPa']:.2f} GPa   Ey = {moduli['E_y_GPa']:.2f} GPa   Gxy = {moduli['G_xy_GPa']:.2f} GPa"
    draw.text((305, footer_y + 24), summary, font=FONT_BODY, fill=hex_to_rgb(CHARCOAL))
    out_png = FIGURES_DIR / "stage_11_composite.png"
    canvas.save(out_png)
    return out_png


def title_card(body: str) -> Image.Image:
    canvas = Image.new("RGB", VIDEO_SIZE, WHITE)
    draw = ImageDraw.Draw(canvas)
    draw.rectangle((0, 0, VIDEO_SIZE[0] - 1, VIDEO_SIZE[1] - 1), outline=hex_to_rgb(CHARCOAL), width=4)
    draw.text((96, 148), TITLE, font=FONT_CARD, fill=hex_to_rgb(CHARCOAL))
    draw.rectangle((96, 246, 344, 312), fill=hex_to_rgb(HORSESHOE))
    draw.text((220, 279), "PASS", font=font(32, True), fill=hex_to_rgb(WHITE), anchor="mm")
    draw_wrapped(draw, (96, 374), body, FONT_CARD_BODY, hex_to_rgb(CHARCOAL), 1620, 10)
    return canvas


def compose_motion_frame(case_pngs: list[pathlib.Path], t: float, progress: float) -> Image.Image:
    canvas = Image.new("RGB", VIDEO_SIZE, WHITE)
    draw = ImageDraw.Draw(canvas)
    draw.rectangle((0, 0, VIDEO_SIZE[0], 80), fill=hex_to_rgb(WHITE))
    draw.text((42, 24), TITLE, font=FONT_SUBTITLE, fill=hex_to_rgb(CHARCOAL))
    draw.text((1535, 24), f"t = {t:.4g} s", font=FONT_SUBTITLE, fill=hex_to_rgb(BLACK70))

    panel_width = 620
    panel_height = 780
    x_positions = [20, 650, 1280]
    labels = ["axial-x", "transverse-y", "shear-xy"]
    for x, label, path in zip(x_positions, labels, case_pngs):
        draw.text((x, 92), label, font=FONT_PANEL, fill=hex_to_rgb(CHARCOAL))
        draw.rectangle((x - 1, 128 - 1, x + panel_width, 128 + panel_height), outline=hex_to_rgb(BLACK30), width=2)
        panel = fit_image(Image.open(path), (panel_width, panel_height))
        canvas.paste(panel, (x, 128))

    y0 = 984
    draw.rectangle((0, 940, VIDEO_SIZE[0], VIDEO_SIZE[1]), fill=hex_to_rgb(BLACK10))
    draw.text((54, 956), "Displacement magnitude shared across all three load cases; moving cursor marks synchronized progress.", font=FONT_SMALL, fill=hex_to_rgb(CHARCOAL))
    x1, x2 = 90, 1830
    draw.line((x1, y0 + 40, x2, y0 + 40), fill=hex_to_rgb(BLACK50), width=3)
    cx = int(x1 + max(0.0, min(1.0, progress)) * (x2 - x1))
    draw.rectangle((x1, y0 + 33, cx, y0 + 47), fill=hex_to_rgb(HORSESHOE))
    draw.line((cx, y0 + 22, cx, y0 + 58), fill=hex_to_rgb(GARNET), width=4)
    draw.text((x1, y0 + 64), "start", font=FONT_SMALL, fill=hex_to_rgb(BLACK70))
    draw.text((x2 - 36, y0 + 64), "end", font=FONT_SMALL, fill=hex_to_rgb(BLACK70))
    return canvas


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


def render_case_scenes(case: str, element_meta: dict[int, dict[str, object]], warp: float, camera, bounds, clim: tuple[float, float], out_dir: pathlib.Path) -> list[pathlib.Path]:
    paths = vtk_series(case)
    rendered: list[pathlib.Path] = []
    for idx, path in enumerate(paths):
        mesh = enrich_grid(pv.read(path), element_meta, warp)
        plotter = pv.Plotter(off_screen=True, window_size=(620, 780))
        setup_plotter(plotter, camera, bounds)
        plotter.add_mesh(
            mesh,
            scalars="displacement_magnitude_mm",
            preference="point",
            cmap=FIELD_RAMP,
            clim=clim,
            show_edges=False,
            scalar_bar_args=None,
        )
        out_png = out_dir / f"{case}_{idx:03d}.png"
        plotter.screenshot(str(out_png))
        plotter.close()
        rendered.append(out_png)
    return rendered


def render_motion(element_meta: dict[int, dict[str, object]], warp: float, camera, bounds, clim: tuple[float, float]) -> pathlib.Path:
    out_mp4 = FIGURES_DIR / "stage_11_motion_hd.mp4"
    case_names = ["stage11_axial_x", "stage11_transverse_y", "stage11_shear_xy"]
    with tempfile.TemporaryDirectory(prefix="stage11_motion_", dir=str(FIGURES_DIR)) as tmp:
        tmpdir = pathlib.Path(tmp)
        scenes = {
            case: render_case_scenes(case, element_meta, warp, camera, bounds, clim, tmpdir)
            for case in case_names
        }
        frame_index = 0
        start = title_card("Three synchronized pubviz reruns: axial-x, transverse-y, and shear-xy displacement evolution.")
        moduli = measured_moduli()
        end = title_card(
            f"Final measured moduli: Ex = {moduli['E_x_GPa']:.2f} GPa, "
            f"Ey = {moduli['E_y_GPa']:.2f} GPa, "
            f"Gxy = {moduli['G_xy_GPa']:.2f} GPa."
        )
        for _ in range(VIDEO_CARD_SECONDS * VIDEO_FPS):
            start.save(tmpdir / f"frame_{frame_index:04d}.png")
            frame_index += 1

        motion_frames = VIDEO_MOTION_SECONDS * VIDEO_FPS
        times = [field_time(pv.read(path)) for path in vtk_series("stage11_axial_x")]
        for i in range(motion_frames):
            scene_idx = round(i * (len(times) - 1) / max(motion_frames - 1, 1))
            case_pngs = [scenes[case][scene_idx] for case in case_names]
            compose_motion_frame(case_pngs, times[scene_idx], i / max(motion_frames - 1, 1)).save(tmpdir / f"frame_{frame_index:04d}.png")
            frame_index += 1

        for _ in range(VIDEO_CARD_SECONDS * VIDEO_FPS):
            end.save(tmpdir / f"frame_{frame_index:04d}.png")
            frame_index += 1

        encode_video(tmpdir, out_mp4)
    return out_mp4


def main() -> int:
    pv.OFF_SCREEN = True
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    element_meta = build_element_meta()
    case_paths = vtk_series("stage11_axial_x") + vtk_series("stage11_transverse_y") + vtk_series("stage11_shear_xy")
    bounds = reference_bounds(case_paths)
    warp = automatic_warp(case_paths, bounds)
    camera = camera_for(bounds)

    chart_png = typst_chart()
    axial_paths = vtk_series("stage11_axial_x")
    disp_clim = scalar_range(axial_paths, "displacement_magnitude_mm", "point", warp, element_meta)
    vm_clim = scalar_range(axial_paths, "von_mises_mpa", "cell", warp, element_meta)
    sigma1_clim = scalar_range(axial_paths, "sigma1_0deg_mpa", "cell", warp, element_meta)

    with tempfile.TemporaryDirectory(prefix="stage11_composite_", dir=str(FIGURES_DIR)) as tmp:
        tmpdir = pathlib.Path(tmp)
        final_axial = axial_paths[-1]
        panels = [
            tmpdir / "disp.png",
            tmpdir / "vm.png",
            tmpdir / "sigma1.png",
        ]
        render_field_panel(final_axial, element_meta, "displacement_magnitude_mm", "point", "Displacement magnitude", "mm", panels[0], warp, camera, bounds, disp_clim)
        render_field_panel(final_axial, element_meta, "von_mises_mpa", "cell", "von Mises stress", "MPa", panels[1], warp, camera, bounds, vm_clim)
        render_field_panel(final_axial, element_meta, "sigma1_0deg_mpa", "cell", "0 deg tow local sigma1", "MPa", panels[2], warp, camera, bounds, sigma1_clim, zero_deg_only=True)
        composite = compose_composite(panels, chart_png)
    motion = render_motion(element_meta, warp, camera, bounds, scalar_range(case_paths, "displacement_magnitude_mm", "point", warp, element_meta))
    print(composite)
    print(motion)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
