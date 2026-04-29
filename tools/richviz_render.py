"""Render rich OpenRadioss animation VTK series to publication figures.

The OpenRadioss ``anim_to_vtk`` converter writes deformed point coordinates
and the displacement vector.  The plotting path reconstructs reference
coordinates from ``points - Displacement`` before applying a controlled visual
warp by ``Displacement``.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
import pathlib
import re
import shutil
import subprocess
import tempfile
from typing import Iterable

import numpy as np
import pyvista as pv


ROOT = pathlib.Path(__file__).resolve().parents[1]
BRAND_RAMP = ["#FFF2E3", "#A49137", "#73000A", "#1F414D"]
GARNET = "#73000A"
ATLANTIC = "#466A9F"
HORSESHOE = "#65780B"
CONGAREE = "#1F414D"
CHARCOAL = "#363636"
BLACK10 = "#ECECEC"
WHITE = "#FFFFFF"
FPS = 24
VIDEO_SECONDS = 5
VIDEO_FRAMES = FPS * VIDEO_SECONDS


@dataclasses.dataclass(frozen=True)
class StageRenderConfig:
    stage: int
    stage_dir: pathlib.Path
    run_name: str
    title: str
    probe_xyz: tuple[float, float, float]
    probe_component: int
    probe_label: str
    load_label: str
    final_force_n: float
    ramp_end_s: float
    xsection_origin: tuple[float, float, float]
    xsection_normal: tuple[float, float, float]
    warp_factor: float
    displacement_reference_m: float
    error_fraction: float
    wall_clock_s: float
    verdict: str
    result_context: str
    camera_position: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]

    @property
    def figures_dir(self) -> pathlib.Path:
        return self.stage_dir / "figures"

    @property
    def run_dir(self) -> pathlib.Path:
        return self.stage_dir / "runs" / self.run_name

    @property
    def prefix(self) -> str:
        return f"stage_{self.stage:02d}"

    @property
    def out_path(self) -> pathlib.Path:
        return self.run_dir / f"{self.run_name}_0001.out"


def parse_stage_config(stage: int) -> StageRenderConfig:
    if stage == 1:
        stage_dir = ROOT / "tests" / "stage_01_beam_bending"
        results = json.loads((stage_dir / "results" / "results.json").read_text(encoding="utf-8"))
        metric = results["metrics"]["4pt_M1"]
        return StageRenderConfig(
            stage=1,
            stage_dir=stage_dir,
            run_name="stage01_4pt_M1_richviz",
            title="Stage 01 Beam Bending, 4-Point M1",
            probe_xyz=(0.100, 0.0, -0.005),
            probe_component=2,
            probe_label="Midspan displacement",
            load_label="Loaded-line resultant",
            final_force_n=1000.0,
            ramp_end_s=0.0045,
            xsection_origin=(0.100, 0.0, 0.0),
            xsection_normal=(1.0, 0.0, 0.0),
            warp_factor=8.0,
            displacement_reference_m=float(metric["delta_eb_m"]),
            error_fraction=float(metric["rel_error_eb"]),
            wall_clock_s=float(metric["wall_clock_s"]),
            verdict=str(results["verdict"]),
            result_context="4pt_M1 gated Euler-Bernoulli check",
            camera_position=((0.255, -0.118, 0.066), (0.100, 0.0, 0.0), (0.0, 0.0, 1.0)),
        )
    if stage == 2:
        stage_dir = ROOT / "tests" / "stage_02_cantilever_large"
        results = json.loads((stage_dir / "results" / "results.json").read_text(encoding="utf-8"))
        metric = results["metrics"]["baseline_alpha_1"]
        err = max(float(metric["err_x"]), float(metric["err_y"]))
        return StageRenderConfig(
            stage=2,
            stage_dir=stage_dir,
            run_name="stage02_implicit_baseline_a1p00_richviz",
            title="Stage 02 Cantilever, alpha=1 Baseline",
            probe_xyz=(1.000, 0.0, 0.0),
            probe_component=1,
            probe_label="Tip sag",
            load_label="Tip-face resultant",
            final_force_n=11.25,
            ramp_end_s=1.0,
            xsection_origin=(0.500, 0.0, 0.0),
            xsection_normal=(1.0, 0.0, 0.0),
            warp_factor=1.0,
            displacement_reference_m=-float(metric["dy_ref_over_L"]),
            error_fraction=err,
            wall_clock_s=float(metric["wall_clock_s"]),
            verdict="alpha=1 PASS",
            result_context="baseline alpha=1 gated elastica check",
            camera_position=((1.38, -0.68, 0.24), (0.50, -0.11, 0.0), (0.0, 0.0, 1.0)),
        )
    raise SystemExit(f"unsupported stage {stage}")


def vtk_paths(config: StageRenderConfig) -> list[pathlib.Path]:
    paths = sorted(config.run_dir.glob(f"{config.run_name}A[0-9][0-9][0-9].vtk"))
    if not paths:
        raise FileNotFoundError(f"no rich VTK frames found in {config.run_dir}")
    return paths


def clean_grid(grid: pv.DataSet, warp_factor: float) -> pv.DataSet:
    grid = grid.copy(deep=True)
    disp = np.asarray(grid.point_data["Displacement"], dtype=float)
    grid.points = np.asarray(grid.points, dtype=float) - disp
    add_fields(grid)
    warped = grid.warp_by_vector("Displacement", factor=warp_factor)
    for key in grid.point_data:
        warped.point_data[key] = grid.point_data[key]
    for key in grid.cell_data:
        warped.cell_data[key] = grid.cell_data[key]
    return warped


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


def add_fields(grid: pv.DataSet) -> None:
    disp = np.asarray(grid.point_data["Displacement"], dtype=float)
    grid.point_data["displacement_magnitude_mm"] = np.linalg.norm(disp, axis=1) * 1000.0
    stress = tensor_stack(grid, "3DELEM_Strs")
    if stress is not None:
        grid.cell_data["von_mises_mpa"] = von_mises_from_tensor(stress) / 1.0e6
    strain = tensor_stack(grid, "3DELEM_Stra")
    if strain is not None:
        grid.cell_data["eq_strain"] = eq_strain_from_tensor(strain)
    plastic_names = [name for name in grid.cell_data.keys() if "Plastic" in name or "EPSP" in name]
    for name in plastic_names:
        arr = np.asarray(grid.cell_data[name], dtype=float)
        if arr.ndim == 1:
            grid.cell_data["plastic_strain"] = arr
            break
        if arr.ndim == 2:
            grid.cell_data["plastic_strain"] = np.linalg.norm(arr, axis=1)
            break


def robust_range(paths: list[pathlib.Path], warp_factor: float) -> dict[str, tuple[float, float]]:
    samples = {"displacement_magnitude_mm": [], "von_mises_mpa": [], "eq_strain": []}
    for path in paths:
        grid = clean_grid(pv.read(path), warp_factor)
        samples["displacement_magnitude_mm"].append(np.asarray(grid.point_data["displacement_magnitude_mm"]))
        if "von_mises_mpa" in grid.cell_data:
            samples["von_mises_mpa"].append(np.asarray(grid.cell_data["von_mises_mpa"]))
        if "eq_strain" in grid.cell_data:
            samples["eq_strain"].append(np.asarray(grid.cell_data["eq_strain"]))
    ranges: dict[str, tuple[float, float]] = {}
    for name, chunks in samples.items():
        if not chunks:
            continue
        values = np.concatenate(chunks)
        high = float(np.nanpercentile(values, 99.5))
        if not math.isfinite(high) or high <= 0.0:
            high = float(np.nanmax(values))
        ranges[name] = (0.0, high if high > 0.0 else 1.0)
    return ranges


def field_time(grid: pv.DataSet) -> float:
    try:
        return float(np.ravel(grid.field_data["TIME"])[0])
    except Exception:
        return 0.0


def add_common_scene(
    plotter: pv.Plotter,
    mesh: pv.DataSet,
    title: str,
    camera_position: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    bounds: bool = True,
) -> None:
    plotter.set_background("white")
    plotter.camera_position = camera_position
    plotter.add_text(title, position="upper_left", font_size=12, color=CHARCOAL, font="courier")
    if bounds:
        plotter.show_bounds(
            grid="front",
            location="outer",
            color=CHARCOAL,
            font_size=8,
            xtitle="x (m)",
            ytitle="y (m)",
            ztitle="z (m)",
        )
    plotter.add_axes(color=CHARCOAL, line_width=2)
    plotter.camera.zoom(1.08)


def scalar_bar(title: str) -> dict[str, object]:
    return {
        "title": title,
        "title_font_size": 12,
        "label_font_size": 10,
        "color": CHARCOAL,
        "vertical": True,
        "fmt": "%.3g",
        "shadow": False,
        "n_labels": 5,
    }


def render_still(
    config: StageRenderConfig,
    path: pathlib.Path,
    scalar: str,
    scalar_location: str,
    clim: tuple[float, float],
    title: str,
    colorbar_title: str,
    out_png: pathlib.Path,
) -> None:
    mesh = clean_grid(pv.read(path), config.warp_factor)
    plotter = pv.Plotter(off_screen=True, window_size=(1800, 1350))
    values = mesh.point_data[scalar] if scalar_location == "point" else mesh.cell_data[scalar]
    plotter.add_mesh(
        mesh,
        scalars=values,
        preference="point" if scalar_location == "point" else "cell",
        cmap=BRAND_RAMP,
        clim=clim,
        show_edges=False,
        scalar_bar_args=scalar_bar(colorbar_title),
    )
    add_common_scene(plotter, mesh, title, config.camera_position)
    annotation = (
        f"{config.result_context}\n"
        f"verdict: {config.verdict} | rel error: {100.0 * config.error_fraction:.3f}% | "
        f"wall clock: {config.wall_clock_s:.1f} s\n"
        f"warp: {config.warp_factor:g}x"
    )
    plotter.add_text(annotation, position="lower_right", font_size=9, color=CHARCOAL, font="courier")
    plotter.screenshot(str(out_png))
    plotter.close()


def render_xsection(
    config: StageRenderConfig,
    path: pathlib.Path,
    clim: tuple[float, float],
    out_png: pathlib.Path,
) -> None:
    mesh = clean_grid(pv.read(path), config.warp_factor)
    section = mesh.slice(normal=config.xsection_normal, origin=config.xsection_origin)
    plotter = pv.Plotter(off_screen=True, window_size=(1800, 1350))
    plotter.add_mesh(
        section,
        scalars="von_mises_mpa",
        preference="cell",
        cmap=BRAND_RAMP,
        clim=clim,
        show_edges=True,
        edge_color=CHARCOAL,
        line_width=0.4,
        scalar_bar_args=scalar_bar("von Mises stress (MPa)"),
    )
    plotter.set_background("white")
    plotter.view_yz()
    plotter.camera.zoom(2.8 if config.stage == 1 else 4.5)
    plotter.show_bounds(
        grid="front",
        location="outer",
        color=CHARCOAL,
        font_size=8,
        xtitle="x (m)",
        ytitle="y (m)",
        ztitle="z (m)",
    )
    plotter.add_text(
        f"{config.title}\nmid-span cross-section, warp {config.warp_factor:g}x",
        position="upper_left",
        font_size=12,
        color=WHITE if config.stage == 1 else CHARCOAL,
        font="courier",
    )
    plotter.screenshot(str(out_png))
    plotter.close()


def render_frame(
    config: StageRenderConfig,
    path: pathlib.Path,
    scalar: str,
    scalar_location: str,
    clim: tuple[float, float],
    out_png: pathlib.Path,
    title: str,
    colorbar_title: str,
) -> None:
    mesh = clean_grid(pv.read(path), config.warp_factor)
    plotter = pv.Plotter(off_screen=True, window_size=(1920, 1080))
    values = mesh.point_data[scalar] if scalar_location == "point" else mesh.cell_data[scalar]
    plotter.add_mesh(
        mesh,
        scalars=values,
        preference="point" if scalar_location == "point" else "cell",
        cmap=BRAND_RAMP,
        clim=clim,
        show_edges=False,
        scalar_bar_args=scalar_bar(colorbar_title),
    )
    add_common_scene(plotter, mesh, title, config.camera_position)
    plotter.add_text(f"t = {field_time(mesh):.4g} s", position="lower_right", font_size=11, color=CHARCOAL, font="courier")
    plotter.screenshot(str(out_png))
    plotter.close()


def render_double_frame(
    config: StageRenderConfig,
    path: pathlib.Path,
    ranges: dict[str, tuple[float, float]],
    out_png: pathlib.Path,
) -> None:
    mesh = clean_grid(pv.read(path), config.warp_factor)
    plotter = pv.Plotter(off_screen=True, window_size=(1920, 1080), shape=(2, 1), border=False)
    for idx, (scalar, location, title, cbtitle) in enumerate(
        [
            ("displacement_magnitude_mm", "point", "Displacement magnitude", "displacement (mm)"),
            ("von_mises_mpa", "cell", "von Mises stress", "von Mises (MPa)"),
        ]
    ):
        plotter.subplot(idx, 0)
        values = mesh.point_data[scalar] if location == "point" else mesh.cell_data[scalar]
        plotter.add_mesh(
            mesh,
            scalars=values,
            preference=location,
            cmap=BRAND_RAMP,
            clim=ranges[scalar],
            show_edges=False,
            scalar_bar_args=scalar_bar(cbtitle),
        )
        add_common_scene(plotter, mesh, f"{config.title} - {title}", config.camera_position, bounds=False)
        if idx == 1:
            plotter.add_text(f"t = {field_time(mesh):.4g} s", position="lower_right", font_size=11, color=CHARCOAL, font="courier")
    plotter.screenshot(str(out_png))
    plotter.close()


def selected_video_paths(paths: list[pathlib.Path]) -> list[pathlib.Path]:
    return [paths[i % len(paths)] for i in range(VIDEO_FRAMES)]


def ffmpeg_encode(frame_dir: pathlib.Path, out_mp4: pathlib.Path) -> None:
    cmd = [
        "ffmpeg",
        "-y",
        "-framerate",
        str(FPS),
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
    ]
    subprocess.run(cmd, check=True, cwd=str(ROOT))


def render_video(
    config: StageRenderConfig,
    paths: list[pathlib.Path],
    ranges: dict[str, tuple[float, float]],
    mode: str,
    out_mp4: pathlib.Path,
) -> None:
    frames = selected_video_paths(paths)
    with tempfile.TemporaryDirectory(prefix=f"{config.prefix}_{mode}_", dir=str(config.figures_dir)) as tmp:
        frame_dir = pathlib.Path(tmp)
        for idx, path in enumerate(frames):
            out_png = frame_dir / f"frame_{idx:04d}.png"
            if mode == "displacement":
                render_frame(
                    config,
                    path,
                    "displacement_magnitude_mm",
                    "point",
                    ranges["displacement_magnitude_mm"],
                    out_png,
                    f"{config.title} - displacement",
                    "displacement (mm)",
                )
            elif mode == "vonmises":
                render_frame(
                    config,
                    path,
                    "von_mises_mpa",
                    "cell",
                    ranges["von_mises_mpa"],
                    out_png,
                    f"{config.title} - von Mises",
                    "von Mises (MPa)",
                )
            elif mode == "doubleview":
                render_double_frame(config, path, ranges, out_png)
            else:
                raise ValueError(mode)
        ffmpeg_encode(frame_dir, out_mp4)


def parse_energy_rows(out_path: pathlib.Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    if not out_path.exists():
        return rows
    for line in out_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = line.split()
        if len(parts) < 10 or not parts[0].isdigit():
            continue
        try:
            time = float(parts[1].replace("D", "E"))
            internal = float(parts[6].replace("D", "E"))
            kinetic_t = float(parts[7].replace("D", "E"))
            kinetic_r = float(parts[8].replace("D", "E"))
        except ValueError:
            continue
        rows.append((time, internal + kinetic_t + kinetic_r))
    return rows


def interpolate_energy(energy_rows: list[tuple[float, float]], t: float) -> float:
    if not energy_rows:
        return 0.0
    times = np.asarray([row[0] for row in energy_rows], dtype=float)
    values = np.asarray([row[1] for row in energy_rows], dtype=float)
    return float(np.interp(t, times, values))


def history_rows(config: StageRenderConfig, paths: list[pathlib.Path]) -> list[dict[str, float]]:
    first = clean_grid(pv.read(paths[0]), config.warp_factor)
    probe = first.find_closest_point(config.probe_xyz)
    energy = parse_energy_rows(config.out_path)
    rows: list[dict[str, float]] = []
    for path in paths:
        grid = clean_grid(pv.read(path), config.warp_factor)
        t = field_time(grid)
        disp = float(np.asarray(grid.point_data["Displacement"])[probe, config.probe_component])
        force = config.final_force_n * min(max(t / config.ramp_end_s, 0.0), 1.0)
        rows.append(
            {
                "time_s": t,
                "probe_displacement_m": disp,
                "force_n": force,
                "total_energy_j": interpolate_energy(energy, t),
            }
        )
    return rows


def write_history_csv(config: StageRenderConfig, rows: list[dict[str, float]]) -> pathlib.Path:
    out = config.figures_dir / f"{config.prefix}_history.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["time_s", "probe_displacement_m", "force_n", "total_energy_j"])
        writer.writeheader()
        writer.writerows(rows)
    return out


def sampled(items: list[dict[str, float]], count: int = 90) -> list[dict[str, float]]:
    if len(items) <= count:
        return items
    idx = np.linspace(0, len(items) - 1, count).round().astype(int)
    return [items[int(i)] for i in idx]


def chart_points(
    rows: list[dict[str, float]],
    x_key: str,
    y_key: str,
    x0: float,
    y0: float,
    width: float,
    height: float,
    y_scale: float = 1.0,
) -> tuple[str, float, float, float]:
    data = sampled(rows)
    xmax = max(row[x_key] for row in rows) or 1.0
    ys = [row[y_key] * y_scale for row in rows]
    ymin = min(0.0, min(ys))
    ymax = max(ys)
    if math.isclose(ymin, ymax):
        ymax = ymin + 1.0
    points = []
    for row in data:
        x = x0 + (row[x_key] / xmax) * width
        yv = row[y_key] * y_scale
        y = y0 + ((yv - ymin) / (ymax - ymin)) * height
        points.append(f"({x:.3f}, {y:.3f})")
    return ", ".join(points), ymin, ymax, xmax


def write_history_typ(config: StageRenderConfig, rows: list[dict[str, float]], csv_path: pathlib.Path) -> pathlib.Path:
    out = config.figures_dir / f"{config.prefix}_history.typ"
    x0, width, height = 1.55, 13.0, 3.15
    panels = [
        ("probe_displacement_m", 1000.0, f"{config.probe_label} (mm)", config.displacement_reference_m * 1000.0),
        ("force_n", 1.0, f"{config.load_label} (N)", config.final_force_n),
        ("total_energy_j", 1000.0, "internal + kinetic energy (mJ)", None),
    ]
    canvas_lines: list[str] = []
    for idx, (key, scale, ylabel, reference) in enumerate(panels):
        y0 = 1.00 + (2 - idx) * 4.10
        points, ymin, ymax, xmax = chart_points(rows, "time_s", key, x0, y0, width, height, scale)
        canvas_lines.extend(
            [
                f"  rect(({x0:.3f}, {y0:.3f}), ({x0 + width:.3f}, {y0 + height:.3f}), fill: white, stroke: charcoal + 0.65pt)",
                f"  content(({x0 - 0.18:.3f}, {y0:.3f}), [{ymin:.3g}], anchor: \"east\")",
                f"  content(({x0 - 0.18:.3f}, {y0 + height:.3f}), [{ymax:.3g}], anchor: \"east\")",
                f"  content(({x0 + width / 2:.3f}, {y0 - 0.42:.3f}), [time (s)], anchor: \"north\")",
                f"  content(({x0 - 1.05:.3f}, {y0 + height / 2:.3f}), [{ylabel}], angle: 90deg)",
                f"  line({points}, stroke: garnet + 1.05pt)",
            ]
        )
        if reference is not None and not math.isclose(ymax, ymin):
            yref = y0 + ((reference - ymin) / (ymax - ymin)) * height
            canvas_lines.extend(
                [
                    f"  line(({x0:.3f}, {yref:.3f}), ({x0 + width:.3f}, {yref:.3f}), stroke: atlantic + 0.8pt)",
                    f"  content(({x0 + width + 0.12:.3f}, {yref:.3f}), [reference], anchor: \"west\")",
                ]
            )
        if key == "probe_displacement_m" and reference is not None and not math.isclose(ymax, ymin):
            gate = 0.01 if config.stage == 1 else 0.02
            for sign in (-1.0, 1.0):
                ygate_val = reference * (1.0 + sign * gate)
                ygate = y0 + ((ygate_val - ymin) / (ymax - ymin)) * height
                canvas_lines.append(f"  line(({x0:.3f}, {ygate:.3f}), ({x0 + width:.3f}, {ygate:.3f}), stroke: horseshoe + 0.55pt)")
        canvas_lines.append(f"  content(({x0 + width:.3f}, {y0 - 0.42:.3f}), [{xmax:.3g}], anchor: \"north\")")
    lines = [
        '#import "@preview/cetz:0.3.4"',
        "",
        '#set page(width: 190mm, height: 160mm, margin: 10mm)',
        '#set text(font: "Libertinus Serif", size: 8.5pt, fill: rgb("#363636"))',
        '#let garnet = rgb("#73000A")',
        '#let atlantic = rgb("#466A9F")',
        '#let horseshoe = rgb("#65780B")',
        '#let charcoal = rgb("#363636")',
        '#let white = rgb("#FFFFFF")',
        "",
        f"#align(center)[#text(size: 11pt, weight: \"bold\")[{config.title} History]]",
        "#v(1.5mm)",
        "#cetz.canvas(length: 1cm, {",
        "  import cetz.draw: *",
        *canvas_lines,
        "})",
        "",
        f"#text(size: 8pt)[Source CSV: #raw(\"{csv_path.name}\"). Garnet = OpenRadioss, atlantic = reference, horseshoe = tolerance band.]",
        "",
    ]
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def compile_typst(typ_path: pathlib.Path) -> pathlib.Path:
    pdf_path = typ_path.with_suffix(".pdf")
    subprocess.run(["typst", "compile", typ_path.name, pdf_path.name], cwd=str(typ_path.parent), check=True)
    return pdf_path


def render_stage(config: StageRenderConfig) -> None:
    pv.OFF_SCREEN = True
    config.figures_dir.mkdir(parents=True, exist_ok=True)
    paths = vtk_paths(config)
    ranges = robust_range(paths, config.warp_factor)
    final = paths[-1]

    render_video(config, paths, ranges, "displacement", config.figures_dir / f"{config.prefix}_motion_displacement.mp4")
    render_video(config, paths, ranges, "vonmises", config.figures_dir / f"{config.prefix}_motion_vonmises.mp4")
    render_video(config, paths, ranges, "doubleview", config.figures_dir / f"{config.prefix}_motion_doubleview.mp4")

    render_still(
        config,
        final,
        "displacement_magnitude_mm",
        "point",
        ranges["displacement_magnitude_mm"],
        f"{config.title} - final displacement",
        "displacement (mm)",
        config.figures_dir / f"{config.prefix}_field_displacement.png",
    )
    render_still(
        config,
        final,
        "von_mises_mpa",
        "cell",
        ranges["von_mises_mpa"],
        f"{config.title} - final von Mises",
        "von Mises (MPa)",
        config.figures_dir / f"{config.prefix}_field_vonmises.png",
    )
    render_still(
        config,
        final,
        "eq_strain",
        "cell",
        ranges["eq_strain"],
        f"{config.title} - final equivalent strain",
        "equiv. strain (-)",
        config.figures_dir / f"{config.prefix}_field_strain.png",
    )
    render_xsection(config, final, ranges["von_mises_mpa"], config.figures_dir / f"{config.prefix}_field_xsection.png")

    rows = history_rows(config, paths)
    csv_path = write_history_csv(config, rows)
    typ_path = write_history_typ(config, rows, csv_path)
    compile_typst(typ_path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=int, choices=[1, 2], required=True)
    args = parser.parse_args(argv)
    if shutil.which("ffmpeg") is None:
        raise SystemExit("ffmpeg not found")
    if shutil.which("typst") is None:
        raise SystemExit("typst not found")
    render_stage(parse_stage_config(args.stage))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
