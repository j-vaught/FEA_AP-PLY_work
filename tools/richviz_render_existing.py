"""Render existing stage VTK/animation outputs for rich visual review.

This renderer is intentionally post-processing only: it does not edit or run
simulation decks.  It consumes VTK files that already exist under a stage's
``runs`` directory and, when only compressed OpenRadioss animation frames are
present, converts those frames in a scratch directory for rendering.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import gzip
import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import tempfile
from typing import Iterable

import numpy as np
import pyvista as pv


ROOT = pathlib.Path(__file__).resolve().parents[1]
OR_ROOT = pathlib.Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
ANIM_TO_VTK = OR_ROOT / "exec" / "anim_to_vtk_linux64_gf"

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
class FrameSpec:
    source: pathlib.Path
    label: str
    run_label: str


@dataclasses.dataclass(frozen=True)
class RenderedFrame:
    vtk_path: pathlib.Path
    label: str
    run_label: str


@dataclasses.dataclass(frozen=True)
class StageConfig:
    stage: int
    stage_dir: pathlib.Path
    title: str
    frames: tuple[FrameSpec, ...]
    run_note: str
    verdict_note: str
    history_kind: str
    final_index: int = -1

    @property
    def figures_dir(self) -> pathlib.Path:
        return self.stage_dir / "figures"

    @property
    def prefix(self) -> str:
        return f"stage_{self.stage:02d}"

    @property
    def results_path(self) -> pathlib.Path:
        return self.stage_dir / "results" / "results.json"

    @property
    def timeseries_path(self) -> pathlib.Path:
        return self.stage_dir / "results" / "timeseries.csv"


@dataclasses.dataclass(frozen=True)
class FieldChoice:
    scalar: str
    location: str
    title: str
    units: str
    description: str


@dataclasses.dataclass(frozen=True)
class HistorySeries:
    label: str
    points: tuple[tuple[float, float], ...]
    color: str


def stage_dir(stage: int) -> pathlib.Path:
    matches = sorted((ROOT / "tests").glob(f"stage_{stage:02d}_*"))
    if not matches:
        raise FileNotFoundError(f"stage {stage:02d} directory not found")
    return matches[0]


def read_results(stage_path: pathlib.Path) -> dict[str, object]:
    return json.loads((stage_path / "results" / "results.json").read_text(encoding="utf-8"))


def stage03_config() -> StageConfig:
    sdir = stage_dir(3)
    frames = []
    for load, tag in [(0.50, "0p50"), (0.80, "0p80"), (1.00, "1p00"), (1.20, "1p20")]:
        run = f"stage03_medium_lf{tag}"
        frames.append(FrameSpec(sdir / "runs" / run / f"{run}A001.vtk", f"load factor {load:.2f}", run))
    return StageConfig(
        stage=3,
        stage_dir=sdir,
        title="Stage 03 Isotropic Dogbone E8",
        frames=tuple(frames),
        run_note="medium-mesh final frames from the load-factor sweep, ordered 0.50 to 1.20",
        verdict_note="Stage 03 verdict is PASS; these figures use the medium mesh that produced the reported metrics.",
        history_kind="stage03",
    )


def stage04_config() -> StageConfig:
    sdir = stage_dir(4)
    frames = []
    for ntheta in [32, 64, 128]:
        run = f"stage04_ntheta{ntheta:03d}"
        frames.append(FrameSpec(sdir / "runs" / run / f"{run}A001.vtk", f"Ntheta {ntheta}", run))
    return StageConfig(
        stage=4,
        stage_dir=sdir,
        title="Stage 04 Open-Hole Kirsch",
        frames=tuple(frames),
        run_note="mesh-refinement final frames for Ntheta 32, 64, and 128; this is a diagnostic sweep, not a dynamic time sequence",
        verdict_note="Stage 04 verdict is PASS; the primary gate is the Ntheta 64 Howland Kt comparison.",
        history_kind="stage04",
    )


def stage05_config() -> StageConfig:
    sdir = stage_dir(5)
    run = "stage05_medium"
    paths = sorted((sdir / "runs" / run).glob(f"{run}A[0-9][0-9][0-9].vtk"))
    frames = tuple(FrameSpec(path, path.stem[-4:], run) for path in paths)
    return StageConfig(
        stage=5,
        stage_dir=sdir,
        title="Stage 05 Dogbone Damage",
        frames=frames,
        run_note="medium-mesh LAW22 damage animation frames A001-A060",
        verdict_note="Stage 05 verdict is PASS; the standard rich figures use the medium mesh while the history plot compares coarse, medium, and fine force histories.",
        history_kind="stage05",
    )


def _stage06_frame_for_index(sdir: pathlib.Path, criterion: str, index: int) -> pathlib.Path | None:
    stem = f"stage06_{criterion}_type6_proxyA{index:03d}"
    vtk = sdir / "runs" / f"{stem}.vtk"
    if vtk.exists():
        return vtk
    gz = sdir / "runs" / f"{stem}.gz"
    if gz.exists():
        return gz
    return None


def stage06_config() -> StageConfig:
    sdir = stage_dir(6)
    frames: list[FrameSpec] = []
    for criterion in ["TSAIWU", "HASHIN", "PUCK"]:
        for index in range(1, 1000):
            path = _stage06_frame_for_index(sdir, criterion, index)
            if path is None:
                if index == 1:
                    continue
                break
            frames.append(FrameSpec(path, f"{criterion} A{index:03d}", f"{criterion} TYPE6 proxy"))
    return StageConfig(
        stage=6,
        stage_dir=sdir,
        title="Stage 06 Composite Failure Criteria",
        frames=tuple(frames),
        run_note="TYPE6/SOL_ORTH proxy animation frames for TSAIWU, HASHIN, and PUCK; canonical LAW25+TYPE14 is blocked",
        verdict_note="Stage 06 verdict is INCONCLUSIVE due to LAW25-on-TYPE14 starter incompatibility; the figures here show /PROP/TYPE6 proxy runs, not the canonical solid path.",
        history_kind="stage06",
    )


def generic_stage_config(stage: int) -> StageConfig | None:
    sdir = stage_dir(stage)
    run_dirs = [p for p in (sdir / "runs").glob("*") if p.is_dir()]
    best: tuple[int, pathlib.Path, list[pathlib.Path]] | None = None
    for run_dir in run_dirs:
        paths = sorted(run_dir.glob("*.vtk"))
        if not paths:
            paths = sorted(run_dir.glob("*.gz"))
        if paths and (best is None or len(paths) > best[0]):
            best = (len(paths), run_dir, paths)
    if best is None:
        loose = sorted((sdir / "runs").glob("*.vtk")) or sorted((sdir / "runs").glob("*.gz"))
        if not loose:
            return None
        frames = tuple(FrameSpec(path, path.stem, path.parent.name) for path in loose)
    else:
        _, run_dir, paths = best
        frames = tuple(FrameSpec(path, path.stem, run_dir.name) for path in paths)
    results = read_results(sdir)
    verdict = str(results.get("verdict", "UNKNOWN"))
    return StageConfig(
        stage=stage,
        stage_dir=sdir,
        title=f"Stage {stage:02d} {str(results.get('name', sdir.name)).replace('_', ' ').title()}",
        frames=frames,
        run_note="automatically selected existing run output with the most frames",
        verdict_note=f"Stage {stage:02d} verdict is {verdict}; figures show existing outputs only.",
        history_kind="generic",
    )


def config_for_stage(stage: int) -> StageConfig | None:
    if stage == 3:
        return stage03_config()
    if stage == 4:
        return stage04_config()
    if stage == 5:
        return stage05_config()
    if stage == 6:
        return stage06_config()
    return generic_stage_config(stage)


def openradioss_env() -> dict[str, str]:
    env = os.environ.copy()
    env["RAD_CFG_PATH"] = str(OR_ROOT / "hm_cfg_files")
    env["RAD_H3D_PATH"] = str(OR_ROOT / "extlib" / "h3d" / "lib" / "linux64")
    reader = str(OR_ROOT / "extlib" / "hm_reader" / "linux64")
    env["LD_LIBRARY_PATH"] = reader + ":" + env.get("LD_LIBRARY_PATH", "")
    return env


def convert_anim_to_vtk(source: pathlib.Path, scratch_dir: pathlib.Path) -> pathlib.Path:
    if source.suffix == ".gz":
        anim_name = source.name[:-3]
        anim_path = scratch_dir / anim_name
        with gzip.open(source, "rb") as src, anim_path.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        anim_path = scratch_dir / source.name
        shutil.copy2(source, anim_path)
    proc = subprocess.run(
        [str(ANIM_TO_VTK), str(anim_path)],
        cwd=str(scratch_dir),
        env=openradioss_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        err = proc.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"anim_to_vtk failed for {source}: {err}")
    vtk_path = scratch_dir / f"{anim_path.name}.vtk"
    if proc.stdout and proc.stdout.lstrip().startswith(b"# vtk"):
        vtk_path.write_bytes(proc.stdout)
    elif not vtk_path.exists():
        candidates = sorted(scratch_dir.glob("*.vtk"))
        if not candidates:
            raise FileNotFoundError(f"anim_to_vtk produced no VTK for {source}")
        candidates[0].replace(vtk_path)
    return vtk_path


def resolve_frames(config: StageConfig, scratch_dir: pathlib.Path) -> list[RenderedFrame]:
    resolved: list[RenderedFrame] = []
    for idx, frame in enumerate(config.frames):
        if frame.source.suffix == ".vtk":
            vtk_path = frame.source
        elif frame.source.suffix == ".gz" or re.search(r"A\d{3}$", frame.source.name):
            if not ANIM_TO_VTK.exists():
                raise FileNotFoundError(f"anim_to_vtk not found at {ANIM_TO_VTK}")
            stage_scratch = scratch_dir / f"{config.prefix}_{idx:03d}"
            stage_scratch.mkdir(parents=True, exist_ok=True)
            vtk_path = convert_anim_to_vtk(frame.source, stage_scratch)
        else:
            raise ValueError(f"unsupported frame input {frame.source}")
        resolved.append(RenderedFrame(vtk_path=vtk_path, label=frame.label, run_label=frame.run_label))
    return resolved


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


def pick_vector_field(grid: pv.DataSet) -> tuple[str | None, str]:
    for name in ["Displacement", "Velocity"]:
        if name in grid.point_data and np.asarray(grid.point_data[name]).ndim == 2 and np.asarray(grid.point_data[name]).shape[1] == 3:
            return name, name
    for name in grid.point_data.keys():
        arr = np.asarray(grid.point_data[name])
        if arr.ndim == 2 and arr.shape[1] == 3:
            return name, name
    return None, "none"


def add_fields(grid: pv.DataSet) -> tuple[str | None, str]:
    vector_name, vector_label = pick_vector_field(grid)
    if vector_name is not None:
        vec = np.asarray(grid.point_data[vector_name], dtype=float)
        grid.point_data["richviz_vector_magnitude_mm"] = np.linalg.norm(vec, axis=1) * 1000.0
    else:
        grid.point_data["richviz_vector_magnitude_mm"] = np.zeros(grid.n_points, dtype=float)

    stress = tensor_stack(grid, "3DELEM_Strs")
    if stress is not None:
        grid.cell_data["von_mises_mpa"] = von_mises_from_tensor(stress) / 1.0e6
    strain = tensor_stack(grid, "3DELEM_Stra")
    if strain is not None:
        grid.cell_data["eq_strain"] = eq_strain_from_tensor(strain)
    for name in grid.cell_data.keys():
        if "Plastic_Strain" in name or "EPSP" in name:
            arr = np.asarray(grid.cell_data[name], dtype=float)
            if arr.ndim == 1:
                grid.cell_data["plastic_strain"] = arr
            elif arr.ndim == 2:
                grid.cell_data["plastic_strain"] = np.linalg.norm(arr, axis=1)
            break
    for name in grid.cell_data.keys():
        if "Damage" in name:
            arr = np.asarray(grid.cell_data[name], dtype=float)
            if arr.ndim == 1:
                grid.cell_data["damage_magnitude"] = arr
            elif arr.ndim == 2:
                grid.cell_data["damage_magnitude"] = np.linalg.norm(arr, axis=1)
            break
    return vector_name, vector_label


def prepare_grid(path: pathlib.Path, warp_factor: float = 0.0) -> tuple[pv.DataSet, str | None, str]:
    grid = pv.read(path).copy(deep=True)
    vector_name, vector_label = add_fields(grid)
    if vector_name == "Displacement":
        disp = np.asarray(grid.point_data[vector_name], dtype=float)
        grid.points = np.asarray(grid.points, dtype=float) - disp
    if vector_name is not None and warp_factor:
        warped = grid.warp_by_vector(vector_name, factor=warp_factor)
        for key in grid.point_data:
            warped.point_data[key] = grid.point_data[key]
        for key in grid.cell_data:
            warped.cell_data[key] = grid.cell_data[key]
        return warped, vector_name, vector_label
    return grid, vector_name, vector_label


def mesh_bounds_union(paths: Iterable[pathlib.Path]) -> tuple[float, float, float, float, float, float]:
    bounds: list[tuple[float, float, float, float, float, float]] = []
    for path in paths:
        mesh, _, _ = prepare_grid(path, 0.0)
        bounds.append(tuple(float(v) for v in mesh.bounds))
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


def max_vector_magnitude(paths: Iterable[pathlib.Path]) -> float:
    maximum = 0.0
    for path in paths:
        grid = pv.read(path)
        vector_name, _ = pick_vector_field(grid)
        if vector_name is None:
            continue
        arr = np.asarray(grid.point_data[vector_name], dtype=float)
        if arr.size:
            maximum = max(maximum, float(np.nanmax(np.linalg.norm(arr, axis=1))))
    return maximum


def automatic_warp_factor(paths: list[pathlib.Path], bounds: tuple[float, float, float, float, float, float]) -> float:
    max_mag = max_vector_magnitude(paths)
    if max_mag <= 0.0:
        return 0.0
    target = 0.06 * characteristic_length(bounds)
    return float(min(max(target / max_mag, 1.0), 50.0))


def camera_position(bounds: tuple[float, float, float, float, float, float]) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    center = np.array([(bounds[0] + bounds[1]) / 2.0, (bounds[2] + bounds[3]) / 2.0, (bounds[4] + bounds[5]) / 2.0])
    radius = characteristic_length(bounds) or 1.0
    direction = np.array([1.7, -1.6, 1.1], dtype=float)
    direction /= np.linalg.norm(direction)
    position = center + direction * radius * 1.75
    return tuple(position.tolist()), tuple(center.tolist()), (0.0, 0.0, 1.0)


def scalar_range(paths: Iterable[pathlib.Path], scalar: str, location: str, warp_factor: float) -> tuple[float, float]:
    chunks: list[np.ndarray] = []
    for path in paths:
        mesh, _, _ = prepare_grid(path, warp_factor)
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
    low = min(0.0, float(np.nanpercentile(values, 0.5)))
    high = float(np.nanpercentile(values, 99.5))
    if not math.isfinite(high) or high <= low:
        high = float(np.nanmax(values))
    if not math.isfinite(high) or high <= low:
        high = low + 1.0
    return (low, high)


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


def add_common_scene(
    plotter: pv.Plotter,
    title: str,
    camera: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    bounds: bool = True,
) -> None:
    plotter.set_background(WHITE)
    plotter.camera_position = camera
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
    plotter.camera.zoom(1.10)


def render_wireframe(config: StageConfig, first: RenderedFrame, camera: object, out_png: pathlib.Path) -> None:
    mesh, _, _ = prepare_grid(first.vtk_path, 0.0)
    plotter = pv.Plotter(off_screen=True, window_size=(1600, 1200))
    plotter.add_mesh(mesh, style="wireframe", color=GARNET, line_width=0.85)
    add_common_scene(plotter, f"{config.title} - first-frame wireframe", camera)
    plotter.add_text(first.label, position="lower_right", font_size=10, color=CHARCOAL, font="courier")
    plotter.screenshot(str(out_png))
    plotter.close()


def render_field(
    config: StageConfig,
    frame: RenderedFrame,
    field: FieldChoice,
    clim: tuple[float, float],
    camera: object,
    warp_factor: float,
    out_png: pathlib.Path,
) -> None:
    mesh, _, vector_label = prepare_grid(frame.vtk_path, warp_factor)
    plotter = pv.Plotter(off_screen=True, window_size=(1600, 1200))
    scalars = mesh.point_data[field.scalar] if field.location == "point" else mesh.cell_data[field.scalar]
    plotter.add_mesh(
        mesh,
        scalars=scalars,
        preference=field.location,
        cmap=BRAND_RAMP,
        clim=clim,
        show_edges=True,
        edge_color=CHARCOAL,
        line_width=0.35,
        scalar_bar_args=scalar_bar(f"{field.title} ({field.units})"),
    )
    add_common_scene(plotter, f"{config.title} - {field.title}", camera)
    t = field_time(mesh)
    plotter.add_text(
        f"{frame.label} | t = {t:.4g} s\nwarp: {warp_factor:.3g}x by {vector_label}\n{config.verdict_note}",
        position="lower_right",
        font_size=8,
        color=CHARCOAL,
        font="courier",
    )
    plotter.screenshot(str(out_png))
    plotter.close()


def render_frame(
    config: StageConfig,
    frame: RenderedFrame,
    field: FieldChoice,
    clim: tuple[float, float],
    camera: object,
    warp_factor: float,
    out_png: pathlib.Path,
) -> None:
    mesh, _, vector_label = prepare_grid(frame.vtk_path, warp_factor)
    plotter = pv.Plotter(off_screen=True, window_size=(1920, 1080))
    scalars = mesh.point_data[field.scalar] if field.location == "point" else mesh.cell_data[field.scalar]
    plotter.add_mesh(
        mesh,
        scalars=scalars,
        preference=field.location,
        cmap=BRAND_RAMP,
        clim=clim,
        show_edges=False,
        scalar_bar_args=scalar_bar(f"{field.title} ({field.units})"),
    )
    add_common_scene(plotter, f"{config.title} - {field.title}", camera)
    plotter.add_text(
        f"{frame.label} | t = {field_time(mesh):.4g} s | warp {warp_factor:.3g}x by {vector_label}",
        position="lower_right",
        font_size=10,
        color=CHARCOAL,
        font="courier",
    )
    plotter.screenshot(str(out_png))
    plotter.close()


def render_double_frame(
    config: StageConfig,
    frame: RenderedFrame,
    disp_field: FieldChoice,
    stress_field: FieldChoice,
    ranges: dict[str, tuple[float, float]],
    camera: object,
    warp_factor: float,
    out_png: pathlib.Path,
) -> None:
    mesh, _, vector_label = prepare_grid(frame.vtk_path, warp_factor)
    plotter = pv.Plotter(off_screen=True, window_size=(1920, 1080), shape=(2, 1), border=False)
    for idx, field in enumerate([disp_field, stress_field]):
        plotter.subplot(idx, 0)
        scalars = mesh.point_data[field.scalar] if field.location == "point" else mesh.cell_data[field.scalar]
        plotter.add_mesh(
            mesh,
            scalars=scalars,
            preference=field.location,
            cmap=BRAND_RAMP,
            clim=ranges[field.scalar],
            show_edges=False,
            scalar_bar_args=scalar_bar(f"{field.title} ({field.units})"),
        )
        add_common_scene(plotter, f"{config.title} - {field.title}", camera, bounds=False)
        if idx == 1:
            plotter.add_text(
                f"{frame.label} | t = {field_time(mesh):.4g} s | warp {warp_factor:.3g}x by {vector_label}",
                position="lower_right",
                font_size=10,
                color=CHARCOAL,
                font="courier",
            )
    plotter.screenshot(str(out_png))
    plotter.close()


def selected_video_frames(frames: list[RenderedFrame]) -> list[RenderedFrame]:
    output_count = max(len(frames), VIDEO_FRAMES)
    return [frames[i % len(frames)] for i in range(output_count)]


def link_or_copy(src: pathlib.Path, dst: pathlib.Path) -> None:
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def ffmpeg_encode(frame_dir: pathlib.Path, out_mp4: pathlib.Path) -> None:
    subprocess.run(
        [
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
        ],
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )


def render_video(
    config: StageConfig,
    frames: list[RenderedFrame],
    field: FieldChoice,
    clim: tuple[float, float],
    camera: object,
    warp_factor: float,
    out_mp4: pathlib.Path,
) -> None:
    with tempfile.TemporaryDirectory(prefix=f"{config.prefix}_{field.scalar}_", dir=str(config.figures_dir)) as tmp:
        frame_dir = pathlib.Path(tmp)
        unique_paths: list[pathlib.Path] = []
        for idx, frame in enumerate(frames):
            out_png = frame_dir / f"unique_{idx:04d}.png"
            render_frame(config, frame, field, clim, camera, warp_factor, out_png)
            unique_paths.append(out_png)
        for idx, _frame in enumerate(selected_video_frames(frames)):
            link_or_copy(unique_paths[idx % len(unique_paths)], frame_dir / f"frame_{idx:04d}.png")
        ffmpeg_encode(frame_dir, out_mp4)


def render_double_video(
    config: StageConfig,
    frames: list[RenderedFrame],
    disp_field: FieldChoice,
    stress_field: FieldChoice,
    ranges: dict[str, tuple[float, float]],
    camera: object,
    warp_factor: float,
    out_mp4: pathlib.Path,
) -> None:
    with tempfile.TemporaryDirectory(prefix=f"{config.prefix}_doubleview_", dir=str(config.figures_dir)) as tmp:
        frame_dir = pathlib.Path(tmp)
        unique_paths: list[pathlib.Path] = []
        for idx, frame in enumerate(frames):
            out_png = frame_dir / f"unique_{idx:04d}.png"
            render_double_frame(config, frame, disp_field, stress_field, ranges, camera, warp_factor, out_png)
            unique_paths.append(out_png)
        for idx, _frame in enumerate(selected_video_frames(frames)):
            link_or_copy(unique_paths[idx % len(unique_paths)], frame_dir / f"frame_{idx:04d}.png")
        ffmpeg_encode(frame_dir, out_mp4)


def available_fields(path: pathlib.Path) -> list[FieldChoice]:
    mesh, _, vector_label = prepare_grid(path, 0.0)
    fields = [
        FieldChoice(
            scalar="richviz_vector_magnitude_mm",
            location="point",
            title=f"{vector_label} magnitude" if vector_label != "Displacement" else "displacement magnitude",
            units="mm",
            description=f"point {vector_label} magnitude in millimetres",
        )
    ]
    if "von_mises_mpa" in mesh.cell_data:
        fields.append(FieldChoice("von_mises_mpa", "cell", "von Mises stress", "MPa", "cell von Mises stress computed from 3D-element stress tensors"))
    if "eq_strain" in mesh.cell_data:
        fields.append(FieldChoice("eq_strain", "cell", "equivalent strain", "-", "cell von-Mises-equivalent strain computed from 3D-element strain tensors"))
    elif "plastic_strain" in mesh.cell_data:
        fields.append(FieldChoice("plastic_strain", "cell", "plastic strain", "-", "cell plastic strain scalar from the VTK output"))
    elif "damage_magnitude" in mesh.cell_data:
        fields.append(FieldChoice("damage_magnitude", "cell", "damage magnitude", "-", "cell damage-array magnitude from the VTK output"))
    return fields


def read_csv_dicts(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def as_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def history_for_stage(config: StageConfig) -> tuple[str, str, list[HistorySeries], str]:
    rows = read_csv_dicts(config.timeseries_path)
    if config.history_kind == "stage03":
        fem = []
        ref = []
        for row in rows:
            eps = as_float(row.get("eps_disp"))
            sigma = as_float(row.get("sigma_mean_Pa"))
            sigma_ref = as_float(row.get("sigma_ref_Pa"))
            if eps is not None and sigma is not None:
                fem.append((eps * 100.0, sigma / 1.0e6))
            if eps is not None and sigma_ref is not None:
                ref.append((eps * 100.0, sigma_ref / 1.0e6))
        return "gauge strain (%)", "gauge stress (MPa)", [
            HistorySeries("OpenRadioss mean stress", tuple(fem), GARNET),
            HistorySeries("closed-form load/area", tuple(ref), ATLANTIC),
        ], "Stress-strain points are extracted from final VTK frames for each load factor."
    if config.history_kind == "stage04":
        fem = []
        howland = []
        kirsch = []
        for row in rows:
            ntheta = as_float(row.get("n_theta"))
            kt = as_float(row.get("kt_fem"))
            if ntheta is not None and kt is not None:
                fem.append((ntheta, kt))
                howland.append((ntheta, 3.035))
                kirsch.append((ntheta, 3.0))
        return "circumferential divisions", "stress concentration Kt", [
            HistorySeries("OpenRadioss recovered Kt", tuple(fem), GARNET),
            HistorySeries("Howland target", tuple(howland), ATLANTIC),
            HistorySeries("Kirsch infinite plate", tuple(kirsch), HORSESHOE),
        ], "The primary gate is Ntheta 64 against the Howland finite-width target."
    if config.history_kind == "stage05":
        series = []
        for key, label, color in [
            ("F_coarse_N", "coarse", ATLANTIC),
            ("F_medium_N", "medium", GARNET),
            ("F_fine_N", "fine", CONGAREE),
        ]:
            pts = []
            for row in rows:
                disp = as_float(row.get("displacement_m"))
                force = as_float(row.get(key))
                if disp is not None and force is not None:
                    pts.append((disp * 1000.0, force / 1000.0))
            series.append(HistorySeries(label, tuple(pts), color))
        return "loaded-end displacement (mm)", "reaction force (kN)", series, "Force-displacement histories compare the coarse, medium, and fine LAW22 mesh sweep."
    if config.history_kind == "stage06":
        series = []
        for criterion, color in [("TSAIWU", GARNET), ("HASHIN", ATLANTIC), ("PUCK", CONGAREE)]:
            pts = []
            for row in rows:
                if row.get("kind") != "analytic_envelope" or row.get("criterion") != criterion:
                    continue
                path_id = as_float(row.get("path_id"))
                reported = as_float(row.get("reported_pa"))
                if path_id is not None and reported is not None:
                    pts.append((path_id, reported / 1.0e6))
            series.append(HistorySeries(criterion, tuple(pts), color))
        return "analytic path id", "reported strength (MPa)", series, "Canonical TYPE14 solver samples were not produced; this plot shows the analytic envelope rows recorded in timeseries.csv."
    numeric_keys = [key for key in rows[0].keys() if key != "stage" and as_float(rows[0].get(key)) is not None]
    if len(numeric_keys) < 2:
        return "row", "value", [HistorySeries("rows", tuple((i, i) for i, _ in enumerate(rows)), GARNET)], "No numeric scalar pair was available."
    x_key, y_key = numeric_keys[0], numeric_keys[1]
    pts = []
    for row in rows:
        x = as_float(row.get(x_key))
        y = as_float(row.get(y_key))
        if x is not None and y is not None:
            pts.append((x, y))
    return x_key, y_key, [HistorySeries(y_key, tuple(pts), GARNET)], "Automatically selected numeric columns from timeseries.csv."


def sample_points(points: tuple[tuple[float, float], ...], count: int = 110) -> tuple[tuple[float, float], ...]:
    if len(points) <= count:
        return points
    idx = np.linspace(0, len(points) - 1, count).round().astype(int)
    return tuple(points[int(i)] for i in idx)


def fmt_num(value: float) -> str:
    if abs(value) >= 1000.0 or (abs(value) < 0.01 and value != 0.0):
        return f"{value:.3g}"
    return f"{value:.4g}"


def chart_line(points: tuple[tuple[float, float], ...], xlim: tuple[float, float], ylim: tuple[float, float], rect: tuple[float, float, float, float]) -> str:
    x0, y0, width, height = rect
    xmin, xmax = xlim
    ymin, ymax = ylim
    pieces = []
    for x, y in sample_points(points):
        px = x0 + ((x - xmin) / (xmax - xmin)) * width if xmax != xmin else x0
        py = y0 + ((y - ymin) / (ymax - ymin)) * height if ymax != ymin else y0
        pieces.append(f"({px:.3f}, {py:.3f})")
    return ", ".join(pieces)


def write_history_typ(config: StageConfig) -> pathlib.Path:
    x_label, y_label, series, note = history_for_stage(config)
    all_points = [pt for item in series for pt in item.points]
    if not all_points:
        all_points = [(0.0, 0.0), (1.0, 1.0)]
    xs = [p[0] for p in all_points]
    ys = [p[1] for p in all_points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(0.0, min(ys)), max(ys)
    if math.isclose(xmin, xmax):
        xmax = xmin + 1.0
    if math.isclose(ymin, ymax):
        ymax = ymin + 1.0
    ypad = 0.05 * (ymax - ymin)
    ymax += ypad
    rect = (1.45, 1.15, 13.5, 6.9)
    lines = [
        '#import "@preview/cetz:0.3.4"',
        "",
        '#set page(width: 185mm, height: 118mm, margin: 10mm)',
        '#set text(font: "Libertinus Serif", size: 8.8pt, fill: rgb("#363636"))',
        '#let garnet = rgb("#73000A")',
        '#let atlantic = rgb("#466A9F")',
        '#let horseshoe = rgb("#65780B")',
        '#let congaree = rgb("#1F414D")',
        '#let charcoal = rgb("#363636")',
        '#let black10 = rgb("#ECECEC")',
        '#let white = rgb("#FFFFFF")',
        "",
        f"#align(center)[#text(size: 11pt, weight: \"bold\")[{config.title} History]]",
        "#v(2mm)",
        "#cetz.canvas(length: 1cm, {",
        "  import cetz.draw: *",
        f"  rect(({rect[0]:.3f}, {rect[1]:.3f}), ({rect[0] + rect[2]:.3f}, {rect[1] + rect[3]:.3f}), fill: white, stroke: charcoal + 0.65pt)",
    ]
    for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
        y = rect[1] + rect[3] * frac
        val = ymin + (ymax - ymin) * frac
        lines.extend(
            [
                f"  line(({rect[0]:.3f}, {y:.3f}), ({rect[0] + rect[2]:.3f}, {y:.3f}), stroke: black10 + 0.45pt)",
                f"  content(({rect[0] - 0.18:.3f}, {y:.3f}), [{fmt_num(val)}], anchor: \"east\")",
            ]
        )
    for item in series:
        if len(item.points) < 2:
            continue
        color_name = {
            GARNET: "garnet",
            ATLANTIC: "atlantic",
            HORSESHOE: "horseshoe",
            CONGAREE: "congaree",
        }.get(item.color, "garnet")
        lines.append(f"  line({chart_line(item.points, (xmin, xmax), (ymin, ymax), rect)}, stroke: {color_name} + 1.05pt)")
    legend_y = rect[1] + rect[3] + 0.60
    legend_x = rect[0] + 0.15
    for idx, item in enumerate(series):
        color_name = {
            GARNET: "garnet",
            ATLANTIC: "atlantic",
            HORSESHOE: "horseshoe",
            CONGAREE: "congaree",
        }.get(item.color, "garnet")
        x = legend_x + idx * 4.25
        lines.extend(
            [
                f"  line(({x:.3f}, {legend_y:.3f}), ({x + 0.45:.3f}, {legend_y:.3f}), stroke: {color_name} + 1.1pt)",
                f"  content(({x + 0.55:.3f}, {legend_y:.3f}), [{item.label}], anchor: \"west\")",
            ]
        )
    lines.extend(
        [
            f"  content(({rect[0]:.3f}, {rect[1] - 0.42:.3f}), [{fmt_num(xmin)}], anchor: \"north\")",
            f"  content(({rect[0] + rect[2]:.3f}, {rect[1] - 0.42:.3f}), [{fmt_num(xmax)}], anchor: \"north\")",
            f"  content(({rect[0] + rect[2] / 2:.3f}, {rect[1] - 0.72:.3f}), [{x_label}], anchor: \"north\")",
            f"  content(({rect[0] - 1.05:.3f}, {rect[1] + rect[3] / 2:.3f}), [{y_label}], angle: 90deg)",
            "})",
            "",
            f"#text(size: 8pt)[Verdict context: {config.verdict_note} {note} Source CSV: #raw(\"{config.timeseries_path.relative_to(config.figures_dir.parent).as_posix()}\").]",
            "",
        ]
    )
    out = config.figures_dir / f"{config.prefix}_history_disp.typ"
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def compile_typst(typ_path: pathlib.Path) -> pathlib.Path:
    pdf_path = typ_path.with_suffix(".pdf")
    subprocess.run(["typst", "compile", typ_path.name, pdf_path.name], cwd=str(typ_path.parent), check=True)
    return pdf_path


def write_readme(
    config: StageConfig,
    frames: list[RenderedFrame],
    final_frame: RenderedFrame,
    fields: list[FieldChoice],
    warp_factor: float,
    produced: set[str],
) -> pathlib.Path:
    results = read_results(config.stage_dir)
    verdict = str(results.get("verdict", "UNKNOWN"))
    final_time = field_time(pv.read(final_frame.vtk_path))
    disp_field = fields[0]
    stress_field = next((field for field in fields if field.scalar == "von_mises_mpa"), None)
    strain_field = next((field for field in fields if field.scalar in {"eq_strain", "plastic_strain", "damage_magnitude"}), None)
    paragraphs = [
        f"# Stage {config.stage:02d} Figure Gallery",
        "",
        f"The current `results/results.json` verdict is {verdict}. {config.verdict_note} The standard rich-visualization set reads existing outputs only: {config.run_note}.",
        "",
        f"`{config.prefix}_mesh_wireframe.png` uses the first available frame (`{frames[0].label}` from `{frames[0].run_label}`) reconstructed to reference coordinates when the `Displacement` vector is present. It is an iso-view garnet wireframe on white with axes in metres; verdict context is {verdict}.",
        "",
        f"`{config.prefix}_field_displacement.png` uses `{final_frame.label}` from `{final_frame.run_label}` at VTK time {final_time:.6g} s, warped {warp_factor:.3g}x by the selected vector field and coloured by {disp_field.description}. The colourbar is in millimetres, the mesh edges are charcoal, the axes are metres, and the verdict context is {verdict}.",
        "",
        f"`{config.prefix}_history_disp.typ` / `{config.prefix}_history_disp.pdf` is a Typst + CeTZ history plot built from `results/timeseries.csv` using brand-colour curves. The plotted scalar is stage-specific: stress-strain for Stage 03, Kt for Stage 04, force-displacement for Stage 05, and the recorded analytic failure envelope for Stage 06; verdict context is {verdict}.",
    ]
    if f"{config.prefix}_motion_displacement.mp4" in produced:
        paragraphs.extend(
            [
                "",
                f"`{config.prefix}_motion_displacement.mp4` is a 1920x1080 H.264 loop at 24 fps using all {len(frames)} rendered frame(s), repeated to 5.0 s. The mesh is warped {warp_factor:.3g}x and coloured by {disp_field.description} with the requested brand ramp; verdict context is {verdict}.",
            ]
        )
    if stress_field is not None and f"{config.prefix}_field_vonmises.png" in produced:
        paragraphs.extend(
            [
                "",
                f"`{config.prefix}_field_vonmises.png` uses `{final_frame.label}` at VTK time {final_time:.6g} s, the same vector warp, and cell von Mises stress in MPa computed from the 3D-element stress tensor. Axes are metres, edges are charcoal, and the verdict context is {verdict}.",
            ]
        )
    if stress_field is not None and f"{config.prefix}_motion_vonmises.mp4" in produced:
        paragraphs.extend(
            [
                "",
                f"`{config.prefix}_motion_vonmises.mp4` uses the same frame sequence and warp as the displacement motion video, coloured by cell von Mises stress in MPa with a fixed colour range over the 5.0 s loop; verdict context is {verdict}.",
            ]
        )
    if f"{config.prefix}_motion_doubleview.mp4" in produced:
        paragraphs.extend(
            [
                "",
                f"`{config.prefix}_motion_doubleview.mp4` is a synchronized 1920x1080 split view of the same frame sequence: vector magnitude in millimetres on the top panel and von Mises stress in MPa on the bottom panel. The lower label reports the source frame and VTK time; verdict context is {verdict}.",
            ]
        )
    if strain_field is not None and f"{config.prefix}_field_strain.png" in produced:
        paragraphs.extend(
            [
                "",
                f"`{config.prefix}_field_strain.png` uses `{final_frame.label}` at VTK time {final_time:.6g} s and the same vector warp, coloured by {strain_field.description}. The scalar is dimensionless, axes are metres, edges are charcoal, and verdict context is {verdict}.",
            ]
        )
    if config.stage == 6:
        paragraphs.extend(
            [
                "",
                "Stage 06 is INCONCLUSIVE because the required LAW25 + TYPE14 canonical solid path is rejected by the starter. These figures intentionally show the TYPE6/SOL_ORTH proxy outputs and analytic diagnostic rows that were actually produced; they are not evidence that the canonical laminate solid path passed.",
            ]
        )
    out = config.figures_dir / "README.md"
    out.write_text("\n".join(paragraphs) + "\n", encoding="utf-8")
    return out


def render_stage(config: StageConfig) -> None:
    if not config.frames:
        raise FileNotFoundError(f"stage {config.stage:02d} has no VTK or animation frames selected")
    pv.OFF_SCREEN = True
    config.figures_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"{config.prefix}_anim_vtk_") as tmp:
        scratch = pathlib.Path(tmp)
        frames = resolve_frames(config, scratch)
        paths = [frame.vtk_path for frame in frames]
        bounds = mesh_bounds_union(paths)
        warp_factor = automatic_warp_factor(paths, bounds)
        camera = camera_position(bounds)
        final_frame = frames[config.final_index]
        fields = available_fields(final_frame.vtk_path)
        disp_field = fields[0]
        stress_field = next((field for field in fields if field.scalar == "von_mises_mpa"), None)
        strain_field = next((field for field in fields if field.scalar in {"eq_strain", "plastic_strain", "damage_magnitude"}), None)
        ranges = {field.scalar: scalar_range(paths, field.scalar, field.location, warp_factor) for field in fields}
        produced: set[str] = set()

        out = config.figures_dir / f"{config.prefix}_mesh_wireframe.png"
        render_wireframe(config, frames[0], camera, out)
        produced.add(out.name)

        out = config.figures_dir / f"{config.prefix}_field_displacement.png"
        render_field(config, final_frame, disp_field, ranges[disp_field.scalar], camera, warp_factor, out)
        produced.add(out.name)

        typ_path = write_history_typ(config)
        compile_typst(typ_path)
        produced.add(typ_path.name)
        produced.add(typ_path.with_suffix(".pdf").name)

        if len(frames) > 1:
            out = config.figures_dir / f"{config.prefix}_motion_displacement.mp4"
            render_video(config, frames, disp_field, ranges[disp_field.scalar], camera, warp_factor, out)
            produced.add(out.name)

            if stress_field is not None:
                out = config.figures_dir / f"{config.prefix}_field_vonmises.png"
                render_field(config, final_frame, stress_field, ranges[stress_field.scalar], camera, warp_factor, out)
                produced.add(out.name)

                out = config.figures_dir / f"{config.prefix}_motion_vonmises.mp4"
                render_video(config, frames, stress_field, ranges[stress_field.scalar], camera, warp_factor, out)
                produced.add(out.name)

                out = config.figures_dir / f"{config.prefix}_motion_doubleview.mp4"
                render_double_video(config, frames, disp_field, stress_field, ranges, camera, warp_factor, out)
                produced.add(out.name)

        if strain_field is not None:
            out = config.figures_dir / f"{config.prefix}_field_strain.png"
            render_field(config, final_frame, strain_field, ranges[strain_field.scalar], camera, warp_factor, out)
            produced.add(out.name)

        write_readme(config, frames, final_frame, fields, warp_factor, produced)
        print(f"rendered stage {config.stage:02d}: {', '.join(sorted(produced))}")


def eligible_stages() -> list[int]:
    stages = []
    for sdir in sorted((ROOT / "tests").glob("stage_[0-9][0-9]_*")):
        try:
            stage = int(sdir.name.split("_")[1])
        except Exception:
            continue
        if stage <= 2:
            continue
        runs = sdir / "runs"
        if not runs.exists():
            continue
        if any(runs.rglob("*.vtk")) or any(runs.rglob("*.gz")):
            stages.append(stage)
    return stages


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=int, action="append", help="stage number to render; may be repeated")
    parser.add_argument("--all", action="store_true", help="render all eligible stages after stage 02")
    args = parser.parse_args(argv)
    if shutil.which("ffmpeg") is None:
        raise SystemExit("ffmpeg not found")
    if shutil.which("typst") is None:
        raise SystemExit("typst not found")
    stages = args.stage or []
    if args.all:
        stages = eligible_stages()
    if not stages:
        raise SystemExit("choose --stage or --all")
    for stage in stages:
        config = config_for_stage(stage)
        if config is None:
            print(f"stage {stage:02d}: no eligible VTK or animation frames")
            continue
        render_stage(config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
