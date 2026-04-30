"""Render a stage 16 V50 shot animation to MP4.

Converts the OpenRadioss A###[.gz] animation frames in a shot run dir to
VTK via anim_to_vtk_linux64_gf, then renders each frame with PyVista
headless OSMesa coloured by displacement magnitude, stitches with ffmpeg.

Usage:
    python tools/stage16_shot_video.py shot_00_m0.80_v80.0
    python tools/stage16_shot_video.py shot_00_m0.80_v80.0 --output /custom/path.mp4
"""
from __future__ import annotations

import argparse
import gzip
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

ROOT = pathlib.Path(__file__).resolve().parents[1]
OR_ROOT = pathlib.Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
ANIM_TO_VTK = OR_ROOT / "exec" / "anim_to_vtk_linux64_gf"
STAGE_DIR = ROOT / "tests" / "stage_16_PW_panel_ballistic"
SWEEP_DIR = STAGE_DIR / "runs" / "stage_16_post_m5" / "phase_c1_v50_sweep"
FIGURES_DIR = STAGE_DIR / "figures"
ANIM_STEM = "stage16_phase_b_single_shot"

GARNET = "#73000A"
ATLANTIC = "#466A9F"
HORSESHOE = "#65780B"
SANDSTORM = "#FFF2E3"
CHARCOAL = "#363636"
WHITE = "#FFFFFF"

FRAME_SIZE = (1920, 1080)
FPS = 15
SIM_DURATION_S = 5.0e-4


def hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for c in candidates:
        if pathlib.Path(c).exists():
            return ImageFont.truetype(c, size)
    return ImageFont.load_default()


def gunzip_keep(src: pathlib.Path, dst: pathlib.Path) -> None:
    if dst.exists() and dst.stat().st_mtime >= src.stat().st_mtime:
        return
    with gzip.open(src, "rb") as f_in, dst.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def radioss_env() -> dict[str, str]:
    env = os.environ.copy()
    env["OR"] = str(OR_ROOT)
    env["RAD_CFG_PATH"] = str(OR_ROOT / "hm_cfg_files")
    env["RAD_H3D_PATH"] = str(OR_ROOT / "extlib" / "h3d" / "lib" / "linux64")
    reader = str(OR_ROOT / "extlib" / "hm_reader" / "linux64")
    env["LD_LIBRARY_PATH"] = reader + ":" + env.get("LD_LIBRARY_PATH", "")
    return env


def convert_one(anim_path: pathlib.Path, env: dict[str, str]) -> pathlib.Path:
    vtk = anim_path.with_suffix(anim_path.suffix + ".vtk")
    if vtk.exists() and vtk.stat().st_mtime >= anim_path.stat().st_mtime:
        return vtk
    proc = subprocess.run(
        [str(ANIM_TO_VTK), str(anim_path)],
        cwd=str(anim_path.parent),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"anim_to_vtk failed for {anim_path}: {proc.stdout[-500:]}")
    if proc.stdout.lstrip().startswith("# vtk"):
        vtk.write_text(proc.stdout, encoding="utf-8")
    elif not vtk.exists():
        raise FileNotFoundError(f"anim_to_vtk produced no VTK for {anim_path}")
    return vtk


def find_frames(run_dir: pathlib.Path) -> list[pathlib.Path]:
    frame_re = re.compile(rf"^{re.escape(ANIM_STEM)}A\d{{3}}(?:\.gz)?$")
    frames: list[pathlib.Path] = []
    for source in sorted(p for p in run_dir.iterdir() if frame_re.match(p.name)):
        if source.suffix == ".gz":
            anim = source.with_suffix("")
            gunzip_keep(source, anim)
        else:
            anim = source
        frames.append(anim)
    return frames


def reference_bounds(vtks: list[pathlib.Path]) -> tuple[float, float, float, float, float, float]:
    bounds = None
    for vtk in vtks[:: max(1, len(vtks) // 8)]:
        mesh = pv.read(str(vtk))
        b = mesh.bounds
        if bounds is None:
            bounds = list(b)
        else:
            for i in range(0, 6, 2):
                bounds[i] = min(bounds[i], b[i])
            for i in range(1, 6, 2):
                bounds[i] = max(bounds[i], b[i])
    return tuple(bounds) if bounds else (-0.05, 0.05, -0.05, 0.05, -0.01, 0.01)


def render_frame(vtk_path: pathlib.Path, frame_idx: int, n_frames: int, out_png: pathlib.Path,
                 disp_clim: tuple[float, float], bounds: tuple[float, float, float, float, float, float]) -> None:
    mesh = pv.read(str(vtk_path))
    if "Displacement" in mesh.point_data:
        disp = np.asarray(mesh.point_data["Displacement"])
        mesh.point_data["disp_mag_mm"] = np.linalg.norm(disp, axis=1) * 1000.0
        scalar_field = "disp_mag_mm"
    elif "Velocity" in mesh.point_data:
        vel = np.asarray(mesh.point_data["Velocity"])
        mesh.point_data["vel_mag_mps"] = np.linalg.norm(vel, axis=1)
        scalar_field = "vel_mag_mps"
    else:
        scalar_field = None

    plotter = pv.Plotter(off_screen=True, window_size=list(FRAME_SIZE))
    plotter.background_color = SANDSTORM
    if scalar_field is not None:
        plotter.add_mesh(
            mesh,
            scalars=scalar_field,
            cmap="inferno",
            clim=disp_clim,
            show_edges=False,
            scalar_bar_args={"title": "Displacement (mm)", "color": CHARCOAL, "title_font_size": 18, "label_font_size": 14, "n_labels": 5},
        )
    else:
        plotter.add_mesh(mesh, color=GARNET, show_edges=False)

    cx = 0.5 * (bounds[0] + bounds[1])
    cy = 0.5 * (bounds[2] + bounds[3])
    cz = 0.5 * (bounds[4] + bounds[5])
    span = max(bounds[1] - bounds[0], bounds[3] - bounds[2])
    cam_dist = span * 1.6
    plotter.camera_position = [(cx + cam_dist * 0.7, cy - cam_dist * 0.7, cz + cam_dist * 0.6),
                                (cx, cy, cz),
                                (0.0, 0.0, 1.0)]
    plotter.screenshot(str(out_png), window_size=list(FRAME_SIZE))
    plotter.close()
    overlay_text(out_png, frame_idx, n_frames)


def overlay_text(png: pathlib.Path, frame_idx: int, n_frames: int) -> None:
    img = Image.open(png).convert("RGB")
    draw = ImageDraw.Draw(img)
    title_font = font(40, bold=True)
    info_font = font(28)
    time_us = (frame_idx / max(1, n_frames - 1)) * SIM_DURATION_S * 1e6
    draw.rectangle((0, 0, FRAME_SIZE[0], 90), fill=hex_to_rgb(CHARCOAL))
    draw.text((30, 18), "Stage 16 — V50 shot 0 — V_i = 80 m/s (0.80 multiplier)", font=title_font, fill=hex_to_rgb(WHITE))
    bar_w = FRAME_SIZE[0] - 60
    draw.rectangle((30, FRAME_SIZE[1] - 70, 30 + bar_w, FRAME_SIZE[1] - 40), outline=hex_to_rgb(CHARCOAL), width=2)
    fill_w = int(bar_w * (frame_idx / max(1, n_frames - 1)))
    draw.rectangle((30, FRAME_SIZE[1] - 70, 30 + fill_w, FRAME_SIZE[1] - 40), fill=hex_to_rgb(GARNET))
    draw.text((30, FRAME_SIZE[1] - 32), f"t = {time_us:6.1f} us   |   frame {frame_idx + 1}/{n_frames}   |   V_r = 0 (arrested)", font=info_font, fill=hex_to_rgb(CHARCOAL))
    img.save(png)


def stitch_video(frames_dir: pathlib.Path, out_mp4: pathlib.Path) -> None:
    pattern = str(frames_dir / "frame_%04d.png")
    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(FPS),
        "-i", pattern,
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        "-crf", "20",
        str(out_mp4),
    ]
    subprocess.run(cmd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("shot_dir", type=str, help="shot subdirectory name under phase_c1_v50_sweep/")
    parser.add_argument("--output", type=pathlib.Path, default=None)
    args = parser.parse_args()

    pv.OFF_SCREEN = True
    run_dir = SWEEP_DIR / args.shot_dir
    if not run_dir.is_dir():
        print(f"shot dir not found: {run_dir}", file=sys.stderr)
        return 1
    out_mp4 = args.output or FIGURES_DIR / f"stage_16_{args.shot_dir}_motion.mp4"
    out_mp4.parent.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] gunzip + index frames in {run_dir}")
    anims = find_frames(run_dir)
    if not anims:
        print(f"no animation frames found in {run_dir}", file=sys.stderr)
        return 1
    print(f"      found {len(anims)} animation frames")

    print(f"[2/4] anim_to_vtk conversion (cached if up to date)")
    env = radioss_env()
    vtks: list[pathlib.Path] = []
    for i, anim in enumerate(anims):
        vtk = convert_one(anim, env)
        vtks.append(vtk)
        if (i + 1) % 10 == 0 or i == len(anims) - 1:
            print(f"      converted {i + 1}/{len(anims)}")

    print(f"[3/4] computing scalar bounds + rendering frames")
    bounds = reference_bounds(vtks)
    disp_clims = []
    for vtk in vtks[:: max(1, len(vtks) // 8)]:
        m = pv.read(str(vtk))
        if "Displacement" in m.point_data:
            disp = np.asarray(m.point_data["Displacement"])
            disp_clims.append(np.linalg.norm(disp, axis=1).max() * 1000.0)
    disp_clim = (0.0, max(disp_clims) if disp_clims else 1.0)
    print(f"      bounds={bounds} disp_clim={disp_clim}")

    with tempfile.TemporaryDirectory(prefix="stage16_shot_frames_") as tmp:
        tmpdir = pathlib.Path(tmp)
        for i, vtk in enumerate(vtks):
            png = tmpdir / f"frame_{i:04d}.png"
            render_frame(vtk, i, len(vtks), png, disp_clim, bounds)
            if (i + 1) % 10 == 0 or i == len(vtks) - 1:
                print(f"      rendered {i + 1}/{len(vtks)}")
        print(f"[4/4] ffmpeg stitch -> {out_mp4}")
        stitch_video(tmpdir, out_mp4)

    print(f"DONE: {out_mp4} ({out_mp4.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
