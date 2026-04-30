"""Clone an existing OpenRadioss run deck for publication animation output.

The helper copies an already-passing starter/engine deck into a sibling
``*_pubviz`` run directory, changes only the job stem and ``/ANIM/DT`` cadence,
runs starter/engine, and converts every animation frame to VTK. It deliberately
does not touch the verification runner outputs or ``results.json`` files.
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


ROOT = pathlib.Path(__file__).resolve().parents[1]
OR_ROOT = pathlib.Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"
ENGINE = OR_ROOT / "exec" / "engine_linux64_gf"
ANIM_TO_VTK = OR_ROOT / "exec" / "anim_to_vtk_linux64_gf"


def fmt_f(*values: float) -> str:
    return "".join(f"{value:20.12g}" for value in values)


def radioss_env() -> dict[str, str]:
    env = os.environ.copy()
    env["OR"] = str(OR_ROOT)
    env["RAD_CFG_PATH"] = str(OR_ROOT / "hm_cfg_files")
    env["RAD_H3D_PATH"] = str(OR_ROOT / "extlib" / "h3d" / "lib" / "linux64")
    reader = str(OR_ROOT / "extlib" / "hm_reader" / "linux64")
    env["LD_LIBRARY_PATH"] = reader + ":" + env.get("LD_LIBRARY_PATH", "")
    return env


def parse_run_time(engine_text: str) -> float:
    lines = engine_text.splitlines()
    for idx, line in enumerate(lines):
        if line.strip().startswith("/RUN/"):
            for follow in lines[idx + 1 :]:
                stripped = follow.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                token = stripped.split()[0].replace("D", "E")
                return float(token)
    raise ValueError("could not find /RUN time in engine deck")


def set_anim_dt(engine_text: str, frame_count: int) -> str:
    run_time = parse_run_time(engine_text)
    lines = engine_text.splitlines()
    for idx, line in enumerate(lines):
        if line.strip() != "/ANIM/DT":
            continue
        insert = idx + 1
        while insert < len(lines) and lines[insert].strip().startswith("#"):
            insert += 1
        if insert >= len(lines):
            raise ValueError("/ANIM/DT has no data line")
        lines[insert] = fmt_f(0.0, run_time / frame_count)
        return "\n".join(lines) + "\n"
    raise ValueError("engine deck has no /ANIM/DT block")


def replace_stem(text: str, source_stem: str, pubviz_stem: str) -> str:
    return text.replace(source_stem, pubviz_stem)


def run_command(argv: list[str], cwd: pathlib.Path, log: pathlib.Path, env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    with log.open("a", encoding="utf-8") as fh:
        fh.write("$ " + " ".join(argv) + "\n")
    proc = subprocess.run(
        argv,
        cwd=str(cwd),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    with log.open("a", encoding="utf-8") as fh:
        fh.write(proc.stdout)
        fh.write(f"\n[exit {proc.returncode}]\n")
    if proc.stdout:
        if len(proc.stdout) <= 4000:
            print(proc.stdout, end="", flush=True)
        else:
            print(f"[captured {len(proc.stdout)} bytes in {log}]", flush=True)
    return proc


def gunzip_keep(src: pathlib.Path, dst: pathlib.Path) -> None:
    with gzip.open(src, "rb") as f_in, dst.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def convert_anim(anim: pathlib.Path, env: dict[str, str], log: pathlib.Path) -> pathlib.Path:
    proc = run_command([str(ANIM_TO_VTK), str(anim)], anim.parent, log, env)
    if proc.returncode != 0:
        raise RuntimeError(f"anim_to_vtk failed for {anim}")
    vtk = anim.with_suffix(anim.suffix + ".vtk")
    if proc.stdout.lstrip().startswith("# vtk"):
        vtk.write_text(proc.stdout, encoding="utf-8")
    elif not vtk.exists():
        candidates = sorted(anim.parent.glob(f"{anim.name}*.vtk")) + sorted(anim.parent.glob("*.vtk"))
        if not candidates:
            raise FileNotFoundError(f"anim_to_vtk produced no VTK for {anim}")
        if candidates[0] != vtk:
            candidates[0].replace(vtk)
    return vtk


def convert_all_animation(workdir: pathlib.Path, stem: str, env: dict[str, str], log: pathlib.Path) -> list[pathlib.Path]:
    converted: list[pathlib.Path] = []
    frame_re = re.compile(rf"^{re.escape(stem)}A\d{{3}}(?:\.gz)?$")
    for source in sorted(p for p in workdir.iterdir() if frame_re.match(p.name)):
        if source.suffix == ".gz":
            anim = source.with_suffix("")
            if not anim.exists() or source.stat().st_mtime > anim.stat().st_mtime:
                gunzip_keep(source, anim)
        else:
            anim = source
        vtk = anim.with_suffix(anim.suffix + ".vtk")
        if vtk.exists() and vtk.stat().st_mtime >= anim.stat().st_mtime:
            converted.append(vtk)
            continue
        converted.append(convert_anim(anim, env, log))
    return converted


def source_decks(stage_dir: pathlib.Path, source_run: str) -> tuple[pathlib.Path, pathlib.Path]:
    source_dir = stage_dir / "runs" / source_run
    if source_dir.is_dir():
        starter = source_dir / f"{source_run}_0000.rad"
        engine = source_dir / f"{source_run}_0001.rad"
    else:
        starter = stage_dir / "runs" / f"{source_run}_0000.rad"
        engine = stage_dir / "runs" / f"{source_run}_0001.rad"
    if not starter.exists() or not engine.exists():
        raise FileNotFoundError(f"missing starter/engine pair for {source_run} under {stage_dir / 'runs'}")
    return starter, engine


def clone_decks(stage_dir: pathlib.Path, source_run: str, output_dir: pathlib.Path, frame_count: int) -> tuple[str, str]:
    source_stem = source_run
    pubviz_stem = output_dir.name
    starter_src, engine_src = source_decks(stage_dir, source_run)
    output_dir.mkdir(parents=True, exist_ok=True)
    starter_text = replace_stem(starter_src.read_text(encoding="utf-8"), source_stem, pubviz_stem)
    engine_text = replace_stem(engine_src.read_text(encoding="utf-8"), source_stem, pubviz_stem)
    engine_text = set_anim_dt(engine_text, frame_count)
    (output_dir / f"{pubviz_stem}_0000.rad").write_text(starter_text, encoding="utf-8")
    (output_dir / f"{pubviz_stem}_0001.rad").write_text(engine_text, encoding="utf-8")
    return source_stem, pubviz_stem


def run_pubviz(stage_dir: pathlib.Path, source_run: str, frame_count: int, threads: int, force: bool) -> pathlib.Path:
    output_dir = stage_dir / "runs" / f"{source_run}_pubviz"
    existing = sorted(output_dir.glob(f"{output_dir.name}A[0-9][0-9][0-9].vtk"))
    if len(existing) >= frame_count and not force:
        print(f"{output_dir}: existing {len(existing)} VTK frames, skipping solver")
        return output_dir
    if output_dir.exists() and force:
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log = output_dir / "runner_pubviz.log"
    log.write_text("", encoding="utf-8")
    _source_stem, pubviz_stem = clone_decks(stage_dir, source_run, output_dir, frame_count)
    env = radioss_env()
    starter = output_dir / f"{pubviz_stem}_0000.rad"
    engine = output_dir / f"{pubviz_stem}_0001.rad"
    proc = run_command([str(STARTER), "-i", str(starter), "-nt", str(threads)], output_dir, log, env)
    if proc.returncode != 0:
        raise RuntimeError(f"starter failed for {pubviz_stem}")
    proc = run_command([str(ENGINE), "-i", str(engine), "-nt", str(threads)], output_dir, log, env)
    if proc.returncode != 0:
        raise RuntimeError(f"engine failed for {pubviz_stem}")
    converted = convert_all_animation(output_dir, pubviz_stem, env, log)
    if len(converted) < frame_count:
        print(f"{pubviz_stem}: solver emitted {len(converted)} VTK frames; requested {frame_count}")
        return output_dir
    print(f"{pubviz_stem}: converted {len(converted)} VTK frames")
    return output_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-dir", required=True, type=pathlib.Path)
    parser.add_argument("--source-run", required=True)
    parser.add_argument("--frames", type=int, default=60)
    parser.add_argument("--threads", type=int, default=int(os.environ.get("RAD_NT", "16")))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    run_pubviz(args.stage_dir.resolve(), args.source_run, args.frames, args.threads, args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
