"""Command-line entry point for the clean-room AP-PLY geometry port."""

from __future__ import annotations

import argparse
from pathlib import Path

from kok_geom.io import generate_mesh


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="kok_geom")
    parser.add_argument("--config", required=True, help="Path to kok_config.json")
    parser.add_argument(
        "--out",
        default=None,
        help="Output .msh path or output directory. Defaults to config.output paths.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    msh_path, orientations_path = generate_mesh(args.config, args.out)
    print(f"wrote {Path(msh_path)}")
    print(f"wrote {Path(orientations_path)}")
    return 0
