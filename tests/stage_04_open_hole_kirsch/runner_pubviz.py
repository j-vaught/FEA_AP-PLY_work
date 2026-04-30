"""Publication animation rerun for Stage 04.

Clones the existing passing Ntheta 64 deck into ``runs/stage04_ntheta064_pubviz``
and changes only ``/ANIM/DT`` so the solver writes 60 animation frames.
"""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.pubviz_regenerate import run_pubviz  # noqa: E402


if __name__ == "__main__":
    run_pubviz(pathlib.Path(__file__).resolve().parent, "stage04_ntheta064", 60, 16, False)
