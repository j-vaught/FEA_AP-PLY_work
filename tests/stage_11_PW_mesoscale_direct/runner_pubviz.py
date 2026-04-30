"""Publication animation reruns for Stage 11.

Clones the three passing Stage 11 decks into sibling ``*_pubviz`` run
directories and changes only ``/ANIM/DT`` so the solver writes 60 animation
frames per load case.
"""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.pubviz_regenerate import run_pubviz  # noqa: E402


if __name__ == "__main__":
    stage_dir = pathlib.Path(__file__).resolve().parent
    for source_run in ("stage11_axial_x", "stage11_transverse_y", "stage11_shear_xy"):
        run_pubviz(stage_dir, source_run, 60, 16, False)
