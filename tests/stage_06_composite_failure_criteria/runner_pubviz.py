"""Publication animation rerun for Stage 06.

Clones the LAW12 + TYPE6 TSAIWU overdrive probe into
``runs/stage06_TSAIWU_law12_type6_pubviz`` and changes only ``/ANIM/DT`` so the
solver writes 60 animation frames.
"""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.pubviz_regenerate import run_pubviz  # noqa: E402


if __name__ == "__main__":
    run_pubviz(pathlib.Path(__file__).resolve().parent, "stage06_TSAIWU_law12_type6", 60, 16, False)
