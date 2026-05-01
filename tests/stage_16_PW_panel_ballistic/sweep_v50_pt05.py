"""Stage 16 focused 3-shot retry with PTHICKFAIL=0.5.

The original Phase C1 sweep at PTHICKFAIL=1.0 returned ALL_ARRESTED for the
0.80-1.20 multiplier range (V_r = 0 at every velocity from 80 to 120 m/s).
Root cause: the strictest erosion criterion (entire element must fail at
all integration points) combined with under-resolved damage-zone meshing
(3 mm) prevents the panel from ever shedding material under impact.

This script runs the upper half of the bracket (V=100, 110, 120 m/s) with
PTHICKFAIL=0.5 (half of the integration points must fail before erosion),
which is the cheapest single-knob change that should let the panel actually
break through. If V=120 perforates and V=100 still arrests, the V50 bracket
is (100, 120) m/s with the new model.

Usage:
    python tests/stage_16_PW_panel_ballistic/sweep_v50_pt05.py
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent.parent
RUNNER = THIS_DIR / "runner.py"
GEOM_DIR = THIS_DIR / "geometry"
RUNS_DIR = THIS_DIR / "runs" / "stage_16_post_m5" / "phase_c2_v50_sweep_pt05"

V50_PUBLISHED_M_S = 100.0
MULTIPLIERS = (1.00, 1.10, 1.20)
THREADS = 32
PTHICKFAIL = 0.5
MESH = GEOM_DIR / "p1_twt_75.msh"
ORIENTATIONS = GEOM_DIR / "orientations.json"


def run_shot(idx: int, multiplier: float) -> dict[str, Any]:
    velocity = multiplier * V50_PUBLISHED_M_S
    run_dir = RUNS_DIR / f"shot_{idx:02d}_m{multiplier:.2f}_v{velocity:.1f}"
    summary_path = run_dir / "phase_b_summary.json"
    if summary_path.exists():
        existing = json.loads(summary_path.read_text())
        if existing.get("engine_rc") == 0 and existing.get("residual_velocity_m_s") is not None:
            print(f"[skip] shot {idx} multiplier {multiplier} already complete")
            return existing
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(RUNNER),
        "--mesh", str(MESH),
        "--orientations", str(ORIENTATIONS),
        "--run-dir", str(run_dir),
        "--velocity", f"{velocity:.4f}",
        "--threads", str(THREADS),
        "--pthickfail", f"{PTHICKFAIL:.4f}",
    ]
    print(f"[shot {idx}] multiplier={multiplier:.2f} velocity={velocity:.2f} m/s pthickfail={PTHICKFAIL} -> {run_dir}")
    start = time.time()
    rc = subprocess.call(cmd, cwd=str(REPO_ROOT))
    elapsed = time.time() - start
    print(f"[shot {idx}] rc={rc} elapsed={elapsed/3600:.2f} h")
    if not summary_path.exists():
        return {"idx": idx, "multiplier": multiplier, "velocity_m_s": velocity, "rc": rc, "summary_missing": True}
    summary = json.loads(summary_path.read_text())
    summary["idx"] = idx
    summary["multiplier"] = multiplier
    summary["wrapper_elapsed_h"] = elapsed / 3600.0
    summary["pthickfail"] = PTHICKFAIL
    return summary


def main() -> int:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"V50 retry sweep starts at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Multipliers: {MULTIPLIERS}  PTHICKFAIL={PTHICKFAIL}")
    shots = []
    for idx, mult in enumerate(MULTIPLIERS):
        shot = run_shot(idx, mult)
        shots.append(shot)
        intermediate = RUNS_DIR / "v50_sweep_intermediate.json"
        intermediate.write_text(json.dumps(shots, indent=2, default=str))
    summary = {
        "stage": 16,
        "name": "PW_panel_ballistic_v50_pt05_retry",
        "pthickfail": PTHICKFAIL,
        "shots": [
            {
                "idx": s.get("idx"),
                "multiplier": s.get("multiplier"),
                "velocity_m_s": s.get("impact_velocity_m_s", s.get("velocity_m_s")),
                "residual_velocity_m_s": s.get("residual_velocity_m_s"),
                "engine_rc": s.get("engine_rc"),
                "wrapper_elapsed_h": s.get("wrapper_elapsed_h"),
            }
            for s in shots
        ],
    }
    summary_path = RUNS_DIR / "v50_sweep_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nRetry sweep complete. Summary at {summary_path}")
    for s in summary["shots"]:
        print(f"  shot {s['idx']:2d}  V_i={s['velocity_m_s']:6.1f} m/s  V_r={s['residual_velocity_m_s']}  rc={s['engine_rc']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
