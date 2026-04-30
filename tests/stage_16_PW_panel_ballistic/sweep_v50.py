"""Stage 16 V50 sweep on the post-Kok-M5 75 mm reduced-section deck.

Drives the existing runner.py for seven shots at multipliers
(0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20) of V50_published (rough 100 m/s).
Each shot runs the full 0.5 ms simulation at -nt 32 with the 207k-elem
75 mm M5 deck, parses residual velocity from T01, fits Recht-Ipson 1963,
computes simulated V50, and writes the aggregate verdict to results.json.

Resume-safe: skips shots whose phase_b_summary.json already records a
completed engine (engine_rc == 0).

Usage:
    python tests/stage_16_PW_panel_ballistic/sweep_v50.py
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
RUNS_DIR = THIS_DIR / "runs" / "stage_16_post_m5" / "phase_c1_v50_sweep"
RESULTS_JSON = THIS_DIR / "results" / "results.json"

V50_PUBLISHED_M_S = 100.0
V50_PUBLISHED_STATUS = "UNVERIFIED"
MULTIPLIERS = (0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20)
THREADS = 32
MESH = GEOM_DIR / "p1_twt_75.msh"
ORIENTATIONS = GEOM_DIR / "orientations.json"
TOLERANCE = 0.07
TOLERANCE_BC_DEGRADED = 0.09


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
    ]
    print(f"[shot {idx}] multiplier={multiplier:.2f} velocity={velocity:.2f} m/s -> {run_dir}")
    print(f"          cmd: {' '.join(cmd)}")
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
    return summary


def recht_ipson_v50(velocities: list[float], residuals: list[float]) -> tuple[float | None, float | None, str]:
    """Fit V_r = a * (V_i^p - V50^p)^(1/p) with p=2 (Recht-Ipson 1963).

    Returns (V50, RMSE, status). For shots where V_i <= V50 the residual is 0
    (projectile arrested); those points are used to bracket V50 from below.
    """
    perforating = [(v, r) for v, r in zip(velocities, residuals) if r is not None and r > 1.0]
    arrested = [(v, r) for v, r in zip(velocities, residuals) if r is None or r <= 1.0]
    if not perforating:
        return None, None, "ALL_ARRESTED"
    if not arrested:
        return None, None, "ALL_PERFORATED"
    v_lo = max(v for v, _ in arrested)
    v_hi = min(v for v, _ in perforating)
    p = 2.0
    best_v50 = None
    best_rmse = float("inf")
    for v50 in [v_lo + (v_hi - v_lo) * t / 200.0 for t in range(201)]:
        if v50 <= 0:
            continue
        residuals_pred = []
        residuals_obs = []
        for v, r in perforating:
            if v <= v50:
                continue
            pred = (v**p - v50**p) ** (1.0 / p)
            residuals_pred.append(pred)
            residuals_obs.append(r)
        if not residuals_pred:
            continue
        rmse = math.sqrt(sum((p_ - o) ** 2 for p_, o in zip(residuals_pred, residuals_obs)) / len(residuals_pred))
        if rmse < best_rmse:
            best_rmse = rmse
            best_v50 = v50
    return best_v50, best_rmse, "FIT_OK"


def main() -> int:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"V50 sweep starts at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"V50_published = {V50_PUBLISHED_M_S} m/s ({V50_PUBLISHED_STATUS})")
    print(f"Mesh: {MESH}")
    print(f"Orientations: {ORIENTATIONS}")
    shots = []
    for idx, mult in enumerate(MULTIPLIERS):
        shot = run_shot(idx, mult)
        shots.append(shot)
        intermediate = RUNS_DIR / "v50_sweep_intermediate.json"
        intermediate.write_text(json.dumps(shots, indent=2, default=str))
    velocities = [s.get("impact_velocity_m_s", s.get("velocity_m_s", float("nan"))) for s in shots]
    residuals = [s.get("residual_velocity_m_s") for s in shots]
    v50_simulated, rmse, fit_status = recht_ipson_v50(velocities, residuals)
    relative_error = None
    verdict = "INCONCLUSIVE"
    if v50_simulated is not None and V50_PUBLISHED_STATUS == "VERIFIED":
        relative_error = abs(v50_simulated - V50_PUBLISHED_M_S) / V50_PUBLISHED_M_S
        verdict = "PASS" if relative_error <= TOLERANCE else ("INCONCLUSIVE" if relative_error <= TOLERANCE_BC_DEGRADED else "FAIL")
    sweep_summary = {
        "stage": 16,
        "name": "PW_panel_ballistic_v50_sweep",
        "fit_status": fit_status,
        "v50_simulated_m_s": v50_simulated,
        "v50_published_m_s": V50_PUBLISHED_M_S,
        "v50_published_status": V50_PUBLISHED_STATUS,
        "relative_error": relative_error,
        "tolerance": TOLERANCE,
        "rmse_residual_velocity_m_s": rmse,
        "verdict_modulo_published_v50": verdict,
        "shots": [
            {
                "idx": s.get("idx"),
                "multiplier": s.get("multiplier"),
                "velocity_m_s": velocities[i],
                "residual_velocity_m_s": residuals[i],
                "engine_rc": s.get("engine_rc"),
                "wrapper_elapsed_h": s.get("wrapper_elapsed_h"),
            }
            for i, s in enumerate(shots)
        ],
    }
    summary_path = RUNS_DIR / "v50_sweep_summary.json"
    summary_path.write_text(json.dumps(sweep_summary, indent=2, default=str))
    print(f"\nSweep complete. Summary at {summary_path}")
    print(f"  V50_simulated: {v50_simulated} m/s  (fit_status={fit_status})")
    print(f"  V50_published: {V50_PUBLISHED_M_S} m/s ({V50_PUBLISHED_STATUS})")
    print(f"  Verdict (modulo published V50): {verdict}")
    if RESULTS_JSON.exists():
        results = json.loads(RESULTS_JSON.read_text())
        metrics = results.setdefault("metrics", {})
        metrics["v50_sweep_reached"] = True
        metrics["v50_simulated_m_s"] = v50_simulated
        metrics["v50_published_m_s"] = V50_PUBLISHED_M_S
        metrics["v50_published_status"] = V50_PUBLISHED_STATUS
        metrics["v50_relative_error"] = relative_error
        metrics["v50_fit_status"] = fit_status
        metrics["v50_sweep_shots"] = sweep_summary["shots"]
        if V50_PUBLISHED_STATUS == "VERIFIED":
            results["verdict"] = verdict
        else:
            results["verdict"] = "INCONCLUSIVE"
            metrics["inconclusive_reason"] = "V50_published is UNVERIFIED; simulated V50 reported but cannot be gated against an authoritative published number"
        RESULTS_JSON.write_text(json.dumps(results, indent=2))
        print(f"  results.json updated at {RESULTS_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
