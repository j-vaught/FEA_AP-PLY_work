# Codex stage 16 V50 sweep brief — post-Kok-M5

You run autonomously on `comech-2422` in tmux session `codex-stage16`. Stage 11 PW_mesoscale_direct just landed PASS at commit `566c7be` with the M5 geometry overhaul (physical undulation solids, interlocked over-and-under, ≤7% resin). Three moduli are inside the 10% gate: `E_x 49.04 GPa (-7.99%)`, `E_y 53.29 GPa (-0.01%)`, `G_xy 22.53 GPa (+9.88%)`. Element count per unit volume is now **lower** than the pre-M5 stage 11 mesh (15,000 TETRA10 in 25 mm³ block vs 25,918 before — roughly 60%), because the resin volume collapsed from 63% to 7%.

Your job: take that corrected geometry pipeline and finish stage 16 — **the project's terminal validation gate** (`tests/stage_16_PW_panel_ballistic/spec.md` §1). Either land PASS / FAIL / INCONCLUSIVE on a real V50 sweep, or land a *measured* (not estimated) compute-budget verdict that we trust.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. No AI provider mentions.

## Required reading

1. `tests/stage_16_PW_panel_ballistic/spec.md` — full stage spec; §1 (terminal goal), §2 (P1-TWT laminate definition + reduced 200×200 mm section), §9 (V50 bracketing sweep), §6 (BC insensitivity check).
2. `tests/stage_16_PW_panel_ballistic/blocker.md` — prior blocker (pre-M5). Wall-clock estimate of 4.19 h/shot was an extrapolation from the *old* 25,918-element stage 11 mesh at `-nt 16`; both inputs to that estimate have changed.
3. `tests/stage_16_PW_panel_ballistic/runner.py` — existing self-contained runner with V50 bracketing, Recht-Ipson fit, BC-insensitivity check. Has dry-run mode.
4. `tests/stage_16_PW_panel_ballistic/geometry/kok_p1_twt_200_config.json` and `kok_p1_twt_100_config.json` — existing config jsons. May need to be regenerated against the new `kok_geom` API after M5.
5. `references/openradioss_law_compatibility_matrix.md` and `openradioss_orientation_convention.md` — recipes (verified): `/MAT/LAW12 + /PROP/TYPE6/SOL_ORTH + /INIBRI/ORTHO`, with `/FAIL/HASHIN` (`IFAIL_SO=1, PTHICKFAIL=1.0`) for erosion and `/INTER/TYPE2` cohesive between adjacent plies.
6. `CODEX_KOK_M5_BRIEF.md` — the brief that produced stage 11 PASS. Geometry decisions you can reuse verbatim.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH
```

**Use `-nt 32` for solver runs.** The 112-core box has been chronically underused. Previous wall-clock estimates were at `-nt 16`; expected ~1.8× speedup from the SMP scaling on this size of problem.

## Phase-gated execution

Do NOT just dive into a full V50 sweep blind. Run three phase-gates and only proceed to the next on a numerical green light. Commit between phases so we have traceable artifacts even if a later phase blocks.

### Phase A — Geometry meshing time-budget probe

Goal: prove the M5 geometry pipeline can mesh the 200×200 mm P1-TWT panel within a tolerable wall-clock (target ≤30 min on `-nt 32`).

1. Update `tests/stage_16_PW_panel_ballistic/geometry/kok_p1_twt_200_config.json` to the post-M5 schema if needed. Mirror the field names used by the working `kok_stage11_config.json`.
2. Invoke the CLI with explicit timing: `time python -m kok_geom --config <200_config.json> --out <out.{msh,inp}>`. Cap at 60 minutes wall-clock.
3. If 200 mm meshes inside 60 min: capture node/element count, commit `stage 16 phase A: 200 mm M5 mesh, <N> nodes / <M> TETRA10 in <T> min`.
4. If 200 mm exceeds 60 min: tile-and-fragment scaling from `plan/kok_port_plan.md` §5 should beat that — diagnose. If still infeasible, downscale to 150×150 mm (Saint-Venant buffer is conservative; 8× tape width = 51 mm minimum from impact center to clamped edge per spec §2.2). Document and proceed.

### Phase B — Single-shot wall-clock measurement

Goal: replace the 4.19 h/shot estimate with a measured number, on the actual M5 geometry, at `-nt 32`, with full erosion + contact stack.

1. Use the ASTM-style mid-velocity shot (`v_impact = V50_published_estimate ≈ 100 m/s` per spec §2 and `runner.py`). If V50_published is still `[UNVERIFIED]`, use 100 m/s as the bracket center; document the digitization gap separately.
2. Build the OpenRadioss deck from the Phase A mesh via `inp2rad` (or the project's verified deck-writer). Apply the LAW12 + TYPE6 + INIBRI/ORTHO + HASHIN + INTER/TYPE2 stack as in stage 11.
3. Run starter + engine with `-nt 32`. Use `/ANIM/DT 5.0e-6` and `/TFILE 1.0e-6` per spec; `/RUN 0.5e-3` simulation duration.
4. Record actual wall-clock. Commit `stage 16 phase B: single-shot measured <T> h on <N>-elem M5 mesh, -nt 32`.

**Decision gate:** if measured wall-clock ≤ 3 h/shot → proceed to Phase C (full sweep is overnight feasible). If 3 < T ≤ 6 h/shot → proceed but only run the **5-shot reduced sweep** (multipliers 0.85, 0.95, 1.00, 1.05, 1.15). If T > 6 h/shot → do NOT launch the sweep; document the measured infeasibility with the *real* numbers and write `INCONCLUSIVE: compute-bound, measured at -nt 32 post-Kok-M5`.

### Phase C — V50 sweep + Recht-Ipson fit

Goal: complete the actual V50 verification.

1. Run the V50 sweep at multipliers `(0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20) · V50_published` (or the reduced 5-shot if Phase B chose that path), all `-nt 32`.
2. Parse residual velocity from each T01; fit Recht-Ipson 1963 model; compute simulated V50.
3. Run the BC-insensitivity check from spec §6 (clamped vs absorbing) on the center-velocity shot only.
4. PASS / FAIL / INCONCLUSIVE per spec §1: 7% tolerance on V50 (loosen to 9% if BC-insensitivity check fails).
5. Update `tests/stage_16_PW_panel_ballistic/results/results.json` and commit `stage 16 PW_panel_ballistic: <verdict>, post-Kok-M5`.

## Acceptance criteria

- **PASS**: V50_simulated within 7% of V50_published (Vakili Rad 2020 P1-TWT) and BC-insensitivity check passes; or within 9% if BC check fails.
- **FAIL**: V50_simulated outside 9% in all BC variants.
- **INCONCLUSIVE — compute-bound**: Phase B measures wall-clock that makes Phase C infeasible; document with concrete `-nt 32`, post-M5 numbers (this is materially different from the prior INCONCLUSIVE).
- **INCONCLUSIVE — V50_published unverified**: if you cannot pin a published V50 from the abstract or thesis citations, document the digitization gap and report the simulated V50 with appropriate caveat.

## Constraints (hard)

- Solid HEXA8 / TETRA10 only.
- Pure Python; gmsh.model.occ.* for any new geometry; MIT-licensed; clean-room (no Kok source).
- LAW12 + TYPE6/SOL_ORTH + INIBRI/ORTHO + FAIL/HASHIN + INTER/TYPE2 — no other constitutive choices.
- Erosion flags `IFAIL_SO=1, PTHICKFAIL=1.0`.
- Hourglass formulation `Iform=24` on TYPE6/TYPE14; `E_hg/E_int ≤ 5%` hard gate.
- No GUI, no closed-source.
- Author = `j-vaught <jvaught@sc.edu>`.

## Iteration discipline

Commit a snapshot **after every phase**, even if it fails or downscales. Messages should look like `stage 16 phase A: <fact>` / `stage 16 phase B: <measurement>` / `stage 16: <verdict>`. We want the traceability of geometry → wall-clock → V50 in git history.

## When to stop

- Stage 16 PASS / FAIL committed → **project terminal goal reached** → exit clean.
- INCONCLUSIVE — compute-bound (with measured numbers, not estimates) → exit clean.
- INCONCLUSIVE — V50_published unverified → exit clean.
- Token / time budget — push everything in flight, exit gracefully.

Begin with Phase A.
