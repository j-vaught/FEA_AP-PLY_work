# Codex Kok geometry M5 brief — physical undulation + interlocked tapes

You run autonomously on `comech-2422` in tmux session `codex-kok`. The previous M1-M4 port produced a *runnable* AP-PLY mesh, but stage 11 FAILed because the geometry is flat-tapes-with-resin (37% tow / 63% resin, no through-thickness interlock) instead of the interlocked AP-PLY architecture. Fix the geometry until **stage 11 PASSes the Kok 2022 stiffness gate within 10%**. **Iterate as long as it takes** — do NOT stop on the first failed attempt.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. No AI provider mentions.

## Required reading (clean-room, do NOT read Kok's source)

1. `tests/stage_11_PW_mesoscale_direct/blocker.md` — diagnoses the three gaps.
2. `tests/stage_11_PW_mesoscale_direct/results/results.json` and `effective_moduli.csv` — current numbers.
3. `plan/kok_port_plan.md` §2 (algorithmic decomposition) and §5 (tile-and-fragment scaling).
4. **Nagelsmit 2013 PhD thesis** Ch. 2 §§2.2 (Concept), 2.3 (Pattern Details). Cite section + figure for every algorithmic decision.
5. **Kok 2022 *Composites Part A***, §3.1 Microscale model. The undulation ratio χ = t/Lu, the over/under interlace pattern, the rotated through-thickness tow axes.
6. `references/openradioss_orientation_convention.md` — `/INIBRI/ORTHO` recipe with rotated axes per element, including the through-thickness component from undulation.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH
```

**Use `-nt 32` for solver runs** (the box has 112 cores; previous runs underused it). For the V50 sweep at stage 11 the runs are cheap (~minutes), but the solver SMP scaling is real on this size of problem.

## What to fix in `src/kok_geom/`

Three concrete geometry bugs, in priority order:

### G1 — Physical undulation solids

In each ply, every tape must follow the **3D undulation curve** specified in Kok 2022 §3.1 Figure 4:
- The neutral axis of a tape ramps up and over the perpendicular tape it's crossing, at undulation ratio **χ = t/L_u ≈ 0.09**, where `t` is the cured ply thickness and `L_u` is the undulation length.
- Build each tape as a swept solid (e.g. `gmsh.model.occ.addThruSections` along a Bspline path with the correct ramp-up / ramp-down regions, or as a Boolean union of straight + ramped segments).
- Each tape's local material axis-1 direction must follow the **tangent of the undulation curve**, NOT the global x-axis. That includes the small through-thickness component during the up/down segments.

### G2 — Interlocked over-and-under pattern

For each pair of crossing tapes (different fiber angles), one goes over and one goes under. The pattern repeats per `tape_spacing s ∈ {1, 2, 3}` per Nagelsmit 2013 §2.3:
- s=1: every crossing alternates one-over-one (the dense interlace).
- s=2: every other crossing.
- s=3: every third.
- Default for stage 11: s=1, four-ply `[0/+45/-45/90]`.

Use the **tile-and-fragment** approach from `plan/kok_port_plan.md` §5: build one unit cell of the over/under pattern as a single OCC compound, then array-replicate. Avoid one giant Boolean cascade.

### G3 — Resin only in actual voids

After all interlocked tapes are placed, the remaining void volume should be **≤ 15-20%**, not 63%. Resin pockets exist only:
- Between adjacent parallel tapes within a ply (small gap if tapes don't fully cover the lane).
- At locations where the ramped tape leaves a triangular pocket against the perpendicular tape it's crossing.

Boolean-fragment the resin domain against the union of all tapes; do not fill all unoccupied space.

## Acceptance criteria (binding, don't relax)

After each geometry change, regenerate the stage 11 mesh and run all three loadings (axial-x, transverse-y, in-plane shear). Compare to Kok 2022 targets at the **single tested architecture** (one-over-one quasi-iso `[0/+45/-45/90]`):

```
E_x  target 53.3 GPa,  PASS within 10% (47.97 - 58.63 GPa)
E_y  target 53.3 GPa,  PASS within 10%
G_xy target 20.5 GPa,  PASS within 10% (18.45 - 22.55 GPa)
```

PASS = **all three** within 10%. INCONCLUSIVE = within 25% (we accept material delta between VTC401 vs IM7/8552). FAIL = beyond 25%.

If you reach **PASS**: flip `tests/stage_11_PW_mesoscale_direct/results/results.json` verdict to PASS, commit `stage 11 PW_mesoscale_direct: PASS, post-Kok-M5`, push.

If you reach **INCONCLUSIVE** but stable: commit `stage 11 PW_mesoscale_direct: INCONCLUSIVE, material delta`, document the VTC401 vs IM7/8552 stiffness arithmetic in `blocker.md`.

If you keep getting >25% off after 5+ geometry iterations: write `blocker.md` with a fully diagnostic "what we tried, what we measured, what's still wrong" report and exit. Do not invent geometric tweaks beyond what the papers actually describe.

## Iteration loop

1. Fix one of G1/G2/G3 in `src/kok_geom/`.
2. Add or update unit tests verifying the fix (volume fractions, undulation curvature, interlace count).
3. Regenerate the stage 11 mesh via `python -m kok_geom --config geometry/kok_stage11_config.json --out tests/stage_11_PW_mesoscale_direct/runs/panel.{msh,inp}` (or whatever the canonical CLI invocation is).
4. Run the three stage 11 loadings.
5. Read the moduli, compare to the targets, write the deltas to `effective_moduli.csv`.
6. Decide: PASS → commit and exit. Otherwise → iterate on the largest-delta cause and go back to (1).

Commit a snapshot **after every iteration** even if it FAILs, with a clear message like `kok M5 iter N: undulation solids, E_x 31.2 GPa (-41%)`. We want to see the geometry-→-stiffness traceability in git history.

## Stage 16 follow-up (only after stage 11 PASS)

If stage 11 PASSes within budget, immediately re-run stage 16 with the corrected geometry on the **reduced 200×200 mm representative section**, `-nt 32`, full V50 sweep at multipliers (0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20)·V50_published. The previous "compute budget exceeded" estimate was at `-nt 16`; with `-nt 32` shots should run in ~2 h each, full sweep overnight.

Commit `stage 16 PW_panel_ballistic: <verdict>, post-Kok-M5` when done.

## Constraints (hard, unchanged)

- Solid HEXA8 / TETRA10 only.
- Pure Python; gmsh.model.occ.* for geometry; MIT-licensed; no Kok source code.
- LAW12 + TYPE6 + /INIBRI/ORTHO recipe — no other constitutive choices.
- No GUI, no closed-source.
- Author = `j-vaught <jvaught@sc.edu>`.

## When to stop

- Stage 11 PASS or INCONCLUSIVE (material delta) committed → if compute permits, run stage 16 + commit → exit.
- 5 iterations without a stiffness improvement of at least 25% per iteration → document in blocker, exit.
- Token / time budget — push everything in flight, exit gracefully.

Begin with G1 (physical undulation solids).
