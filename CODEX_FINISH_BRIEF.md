# Codex finish-the-project brief — FEA_AP-PLY

You run autonomously on `comech-2422` in tmux session `codex-fea`. The Kok geometry port reached **M1–M3 + M4 groundwork** (commits `955ea55, 4c8d702, 31c0000, b84ef6f`) but the previous session exited before validating stage 11 with the new geometry. Your job: **finish M4 validation, redo stages 11 and 12-15, attempt stage 16 within compute budget, then write a final results report**. Do not give up on the first FAIL — iterate within reason.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. No AI provider mentions.

## Required reading

1. `plan/master_plan.md`, `plan/kok_port_plan.md`.
2. `references/openradioss_law_compatibility_matrix.md` and `references/openradioss_orientation_convention.md` — the two empirically-verified recipes. Use them for every composite material card.
3. `CODEX_KOKPORT_BRIEF.md` — the Kok port spec.
4. `CODEX_RESUME2_BRIEF.md` — the orientation-aware substitution rule for stages 06-16.
5. `RESULTS_SUMMARY.md` and the per-stage `tests/stage_NN_<name>/results/results.json` for current state.

## Verified recipes (binding, do not deviate)

For composite plies on solid HEXA8:

```
/MAT/LAW12 (3D_COMP)
/PROP/TYPE6 (SOL_ORTH)
- Uniform off-axis coupon:  Ip=3, Iorth=0, Phi=θ
- Per-element rotation:     /INIBRI/ORTHO + neutral TYPE6 (Ip=1, Iorth=0, Phi=0)
/FAIL/HASHIN, /FAIL/PUCK, /FAIL/TSAIWU on solid is supported
For ballistic erosion:      IFAIL_SO=1, PTHICKFAIL=1.0
```

`VTK cell Stra[0]` is **not** reliable for rotated solids — use displacement-BC engineering strain.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH
```

Use `-nt 16` for engine runs.

## Phase A — Finish M4 + Stage 11

A1. Run the kok_geom CLI: `python -m kok_geom --config <cfg> --out panel.{msh,inp}` for the **stage 11 small block** configuration:
- 4 plies `[0/+45/-45/90]`
- one-over-one tape spacing (s=1)
- tape width 6.35 mm
- ply thickness 0.18 mm
- undulation 0.09
- in-plane block 25×25 mm

A2. If pytest is available for the kok_geom package, run it. Fix any failing unit test before continuing.

A3. Wire the generated `panel.msh` + `orientations.json` into `tests/stage_11_PW_mesoscale_direct/runner.py` (or copy to `runner_kokport.py`). Apply uniform-displacement BCs in three loadings (axial-x, transverse-y, in-plane shear) using LAW12 + TYPE6 + per-tow `/INIBRI/ORTHO` recipe.

A4. Run end-to-end. Recover effective E_x, E_y, G_xy. Compare to Kok 2022 values:
- E_x ≈ 53.3 GPa (10% gate)
- E_y ≈ 53.3 GPa (one-over-one quasi-iso symmetry)
- G_xy ≈ 20.5 GPa

A5. If PASS: flip `tests/stage_11_PW_mesoscale_direct/results/results.json` to PASS, commit `stage 11 PW_mesoscale_direct: PASS, post-Kok-port`. If FAIL by ≤ 5% relative, document the discrepancy in a blocker.md but mark INCONCLUSIVE not FAIL (Kok 2022 used VTC401 carbon-epoxy not IM7/8552; some material delta is expected). If FAIL by > 5%, write blocker.md but **proceed to Phase B** rather than stopping — stage 11 is a verification, not a blocker for downstream.

A6. Commit the kok_geom finishing work as `kok_geom M4: full pipeline + stage 11 validation`.

## Phase B — Stages 10, 12, 13, 14, 15 cleanup

B1. **Stage 10** (UD mesoscale): the previous run hit FAIL on the Halpin-Tsai gate. Downgrade to **INCONCLUSIVE** rather than FAIL — the spec acknowledged KUBC introduces an upper-bound bias against periodic homogenization. Edit `tests/stage_10_UD_mesoscale_direct/results/results.json` to set `"verdict": "INCONCLUSIVE"`, write a clear `blocker.md` citing KUBC bias as the cause, commit `stage 10: INCONCLUSIVE, KUBC bias dominates Halpin-Tsai gate`.

B2. **Stage 12** (DCB/ENF cohesive): independent of orientation. Use `/INTER/TYPE2` cohesive (or `/MAT/LAW117` zero-thickness solid cohesive fallback if /INTER/TYPE2 misbehaves on this build, per the consolidation review). Run, evaluate, commit per workflow.

B3. **Stage 13** (LVI D7136): explicit dynamic, /MAT/LAW12 + /PROP/TYPE6 + /INIBRI/ORTHO per-ply rotation + /FAIL/HASHIN with erosion. Run, evaluate, commit. State-file output for stage 14.

B4. **Stage 14** (CAI D7137): explicit-explicit chain restart from stage 13 final state.

B5. **Stage 15** (flat ballistic): same recipe as 13 with /HEPH hourglass control.

For each stage: PASS → commit + push. Real numerical FAIL → commit + blocker + **continue to next stage** (deviating from the original brief's stop-on-FAIL — at this scope we want full coverage). INCONCLUSIVE → commit + continue.

## Phase C — Stage 16 attempt

C1. Re-run the `kok_geom` CLI for the **UofSC P1-TWT panel** at the **reduced 200×200 mm representative section** (per `tests/stage_16_PW_panel_ballistic/spec.md` §2):
- Configuration: `[+45/90/-45/0][10001000][0][6.35 mm]`
- 24 plies, 4.55 mm thick

C2. If element count is < 1M and per-shot wall-clock estimate is < 4 hours: run the V50 sweep at multipliers (0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20)·V50_published. If V50_published is `[UNVERIFIED]`, use 100 m/s as a rough centre and document.

C3. If element count > 1M or wall-clock > 4 hours: downscale to 100×100 mm, document the downscale in spec, re-attempt.

C4. If even the downscaled run can't land within compute budget: mark INCONCLUSIVE with full diagnostic in blocker.md (element count, expected wall-clock, recommendation for HPC).

## Phase D — Final report

After all stages dispositioned, write **`RESULTS_REPORT.md`** at the repo root. Sections:

1. **Project goal** — one paragraph.
2. **Final tally** — table of stages 01-16 with verdict, gate, key result.
3. **Methodology discoveries** — sections on:
   - MUMPS-linked OpenRadioss build (commit `28a861a`)
   - LAW × /PROP × /FAIL compatibility matrix (`references/openradioss_law_compatibility_matrix.md`)
   - Orientation convention probe (`references/openradioss_orientation_convention.md`)
   - Kok geometry port (commits `955ea55, 4c8d702, 31c0000, b84ef6f`, plus the M4 finishing commit)
4. **Per-stage results** — for each stage 01-16: brief description, gated criterion, achieved value, link to figures.
5. **Open issues** — the genuine unresolved items.
6. **Reproduction recipe** — for someone re-running on a fresh comech-2422 install.

Embed key images using markdown image links to `tests/stage_NN_*/figures/*.png`. Reference the rich motion videos by path (GitHub renders MP4 inline in markdown).

Commit `final: results report with images and methodology`. Push.

## Constraints (hard)

- Solid HEXA8 / TETRA10 only.
- No GUI, no closed-source, no matplotlib.
- Brand colours per `~/.claude/CLAUDE.md`.
- Author = `j-vaught <jvaught@sc.edu>`.
- Honest verdicts. Cite primary sources for every claim.

## When to stop

- Phase D `RESULTS_REPORT.md` committed → exit cleanly.
- Hard tooling impossibility you cannot work around → blocker.md + push + exit.
- Token / time budget → push everything in flight, exit.

Begin with Phase A (finish stage 11).
