# Codex resume brief — FEA_AP-PLY (post-LAW-matrix)

You run autonomously on `comech-2422` in tmux session `codex-fea`. The previous run hit a `LAW25 + /PROP/TYPE14` cascade across stages 06-09 (all INCONCLUSIVE). The empirical compatibility matrix at `references/openradioss_law_compatibility_matrix.md` has revealed the real fix: **swap `/PROP/TYPE14` → `/PROP/TYPE6` (`SOL_ORTH`) for composite materials**. `TYPE6` is an orthotropic *solid* property, NOT a shell — the user's "all solid" constraint is preserved.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. Do not mention any AI provider.

## Required reading (in order)

1. `README.md`, `plan/master_plan.md`, `CODEX_BRIEF.md` (the original autonomous brief — workflow rules still apply).
2. **`references/openradioss_law_compatibility_matrix.md`** — the verified substitution table. Bind every composite-material change in the runners to a row of this matrix. Do NOT use unverified combinations.
3. `RESULTS_SUMMARY.md` and the existing per-stage `tests/stage_NN_<name>/results/results.json` — current state.
4. The four stage-06-09 `blocker.md` files — to confirm the root cause matches the matrix's Error 3047 finding.

## What changed since last run

- Stages 01, 03, 04, **05 PASS** (don't redo).
- Stage 02 INCONCLUSIVE (Newton-Raphson divergence at α=3 — *not* a LAW/PROP issue; leave the existing INCONCLUSIVE commit, do not retry stage 02 here).
- Stages 06, 07, 08, 09 INCONCLUSIVE — **redo with the matrix substitutions below**.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss     # MUMPS-linked binaries
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH
```

Use `-nt 16` for solver runs (the box has 112 cores; previous runs at `-nt 8/16` ran cleanly).

## The substitution rule (binding)

For every stage with composite materials (06, 07, 08, 09, 11, 13, 14, 15, 16):

| Replace | With |
|---|---|
| `/PROP/TYPE14` (referencing a composite LAW) | `/PROP/TYPE6` (`SOL_ORTH`) |
| `/MAT/LAW25` (COMPSH — shell only) | `/MAT/LAW12` (3D_COMP — solid composite) |
| /FAIL/HASHIN, /FAIL/PUCK, /FAIL/TSAIWU | unchanged (they parse on LAW12 + TYPE6 per matrix) |

The matrix verified that `LAW12 + TYPE6 + {HASHIN, PUCK, TSAIWU}` all parse cleanly via `starter -check`. Do NOT use `LAW25` or `TYPE14` for composite solids in this build. Isotropic stages (01, 03, 04) remain `LAW1/LAW2/LAW36 + TYPE14` — those are unchanged.

For ballistic erosion (stages 13-16), the matrix recommendation: `LAW12 + TYPE6 + /FAIL/HASHIN` with `IFAIL_SO=1`, `PTHICKFAIL=1.0`.

## Per-stage workflow

For each stage `NN ∈ {06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 16}`, in order:

1. `git pull --rebase`.
2. Read `tests/stage_NN_<name>/spec.md` for the test goal.
3. Edit `tests/stage_NN_<name>/runner.py` deck templates: replace `/PROP/TYPE14` → `/PROP/TYPE6`, replace `/MAT/LAW25` → `/MAT/LAW12`. Keep all `/FAIL/*` cards as specified.
4. Re-run end-to-end (starter + engine + anim_to_vtk + PyVista + verdict).
5. If verdict = PASS: commit `stage NN <name>: PASS, post-matrix fix`, push.
6. If verdict = FAIL (real numerical mismatch, not a parsing issue): commit `stage NN <name>: FAIL, <reason>`, write `blocker.md`, **stop** (per original brief).
7. If verdict = INCONCLUSIVE (toolchain still blocked despite matrix substitution): commit, write `blocker.md`, **proceed** to next stage.
8. Update `RESULTS_SUMMARY.md` after each stage.

## Stage-specific notes

- **Stages 06-09**: redo only with the matrix substitution. These were INCONCLUSIVE solely due to TYPE14; should now PASS.
- **Stage 10** (UD mesoscale): **isotropic** fiber + matrix at this scale. Uses LAW1, no composite cards. Should NOT need the substitution. Just run.
- **Stage 11** (PW mesoscale): **requires the Kok geometry port** (`src/kok_geom/` — currently empty). Per `plan/kok_port_plan.md`, the port is a 4-week milestoned effort. **You will not finish stage 11 in this session.** Mark INCONCLUSIVE with a note that the Kok port is the next major work item, write a clear blocker.md, proceed to stage 12.
- **Stage 12** (DCB/ENF): cohesive interface. May still hit the `/INTER/TYPE2` cohesive vs LAW117 fallback documented in consolidation review. The matrix doesn't directly address cohesive cards — stage 12's runner already has the auto-fallback.
- **Stages 13-14** (LVI + CAI): explicit dynamic + state-restart chain. Use the matrix substitution. Stage 14 chains off stage 13.
- **Stages 15-16** (ballistic): same substitution + element erosion flags. Stage 16 needs the Kok geometry; if Kok port isn't done, mark INCONCLUSIVE with the same note as stage 11 and exit.

## Constraints (hard, unchanged)

- Solid HEXA8 / TETRA10 only. (TYPE6 is still solid — orthotropic solid property.)
- No GUI, no closed-source, no matplotlib.
- Brand colours per `~/.claude/CLAUDE.md`.
- Author = `j-vaught <jvaught@sc.edu>`.

## When to stop

- All redoable stages (06-16) processed (PASS / FAIL+stop / INCONCLUSIVE+continue per the rules) → write final `RESULTS_SUMMARY.md`, commit, push, exit.
- Genuine numerical FAIL on any stage → blocker.md + stop.
- Token / time budget → push everything and exit gracefully.

Begin with stage 06.
