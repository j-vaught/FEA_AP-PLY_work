# Codex LAW × /PROP compatibility matrix brief — FEA_AP-PLY

You run autonomously on `comech-2422` in tmux session `codex-lawmatrix`. The previous codex-fea hit a cascade of `LAW25 TYPE14 starter blocked` INCONCLUSIVE verdicts on stages 06-09. The root cause is almost certainly that `/MAT/LAW25` is **CRASURV / COMPSH** (composite *shell* card by design) and is not accepted on `/PROP/TYPE14` solid bricks. The fix is a different `/MAT/LAW*` for solid composites — but instead of guessing, **discover the compatibility matrix empirically**.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. Do not mention any AI provider.

## Goal (one shot)

Produce `references/openradioss_law_compatibility_matrix.md` — a definitive markdown matrix telling future codex sessions exactly which LAW × /PROP × /FAIL combinations parse cleanly in this OpenRadioss build. Then exit.

## Required reading (short)

1. `README.md`, `plan/master_plan.md` — for the constraint context (solid HEXA only).
2. `tests/stage_06_composite_failure_criteria/blocker.md`, `tests/stage_07_UD_tow_D3039/blocker.md` — see what error codex-fea hit.
3. `references/openradioss_endtoend_audit.md` — for the LAW recommendations cited there.
4. `/mnt/storage/j-vaught/openradioss/OpenRadioss-src/hm_cfg_files/config/CFG/` — *the source-of-truth keyword definitions for every /MAT/LAW*. Look for `LAW000.cfg`, `LAW001.cfg`, etc. (or whatever the layout is) to confirm which LAWs exist in this build and their canonical solid-vs-shell applicability.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH
```

Use the **MUMPS-linked** binaries already in `$OR/exec/` (commit `28a861a`).

## Methodology

### A. Enumerate candidate LAWs

From the OpenRadioss CFG directory `hm_cfg_files/config/CFG/`, list every `/MAT/LAW*.cfg` (or equivalent) that mentions in its name or header any of: composite, orthotropic, anisotropic, transversely isotropic, tsai, hashin, puck, ply, fiber, fibre. Likely candidates (verify each exists in the build):

- LAW1 (linear elastic isotropic) — control / sanity
- LAW2 (PLAS_JOHNS, Johnson-Cook)
- LAW12 (anisotropic Tsai-Wu solid?)
- LAW14 (composite solid orthotropic)
- LAW15 (Drucker-Prager / cohesive — check)
- LAW17 (composite solid)
- LAW25 (COMPSH composite shell — known failure on TYPE14, included for matrix completeness)
- LAW28 (composite solid)
- LAW42 (Ogden hyperelastic)
- LAW53 (Tsai-Wu solid)
- LAW55 (Hashin)
- LAW58 (anisotropic)
- LAW62 (composite solid)
- any other LAW*.cfg present that fits the criteria

### B. Build a minimal probe deck

For each candidate LAW, build the **simplest possible starter deck**:

- 1 HEXA8 solid brick (TYPE14) — for the solid track.
- 1 SHELL element (TYPE1 or TYPE6) — for the shell-control track (some LAWs only work on shells; we need the contrast).
- /MAT/LAW{candidate} with the simplest representative parameter set (look up `/MAT/LAW{N}` in `hm_cfg_files/config/CFG/...` for the keyword stub).
- /PROP/TYPE14 referencing the LAW for the solid; /PROP/TYPE1 or TYPE6 for the shell.
- A simple /BCS clamp + /CLOAD or /IMPDISP load.
- /END.

Generate, then run `$OR/exec/starter_linux64_gf -i deck_<LAWN>_solid.rad -check` and capture exit code + ERROR ID.

Repeat with `/FAIL/HASHIN`, `/FAIL/PUCK`, `/FAIL/TSAIWU`, `/FAIL/CHANG`, `/FAIL/USER1`, `/FAIL/JOHNSON_COOK` attached to the LAW. Run starter -check each time, record outcomes.

### C. Build the matrix

Schema for the markdown matrix at `references/openradioss_law_compatibility_matrix.md`:

```
# OpenRadioss LAW × /PROP × /FAIL compatibility matrix

## Verified on
- OpenRadioss build: <starter -version output>
- MUMPS-linked: yes
- Date: 2026-04-29

## Matrix

| LAW       | TYPE14 (solid) | TYPE1/6 (shell) | /FAIL/HASHIN | /FAIL/PUCK | /FAIL/TSAIWU | Other /FAIL | Notes |
|-----------|---------------:|----------------:|-------------:|-----------:|-------------:|-------------|-------|
| LAW1      | OK            | OK             | n/a          | n/a        | n/a          | n/a         | isotropic baseline |
| LAW2      | OK            | OK             | -            | -          | -            | JOHNSON_COOK | iso plastic, no composite damage |
| LAW12     | <verdict>     | <verdict>      | <verdict>    | <verdict>  | <verdict>    | <verdict>   | <error id if blocked> |
| ...       | ...           | ...            | ...          | ...        | ...          | ...         | ... |
```

Verdicts: `OK` (starter -check returns 0), `BLOCKED:<error_id>` (starter rejects with specific error ID), `n/a` (the FAIL/LAW combo doesn't make physical sense).

## Recommendations section

After the matrix, write a **"Recommended substitutions"** section. For each downstream need:

- "Linear orthotropic solid (no damage)" → which LAW
- "Solid composite with Hashin failure" → which LAW + which /FAIL
- "Solid composite with Tsai-Wu failure" → which LAW + which /FAIL
- "Solid composite with Puck failure" → which LAW + which /FAIL
- "Solid composite with element erosion (for ballistic)" → which LAW + flags

Each recommendation must be **directly verified** in the matrix above (cite the row) — do not extrapolate.

## Output

1. `references/openradioss_law_compatibility_matrix.md` — the matrix + recommendations.
2. `references/lawmatrix_probe_decks/` — the small probe decks (`deck_LAW12_solid.rad`, etc.) and their `*.out` listing files for traceability. Commit them, but they're tiny.
3. Single commit: `lawmatrix: empirical /MAT/LAW × /PROP × /FAIL compatibility for solid composites`.
4. Push.
5. Exit cleanly.

## Constraints (hard)

- Solid elements only for the actual project (this matrix exercise IS allowed to test shells for the contrast — that is the whole point).
- No GUI, no matplotlib.
- No closed-source software.
- Author = `j-vaught <jvaught@sc.edu>`.

## When to stop

- Matrix written, recommendations section written, committed, pushed — exit.
- An OpenRadioss bug or environment issue you can't work around — write `references/lawmatrix_blocker.md`, push, exit.
- Time budget — push partial matrix and exit gracefully.

Begin.
