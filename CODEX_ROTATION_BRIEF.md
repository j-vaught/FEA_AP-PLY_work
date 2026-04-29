# Codex rotation-correctness brief — FEA_AP-PLY

You run autonomously on `comech-2422` in tmux session `codex-rotation`. The previous codex-fea hit a real numerical FAIL on stage 07: `LAW12 + /PROP/TYPE6` parses and runs, but the off-axis 45° rotation does not produce the analytic E_x. Four orientation conventions were tried (Phi, /SKEW, Vx/Vy/Vz vector, ip/iorth flags); none matched. Your job is **a definitive empirical answer on how to apply per-element fiber orientation correctly** in this OpenRadioss build, OR a documented impossibility statement if none of the available conventions can match the analytic transformation.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. No AI provider mentions.

## Required reading

1. `references/openradioss_law_compatibility_matrix.md` — confirms LAW12 + TYPE6 parses; doesn't address rotation.
2. `tests/stage_07_UD_tow_D3039/blocker.md` — the failure description.
3. `tests/stage_07_UD_tow_D3039/runner.py` — see what the previous codex tried.
4. The OpenRadioss source CFGs at `/mnt/storage/j-vaught/openradioss/OpenRadioss-src/hm_cfg_files/config/CFG/RADIOSS/MAT/LAW012*.cfg` and `LAW014*.cfg` — these are the **canonical input grammar** for those LAWs. Read them. They tell you exactly which keyword fields exist and what each does.
5. The `/PROP/TYPE6` (`SOL_ORTH`) CFG at `hm_cfg_files/config/CFG/RADIOSS/PROP/TYPE006*.cfg` — same.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
export OR_SRC=/mnt/storage/j-vaught/openradioss/OpenRadioss-src
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH
```

## Methodology

### Phase A — Doc and example reconnaissance

1. **Read the LAW12 CFG and LAW14 CFG** files identified above. Tabulate every keyword field. Distinguish:
   - Material constants (E1, E2, ν12, G12, etc.)
   - Strength values
   - Orientation slots (skew_ID, Phi, Vx/Vy/Vz, Px/Py/Pz, Iorth, Ip)
   - Numerical flags

2. **Read the /PROP/TYPE6 CFG.** Note especially the `Skew_ID`, `Phi`, vector fields, and `Iorth`/`Ip` columns. Confirm whether the rotation is set in the LAW or in the PROP (or both).

3. **Search the OpenRadioss `qa-tests`** under `/mnt/storage/j-vaught/openradioss/OpenRadioss-src/qa-tests/` for any deck using LAW12 or LAW14 with non-trivial orientation (anything that grep matches for `LAW12` or `LAW14` plus `Phi` ≠ 0 or `/SKEW` ≠ default). Cite up to 5 examples. For each, dump the LAW12/14 + /PROP/TYPE6 + /SKEW block to confirm the canonical syntax.

Save these findings to `references/openradioss_orientation_recon.md` with cited examples and a "documented input grammar" table.

### Phase B — Probe sweep

Build a single probe driver that, for each test combination:

1. Defines the IM7/8552 LAW12 card (use the same elastic constants as `references/material_cards/im7_8552.json`):
   - E1 = 171.4 GPa, E2 = E3 = 9.08 GPa
   - ν12 = ν13 = 0.32, ν23 = 0.50
   - G12 = G13 = 5.29 GPa, G23 = 3.08 GPa
2. Defines a 1-element HEXA8 cube (or 4×4×1 small block for averaging) on /PROP/TYPE6.
3. Applies the orientation per the candidate convention.
4. Applies a uniaxial extensional displacement BC.
5. Runs `starter -check` first (must pass), then the engine.
6. Extracts the area-averaged σ and ε from the VTK; computes the apparent E_x.
7. Compares to the analytic transformation:
   ```
   1/E_x(θ) = c⁴/E1 + s⁴/E2 + (1/G12 − 2ν12/E1) c² s²
   ```
   Reference values for IM7/8552:
   - θ=0°  → E_x = E1 = 171.4 GPa
   - θ=30° → E_x ≈ 22.94 GPa
   - θ=45° → E_x ≈ 13.28 GPa
   - θ=60° → E_x ≈ 10.07 GPa
   - θ=90° → E_x = E2 = 9.08 GPa

### Combinations to sweep

Cross product of:

- **LAW**: {LAW12, LAW14}
- **Orientation mechanism**: 
  - M1: `Phi` field on `/PROP/TYPE6` (in-plane angle)
  - M2: `/SKEW/FIX` referenced via `Skew_ID` on `/PROP/TYPE6`, fiber along skew x-axis
  - M3: `Vx, Vy, Vz` reference vector on `/PROP/TYPE6` (orientation vector for ply 1-direction)
  - M4: `Vx, Vy, Vz` + `Px, Py, Pz` (full triad, second vector for transverse direction)
  - M5: per-element orientation via `/INIBRI/ORTHO` if available (look up — possibly the right answer)
  - M6: any other mechanism revealed by the CFG read in Phase A
- **Iorth flag**: {0, 1}
- **Ip flag**: {1, 3}
- **Test angles**: {0°, 30°, 45°, 60°, 90°}

For each (LAW × Mechanism × Iorth × Ip × angle) cell, record:
- Starter exit code
- Engine exit code  
- Recovered E_x (GPa)
- Analytic E_x (GPa)
- Relative error
- Verdict: PASS (≤ 2% on at least 3 angles), FAIL (numerics off), BLOCKED (parsing/engine error)

### Phase C — Output

Write `references/openradioss_orientation_convention.md` with:

1. **Recipe section**: the recommended (LAW, mechanism, flags) for arbitrary fiber rotation, citing the row of the matrix that supports it.
2. **Example deck snippet**: a copy-pastable `/MAT/LAW12 + /PROP/TYPE6 + /SKEW + ...` block that reproduces analytic E_x at θ ∈ {0, 45, 90} within 2%.
3. **Comparison table**: full results from Phase B.
4. **Honest impossibility statement** if no combination works: list every variant tried, the lowest off-axis error achieved, and recommend a path forward (e.g. switch to LAW28, accept off-axis bias for stage 16, or rebuild OpenRadioss with a different config).

Commit: `rotation: empirical /MAT/LAW × /PROP/TYPE6 orientation convention probe`. Push.

## Constraints (hard)

- Solid HEXA8 only.
- No GUI, no closed-source.
- Cite primary sources (CFG files, qa-tests).
- Author = `j-vaught <jvaught@sc.edu>`.
- Do not modify the existing stage runners. This is a focused side investigation.

## When to stop

- Found a working convention (≤ 2% on 3+ angles) — commit recipe + matrix + exit.
- Exhausted all reasonable variants without success — commit impossibility statement + recommendation + exit.
- Token / time budget — commit partial, push, exit.

Begin with Phase A (doc and example recon).
