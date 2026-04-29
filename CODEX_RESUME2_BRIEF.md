# Codex resume brief #2 — FEA_AP-PLY (post-rotation)

You run autonomously on `comech-2422` in tmux session `codex-fea`. Stage 07 hit a real numerical FAIL because of an incorrect orientation convention (`Ip=1, Iorth=1, skew_ID=1, Phi=0` does not consume the skew frame). The empirical rotation probe at **`references/openradioss_orientation_convention.md`** (commit `b9ec4b0`) discovered three working recipes that match the analytic transformation across angles {0, 30, 45, 60, 90}° within **0.24%**. Apply those recipes and redo stages 07-09, then continue 10, 12, 13-16.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. No AI provider mentions.

## Required reading

1. `README.md`, `plan/master_plan.md`.
2. **`references/openradioss_orientation_convention.md`** — the verified orientation recipes. Bind every composite material card change to a row of this doc.
3. `references/openradioss_law_compatibility_matrix.md` — confirms `/MAT/LAW12` + `/PROP/TYPE6` parses; do NOT use `/PROP/TYPE14` for composite LAWs.
4. `tests/stage_07_UD_tow_D3039/blocker.md` — the FAIL diagnosis to confirm the orientation cause.
5. `RESULTS_SUMMARY.md` and per-stage `tests/stage_NN_<name>/results/results.json` — current state.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH
```

Use `-nt 16`.

## The substitution + orientation rules (binding)

For composite materials in stages 06-16, use:

| Layer | Setting |
|---|---|
| Material | `/MAT/LAW12` (or `/MAT/LAW14`) |
| Property | `/PROP/TYPE6` (`SOL_ORTH`) — solid HEXA8, NOT TYPE14 |
| Uniform off-axis (single-angle coupon) | `Ip=3, Iorth=0, Phi=θ` on `/PROP/TYPE6` |
| Per-element rotation (laminate, woven) | `/INIBRI/ORTHO` per brick + neutral `/PROP/TYPE6` (`Ip=1, Iorth=0, Phi=0`) |
| Vector alternative | `/PROP/TYPE6` with `Vx,Vy,Vz=(cos θ, sin θ, 0)`, `Ip=13, Iorth=0` |
| /FAIL cards | unchanged (HASHIN, PUCK, TSAIWU all parse on LAW12+TYPE6 per matrix) |

**Strain measurement caveat** (also from the convention doc): for rotated TYPE6 solids, **VTK cell `Stra[0]` is NOT a reliable global axial engineering strain**. Compute strain from displacement-BC engineering values instead.

## Per-stage workflow

For each stage in {07, 08, 09, 10, 12, 13, 14, 15} (skip 11 and 16 — blocked by Kok geometry port):

1. `git pull --rebase`.
2. Edit `tests/stage_NN_<name>/runner.py` deck templates per the rule above:
   - Replace any `/PROP/TYPE14 + LAW{12,14,25}` with `/PROP/TYPE6 + LAW12`.
   - For uniform off-axis tests (07, single-angle parts of 08), set `Ip=3, Iorth=0, Phi=θ` on TYPE6.
   - For laminate stacks (09), set per-ply `/PROP/TYPE6` with `Ip=3, Iorth=0, Phi=θ_ply`.
   - For impact stages with rotated tows (13-15), use `/INIBRI/ORTHO` per element.
3. Run end-to-end (starter + engine + anim_to_vtk + verdict).
4. PASS → commit `stage NN <name>: PASS, post-rotation fix`, push.
5. FAIL (real numerical mismatch) → commit FAIL + blocker.md, **stop**.
6. INCONCLUSIVE (toolchain-blocked) → commit + blocker.md, **proceed**.
7. Update `RESULTS_SUMMARY.md` after each stage.

## Stage-specific notes

- **Stage 07**: redo with `Ip=3, Iorth=0, Phi=45` for the off-axis case. Expected PASS based on the convention doc's verified 0.24% match at 45°.
- **Stage 08** (ply rotation sweep at 0,15,30,45,60,75,90): use `Ip=3, Iorth=0, Phi=θ` for each angle. Expected PASS.
- **Stage 09** (cross-ply `[0/90]_s` and quasi-iso `[0/+45/-45/90]_s`): use per-ply `/PROP/TYPE6` with the right Phi for each ply. Expected PASS against CLT.
- **Stage 10** (UD mesoscale): isotropic fiber + matrix at this scale; no orientation issue. Run as-is.
- **Stage 11** (PW mesoscale): **requires Kok geometry port (`src/kok_geom/` empty)**. Mark INCONCLUSIVE with a clear blocker.md citing `plan/kok_port_plan.md` as the next major work item, proceed to 12.
- **Stage 12** (DCB/ENF cohesive): independent of orientation; should run with /INTER/TYPE2 or /MAT/LAW117 fallback as the previous brief specified.
- **Stage 13** (LVI): `/INIBRI/ORTHO` for rotated layups + /FAIL/HASHIN with erosion.
- **Stage 14** (CAI): chained restart from stage 13.
- **Stage 15** (flat ballistic): same /INIBRI/ORTHO recipe.
- **Stage 16** (PW panel ballistic): **requires Kok port**. Mark INCONCLUSIVE with same note as stage 11 and exit.

## Constraints (hard, unchanged)

- Solid HEXA8 / TETRA10 only (TYPE6 is solid orthotropic).
- No GUI, no closed source, no matplotlib.
- Brand colours per `~/.claude/CLAUDE.md`.
- Author = `j-vaught <jvaught@sc.edu>`.

## When to stop

- All redoable stages processed → write final `RESULTS_SUMMARY.md`, commit, push, exit.
- Genuine numerical FAIL on any stage → blocker.md + stop.
- Token / time budget → push everything and exit gracefully.

Begin with stage 07.
