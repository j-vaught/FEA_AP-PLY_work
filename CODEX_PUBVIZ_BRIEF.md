# Codex publication-grade visualization brief — FEA_AP-PLY

You run autonomously on `comech-2422` in tmux session `codex-pubviz`. This is the **final visualization pass** — produce publication-grade figures and motion videos that exceed the existing first-pass output. The previous richviz already wrote the basic motion videos for stages 01-06; your job is to make them better and to add equivalents for stages 07-16 where solver output exists. After this you exit and the project wraps.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. No AI provider mentions.

## Required reading

1. `RESULTS_REPORT.md` (commit `aad2d94`) — the final deliverable. Each section's images should be cross-linked from the report.
2. `~/.claude/CLAUDE.md` — brand-colour palette, no rounded corners, no matplotlib.
3. `references/openradioss_orientation_convention.md` — strain measurement caveat (use displacement-BC engineering strain, not VTK Stra[0], for rotated solids).

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply  # mesalib already installed → OSMesa headless render works
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
```

## What "better than now" means concretely

Current figures (commits `213e6c3`, `6f53873`, `8a6f841`, `ebbfa41`, `8937061`, `331e12d`, `6519795`) have:
- ~5 MP4s with 1-15 frames (mostly looped from limited anim output).
- ~25 PNG stills, often single-field colormaps.
- Basic Typst plots from CSV.

Upgrade to:

### A. Multi-panel composite figures
For each PASS stage 01-09 with VTK output, produce **`stage_NN_composite.png`** (2400×1600, 4-panel):
- Top-left: deformed mesh coloured by displacement magnitude (mm).
- Top-right: same warp coloured by von Mises stress (MPa).
- Bottom-left: same warp coloured by equivalent strain (%) or fiber-direction stress (MPa) — pick whichever is most informative for the stage.
- Bottom-right: 2D Typst+CeTZ time-history (load-displacement, modulus-vs-angle, etc., using the existing `stage_NN_*_probe.typ` data).
Each panel has its own colourbar with brand colors and units. Stage label and verdict in a small footer band.

### B. Higher-quality motion videos
For each PASS stage with multi-frame anim output, produce **`stage_NN_motion_hd.mp4`** (1920×1080, 30 fps, 8-12 s, H.264):
- If the existing run only has 1-2 frames, **regenerate the deck with `/ANIM/DT` set to write 50+ frames** (mirror the richviz pattern). DO NOT touch already-passing verification artifacts; write to `tests/stage_NN_*/runs/<rootname>_pubviz/` so they're separate.
- Synchronized split-screen: left half displacement-coloured, right half von Mises-coloured, with a moving time cursor in the lower bar.
- Title card at the start (1 s static), end card (1 s) with verdict + key result.

### C. Cross-stage summary figure
**`figures/cross_stage_summary.typ`** + compiled `cross_stage_summary.pdf`:
- Single Typst+CeTZ figure that visualises the 16-stage final tally as a colored grid.
- Columns: stage number; rows: gate result vs criterion.
- Colour: green for PASS, amber for INCONCLUSIVE (toolchain-blocked or downgrade), red for FAIL, grey for not-attempted.
- Brand colours: PASS = `#65780B` (horseshoe), INCONCLUSIVE = `#A49137` (honeycomb), FAIL = `#73000A` (garnet), N/A = `#A2A2A2` (50% black).
- Bottom row: methodology chips for the three discoveries (MUMPS rebuild, LAW matrix, orientation convention) + Kok port progress (M1-M4).

### D. Per-method explainer figures
Three small Typst+CeTZ figures, one per discovery:
- `figures/discovery_law_matrix.typ` — visualisation of which LAW × PROP × FAIL combos parse, abstracted from `references/openradioss_law_compatibility_matrix.md`.
- `figures/discovery_orientation_convention.typ` — the working recipe schematic (Phi vs /SKEW vs /INIBRI/ORTHO) with the verified 0.24% match noted.
- `figures/discovery_kok_port.typ` — the M1-M4 milestone timeline plus the stage-11 stiffness gap (Ex 17 vs target 53 GPa).

## Per-stage scope

| Stage | Verdict | Composite PNG (A) | HD MP4 (B) | Notes |
|---|---|---|---|---|
| 01 | PASS | yes | yes (regenerate 50 frames) | beam bending |
| 02 | INCONCLUSIVE | yes (alpha=1 only) | yes | cantilever, alpha=1 only |
| 03 | PASS | yes | yes | dogbone |
| 04 | PASS | yes | yes | open hole |
| 05 | PASS | yes | yes | damage |
| 06 | PASS | yes | yes | failure-card probe |
| 07 | PASS | yes | skip MP4 if no anim frames | UD tow |
| 08 | PASS | yes | skip MP4 | ply rotation sweep — Typst plot is primary |
| 09 | PASS | yes | skip MP4 | laminate CLT |
| 10-16 | INCONCLUSIVE / FAIL | wireframe only if VTK exists | skip | document blocker, do not rerun the failed solver attempt |

## Workflow

1. `git pull --rebase`.
2. Process stages in order 01 → 16. Each stage:
   - Check existing artifacts; do not overwrite if they're already at HD quality.
   - Build the new deliverables under `tests/stage_NN_*/figures/` (composite PNG and HD MP4 with `_hd` or `_composite` suffixes).
   - For stages requiring deck regeneration (multi-frame anim output): write a separate `runner_pubviz.py` that copies the pristine runner, modifies only `/ANIM/DT`, runs in `runs/<rootname>_pubviz/`. Do not pollute verification artifacts.
3. Commit per-stage `pubviz stage NN: composite + HD motion`.
4. After all stages: build the cross-stage summary and discovery explainers (C, D), commit `pubviz: cross-stage summary and methodology figures`.
5. Update `RESULTS_REPORT.md` to reference the new images at the top of each stage row (replace the existing `[PNG]` and `[MP4]` links with the `_composite.png` and `_hd.mp4` ones). Commit `final: results report figures upgrade to publication-grade`.
6. Push, exit cleanly.

## Constraints (hard)

- Solid HEXA8 / TETRA10 only in any deck regeneration.
- Brand colours strict; no rounded corners; no matplotlib.
- Typst + CeTZ for 2D plots; PyVista headless OSMesa for 3D renders; ffmpeg for MP4 muxing.
- No GUI, no closed source.
- Author = `j-vaught <jvaught@sc.edu>`.
- Do not change verification verdicts or stage results.json. Only add figures and update the report's figure links.

## When to stop

- All stages have composite/HD figures (where applicable) + summary + discovery explainers committed → exit.
- ffmpeg/PyVista/Typst tooling fault you cannot work around → write `figures/_pubviz_blocker.md`, commit, push, exit.
- Token / time budget — push everything, exit gracefully.

Begin.
