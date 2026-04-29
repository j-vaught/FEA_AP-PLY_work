# Codex rich-visualisation brief — FEA_AP-PLY

You run autonomously on `comech-2422` in tmux session `codex-richviz`. Your job is to produce **publication-quality, motion-rich, multi-field** figures and videos for stages that have already been verified (currently stage 01 PASS, stage 02 α=1 PASS / α=3 blocked). Operate in your own lane; do not modify simulation runners belonging to `codex-fea`'s active stage work.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. Do not mention any AI provider.

## Why this exists

The first viz pass produced single-frame MP4s because the verification decks only wrote one animation frame (`/ANIM/DT` was set sparse to keep file sizes small for the gate check). The figures show stress/strain/displacement only via the inadequate default `Velocity` field. We need real motion + real fields, on the gated cases.

## Required reading

1. `README.md`, `plan/master_plan.md`, `CODEX_BRIEF.md`.
2. `tests/stage_NN_<name>/spec.md` — for what each stage is showing.
3. `tests/stage_NN_<name>/runner.py` — only to understand how the deck is templated; you may COPY (do not edit) it into a `runner_richviz.py` sibling that writes a richer deck.
4. The previously-committed figures under `tests/stage_NN_<name>/figures/` — keep them; add new ones with `_motion_*` and `_field_*` prefixes.
5. `~/.claude/CLAUDE.md` — brand-colour palette and no-rounded-corner rules.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
# OSMesa headless render is enabled (mesalib package).
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss     # MUMPS-linked binaries live here
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH
```

Use `pgrep -a engine_linux64` to detect when `codex-fea` is running its own engine — back off with `time.sleep(60)` rather than competing for the same node. The host has 112 cores; one solver run uses 1 thread, but staggering avoids cache thrash.

## What to do per stage

For every stage `NN` with a `tests/stage_NN_<name>/results/results.json` whose `verdict` is `PASS` or where at least one gated metric is PASS (stage 02 has α=1 PASS):

### A. Rich rerun

1. Make a copy of the existing runner: `cp runner.py runner_richviz.py`.
2. In `runner_richviz.py`, expand the engine-deck template's animation cards to:
   ```text
   /ANIM/DT
   <run_time / 60>
   /ANIM/BRICK/STRS                # 3D-element Cauchy stress, all components
   /ANIM/BRICK/STRA/STRAIN/ALL     # 3D-element strain, all components
   /ANIM/BRICK/EPSP                # plastic strain (where defined)
   /ANIM/BRICK/DAMA                # damage variable (where defined)
   /ANIM/NODA/DIS                  # nodal displacement vector
   /ANIM/NODA/VEL                  # nodal velocity vector
   /ANIM/GZIP                      # default — keep .gz
   ```
   Goal: ~50-60 evenly-spaced frames over the gated alpha (or per the existing run-time).
3. Re-run only the gated cases. Stage 01: 3pt and 4pt at gated mesh M1 (already passing). Stage 02: α=1 baseline (already passing). DO NOT attempt α=3, α=5 — those belong to `codex-fea`.
4. Each rich rerun lives in `tests/stage_NN_<name>/runs/<rootname>_richviz/` to keep verification artifacts pristine.

### B. Render

For each `<rootname>_richviz/<rootname>A0*.gz` series:

1. `gunzip` and convert with `$OR/exec/anim_to_vtk_linux64_gf` to `.vtk` per frame.
2. Load all frames into PyVista. Compute per-frame:
   - `displacement_magnitude` from the `Displacement` vector (in mm or m, label clearly).
   - `von_mises` from the 3D-element Cauchy stress components (Tr σ → σ_dev → von Mises). PyVista has helpers; cite the formula in code comments only if needed.
   - `eq_strain` (von-Mises-equivalent strain) from the strain tensor for stages with damage.
   - `plastic_strain` directly from the array if present.

3. Produce these MP4s (1080p, 24 fps, 4-8 s; loop the sequence until reaching the duration; H.264 via ffmpeg):

   - `tests/stage_NN_<name>/figures/stage_NN_motion_displacement.mp4` — mesh **warped by `Displacement`**, coloured by displacement magnitude. Colourbar in mm with brand-colour ramp `#FFF2E3 → #A49137 → #73000A → #1F414D`.
   - `tests/stage_NN_<name>/figures/stage_NN_motion_vonmises.mp4` — same warp, coloured by von Mises stress. Colourbar in MPa.
   - `tests/stage_NN_<name>/figures/stage_NN_motion_doubleview.mp4` — split-screen 1920×1080: top half displacement, bottom half von Mises, synchronized timestep, time-of-frame label in lower-right.

4. Produce these stills (PNG, ≥1600×1200):

   - `tests/stage_NN_<name>/figures/stage_NN_field_displacement.png` — final frame, warped, displacement-coloured, with axes (m), colourbar with units, and a small inline annotation of `verdict, rel_error, wall_clock` from `results.json`.
   - `tests/stage_NN_<name>/figures/stage_NN_field_vonmises.png` — final frame, von-Mises-coloured.
   - `tests/stage_NN_<name>/figures/stage_NN_field_strain.png` — final frame, strain-coloured (eq strain or principal strain).
   - `tests/stage_NN_<name>/figures/stage_NN_field_xsection.png` — mid-span cross-section showing the through-thickness gradient of σ or ε. Use a clipping plane.

### C. Time-series, in Typst + CeTZ

For each stage, augment the existing `stage_NN_timeseries.typ` (or write a new `stage_NN_history.typ`) with:

1. Midspan / probe-node displacement vs simulation time.
2. Reaction force at the loaded face vs simulation time.
3. Total energy = internal + kinetic vs time (sanity check, especially for explicit DR runs).

Brand colours: garnet for FEM curves, atlantic for closed-form / reference, horseshoe for the spec gate. No matplotlib. Compile to `.pdf` and commit both `.typ` and `.pdf`.

### D. README update

Update each stage's `tests/stage_NN_<name>/figures/README.md` with one paragraph per new figure: which timestep, which colour mapping, which physical units, and the verdict context.

## Workflow

1. `git pull --rebase` at the start.
2. Process stages in numerical order (01 first, then 02 α=1).
3. After each stage:
   - `git add tests/stage_NN_*/runs/*_richviz/ tests/stage_NN_*/figures/`.
   - Commit `richviz stage NN: motion videos + multi-field figures`.
   - Push.
4. After both stages are done, write `RICHVIZ_SUMMARY.md` listing every artifact produced and exit.

## Constraints (hard)

- Solid HEXA8 / TETRA10 only. Do not change element types.
- No matplotlib. No GUI tooling. No closed-source software.
- Brand colours strict; no rounded corners on any plot frame.
- Keep existing committed figures intact; only add `_motion_*`, `_field_*`, `_history_*` siblings.
- Author = `j-vaught <jvaught@sc.edu>`.
- Do NOT touch simulation runners belonging to `codex-fea`'s active stages 03-16. Only rich-rerun cases that have already passed.
- Do not run two engines at the same time as `codex-fea`'s solver — back off with sleep when you see one running.
- Do not commit `.gz` raw anim files (keep them under `_richviz/`); only commit derived `.vtk` (one or two per stage as small representative samples), `.png`, `.mp4`, `.typ`, `.pdf`.

## When to stop

- Both stages (01 PASS, 02 α=1 PASS) have full rich-viz artifact sets — write `RICHVIZ_SUMMARY.md`, commit, push, exit.
- Hard tooling error (PyVista / ffmpeg / OSMesa) you cannot work around — write `figures/_richviz_blocker.md`, push, exit.
- Token / time budget — push everything and exit gracefully.

Begin.
