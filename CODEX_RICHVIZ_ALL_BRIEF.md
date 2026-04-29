# Codex rich-visualisation ALL-STAGES brief — FEA_AP-PLY

Replaces `CODEX_RICHVIZ_BRIEF.md` for this run. **Render figures for every stage that has VTK output, regardless of verdict** — not just PASS / partial-PASS stages. The user wants to see what was actually run for INCONCLUSIVE stages too (the proxy runs, the diagnostic runs, the not-quite-passing iterations).

You run autonomously on `comech-2422` in tmux session `codex-richviz`. Operate in your own lane; do not modify simulation runners belonging to active codex-lawmatrix or codex-fea sessions.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. Do not mention any AI provider.

## What's already done (skip these — pristine)

- `tests/stage_01_beam_bending/figures/` — wireframe + deformed + motion videos committed at `213e6c3`.
- `tests/stage_02_cantilever_large/figures/` — committed at `6f53873`.

## What to do for stages 03-09 (currently exist with verdicts on disk)

For **every** stage `NN ∈ {03, 04, 05, 06, 07, 08, 09}` (and any later stages when they appear) with at least one `.vtk` or `.gz` anim file under `tests/stage_NN_<name>/runs/<rootname>/`:

1. Skip if `tests/stage_NN_<name>/figures/stage_NN_motion_displacement.mp4` already exists (already done).
2. Otherwise produce, in order of priority:

   ### Priority A — Fast, always do these
   - `tests/stage_NN_<name>/figures/stage_NN_mesh_wireframe.png` — iso-view wireframe in garnet `#73000A` on white, units in metres on axes. From the FIRST anim frame's geometry.
   - `tests/stage_NN_<name>/figures/stage_NN_field_displacement.png` — final-frame mesh warped by `Displacement`, coloured by displacement magnitude in mm with the brand-colour ramp `#FFF2E3 → #A49137 → #73000A → #1F414D`. If `Displacement` array isn't present, use whatever vector field IS in the VTK and label honestly (Velocity, etc.). 1600×1200, edges in `#363636`.
   - `tests/stage_NN_<name>/figures/stage_NN_history_disp.typ` + compiled `.pdf` — Typst+CeTZ time-history of the primary scalar (midspan deflection, K_t, peak load, etc., depending on stage) from `results/timeseries.csv`. Brand colours.

   ### Priority B — If multiple frames exist
   - `tests/stage_NN_<name>/figures/stage_NN_motion_displacement.mp4` — 1080p / 24 fps, all frames, warped by Displacement, coloured by displacement magnitude. 4-8 s minimum (loop frames if needed).
   - `tests/stage_NN_<name>/figures/stage_NN_field_vonmises.png` — final-frame coloured by von Mises stress in MPa, computed from `3DELEM_Strs_Intg_Point111__` cell-array if present.
   - `tests/stage_NN_<name>/figures/stage_NN_motion_vonmises.mp4` — same warp, von Mises coloured.

   ### Priority C — Only if Priority A+B succeed and time budget allows
   - `tests/stage_NN_<name>/figures/stage_NN_motion_doubleview.mp4` — split-screen 1920×1080, top half displacement / bottom half von Mises.
   - `tests/stage_NN_<name>/figures/stage_NN_field_strain.png` — final-frame strain-coloured.

3. Update or write `tests/stage_NN_<name>/figures/README.md` with one paragraph per figure: which timestep, which colour mapping, which physical units, the verdict from `results.json`, and an honest note about INCONCLUSIVE stages (e.g. *"Stage 06 verdict is INCONCLUSIVE due to LAW25-on-TYPE14 starter incompatibility; the figures here show the /PROP/TYPE6 proxy runs, not the canonical solid path."*).

4. Commit `richviz stage NN: figures for <verdict> outputs` after each stage. Push.

## Constraints (hard)

- Solid HEXA8 / TETRA10 only in any rich-rerun deck modifications.
- **Do NOT re-run failed/blocked stages**. For stages 06-09 (LAW25+TYPE14 blocked), **render whatever VTK output already exists** under `runs/`, even if it's the proxy run. Document honestly.
- No matplotlib, no GUI. Typst+CeTZ for plots. PyVista for 3D.
- Brand colours strict; no rounded corners.
- Author = `j-vaught <jvaught@sc.edu>`.
- Do NOT touch simulation runners belonging to `codex-fea` (when it relaunches) or `codex-lawmatrix`.

## When to stop

- Stages 01-09 all have figure sets — write `RICHVIZ_ALL_SUMMARY.md` and exit.
- New stages 10+ land while you're running — process them too.
- After 2 hours total wall-clock, push everything and exit gracefully.

Begin.
