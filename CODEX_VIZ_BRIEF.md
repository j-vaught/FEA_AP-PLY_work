# Codex visualization brief — FEA_AP-PLY

You are a dedicated visualization codex running in tmux session `codex-viz` on `comech-2422`. Your job is to produce publication-grade figures for each stage of the FEA test ladder as the main codex (running in `codex-fea`) produces VTK / CSV outputs. Operate independently — do not modify simulation code or reruns.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. Do not mention any AI provider in commits.

## Required reading

1. `README.md` — repo overview.
2. `plan/master_plan.md` — committed direction and the role of figures.
3. `~/.claude/CLAUDE.md` (user instructions) — brand-colour palette and figure-style rules. **Brand colours only**, no rounded edges, no matplotlib, all quantitative plots in **Typst + CeTZ**. PNG / MP4 for 3D renders is fine via PyVista.
4. `tests/stage_NN_<name>/spec.md` — read each stage's spec for what to visualize before rendering.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
# OSMesa headless render is enabled (mesalib package is installed in feaapply env).
# vtkRenderWindow().OffScreenRenderingOn() works without an X server.
```

Test rig:
```python
import vtk
rw = vtk.vtkRenderWindow(); rw.OffScreenRenderingOn(); rw.Render()
# expected: prints "VTK headless render OK", may emit a benign DISPLAY warning
```

## Per-stage deliverables

For each stage `NN` that has at least one `.vtk` file in `tests/stage_NN_<name>/runs/<rootname>/`:

1. **Mesh wireframe** (iso view) — `tests/stage_NN_<name>/figures/stage_NN_mesh_wireframe.png`
   - Stroke colour: garnet `#73000A` (UofSC primary).
   - Background: white. Axes legend bottom-left. Bounds box thin charcoal `#363636`.
   - Square aspect / no rounded corners.

2. **Deformed view, colour-mapped by displacement magnitude** — `tests/stage_NN_<name>/figures/stage_NN_deformed.png`
   - If a `Displacement` or `U` array exists in the VTK, warp the mesh by it (PyVista `mesh.warp_by_vector('Displacement')`).
   - If only `Velocity` is available (which is the OpenRadioss anim default), DO NOT warp — instead colour the undeformed mesh by velocity magnitude and label the colorbar accordingly. Do not invent a displacement field.
   - Colormap: monotone garnet→atlantic ramp (colors `#FFF2E3`, `#A49137`, `#73000A`, `#1F414D` interpolated). Colourbar with units (mm or m/s).
   - Edge lines drawn at stroke `#363636`, line width 0.4.

3. **Multi-frame animation** — `tests/stage_NN_<name>/figures/stage_NN_animation.mp4` (if more than one anim frame exists)
   - Detect `<rootname>A001.gz`, `A002.gz`, … — gunzip each, run `$OR/exec/anim_to_vtk_linux64_gf` to get `<rootname>A00N.vtk`, then build a frame sequence in PyVista at 24 fps.
   - Same colourmap as the deformed view.
   - 4 s minimum (loop the sequence if necessary). 1080p (1920×1080).
   - **If only one frame exists, skip the MP4** — write a one-line `figures/stage_NN_no_animation.txt` explaining only A001 was written.

4. **Time-series plot** — `tests/stage_NN_<name>/figures/stage_NN_timeseries.typ`
   - Source Typst + CeTZ from the stage's `results/timeseries.csv`.
   - Plot the primary scalar quantity vs. mesh-or-time (for stage 01: rel_error_eb vs. mesh, log scale on Y if helpful).
   - Use brand colours per the global palette. Garnet for mismatch / fail, horseshoe for spec gate, atlantic for closed-form reference.
   - Render the Typst source with `typst compile` to a sibling PDF, but commit both `.typ` source and `.pdf`.

5. **Stage gallery summary** — `tests/stage_NN_<name>/figures/README.md`
   - One paragraph stating what each figure shows, which timestep the snapshot represents, and the verdict from `results.json`.

## Workflow

1. `git pull --rebase` (sync with main codex's commits).
2. List `tests/stage_*/runs/*/` for any stage that has a `.vtk` file.
3. For each such stage, check whether `tests/stage_NN_*/figures/stage_NN_mesh_wireframe.png` already exists. **Skip** stages that already have the wireframe (figures are already done).
4. For new stages, produce all five deliverables in order (1, 2, 3, 4, 5).
5. After each stage's figures are written: `git add tests/stage_NN_*/figures/`, commit `"viz stage NN: wireframe, deformed, anim, timeseries"`, push.
6. After all currently-available stages are figured: poll `tests/stage_*/runs/` every 5 minutes for 1 hour. If new VTK output appears for a stage that lacks figures, render it. Otherwise exit gracefully after the polling window.

## Constraints (hard)

- No matplotlib. No seaborn. No plotly.
- Typst + CeTZ only for quantitative plots. PyVista headless for 3D renders. Brand colours.
- No GUI — the env is headless. Use `pyvista.start_xvfb()` only if `mesalib`-based OSMesa fails (it shouldn't — but if it does, fall back to writing OBJ + glTF and document as `figures/stage_NN_render_failed.md`).
- Do not modify simulation code, runners, deck templates, or anything outside `tests/stage_*/figures/`.
- Do not duplicate work in `codex-fea`'s lane.
- Author identity = `j-vaught <jvaught@sc.edu>`.
- All units explicit on every figure.

## When to stop

- Polling window expired (1 h with no new VTK in stages that lack figures) — write `RESULTS_SUMMARY_VIZ.md` with the list of figures rendered and exit.
- Hard error on PyVista / VTK that you cannot work around — write `figures/_blocker.md` and exit.
- Stages 1-16 all have complete figure sets — write summary and exit.

Begin.
