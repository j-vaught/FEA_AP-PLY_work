# Rich Visualization All-Stages Summary

Generated on 2026-04-29 from existing VTK and compressed OpenRadioss animation outputs. No simulation decks were edited and no blocked stage was rerun.

## Completed Figure Sets

Stage 01 and Stage 02 were already complete before this pass, with figure sets under `tests/stage_01_beam_bending/figures/` and `tests/stage_02_cantilever_large/figures/`.

Stage 03 `iso_dogbone_E8` was rendered from the existing medium-mesh load-factor VTK sweep. Commit: `6519795` (`richviz stage 03: figures for PASS outputs`). Artifacts include `stage_03_mesh_wireframe.png`, `stage_03_field_displacement.png`, `stage_03_history_disp.typ/.pdf`, `stage_03_motion_displacement.mp4`, `stage_03_field_vonmises.png`, `stage_03_motion_vonmises.mp4`, `stage_03_motion_doubleview.mp4`, `stage_03_field_strain.png`, and `figures/README.md`.

Stage 04 `open_hole_kirsch` was rendered from the existing Ntheta 32/64/128 diagnostic final-frame sweep. Commit: `331e12d` (`richviz stage 04: figures for PASS outputs`). The MP4s are documented as a mesh-refinement diagnostic sequence, not a dynamic time animation.

Stage 05 `dogbone_damage` was rendered from the existing 60-frame medium-mesh LAW22 sequence. Commit: `8a6f841` (`richviz stage 05: figures for PASS outputs`). The history plot compares the coarse, medium, and fine force-displacement histories from `results/timeseries.csv`.

Stage 06 `composite_failure_criteria` was rendered from the existing TYPE6/SOL_ORTH proxy outputs for TSAIWU, HASHIN, and PUCK. Commit: `ebbfa41` (`richviz stage 06: figures for INCONCLUSIVE outputs`). The README states that Stage 06 remains INCONCLUSIVE because LAW25 + TYPE14 is rejected by the starter, and these figures are proxy outputs rather than the canonical solid path.

## Not Rendered

Stages 07, 08, and 09 currently have INCONCLUSIVE verdicts on disk but no `.vtk` or `.gz` animation files under their `runs/` directories, so there was no existing field output to render in this pass.

Stages 10 through 16 also had no `.vtk` or `.gz` animation files under `runs/` at the final scan.

## Renderer

The reusable post-processing entry point added in this pass is `tools/richviz_render_existing.py`. It uses PyVista off-screen rendering for 3D figures and videos, Typst + CeTZ for history plots, and ffmpeg for H.264 MP4 encoding. Missing VTK frames from compressed animation files are converted in a temporary scratch directory only; no new solver runs are launched.
