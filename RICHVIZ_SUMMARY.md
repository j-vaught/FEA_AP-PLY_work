# Rich Visualization Summary

Generated on 2026-04-29 for the already gated cases only: Stage 01 M1 3-point / 4-point bending and Stage 02 baseline alpha=1. Stage 02 alpha=3 and alpha=5 were not rerun.

## Stage 01 Beam Bending

Context: `results/results.json` verdict is PASS. The gated M1 relative displacement error is 0.968% for both 3-point and 4-point bending against Euler-Bernoulli.

Rich reruns:

- `tests/stage_01_beam_bending/runner_richviz.py`
- `tests/stage_01_beam_bending/runs/stage01_3pt_M1_richviz/` generated 61 animation/VTK frames over 0.0-0.00600004 s
- `tests/stage_01_beam_bending/runs/stage01_4pt_M1_richviz/` generated 61 animation/VTK frames over 0.0-0.00600001 s
- committed representative VTKs: `stage01_3pt_M1_richvizA061.vtk`, `stage01_4pt_M1_richvizA061.vtk`

Published artifacts:

- `tests/stage_01_beam_bending/figures/stage_01_motion_displacement.mp4`
- `tests/stage_01_beam_bending/figures/stage_01_motion_vonmises.mp4`
- `tests/stage_01_beam_bending/figures/stage_01_motion_doubleview.mp4`
- `tests/stage_01_beam_bending/figures/stage_01_field_displacement.png`
- `tests/stage_01_beam_bending/figures/stage_01_field_vonmises.png`
- `tests/stage_01_beam_bending/figures/stage_01_field_strain.png`
- `tests/stage_01_beam_bending/figures/stage_01_field_xsection.png`
- `tests/stage_01_beam_bending/figures/stage_01_history.typ`
- `tests/stage_01_beam_bending/figures/stage_01_history.pdf`
- `tests/stage_01_beam_bending/figures/stage_01_history.csv`

## Stage 02 Cantilever Large Deflection

Context: `results/results.json` remains INCONCLUSIVE because alpha=3 and alpha=5 are blocked, but the baseline alpha=1 gated metric is PASS (`dx/L` error 1.574%, `dy/L` error 0.871%).

Rich rerun:

- `tests/stage_02_cantilever_large/runner_richviz.py`
- `tests/stage_02_cantilever_large/runs/stage02_implicit_baseline_a1p00_richviz/` generated 61 animation/VTK frames over 0.0-1.0 s
- committed representative VTK: `stage02_implicit_baseline_a1p00_richvizA061.vtk`

Published artifacts:

- `tests/stage_02_cantilever_large/figures/stage_02_motion_displacement.mp4`
- `tests/stage_02_cantilever_large/figures/stage_02_motion_vonmises.mp4`
- `tests/stage_02_cantilever_large/figures/stage_02_motion_doubleview.mp4`
- `tests/stage_02_cantilever_large/figures/stage_02_field_displacement.png`
- `tests/stage_02_cantilever_large/figures/stage_02_field_vonmises.png`
- `tests/stage_02_cantilever_large/figures/stage_02_field_strain.png`
- `tests/stage_02_cantilever_large/figures/stage_02_field_xsection.png`
- `tests/stage_02_cantilever_large/figures/stage_02_history.typ`
- `tests/stage_02_cantilever_large/figures/stage_02_history.pdf`
- `tests/stage_02_cantilever_large/figures/stage_02_history.csv`

## Rendering Notes

All motion videos are H.264 MP4 at 1920x1080, 24 fps, 5.0 s. All field PNGs are 1800x1350. Colour maps use the requested brand ramp `#FFF2E3 -> #A49137 -> #73000A -> #1F414D`; Typst + CeTZ histories use Garnet for OpenRadioss, Atlantic for reference/final targets, and Horseshoe for tolerance bands.

Raw `.gz` animation files and solver restart/output files remain local under `_richviz` run directories and were not committed. The renderer used for both stages is `tools/richviz_render.py`.
