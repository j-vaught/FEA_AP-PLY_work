# Stage 06 Figure Gallery

The current `results/results.json` verdict is INCONCLUSIVE. Stage 06 verdict is INCONCLUSIVE due to LAW25-on-TYPE14 starter incompatibility; the figures here show /PROP/TYPE6 proxy runs, not the canonical solid path. The standard rich-visualization set reads existing outputs only: TYPE6/SOL_ORTH proxy animation frames for TSAIWU, HASHIN, and PUCK; canonical LAW25+TYPE14 is blocked.

`stage_06_mesh_wireframe.png` uses the first available frame (`TSAIWU A001` from `TSAIWU TYPE6 proxy`) reconstructed to reference coordinates when the `Displacement` vector is present. It is an iso-view garnet wireframe on white with axes in metres; verdict context is INCONCLUSIVE.

`stage_06_field_displacement.png` uses `PUCK A004` from `PUCK TYPE6 proxy` at VTK time 0.000819745 s, warped 6.34x by the selected vector field and coloured by point Displacement magnitude in millimetres. The colourbar is in millimetres, the mesh edges are charcoal, the axes are metres, and the verdict context is INCONCLUSIVE.

`stage_06_history_disp.typ` / `stage_06_history_disp.pdf` is a Typst + CeTZ history plot built from `results/timeseries.csv` using brand-colour curves. The plotted scalar is stage-specific: stress-strain for Stage 03, Kt for Stage 04, force-displacement for Stage 05, and the recorded analytic failure envelope for Stage 06; verdict context is INCONCLUSIVE.

`stage_06_motion_displacement.mp4` is a 1920x1080 H.264 loop at 24 fps using all 12 rendered frame(s), repeated to 5.0 s. The mesh is warped 6.34x and coloured by point Displacement magnitude in millimetres with the requested brand ramp; verdict context is INCONCLUSIVE.

`stage_06_field_vonmises.png` uses `PUCK A004` at VTK time 0.000819745 s, the same vector warp, and cell von Mises stress in MPa computed from the 3D-element stress tensor. Axes are metres, edges are charcoal, and the verdict context is INCONCLUSIVE.

`stage_06_motion_vonmises.mp4` uses the same frame sequence and warp as the displacement motion video, coloured by cell von Mises stress in MPa with a fixed colour range over the 5.0 s loop; verdict context is INCONCLUSIVE.

`stage_06_motion_doubleview.mp4` is a synchronized 1920x1080 split view of the same frame sequence: vector magnitude in millimetres on the top panel and von Mises stress in MPa on the bottom panel. The lower label reports the source frame and VTK time; verdict context is INCONCLUSIVE.

`stage_06_field_strain.png` uses `PUCK A004` at VTK time 0.000819745 s and the same vector warp, coloured by cell plastic strain scalar from the VTK output. The scalar is dimensionless, axes are metres, edges are charcoal, and verdict context is INCONCLUSIVE.

Stage 06 is INCONCLUSIVE because the required LAW25 + TYPE14 canonical solid path is rejected by the starter. These figures intentionally show the TYPE6/SOL_ORTH proxy outputs and analytic diagnostic rows that were actually produced; they are not evidence that the canonical laminate solid path passed.
