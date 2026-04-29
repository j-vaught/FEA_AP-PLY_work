# Stage 04 Figure Gallery

The current `results/results.json` verdict is PASS. Stage 04 verdict is PASS; the primary gate is the Ntheta 64 Howland Kt comparison. The standard rich-visualization set reads existing outputs only: mesh-refinement final frames for Ntheta 32, 64, and 128; this is a diagnostic sweep, not a dynamic time sequence.

`stage_04_mesh_wireframe.png` uses the first available frame (`Ntheta 32` from `stage04_ntheta032`) reconstructed to reference coordinates when the `Displacement` vector is present. It is an iso-view garnet wireframe on white with axes in metres; verdict context is PASS.

`stage_04_field_displacement.png` uses `Ntheta 128` from `stage04_ntheta128` at VTK time 1 s, warped 50x by the selected vector field and coloured by point Displacement magnitude in millimetres. The colourbar is in millimetres, the mesh edges are charcoal, the axes are metres, and the verdict context is PASS.

`stage_04_history_disp.typ` / `stage_04_history_disp.pdf` is a Typst + CeTZ history plot built from `results/timeseries.csv` using brand-colour curves. The plotted scalar is stage-specific: stress-strain for Stage 03, Kt for Stage 04, force-displacement for Stage 05, and the recorded analytic failure envelope for Stage 06; verdict context is PASS.

`stage_04_motion_displacement.mp4` is a 1920x1080 H.264 loop at 24 fps using all 3 rendered frame(s), repeated to 5.0 s. The mesh is warped 50x and coloured by point Displacement magnitude in millimetres with the requested brand ramp; verdict context is PASS.

`stage_04_field_vonmises.png` uses `Ntheta 128` at VTK time 1 s, the same vector warp, and cell von Mises stress in MPa computed from the 3D-element stress tensor. Axes are metres, edges are charcoal, and the verdict context is PASS.

`stage_04_motion_vonmises.mp4` uses the same frame sequence and warp as the displacement motion video, coloured by cell von Mises stress in MPa with a fixed colour range over the 5.0 s loop; verdict context is PASS.

`stage_04_motion_doubleview.mp4` is a synchronized 1920x1080 split view of the same frame sequence: vector magnitude in millimetres on the top panel and von Mises stress in MPa on the bottom panel. The lower label reports the source frame and VTK time; verdict context is PASS.

`stage_04_field_strain.png` uses `Ntheta 128` at VTK time 1 s and the same vector warp, coloured by cell von-Mises-equivalent strain computed from 3D-element strain tensors. The scalar is dimensionless, axes are metres, edges are charcoal, and verdict context is PASS.
