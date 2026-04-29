# Stage 02 Figure Gallery

`stage_02_mesh_wireframe.png`, `stage_02_deformed.png`, `stage_02_no_animation.txt`, and `stage_02_timeseries.typ` / `stage_02_timeseries.pdf` are the original verification-view figures retained for provenance. The current `results/results.json` verdict is INCONCLUSIVE only because alpha=3 and alpha=5 remain blocked; the baseline alpha=1 row is PASS with `dx/L` error 1.574% and `dy/L` error 0.871% against the Bisshopp-Drucker elastica reference.

`stage_02_motion_displacement.mp4` uses the gated baseline alpha=1 rich series from `runs/stage02_implicit_baseline_a1p00_richviz/`, frames A001-A061 over 0.0-1.0 s. The mesh is reconstructed to reference coordinates, warped 1x by the `Displacement` vector, and coloured by displacement magnitude in millimetres with the Sandstorm-Honeycomb-Garnet-Congaree ramp; context is the alpha=1 PASS elastica check.

`stage_02_motion_vonmises.mp4` uses the same alpha=1 rich frames and 1x displacement warp, coloured by cell von Mises stress in MPa computed from the 3D Cauchy stress tensor. The colour range is held constant through the 5 s H.264 loop, so stress growth remains comparable over the implicit load ramp.

`stage_02_motion_doubleview.mp4` is a synchronized 1920x1080 split view of the alpha=1 rich series: displacement magnitude in millimetres on the top panel and von Mises stress in MPa on the bottom panel. The lower-right time label reports the OpenRadioss frame time; alpha=3 and alpha=5 are intentionally absent from this visualization set.

`stage_02_field_displacement.png` is the final alpha=1 rich frame A061 at 1.0 s, reconstructed to reference coordinates, warped 1x by `Displacement`, and coloured by displacement magnitude in millimetres. Axes are in metres and the inline annotation records the alpha=1 PASS context, 1.574% maximum gated relative error, and 84.4 s wall clock from the verification result.

`stage_02_field_vonmises.png` is the same final frame and warp, coloured by von Mises stress in MPa from the averaged 3D-element Cauchy stress tensors. The plotted state is the alpha=1 baseline only.

`stage_02_field_strain.png` is the same final frame and warp, coloured by von-Mises-equivalent strain computed from the 3D-element strain tensors. Strain is dimensionless, and the verdict context remains alpha=1 PASS.

`stage_02_field_xsection.png` slices the final alpha=1 rich frame at mid-span (`x = 0.500 m`) and colours the through-section by von Mises stress in MPa. It shows the bending stress gradient in the large-deflection cantilever while avoiding any alpha=3 or alpha=5 rerun.

`stage_02_history.typ` / `stage_02_history.pdf` and `stage_02_history.csv` plot the alpha=1 rich frame history in Typst + CeTZ: tip sag in millimetres versus simulation time, tip-face resultant force in newtons, and total energy (`internal + kinetic`) in millijoules. Garnet denotes OpenRadioss histories, Atlantic denotes the reference/final target, and Horseshoe marks the 2% displacement tolerance band.
