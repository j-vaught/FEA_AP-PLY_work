# Stage 01 Figure Gallery

`stage_01_mesh_wireframe.png`, `stage_01_deformed.png`, `stage_01_animation.mp4`, and `stage_01_timeseries.typ` / `stage_01_timeseries.pdf` are the original verification-view figures retained for provenance. The current `results/results.json` verdict is PASS: gated M1 relative displacement errors are 0.968% for both 3-point and 4-point bending against the Euler-Bernoulli reference.

`stage_01_motion_displacement.mp4` uses the gated 4-point M1 rich series from `runs/stage01_4pt_M1_richviz/`, frames A001-A061 over 0.0-0.00600001 s. The mesh is reconstructed to reference coordinates, warped 8x by the `Displacement` vector, and coloured by displacement magnitude in millimetres with the Sandstorm-Honeycomb-Garnet-Congaree ramp; context is the PASS 4pt_M1 displacement gate with 0.968% relative error.

`stage_01_motion_vonmises.mp4` uses the same 4-point M1 rich frames and warp, coloured by cell von Mises stress in MPa computed from the 3D Cauchy stress tensor. The colour range is held constant through the 5 s H.264 loop, so the stress build-up remains comparable from frame to frame under the PASS 4pt_M1 verdict context.

`stage_01_motion_doubleview.mp4` is a synchronized 1920x1080 split view of the same rich series: displacement magnitude in millimetres on the top panel and von Mises stress in MPa on the bottom panel. The lower-right time label reports the OpenRadioss frame time, with both panels using the same 8x displacement warp and PASS 4pt_M1 context.

`stage_01_field_displacement.png` is the final rich frame A061 at 0.00600001 s, reconstructed to reference coordinates, warped 8x by `Displacement`, and coloured by displacement magnitude in millimetres. Axes are in metres and the inline annotation records the 4pt_M1 gated verdict, 0.968% relative error, and 451.6 s verification wall clock.

`stage_01_field_vonmises.png` is the same final frame and warp, coloured by von Mises stress in MPa from the 3D-element Cauchy stress tensor. The plotted state corresponds to the PASS 4pt_M1 Euler-Bernoulli gate and uses the same brand ramp as the motion stress video.

`stage_01_field_strain.png` is the same final frame and warp, coloured by von-Mises-equivalent strain computed from the 3D-element strain tensor. Strain is dimensionless, and the verdict context is the PASS 4pt_M1 displacement gate.

`stage_01_field_xsection.png` slices the final 4-point M1 rich frame at mid-span (`x = 0.100 m`) and colours the through-thickness section by von Mises stress in MPa. It highlights the expected bending gradient across the 10 mm thickness under the same PASS 4pt_M1 context.

`stage_01_history.typ` / `stage_01_history.pdf` and `stage_01_history.csv` plot the rich 4-point M1 frame history in Typst + CeTZ: midspan displacement in millimetres versus simulation time, loaded-line resultant force in newtons, and total energy (`internal + kinetic`) in millijoules. Garnet denotes OpenRadioss histories, Atlantic denotes the reference/final target, and Horseshoe marks the 1% displacement tolerance band.
