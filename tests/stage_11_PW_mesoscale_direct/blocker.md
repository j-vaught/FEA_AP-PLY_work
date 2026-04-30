# Stage 11 Blocker - AP-PLY stiffness gate

Author: J.C. Vaught

The Kok geometry port generated the 25 mm x 25 mm, four-ply `[0/+45/-45/90]` AP-PLY block and the OpenRadioss decks used `LAW12 + TYPE6/SOL_ORTH` with per-tow `/INIBRI/ORTHO`. The initial TETRA10 formulation `Itetra10=2` was incompatible with imposed kinematic conditions in this build (`ERROR ID : 1208`), so the validation decks were rerun with the documented alternate TETRA10 formulation `Itetra10=1000`. Starter and engine then completed for all three load cases.

The solver completed status and numerical comparison are:

- E_x_GPa: measured `16.705859039432056` GPa, target `53.3` GPa, error `68.65692487911434` percent.
- E_y_GPa: measured `16.96796806103866` GPa, target `53.3` GPa, error `68.16516311249782` percent.
- G_xy_GPa: measured `0.02289196065608415` GPa, target `20.5` GPa, error `99.88833189923861` percent.

This is not a near-tolerance material-system discrepancy. The current clean-room port represents each ply as discrete straight tows plus pure-resin remainder. For the stage-11 one-over-one 25 mm block this creates a high resin volume (`16451` resin TETRA10 vs `9467` tow TETRA10), and it does not yet reproduce the Kok 2022 mixed tow/resin unit-cell stiffness. The geometry metadata records undulation crossings, but the current mesh does not yet partition separate undulation solids with rotated through-thickness tow axes. Those M4 gaps dominate the stiffness miss.

Recommended next fix: implement physical undulation volumes and Kok-style resin/tow mixed unit cells in `kok_geom`, then rerun this same runner. The solver-side recipe itself is verified here: LAW12 + TYPE6 + `/INIBRI/ORTHO` starts and completes on the generated TETRA10 mesh.

Verdict: `FAIL`.
