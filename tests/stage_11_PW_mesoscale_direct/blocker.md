# Stage 11 Blocker - AP-PLY stiffness gate

Author: J.C. Vaught

The Kok geometry port generated the 25 mm x 25 mm, four-ply `[0/+45/-45/90]` AP-PLY block and the OpenRadioss decks used `LAW12 + TYPE6/SOL_ORTH` with per-tow `/INIBRI/ORTHO`.

The solver completed status and numerical comparison are:

- E_x_GPa: measured `53.76192024927913` GPa, target `53.3` GPa, error `0.8666421187225736` percent.
- E_y_GPa: measured `59.306453208261246` GPa, target `53.3` GPa, error `11.26914297985225` percent.
- G_xy_GPa: measured `14.383728664547622` GPa, target `20.5` GPa, error `29.83546992903599` percent.

Kok 2022 used VTC401 carbon/epoxy target data, while this project stage uses the canonical IM7/8552 material card. Material delta and the coarsened TETRA10 validation mesh are the leading suspected causes if the miss is near the tolerance band.

Verdict: `FAIL`.
