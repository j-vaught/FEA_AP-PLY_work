# Stage 11 Blocker - AP-PLY stiffness gate

Author: J.C. Vaught

The Kok geometry port generated the 25 mm x 25 mm, four-ply `[0/+45/-45/90]` AP-PLY block and the OpenRadioss decks used `LAW12 + TYPE6/SOL_ORTH` with per-tow `/INIBRI/ORTHO`.

The solver completed status and numerical comparison are:

- E_x_GPa: measured `53.9154639276` GPa, target `53.3` GPa, error `1.154716562101321` percent.
- E_y_GPa: measured `58.435006723866664` GPa, target `53.3` GPa, error `9.634158956597874` percent.
- G_xy_GPa: measured `14.428740611936254` GPa, target `20.5` GPa, error `29.61589945396949` percent.

Kok 2022 used VTC401 carbon/epoxy target data, while this project stage uses the canonical IM7/8552 material card. Material delta and the coarsened TETRA10 validation mesh are the leading suspected causes if the miss is near the tolerance band.

Verdict: `FAIL`.
