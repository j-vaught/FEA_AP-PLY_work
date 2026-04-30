# Stage 11 Blocker - AP-PLY stiffness gate

Author: J.C. Vaught

The Kok geometry port generated the 25 mm x 25 mm, four-ply `[0/+45/-45/90]` AP-PLY block and the OpenRadioss decks used `LAW12 + TYPE6/SOL_ORTH` with per-tow `/INIBRI/ORTHO`.

The solver completed status and numerical comparison are:

- E_x_GPa: measured `43.53917679364661` GPa, target `53.3` GPa, error `18.312989130118925` percent.
- E_y_GPa: measured `46.864719528066665` GPa, target `53.3` GPa, error `12.073696945465914` percent.
- G_xy_GPa: measured `11.362530276086744` GPa, target `20.5` GPa, error `44.573023043479296` percent.

Kok 2022 used VTC401 carbon/epoxy target data, while this project stage uses the canonical IM7/8552 material card. Material delta and the coarsened TETRA10 validation mesh are the leading suspected causes if the miss is near the tolerance band.

Verdict: `FAIL`.
