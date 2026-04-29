# Stage 08 Blocker - LAW25 + TYPE14 ply-rotation sweep does not start

Author: J.C. Vaught

The canonical stage 08 sweep requires `/MAT/LAW25` on `/PROP/TYPE14` solid HEXA8 coupons with `/SKEW/FIX` frames at seven ply angles. Every canonical deck is rejected by the installed OpenRadioss starter with `ERROR ID : 3047` material/property compatibility before the engine can run.

The runner also generated `/PROP/TYPE6` (`/PROP/SOL_ORTH`) proxy decks at the same seven angles. Those proxy starters succeed, which isolates the issue to the required TYPE14 path rather than to the material card or skew definitions.

Verdict: `INCONCLUSIVE` due to toolchain/material-property compatibility, not a numerical off-axis stiffness failure.
