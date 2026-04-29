# Stage 06 Blocker - LAW25 + TYPE14 rejected by starter

Author: J.C. Vaught

The required canonical deck uses `/MAT/LAW25` on `/PROP/TYPE14` with solid HEXA8 elements. The installed OpenRadioss starter rejects that material/property pair before the engine can run.

Observed starter finding:

- `ERROR ID : 3047`
- `ERROR IN MATERIAL/PROPERTY COMPATIBILITY`
- `PROPERTY ID 1 OF TYPE 14 IS NOT COMPATIBLE WITH MATERIAL ID 1 OF TYPE 25`

A separate all-solid `/PROP/TYPE6` (`/PROP/SOL_ORTH`) proxy using the same `/MAT/LAW25` card successfully starts and the engine reports solid-element failure for `/FAIL/TSAIWU`, `/FAIL/HASHIN`, and `/FAIL/PUCK`. That proves the failure cards are present in this binary, but it does not satisfy the stage 06 `/PROP/TYPE14` requirement or consolidation item B.

Verdict: `INCONCLUSIVE` due to toolchain/material-property compatibility, not a numerical failure against the WWFE-II principal-axis strengths.
