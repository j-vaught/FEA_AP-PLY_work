# Stage 09 Blocker - LAW25 + TYPE14 laminate stack does not start

Author: J.C. Vaught

The canonical stage 09 CLT verification requires one `/MAT/LAW25` and one `/PROP/TYPE14` per solid ply. The installed OpenRadioss starter rejects the LAW25 + TYPE14 pairing with `ERROR ID : 3047` material/property compatibility for both the `[0/90]s` and `[0/+45/-45/90]s` stacks.

Analytic CLT A-matrices are written to `results/timeseries.csv`. `/PROP/TYPE6` (`/PROP/SOL_ORTH`) proxy stacks start successfully, so the blocker is the required TYPE14 path rather than the material card or layup bookkeeping.

Verdict: `INCONCLUSIVE` due to toolchain/material-property compatibility, not a numerical CLT failure.
