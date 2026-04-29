# Stage 07 Blocker - LAW25 + TYPE14 D3039 coupons do not start

Author: J.C. Vaught

Both canonical D3039 coupons require `/MAT/LAW25` on `/PROP/TYPE14` solid HEXA8 elements. The installed OpenRadioss starter rejects that pairing with `ERROR ID : 3047` material/property compatibility before any modulus extraction can be run.

The runner also generated `/PROP/TYPE6` (`/PROP/SOL_ORTH`) proxy decks with `/SKEW/FIX` ply frames. Those proxy decks start successfully, confirming the material card itself is accepted on a solid orthotropic property, but the proxy does not satisfy the stage 07 `/PROP/TYPE14` requirement.

Verdict: `INCONCLUSIVE` due to toolchain/material-property compatibility, not a numerical D3039 modulus failure.
