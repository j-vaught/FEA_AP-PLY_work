# Stage 10 Blocker - UD mesoscale modulus mismatch

Author: J.C. Vaught

The checked-in stage 10 wrapper was not runnable on this host because it expected an absent `UD_mesoscale_hex.inp`, an absent Lima `or_paths.json`, and a VTKHDF converter path. I replaced it with a self-contained OpenRadioss runner that generates a 20x20x4 HEXA8 voxel cell with a 0.600 fiber volume fraction.

The isotropic fiber/matrix mesoscale cell starts, runs, converts to VTK, and post-processes, but two effective moduli miss the 5 percent gate.

- transverse_x: FEM `22.516 GPa`, target `20.250 GPa`, error `11.189%`.
- shear_xz: FEM `0.004 GPa`, target `5.553 GPa`, error `99.934%`.

Verdict: `FAIL`.
