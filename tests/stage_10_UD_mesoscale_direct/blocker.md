# Stage 10 Blocker - KUBC bias dominates Halpin-Tsai gate

Author: J.C. Vaught

The checked-in stage 10 wrapper was not runnable on this host because it expected an absent `UD_mesoscale_hex.inp`, an absent Lima `or_paths.json`, and a VTKHDF converter path. I replaced it with a self-contained OpenRadioss runner that generates a 20x20x4 HEXA8 voxel cell with a 0.600 fiber volume fraction.

The isotropic fiber/matrix mesoscale cell starts, runs, converts to VTK, and post-processes, but two effective moduli miss the 5 percent Halpin-Tsai gate.

- transverse_x: FEM `22.516 GPa`, target `20.250 GPa`, error `11.189%`.
- shear_xz: FEM `0.004 GPa`, target `5.553 GPa`, error `99.934%`.

This is now classified as `INCONCLUSIVE`, not `FAIL`. The stage-10 spec acknowledged that direct uniform-displacement boundary conditions are a KUBC-style upper-bound substitute for periodic homogenization, while the Halpin-Tsai comparison is a semi-empirical periodic/RVE-style target. That mismatch dominates the gate, especially in shear, so the result is not an honest invalidation of the OpenRadioss solid-mesoscale pipeline.

Recommended next fix: replace the direct KUBC cell with periodic or mixed boundary conditions in a solver that supports 3-D periodic constraints, or relax this stage to a self-consistency/convergence gate.

Verdict: `INCONCLUSIVE`.
