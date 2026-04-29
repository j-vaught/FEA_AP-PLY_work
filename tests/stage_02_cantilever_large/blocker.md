# Stage 02 Blocker: Nonlinear Implicit Requires MUMPS

**Author.** J.C. Vaught
**Date.** 2026-04-29

## Verdict

`INCONCLUSIVE`

## What was attempted

The Stage 02 runner now builds the solid HEXA8 cantilever deck directly, writes OpenRadioss starter/engine files, and probes the nonlinear implicit path before falling back to explicit dynamic relaxation.

The implicit starter completes, but the engine terminates immediately with:

```text
Fatal error: MUMPS required
```

This confirms that the installed OpenRadioss prebuilt at `/mnt/storage/j-vaught/openradioss/OpenRadioss` is not MUMPS-linked.

The automatic explicit fallback was probed on the coarse mesh at `alpha=1.0` with `/KEREL`, `/DAMP`, and `/DT/NODA/CST`. It runs end-to-end, emits VTK, and is recorded in `results/timeseries.csv`, but the one-second mass-scaled run remains far from quasi-static equilibrium:

```text
dy/L FEM = 3.400685e-05
dy/L ref = 3.017208e-01
dx/L FEM = 1.033401e-08
dx/L ref = 5.643324e-02
```

Because this is not the gated baseline mesh and because the intended implicit path is unavailable, this is a toolchain inconclusive result rather than a numerical verification failure.

## Required resolution

Either provide a MUMPS-linked OpenRadioss build for `/IMPL/NONLIN`, or allocate a much longer explicit dynamic-relaxation run plan with kinetic-energy checks before treating the fallback as a verification result.
