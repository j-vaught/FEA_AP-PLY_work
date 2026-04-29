# Stage 02 Blocker: Large-Deflection Implicit Convergence

**Author.** J.C. Vaught
**Date.** 2026-04-29

## Verdict

`INCONCLUSIVE`

## Status

The earlier MUMPS blocker is resolved. OpenRadioss was rebuilt from source with the conda-forge MUMPS MPI libraries, and the installed engine now reaches `DMUMPS 5.8.2` under `/IMPL/NONLIN`.

The refined baseline `80 x 4 x 6` solid mesh with `Isolid=14`, `Ismstr=11`, and `Icpre=1` completes the gated `alpha=1` case:

```text
dx/L FEM = 5.732170e-02, ref = 5.643324e-02, error = 1.574%
dy/L FEM = 3.043490e-01, ref = 3.017208e-01, error = 0.871%
```

The required `alpha=3` and `alpha=5` cases do not yet complete as a valid verification. The blocker is now nonlinear path convergence, not missing MUMPS.

## Failed follow-up paths

- Load control, `/IMPL/DT/2`, force norm (`NITOL=2`) reaches only about `T=0.066` for `alpha=3` before timestep-limit termination.
- Load control, displacement norm (`NITOL=3`) advances farther, to about `T=0.271`, then terminates on timestep limit.
- Riks/arc-length (`/IMPL/DT/3`) can traverse close to the target load, but either exits with timestep-limit failure at the final step before writing a usable animation frame or reverses into negative arc-length steps during a held-load variant.
- `/IMPL/QSTAT` stalls almost immediately near `T=8.6e-4`.
- Displacement-control variants are not acceptable as verification: whole-face `IMPDISP` overconstrains tip rotation, single-node `IMPDISP` pulls a local solid node, and a trial `/RBODY` tip control over-stiffens the end condition.
- `Isolid=24` HEPH is much slower and stalls early for this case.

## Required resolution

Stage 02 needs a robust OpenRadioss nonlinear static path for the `alpha=3` and `alpha=5` dead-load cantilever cases, or a documented explicit/quasi-static fallback with kinetic-energy and mesh-convergence checks. Until then, stages 03-16 remain unattempted because Stage 02 is the geometric-nonlinearity gate.
