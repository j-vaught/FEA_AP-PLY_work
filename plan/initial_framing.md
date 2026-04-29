# FEA_AP-PLY — initial framing (pre-research) [SUPERSEDED by master_plan.md]

> Kept for the record. The hybrid FEniCSx+OpenRadioss recommendation in this document was rejected by the user (single-tool, all-solid constraints). The committed direction is single-tool OpenRadioss with reframed stages 10-11 — see `master_plan.md`.



This document captures my honest framing of the problem *before* the five parallel research agents return. The master plan in `master_plan.md` will be written *after* the research and may overrule anything below.

## What you asked for

A pure open-source, no-GUI, code-only pipeline that can simulate, in increasing order of difficulty,

1. tensile of a single composite tow,
2. tensile of a full layup with epoxy,
3. impact of composites,
4. how a pseudo-woven preform deforms to fit a tool,
5. and ultimately a high-velocity projectile penetrating a fully cured UofSC pseudo-woven panel.

Strong preference for FEniCSx alone; willingness to fall back to Code_Aster, CalculiX, Kratos, Akantu, or MPM only if FEniCSx genuinely cannot do something.

## My honest a-priori position (will be checked against research)

**FEniCSx is excellent for stages 1-3 of the test progression** (linear elastic, anisotropic ply, laminate, mesoscale RVE). It is a finite-element library with first-class Python, excellent anisotropic constitutive support via UFL, periodic BCs via `dolfinx_mpc`, and clean PyVista / XDMF output that satisfies the no-GUI constraint.

**FEniCSx alone is a poor fit for the final ballistic stage.** High-velocity projectile penetration needs explicit time integration with very small time steps, robust frictional contact between many surfaces, large element distortions, and either element erosion or a particle method (MPM/SPH) to handle fragmentation. None of these are first-class in dolfinx today. Forcing FEniCSx to do this would mean writing a whole explicit dynamics + contact + erosion stack from scratch — months of work that duplicates what CalculiX, Kratos, or Akantu already do.

The honest answer is therefore likely a **hybrid**: FEniCSx for stages 1-12 (everything quasi-static or low-rate), and a second tool for the impact stages. The research wave will confirm whether this is right.

## Test progression — initial sketch (will be expanded by per-test agents)

This is the staged sequence I'd want to verify, each isolating one capability:

1. Linear-elastic beam (3-pt and 4-pt bending) — closed-form check.
2. Cantilever with geometric nonlinearity — large displacement.
3. Isotropic dogbone (ASTM E8) — uniaxial stress baseline.
4. Open-hole tension — stress concentration.
5. Dogbone with isotropic damage (Mazars-like) — softening branch.
6. Dogbone with multiple failure criteria (Tsai-Wu, Hashin, Puck) on the same coupon — see how each predicts failure.
7. Single-tow UD tension (ASTM D3039), on-axis and off-axis — anisotropy check.
8. Ply with rotation — verify stiffness transformation.
9. Cross-ply and quasi-isotropic laminate — classical lamination theory check.
10. Mesoscale RVE of a UD ply (fiber + matrix) with periodic BCs — homogenized stiffness vs. micromechanics.
11. Mesoscale RVE of a pseudo-woven cell (tow + matrix) — homogenized properties.
12. Cohesive-zone delamination tests — DCB (mode I, ASTM D5528) and ENF (mode II, ASTM D7905).
13. Drop-tower impact (ASTM D7136) — low-velocity.
14. Compression after impact (ASTM D7137).
15. High-velocity / ballistic impact — V50 and residual velocity on a flat coupon.
16. Final goal — projectile through a UofSC pseudo-woven panel.

Each stage has a single explicit success criterion (e.g. stiffness within 2% of analytic, load-displacement RMS within 5% of reference) so we know "this stage works" before adding complexity.

## Visualization

Constraint is no GUI. Two compatible options:

- **PyVista headless** in the simulation script — write PNGs / MP4 directly, fully Pythonic, fits the no-GUI constraint cleanly.
- **ParaView's `pvbatch`** — scriptable from Python, reads XDMF/HDF5 from FEniCSx, can produce publication-grade images and animations.

Plots and figure overlays must be Typst + CeTZ per global preferences (export numerical data from FEniCSx, draw the plot in Typst).

## What the research wave will resolve

- Whether `dolfinx_mpc` / `dolfiny` / `MFront-MGIS` together can carry stages 1-12 cleanly.
- Whether any open-source code can do composite ballistic impact end-to-end (Akantu and Kratos are the most likely candidates).
- Exactly what the UofSC pseudo-woven and Delft AP-PLY architectures look like geometrically — needed before we can mesh them.
- Whether there is open data (micro-CT scans, mechanical raw data, ply layouts) we can use as ground truth.

## Status

Five opus research agents are running in parallel. When they return I will write `master_plan.md` and *then* spawn one further opus agent per test stage to flesh out the specification. No code will be written before the master plan is approved.
