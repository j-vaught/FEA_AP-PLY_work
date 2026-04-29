# Open-Source FEA Alternatives for High-Velocity Composite Impact

**Author:** J.C. Vaught
**Date:** 2026-04-29
**Goal:** Identify open-source / academically-licensed code paths that can simulate
high-velocity projectile penetration of an AP-PLY pseudo-woven composite laminate,
covering coupon-scale tensile validation through ballistic-scale fragmentation.

The user's hardest requirement is *composite ballistic impact with progressive damage
plus element erosion / fragmentation* in a fully scriptable (no-GUI) workflow.
This document compares fifteen candidates against six capabilities:

1. (a) Tensile coupon (anisotropic small-strain elasticity)
2. (b) Anisotropic composite layup (laminate / lamina-level)
3. (c) Progressive damage (Hashin / Puck / CDM / cohesive)
4. (d) Contact mechanics (penalty / mortar / node-to-surface)
5. (e) Explicit dynamics (central difference / forward Euler)
6. (f) High-velocity impact with element erosion and fragmentation

A "Yes" requires a solver-native implementation; a "Yes (UMAT)" means the user must
write the constitutive law; "No" means the capability is absent or experimental.

---

## 1. CalculiX (ccx)

- **License:** GPL v2.
- **Composites:** Composite shell elements (8-node S8 / 6-node S6) with arbitrary
  number of layers, orthotropic plies, fiber-direction ORIENTATION cards.
  No native woven / pseudo-woven elements.
- **Damage:** No first-class composite damage model in ccx 2.23. Damage must be
  written as a user material in Fortran (`umat_user.f`) inside ccx itself.
  No SDV-driven element-deletion API; community forum confirms element removal is
  not natively supported and requires source patching.
- **Contact:** Node-to-face penalty (small + large sliding), face-to-face penalty
  (small sliding), and recently massless contact for explicit dynamics.
- **Dynamics:** Both implicit (HHT-alpha) and explicit central difference.
  Explicit supports selective mass scaling.
- **Erosion / fragmentation:** *Not natively.* Requires ccx source modification.
- **Scripting / GUI:** Pure ASCII Abaqus-syntax `.inp`. No GUI required;
  `ccx2.21 -i job` runs head-less. Python wrappers exist (e.g. `pycalculix`).
- **Impact case studies:** Limited. The NASA TM 2019 multi-scale paper (Pineda et
  al.) wraps ccx with FEAMAC for multi-scale composites, but does not address
  ballistic impact. No published projectile penetration with ccx.
- **Maintenance:** Active, version 2.23 (early 2025), single-maintainer
  (G. Dhondt, MTU Aero Engines).
- **Best fit for:** Implicit-dynamic and quasi-static composite stress recovery,
  tensile coupons, low-velocity drop-weight at small scale. **Not** high-velocity
  ballistic impact.

## 2. Code_Aster

- **License:** GPL v3 (until v17) / LGPL portions.
- **Composites:** `DEFI_COMPOSITE` defines arbitrary layered shells with
  orthotropic plies; post-process Hashin / Tsai-Wu via field operators.
- **Damage:** Many CDM models (`ENDO_ISOT_BETON`, `ENDO_ORTH_BETON`, gradient
  damage, GTN, Lemaître-Chaboche). No production-grade Hashin/Puck damage law for
  CFRP; Hashin appears only as a *post-processing failure index*.
- **Contact:** Mature Lagrangian, augmented-Lagrangian, and continuous mortar
  contact. THM coupling.
- **Dynamics:** `DYNA_NON_LINE` operator supports both implicit (HHT, Newmark)
  and explicit central-difference schemes. Explicit is officially documented but
  rarely used at high velocity in published literature.
- **Erosion / fragmentation:** Element death not standard. The XFEM/cohesive
  modules can grow cracks but are *implicit-quasistatic* in practice.
- **Scripting / GUI:** Python `.comm` command files; `as_run` runs CLI. Salome-Meca
  GUI is optional.
- **Impact case studies:** Coupled Eulerian-Lagrangian explicit fluid-structure
  cases exist on the forum, but no peer-reviewed composite projectile penetration
  published with Code_Aster.
- **Maintenance:** Active, EDF-supported. Latest stable v17 (2024).
- **Best fit for:** Implicit nonlinear composite analysis, low-velocity impact,
  thermo-mechanical. **Not** high-velocity erosion-driven penetration.

## 3. Kratos Multiphysics

- **License:** BSD-4 (most permissive on this list).
- **Composites:** `StructuralMechanicsApplication` `SerialParallelRuleOfMixtures`
  and `ParallelRuleOfMixtures` constitutive laws stack fiber+matrix laws. Layer-by-
  layer composite shell modelling supported. `combination_factors` parameterise
  volumetric participation.
- **Damage:** Multiple isotropic, orthotropic, and traction-only damage laws;
  parabolic hardening, exponential softening, asymmetric tension/compression. Can
  combine with rule-of-mixtures for ply-level CDM.
- **Contact:** Dedicated `ContactStructuralMechanicsApplication` (mortar +
  augmented-Lagrangian), plus `DEMApplication` for sphere/non-sphere DEM.
- **Dynamics:** `MPMExplicitSolver` (central difference) and FEM explicit solver
  exist. The `ParticleMechanicsApplication` (a.k.a. `MPMApplication`) is the most
  relevant component, explicitly designed for "extreme events involving impact,
  penetration, fragmentation, blast, multi-phase interaction, and failure
  evolution" per the application README.
- **Erosion / fragmentation:** MPM is a meshfree particle method, so
  fragmentation arises naturally from particle separation. Native MPM-FEM and
  MPM-DEM coupling.
- **Scripting / GUI:** Pure Python driver via `KratosMultiphysics` PyPI; JSON
  ProjectParameters; head-less. GiD GUI is optional.
- **Impact case studies:** Several PFEM/MPM papers from CIMNE on slope failure
  and granular impact; few peer-reviewed *composite ballistic* applications, but
  the MPM application is actively used for impact research.
- **Maintenance:** Very active, v10.4.0 (Dec 2025), CIMNE-led, large international
  contributor base.
- **Best fit for:** **Strongest open-source candidate for ballistic-scale impact
  with progressive damage in a Python pipeline.** Combines explicit MPM,
  composite rule-of-mixtures damage, and contact in one BSD-licensed tool.

## 4. Akantu (EPFL LSMS)

- **License:** LGPL v3.
- **Composites:** No native composite/laminate constitutive library. Anisotropic
  elasticity is supported; ply-level damage requires user material.
- **Damage:** Continuum damage (Mazars, Marigo, non-local), phase-field brittle
  fracture, and *extrinsic / intrinsic cohesive zone elements with parallel
  dynamic insertion* — the central novel feature of Akantu.
- **Contact:** Dedicated contact-mechanics module (penalty, search algorithms),
  designed for fragment-fragment interaction.
- **Dynamics:** Explicit central-difference is a primary use case; the library
  was designed around dynamic fracture and fragmentation.
- **Erosion / fragmentation:** **Best-in-class open-source fragmentation.**
  Vocialta-Richart-Molinari (2017, IJNME) and Molinari et al. (2007 IJNME) on
  cohesive-element insertion underpin the library. Fragments are tracked
  explicitly through cohesive failure, not erosion.
- **Scripting / GUI:** C++ first; Python bindings since v4. Driven from Python
  scripts or C++ mains. No GUI.
- **Impact case studies:** Tempered-glass fragmentation (Vocialta 2017),
  spall and dynamic-fragmentation studies (LSMS lab). No published *composite
  ballistic plate* paper located, but underlying capability exists.
- **Maintenance:** Active, v5.0 (C++17, 2024), GitLab, EPFL LSMS.
- **Best fit for:** **Fragmentation physics + brittle/quasi-brittle dynamic
  fracture.** A strong candidate if AP-PLY damage is modelled as inter-tow
  cohesive interfaces. Composite ply CDM must be self-implemented.

## 5. MoFEM (University of Glasgow)

- **License:** MIT (core); some extensions LGPL.
- **Composites:** Hierarchical / hp-FEM oriented. Anisotropic elasticity yes;
  no off-the-shelf laminate UD ply damage library.
- **Damage:** Cohesive interface elements available (`CohesiveInterfaceElement`
  documented in API), gradient damage, and phase-field fracture modules.
- **Contact:** Mortar contact via `SimpleContact` module.
- **Dynamics:** Implicit transient available; explicit dynamics not advertised
  as production-ready.
- **Erosion / fragmentation:** No element-erosion API; fracture via cohesive
  / phase-field is implicit.
- **Scripting / GUI:** C++; Python bindings exist but documentation is thin.
  No GUI.
- **Impact case studies:** No composite ballistic-impact case study located.
- **Maintenance:** Active, Glasgow research group, smaller community than Kratos
  or Code_Aster.
- **Best fit for:** Hierarchical hp-adaptive fracture analysis, research-grade.
  **Not** ballistic impact.

## 6. deal.II

- **License:** LGPL v2.1.
- **Composites:** Generic FE library; user must implement laminate logic from
  scratch.
- **Damage:** None native; user-implemented.
- **Contact:** Some user examples; no production contact module.
- **Dynamics:** User implements time integration. Explicit central-difference is
  feasible but DIY.
- **Erosion / fragmentation:** Not natively.
- **Scripting / GUI:** C++; Python bindings via PyDealII.
- **Impact case studies:** None directly.
- **Maintenance:** Very active, v9.6 (2024).
- **Best fit for:** Method development, custom solvers. **Not** a turnkey
  ballistic code.

## 7. FEAP (R.L. Taylor, UC Berkeley)

- **License:** Custom academic (per-user fee for full version; free FEAPpv).
- **Composites:** Layered solid / shell elements with orthotropic plies.
  Composite damage models in research user-elements (notably AceGen-generated
  user elements) but not in the standard distribution.
- **Damage:** Visco-elasticity with damage, elasto-plasticity, user materials.
- **Contact:** Tied / sliding via user elements; not a strong area.
- **Dynamics:** Newmark implicit; FORWard explicit (1st-order ODE) added in 8.6.
  Explicit central-difference for elasto-dynamics is supported but not the
  primary use case.
- **Erosion / fragmentation:** None native.
- **Scripting / GUI:** Plain text command files; no GUI required.
- **Impact case studies:** Many in literature for soil dynamics; sparse for
  composite ballistic.
- **Maintenance:** Slow, source distributed by Berkeley, small community.
- **Best fit for:** Research user-element development, education. **Not**
  high-velocity composite impact out of the box.

## 8. Project Chrono

- **License:** BSD-3.
- **Composites:** Beam/shell/solid FEA modules; shell composite layups exist.
  Not the strong suit.
- **Damage:** None native.
- **Contact:** Excellent multibody contact (penalty + complementarity), but
  not designed for fine-mesh continuum failure.
- **Dynamics:** Real-time multibody + flexible-body dynamics; not classical
  high-rate explicit FEA.
- **Erosion / fragmentation:** Can simulate rigid-body fragmentation but not
  continuum erosion.
- **Scripting / GUI:** C++ + PyChrono; head-less.
- **Impact case studies:** Multibody / vehicle / granular; no composite ballistic.
- **Maintenance:** Active, UW-Madison / UniParma.
- **Best fit for:** Coupled rigid-body / DEM environments around an impact event
  rather than the impact itself.

## 9. Karamelo (Adelaide / Monash)

- **License:** GPL v2.
- **Composites:** No native composite constitutive library; user must add.
- **Damage:** Limited; framework supports adding constitutive models.
- **Contact:** MPM grid-based contact (Total Lagrangian MPM contact).
- **Dynamics:** Explicit MPM (USL/MUSL/CPDI), B-spline shape functions, GPU via
  Kokkos branch.
- **Erosion / fragmentation:** Native through MPM particle separation.
- **Scripting / GUI:** LAMMPS-style text input script; C++ source. No Python API.
- **Impact case studies:** General extreme-deformation demos in De Vaucorbeil et
  al. (CPM 2021); no composite ballistic publications located.
- **Maintenance:** Moderate — 770+ commits, no formal releases, ~94 stars.
- **Best fit for:** MPM research where users will write their own constitutive
  laws. Smaller community than Kratos MPM.

## 10. Uintah / MPM (University of Utah)

- **License:** MIT.
- **Composites:** Anisotropic elastic / plastic models; some ply-level
  implementations in research forks. AnisoMPM (Wolper et al. 2020) is a closely
  related anisotropic-damage MPM, but lives in `ziran2020`, not Uintah core.
- **Damage:** Continuum damage, Johnson-Cook with damage, user-defined.
- **Contact:** MPM grid contact (frictional, no-slip, separable) — well
  established at Utah.
- **Dynamics:** Explicit MPM (the canonical Sulsky-Schreyer MPM lineage).
- **Erosion / fragmentation:** Native; particle separation gives fragmentation
  without remeshing.
- **Scripting / GUI:** XML problem-spec (`.ups`); C++ source; head-less; no
  Python API.
- **Impact case studies:** Banerjee, Guilkey, et al. on shaped-charge jet
  penetration of aluminum (peer-reviewed); ample literature on MPM penetration
  in metals. Composite plate ballistic impact in Uintah is rarer.
- **Maintenance:** Active, large code base, SCI Institute (Utah).
- **Best fit for:** Large-scale MPM penetration research where user is willing
  to write composite constitutive models in C++.

## 11. Anura3D

- **License:** Open-source (custom — released by Anura3D Research Community
  April 2021, on GitHub, GPL-style).
- **Composites:** Geomechanics-focused. Not designed for fiber composites.
- **Damage:** Soil constitutive models, UMAT-style interface for user models.
- **Contact:** MPM grid contact for soil-structure interaction.
- **Dynamics:** Explicit MPM.
- **Erosion / fragmentation:** Through MPM particle separation; not really
  intended for projectile fragmentation.
- **Scripting / GUI:** GiD-based pre/post (optional GUI). Fortran source, no
  Python API.
- **Impact case studies:** Saturated/unsaturated soil penetration, landslide-
  structure interaction; no composite ballistic.
- **Maintenance:** Active community, CIMNE workshop 2025.
- **Best fit for:** Soil-water-structure MPM. **Not relevant** to AP-PLY
  ballistic impact.

## 12. SU2

Skipped per user instruction (CFD-only).

## 13. GetFEM

- **License:** LGPL v3.
- **Composites:** Composite finite-element basis options
  (`FEM_STRUCTURED_COMPOSITE`, hierarchical composite); not laminate-composite
  per se. Anisotropic elasticity yes.
- **Damage:** Generic-assembly framework allows user-defined CDM; nothing
  off-the-shelf for fiber composites.
- **Contact:** Mature contact / friction module (Alart-Curnier, augmented
  Lagrangian).
- **Dynamics:** Generic explicit RHS support; central-difference time stepping
  is supported but rarely used at high velocity.
- **Erosion / fragmentation:** None native.
- **Scripting / GUI:** Python, MATLAB, Octave, Scilab; no GUI.
- **Impact case studies:** Contact mechanics / forming; no composite ballistic.
- **Maintenance:** Active, Yves Renard (INSA Lyon).
- **Best fit for:** Implicit contact-mechanics research with custom
  constitutive laws.

## 14. NGSolve / Netgen

- **License:** LGPL.
- **Composites:** Generic FE; high-order. No off-the-shelf laminate / ply damage.
- **Damage:** Research papers on finite-strain damage with NGSolve exist
  (Bartels-Mielke 2020) but are user-implemented.
- **Contact:** Limited. No production contact module.
- **Dynamics:** Implicit transient; explicit central-difference can be assembled
  by user but not turnkey.
- **Erosion / fragmentation:** None.
- **Scripting / GUI:** Python-first (NGS-Py); head-less.
- **Impact case studies:** None for composite ballistic.
- **Maintenance:** Very active, TU Wien.
- **Best fit for:** High-order spectral / hp-elasticity research, NOT impact.

## 15. OOFEM (CTU Prague)

- **License:** LGPL v2.1.
- **Composites:** Layered shells, orthotropic plies, fiber-reinforced concrete
  models. Anisotropic damage models (continuum-damage with ortho-tropic
  evolution) — more developed than CalculiX.
- **Damage:** Many CDM models (Mazars, Mazars-Pijaudier-Cabot, Microplane M4/M7,
  isotropic and orthotropic damage). Tsai-Wu-style failure criteria can be added
  via material registration.
- **Contact:** Basic Lagrange-multiplier and penalty; not as mature as Kratos.
- **Dynamics:** Implicit Newmark and explicit central-difference engines exist
  in the SM (structural-mechanics) module.
- **Erosion / fragmentation:** No element-erosion API in mainline; cohesive
  surfaces possible.
- **Scripting / GUI:** Python bindings (`oofem` on PyPI, v3.0 Dec 2025) allow
  custom constitutive models, elements, time functions in Python. Plain `.in`
  text input also supported. No GUI required.
- **Impact case studies:** Quasi-static / cyclic CFRP, concrete impact loading;
  no peer-reviewed composite ballistic-impact paper located.
- **Maintenance:** Active, version 3.0 released Dec 2025, two universities
  (CTU + Chalmers).
- **Best fit for:** Coupon-scale composite progressive damage, low-velocity
  impact. Excellent Python extensibility. Mid-range candidate for ballistic.

---

## Comparison Matrix

Legend: `Y` = native, production-ready; `~` = possible via user material /
extension; `N` = absent. "Script" = Python or text input only, head-less.
"Free?" = no license fee.

| Code | Lic. | (a) Coupon | (b) Layup | (c) Damage | (d) Contact | (e) Explicit | (f) Impact / Erosion | Script | GUI? | Impact lit. |
|---|---|---|---|---|---|---|---|---|---|---|
| CalculiX | GPL | Y | Y | ~ (UMAT) | Y | Y | ~ (patch) | Y (.inp) | none | weak |
| Code_Aster | GPL | Y | Y | ~ (post) | Y | Y | N | Y (.comm) | optional | weak |
| **Kratos** | **BSD** | **Y** | **Y** | **Y** | **Y** | **Y (FEM+MPM)** | **Y (MPM)** | **Y (Python)** | **none** | **moderate** |
| **Akantu** | **LGPL** | **Y** | **~** | **Y (CZM/PF)** | **Y** | **Y** | **Y (cohesive)** | **Y (Py/C++)** | **none** | **fragmentation** |
| MoFEM | MIT | Y | ~ | ~ | Y | ~ | N | Y (C++/Py) | none | none |
| deal.II | LGPL | Y | DIY | DIY | DIY | DIY | DIY | Y (C++/Py) | none | none |
| FEAP | acad. | Y | Y | ~ | ~ | Y | N | Y (text) | none | sparse |
| Chrono | BSD | Y | ~ | N | Y (MBD) | Y | ~ (MBD) | Y (PyChrono) | none | MBD-only |
| Karamelo | GPL | ~ | DIY | DIY | Y (MPM) | Y | Y (MPM) | Y (text) | none | extreme deform. |
| **Uintah** | **MIT** | **Y** | **~** | **Y (JC, CDM)** | **Y (MPM)** | **Y (MPM)** | **Y (MPM)** | **Y (XML)** | **none** | **strong (metals)** |
| Anura3D | GPL-like | Y | N | ~ (soil) | Y (MPM) | Y | Y (MPM) | Y (GiD/text) | optional | soils |
| GetFEM | LGPL | Y | ~ | DIY | Y | ~ | N | Y (Python) | none | none |
| NGSolve | LGPL | Y | DIY | DIY | DIY | DIY | N | Y (Python) | none | none |
| **OOFEM** | **LGPL** | **Y** | **Y** | **Y (multi)** | **Y** | **Y** | **~** | **Y (Py/text)** | **none** | **weak** |
| FEniCSx | LGPL | Y | ~ | DIY | DIY | DIY (UL) | DIY | Y (Python) | none | none for impact |

The four codes in **bold** are the only ones that simultaneously offer
explicit dynamics, contact, and a fragmentation/erosion mechanism; FEniCSx is
listed at bottom as the user's reference baseline.

---

## Honest Assessment of the Hard Requirement

No single open-source code is a drop-in replacement for LS-DYNA / Abaqus-Explicit
+ MAT_162 (composite damage with element erosion) for high-velocity penetration of
fiber composites. The closest candidates are:

1. **Kratos `ParticleMechanicsApplication` (MPM)** — explicit, BSD, Python-driven,
   ParallelRuleOfMixtures composite damage, fragmentation by particle separation.
   *Realistic path with substantial user effort to validate the AP-PLY ply law.*

2. **Akantu** — explicit, LGPL, world-class cohesive-element fragmentation,
   Python bindings. *Realistic path if AP-PLY failure is recast as cohesive
   inter-tow / inter-ply interfaces.* User must add ply CDM.

3. **Uintah / MPM** — explicit, MIT, strong penetration literature in metals;
   composite ply CDM must be implemented in C++. *Most heavyweight, best
   for HPC.*

4. **OOFEM + cohesive surfaces** — best off-the-shelf composite damage library
   with Python bindings, but explicit-impact track record is limited.

FEniCSx alone cannot do this. There is *no* turnkey FEniCSx contact + explicit
+ erosion stack as of v0.9 (2025). Custom updated-Lagrangian explicit dynamics
is feasible (Bouclier et al. 2025) but high-velocity contact and erosion are
research projects, not features.

---

## Recommended Hybrid Pipeline

Given the user's preference for FEniCSx and the scriptability requirement, the
recommended split is:

### Stage 1: Coupon and Lamina Validation (FEniCSx)

- Build the AP-PLY mesoscale RVE in Gmsh (`gmsh -python`).
- Solve quasi-static tension / shear / off-axis tests in DOLFINx with a custom
  Hashin-style continuum-damage UMAT-equivalent (via `mgis.fenics` /
  `dolfinx_external_operator` / `MFront`).
- Calibrate ply elastic, strength, and fracture-energy parameters against
  experimental coupon data.
- Output a `material_card.json` with calibrated parameters.

### Stage 2: Low-Velocity Drop-Weight (FEniCSx OR OOFEM)

- For drop-weight at <10 m/s, OOFEM's implicit Newmark + ortho-tropic damage +
  cohesive interfaces is adequate; driveable from Python.
- Alternatively, FEniCSx-explicit (updated-Lagrangian central difference) with
  penalty contact, following the 2025 Bouclier-style implementation. Validate
  against the coupon-calibrated material card.

### Stage 3: High-Velocity Ballistic Impact (Kratos MPM, primary; Akantu, alt.)

- **Primary: Kratos `MPMApplication`.** The composite is meshed as MPM particles
  with the calibrated ply law translated into a SerialParallelRuleOfMixtures
  constitutive law. Projectile is meshed with the same MPM or as a rigid FEM
  contact body. Explicit central difference; particle separation provides
  fragmentation. All driven from Python; ProjectParameters.json + custom
  constitutive law in C++ (one-time investment).
- **Alternative: Akantu.** Use cohesive-element insertion at every ply
  interface and at every tow boundary (AP-PLY pseudo-weave). Bulk plies use a
  Mazars-orthotropic or user damage law. Explicit central-difference, parallel.
  Python or C++ driver; output to Paraview.
- **Fallback: Uintah.** If MPI scaling and HPC throughput dominate (cluster-class
  runs), Uintah's MPM is more battle-tested for penetration than Kratos's.

### Stage 4: Post-Processing (FEniCSx + Typst/CeTZ)

- Read VTK / XDMF output from Stage 3 with `meshio` / `pyvista`.
- Extract penetration depth vs. time, residual velocity, energy absorption,
  fragment count.
- Plot via Typst/CeTZ per project visualization standards.

---

## Bottom Line

- **One-code answer:** **Kratos Multiphysics**, BSD, Python-driven, with the
  MPMApplication for ballistic and the StructuralMechanicsApplication for the
  coupon scale. It is the only candidate that natively covers all six
  capabilities under a permissive license and a no-GUI Python pipeline.

- **Two-code answer (preferred):** **FEniCSx for coupon + lamina calibration,
  Kratos MPM for ballistic.** This keeps the user inside FEniCSx where they
  prefer, while delegating the explicit-impact problem to a code that was
  built for it.

- **Fragmentation-physics answer:** **FEniCSx + Akantu.** If fragment
  size/velocity distributions are the deliverable rather than penetration
  depth, Akantu's cohesive insertion is the most rigorous open-source method.

The HPC answer is Uintah-MPM but with the largest custom constitutive
development effort.

---

## References

See `impact_refs.bib` for full BibTeX entries.

- Vocialta, Richart, Molinari (2017) — Akantu cohesive fragmentation.
- Banerjee et al. (2009) — Uintah-MPM shaped-charge penetration of Al.
- Wolper et al. (2020) — AnisoMPM anisotropic damage.
- De Vaucorbeil et al. (2021) — Karamelo MPM package.
- Patzák & Bittnar (2001) — OOFEM design.
- Pineda et al. (2019) — CalculiX + FEAMAC multi-scale composites (NASA).
- Bouclier et al. (2025) — Custom constitutive in FEniCSx.
- Lau, Belnoue, Hallett (2023) — AP-PLY automatic ply detection.
- Belnoue, Pinho et al. — AP-PLY delamination prediction.
