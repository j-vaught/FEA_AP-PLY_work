# Stage 05 Specification: Notched Dogbone with Isotropic Ductile Damage (Lemaitre / Mazars-class CDM Demonstrator)

Author: J.C. Vaught
Date: 2026-04-29
Stage tier: E (Validation, capability demonstrator with documented limitation)
OpenRadioss audit verdict: MARGINAL (`references/openradioss_endtoend_audit.md`, row 5)

This specification operationalizes Stage 5 of the master plan (`plan/master_plan.md` table in section 3): introduce continuum damage softening on a metallic dogbone, prove load-displacement mesh objectivity through the elastic plus plastic plus damage-onset regions, and honestly document the post-peak local-CDM mesh sensitivity that the OpenRadioss native cards inherit from the absence of nonlocal regularization.

---

## 1. Goal and capability isolation

The capability isolated by this stage is *strain-softening continuum damage on a solid mesh, with mesh objectivity demonstrated through damage onset*. Earlier stages cover linear elasticity (Stage 1), geometric nonlinearity (Stage 2), J2 plasticity in a smooth dogbone (Stage 3), and stress-concentration mesh refinement (Stage 4). Stage 5 adds exactly one new physics ingredient. The stiffness reduction is driven by an internal scalar damage variable $D \in [0, 1]$ coupled to the equivalent plastic strain $\bar{\varepsilon}^{p}$, which is the textbook Lemaitre 1985 isotropic ductile damage formulation (Lemaitre 1985, *Journal of Engineering Materials and Technology*).

Two things are explicitly *not* isolated here. Anisotropic damage (fiber and matrix damage) is reserved for Stages 6 and 13. Cohesive interfacial separation is reserved for Stage 12. This stage stays inside an isotropic continuum on a single bulk metal.

---

## 2. Standard, references, and the LAW22 versus Mazars substitution

There is no dedicated ASTM standard for the softening branch. Verification is against the canonical CDM benchmarks. Lemaitre 1985 gives the underlying continuum damage mechanics (CDM) theory for ductile metals. Mazars 1986 (*Engineering Fracture Mechanics*) gives the original isotropic damage formulation, originally posed for concrete and conceptually identical at the kinematic level (an effective-stress hypothesis with a scalar $D$ that grows with an equivalent strain measure). Pijaudier-Cabot and Mazars 1989 (*Journal of Engineering Mechanics*, formalized in Pijaudier-Cabot and Mazars 2001 *Handbook of Materials Behavior Models*) introduced the nonlocal-integral and gradient-enhanced regularization that closes the post-peak mesh-pathology of any local damage model. The audit row for Stage 5 (`openradioss_endtoend_audit.md` row 5) records that *OpenRadioss does not ship a Mazars card by name and DOCUMENTATION NOT LOCATED for a true nonlocal or gradient-regularized Mazars*; the closest native substitute on solid bricks is `/MAT/LAW22` (ductile damage with strain criterion) or `/MAT/LAW23` (visco-damage variant). LAW22 is therefore the substrate for this stage.

The substitution is acceptable as an isotropic damage *capability demonstrator* because the underlying CDM theory of LAW22 is the same effective-stress, scalar-$D$ formulation as Lemaitre and Mazars. Its limitation relative to Mazars-Pijaudier-Cabot is the absence of nonlocal regularization, which is why the pass criterion below is split into a pre-peak region (where the result must be mesh-objective) and a post-peak region (where mesh divergence is expected and accepted, with a written diagnostic). This is consistent with the master-plan §10 risk register entry on regularization length and with the audit's MARGINAL verdict.

Citations.

- Lemaitre, J. (1985), "A continuous damage mechanics model for ductile fracture", *J. Eng. Mater. Tech.* 107(1), 83-89. (`Lemaitre1985` in `references/test_progression_refs.bib`.)
- Mazars, J. (1986), "A description of micro- and macroscale damage of concrete structures", *Eng. Fract. Mech.* 25(5-6), 729-737. (`Mazars1986`.)
- Pijaudier-Cabot, G., and Mazars, J. (2001), "Damage models for concrete", in *Handbook of Materials Behavior Models*, Lemaitre (ed.), Academic Press, 500-512. (`PijaudierCabotMazars2001`; original 1989 nonlocal paper cited inside.)
- Vassaux, M., Oliver-Leblond, C., Richard, B., and Ragueneau, F. (2022), "A modified Mazars damage model with energy regularization", *Eng. Fract. Mech.* 261, 108213. (`Vassaux2022ModifiedMazars`.)
- Altair Engineering (2025), "Ductile damage model: LAW22 and LAW23 for solids and shells". (`AltairRadiossDuctileDamage`.)
- Altair Engineering (2025), "/PROP/TYPE14 (SOLID) general solid property reference". (`AltairRadiossPropType14`.)
- Altair Engineering (2025), "Implicit Analysis Activation in Radioss / OpenRadioss". (`AltairRadiossImplActivation`.)
- de Souza Neto, E.A., Peric, D., and Owen, D.R.J. (2008), *Computational Methods for Plasticity*, Wiley, Ch. 12-13 (cited as the canonical 1D Lemaitre bar verification problem in `test_progression_literature.md` Stage 5).

---

## 3. Geometry

The geometry is the ASTM E8 sub-size rectangular dogbone outline of Stage 3, modified by adding a symmetric pair of semicircular notches in the gauge section to localize damage at a known location. Localization at the notch root removes the geometric-imperfection ambiguity that the smooth Stage 3 dogbone would otherwise introduce, and it gives the runner a single, well-posed location at which to evaluate damage onset.

All dimensions are in SI millimetres for layout convenience; the deck is converted to consistent SI metres in the runner (length in m, mass in kg, time in s, force in N, stress in Pa). Numerical fields below are exact.

| Symbol | Description | Value |
|---|---|---|
| $L_\text{tot}$ | total specimen length | 200.0 mm |
| $L_g$ | reduced gauge section length (between fillets) | 60.0 mm |
| $w_g$ | reduced section nominal width (before notch) | 12.5 mm |
| $w_n$ | net width at notch root | 8.0 mm |
| $r_n$ | notch radius (semicircular) | 2.25 mm |
| $w_\text{grip}$ | grip width | 25.0 mm |
| $r_f$ | shoulder fillet radius (grip to gauge) | 12.5 mm |
| $t$ | thickness | 6.0 mm |

The notch radius $r_n = (w_g - w_n)/2 = 2.25$ mm gives a smooth semicircular cutout. The net cross section at the notch is $A_n = w_n \cdot t = 8.0 \times 6.0 = 48.0\ \text{mm}^2$. The notch root is centred axially at $x = L_\text{tot}/2$ and located at $y = \pm w_g/2$, mirrored across the centerline. Both notches are nominally identical so that the symmetric stress field puts the damage onset on the midline ligament.

The geometry is built in GMSH via the Python API in the runner. Two-fold symmetry could be exploited to halve the mesh, but the stage runs the full geometry to make the boundary-condition write-up identical to Stage 3 and to keep the mesh-objectivity comparison free of any symmetry-plane artefacts.

---

## 4. Mesh

Solid HEXA8 (`/PROP/TYPE14`). Three meshes are built for the objectivity study.

| Mesh | Target element size $h$ in notch ligament | Approximate element count | Through-thickness divisions |
|---|---|---|---|
| coarse | 1.00 mm | $\sim 12{,}000$ | 6 |
| medium | 0.50 mm | $\sim 90{,}000$ | 12 |
| fine | 0.25 mm | $\sim 700{,}000$ | 24 |

The grip and shoulder regions are biased to coarser elements ($h \approx 2$ mm) to keep the global element count manageable; only the gauge section ligament between notches is refined. Element size in the ligament is the nominal $h$ in the table. Through-thickness divisions are set to keep aspect ratio near unity in the most refined region. The transition from grip to gauge to ligament is built with a transfinite frame on the grip and an unstructured but sized hex mesh in the gauge, recombined to all-hex via the GMSH `Recombine Surface` directive applied to the extrusion source.

Stage 4 (Kirsch, ASTM D5766) already established the mesh-refinement protocol for stress concentration; the same refinement-zone discipline applies here.

---

## 5. Material card (the core of the stage; LAW22 versus Mazars rationale)

Material is a generic ductile structural metal. Either A36 mild steel or DP780 dual-phase steel is acceptable; the stage chooses A36 for the canonical dogbone tradition and to keep the elastic constants cleanly recognizable, with DP780 retained as an optional alternate in the deck template.

Constitutive law family. Linear isotropic elasticity for $\bar{\varepsilon}^p \le 0$, von Mises J2 plasticity with a tabulated isotropic hardening curve, and Lemaitre-style scalar isotropic damage that activates at a damage-onset plastic strain $\bar{\varepsilon}^p_D$ and saturates at a critical plastic strain $\bar{\varepsilon}^p_R$ where the element is deleted. In OpenRadioss this is `/MAT/LAW22` (DAMA), or `/MAT/LAW23` (DAMAGE) if rate sensitivity is desired; this stage selects LAW22 because rate effects are not part of the capability under test.

Numerical values.

| Symbol | Description | Value | Source |
|---|---|---|---|
| $\rho$ | density | $7850\ \mathrm{kg/m^3}$ | A36 nominal |
| $E$ | Young's modulus | $200\ \mathrm{GPa} = 2.00\times 10^{11}\ \mathrm{Pa}$ | A36 nominal |
| $\nu$ | Poisson's ratio | 0.30 | A36 nominal |
| $\sigma_y$ | initial yield stress | $250\ \mathrm{MPa} = 2.50\times 10^{8}\ \mathrm{Pa}$ | A36 nominal |
| hardening | tabulated $\sigma$ vs. $\bar{\varepsilon}^p$ | linear from 250 MPa at $\bar{\varepsilon}^p = 0$ to 450 MPa at $\bar{\varepsilon}^p = 0.20$, plateau to $\bar{\varepsilon}^p = 0.5$ | A36 / DP780 representative |
| $\bar{\varepsilon}^p_D$ | damage onset plastic strain | 0.05 | per stage brief |
| $\bar{\varepsilon}^p_R$ | full damage / element deletion plastic strain | 0.50 | per stage brief |
| $D_{\max}$ | damage saturation value | 0.999 | LAW22 default safe value |

Why LAW22 is the right native substitute for a Mazars-class isotropic damage demonstrator. The Mazars 1986 model and the Lemaitre 1985 model share the same kinematic core: an effective stress $\tilde{\sigma} = \sigma / (1 - D)$ acts on the undamaged matrix, while a scalar $D \in [0, 1]$ grows with a strain-like internal variable. Mazars uses an equivalent positive-strain norm and a tension-compression decomposition appropriate for concrete; Lemaitre uses an equivalent plastic strain and a damage threshold appropriate for ductile metals. OpenRadioss `/MAT/LAW22` implements the Lemaitre-style equivalent-plastic-strain trigger plus a separate strain criterion in tension, compression, and shear, which is exactly the *engineering-metal* member of this CDM family. It is *not* a literal Mazars implementation, and it is *not* a nonlocal model; the audit row records this as a MARGINAL verdict precisely because nonlocal regularization is absent. The substitution is honest because (i) the underlying CDM kinematics are identical, (ii) the parameter set $(\bar{\varepsilon}^p_D, \bar{\varepsilon}^p_R)$ maps one-to-one to a Mazars onset and rupture pair, and (iii) the demonstrated capability (softening branch with mesh objectivity through damage onset) is the same physics question independent of the specific damage-driver functional form.

What LAW22 cannot do, and what we therefore concede in writing. LAW22 is *local* in the sense of de Souza Neto-Peric-Owen 2008 Chapter 12: damage at a Gauss point is driven only by the local equivalent strain at that Gauss point, with no characteristic length built in. Once the post-peak softening regime begins, the deformation localizes into a band whose width scales with the element size, which means the dissipated energy per unit fracture area, and therefore the load-displacement curve, is mesh-dependent in the fully softened regime. This is a known and documented property of any local CDM model and is the reason Pijaudier-Cabot-Mazars 1989 introduced nonlocal regularization. OpenRadioss has no documented nonlocal or gradient-regularization keyword (audit row 5, "DOCUMENTATION NOT LOCATED for a true non-local / gradient-regularized Mazars"). Stage 5's pass criterion in section 8 below honours this fact: the mesh-objectivity tolerance is enforced from the elastic regime through damage onset, and divergence in the fully softened regime is documented as expected, not failed.

DP780 alternate. If the user wishes to demonstrate a higher-strength steel: $E = 200$ GPa, $\nu = 0.30$, $\sigma_y = 500$ MPa, hardening to 780 MPa at $\bar{\varepsilon}^p = 0.10$, $\bar{\varepsilon}^p_D = 0.04$, $\bar{\varepsilon}^p_R = 0.20$. The deck template carries this as a swappable parameter block.

---

## 6. Boundary conditions and loading

The boundary conditions are identical in spirit to Stage 3 (`tests/stage_03_dogbone_tension/spec.md`), so the only physics added on top of Stage 3 is the LAW22 damage variable.

Reference frame and axes. Axis $x$ runs along the specimen length, $y$ across the width, $z$ through the thickness. The origin is at the centroid of the left grip face.

Constraints.

- Left grip face ($x = 0$): all nodes constrained to $u_x = 0$, $u_y = 0$, $u_z = 0$ (encastre).
- Right grip face ($x = L_\text{tot}$): all nodes prescribed $u_x = u_x^{\text{imposed}}(t)$, $u_y$ and $u_z$ free except a single corner node fixed in $y$ and $z$ to suppress rigid-body rotation about $x$.

Loading. Quasi-static stretch. The right grip is displaced from $u_x = 0$ to $u_x = 6.0$ mm over the full simulation, applied through a `/IMPDISP` linear function. Six millimetres at a 60 mm gauge length corresponds to nominal strain $\sim 10\%$ in the gauge, well past the damage-onset plastic strain of $\bar{\varepsilon}^p_D = 0.05$ and into the softening regime, but short of full rupture so the runner can sample the post-peak comparison.

Rate of loading. The capability is a quasi-static one. Two equally documented options exist in OpenRadioss for this stage.

Option A (preferred, primary path). `/IMPL/QSTAT` implicit quasi-static with `/IMPL/SOLVER` MUMPS-linked. Required if the engine binary was built with the implicit branch. Convergence is non-trivial in the post-damage-onset regime because the tangent stiffness loses positive-definiteness; the solver must therefore be configured with arc-length continuation enabled (`/IMPL/NONLIN/SMDISP`).

Option B (fallback, also documented). Explicit dynamic relaxation. The displacement is ramped over a long pseudo-time (1.0 s) with a smooth tanh ramp, mass scaling enabled with a target of stable time step $\Delta t_\text{stab} \ge 1\times 10^{-7}$ s, and `/DAMP` global damping with a Rayleigh coefficient sufficient to suppress kinetic energy below 5% of internal energy throughout the run (the standard quasi-static-explicit acceptance test). The runner detects whether MUMPS implicit is available; if not, it falls back to Option B and writes the choice into the run log. This matches master-plan §10 risk 2 (MUMPS-linked OpenRadioss build) and audit row 1 (implicit MARGINAL verdict).

Output requested. Time history `T01` of the right-grip reaction force in $x$, displacement of the right-grip control node in $x$, internal energy, kinetic energy, hourglass energy. Animation `.anim` at 50 frames over the run, plus `/H3D/ELEM/DAMA` so that damage variable $D$ is recoverable per element via the Vortex-Radioss reader.

---

## 7. OpenRadioss deck skeleton

The `_0000.rad` (starter) deck declares the geometry, materials, properties, boundary conditions, and the implicit-or-explicit time integration scheme. The `_0001.rad` (engine) deck declares the run controls. Cards listed below by name; values templated at runtime by `runner.py`.

Starter `_0000.rad`.

```
/BEGIN
  Stage 05 dogbone with isotropic ductile damage (LAW22)
/UNIT/1
  Mg mm s
/MAT/LAW22/1
  density rho, E, nu, sigma_y, hardening table id, eps_d, eps_r, D_max
/PROP/TYPE14/1
  general solid (HEXA8) property; isotropic frame
/PART/1
  material 1, property 1
/SUBDOMAIN
  one subdomain wrapping the gauge section ligament (so that /H3D/ELEM/DAMA writes only there)
/NODE
  ... (from GMSH .inp via inp2rad)
/BRICK
  ... (HEXA8 elements; element block from inp2rad)
/GRNOD/PART, /GRNOD/BOX
  group definitions for left grip face, right grip face, control node
/BCS/1
  full encastre on left grip face group
/BCS/2
  y, z constrained on right grip face control node only
/IMPDISP/1
  imposed x-displacement on right grip face group, function 1, ramp from 0 to 6 mm
/FUNCT/1
  ramp function (0,0) to (1.0, 1.0)
/IMPL/QSTAT          # if Option A
/IMPL/SOLVER/MUMPS
/IMPL/NONLIN/SMDISP  # arc-length continuation
/DT/BRICK/1          # explicit time-step targets if Option B
/DAMP/1              # Rayleigh damping if Option B
/MASS/SCAL           # if Option B
/TH/PART, /TH/NODE   # time-history requests
/ANIM/ELEM/DAMA, /H3D/ELEM/DAMA
/END
```

Engine `_0001.rad`.

```
/RUN/STAGE05/1
/TFILE/0.001         # time-history sampling
/ANIM/DT/0.02        # animation cadence, 50 frames over 1 s
/STOP                # at t = 1.0 s (reached imposed-disp end)
/END
```

The exact card field syntax is left to the templater; OpenRadioss reads the deck in fixed-format columns, and the templater enforces those columns. The starter is exercised by `starter_linux64_gf -i job_0000.rad`, and the engine by `engine_linux64_gf -i job_0001.rad -nt N` per the audit-confirmed CLI workflow.

Mesh import. The deck is built from a GMSH-produced `.inp` (Abaqus deck) converted via OpenRadioss `inp2rad` (audit row C: PASS). The runner orchestrates GMSH then `inp2rad` then `starter` then `engine`, all from a single Python entry point.

---

## 8. Reference solution and validation success criterion

There is no closed-form reference for the full notched-dogbone load-displacement curve; the reference solution for Stage 5 is the *self-consistency of the three meshes through damage onset*, plus a sanity check on the elastic regime against analytic axial stiffness.

Elastic-regime sanity check. The notched specimen's nominal effective axial stiffness, treating the gauge ligament as a uniaxial bar of length $L_g$ and cross section $A_n$ in series with a stiffer grip of length $L_\text{tot} - L_g$ and cross section $w_\text{grip} \cdot t$, gives an axial compliance prediction good to a few percent. The runner extracts the slope of the load-displacement curve in the first 10% of the loading and compares to this analytic estimate within 10% (loose tolerance because the fillet contributes a non-trivial Saint-Venant correction).

Mesh-objectivity criterion (the stage's primary pass criterion).

Define the load-displacement curves $F_h(u)$ for each of the three meshes $h \in \{1.00, 0.50, 0.25\}\ \text{mm}$. Define damage-onset displacement $u_D$ as the displacement at which the damage variable first exceeds $D = 0.01$ at any element in the medium mesh; this is a stable definition because the pre-damage response is already mesh-objective. Compute the windowed RMS difference

$$\mathrm{RMSE}_{h_1, h_2} = \sqrt{\frac{1}{u_D} \int_{0}^{u_D} \left[ \frac{F_{h_1}(u) - F_{h_2}(u)}{F_\text{peak}} \right]^2 du}$$

normalized by the peak load $F_\text{peak} = \max F_h(u)$ on the medium mesh.

Pass: $\mathrm{RMSE}_{h_1, h_2} \le 0.05$ for both pairs (coarse, medium) and (medium, fine), evaluated on the displacement window $[0, u_D]$.

Documented and accepted divergence: in the displacement window $[u_D, 6\ \text{mm}]$, the three curves are *not* required to coincide. Mesh-dependent post-peak softening is the documented limitation of LAW22 inherited from its absence of nonlocal regularization. The runner reports the post-peak RMS for transparency, prints both numbers, and labels the post-peak figure with "expected mesh-dependent softening, local CDM, no nonlocal regularization (Pijaudier-Cabot-Mazars 1989)".

This split-criterion design is the honest mapping of the master-plan Stage 5 pass criterion ("mesh-objective load-deflection within 5% RMS") onto the audited capability of OpenRadioss: 5% RMS through damage onset is achievable; 5% RMS in the fully softened regime is not, and would require a tool with nonlocal regularization (FEniCSx + MFront per the audit's recommended hybrid for Stage 5, or a userlib_sdk Fortran nonlocal extension).

Reporting. The runner writes a CSV `mesh_objectivity.csv` with columns `(displacement_m, F_coarse_N, F_medium_N, F_fine_N, D_coarse, D_medium, D_fine)`, plus a summary JSON `mesh_objectivity_summary.json` with the two RMSE numbers and the verdict (PASS / FAIL).

---

## 9. Per-stage runner script overview

The runner `runner.py` does the following, in order, for each of the three meshes.

1. Build the GMSH model from `geometry.geo` (or directly via the `gmsh` Python API), parameterized by a single element-size argument $h \in \{1.0, 0.5, 0.25\}\ \text{mm}$.
2. Export an Abaqus `.inp` via `gmsh.write(...)` then convert to OpenRadioss `.rad` via `inp2rad`.
3. Template the LAW22 material card and the boundary condition cards into the starter deck.
4. Detect implicit-MUMPS availability; pick Option A (`/IMPL/QSTAT`) or Option B (explicit dynamic relaxation) and write that decision to the run log.
5. Invoke OpenRadioss `starter` then `engine` inside Lima per the master-plan toolchain.
6. Read the `T01` time history with `vortex-radioss` and extract the right-grip reaction force vs. displacement.
7. Read the `.anim` (or `.h3d`) and extract per-element damage variable $D$ for the gauge ligament; compute $u_D$ as the displacement at which $\max_e D > 0.01$ on the medium mesh.
8. Compute pre-onset and post-onset windowed RMSE per section 8.
9. Write `mesh_objectivity.csv` and `mesh_objectivity_summary.json`; print a PASS or FAIL verdict.
10. (Optional) emit a CeTZ-ready CSV for the load-displacement plot and a CeTZ-ready CSV for damage contour at $u = u_D$, both consumed by the Typst report later.

The runner does not generate matplotlib figures; per global preferences, plotting is exclusively Typst plus CeTZ from CSV inputs.

---

## 10. Known limitations, risks, and what the next stage assumes

Limitations.

1. LAW22 is local. Post-peak mesh dependence is expected and is honestly reported, not failed (section 8). A reviewer asking for a strictly mesh-objective post-peak softening must be referred to a nonlocal or gradient model (Pijaudier-Cabot-Mazars 1989, Vassaux 2022 with energy regularization), which OpenRadioss does not natively offer; the documented hybrid alternative in the audit is FEniCSx plus MFront with a Mazars primitive (Helfer 2015).
2. Notch geometry localizes damage at a known place but it does not provide a closed-form softening curve; verification is mesh-self-consistency, not closed-form match.
3. A36 versus DP780 is a parameter swap, not a re-derivation; if the user replaces the metal, the parameter block changes but the deck does not.
4. Implicit static convergence in the softening regime requires arc-length continuation. If `/IMPL/NONLIN/SMDISP` arc-length is unavailable in the build (audit row 1 mentions "build issues for implicit-only static"), the explicit-dynamic-relaxation fallback is documented and acceptable.
5. The MUMPS dependency for `/IMPL/QSTAT` is the same risk recorded in master-plan §10 risk 2; the runner's automatic fallback to explicit dynamic relaxation neutralizes the operational impact.

Risks (per master-plan §10 nomenclature).

- macOS Lima file-I/O performance on the fine mesh (700k elements) could be slow; the runner allows `--engine-cores` to be passed to OpenRadioss `engine` to mitigate.
- Hourglass energy on HEXA8 with under-integration can become a confounder; the deck enforces $E_\text{hg} / E_\text{int} \le 5\%$ as a sanity gate (printed by the runner, separate from the mesh-objectivity verdict).
- inp2rad is "beta" (audit row C); the runner sanity-checks element count and node count consistency between the GMSH `.inp` and the converted `.rad`.

What the next stage (Stage 6, Tsai-Wu / Hashin / Puck side-by-side) assumes from this stage.

- Solid HEXA8 dogbone meshing pipeline through GMSH plus inp2rad is reliable.
- BCs and time-history extraction infrastructure is reusable; Stage 6 swaps `/MAT/LAW22` for `/MAT/LAW25` plus `/FAIL/HASHIN`, `/FAIL/PUCK`, `/FAIL/TSAIWU` cards on the same outline, without touching geometry, BCs, or runner harness.
- The mesh-objectivity reporting pattern (CSV plus summary JSON plus PASS/FAIL printout) is the template for every later stage.
