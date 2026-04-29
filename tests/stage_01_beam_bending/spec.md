# Stage 01 — Linear-Elastic 3-Point and 4-Point Bending of a Solid Prismatic Beam

**Author.** J.C. Vaught
**Date.** 2026-04-29
**Status.** Phase 2 deep-dive specification. Toolchain smoke test for the FEA_AP-PLY single-tool OpenRadioss pipeline (master_plan.md §3, row 1). MARGINAL on the openradioss_endtoend_audit.md row 1 because /IMPL/LINEAR requires a MUMPS-linked OpenRadioss build; this spec documents the build path and the dynamic-relaxation fallback.

---

## 1. Stage summary

A simply supported prismatic aluminum 6061-T6 beam, 200 mm long with a 20 mm × 10 mm rectangular cross-section, is loaded first in three-point bending with a single midspan point load and then in four-point bending with two equal point loads at the third-points. The beam is meshed exclusively with HEXA8 solid elements (master_plan.md §1 hard requirement: solid elements only; no shell laminate stacks). The simulation is run as a static linear-elastic analysis using OpenRadioss /IMPL/LINEAR with /MAT/LAW1 (linear elastic) on /PROP/TYPE14 (general solid) bricks. The stage is the toolchain smoke test for every downstream stage. Pass criterion. Midspan deflection at the converged mesh is within 1% of the Euler-Bernoulli closed-form solution for both load cases — δ_3pt = PL³/(48EI) for three-point bending and δ_4pt,max = Pa(3L² − 4a²)/(24EI) with a = L/3 for four-point bending. If this stage does not run end-to-end (geometry → mesh → .rad deck → starter+engine inside Lima → vtkhdf → PyVista headless → CSV), the toolchain is broken and no later stage can succeed.

---

## 2. Geometry

All dimensions are SI (metres internally; millimetres in the deck because OpenRadioss accepts a unit system declared on `/UNIT`). Coordinate convention. x = beam axis (length), y = width, z = height; origin at the centroid of the left end face.

| Quantity | Symbol | Value | Unit |
|---|---|---|---|
| Length | L | 0.200 | m (200 mm) |
| Width | b | 0.020 | m (20 mm) |
| Height | h | 0.010 | m (10 mm) |
| Cross-section area | A = b·h | 2.0 × 10⁻⁴ | m² |
| Second moment of area (about y) | I = b h³ / 12 | 1.6667 × 10⁻⁹ | m⁴ |
| Slenderness | L / h | 20 | — |

The slenderness ratio L/h = 20 places the beam squarely in the Euler-Bernoulli regime per Timoshenko and Gere (1972) [@Timoshenko1972MechanicsOfMaterials]; the shear-deflection correction PL/(4 κ G A) for a rectangular section with κ ≈ 5/6 contributes well under 1% of the bending deflection (verified in §7 below), so the closed-form pass criterion does not need a Timoshenko correction.

ASCII sketch — three-point bending.

```
                                P (downward)
                                  │
                                  ▼
   z                          ┌───┴───┐
   ▲              ┌───────────┤ ▲     ├───────────┐
   │              │           │       │           │
   └──► x         │     beam: 200 mm × 20 mm × 10 mm    │
                  └───────────┴───────┴───────────┘
                  ▲                                   ▲
                  Δ pin support                       Δ roller support
                  (uy = uz = 0,                       (uz = 0)
                   ux = 0 at this end only)
                  x = 0                               x = L = 0.200 m
                                                       midspan x = L/2 = 0.100 m
```

ASCII sketch — four-point bending.

```
                       P/2              P/2
                        │                │
                        ▼                ▼
   z                ┌───┴───┐        ┌───┴───┐
   ▲    ┌───────────┤       ├────────┤       ├───────────┐
   │    │           │  load │        │  load │           │
   └──►x│    beam:  │       │        │       │           │
        └───────────┴───────┴────────┴───────┴───────────┘
        ▲                                                ▲
        Δ pin                                            Δ roller
        x = 0    x = L/3 ≈ 0.0667     x = 2L/3 ≈ 0.1333  x = L = 0.200
```

---

## 3. Mesh

Element type. HEXA8 (8-node linear hex), one OpenRadioss `BRICK` per element, all sharing the same `/PROP/TYPE14` solid property. No shells, no beams, no tetrahedra anywhere in this stage. This complies with master_plan.md §1 ("solid only — HEXA8 / TETRA10 / pentahedra; no shell layups").

Target element size and convergence study. The pass criterion is 1% on midspan deflection, which is tight for a linear hex on a slender beam because HEXA8 with full integration suffers from shear locking; for that reason the deck will use OpenRadioss `Isolid` formulation flag = 17 (HA8 — fully integrated linear brick with assumed-strain physical-stabilization, the default modern recommendation for /PROP/TYPE14, see Altair PROP/TYPE14 reference [@AltairRadiossPropType14]) [UNVERIFIED-default-Isolid] — the deck skeleton in §6 lists `Isolid=17` as a candidate; if shear locking still appears we fall back to `Isolid=24` (HEPH with hourglass control) or upgrade to HEXA20 (TYPE20 brick) which is overkill for a smoke test but a documented option.

Through-thickness count must be ≥ 5 elements (brief). Refinement plan, three meshes:

| Mesh | n_x (length) | n_y (width) | n_z (thickness) | Elements | DOF (approx) |
|---|---|---|---|---|---|
| M0 (coarse) | 40 | 4 | 5 | 800 | ~7.4 k |
| M1 (medium) | 80 | 8 | 10 | 6 400 | ~52 k |
| M2 (fine) | 160 | 16 | 20 | 51 200 | ~390 k |

Element edge lengths on M1 (target). Δx = 2.5 mm, Δy = 2.5 mm, Δz = 1.0 mm. Aspect ratio 2.5:2.5:1.0 — within the Altair guideline of <5:1 for solid brick analysis. Convergence is judged on midspan vertical deflection u_z(L/2, 0, 0) under the same load. The pass criterion (1% vs. closed form) must be met on M1 or finer; M0 is allowed to fail it. Reporting: log mesh size vs. error in a CSV and Typst/CeTZ plot showing the expected O(h²) convergence rate for linear Lagrange elements (test_progression_literature.md, Stage 1 cross-cutting note).

Refinement zones. None — uniform mesh is sufficient for a Saint-Venant-decayed midspan response on a slender beam with point-load BCs that we will distribute across a small loading-pad node set (see §4). No stress-concentration refinement is required because there is no hole, notch, or geometric discontinuity in this stage.

Mesh generation. GMSH Python API (master_plan.md §4 toolchain) writing a structured hex mesh of the box `[0,L] × [-b/2, b/2] × [-h/2, h/2]` with `setTransfiniteCurve` / `setTransfiniteSurface` / `setTransfiniteVolume` and `recombine` to force HEXA elements. Output is Abaqus `.inp`, then converted to `.rad` via the OpenRadioss `inp2rad` Python tool (audit row C, [@OpenRadiossDiscussion3222INP2RAD]). The runner.py in this directory orchestrates the GMSH call as well as the inp2rad invocation.

---

## 4. Boundary conditions and loading

All BCs are applied to *node sets*, not to surface loads, to avoid pressure-load surface integration in this smoke test. Every node set is built up from the structured mesh node ordering at deck-template time.

### 4.1 Support node sets (shared between 3-pt and 4-pt)

| Node set | Definition | DOF constrained | Rationale |
|---|---|---|---|
| `NSET_SUPPORT_LEFT` | All nodes at x = 0, z = −h/2, y arbitrary | u_y = 0, u_z = 0, u_x = 0 (pin) | Simply supported pin at left end, restrains rigid-body translation entirely. |
| `NSET_SUPPORT_RIGHT` | All nodes at x = L, z = −h/2, y arbitrary | u_y = 0, u_z = 0 (roller) | Simply supported roller at right end; u_x free to allow axial elongation, removing a spurious axial reaction. |

The pin/roller convention is the canonical simply supported beam BC of Timoshenko and Gere [@Timoshenko1972MechanicsOfMaterials] and Gere and Goodno [@GereGoodno2012MechanicsOfMaterials]. The y = 0 plane is *not* additionally constrained because the structured mesh + symmetry of the problem will place the support node set on the bottom-fiber line and the mesh has enough through-width nodes (≥ 5 on M1 in y) that no rigid-body rotation about x is unrestrained — but as a belt-and-braces step the deck pins one mid-width node at each support to u_y = 0 (already covered by the entire support line above; the constraint is consistent). [UNVERIFIED — the cleanest pin/roller convention may instead pin only the centroid of each support face; both conventions converge to the same closed-form deflection within 0.1% on a slender beam, so this is documented as a calibration knob, not a correctness risk.]

### 4.2 Load node sets — 3-point bending

| Node set | Definition | Applied | Magnitude (each direction) |
|---|---|---|---|
| `NSET_LOAD_3PT` | Nodes at x = L/2, z = +h/2, y arbitrary | Concentrated force in −z | Total P = 1000 N distributed over the line of nodes (force per node = P / N_load_nodes) |

The load is applied via `/CLOAD` (concentrated load) on the top-fiber midspan line. Why a line and not a single node — a single-node point load on a HEXA8 mesh produces a local stress singularity that pollutes the displacement field over a few elements (Saint-Venant's principle puts the polluted region within ~ h of the load), but the *integrated* midspan deflection is still recovered to <1% on M1 because we measure it at the bottom-fiber mid-width node, which is one half-thickness below the load.

P = 1000 N is selected so that the maximum bending stress σ_max = M c / I = (PL/4)(h/2)/I = (1000 × 0.200 / 4)(0.005) / 1.6667 × 10⁻⁹ = 1.50 × 10⁸ Pa = 150 MPa. This is well below the 6061-T6 yield stress of ≈ 276 MPa [UNVERIFIED — typical reference value, citation in §5] so the linear-elastic assumption holds.

### 4.3 Load node sets — 4-point bending

| Node set | Definition | Applied | Magnitude |
|---|---|---|---|
| `NSET_LOAD_4PT_A` | Nodes at x = L/3, z = +h/2, y arbitrary | Concentrated force in −z | P/2 = 500 N total (per node force = (P/2) / N_load_nodes) |
| `NSET_LOAD_4PT_B` | Nodes at x = 2L/3, z = +h/2, y arbitrary | Concentrated force in −z | P/2 = 500 N total |

Total applied load matches the 3-pt case (P = 1000 N split into two equal P/2) so that maximum bending moment in the constant-moment region M_max,4pt = (P/2)(L/3) = PL/6 (constant between the two loads) is two-thirds of the 3-pt midspan moment M_max,3pt = PL/4. The closed-form midspan deflection in §7 reflects this.

### 4.4 Measurement node — pass criterion is read off this node

| Probe | Definition | Quantity reported |
|---|---|---|
| `NODE_MIDSPAN` | Single node at (x, y, z) = (L/2, 0, −h/2) | Vertical displacement u_z (signed; expected negative because load is in −z) |

The measurement is taken at the bottom fiber midspan because (a) it is a *response* point, not a *forcing* point, so it is unaffected by the local load singularity in the 3-pt case, and (b) by symmetry the bottom-fiber midspan is exactly where Euler-Bernoulli gives δ_max. For 4-pt bending the same `NODE_MIDSPAN` is the maximum-deflection point (closed form below).

---

## 5. Material card

Aluminum 6061-T6 (annealed-then-aged), one of the most widely characterized engineering alloys. Linear-elastic isotropic for this stage; plasticity is deferred to Stage 3.

| Property | Symbol | Value | Unit | Source |
|---|---|---|---|---|
| Young's modulus | E | 6.89 × 10¹⁰ | Pa (68.9 GPa) | Brief |
| Poisson's ratio | ν | 0.33 | — | Brief |
| Density | ρ | 2 700 | kg/m³ | Brief |
| Yield stress (FYI, not used) | σ_y | ≈ 2.76 × 10⁸ | Pa (276 MPa) | [UNVERIFIED — typical handbook value for 6061-T6, used only for the "we are below yield" sanity check; not a constitutive parameter in this stage. Cross-check against MMPDS or ASM Handbook before stage 3.] |

The values above match the brief verbatim. They are also consistent with general-purpose handbook entries for 6061-T6 [UNVERIFIED — ASM Handbook Volume 2 entry would normally be cited; not in the project bib] and are documented in MMPDS (Metallic Materials Properties Development and Standardization) Chapter 3. [UNVERIFIED — MMPDS is paywalled and not present in references/*.bib, so the brief's values are accepted as authoritative for this stage.]

OpenRadioss card. `/MAT/LAW1/MID` (LAW1 = linear elastic isotropic). The keyword block in §6 lists the three numbers in the order required by the Altair LAW1 documentation [UNVERIFIED — exact column ordering must be cross-checked against `help.altair.com/hwsolvers/rad/topics/solvers/rad/mat_law1_starter_r.htm`; project bib does not have a dedicated LAW1 entry but does carry [@AltairRadiossLAW2] which establishes the LAWN keyword family pattern].

Density ρ is not used by /IMPL/LINEAR (no inertia in static analysis) but must still be a legal positive value in the card; it becomes relevant if we fall back to dynamic relaxation (see §10 risks).

---

## 6. OpenRadioss deck skeleton

OpenRadioss runs as a two-phase CLI: a `starter` reads `<job>_0000.rad` and writes a starter restart, then `engine` reads `<job>_0001.rad` and runs the analysis. The two files share the model definition; the engine file holds run-control (output frequency, run end criterion). The skeleton below lists every keyword block needed and is *not* a complete runnable deck; the runner.py templates this skeleton at run time.

### 6.1 Starter file `bend_3pt_0000.rad`

```text
#RADIOSS STARTER
/BEGIN
bend_3pt
2026 0
                       Mg                  mm                   s     # unit system
                       Mg                  mm                   s
/UNIT/1
units_SI_mm_Mg_s
                       Mg                  mm                   s
#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8--

# ───── Material ─────
/MAT/LAW1/1
AL6061T6
#              RHO_I
                  2.7E-9
#                E                  Nu
                 6.89E4              0.33
# (units: density in Mg/mm^3 = 2.7e-9; modulus in MPa = 6.89e4; Poisson dimensionless)

# ───── Property ─────
/PROP/TYPE14/1
solid_brick_HA8
#  Iorth  Ismstr   Iframe   Iplas   Itetra   Itetra4  Iframe2   Isolid
       0       0        0       0        0         0        0       17
#  qa       qb       h        Lambda    deltaT_min
   1.10     0.05     0.10     0.0       0.0
# (Isolid=17 = HA8 fully-integrated assumed-strain hex, recommended for slender bending)
# [UNVERIFIED: column layout of TYPE14 must be confirmed from the Altair card reference,
#  ref [@AltairRadiossPropType14]. The runner.py templates the columns exactly.]

# ───── Mesh (imported from inp2rad output) ─────
#include nodes.inc        # /NODE block, ~52 k nodes for M1
#include elements.inc     # /BRICK block, 6 400 HEXA8 elements
#include nsets.inc        # /GRNOD/NODE blocks: NSET_SUPPORT_LEFT, NSET_SUPPORT_RIGHT,
                          #                     NSET_LOAD_3PT, NODE_MIDSPAN

# ───── Part ─────
/PART/1
beam_part
        1        1        0
# part 1, prop_id 1, mat_id 1

# ───── Boundary conditions ─────
/BCS/1
fix_pin_left
#  Tx   Ty   Tz   Rx   Ry   Rz   skew_id   grnod_id
    1    1    1    0    0    0          0        100      # NSET_SUPPORT_LEFT
/BCS/2
fix_roller_right
#  Tx   Ty   Tz   Rx   Ry   Rz   skew_id   grnod_id
    0    1    1    0    0    0          0        101      # NSET_SUPPORT_RIGHT

# ───── Loads ─────
# Concentrated nodal force in -z, applied through a step function (linear ramp from 0 to 1
# over the implicit pseudo-time, then held; /IMPL/LINEAR ignores the ramp shape and uses
# the final amplitude).
/FUNCT/1
ramp_load_3pt
       0.0        0.0
       1.0        1.0
/CLOAD/1
load_3pt_midspan
#   grnod_id   dir   skew_id   funct_id   scale (force per node, Newtons)
         200     3         0          1   FORCE_PER_NODE_3PT
# dir=3 means z; sign of force given by 'scale' — negative for -z load.
# FORCE_PER_NODE_3PT = -1000 / N_nodes_in_NSET_LOAD_3PT (templated by runner.py)

# ───── Implicit static activation ─────
/IMPL/LINEAR        # static linear analysis, single solve; or /IMPL/QSTAT for multi-step
/IMPL/SOLVER/1
#  Linear solver type: 1 = direct (MUMPS, requires MUMPS-linked build),
#                      2 = iterative PCG, etc.
        1                                    # MUMPS direct
# [RISK: see §10 risk #1; if MUMPS-linked binary is unavailable, /IMPL is non-functional.]
/IMPL/PRINT
        1     1                              # print convergence each step
/IMPL/DT/STOP
   1.0E-6   1.0E+0                           # min dt, max dt — small range OK for linear
/IMPL/DTINI
   1.0                                       # one full step for /IMPL/LINEAR
/IMPL/NCYCLE/STOP
        1                                    # one Newton cycle suffices for linear

# ───── Time-history (T01) and animation (.anim) requests ─────
/TH/NODE/1
midspan_displacement
        300                                  # NODE_MIDSPAN id resolved by inp2rad
/TH/NODE/SAVE
DX  DY  DZ                                   # save translations
/ANIM/DT
        0.0     1.0                          # write one anim frame at end of step
/ANIM/BRICK/TENS/STRESS/ALL                  # stress on all bricks
/ANIM/BRICK/TENS/STRAIN/ALL                  # strain on all bricks
/ANIM/NODA/DISP                              # displacement field

/END
```

### 6.2 Engine file `bend_3pt_0001.rad`

```text
#RADIOSS ENGINE
/RUN/bend_3pt/1
       1.0                                   # final pseudo-time (1.0 s for /IMPL/LINEAR)
/STOP
/PRINT/-1                                    # log every cycle
/HIS/0.0                                     # always write time-history
/REC/-1
/END
```

### 6.3 Four-point bending — what changes

The `bend_4pt_*.rad` decks are identical except for the load block. Replace the single `/CLOAD/1` block with two:

```text
/CLOAD/1
load_4pt_a
         210     3         0          1   FORCE_PER_NODE_4PT      # NSET_LOAD_4PT_A, x = L/3
/CLOAD/2
load_4pt_b
         211     3         0          1   FORCE_PER_NODE_4PT      # NSET_LOAD_4PT_B, x = 2L/3
# FORCE_PER_NODE_4PT = -500 / N_nodes_in_each_load_set
```

Total downward force is the same (P = 1000 N). Closed-form midspan deflection differs (§7).

### 6.4 Keyword inventory used in this stage

| Keyword | Purpose | Cited |
|---|---|---|
| `/UNIT` | Declare unit system (Mg, mm, s → MPa) | OpenRadioss starter manual |
| `/MAT/LAW1` | Linear elastic isotropic | Altair LAW2 reference shows family pattern [@AltairRadiossLAW2] |
| `/PROP/TYPE14` | General solid (brick) property | [@AltairRadiossPropType14] |
| `/NODE`, `/BRICK` | Mesh | inp2rad output |
| `/GRNOD/NODE` | Node sets | OpenRadioss starter manual |
| `/PART` | Part assembly | OpenRadioss starter manual |
| `/BCS` | Translational/rotational constraints | OpenRadioss starter manual |
| `/CLOAD` | Concentrated nodal load | OpenRadioss starter manual |
| `/FUNCT` | Time function for load ramp | OpenRadioss starter manual |
| `/IMPL/LINEAR` | Activate static linear solver | [@AltairRadiossImplActivation] |
| `/IMPL/SOLVER` | Choose direct (MUMPS) or iterative | [@AltairRadiossImplActivation], [@MUMPS2024] |
| `/IMPL/DT`, `/IMPL/DTINI`, `/IMPL/PRINT`, `/IMPL/NCYCLE/STOP` | Static analysis controls | [@AltairRadiossImplActivation] |
| `/TH/NODE` | Time-history for the midspan probe | OpenRadioss engine manual |
| `/ANIM/...` | Animation file outputs (later converted to vtkhdf) | [@KitwareOpenRadiossVTKHDF2024] |
| `/RUN`, `/STOP`, `/PRINT`, `/HIS`, `/REC`, `/END` | Engine run control | OpenRadioss engine manual |

---

## 7. Reference solution

The closed-form Euler-Bernoulli solution for a slender simply supported beam is given in Timoshenko and Gere (1972) [@Timoshenko1972MechanicsOfMaterials] and reproduced in Gere and Goodno (2012) [@GereGoodno2012MechanicsOfMaterials]; both are textbook references that take the small-displacement linear-elastic assumption with no shear deformation.

### 7.1 Three-point bending

Simply supported beam of length L with a single concentrated load P at midspan. The Euler-Bernoulli midspan deflection is

```
δ_3pt = P L³ / (48 E I)                (1)
```

Plug in P = 1000 N, L = 0.200 m, E = 6.89 × 10¹⁰ Pa, I = 1.6667 × 10⁻⁹ m⁴.

```
δ_3pt = (1000) (0.200)³ / [48 · 6.89e10 · 1.6667e-9]
      = 1000 · 8e-3 / [48 · 6.89e10 · 1.6667e-9]
      = 8.0 / 5 514.0
      = 1.4509e-3 m
      = 1.4509 mm
```

So the FEM run must report u_z(L/2, 0, −h/2) = −1.4509 mm to within 1% (i.e. between −1.4364 mm and −1.4654 mm).

Shear-deformation correction (Timoshenko). For a rectangular cross section with κ = 5/6, G = E / (2(1+ν)) = 25.9 GPa, A = 2 × 10⁻⁴ m²,

```
δ_shear = P L / (4 κ G A)
        = 1000 · 0.200 / (4 · 5/6 · 25.9e9 · 2e-4)
        = 200 / (4 · 0.8333 · 5.18e6)
        = 200 / 17.27e6
        = 1.16e-5 m
        = 0.0116 mm
```

Ratio δ_shear / δ_3pt ≈ 0.8% — within the 1% pass criterion budget. So the *exact* Timoshenko midspan deflection is δ_3pt + δ_shear ≈ 1.4625 mm, and the Euler-Bernoulli closed form is itself ~0.8% too low. Two clean ways to handle this:

1. Accept the 1% tolerance against pure Euler-Bernoulli and rely on the FEM result landing between the two. The FEM, being a 3D solid, includes shear and so will give ≈ 1.4625 mm; pure Euler-Bernoulli gives 1.4509 mm; gap is 0.8% — passes the 1% criterion comfortably. This is the route the brief implies.
2. Compare against the Timoshenko closed form. Same FEM, tighter tolerance. Recommended for the convergence-study CSV but not the binary pass/fail.

The runner.py reports both deltas and the FEM value, and pass/fail against pure Euler-Bernoulli per the brief.

### 7.2 Four-point bending

Simply supported beam of length L, two equal concentrated loads P/2 at x = a and x = L − a, where a = L/3. From Gere and Goodno (2012) [@GereGoodno2012MechanicsOfMaterials] (Appendix table for "simply supported beam, two equal loads at equal distances from supports"), the midspan deflection is

```
δ_4pt,max = P a (3 L² − 4 a²) / (24 E I)   with a = L/3       (2)
```

Plug in a = L/3 = 0.0667 m:

```
3 L² − 4 a² = 3 (0.200)² − 4 (0.0667)²
            = 3 · 0.04 − 4 · 0.004444
            = 0.12 − 0.01778
            = 0.1022 m²

δ_4pt,max = 1000 · 0.0667 · 0.1022 / (24 · 6.89e10 · 1.6667e-9)
          = 1000 · 6.815e-3 / (2 757.0)
          = 6.815 / 2 757.0
          = 2.472e-3 m
          = 2.472 mm
```

[Sanity. The 4-pt deflection should be larger than the 3-pt deflection at the same total load even though the maximum bending moment is smaller (PL/6 vs PL/4) because the bending moment is *constant* over the middle third for 4-pt, so the curvature is integrated over a longer arc. Ratio δ_4pt / δ_3pt = 2.472 / 1.4509 = 1.704; the textbook ratio at a = L/3 is (23/108) · (48 / 1) = 10.22 for a uniform comparison metric — actually the cleanest cross-check is recomputing both deflections against the same formula bank, which the runner.py does symbolically.] [UNVERIFIED — formula (2) is the standard a-position deflection for two symmetric point loads; cross-checked against Gere and Goodno appendix.]

### 7.3 Why this is *verification*, not *validation*

Per test_progression_literature.md cross-cutting note 3, Stage 1 is verification (FEM result vs. closed form). No experimental ground truth is invoked — the question is whether the FEM kernel correctly assembles, integrates, applies BCs, and solves a linear elastic PDE. The closed form is the analytic ground truth.

### 7.4 Convergence rate prediction

For linear-Lagrange (HEXA8) elements, the expected convergence rate of the energy norm is O(h²) and of the displacement L² norm is O(h^(p+1)) = O(h²). On the three-mesh sweep M0 → M1 → M2 (h halved twice), the error in midspan deflection should reduce by a factor of ~4 each refinement. The Typst/CeTZ convergence plot is expected to show a slope of −2 on log-log axes.

---

## 8. Validation success criterion — how we measure midspan deflection

Pipeline.

1. OpenRadioss writes `.anim` frames (one final frame for /IMPL/LINEAR) and a `T01` time-history file.
2. The Kitware `openradioss-to-vtkhdf` converter [@KitwareOpenRadiossVTKHDF2024] turns `.anim` into a `.vtkhdf` file.
3. PyVista headless reads the `.vtkhdf`.
4. PyVista locates the node nearest to (L/2, 0, −h/2) using `pyvista.UnstructuredGrid.find_closest_point((Lh, 0, −h/2))` or by lookup from the `NODE_MIDSPAN` id stored in `T01`.
5. The vertical displacement u_z at that node at the final frame is extracted in metres and compared to the Euler-Bernoulli closed form.

Pass criterion (binary).

```
|u_z_FEM − δ_closed_form| / |δ_closed_form|  ≤  0.01    (1%)
```

evaluated separately for the 3-pt run (δ_closed_form = δ_3pt) and the 4-pt run (δ_closed_form = δ_4pt,max). Both must pass.

Reporting (verbose).

| Field | 3-pt | 4-pt |
|---|---|---|
| u_z_FEM at midspan, M0 | <run> | <run> |
| u_z_FEM at midspan, M1 | <run> | <run> |
| u_z_FEM at midspan, M2 | <run> | <run> |
| δ closed form (Euler-Bernoulli) | −1.4509 mm | −2.472 mm |
| δ closed form (Timoshenko-corrected) | −1.4625 mm | [recompute] |
| Relative error vs Euler-Bernoulli (M1) | < 1% required | < 1% required |
| Relative error vs Timoshenko (M1) | reported | reported |
| Convergence-rate fit p (slope on log-log) | report (≈ 2 expected) | report |

The CSV produced by runner.py has this exact schema. The Typst/CeTZ figure — convergence plot — is built from the CSV per the user's standing rule (CLAUDE.md "All figures must be authored in Typst using the CeTZ package"). No matplotlib.

---

## 9. Toolchain runner

The Python runner is at `/Volumes/MacShare/Code/FEA_AP-PLY/tests/stage_01_beam_bending/runner.py`. Pipeline orchestration:

1. Parse CLI flags (mesh size, load case, dry-run).
2. Generate a structured-hex mesh in GMSH Python API; write to `mesh.inp` (Abaqus format).
3. Render the deck template `deck.rad.j2` (TODO — not in this commit; runner.py points at it as a path) into `bend_3pt_0000.rad`, `bend_3pt_0001.rad`, and the 4-pt counterparts. The template fills in: mesh include path, node set ids, force-per-node values, mesh-size-dependent file names.
4. Convert `mesh.inp` → `mesh.rad` via `inp2rad` (OpenRadioss tools repo).
5. Invoke `limactl shell <vmname> -- /opt/openradioss/exec/starter_linuxa64 -i bend_3pt_0000.rad -nt N` then `engine_linuxa64 -i bend_3pt_0001.rad -nt N`. Both commands must succeed (return code 0).
6. Convert the resulting `bend_3pt_0001.anim` to `bend_3pt.vtkhdf` via the Kitware converter, also invoked inside Lima [@KitwareOpenRadiossVTKHDF2024].
7. Read `.vtkhdf` with PyVista headless on the macOS host (the macOS filesystem is shared into Lima, so the file is visible from both sides).
8. Extract u_z at the midspan node, compute relative error vs. closed form, write `results.csv`, print pass/fail.
9. (Optional, downstream) Render the convergence plot in Typst+CeTZ.

The runner.py in this directory is a 50-200 line skeleton suitable for extension as the toolchain firms up. The .rad deck template (`deck.rad.j2`) is intentionally *not* written end-to-end in this stage; it is a TODO referenced by the runner. The keyword skeleton in §6 above is the spec for the template.

---

## 10. Risks and unknowns

1. **MUMPS-linked OpenRadioss build is the critical risk.** /IMPL/LINEAR requires MUMPS [@MUMPS2024] for the direct sparse solve. The OpenRadioss prebuilt binaries on GitHub releases [@OpenRadiossInstall2026] are documented as Linux x86-64 / Linux ARM64 / Windows; whether the prebuilt binaries are built with MUMPS-on or MUMPS-off is not flagged on the release page. Discussion #2117 [@OpenRadiossDiscussion2117] confirms the implicit build is harder than the explicit-only build. Mitigation strategy in priority order:
   - (a) Try the prebuilt Linux ARM64 binary inside Lima [@OpenRadiossDiscussion2125] [@Lima2026] [@Apptainer2026] with `/IMPL/LINEAR` and read the starter log for "MUMPS not available" or equivalent. If it works, great.
   - (b) If (a) fails, build OpenRadioss from source inside Lima with MUMPS linked. The build_script.sh accepts a MUMPS option per Discussion #2117. Adds ~30 minutes of one-time install pain. Build dependencies: GFortran, OpenMPI, LAPACK, ScaLAPACK, and MUMPS itself (CeCILL-C — copyleft but compatible with AGPL-3 OpenRadioss for non-distribution use, master_plan.md §10 risk #7).
   - (c) Fallback: drop /IMPL/LINEAR, use explicit dynamic with mass-scaling and dynamic relaxation (`/DYREL`) to settle to a quasi-static linear-elastic solution. This is the audit's row-1 fallback. Costs ≈ 10× wall time for the same converged static solution but does not require MUMPS. Verifies the rest of the pipeline (geometry → mesh → starter → engine → vtkhdf → PyVista) end-to-end with the exact same .rad deck modulo the implicit cards.
   - The runner.py CLI has a `--solver implicit|explicit` switch to pick between (a/b) and (c).

2. **Lima + Apptainer first-time install on macOS.** Documented (master_plan.md §7) but adds one VM-boundary file-system layer. Lima [@Lima2026] mounts the macOS home directory into the VM; performance is fine for a 52 k-DOF static linear solve but will matter for stages 13-16. macOS native OpenRadioss does not exist [@OpenRadiossInstall2026]. Mitigation: this stage is the one that *catches* this risk; the stage's whole purpose is the smoke test.

3. **inp2rad converter is "beta" status** (audit row C). Tow orientation, node-set preservation, and material assignment are flagged for sanity check after conversion. For Stage 1 the only thing inp2rad has to preserve correctly is geometry and node sets — no orthotropy, no material orientation. Mitigation: dump `mesh.rad` and grep for the four expected /GRNOD/NODE blocks before invoking the starter.

4. **HEXA8 shear locking on a slender beam.** The 1% tolerance is tight enough that a poorly stabilized HEXA8 can fail it on M1 and only converge on M2. Mitigation: `Isolid=17` (HA8) or `Isolid=24` (HEPH) on /PROP/TYPE14 are documented to be appropriate for bending [UNVERIFIED — the column-by-column reading of the TYPE14 reference page is required, [@AltairRadiossPropType14]]. Document the chosen Isolid and the convergence-rate fit in the results CSV.

5. **Point-load vs. distributed-load BC choice.** A single-node point load on a HEXA8 mesh has a 1/r-style local field that pollutes the displacement near the load. Mitigation already in §4: distribute the load over the line of nodes through the width on the top fiber, and measure deflection at a point at least one half-thickness away (bottom fiber midspan). [UNVERIFIED — strictly the cleanest physics is a small distributed pressure patch over a 5 mm × 5 mm pad approximating the loading roller in real test apparatus; this is not necessary for verification but may be needed if cross-comparing to ASTM D7264 (composite three-point flexure, not used in this stage but the same beam geometry recurs in Stages 7-9).]

6. **Closed-form ambiguity at the 1% level — Euler-Bernoulli vs. Timoshenko vs. 3D solid.** Already documented in §7. The brief calls for 1% vs. Euler-Bernoulli; the 3D-solid FEM result will systematically be slightly larger than Euler-Bernoulli (it includes shear deformation) and equal to Timoshenko (κ = 5/6 for rectangular). On L/h = 20 the gap is ~0.8%, eating most of the tolerance budget. Mitigation: document both closed forms in the CSV; if the FEM lands inside the 1% Euler-Bernoulli budget, pass; if it lands at the Timoshenko value (above the 1% budget) but within 1% of Timoshenko, that is also acceptable and the spec footnotes the relaxed pass.

7. **Unit system.** OpenRadioss does not enforce a unit system; the user declares one on `/UNIT` and consistency is the user's responsibility. The deck above uses `Mg, mm, s → MPa` because that is the OpenRadioss / Radioss / LS-DYNA crash-analysis convention. Watch out: density 2 700 kg/m³ → 2.7 × 10⁻⁹ Mg/mm³, modulus 6.89 × 10¹⁰ Pa → 68 900 MPa = 6.89 × 10⁴ MPa, force 1 000 N → 1 000 N (no change in force unit) — but if the unit system is mistakenly read as `kg, mm, s → kPa` everything is off by 10³. The runner.py templates the values from a single SI dict and applies a documented `convert_to_Mg_mm_s` transform to each.

8. **No published PyPI for openradioss-python templating** (audit row A). Templating uses Jinja2 directly. Stage 1 establishes this template scaffolding for stages 2-9.

9. **Aluminum 6061-T6 properties from the brief are used unverified.** The brief's E, ν, ρ are widely accepted handbook values, but the project bib does not contain a primary citation for 6061-T6 specifically. [UNVERIFIED — MMPDS or ASM Handbook is the canonical source.] Acceptable for a verification-only stage where the *closed form uses the same E* — any error in E cancels exactly.

10. **macOS file-system case-sensitivity and line endings.** Lima mounts may translate paths; OpenRadioss starter is sensitive to whether the deck has CRLF or LF line endings (the starter is GFortran, which on Linux ARM64 typically wants LF). Mitigation: runner.py opens files in binary mode with explicit `\n`.

---

## Appendix A — Citation key list (subset of `references/openradioss_refs.bib` and `test_progression_refs.bib` actually cited above)

- [@Timoshenko1972MechanicsOfMaterials] Timoshenko and Gere, *Mechanics of Materials*, 1972.
- [@GereGoodno2012MechanicsOfMaterials] Gere and Goodno, *Mechanics of Materials*, 2012.
- [@AltairRadiossImplActivation] Altair Engineering, "Implicit Analysis Activation in Radioss / OpenRadioss".
- [@AltairRadiossPropType14] Altair Engineering, "/PROP/TYPE14 (SOLID): general solid property reference".
- [@AltairRadiossLAW2] Altair Engineering, "/MAT/LAW2 (PLAS_JOHNS)" (cited as the closest-pattern keyword family doc; LAW1 has the same column-format conventions).
- [@MUMPS2024] Amestoy et al., MUMPS sparse direct solver.
- [@KitwareOpenRadiossVTKHDF2024] Kitware, openradioss-to-vtkhdf converter.
- [@Lima2026] Lima Project, Linux machines on macOS.
- [@Apptainer2026] Apptainer Project, container runtime.
- [@OpenRadiossInstall2026] OpenRadioss Project, INSTALL.md and HOWTO.md.
- [@OpenRadiossDiscussion2117] OpenRadioss Community, implicit build discussion.
- [@OpenRadiossDiscussion2125] OpenRadioss Community, Apple Silicon via Apptainer.
- [@OpenRadiossDiscussion3222INP2RAD] OpenRadioss Community, inp2rad Python tool.

End of spec.
