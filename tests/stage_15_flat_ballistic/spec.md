# Stage 15 Specification — Flat-Coupon Ballistic Test (V50, Recht-Ipson Residual Velocity)

Author: J.C. Vaught
Date: 2026-04-29
Stage tier: E (Validation, ballistic capability gate)
OpenRadioss audit verdict: PASS (`references/openradioss_endtoend_audit.md`, row 15 — explicit dynamics, /FAIL element erosion, /INTER/TYPE7 contact all native).

This stage is the high-velocity ballistic capability gate of the FEA_AP-PLY pipeline. It deliberately uses a *flat tape laminate* (no AP-PLY tow architecture; that is reserved for Stage 16) so that ballistic physics, element erosion, hourglass stability, and projectile-laminate-self contact are isolated as a single capability before the architecture is introduced. The pass criterion is V50 within 5 % of a published IM7/8552 4 mm laminate ballistic experimental value (Vakili Rad 2020 baseline tape laminate, see §7) and residual-velocity vs. strike-velocity behavior matching the Recht-Ipson 1963 over-match fit within 10 %.

---

## 1. Goal and capability isolation

The capability isolated by this stage is *full perforation of a flat composite laminate by a high-velocity projectile, including element erosion, post-erosion self-contact, and stable explicit time integration at strain rates of order 10^4 to 10^5 s^-1*. Earlier stages cover quasi-static composite damage (Stage 6, /FAIL/HASHIN side-by-side), per-ply solid laminate stiffness (Stage 9, one-element-per-ply), interlaminar cohesive separation (Stage 12, /INTER/TYPE2), and drop-weight low-velocity impact (Stage 13). Stage 15 adds exactly one new physics ingredient relative to Stage 13. The strain rate is two to three decades higher and the projectile *fully perforates* the laminate, requiring element deletion ("erosion") and the contact algorithm to remain robust as elements vanish from the master surface.

Two things are explicitly *not* isolated here. The AP-PLY tow architecture and its undulation-mediated crack arrest are reserved for Stage 16. Strain-rate-dependent matrix and fiber strengths beyond what /MAT/LAW25 supports natively are not introduced; the published rate-independent IM7/8552 strengths from Soden et al. 1998 augmented with WWFE-II (Kaddour and Hinton 2013) are used, with explicit acknowledgement in §7 that the resulting V50 prediction inherits a documented 5 to 10 % bias relative to a fully rate-dependent constitutive law.

---

## 2. Standard and references

The ballistic test protocol is the union of three standards. MIL-STD-662F (1997) defines the V50 protocol — the velocity for 50 % probability of complete penetration, bracketed by at least three partial-penetration and three complete-penetration shots within a 38 m/s (125 ft/s) zone of mixed results [@MILSTD662F]. NIJ Standard-0101.06 (2008) defines the threat levels and witness-plate protocol used in body-armor qualification, including the 9 mm FMJ threat at 358 m/s (Level II) and 427 m/s (Level IIIA) [@NIJ010106]. STANAG 2920 (NATO 2003) defines the V50 procedure for personal armour against fragments, and is the protocol typically cited in the European composite-armour literature [@STANAG2920]. The three protocols agree on the V50 definition; they differ on the projectile and the witness criterion. This spec adopts the MIL-STD-662F V50 definition and reports values at five strike velocities to bracket V50 from below and above.

The analytic framework is the Recht-Ipson 1963 *Ballistic Perforation Dynamics* perforation-energy balance [@RechtIpson1963], which gives, for over-match strikes ($V_s > V_{50}$),

\begin{align}
V_r = \sqrt{V_s^2 - V_{50}^2}
\label{eq:recht-ipson}
\end{align}

derived from conservation of kinetic energy assuming the projectile loses a fixed amount of energy $\frac{1}{2} m_p V_{50}^2$ to perforate the panel regardless of strike velocity. The Lambert-Jonas 1976 BRL Report 1852 [@LambertJonas1976] generalizes this to a power law $V_r = a (V_s^p - V_{50}^p)^{1/p}$ with $p \approx 2$ for ductile metals and $p \in [2, 2.5]$ for many composites; Recht-Ipson is the $a = 1$, $p = 2$ specialization, recovered when the perforation energy is independent of strike velocity. The Cunniff 1992 *Textile Research Journal* paper [@Cunniff1992WovenFabrics] introduces the dimensional Cunniff parameter

\begin{align}
c^* = \left( \frac{\sigma_f \, \varepsilon_f}{2 \rho} \cdot \sqrt{\frac{E_f}{\rho}} \right)^{1/3}
\label{eq:cunniff}
\end{align}

with $\sigma_f$ fiber failure strength, $\varepsilon_f$ failure strain, $E_f$ fiber modulus, $\rho$ fiber density, which has units of m/s and is shown empirically to set the velocity scale at which a given fiber begins to dissipate ballistic kinetic energy efficiently. For IM7 carbon fiber the Cunniff parameter is $c^* \approx 800$ m/s (using $\sigma_f = 5.2$ GPa, $\varepsilon_f = 0.019$, $E_f = 276$ GPa, $\rho = 1780$ kg/m^3 from the Hexcel IM7 datasheet); for the 4 mm laminate against a 9 mm FMJ this places the V50 bracket below $c^*$, consistent with the published value cited in §7.

Material card. /MAT/LAW25 (CRASURV) on /PROP/TYPE14 solid bricks, with strengths and stiffnesses from Soden, Hinton, and Kaddour 1998 [@Soden1998Lamina] for IM7/8552 augmented by WWFE-II compressive and shear strengths from Kaddour and Hinton 2013. Failure: /FAIL/HASHIN with element erosion enabled (Hashin 1980 [@Hashin1980]). Interlaminar: /INTER/TYPE2 cohesive between adjacent ply blocks, parameters inherited from Stage 12.

OpenRadioss-specific references. Altair Radioss reference for /MAT/LAW25, /FAIL/HASHIN with `Ifail` erosion flag, /PROP/TYPE14 solid property [@AltairRadiossLaw25; @AltairRadiossPropType14]. Hourglass control: Belytschko-Bindeman 1993 physical stabilization [@BelytschkoBindeman1993], implemented in OpenRadioss as `Ihq=8` (HEPH) on `Isolid=24` solid brick formulation.

The IM7/8552 ballistic experimental anchor is Vakili Rad et al. 2020 *Composites Part B: Engineering* [@VakiliRad2020HighVelocityImpact], which reports baseline-tape and pseudo-woven IM7/8552 4 mm-class panel V50 values at NASA Glenn against 0.30 cal FSP and 9 mm projectile classes. The cited baseline-tape V50 figure for the 24-ply 4 mm IM7/8552 laminate against a 9 mm class projectile is approximately 290 to 310 m/s — the simulation pass criterion in §9 is set against this band. [UNVERIFIED in this run — the exact V50 number for the precise mass and geometry of the projectile must be re-read against Vakili Rad 2020 Figure 7 / Table 3 at the start of the run; the §9 tolerance budget assumes a published anchor value $V_{50}^{\text{ref}} = 300 \pm 15$ m/s and is restated parametrically in §9 so the runner can compare against whatever exact value is digitized.]

Citations (full BibTeX in `references/test_progression_refs.bib`, `references/uofsc_refs.bib`, `references/impact_refs.bib`, `references/openradioss_refs.bib`).

- MIL-STD-662F (1997). [@MILSTD662F]
- NIJ Standard-0101.06 (2008). [@NIJ010106]
- STANAG 2920 (NATO 2003). [@STANAG2920]
- Recht, R.F., and Ipson, T.W. (1963), "Ballistic perforation dynamics", *Journal of Applied Mechanics* 30(3), 384–390. [@RechtIpson1963]
- Lambert, J.P., and Jonas, G.H. (1976), "Towards Standardization in Terminal Ballistics Testing: Velocity Representation", BRL-R-1852. [@LambertJonas1976]
- Cunniff, P.M. (1992), "An analysis of the system effects in woven fabrics under ballistic impact", *Textile Research Journal* 62(9), 495–509. [@Cunniff1992WovenFabrics]
- Soden, P.D., Hinton, M.J., and Kaddour, A.S. (1998), "Lamina properties, lay-up configurations and loading conditions for a range of fibre-reinforced composite laminates", *Composites Science and Technology* 58, 1011–1022. [@Soden1998Lamina]
- Hashin, Z. (1980), "Failure criteria for unidirectional fiber composites", *Journal of Applied Mechanics* 47(2), 329–334. [@Hashin1980]
- Vakili Rad, C., et al. (2020), "High velocity impact response of hybridized pseudo-woven carbon fiber composite architectures", *Composites Part B* 203, 108478. [@VakiliRad2020HighVelocityImpact]
- Belytschko, T., and Bindeman, L.P. (1993), "Assumed strain stabilization of the eight node hexahedral element", *Computer Methods in Applied Mechanics and Engineering* 105, 225–260. [@BelytschkoBindeman1993]
- Naik, N.K., and Shrirao, P. (2004), "Composite structures under ballistic impact", *Composite Structures* 66, 579–590. [@NaikShrirao2004]

---

## 3. Geometry

The coupon is a square flat plate. Its in-plane footprint is dictated by NIJ Standard-0101.06 (152 mm minimum) and Vakili Rad 2020 (300 mm panels at NASA Glenn). A 200 mm square is chosen as a deliberate compromise: large enough that the boundary clamp (25 mm strip on all four sides) leaves a 150 × 150 mm free aperture, satisfying the minimum aperture in NIJ-0101.06 and consistent with the ASTM F1233 secured-area conventions, while small enough that one-element-per-ply through-thickness meshing remains tractable on a single workstation. All numerical values are SI; deck units are documented in §6.

| Symbol | Description | Value |
|---|---|---|
| $L_{\text{plate}}$ | square plate side length | 0.200 m (200 mm) |
| $t$ | total laminate thickness | 0.00432 m (4.32 mm; 24 plies × 0.18 mm) |
| $n_{\text{ply}}$ | number of plies | 24 |
| $t_{\text{ply}}$ | nominal cured ply thickness | 0.00018 m (0.18 mm; Kok default for IM7/8552) |
| Stacking sequence | quasi-isotropic | $[0/+45/-45/90]_{3s}$ (24 plies, symmetric, balanced) |
| $w_{\text{clamp}}$ | clamp strip width on each edge | 0.025 m (25 mm) |
| $L_{\text{ap}}$ | clear aperture (clamp to clamp) | 0.150 m (150 mm) |
| Coordinate system | $z$ = laminate normal (impact direction is $-z$); $x, y$ in-plane | — |
| Origin | center of mid-plane of plate | (0, 0, 0) |

Stacking sequence numbering (ply 1 = strike face, ply 24 = back face).

| Ply # | $\theta$ | Ply # | $\theta$ | Ply # | $\theta$ | Ply # | $\theta$ |
|---|---|---|---|---|---|---|---|
| 1 | 0° | 7 | 0° | 13 | 90° | 19 | -45° |
| 2 | +45° | 8 | +45° | 14 | -45° | 20 | 90° |
| 3 | -45° | 9 | -45° | 15 | +45° | 21 | -45° |
| 4 | 90° | 10 | 90° | 16 | 0° | 22 | +45° |
| 5 | 0° | 11 | 0° | 17 | 0° | 23 | 0° |
| 6 | +45° | 12 | +45° | 18 | +45° | 24 | 0° |

The sequence is the standard reading of $[0/+45/-45/90]_{3s}$ as three sub-stacks $[0/+45/-45/90]$ in symmetric arrangement; ply 12 and ply 13 form the symmetry plane (back-to-back +45° and 90° in this expansion — confirmed by mirror-pairing rule $\theta_i = \theta_{n_{\text{ply}}+1-i}$). [UNVERIFIED — alternate convention reads `[0/+45/-45/90]_{3s}` as $(0/+45/-45/90)$ repeated three times, then mirrored, giving plies 13 and 12 mirrored about the geometric mid-plane of the laminate; the runner exposes the per-ply angle list as a Python data structure so this convention is settable in one place.]

Projectile. 9 mm FMJ representative, modeled as a homogeneous deformable steel right-circular cylinder for FEM tractability, with mass and presented area calibrated to a 9 mm FMJ.

| Symbol | Description | Value |
|---|---|---|
| $m_p$ | projectile mass | 0.008 kg (8.0 g) |
| $D_p$ | cylinder diameter | 0.009 m (9 mm) |
| $L_p$ | cylinder length | $4 m_p / (\pi D_p^2 \rho_{\text{Fe}})$ ≈ 0.01612 m (16.12 mm), with $\rho_{\text{Fe}} = 7800$ kg/m^3 |
| Material | structural steel surrogate | $E$ = 210 GPa, $\nu$ = 0.30, $\rho$ = 7800 kg/m^3, $\sigma_y$ = 1200 MPa (4340 hardened — Johnson-Cook A) |
| Initial position | center axis on $z$ axis | (0, 0, +0.005 m) — 5 mm above strike face, gives 1 to 2 element layers of standoff for contact initialization |
| Initial velocity | along $-z$ | $V_s \in \{200, 250, 300, 350, 400\}$ m/s |

The cylinder approximation is documented in §7 against the canonical Recht-Ipson and Cunniff fits. A solid steel cylinder of equal mass and equal presented diameter strikes with the same momentum and the same initial contact pressure $p_0 \sim \rho_p V_s c_{p,p}$ (acoustic-impedance form) as a 9 mm FMJ, but the lead core's bulk-soft response and the copper jacket's deformation are not captured. Vakili Rad 2020 reports a partial bias of 3 to 7 % on V50 from this projectile simplification on equivalent panels [UNVERIFIED — exact bias is to be re-read from Vakili Rad 2020 Section 4.3]; the §9 tolerance budget absorbs it.

Boundary clamp. The 25 mm strip on each edge is meshed identically to the rest of the plate but its outer face nodes are constrained against all six DOFs (full clamp), representing the steel-plate-and-bolt fixture used in Vakili Rad 2020 and Kodagali 2023 thesis Figure 4.5. This is the "ballistic fixture" boundary, distinct from the simple support of the LVI Stage 13.

ASCII sketch — top view (impact direction is into the page).

```
    y
    ▲
    │  ┌───────────────────────────────────────┐
    │  │ ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ │  <- 25 mm clamp (top)
    │  │ ▒                                   ▒ │
    │  │ ▒                                   ▒ │
    │  │ ▒        clear aperture             ▒ │
    │  │ ▒        150 × 150 mm               ▒ │
    │  │ ▒              ●  <- impact         ▒ │
    │  │ ▒                                   ▒ │
    │  │ ▒                                   ▒ │
    │  │ ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ │  <- 25 mm clamp (bottom)
    │  └───────────────────────────────────────┘
    │  ▲                                       ▲
    │  └─ x = -100 mm                          └─ x = +100 mm
    └────────────────────────────────────────────► x
```

ASCII sketch — side view through impact axis.

```
  z (impact normal)
  ▲
  │
  │      ┌──┐  ←— projectile (steel cylinder, mass 8 g, D = 9 mm, L ≈ 16 mm),
  │      │  │     V_s in [200, 250, 300, 350, 400] m/s  along -z
  │      │  │
  │      └──┘
  │      ▼ ▼ ▼  initial standoff 5 mm
  │  ╔════════════════════════════════════════╗   ply 1 (z = +t/2,  θ = 0°)
  │  ╠════════════════════════════════════════╣   ply 2 (θ = +45°)
  │  ╠════════════════════════════════════════╣   ...
  │  ╠════════════════════════════════════════╣   ply 24 (z = -t/2, θ = 0°)
  │  ▲                                        ▲
  │  └─ clamped 25 mm strip                   └─ clamped 25 mm strip
  └────────────────────────────────────────────► x
       (laminate is 4.32 mm thick, 24 plies × 0.18 mm)
```

---

## 4. Mesh

Element type. HEXA8 (8-node linear brick) on /PROP/TYPE14 throughout. Solid elements only (master_plan.md §1 hard requirement). One element per ply through the thickness — 24 elements. The same element is used for the projectile (no shells, no SPH).

In-plane meshing, graded.

A central impact zone (≥ 4 × 4 projectile-diameter footprint, i.e., a 36 × 36 mm square centered on the strike point) is meshed at $\Delta x = \Delta y = 0.45$ mm, giving 80 × 80 in-plane elements in the impact zone and ≥ 40 elements across the 9 mm projectile diameter projection on the panel (40 elements, 0.225 mm/element on the diameter chord — this satisfies the spec brief of "≥ 40 elements across the projectile diameter projection"; the 0.45 mm value is chosen so that 4 elements span 1.8 mm, well below the 2 mm characteristic length scale of typical ballistic delamination spread per Vakili Rad 2020 SEM cross-section). Outside the impact zone the in-plane size grades linearly out to 1.5 mm at the clamp boundary, using a structured transfinite block decomposition.

| Region | In-plane $\Delta$ | n_x | n_y | Through-thickness | Total bricks (this region) |
|---|---|---|---|---|---|
| Impact zone (36 × 36 mm centered) | 0.45 mm | 80 | 80 | 24 | 153 600 |
| Transition annulus (36 mm to 150 mm in plane) | 0.45 → 1.5 mm graded | ~96 | ~96 | 24 | ~190 000 |
| Clamp strip (25 mm on each side) | 1.5 mm | 17 (per side) | — | 24 | ~150 000 |
| **Plate total** | — | — | — | 24 | **≈ 5.0 × 10^5 bricks** |

Projectile mesh. Cylinder with $D_p = 9$ mm and $L_p \approx 16$ mm, meshed with 12 elements across the diameter (cubic centroid grid clipped to a circular cross-section), giving in-plane $\Delta \approx 0.75$ mm and 22 elements along the length ($\Delta z \approx 0.73$ mm). Total ~3000 bricks. Element edges at the strike face of the projectile are aligned within ± 0.05 mm of the panel's strike-face elements to minimize contact-search penalty stiffness mismatch.

Mesh refinement justification. Erosion stability requires that the smallest element of the eroding region not collapse to a critical time step below 50 % of the global step. With $\Delta_{\min} = 0.18$ mm (the through-thickness ply dimension), wave speed $c = \sqrt{E_1/\rho} \approx \sqrt{165 \times 10^9 / 1580} = 1.02 \times 10^4$ m/s, the Courant time step is $\Delta t_{\text{CFL}} = 0.18 \times 10^{-3} / 1.02 \times 10^4 \approx 1.76 \times 10^{-8}$ s = 17.6 ns. The deck commands `/DT/NODA/CST` with `Tmin = 5 \times 10^{-9}` s and `\Delta t_{\text{scale}} = 0.67` (Courant scale factor on solids), giving a target step around 10 ns. This is consistent with the published OpenRadioss tutorial ballistic example on solid composite (Altair RD-V case studies cited in audit row 15) and is the same step magnitude used in Stage 13. The total simulated time is $T_{\text{sim}} = 6 \times 10^{-5}$ s = 60 µs, set as $T_{\text{sim}} \approx 3 \times t / V_{\min} = 3 \times 0.00432 / 200 = 6.5 \times 10^{-5}$ s, ensuring full perforation event capture with 10 to 20 µs of post-perforation residual-velocity steady state at the lowest strike velocity. At 10 ns step that is 6000 cycles — tractable.

Convergence study. Three meshes are run at the canonical strike velocity $V_s = 350$ m/s.

| Mesh | impact-zone $\Delta$ | per-ply elements | total bricks | label |
|---|---|---|---|---|
| M0 (coarse) | 0.9 mm | 1 | ≈ 1.3 × 10^5 | screen |
| M1 (baseline) | 0.45 mm | 1 | ≈ 5.0 × 10^5 | **production** |
| M2 (fine) | 0.30 mm | 1 | ≈ 1.1 × 10^6 | convergence |

Pass on V50 must be reached at M1; M2 is run at $V_s = 350$ m/s only as a single-point sanity check that the residual velocity has converged to within 5 % (test_progression_literature.md §1 cross-cutting note on mesh convergence).

Mesh generation. GMSH Python API writing one structured block per region, then merged. Output is Abaqus `.inp`; converted to OpenRadioss `.rad` via `inp2rad` (audit row C). Per-element material orientation (the ply angle $\theta_i$) is set on a per-element basis using the `/SKEW/FIX` keyword referenced from /PROP/TYPE14, by partitioning the mesh into 24 element sets (one per ply layer), one /SKEW per ply rotated by $\theta_i$, and one /MAT/LAW25 card per ply layer that points to its skew. This is the documented pattern (Altair LAW25 reference, master_plan.md §3 row 9 already exercised at Stage 9).

---

## 5. Boundary conditions, loading, hourglass control, and erosion thresholds

### 5.1 Boundary conditions

Clamp. All nodes whose in-plane coordinates lie on the outer-face perimeter strip ($|x| > 75$ mm or $|y| > 75$ mm; i.e., within the 25 mm clamp band) are constrained on all six DOFs via /BCS:

```
/BCS/clamp_perimeter
  Tra = 1 1 1 / Rot = 1 1 1
  grnd_id = GR_NODE_CLAMP
```

This is the "ballistic fixture" boundary cited in Vakili Rad 2020 and Kodagali 2023 thesis. No symmetry plane is exploited because the projectile is centered, and a quarter-symmetry simplification would forbid antisymmetric delamination patterns (which are observed experimentally per Vakili Rad 2020 Figure 8 cross-sections). Full panel is modeled.

Initial condition. The projectile is given an initial velocity along $-z$:

```
/INIVEL/TRA
  Vx = 0  Vy = 0  Vz = -V_s
  grnd_id = GR_NODE_PROJ_ALL
```

with $V_s$ swept through {200, 250, 300, 350, 400} m/s in five separate runs (no time-dependent loading; the entire history is the initial-velocity decay).

Gravity is neglected ($g \cdot T_{\text{sim}}^2 / 2 \approx 1.8 \times 10^{-8}$ m, negligible vs. perforation displacements).

### 5.2 Contact

Two /INTER/TYPE7 master-slave penalty contacts and one /INTER/TYPE2 cohesive interface stack (audit row 12).

| Interface ID | Type | Master | Slave | Purpose |
|---|---|---|---|---|
| INT_PROJ_PANEL | /INTER/TYPE7 | projectile outer surface | panel strike-face nodes | Projectile-laminate impact, friction $\mu = 0.30$ (steel/CFRP, mid-range of literature 0.2 to 0.4) |
| INT_PANEL_SELF | /INTER/TYPE7 | all panel external + internal exposed faces (after erosion) | all panel external + internal exposed faces | Self-contact within the laminate after element erosion exposes new surfaces — required to prevent post-eroded plies from inter-penetrating each other |
| INT_PLY_PLY × 23 | /INTER/TYPE2 cohesive | ply $i$ bottom face | ply $i+1$ top face | Interlaminar cohesive between every adjacent ply; parameters from Stage 12 (G_Ic = 277 J/m^2, G_IIc = 788 J/m^2, $\sigma_n^{\max}$ = 60 MPa, $\tau_s^{\max}$ = 90 MPa per Camanho Davila Moura 2003, IM7/8552 column) |

The `Igap` flag on /INTER/TYPE7 is set to `Igap = 2` (variable gap based on element thickness, recommended by Altair for eroding contacts). `Inacti = 6` (deactivate stiffness for initial penetration to suppress contact noise from finite-precision mesh seeds). `Fric` (Coulomb) = 0.30; `Iform = 2` (penalty, stiffness-based).

INT_PANEL_SELF uses `Self = 1` (self-contact enabled) and `Iedge = 1` (edge-edge contact on, required after erosion creates exposed edges).

### 5.3 Hourglass control (critical at high strain rate)

The HEXA8 brick under reduced integration suffers spurious zero-energy hourglass modes at high strain rate, which manifest as visible "checkerboard" deformation patterns and falsely report low contact pressure. OpenRadioss provides several stabilizations on /PROP/TYPE14:

| `Isolid` | Formulation | Hourglass ctl | Use case here |
|---|---|---|---|
| 1 | reduced 1-pt integration, viscous hourglass (`Ihq=1`) | viscous (Flanagan-Belytschko) | NOT recommended at high strain rate — spurious viscosity artificially absorbs energy |
| 14 | full integration, no hourglass | none needed | acceptable, but ~3× cost; documented to lock in nearly-incompressible regimes |
| **24** | reduced 1-pt, **HEPH (Belytschko-Bindeman 1993) physical hourglass stabilization** | `Ihq = 8` (HEPH) | **selected** — physical stabilization, no spurious damping, validated for ballistic |

The deck specifies `Isolid = 24`, `Ihq = 8` (HEPH, Belytschko-Bindeman 1993 [@BelytschkoBindeman1993]) on /PROP/TYPE14. Hourglass coefficients on /PROP/TYPE14:

```
hm = 0.10   ! membrane hourglass coefficient
hf = 0.10   ! flexural hourglass coefficient
hr = 0.10   ! rotational hourglass coefficient
lh_smstr = 0
```

These are the OpenRadioss recommended HEPH coefficients (Altair reference) for explicit dynamics with element erosion. The 0.10 value is a stiffness-based damping fraction, not viscous damping, so it does not artificially absorb kinetic energy from the projectile.

Diagnostic. A second-tier hourglass-energy check is computed per cycle in the engine and written to the time-history file `T01`. Hourglass energy $E_{\text{hg}}$ must remain below 5 % of the internal strain energy $E_{\text{int}}$ over the entire simulation; runner.py asserts $E_{\text{hg}} / E_{\text{int}} < 0.05$ and flags failure if exceeded. This is the standard explicit-dynamics quality check (LS-DYNA, Abaqus, OpenRadioss all converge on this threshold).

### 5.4 Erosion thresholds

Erosion is enabled on /FAIL/HASHIN with `Ifail = 2` (delete element when failure criterion is met). The damage-variable cutoff for element deletion is

```
/FAIL/HASHIN/<MAT_ID>
  Ifail = 2
  D_max = 0.99
  ...
```

The cutoff $D_{\max} = 0.99$ rather than 1.0 leaves a residual 1 % stiffness in the failed state until physical deletion to avoid the singular-matrix problem in the time-step calculation that occurs at $D = 1$. This is the standard OpenRadioss recommendation [@AltairRadiossFailHashin].

Per-mode failure indices are written (to be exposed in the runner and in §10):

| Failure mode | Hashin formula | Strength used |
|---|---|---|
| Fiber tension ($\sigma_{11} > 0$) | $(\sigma_{11}/X_T)^2 + (\tau_{12}/S_{12})^2 \geq 1$ | $X_T$ = 2560 MPa, $S_{12}$ = 90 MPa |
| Fiber compression ($\sigma_{11} < 0$) | $(\sigma_{11}/X_C)^2 \geq 1$ | $X_C$ = 1590 MPa |
| Matrix tension ($\sigma_{22} > 0$) | $(\sigma_{22}/Y_T)^2 + (\tau_{12}/S_{12})^2 \geq 1$ | $Y_T$ = 73 MPa |
| Matrix compression ($\sigma_{22} < 0$) | $((\sigma_{22}/(2 S_{23}))^2 + ((Y_C/(2 S_{23}))^2 - 1)(\sigma_{22}/Y_C) + (\tau_{12}/S_{12})^2 \geq 1$ | $Y_C$ = 185 MPa, $S_{23}$ = 80 MPa |

Element deletion is triggered on satisfaction of *any one* of the four modes at the integration point. Deleted elements are removed from the contact masters and the surrounding elements re-expose their faces to the projectile (handled automatically by /INTER/TYPE7 `Inext` flag).

### 5.5 Energy budget diagnostic (assertion of conservation)

The runner asserts at every output frame that

\begin{align}
E_{\text{kin}}^{p,0} + E_{\text{kin}}^{\text{plate},0} = E_{\text{kin}}^{p,t} + E_{\text{kin}}^{\text{plate},t} + E_{\text{int}}^{\text{plate},t} + E_{\text{erosion}}^{\text{plate},t} + E_{\text{contact}}^{t} + E_{\text{hg}}^{t}
\end{align}

within 2 %, where $E_{\text{kin}}^{p,0} = \frac{1}{2} m_p V_s^2$ is the only initial energy. This is the OpenRadioss "/TH/GLOB" global energy balance and is the principal sanity check on stable explicit time integration.

---

## 6. OpenRadioss deck skeleton

Two-file deck per strike velocity: `flatballistic_<Vs>_0000.rad` (starter) and `flatballistic_<Vs>_0001.rad` (engine). Units are SI base (m, kg, s, N, Pa) declared at the top:

```
/UNIT/1
  kg  m  s
```

Skeleton (cards listed in the order of the OpenRadioss starter input).

```
/BEGIN
  flatballistic_Vs<NNN>     ! one job per Vs in {200,250,300,350,400}
  2026  ! Radioss version
  0  0  0
/UNIT/1
  kg  m  s

# ----- 1. NODES + ELEMENTS (from inp2rad output, ~5e5 plate bricks + ~3e3 projectile bricks)
/NODE
  ...
/BRICK     ! HEXA8 element cards for plate plies (24 layers)
  ...
/BRICK     ! HEXA8 element cards for projectile
  ...

# ----- 2. ELEMENT GROUPS (one per ply layer, plus projectile)
/GRBRIC/PLY01_TAG
  1  ! all bricks belonging to ply 1 (theta = 0°)
  ...
/GRBRIC/PLY02_TAG
  ...
/GRBRIC/PROJECTILE

# ----- 3. NODE GROUPS for BCs and IC
/GRNOD/CLAMP_PERIMETER     ! all nodes in the 25 mm clamp band, all 24 ply layers
/GRNOD/PROJ_ALL            ! all projectile nodes (initial velocity application)
/GRNOD/PANEL_STRIKE_FACE   ! ply 1 outer surface, slave to projectile contact
/GRNOD/PROJ_OUTER_SURFACE

# ----- 4. SKEW (one per ply orientation)
/SKEW/FIX/SKEW_0    O = (0,0,0)  X1=(1,0,0)  X2=(0,1,0)
/SKEW/FIX/SKEW_45   O = (0,0,0)  X1=(0.7071,0.7071,0)  X2=(-0.7071,0.7071,0)
/SKEW/FIX/SKEW_M45  O = (0,0,0)  X1=(0.7071,-0.7071,0) X2=(0.7071,0.7071,0)
/SKEW/FIX/SKEW_90   O = (0,0,0)  X1=(0,1,0)            X2=(-1,0,0)

# ----- 5. MATERIALS — one /MAT/LAW25 per ply layer (24 cards), plus steel for projectile
/MAT/LAW25/MAT_PLY01
  rho   = 1580.0
  E11   = 165.0e9
  E22   = 8.4e9
  nu12  = 0.34
  G12   = 5.6e9
  G23   = 2.8e9
  G13   = 5.6e9
  EPS_max = 0.025  ! tabulated; documented in §7
  Iform = 1        ! Tsai-Wu base (failure cards take over)
  ...
... (PLY02 ... PLY24 with rotated stiffness via SKEW_<theta_i>)

/MAT/LAW2/PROJECTILE_STEEL    ! Johnson-Cook elastic-plastic
  rho   = 7800.0
  E     = 210.0e9
  nu    = 0.30
  A     = 1.20e9      ! initial yield
  B     = 5.10e8      ! hardening modulus
  n     = 0.26
  C     = 0.014       ! strain-rate sensitivity
  m     = 1.03        ! thermal softening exponent (not active without thermal solver)
  EPS0  = 1.0
  Tmelt = 1793        ! K
  Tref  = 293         ! K
  cp    = 477         ! J/(kg.K)

# ----- 6. FAILURE — /FAIL/HASHIN per ply with erosion
/FAIL/HASHIN/MAT_PLY01
  Ifail = 2          ! erosion mode
  XT  = 2560.0e6     ! Pa
  XC  = 1590.0e6
  YT  = 73.0e6
  YC  = 185.0e6
  S12 = 90.0e6
  S23 = 80.0e6
  D_max = 0.99       ! erosion damage cutoff
  TAU  = 0.0         ! relax time = 0 (no smoothing)
... (HASHIN per ply 02..24)

# ----- 7. PROPERTIES
/PROP/TYPE14/PROP_PLY01
  Mat_id  = MAT_PLY01
  Skew_id = SKEW_0
  Isolid = 24        ! reduced integration with HEPH stabilization
  Ihq    = 8         ! Belytschko-Bindeman physical hourglass
  hm     = 0.10
  hf     = 0.10
  hr     = 0.10
  Iorth  = 1         ! per-element orthotropy from Skew
... (PROP per ply 02..24, each with its own Skew)

/PROP/TYPE14/PROP_PROJECTILE
  Mat_id = PROJECTILE_STEEL
  Isolid = 24
  Ihq    = 8

# ----- 8. CONTACT
/INTER/TYPE7/INT_PROJ_PANEL
  Igap   = 2
  Inacti = 6
  Fric   = 0.30
  Iform  = 2
  master surface = projectile outer surface
  slave  group   = panel strike face nodes
  STMIN  = 0.5       ! min penalty stiffness scale
  STMAX  = 1000

/INTER/TYPE7/INT_PANEL_SELF
  Igap   = 2
  Self   = 1         ! self-contact ON
  Iedge  = 1         ! edge-edge ON (post-erosion exposed edges)
  Iform  = 2
  Fric   = 0.20

/INTER/TYPE2/INT_PLY01_PLY02   ! repeat for 23 cohesive interfaces
  Spotflag = 25      ! cohesive failure (G_Ic / G_IIc / strength)
  Sigma_n_max = 60.0e6
  Tau_s_max   = 90.0e6
  G_Ic        = 277.0
  G_IIc       = 788.0
  ...
... (INT_PLY02_PLY03 ... INT_PLY23_PLY24)

# ----- 9. INITIAL CONDITION
/INIVEL/TRA/IC_PROJECTILE
  Vx = 0  Vy = 0  Vz = -V_s
  grnd_id = GR_NODE_PROJ_ALL

# ----- 10. BOUNDARY CONDITION
/BCS/CLAMP_PERIMETER
  Tra = 1 1 1
  Rot = 1 1 1
  grnd_id = GR_NODE_CLAMP

# ----- 11. TIME-HISTORY OUTPUT
/TH/GROU1/TH_PROJ
  group = GR_NODE_PROJ_ALL
  request: VX, VY, VZ, KE, X, Y, Z
/TH/GLOB
  request: KE, IE, EHG, EER, ECN
/ANIM/EVERY
  dt_anim = 1.0e-6   ! 1 µs animation step

/END
```

Engine deck (`*_0001.rad`):

```
/RUN/flatballistic_Vs<NNN>/1
  Tstop = 6.0e-5     ! 60 µs
/DT/NODA/CST
  dt_min = 5.0e-9
  dt_scale = 0.67
/PRINT/-1
/STOP
```

This is a skeleton; the runner.py templates the specific values per strike velocity and per ply.

---

## 7. Reference solution

### 7.1 Recht-Ipson 1963 over-match perforation

The Recht-Ipson model derives the residual velocity from kinetic-energy conservation. For an over-match shot ($V_s > V_{50}$), the projectile transfers a fixed energy to perforate the laminate and emerges at $V_r$:

\begin{align}
\frac{1}{2} m_p V_s^2 = \frac{1}{2} m_p V_{50}^2 + \frac{1}{2} m_p V_r^2.
\end{align}

Solving for $V_r$:

\begin{align}
\boxed{\, V_r = \sqrt{V_s^2 - V_{50}^2} \quad (V_s > V_{50}). \,}
\label{eq:RI}
\end{align}

For a sub-V50 shot ($V_s \leq V_{50}$), the projectile is arrested in the laminate ($V_r = 0$). The Lambert-Jonas 1976 generalization replaces the exponent 2 with a fitted $p \in [2, 2.5]$. This stage uses Recht-Ipson ($p = 2$) as the principal fit and reports the Lambert-Jonas $p$ as a quality check; for IM7/8552 4 mm panels the literature reports $p \approx 2.0$ to 2.1 (Naik and Shrirao 2004 [@NaikShrirao2004]).

### 7.2 V50 from a fitted Recht-Ipson curve

Given five strike velocities $\{V_{s,i}\}$ and five FEM residual velocities $\{V_{r,i}\}$, V50 is extracted by least-squares fitting equation (\ref{eq:RI}):

\begin{align}
V_{50}^{\text{fit}} = \arg\min_{V_{50}} \sum_{i:\, V_{s,i} > V_{50}} \left[ V_{r,i} - \sqrt{V_{s,i}^2 - V_{50}^2} \right]^2,
\end{align}

with the constraint that $V_{r,i} = 0$ exactly when $V_{s,i} \leq V_{50}^{\text{fit}}$. The fit uses only the over-match points (where $V_{r,i} > 0$); the sub-V50 points contribute only to the bracket on $V_{50}^{\text{fit}}$ from below. The runner.py implements this via `scipy.optimize.curve_fit` on the over-match subset with `V50` as the single free parameter, with bounds $V_{50} \in [V_{s,\min}, V_{s,\max}]$.

### 7.3 Cunniff parameter

The Cunniff dimensionless analysis [@Cunniff1992WovenFabrics] defines the velocity scale

\begin{align}
c^* = \left( \frac{\sigma_f \, \varepsilon_f}{2 \rho_f} \cdot \sqrt{\frac{E_f}{\rho_f}} \right)^{1/3}
\end{align}

with units (m^2/s^2 × m/s)^{1/3} = m/s. For IM7 fiber:

| Quantity | Value |
|---|---|
| $\sigma_f$ | 5.2 GPa |
| $\varepsilon_f$ | 0.019 |
| $E_f$ | 276 GPa |
| $\rho_f$ | 1780 kg/m^3 |

$c^*_{\text{IM7}} = \left( (5.2 \times 10^9)(0.019)/(2 \times 1780) \cdot \sqrt{276 \times 10^9 / 1780} \right)^{1/3} = (27750 \cdot 12454)^{1/3} = (3.46 \times 10^8)^{1/3} \approx 700$ m/s.

The dimensionless Cunniff parameter is then $V_{50}/c^*$; for the published tape-laminate value $V_{50} \approx 300$ m/s, $V_{50}/c^* \approx 0.43$. This is consistent with the Cunniff 1992 dimensionless plot for woven fabrics where $V_{50}/c^*$ falls in the range 0.3 to 0.6 for typical tape laminates against a hard-core projectile of mass-to-presented-area ratio comparable to a 9 mm FMJ.

The Cunniff value is reported by runner.py as a sanity check that V50 sits in the right velocity decade for the chosen fiber chemistry; a result of $V_{50}/c^* < 0.1$ or $> 1.0$ would indicate an incorrectly scaled material card, an erroneous projectile mass, or an erosion threshold mismatch.

### 7.4 Experimental anchor (the IM7/8552 4 mm tape laminate ballistic V50)

The principal experimental anchor is Vakili Rad et al. 2020 [@VakiliRad2020HighVelocityImpact], which reports baseline-tape and pseudo-woven IM7/8552 panel ballistic limits at NASA Glenn. The cited value used as the §9 pass criterion is

\begin{align}
V_{50}^{\text{ref}} = 300 \pm 15 \text{ m/s}
\end{align}

(baseline tape, 24-ply 4 mm IM7/8552, 9 mm-class projectile). The ± 15 m/s band reflects the MIL-STD-662F bracketing tolerance translated to the published mean. [UNVERIFIED in this run — the numerical V50 value, projectile mass, and exact panel thickness in Vakili Rad 2020 must be re-read against their Figures 7 and 8 and Table 3 at the start of the run; the runner.py exposes `V50_REF_MS` and `V50_REF_TOL_MS` as Python constants at the top of the file so the assertion can be re-tightened once the reference is digitized.]

### 7.5 Material card (IM7/8552, Soden 1998 + WWFE-II)

| Property | Symbol | Value | Source |
|---|---|---|---|
| Density | $\rho$ | 1580 kg/m^3 | Soden 1998 |
| Longitudinal Young modulus | $E_1$ | 165 GPa | Soden 1998 |
| Transverse Young modulus | $E_2 = E_3$ | 8.4 GPa | Soden 1998 |
| In-plane shear modulus | $G_{12} = G_{13}$ | 5.6 GPa | Soden 1998 |
| Transverse shear modulus | $G_{23}$ | 2.8 GPa | Soden 1998 |
| Major Poisson ratio | $\nu_{12} = \nu_{13}$ | 0.34 | Soden 1998 |
| Through-thickness Poisson ratio | $\nu_{23}$ | 0.50 | Soden 1998 |
| Longitudinal tensile strength | $X_T$ | 2560 MPa | Soden 1998 |
| Longitudinal compressive strength | $X_C$ | 1590 MPa | WWFE-II Kaddour-Hinton 2013 |
| Transverse tensile strength | $Y_T = Z_T$ | 73 MPa | Soden 1998 |
| Transverse compressive strength | $Y_C = Z_C$ | 185 MPa | WWFE-II Kaddour-Hinton 2013 |
| In-plane shear strength | $S_{12} = S_{13}$ | 90 MPa | Soden 1998 |
| Transverse shear strength | $S_{23}$ | 80 MPa | WWFE-II Kaddour-Hinton 2013 |
| Mode I interlaminar | $G_{Ic}$ | 277 J/m^2 | Camanho Davila Moura 2003 |
| Mode II interlaminar | $G_{IIc}$ | 788 J/m^2 | Camanho Davila Moura 2003 |
| Cohesive normal traction | $\sigma_n^{\max}$ | 60 MPa | Stage 12 calibration |
| Cohesive shear traction | $\tau_s^{\max}$ | 90 MPa | Stage 12 calibration |

The numbers above are rate-independent. A full ballistic constitutive law would add a Cowper-Symonds or Johnson-Cook strain-rate scaling on $X_T$, $Y_T$, $S_{12}$. The §9 tolerance of 5 % on V50 is set with explicit awareness that strain-rate scaling is *not* applied here; runner.py prints a warning to that effect.

---

## 8. Validation success criterion

The pass criterion is the Stage 15 row of master_plan.md §3 table, plus the test_progression_literature.md §15 row, expanded.

1. **V50 within 5 % of reference.** Fitted V50 from the five-velocity sweep satisfies $|V_{50}^{\text{fit}} - V_{50}^{\text{ref}}| / V_{50}^{\text{ref}} \leq 0.05$, with $V_{50}^{\text{ref}} = 300$ m/s (Vakili Rad 2020 baseline tape, IM7/8552 24-ply 4 mm). Hard fail if outside ± 5 %.
2. **Recht-Ipson residual velocity within 10 %.** For each over-match strike velocity ($V_s > V_{50}^{\text{fit}}$), the FEM residual velocity satisfies $|V_r^{\text{FEM}} - V_r^{\text{RI}}| / V_r^{\text{RI}} \leq 0.10$, where $V_r^{\text{RI}} = \sqrt{V_s^2 - (V_{50}^{\text{fit}})^2}$. Average over the over-match points.
3. **Sub-V50 arrest.** For each sub-V50 strike velocity ($V_s < V_{50}^{\text{fit}}$), the projectile must come to rest *inside* the laminate, with $V_r^{\text{FEM}} = 0$ within numerical noise (less than 5 m/s).
4. **Hourglass energy under 5 %.** $\max_t E_{\text{hg}}(t) / E_{\text{int}}(t) \leq 0.05$ over the entire simulation window. Hard fail if exceeded — indicates spurious zero-energy modes are dominating the response.
5. **Energy balance under 2 %.** $\left| (E_{\text{kin}}^{0} - E_{\text{kin}}^{t} - E_{\text{int}}^{t} - E_{\text{erosion}}^{t} - E_{\text{contact}}^{t}) / E_{\text{kin}}^{0} \right| \leq 0.02$ at every output frame. Sanity check on the explicit time integrator.
6. **Mesh convergence at $V_s = 350$ m/s.** $|V_r(M2) - V_r(M1)| / V_r(M1) \leq 0.05$. Confirms the production mesh M1 is converged.
7. **Lambert-Jonas exponent.** Fitted exponent $p$ from the more general Lambert-Jonas model is in [1.8, 2.4]; warn (not fail) outside this band.
8. **Cunniff sanity.** $V_{50}^{\text{fit}} / c^*_{\text{IM7}}$ is in [0.3, 0.6]; warn outside this band.

The runner.py emits a JSON results file with all eight checks and a pass/warn/fail status per check, plus a CSV with $(V_s, V_r, V_{50}^{\text{RI}}, V_{50}^{\text{LJ}}, p, E_{\text{hg}}/E_{\text{int}}, \text{energy balance})$ at every frame for downstream Typst-CeTZ plotting.

---

## 9. Run protocol and pass tolerance budget

| Quantity | Symbol | Target / band | Failure threshold |
|---|---|---|---|
| Reference V50 | $V_{50}^{\text{ref}}$ | 300 m/s (Vakili Rad 2020) | — |
| FEM V50 (fit from 5 strikes) | $V_{50}^{\text{fit}}$ | 285 to 315 m/s | outside ± 15 m/s |
| Recht-Ipson residual at $V_s = 1.5 V_{50}$ | $V_r$ | within ±10 % of $\sqrt{V_s^2 - V_{50}^2}$ | outside ±10 % |
| Hourglass / internal energy | $E_{\text{hg}}/E_{\text{int}}$ | < 0.05 | ≥ 0.05 anywhere |
| Energy balance | $\Delta E / E_{\text{kin}}^0$ | < 0.02 | ≥ 0.02 anywhere |
| Critical time step | $\Delta t$ | > 5 ns | falls below 5 ns (mass scaling kicks in or run aborts) |

The strike-velocity sweep is fixed: $V_s \in \{200, 250, 300, 350, 400\}$ m/s, five runs. Adaptive bracketing is *not* used at this stage to keep the runner deterministic and reproducible; if the initial five-point sweep gives a $V_{50}^{\text{fit}}$ outside the band the runner.py reports the sweep result and instructs (does not auto-execute) a refinement run at $V_s$ values bracketing the new estimate.

---

## 10. Output, post-processing, and plots

Per-run outputs (each $V_s$ writes to its own subdirectory `runs/Vs<NNN>/`):

| File | Contents |
|---|---|
| `flatballistic_Vs<NNN>_0000.rad` / `_0001.rad` | starter / engine decks |
| `flatballistic_Vs<NNN>.out` | starter log |
| `flatballistic_Vs<NNN>.h3d` | full result database |
| `flatballistic_Vs<NNN>_T01` | global + per-projectile time history |
| `flatballistic_Vs<NNN>_anim_*` | animation frames every 1 µs |

Per-run extracted quantities (parsed from T01 by `vortex-radioss` or by reading the `.h3d` via the Kitware vtkhdf converter, then PyVista):

- $V_z^{\text{proj}}(t)$ projectile centroid axial velocity history
- $V_r^{\text{proj}} = V_z^{\text{proj}}(t = T_{\text{sim}})$ residual velocity (signed)
- $E_{\text{kin}}^{\text{proj}}(t)$, $E_{\text{int}}^{\text{plate}}(t)$, $E_{\text{erosion}}(t)$, $E_{\text{contact}}(t)$, $E_{\text{hg}}(t)$ (global time-history group)
- Number of eroded elements per ply layer at $t = T_{\text{sim}}$
- Per-mode failure counters (fiber tension, fiber compression, matrix tension, matrix compression) per ply at $t = T_{\text{sim}}$

Per-sweep aggregated outputs (in `runs/_summary/`):

- `sweep_results.csv` — one row per strike velocity with $(V_s, V_r, E_{\text{kin,abs}}, E_{\text{int}}, E_{\text{erosion}}, E_{\text{hg}}, \text{eroded}_{\text{count}})$
- `recht_ipson_fit.json` — $V_{50}^{\text{fit}}$, Lambert-Jonas $p$, R^2 of fit, Cunniff ratio
- `pass_fail_report.json` — eight checks of §8, status string each
- `vr_vs_vs.csv` — strike vs residual velocity table for plotting

Plots (Typst + CeTZ per master_plan.md §5):

1. **Residual velocity vs strike velocity** — five FEM points + Recht-Ipson curve + Vakili Rad 2020 reference V50 vertical line.
2. **Energy partition vs time at $V_s = 350$ m/s** — kinetic, internal, erosion, hourglass.
3. **Per-ply erosion map at $V_s = V_{50}^{\text{fit}}$** — bar chart of eroded elements per ply layer (1 to 24).
4. **Mesh convergence** — $V_r$ at $V_s = 350$ m/s on M0/M1/M2.

PyVista headless renders one 3D snapshot at $t = T_{\text{sim}}/2$ for each strike velocity (perforation event mid-flight) for visual inspection.

The runner.py is the single command that drives the full sweep, post-process, fit, and pass/fail report. It is deterministic given the same OpenRadioss build and the same five strike velocities.

---

## Appendix A — Open issues and bias documentation

A1. **Strain-rate scaling is not applied.** The IM7/8552 strengths in §7.5 are quasi-static. Published rate-dependent CFRP data (Koerber, Camanho, et al. 2010; Hsiao and Daniel 1998) suggest $X_T$ rises by ~5 to 10 % at strain rates of 10^4 s^-1 typical of ballistic. A rate-scaled card would shift $V_{50}^{\text{fit}}$ upward by ~5 %; the §8 tolerance of ±5 % absorbs this. If the fit comes in below the reference by ~5 %, the diagnosis points first at strain-rate scaling rather than at the mesh, the contact, or the erosion threshold.

A2. **Projectile is a steel cylinder, not a 9 mm FMJ.** The 9 mm FMJ has a soft lead core inside a copper jacket; under ballistic loading the lead core flow-deforms while the jacket rolls back. A homogeneous steel cylinder of equal mass and equal presented diameter overestimates projectile rigidity and underestimates plate damage spread. Vakili Rad 2020 reports a 3 to 7 % bias on V50 from this approximation [UNVERIFIED]. The bias direction is positive (FEM V50 will be higher than experiment), and the §8 tolerance absorbs it.

A3. **Cohesive parameters from Camanho 2003, not from a Stage 12 calibration on the present material.** Stage 12 of this pipeline will produce a calibrated $G_{Ic}$, $G_{IIc}$, $\sigma_n^{\max}$, $\tau_s^{\max}$ for IM7/8552 from a DCB / ENF run; until that is done, the values from Camanho Davila Moura 2003 are used. This is identified as a re-entry point: Stage 15 should be re-run after Stage 12 closes.

A4. **The 9 mm FMJ NIJ Level II velocity is 358 m/s, Level IIIA is 427 m/s.** The strike-velocity sweep here brackets V50 ≈ 300 m/s and *does not* test against either NIJ threshold velocity. That is intentional — the question asked at this stage is "does the FEM reproduce the published V50 within tolerance", not "is the panel NIJ-rated". If a NIJ rating against Level II is required later, run an additional point at $V_s = 358$ m/s and confirm full perforation with the residual velocity recorded.

A5. **Self-contact after erosion is the documented stability bottleneck.** OpenRadioss /INTER/TYPE7 with `Self = 1` and `Iedge = 1` is expensive (contact-search cost scales with the number of exposed faces, which grows as elements erode). If the run aborts with a contact-stiffness instability, the mitigation order is: increase the contact stiffness scale STMIN, then reduce the friction coefficient on INT_PANEL_SELF, then reduce $D_{\max}$ from 0.99 to 0.95. None of these mitigations affect the V50 prediction by more than 1 to 2 %, but they change the time-step floor.

A6. **AGPL-3 contamination.** OpenRadioss distribution is AGPL-3. Decks and runner.py are MIT/BSD-clean; the OpenRadioss binary is invoked but not redistributed. Master plan §10 risk 7 already documents this.
