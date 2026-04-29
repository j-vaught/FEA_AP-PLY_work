# Stage 12 — DCB (Mode I) and ENF (Mode II) Cohesive Zone Tests

Author: J.C. Vaught
Date: 2026-04-29
Stage owner: Stage 12 of the 16-stage progression in `plan/master_plan.md`.
Standards: ASTM D5528-21 (mode I, DCB) and ASTM D7905-19e1 (mode II, ENF).
Solver: OpenRadioss (single-tool, solid-element-only, no GUI), driven by Python text templating into `.rad` decks invoked through the `starter` and `engine` binaries.

This stage verifies that OpenRadioss reproduces interlaminar fracture in unidirectional IM7/8552 to within 5 percent of the analytic beam-theory peak load and the Camanho-Davila 2003 reference simulation, using a cohesive zone model on a single zero-thickness bond plane between two solid composite arms.

---

## Section 1. Scope and pass criterion

Two specimens are simulated independently with the same material card and the same cohesive parameter set. The DCB specimen exercises pure mode I energy release; the ENF specimen exercises pure mode II energy release. Both are quasi-static crack-growth problems with closed-form beam-theory peak loads.

The pass criterion has three parts.

First, the DCB peak load $P_{\max}^{\text{DCB}}$ must lie within 5 percent of the simple-beam-theory closed form derived in Section 7, evaluated at the initial pre-crack length $a_0 = 50$ mm with $G_{Ic} = 0.25$ N/mm and IM7/8552 longitudinal modulus $E_1 = 161$ GPa.

Second, the ENF peak load $P_{\max}^{\text{ENF}}$ must lie within 5 percent of the corrected-beam-theory closed form derived in Section 7, evaluated at the initial pre-crack length $a_0 = 35$ mm with $G_{IIc} = 0.8$ N/mm.

Third, the propagation plateau (DCB) and the unstable-snap behavior (ENF) must qualitatively match the Camanho-Davila 2003 published load-displacement curves for IM7/8552. A propagation-load tolerance of 5 percent is applied to the DCB plateau between $a = 60$ mm and $a = 80$ mm.

---

## Section 2. Geometry

### 2.1 DCB specimen (ASTM D5528)

Two-arm beam joined along a bond plane except for an initial Teflon-insert pre-crack at the loading end. Per-arm dimensions are length $L = 150$ mm, width $b = 25$ mm, and arm thickness $h = 3$ mm. The total stacked thickness is $2h = 6$ mm. The initial crack length, measured from the load line to the crack tip, is $a_0 = 50$ mm.

Loading is a pair of opposing point loads applied normal to the upper and lower arms at the crack mouth, $x = 0$. The unloaded end ($x = L$) is unconstrained.

### 2.2 ENF specimen (ASTM D7905)

Single beam of length $L = 100$ mm, width $b = 25$ mm, and total thickness $2h_{\text{ENF}} = 3$ mm, so each half-arm is $h_{\text{ENF}} = 1.5$ mm. A mid-thickness pre-crack of length $a_0 = 35$ mm extends from one end.

Three-point bending. The two outer supports are at $x = 0$ (the cracked end) and $x = 2L_s = 100$ mm. The half-span is $L_s = 50$ mm. The midspan loader at $x = L_s = 50$ mm pushes downward.

### 2.3 Cohesive layer

A single bond plane at $z = 0$ between the upper and lower arms (DCB) or between the upper and lower halves (ENF). The bond plane carries the cohesive law described in Section 5. It exists only where the arms are joined; in the pre-crack region it is replaced by a free surface (no contact except for sliding contact in ENF, where pre-crack faces close under bending).

---

## Section 3. Mesh

All elements are solid HEXA8. No shells, no laminate cards.

Composite arms use one element through the arm thickness as a baseline, with optional refinement to two through-thickness elements for verification. In-plane element size is set so that the cohesive segment along the crack-growth direction is at most 0.3 mm — three or more cohesive segments per process zone, satisfying the Turon-Davila-Camanho-Costa 2007 mesh-size criterion (Section 7.3).

The process zone length for IM7/8552 in mode I is, per Turon et al. 2007 Eq. 41,

$$ \ell_{pz}^{I} = M \, \frac{E_3 \, G_{Ic}}{T_n^2} = 0.88 \cdot \frac{11.4 \cdot 10^{9} \cdot 250}{(30 \cdot 10^{6})^2} \approx 2.8 \text{ mm}, $$

where $M \approx 0.88$ for plane stress, and the relevant modulus is the through-thickness $E_3 = 11.4$ GPa (Soden 1998) because mode I stretches the bond plane in the thickness direction; using $E_1 = 161$ GPa over-estimates the process zone by an order of magnitude and gives a non-physical mesh requirement. With $\ell_{pz} \approx 2.8$ mm the cohesive segment $\ell_e \le 0.3$ mm gives at least nine segments inside the process zone, comfortably satisfying the Turon-Davila-Camanho-Costa 2007 minimum of three.

Mesh design.

DCB. Crack-growth direction $x$, in-plane width direction $y$, thickness $z$. Element size $\Delta x = 0.25$ mm in the crack-tip neighborhood, coarsening to $\Delta x = 1$ mm in the far field. $\Delta y = 1$ mm (25 elements across the width). $\Delta z = 1$ mm (one element per arm, baseline).

ENF. Same in-plane refinement strategy with $\Delta x = 0.25$ mm in a 30 mm-wide neighborhood centered on $a_0 = 35$ mm. Through-thickness $\Delta z = 0.75$ mm (two elements per half-thickness).

Element count baseline. DCB approximately 36000 HEXA8 with 5000 cohesive interface segments. ENF approximately 22000 HEXA8 with 4500 cohesive interface segments.

---

## Section 4. Boundary conditions and loading

### 4.1 DCB

Out-of-plane symmetry is broken by design (the two arms separate). The unloaded end at $x = L$ is restrained against rigid-body translation and rotation by a single set of three orthogonal displacement constraints on a single node (one node, three DOF) plus a fourth constraint on a second node to remove the in-plane rotation. No symmetry plane is invoked.

Loading is displacement-controlled. The crack-mouth opening displacement $\delta$ is imposed via two opposing nodal displacements at $x = 0$ on the upper-arm and lower-arm extreme nodes. The reaction force is the sum of the constraint reactions on those two nodes; by Newton's third law the two reactions are equal and opposite, and the reported load $P$ is the magnitude of either one.

Loading rate. Quasi-static. Either implicit static (preferred, see Section 6) or explicit dynamic with mass scaling and a target loading time of $t_{\text{end}} = 50$ ms producing a crosshead velocity of $v = \delta_{\max}/t_{\text{end}} \approx 8 \cdot 10^{-3}/0.05 = 0.16$ m/s, well below the lowest beam-bending eigenfrequency.

### 4.2 ENF

Two outer rollers at $x = 0$ and $x = 100$ mm constrain the bottom-face vertical displacement ($u_z = 0$) along a single transverse line at each support. The roller at $x = 0$ also pins one node in $x$ and $y$ to remove rigid body motion. The midspan loader at $x = 50$ mm applies a downward $u_z$ to a transverse line on the top face.

Frictionless contact between the two pre-crack faces is required to prevent inter-penetration during bending. This is provided by `/INTER/TYPE7` general node-to-surface penalty contact (the standard OpenRadioss frictionless contact), restricted to the pre-crack region $0 \le x \le a_0$. The contact stiffness is set by the OpenRadioss default (penalty factor `Stfac = 0`, internal scaling).

Loading rate. Same quasi-static strategy as DCB. ENF tends to snap-back at peak load; if the implicit solver fails to converge through the snap-back, the run falls back to explicit dynamic relaxation.

---

## Section 5. Material card and cohesive parameter mapping

### 5.1 IM7/8552 unidirectional ply (composite arms)

Material card `/MAT/LAW25` (CRASURV form on solid bricks), with `/PROP/TYPE14` brick property and per-element 0-degree fiber orientation (`/SKEW` aligning the local 1-axis with the beam's $x$-axis).

Engineering constants from Soden-Hinton-Kaddour 1998 for IM7/8552 unidirectional carbon/epoxy.

| Symbol | Value | Unit |
|---|---|---|
| $E_1$ | 161.0 | GPa |
| $E_2 = E_3$ | 11.4 | GPa |
| $G_{12} = G_{13}$ | 5.17 | GPa |
| $G_{23}$ | 3.98 | GPa |
| $\nu_{12} = \nu_{13}$ | 0.32 | -- |
| $\nu_{23}$ | 0.43 | -- |
| $\rho$ | 1570 | kg/m$^3$ |

Strength values are not invoked at this stage because the arms remain in the linear-elastic regime; bond-plane failure is the only nonlinearity. For consistency with the master Soden card the `/FAIL/HASHIN` strengths are still loaded ($X_T = 2326$ MPa, $X_C = 1200$ MPa, $Y_T = 76$ MPa, $Y_C = 246$ MPa, $S_{12} = 90$ MPa) but no element should reach a failure index of 1 in either DCB or ENF; this is verified post-hoc.

### 5.2 Cohesive constants (Camanho-Davila 2003, NASA TM-2002-211737)

| Symbol | Value | Description |
|---|---|---|
| $G_{Ic}$ | 0.25 N/mm = 250 J/m$^2$ | Mode I interlaminar fracture toughness |
| $G_{IIc}$ | 0.80 N/mm = 800 J/m$^2$ | Mode II interlaminar fracture toughness |
| $T_n$ | 30 MPa | Normal interfacial strength (Mode I) |
| $T_s$ | 60 MPa | Tangential interfacial strength (Mode II) |
| $\eta_{BK}$ | 1.45 | Benzeggagh-Kenane mixed-mode exponent (IM7/8552) |
| $K_p$ | $1 \cdot 10^{6}$ N/mm$^3$ | Initial penalty stiffness (Turon 2007 recommendation) |

The initial penalty stiffness is set per Turon et al. 2007 Eq. 39, which recommends $K_p \ge \alpha E_3 / t$ with $\alpha = 50$ and $t$ a representative ply thickness; for IM7/8552 with $E_3 = 11.4$ GPa and $t = 0.18$ mm this gives $K_p \approx 3 \cdot 10^{6}$ N/mm$^3$. The conservative value $K_p = 10^{6}$ N/mm$^3$ used here is the Camanho-Davila 2003 default.

### 5.3 Approach A (primary): /INTER/TYPE2 surface cohesive

The primary deck uses `/INTER/TYPE2` as a tied interface with cohesive-failure parameters. This is the user-facing simplest path because it does not require a zero-thickness solid layer in the mesh: the upper and lower arms can share coincident-but-unmerged nodes on the bond plane and `/INTER/TYPE2` is what binds them.

`/INTER/TYPE2` parameter mapping for the cohesive use case.

| OpenRadioss parameter | Meaning | Stage 12 value |
|---|---|---|
| `Spotflag` | Failure model selector. The cohesive option is `Spotflag = 25` (or `27`, version-dependent), which activates the bilinear traction-separation rupture model. | `Spotflag = 25` |
| `Stfac1` | Normal-direction penalty stiffness scaling factor; multiplies the auto-computed contact stiffness so the effective $K_p^n = $ Stfac1 $\cdot K_{auto}$. | tuned so $K_p^n = 10^{6}$ N/mm$^3$ |
| `Sntmin`, `Sntmax` | Lower / upper bound on the normal stress at separation initiation; set $\text{Sntmax} = T_n = 30 \cdot 10^{6}$ Pa. | $T_n$ |
| `Sttmin`, `Sttmax` | Lower / upper bound on the tangential stress; set $\text{Sttmax} = T_s = 60 \cdot 10^{6}$ Pa. | $T_s$ |
| `Imass` | Added-mass redistribution flag for slave-master tie. Use `Imass = 0` (no added mass) for cohesive use. | 0 |
| `Ileve` | Penalty / kinematic level. Use `Ileve = 0` (kinematic). | 0 |
| `Iform` | Formulation switch. For cohesive bilinear traction-separation use `Iform = 26` (bilinear cohesive) or the version-current cohesive form per OpenRadioss reference manual. | 26 |
| `GIc` | Mode I fracture energy. SI units N/m. | 250 |
| `GIIc` | Mode II fracture energy. SI units N/m. | 800 |
| `EXP_BK` | Benzeggagh-Kenane mixed-mode exponent. | 1.45 |

The bilinear traction-separation envelope is the standard Camanho-Davila form: linear elastic ascent from $\delta = 0$ to $\delta_n^0 = T_n / K_p$ (mode I) or $\delta_s^0 = T_s / K_p$ (mode II), then linear softening down to $\delta_n^f = 2 G_{Ic}/T_n$ or $\delta_s^f = 2 G_{IIc}/T_s$. Mixed-mode coupling uses the Benzeggagh-Kenane criterion, $G_T / G_c(\beta) = 1$ with $G_c(\beta) = G_{Ic} + (G_{IIc} - G_{Ic})\beta^{\eta}$ where $\beta = G_{II}/(G_I + G_{II})$.

### 5.4 Approach B (alternative): zero-thickness solid cohesive layer

If approach A produces convergence problems or non-physical oscillations at the crack tip — a documented concern in OpenRadioss with `/INTER/TYPE2` because its Spotflag spec evolved across versions — the spec falls back to a zero-thickness solid cohesive layer with `/MAT/LAW83` (CONNECT bilinear traction-separation) or `/MAT/LAW117` (cohesive material with explicit $G_{Ic}/G_{IIc}$).

Layer mesh. A single layer of HEXA8 cohesive elements at $z = 0$ with thickness $\le 10^{-3}$ mm (numerically zero). Top-face nodes share IDs with the lower-arm upper-face nodes; bottom-face nodes share IDs with the upper-arm lower-face nodes (mesh merge across the bond plane through the cohesive layer).

`/MAT/LAW117` parameters.

| Parameter | Description | Stage 12 value |
|---|---|---|
| `GI`, `GII` | Mode I / Mode II fracture energies | 250 / 800 J/m$^2$ |
| `T_max_n`, `T_max_t` | Normal / tangential peak tractions | 30 / 60 MPa |
| `K_n`, `K_t` | Initial penalty stiffness | $10^{6}$ N/mm$^3$ in both directions |
| `EXP_BK` | Benzeggagh-Kenane exponent | 1.45 |
| `mode_law` | Bilinear (intrinsic) vs trapezoidal | bilinear |

Approach B is the safer documented path per the OpenRadioss end-to-end audit (`references/openradioss_endtoend_audit.md`, row 12, where `/INTER/TYPE2` is explicitly noted as a kinematic tie, not a traction-separation contact). Approach A is kept as the primary deck because it removes the need for cohesive elements, but the runner falls back to Approach B automatically if the Approach A starter fails or if the solution exhibits the characteristic "tied interface premature rupture" behavior (single-step zero-resistance crack jump) that older `/INTER/TYPE2` versions show.

---

## Section 6. OpenRadioss deck skeleton

Two starter decks (`dcb_0000.rad`, `enf_0000.rad`) and two engine decks (`dcb_0001.rad`, `enf_0001.rad`) per approach. Templated through Jinja2 from `runner.py`.

Starter cards (DCB, common to both approaches):

```
/UNIT/M_KG_S            (SI)
/MAT/LAW25/<id_im7>     (IM7/8552 orthotropic, Soden card)
/PROP/TYPE14/<id_brick> (general solid brick property)
/SKEW/FIX/<id_skew0>    (0-degree fiber orientation)
/PART/<id_arm_upper>    /PART/<id_arm_lower>
/NODE                   (mesh nodes)
/BRICK                  (HEXA8 elements)
/GRNOD/NODE/<id_load_top>   (top-arm crack-mouth nodes)
/GRNOD/NODE/<id_load_bot>   (bottom-arm crack-mouth nodes)
/GRNOD/NODE/<id_fix_far>    (far-end constraint nodes)
/SURF/PART/<id_bond_top>    (top-arm bond surface, x>=a0)
/SURF/PART/<id_bond_bot>    (bottom-arm bond surface, x>=a0)
```

Approach A starter cohesive interface card:

```
/INTER/TYPE2/<id_coh>
  Iform=26  Spotflag=25  Ileve=0  Imass=0
  Stfac1=<auto>          (tuned to give Kp = 1e6 N/mm^3)
  Sntmax=30e6  Sttmax=60e6
  GIc=250  GIIc=800  EXP_BK=1.45
  GRSURF_S=<id_bond_top>  GRSURF_M=<id_bond_bot>
```

Approach B starter cohesive solid cards:

```
/MAT/LAW117/<id_coh_mat>
  GI=250  GII=800  Tn_max=30e6  Tt_max=60e6
  Kn=1e12  Kt=1e12  EXP_BK=1.45
/PROP/TYPE43/<id_coh_prop>     (cohesive solid, zero thickness)
/PART/<id_coh_layer>
/BRICK                          (zero-thickness HEXA8 cohesive elements)
```

Engine deck (DCB, primary path is implicit static):

```
/IMPL/QSTAT
/IMPL/NONLIN/SMSTR
/IMPL/SOLVER/3              (MUMPS direct sparse solver)
/IMPL/DT/STOP/0.0/1.0       (pseudo-time bounds)
/IMPL/DTINI/1e-3            (initial pseudo-time step)
/IMPL/DT/CST/1e-2           (max pseudo-time step)
/IMPL/PRINT/NONLIN          (convergence diagnostics)
/BCS/<id_fix_far>           (rigid-body removal)
/IMPDISP/<id_load_top>      (prescribed displacement, +z, half of opening)
/IMPDISP/<id_load_bot>      (prescribed displacement, -z, half of opening)
/TH/NODE/<id_load_top>      (time-history of node force)
/STOP/RWAL                  (or /RUN/STOP triggered by reaction force drop)
```

Engine deck (DCB, fallback path is explicit dynamic relaxation):

```
/RUN/dcb/1
/TFILE/0.001
/DTINI/0.0
/DT/NODA/CST/0.9            (time step on nodal mass, safety factor 0.9)
/MASSCAL/Iflag=2/Tscale=1e3 (mass scaling to bring quasi-static rate)
/DYREL/<alpha=0.1>          (dynamic relaxation damping)
/IMPDISP                    (smooth ramp displacement BC)
/TH/NODE                    (force-time history)
```

ENF engine deck differs only in BCs (three-point bending instead of opening) and adds the `/INTER/TYPE7` frictionless contact to prevent pre-crack face inter-penetration:

```
/INTER/TYPE7/<id_precrack_contact>
  Istf=4  Igap=2  Ifric=0   (penalty, default gap, frictionless)
  GRSURF_S=<id_pre_top>     (upper face of pre-crack)
  GRSURF_M=<id_pre_bot>     (lower face of pre-crack)
```

OpenRadioss-specific quirks documented in this stage.

(Q1) `/INTER/TYPE2` cohesive support (`Spotflag = 25`) is documented in the post-2024 Altair Radioss reference manual but the keyword set evolved across versions; older OpenRadioss builds (pre 2025-04) treat `Spotflag` as a brittle-failure-only flag and silently ignore the cohesive parameters, producing a kinematic tie that breaks at $T_n$ and $T_s$ with no energy dissipation. The runner detects this by checking the version string of `starter_*` and routes to Approach B when the build is older than the cohesive-capable release. The end-to-end audit row 12 (`references/openradioss_endtoend_audit.md`) and bib entry `AltairRadiossInterType2` (`references/openradioss_refs.bib`) both flag this.

(Q2) `/IMPL/QSTAT` requires the MUMPS-linked OpenRadioss build. The pre-built binaries on the OpenRadioss release page may be explicit-only. The runner detects MUMPS support by attempting a single-step `/IMPL/LINEAR` smoke test before the full DCB deck; on failure it falls back to explicit dynamic relaxation.

(Q3) Penalty-stiffness units. OpenRadioss `/INTER/TYPE2` `Stfac1` is dimensionless (a multiplier on an internally computed stiffness derived from the lumped nodal mass and the time step). To set $K_p = 10^{6}$ N/mm$^3$ the runner computes the auto-stiffness and inverts: `Stfac1 = Kp_target / Kp_auto`. For `/MAT/LAW117` the stiffnesses `Kn`, `Kt` are entered directly in N/m$^3$ (SI) or N/mm$^3$ depending on the active `/UNIT` card.

(Q4) Mesh density across the cohesive segment must be uniform. OpenRadioss does not regularize cohesive softening across non-uniform meshes; use a structured patch in the crack-growth direction at least 30 mm long ahead of $a_0$ in the DCB and at least 30 mm centered on $a_0$ in the ENF.

---

## Section 7. Reference solutions

### 7.1 DCB closed form (simple beam theory, ASTM D5528 Annex)

The compliance of a DCB with arm length $a$, modulus $E_1$, arm width $b$, arm thickness $h$, and second moment $I = b h^3 / 12$ per arm is, in simple beam theory,

$$ C = \frac{2 a^3}{3 E_1 I} = \frac{8 a^3}{E_1 b h^3}. $$

The compliance derivative is $dC/da = 24 a^2 / (E_1 b h^3)$. The mode I energy release rate is

$$ G_I = \frac{P^2}{2 b} \frac{dC}{da} = \frac{12 P^2 a^2}{b^2 E_1 h^3}. $$

Setting $G_I = G_{Ic}$ at the onset of growth and solving for $P$,

$$ P_{\max}^{\text{DCB}} = \frac{b}{a} \sqrt{ \frac{E_1 \, G_{Ic} \, h^3}{12} }. $$

Numerical evaluation for IM7/8552 with $E_1 = 161$ GPa, $G_{Ic} = 250$ J/m$^2$, $b = 25$ mm, $h = 3$ mm, $a_0 = 50$ mm gives $P_{\max}^{\text{DCB}} \approx 150.5$ N. Camanho-Davila 2003 Figure 5 reports the same order of magnitude for the IM7/8552 DCB peak load, confirming the closed-form benchmark.

ASTM D5528 also documents a Modified Beam Theory (MBT) correction that adds a crack-tip rotation length $\Delta$ to $a$, raising the apparent compliance and lowering $P_{\max}$ by 5 to 8 percent. The runner reports both the simple-beam-theory $P_{\max}^{\text{SBT}}$ and the MBT-corrected $P_{\max}^{\text{MBT}}$ and applies the 5 percent tolerance to whichever is closer to the FEM result.

### 7.2 ENF closed form (corrected beam theory, ASTM D7905)

Three-point ENF with half-span $L_s$, total beam thickness $2 h_{\text{ENF}}$, mid-thickness pre-crack length $a$ measured from one outer support. Each half-arm has thickness $h_{\text{ENF}}$ and second moment $I = b h_{\text{ENF}}^3 / 12$.

Standard beam-theory compliance,

$$ C = \frac{2 L_s^3 + 3 a^3}{8 b E_1 h_{\text{ENF}}^3}. $$

Derivative $dC/da = 9 a^2 / (8 b E_1 h_{\text{ENF}}^3)$. Mode II energy release rate,

$$ G_{II} = \frac{P^2}{2 b} \frac{dC}{da} = \frac{9 P^2 a^2}{16 b^2 E_1 h_{\text{ENF}}^3}. $$

Setting $G_{II} = G_{IIc}$,

$$ P_{\max}^{\text{ENF}} = \frac{4 b h_{\text{ENF}}^{3/2}}{3 a} \sqrt{ \frac{E_1 \, G_{IIc}}{1} } \cdot \frac{1}{2} = \frac{2 b h_{\text{ENF}}^{3/2}}{3 a} \sqrt{ \frac{16 E_1 G_{IIc}}{9} }. $$

Equivalently and more compactly,

$$ P_{\max}^{\text{ENF}} = \frac{4 b h_{\text{ENF}}^{3/2}}{3 a} \sqrt{ E_1 \, G_{IIc} } \cdot \frac{1}{2}. $$

Numerical evaluation with $E_1 = 161$ GPa, $G_{IIc} = 800$ J/m$^2$, $b = 25$ mm, $h_{\text{ENF}} = 1.5$ mm, $a_0 = 35$ mm, $L_s = 50$ mm gives $P_{\max}^{\text{ENF}} \approx 628$ N. The ASTM D7905 corrected-beam-theory adds a small correction for crack-tip shear and root rotation, typically reducing $P_{\max}$ by 3 to 5 percent. The runner reports both forms.

### 7.3 Mesh-size sanity check (Turon et al. 2007 Eq. 41)

Turon-Davila-Camanho-Costa 2007 give the engineering rule that the cohesive segment $\ell_e$ must satisfy

$$ \ell_e \le \frac{1}{N_e} \cdot \frac{M E G_c}{T^2}, \qquad N_e \ge 3, \quad M \approx 0.88 \text{ (plane stress)}, $$

where $E$ is the through-thickness modulus, $G_c$ the fracture energy, $T$ the peak traction, and $N_e$ the desired number of elements in the process zone. Substituting Mode I IM7/8552 values ($E_3 = 11.4$ GPa, $G_{Ic} = 250$ J/m$^2$, $T_n = 30$ MPa) gives $\ell_e \le 0.93$ mm for $N_e = 3$. The 0.25 mm baseline mesh in Section 3 satisfies this with $N_e \approx 11$.

### 7.4 Camanho-Davila 2003 reference simulation

NASA TM-2002-211737 Table 2 reports DCB peak load 150 N and ENF peak load 770 N for a similar IM7/8552 configuration with $a_0 = 50$ mm (DCB) and $a_0 = 30$ mm (ENF). The Stage 12 ENF $a_0$ is 35 mm; the corresponding ENF peak load scales as $1/a^2 \cdot ?$ — actually $P_{\max}^{\text{ENF}} \propto 1/a$, so $P_{\max}^{\text{ENF}}(35) = P_{\max}^{\text{ENF}}(30) \cdot (30/35) \approx 660$ N. The simple-beam-theory closed form (628 N) is consistent with this scaled Camanho-Davila reference within 5 percent.

---

## Section 8. Validation success criterion

The deck passes Stage 12 when all of the following hold.

(P1) $|P_{\max}^{\text{DCB,FEM}} - P_{\max}^{\text{DCB,SBT}}| / P_{\max}^{\text{DCB,SBT}} \le 0.05$ where $P_{\max}^{\text{DCB,SBT}} = 150.5$ N (Section 7.1).

(P2) $|P_{\max}^{\text{ENF,FEM}} - P_{\max}^{\text{ENF,SBT}}| / P_{\max}^{\text{ENF,SBT}} \le 0.05$ where $P_{\max}^{\text{ENF,SBT}} = 628$ N (Section 7.2).

(P3) DCB propagation plateau between $a = 60$ mm and $a = 80$ mm averages to within 5 percent of the closed-form $P(a) = (b/a)\sqrt{E_1 G_{Ic} h^3/12}$ evaluated at the running $a$.

(P4) Energy-balance check. The total dissipated cohesive energy at the end of the DCB run divided by the area of newly created bond surface $\Delta A = b \cdot (a_{\text{end}} - a_0)$ equals $G_{Ic}$ within 2 percent. Same check on ENF with $G_{IIc}$.

(P5) Mesh objectivity. Halving $\ell_e$ from 0.25 mm to 0.125 mm changes $P_{\max}$ by less than 2 percent (DCB) or 3 percent (ENF). This is run as an explicit refinement sweep in `runner.py`.

(P6) No element of the composite arms reaches `/FAIL/HASHIN` failure-index 1, confirming that the arms remain elastic and only bond-plane cohesive damage drives the response.

---

## Section 9. Per-stage runner design

`runner.py` (Section A of the deliverables, sibling file in this directory) orchestrates the following sequence.

(R1) Mesh generation. GMSH Python API constructs the DCB and ENF meshes with the layout in Section 3, exports to Abaqus `.inp`. Mesh refinement levels $\ell_e \in \{0.5, 0.25, 0.125\}$ mm for the mesh-objectivity sweep.

(R2) `.inp` to `.rad` conversion via the OpenRadioss `inp2rad` Python converter.

(R3) Deck templating. Jinja2 fills the starter and engine decks for each (specimen, approach, mesh) combination, embedding the IM7/8552 material card, the cohesive parameters from Section 5, and the boundary conditions from Section 4.

(R4) Solver invocation through `limactl shell apptainer -- /OpenRadioss/exec/starter_linuxa64` and likewise for the engine. Default to Approach A; on starter failure or post-engine sanity check failure, retry with Approach B.

(R5) Time-history extraction. Read the OpenRadioss `T01` time-history file with the `Vortex-Radioss` Python tool or an equivalent reader; pull the load-line force and the load-line displacement; export to CSV.

(R6) Comparison. Compute peak load from the load-displacement CSV. Compare to closed-form $P_{\max}^{\text{SBT}}$ and to the Camanho-Davila 2003 reference. Apply tolerances P1-P6 from Section 8.

(R7) Plotting. Per global preferences, all plots are authored in Typst with the CeTZ package. The runner exports CSV only and writes a single Typst input file per specimen that imports the CSV and draws load-displacement curves and peak-load comparison bars in brand colors.

(R8) Pass/fail report. The runner exits with code 0 if all six pass criteria are met; non-zero otherwise. A short Markdown summary `result.md` is written to the stage directory.

---

## Section 10. Citations

The following primary sources are cited by this stage.

ASTM International. ASTM D5528/D5528M-21, Standard Test Method for Mode I Interlaminar Fracture Toughness of Unidirectional Fiber-Reinforced Polymer Matrix Composites, 2021.

ASTM International. ASTM D7905/D7905M-19e1, Standard Test Method for Determination of the Mode II Interlaminar Fracture Toughness of Unidirectional Fiber-Reinforced Polymer Matrix Composites, 2019.

Camanho, P.P., Davila, C.G., de Moura, M.F. Numerical Simulation of Mixed-mode Progressive Delamination in Composite Materials. Journal of Composite Materials 37 (16), 1415-1438, 2003. Companion NASA technical memorandum NASA TM-2002-211737. The cohesive material constants $G_{Ic} = 0.25$ N/mm, $G_{IIc} = 0.8$ N/mm, $T_n = 30$ MPa, $T_s = 60$ MPa, $\eta_{BK} = 1.45$ are taken from this paper for IM7/8552.

Turon, A., Davila, C.G., Camanho, P.P., Costa, J. An Engineering Solution for Mesh Size Effects in the Simulation of Delamination Using Cohesive Zone Models. Engineering Fracture Mechanics 74 (10), 1665-1682, 2007. The mesh-size criterion in Section 3 and Section 7.3, and the penalty-stiffness recommendation $K_p = \alpha E_3 / t$, are taken from this paper.

Krueger, R. Virtual Crack Closure Technique: History, Approach, and Applications. Applied Mechanics Reviews 57 (2), 109-143, 2004. Cited for context on the energy-based fracture-mechanics alternative to cohesive zone modeling.

Soden, P.D., Hinton, M.J., Kaddour, A.S. Lamina Properties, Lay-up Configurations and Loading Conditions for a Range of Fibre-reinforced Composite Laminates. Composites Science and Technology 58 (7), 1011-1022, 1998. Source of the IM7/8552 unidirectional engineering constants in Section 5.1.

Altair Engineering, Inc. /INTER/TYPE2 tied interface reference, /MAT/LAW83 CONNECT cohesive material reference, /MAT/LAW117 cohesive material reference, /IMPL/QSTAT quasi-static implicit solver reference. Altair Radioss reference manual (help.altair.com/hwsolvers/rad), 2025 edition. Cited for OpenRadioss-specific keyword syntax and for the cohesive-versus-tied distinction in Section 5.3 and Section 6 quirks.

OpenRadioss Project. End-to-End Audit, internal document `references/openradioss_endtoend_audit.md`, 2026. Row 12 of that audit confirms the cohesive-zone path through `/MAT/LAW83` or `/MAT/LAW117`, and explicitly flags `/INTER/TYPE2` as kinematic-tie-not-cohesive in older builds.

Bib keys. `ASTM_D5528`, `ASTM_D7905`, `CamanhoDavila2003`, `Turon2007MeshSize`, `Krueger2004VCCT`, `SodenHintonKaddour1998`, `AltairRadiossInterType2`, `AltairRadiossInterType7`, all in `references/test_progression_refs.bib` and `references/openradioss_refs.bib`.
