# Stage 03 - Isotropic Dogbone Tensile Coupon (ASTM E8/E8M)

Author. J.C. Vaught
Date. 2026-04-29
Status. Specification, pre-run.

---

## 1. Stage summary

This stage isolates the elasto-plastic constitutive integration on a solid
hexahedral mesh and verifies that the FEM kernel reproduces uniform uniaxial
stress and strain in the gauge section of an ASTM E8/E8M-22 standard
rectangular plate-type tensile coupon, that the apparent Young's modulus
matches the input, and that yield onset occurs at the input yield stress.
It is the first stage that exercises a non-linear material law and the
first stage that imposes displacement control through a grip face on a
solid mesh. It builds on the verified linear-elastic kernel from stage 1
(linear-elastic beam) and the verified geometric-nonlinear kernel from
stage 2 (large-deflection cantilever). All elements are solid HEXA8 per
the master plan's hard solid-only constraint. The pass criterion is
quantitative on three independent measurements. The axial Cauchy stress
$\sigma_{xx}$ in the gauge section before yield must match the
applied-load-over-area ratio $P/A_{\text{gauge}}$ within $1\%$. The axial
engineering strain $\varepsilon_{xx}$ in the gauge section must match the
gauge-length-displacement ratio $\delta/L_{\text{gauge}}$ within $1\%$.
Yield onset, defined as the strain at which the apparent stiffness drops
below $0.99 E$, must occur at $\sigma_{xx} \in [0.99, 1.01] \sigma_y$.
Saint-Venant decay past the fillets and uniform stress-strain in the
constant-cross-section gauge are the qualitative observables that
underpin those quantitative checks. The reference solution is closed-form
in the elastic regime ($\sigma = E \varepsilon$, $\sigma = P/A$,
$\varepsilon = \delta/L$) and is taken from the input material card at
yield onset. Material data are A36 low-carbon steel from the master plan
brief. Standards citation throughout is ASTM E8/E8M-22 \cite{ASTM_E8}.

---

## 2. Geometry (ASTM E8/E8M-22, rectangular plate-type "Specimen 1")

ASTM E8/E8M-22 \cite{ASTM_E8} Table 1 defines six standard rectangular
plate-type tension specimens. The geometry below is the standard "plate
type, $40$ mm gauge length" specimen with the gauge length adjusted to
the canonical $L_0 = 50$ mm value used in the master plan brief and in
the SI-aligned ISO 6892-1 dogbone literature. Thickness is taken as the
master-plan-specified $t = 6$ mm, which is inside the E8/E8M plate-type
applicable range ($5 \le t \le 19$ mm for the standard plate type per
section 6.4 of E8/E8M-22) \cite{ASTM_E8}.

### 2.1 Dimension table

| Symbol | Quantity | Value | Source |
|---|---|---|---|
| $L_0$ | Gauge length | $50.0$ mm | master plan brief |
| $W$ | Reduced-section (gauge) width | $12.5$ mm | master plan brief |
| $t$ | Thickness | $6.0$ mm | master plan brief |
| $L_r$ | Reduced-section length | $60.0$ mm | $L_0 + 10$ mm per E8/E8M-22 sec.~6.4 |
| $R$ | Fillet radius (gauge to grip) | $12.5$ mm | E8/E8M-22 sec.~6.4 (min. $12.5$ mm) |
| $W_g$ | Grip-section width | $20.0$ mm | E8/E8M-22 sec.~6.4 (gauge $+$ $\ge 6$ mm) |
| $L_g$ | Grip-section length (each end) | $50.0$ mm | E8/E8M-22 sec.~6.4 |
| $L_T$ | Total specimen length | $L_r + 2(L_g + R) = 185.0$ mm | derived |
| $A_{\text{gauge}}$ | Gauge-section cross-section area | $W \cdot t = 75.0$ mm$^2$ | derived |

The transition from gauge to grip is a tangent circular fillet of radius
$R = 12.5$ mm, satisfying E8/E8M-22 section 6.4 minimum-radius
requirement that the fillet radius equal or exceed the gauge half-width.
The fillet is tangent to the parallel-section side edges and to the grip
side edges; this tangency is the standard E8/E8M sketch in Annex A1
Figure 1 of the document \cite{ASTM_E8}.

### 2.2 ASCII sketch (top view, in-plane, mm; not to scale)

```
                     (+y)
                      |
                      |
        x = -L_T/2    |      x = +L_T/2
            v         |          v
        +-------+     |      +-------+
        |       |\    |     /|       |
   W_g/2 -      | \   |    / |       - W_g/2
        | grip  |  \  |   /  | grip  |     <- y = +W_g/2
        |  L    |   \ |  /   |   R   |
        |       |  R \| /    |       |
        +-------+----/+\-----+-------+     <- y = +W/2  (gauge edge)
                | gauge L_0  |
                |  width W   |
                |  ----------|------> (+x)  axis of loading
                |            |
        +-------+----\+/-----+-------+     <- y = -W/2  (gauge edge)
        |       |   / | \    |       |
        |  L    |  / R|R \   |   R   |
        | grip  | /   |   \  | grip  |     <- y = -W_g/2
   W_g/2 -      |/    |    \ |       - W_g/2
        |       |     |     \|       |
        +-------+     |      +-------+
            ^         |          ^
        x = -L_T/2    |      x = +L_T/2
                      |
                     (-y)

         |<--L_g-->|<-R->|<--L_r-->|<-R->|<--L_g-->|
         |                                         |
         |<--------------- L_T = 185 mm -----------|

      Thickness t = 6 mm extruded out of the page along (+z).
```

The $x$-axis is the loading axis, $y$ is in-plane transverse, and $z$ is
through-thickness. The gauge length $L_0 = 50$ mm is centered on the
origin and extends from $x = -25$ mm to $x = +25$ mm. Gauge stress and
strain are measured in this region by the post-processing layer.

---

## 3. Mesh

The mesh uses solid HEXA8 elements (linear hexahedra, eight nodes,
trilinear shape functions, reduced or full integration per the property
card setting; default is full $2 \times 2 \times 2$ Gauss integration on
\texttt{/PROP/TYPE14}) \cite{AltairRadiossPropType14}. Zoning is graded so
that the gauge section is the most refined region and the grip sections
are coarsest, with a smooth transition through the fillets.

### 3.1 Target element counts

| Region | Through-thickness | Across width | Along length | Element size |
|---|---|---|---|---|
| Gauge section ($\|x\| \le L_0/2$) | $\ge 3$ (default $4$) | $\ge 6$ (default $8$) | $\ge 40$ (default $50$) | $\sim 1.0$ mm |
| Fillet transition | $4$ | $8$ tapering to $12$ | $\sim 8$ along arc | $\sim 1.5$ mm |
| Grip section ($L_0/2 + R \le \|x\| \le L_T/2$) | $4$ | $12$ | $\sim 25$ | $\sim 2.0$ mm |

This satisfies the master-plan-brief minima ($\ge 3$ through-thickness,
$\ge 6$ across width, $\ge 40$ along length). The $\ge 3$ through-thickness
requirement is what distinguishes a true 3D solid-element verification
from a plane-stress shortcut and is the load-bearing reason the master
plan permits no shell elements.

### 3.2 Mesh-convergence study

Three meshes are run at $h$, $h/\sqrt{2}$, and $h/2$ uniform refinement
in the gauge section. The coarse-mesh element count is approximately
$4 \times 8 \times 50 = 1600$ HEXA8 elements in the gauge alone, plus
about $5000$ elements in the fillet and grip regions, total
$\sim 6600$ elements coarse. The medium and fine meshes are
$\sim 13\,000$ and $\sim 26\,000$ elements respectively. All three meshes
must satisfy the pass criterion individually; the convergence rate of
the apparent modulus and gauge-stress error must be $O(h^2)$ for linear
hexahedra.

### 3.3 Mesh generation toolchain

GMSH (Python API) builds the geometry from the dimension table above
using a sequence of \texttt{addPoint}, \texttt{addLine},
\texttt{addCircleArc}, \texttt{addPlaneSurface},
\texttt{extrude} (along $z$ to thickness $t$) operations. The meshing
algorithm is \texttt{Frontal-Delaunay for Quads} on the in-plane faces,
followed by \texttt{Recombine} to produce a structured quad mesh,
followed by extrusion to produce HEXA8. The resulting mesh is exported
as Abaqus \texttt{.inp} via GMSH's native Abaqus writer; the
\texttt{.inp} is converted to OpenRadioss \texttt{.rad} via the
\texttt{inp2rad} Python converter \cite{OpenRadiossDiscussion3222INP2RAD}
shipped in the \texttt{OpenRadioss/Tools} repository
\cite{OpenRadiossTools2026}. Element sets (\texttt{*ELSET}) are written
for the gauge section, fillet, and grip sections so that material and
property assignments survive the conversion.

---

## 4. Boundary conditions and loading

Displacement control on the loaded grip face. The opposite grip face is
fully clamped. Lateral DOFs on the loaded grip face are suppressed to
avoid rigid-body rotation and twist; only the axial DOF moves. This is
the standard FEM idealization of a rigid wedge-grip in a screw-driven
test machine.

### 4.1 Boundary conditions

| Face | $x$-extent | $u_x$ | $u_y$ | $u_z$ | $\theta_{x,y,z}$ |
|---|---|---|---|---|---|
| Fixed grip (left) | $x = -L_T/2$ | $0$ | $0$ | $0$ | n/a (solid) |
| Loaded grip (right) | $x = +L_T/2$ | $\delta(t)$ (imposed) | $0$ | $0$ | n/a (solid) |
| All other faces | interior or external | free | free | free | n/a |

Solid HEXA8 elements have only translational DOFs, so rotational
restraints are not applicable. The lateral $u_y = u_z = 0$ constraints
on the loaded grip face are applied via \texttt{/BCS}
\cite{AltairRadiossBCS}. The clamp on the fixed grip face is also
\texttt{/BCS}. The axial displacement on the loaded face is applied via
\texttt{/IMPDISP} \cite{AltairRadiossImpDisp} driven by a piecewise-linear
function of time defined in \texttt{/FUNCT}.

### 4.2 Load history

The displacement is ramped linearly from $0$ to $\delta_{\max}$ over a
time interval $T_{\text{end}}$. Strain at the end of the ramp targets
$\varepsilon_{xx}^{\text{end}} = 0.005 = 0.5\%$, comfortably above the
elastic limit $\varepsilon_y = \sigma_y/E = 250 / 200{,}000 = 0.00125 =
0.125\%$, so the simulation passes through yield onset and a small
amount of post-yield plastic flow. Total imposed displacement is
$\delta_{\max} = \varepsilon^{\text{end}} \cdot L_{\text{eff}}$, where
$L_{\text{eff}}$ is the effective compliant length between grips
including fillet and partial-grip compliance; for a first pass we take
$L_{\text{eff}} \approx L_0 + 2 R = 75$ mm and so
$\delta_{\max} = 0.005 \cdot 75 = 0.375$ mm, then iterate once on the
basis of the first run if the actual gauge strain at $T_{\text{end}}$
overshoots or undershoots $0.5\%$. For the implicit \texttt{/IMPL/QSTAT}
solver, the time variable is pseudo-time and $T_{\text{end}}$ is set to
$1.0$ s.

---

## 5. Material card (LAW2 vs LAW36 trade-off)

Material is ASTM A36 low-carbon hot-rolled structural steel. Properties
from the master-plan brief are summarized below. The choice between
\texttt{/MAT/LAW2} (Johnson-Cook with rate effects off)
\cite{AltairRadiossLAW2} and \texttt{/MAT/LAW36} (PLAS\_TAB tabulated
piecewise-linear) \cite{AltairRadiossLAW36} is made deliberately and
documented here.

### 5.1 Input properties (A36 steel)

| Symbol | Quantity | Value |
|---|---|---|
| $E$ | Young's modulus | $200$ GPa $= 2.0 \times 10^{11}$ Pa |
| $\nu$ | Poisson's ratio | $0.30$ |
| $\rho$ | Density | $7850$ kg/m$^3$ |
| $\sigma_y$ | Yield stress (initial) | $250$ MPa |
| $\sigma_u$ | Ultimate tensile strength | $400$ MPa |
| $\varepsilon^{p}_u$ | Plastic strain at UTS (estimate) | $0.18$ |

Strain-rate effects are off because the loading is quasi-static and
implicit. Temperature effects are off because the test is isothermal.
Damage and failure are off because stage 3's pass criterion ends at
small post-yield plasticity, well below necking instability. Damage is
introduced in stage 5.

### 5.2 LAW2 (PLAS\_JOHNS) form

Johnson-Cook stress-strain law with rate and thermal terms switched
off reduces to the power-law hardening
$\sigma_y(\bar{\varepsilon}^p) = a + b (\bar{\varepsilon}^p)^n$, where $a$
is the initial yield, $b$ is the hardening modulus, and $n$ is the
hardening exponent. For A36 a standard fit is $a = 250$ MPa,
$b = 275$ MPa, $n = 0.36$, calibrated so that
$\sigma(\bar{\varepsilon}^p = 0.18) = 250 + 275 \cdot 0.18^{0.36}
\approx 400$ MPa, recovering the input UTS at the input plastic strain
at UTS. The rate sensitivity $c$ and reference plastic strain rate
$\dot{\varepsilon}_0$ are set so that
$1 + c \ln(\dot{\bar{\varepsilon}}^p / \dot{\varepsilon}_0) = 1$ at the
quasi-static rate of the test, equivalent to setting $c = 0$ which the
Altair documentation explicitly supports as the rate-effects-off branch
\cite{AltairRadiossLAW2}.

### 5.3 LAW36 (PLAS\_TAB) form

A user-defined function $\sigma_y(\bar{\varepsilon}^p)$ supplied as a
table of $(\bar{\varepsilon}^p, \sigma_y)$ pairs. For A36 the table is
$(0, 250)$, $(0.02, 290)$, $(0.05, 320)$, $(0.10, 360)$, $(0.18, 400)$
MPa, linearly interpolated between knots and held constant beyond
$\bar{\varepsilon}^p = 0.18$ for the small post-yield excursion of this
stage \cite{AltairRadiossLAW36}.

### 5.4 Trade-off and choice

| Criterion | LAW2 (Johnson-Cook, rate off) | LAW36 (PLAS\_TAB) |
|---|---|---|
| Form | Smooth power law $\sigma = a + b \varepsilon^n$ | Piecewise linear from a table |
| Number of parameters | $3$ ($a, b, n$) | $\ge 2$ knots in a table |
| Calibration | Three-parameter fit to handbook curve | Direct copy of handbook curve points |
| Reuse in stages 4-5 | Drop-in for any metallic plasticity | Same; table can be expanded |
| Failure / damage | Built-in Johnson-Cook damage available later | Requires \texttt{/FAIL} card |
| Rate dependence later | Re-enabled by setting $c \ne 0$, $\dot{\varepsilon}_0$ | Requires multiple curves indexed by rate |
| Documentation | \cite{AltairRadiossLAW2} | \cite{AltairRadiossLAW36} |

**Choice. \texttt{/MAT/LAW2} (Johnson-Cook with $c = 0$).** Reasoning.
First, LAW2 is a smooth analytic curve, so the elastic-to-plastic
transition is $C^0$ in stress and $C^{-1}$ only at the single point
$\bar{\varepsilon}^p = 0$, whereas LAW36 has a $C^{-1}$ slope
discontinuity at every interior knot of the table; the implicit Newton
iteration converges more reliably on the smoother LAW2 form for a
first-pass verification. Second, LAW2's three parameters $(a, b, n)$
correspond directly to the three numbers we have in the master-plan
brief ($\sigma_y$, $\sigma_u$, and an estimated $\varepsilon^p_u$),
whereas LAW36 begs the question of how many knots to put in the table
and where to place them. Third, when later stages (4, 5) reuse this
material card, LAW2's smooth analytic form is more portable. Fourth,
both cards are compatible with \texttt{/PROP/TYPE14} solid bricks per
the OpenRadioss reference manual \cite{AltairRadiossPropType14}, so the
choice has no impact on the property card. The runner.py is built
parameterized so that the material card can be flipped to LAW36 in a
single config flag, and the validation script is form-agnostic; this
preserves the option to re-run with LAW36 if a regression is found.

### 5.5 Material card excerpt (LAW2, SI-meter-second-kg unit set)

```
/MAT/LAW2/1
A36 steel quasi-static, rate effects off
# rho   E         nu     a        b        n      eps_max  sig_max
7850.0  2.0e11    0.3
2.5e8   2.75e8    0.36   0.0      0.0
# c     eps_dot_0  ICC  Fsmooth  Fcut  Chard  Tmelt  rho_Cp  Tref  ...
0.0     1.0        0    0        0     0      0      0       0
```

(SI base units: kg, m, s; stresses therefore in Pa.) The full deck
template lives in \texttt{runner.py} and uses Jinja2-style
substitutions to expand the values from the master-plan brief.

---

## 6. OpenRadioss deck skeleton

Two-file deck per the OpenRadioss starter / engine convention, with
filenames \texttt{stage03\_dogbone\_E8\_0000.rad} (starter) and
\texttt{stage03\_dogbone\_E8\_0001.rad} (engine).

### 6.1 Starter deck (\texttt{\_0000.rad}) keyword cards

| Keyword | Purpose | Reference |
|---|---|---|
| \texttt{/BEGIN} | Deck header, unit system, run title | OpenRadioss reference |
| \texttt{/MAT/LAW2/1} | A36 steel Johnson-Cook (rate off) | \cite{AltairRadiossLAW2} |
| \texttt{/PROP/TYPE14/1} | General solid property, no orthotropic frame | \cite{AltairRadiossPropType14} |
| \texttt{/PART/1} | Couples material to property to element set | OpenRadioss reference |
| \texttt{/NODE} | Nodal coordinates from inp2rad-converted GMSH mesh | \cite{OpenRadiossDiscussion3222INP2RAD} |
| \texttt{/BRICK} | HEXA8 connectivity | OpenRadioss reference |
| \texttt{/GRNOD/NODE/1} | Node group: fixed-grip face nodes | OpenRadioss reference |
| \texttt{/GRNOD/NODE/2} | Node group: loaded-grip face nodes | OpenRadioss reference |
| \texttt{/GRNOD/NODE/3} | Node group: gauge-section reference nodes (for $\varepsilon$ gauge) | OpenRadioss reference |
| \texttt{/GRBRIC/BRIC/1} | Element group: gauge-section bricks (for $\sigma$ gauge) | OpenRadioss reference |
| \texttt{/BCS/1} | Clamp on \texttt{GRNOD/1}: $u_x = u_y = u_z = 0$ | \cite{AltairRadiossBCS} |
| \texttt{/BCS/2} | Lateral on \texttt{GRNOD/2}: $u_y = u_z = 0$ | \cite{AltairRadiossBCS} |
| \texttt{/FUNCT/1} | Linear ramp $f(t) = t / T_{\text{end}}$ | OpenRadioss reference |
| \texttt{/IMPDISP/1} | Imposed displacement $u_x = \delta_{\max} \cdot f(t)$ on \texttt{GRNOD/2} | \cite{AltairRadiossImpDisp} |
| \texttt{/TH/NODE/1} | Output history on \texttt{GRNOD/1}: reaction $F_x$ | \cite{AltairRadiossTH} |
| \texttt{/TH/NODE/2} | Output history on \texttt{GRNOD/3}: gauge endpoint $u_x$ | \cite{AltairRadiossTH} |
| \texttt{/TH/BRICK/1} | Output history on \texttt{GRBRIC/1}: gauge $\sigma_{xx}$, $\bar{\varepsilon}^p$ | \cite{AltairRadiossTH} |
| \texttt{/END} | Deck terminator | OpenRadioss reference |

### 6.2 Engine deck (\texttt{\_0001.rad}) keyword cards

| Keyword | Purpose | Reference |
|---|---|---|
| \texttt{/RUN/...} | Run title and duration ($T_{\text{end}} = 1.0$ s) | OpenRadioss reference |
| \texttt{/IMPL/QSTAT} | Quasi-static implicit solver activation | \cite{AltairRadiossImplQstat} |
| \texttt{/IMPL/SOLVER} | MUMPS direct solver selection | \cite{AltairRadiossImplActivation} |
| \texttt{/IMPL/DT/STOP} | Initial pseudo-time step and stopping rule | \cite{AltairRadiossImplActivation} |
| \texttt{/IMPL/NONLIN} | Newton-Raphson nonlinear iteration controls | \cite{AltairRadiossImplActivation} |
| \texttt{/IMPL/PRINT} | Iteration log verbosity | \cite{AltairRadiossImplActivation} |
| \texttt{/ANIM/DT} | Animation output cadence | OpenRadioss reference |
| \texttt{/ANIM/BRICK/TENS/STRESS} | Per-element Cauchy stress tensor in animation | OpenRadioss reference |
| \texttt{/ANIM/BRICK/EPSP} | Per-element equivalent plastic strain in animation | OpenRadioss reference |
| \texttt{/ANIM/NODA/DISP} | Per-node displacement in animation | OpenRadioss reference |
| \texttt{/STOP} | Engine terminator | OpenRadioss reference |

This deck skeleton follows the OpenRadioss tensile-test Confluence
tutorial \cite{OpenRadiossConfluenceTensile2024} for keyword choice and
ordering, with the substitutions documented above for the A36 material,
the SI unit set, and the explicit \texttt{/IMPL/QSTAT} solver instead of
the tutorial's explicit-with-mass-scaling path.

---

## 7. Reference solution

### 7.1 Closed-form (elastic regime)

Before yield, the gauge section sees a uniaxial Cauchy stress state
$\sigma_{xx} = P / A_{\text{gauge}}$, $\sigma_{yy} = \sigma_{zz} = 0$,
shear stresses zero, with strain
$\varepsilon_{xx} = \sigma_{xx} / E$ and lateral strains
$\varepsilon_{yy} = \varepsilon_{zz} = -\nu \varepsilon_{xx}$. The
load-displacement relation in the elastic regime is
$P(\delta) = E A_{\text{gauge}} \delta / L_{\text{eff}}$, with the
caveat that $L_{\text{eff}}$ depends on the compliance of the fillet
and grip-section material outside the gauge. For a perfect specimen
clamped at $x = \pm L_T/2$ with rigid grips, classical Saint-Venant
analysis gives that the gauge-section stress is uniform within $0.5\%$
of $P/A_{\text{gauge}}$ at distances greater than approximately one
gauge width $W = 12.5$ mm from the fillet \cite{ASTM_E8}, well within
the gauge length of $L_0 = 50$ mm.

The yield onset occurs at
$\sigma_{xx}^{\text{yield}} = \sigma_y = 250$ MPa, equivalently at
$\varepsilon_{xx}^{\text{yield}} = \sigma_y / E = 1.25 \times 10^{-3}$,
equivalently at $P^{\text{yield}} = \sigma_y A_{\text{gauge}} = 250
\times 75 = 18\,750$ N $= 18.75$ kN.

Apparent gauge-section axial stiffness from the integrated
stress-displacement curve is
$K_{\text{app}} = \mathrm{d}P / \mathrm{d}\delta = E A_{\text{gauge}} /
L_{\text{eff}}$. For displacement control with $\delta$ measured grip-
to-grip and load $P$ measured as the reaction on the fixed grip face,
this gives an effective length $L_{\text{eff}}$ rather than the gauge
length $L_0$; the proper apparent-modulus extraction therefore uses
the gauge-section nodal-extensometer signal $\delta_{L_0}$ between
$x = -L_0/2$ and $x = +L_0/2$, not the grip-to-grip $\delta$.

### 7.2 Experimental $\sigma$-$\varepsilon$ curve

Standard A36 stress-strain curves are tabulated in Boresi and Schmidt,
*Advanced Mechanics of Materials*, in MIL-HDBK-5J for hot-rolled
structural steel, and in the AISI Steel Construction Manual. A
representative engineering-stress-engineering-strain curve for A36 has
the following anchor points: $\sigma_y \approx 250$ MPa at
$\varepsilon \approx 0.0012$, a Lueders-band plateau from $\varepsilon
\approx 0.0012$ to $\varepsilon \approx 0.015$ at approximately
constant $\sigma$, then strain-hardening up to $\sigma_u \approx 400$
MPa at $\varepsilon \approx 0.20$, then necking down to fracture at
$\varepsilon \approx 0.30$. Stage 3 only exercises the curve up to
$\varepsilon \approx 0.005$, just past the elastic limit and before
the Lueders plateau is fully resolved, so the LAW2 fit's smoothing of
the Lueders plateau is acceptable here.

The reference solution for the validation check is constructed
analytically from the input material parameters; the experimental
curve is invoked only as a sanity-check overlay in the diagnostic
plots produced by the runner.

---

## 8. Validation success criterion

The run passes if and only if all four numerical checks below evaluate
true on the medium mesh, and all four converge with mesh refinement
toward zero residual at the expected $O(h^2)$ rate.

| Check | Quantity | Reference | Tolerance |
|---|---|---|---|
| C1 (uniformity) | std.~dev.~of $\sigma_{xx}$ across gauge bricks at $P = 0.8 P^{\text{yield}}$ | $0$ (perfectly uniform) | $\le 1\%$ of mean |
| C2 (stress vs $P/A$) | gauge-mean $\sigma_{xx}$ at $P = 0.8 P^{\text{yield}}$ | $P / A_{\text{gauge}}$ | $\le 1\%$ relative error |
| C3 (strain vs $\delta/L$) | gauge-mean $\varepsilon_{xx}$ at the same load | $\delta_{L_0} / L_0$ | $\le 1\%$ relative error |
| C4 (yield onset) | $\sigma_{xx}$ at the strain where $K_{\text{app}}$ first drops below $0.99 E$ | $\sigma_y = 250$ MPa | $\in [0.99, 1.01] \sigma_y$ |

Optional secondary checks (informational only, not pass / fail).

| Check | Quantity | Reference | Tolerance |
|---|---|---|---|
| S1 | apparent $E$ from gauge $\sigma$ vs $\varepsilon$ slope below yield | $E = 200$ GPa | $\le 1\%$ relative error |
| S2 | gauge-mean $\varepsilon_{yy}/\varepsilon_{xx}$ below yield | $-\nu = -0.30$ | $\le 1\%$ relative error |
| S3 | post-yield gauge-mean $\sigma$ at $\bar{\varepsilon}^p = 0.001$ | $a + b (0.001)^n$ from LAW2 | $\le 2\%$ relative error |

C1 demonstrates Saint-Venant decay past the fillets. C2 and C3 are the
master-plan-brief pass criteria. C4 is the master-plan-brief yield-
onset check. S1 reproduces the master plan §3 row-3 "ε within 1%"
detailed verification. S2 cross-checks Poisson contraction. S3 confirms
LAW2 is integrated correctly into the post-yield branch.

---

## 9. Toolchain runner

The runner is \texttt{runner.py} in this directory. It executes the
following pipeline.

1. Generate the GMSH mesh using the in-plane CAD construction described
   in §3.3, extruded through-thickness to HEXA8 elements. Element sets
   are written for fixed grip, loaded grip, gauge nodes, and gauge
   bricks. The output is \texttt{mesh.inp} (Abaqus format).
2. Convert \texttt{mesh.inp} to OpenRadioss \texttt{.rad} via
   \texttt{inp2rad} \cite{OpenRadiossDiscussion3222INP2RAD}.
3. Render the starter deck \texttt{stage03\_dogbone\_E8\_0000.rad} from
   a Jinja2-style template using the dimension table from §2 and the
   material card from §5. The mesh nodes and elements from step 2 are
   embedded by direct text inclusion of the converted node and brick
   blocks.
4. Render the engine deck \texttt{stage03\_dogbone\_E8\_0001.rad} from
   a Jinja2-style template with the implicit-static keyword cards from
   §6.2.
5. Invoke OpenRadioss starter and engine inside Lima
   \cite{Lima2026, Apptainer2026, OpenRadiossDiscussion2125, OpenRadiossInstall2026}
   via \texttt{limactl shell} subprocess calls. The Lima mount points
   the macOS working directory at \texttt{/mnt/host} inside the VM, so
   file I/O is in-place.
6. After the engine completes, parse the \texttt{T01} time-history file
   using the \texttt{Vortex-Radioss} reader \cite{VortexRadioss2024} to
   extract per-time-step nodal reaction force on the fixed grip,
   gauge-section nodal displacement, and gauge-section per-element
   $\sigma_{xx}$ and $\bar{\varepsilon}^p$.
7. Compute the four pass-fail checks C1-C4 of §8 and the three
   secondary checks S1-S3. Emit a \texttt{results.json} dictionary
   with the booleans, the residuals, the tolerance, and a pass / fail
   summary.
8. Export a CSV of the load-displacement curve and a CSV of the
   gauge-section $\sigma$-$\varepsilon$ curve.
9. (Optional, if \texttt{--plot} flag is set.) Generate Typst + CeTZ
   source for the load-displacement and stress-strain plots, compile
   to PDF and PNG using the Typst CLI per the user's global preferences.
   Brand colors are the user's standard palette (Garnet, Atlantic, et
   al.), no rounded edges, high contrast.

The runner is parameterized through a small Python dataclass at the
top of the file so that the material law (LAW2 vs LAW36), the mesh
density, the imposed displacement, and the Lima VM name can all be
overridden by command-line argument or environment variable without
editing the templates. Default values reproduce the configuration
documented in this spec.

---

## 10. Risks and unknowns

1. **Implicit-build availability inside Lima.** \texttt{/IMPL/QSTAT}
   needs the MUMPS-linked OpenRadioss build. Per
   \cite{OpenRadiossDiscussion2117} this is harder to build than the
   explicit-only path. Mitigation: if MUMPS is not available, fall
   back to explicit quasi-static with mass scaling per the Confluence
   tensile-test tutorial \cite{OpenRadiossConfluenceTensile2024}; the
   pass criterion is unchanged but the run takes longer.

2. **inp2rad beta status.** Element-set membership and node-set
   membership for the gauge / fillet / grip subdivisions could be
   dropped or renumbered by the converter
   \cite{OpenRadiossDiscussion3222INP2RAD}. Mitigation: after
   conversion, validate that \texttt{/GRNOD} and \texttt{/GRBRIC}
   blocks have the expected node and element counts, and abort with
   a diagnostic if they do not.

3. **Lueders plateau is not in LAW2.** A36 has a real Lueders plateau
   (load drop and reload at constant stress over $\varepsilon \in
   [0.0012, 0.015]$). LAW2's smooth power-law hardening misses this.
   Mitigation: stage 3's pass criterion ends at $\varepsilon =
   0.005$, before the Lueders plateau is over, so the smoothing only
   matters for the secondary check S3. If S3 fails, switch to LAW36
   and supply a tabulated curve that includes the plateau.

4. **Effective-length $L_{\text{eff}}$ is mesh-dependent.** The
   compliance of the fillet and partial grip is mesh-resolution-
   dependent at the percent level. Mitigation: take all stress and
   strain measurements from the gauge-section interior using the
   gauge-section nodal extensometer (gauge-endpoint nodes at
   $x = \pm L_0/2$), not from the grip-to-grip displacement. This
   makes C3 mesh-independent to the same order as the elastic kernel.

5. **Reduced vs full integration on HEXA8.** Reduced one-point
   integration on \texttt{/PROP/TYPE14} is faster but admits hourglass
   modes; full $2 \times 2 \times 2$ integration is the default and is
   used here. Mitigation: the master plan calls for the default; if a
   future study moves to reduced integration, the hourglass control
   parameters \texttt{Ihq} and \texttt{qh} must be set non-zero.

6. **Compute time on Apple Silicon via Lima.** All work crosses the
   Lima VM boundary, which adds I/O overhead. For the medium mesh
   ($\sim 13\,000$ HEXA8) the implicit run is expected to take
   minutes, well within practical bounds; for the fine mesh
   ($\sim 26\,000$) it may take tens of minutes. Mitigation: run
   coarse first, use it to commission the toolchain, then queue the
   medium and fine meshes overnight.

7. **Definition of "yield onset" depends on probe granularity.**
   Detecting when the apparent stiffness first drops below $0.99 E$
   requires a fine enough output cadence on the load-displacement
   curve. Mitigation: \texttt{/ANIM/DT} and \texttt{/TH} cadence are
   set so that there are at least $50$ time-history samples in the
   elastic regime, providing a smooth numerical derivative.

8. **GMSH structured-quad-extrusion fragility.** The fillet region's
   structured-quad recombine may fail for certain edge mesh densities
   on the fillet arc and produce mixed quad / triangle regions which
   then yield mixed HEXA / wedge in the extrusion. Mitigation: GMSH
   options pin the fillet arc to a fixed number of subdivisions
   ($8$ on the coarse, $12$ on medium, $16$ on fine) and the script
   asserts on output that no triangles or wedges are present.

9. **Possible documented gap on per-element gauge-stress output.**
   \texttt{/TH/BRICK} writes integration-point or element-averaged
   stress depending on the option flag. Mitigation: the runner
   queries the resulting time history for the expected output
   variable name and field shape and aborts with a diagnostic if the
   mapping is wrong.
