# Stage 07 - ASTM D3039 Unidirectional Tow Tension, On-Axis and Off-Axis

Author: J.C. Vaught
Date: 2026-04-29
Standard: ASTM D3039/D3039M-17
Material card: IM7/8552 from Soden, Hinton, Kaddour 1998
Solver: OpenRadioss (implicit static, MUMPS-linked build)
Element family: solid HEXA8 only

This stage is the cleanest verification of solid orthotropic material in
OpenRadioss. It isolates one capability, the assembly and rotation of the
3D orthotropic stiffness tensor on a /PROP/TYPE14 brick driven by /MAT/LAW25,
and proves that OpenRadioss recovers the exact closed-form on-axis modulus
$E_1$ and the closed-form 45 degree off-axis modulus $E_x$. If this stage
fails, every downstream composite stage (8 cross-ply, 9 quasi-iso, 11 AP-PLY
mesoscale, 12 DCB/ENF, 13 LVI, 14 CAI, 15-16 ballistic) is unreliable
because the orthotropic stiffness assembly is shared by all of them.

---

## 1. Goal and Scope

Two FEM runs on geometrically identical-material-class coupons under
displacement-controlled uniaxial tension.

Run 7A. On-axis 0 degree UD tension. Fibers aligned with the loading axis.
Recover apparent longitudinal modulus $E_1^{\text{FEM}}$ from the linear
slope of $\sigma_{xx}$ versus $\varepsilon_{xx}$ in the gauge section.
Compare to the LAW25 input value $E_1 = 171$ GPa. Pass if
$|E_1^{\text{FEM}} - E_1| / E_1 \le 0.02$.

Run 7B. Off-axis 45 degree UD tension. Fibers rotated 45 degrees in plane
about the through-thickness axis $z$ relative to the loading axis $x$.
Recover apparent uniaxial modulus $E_x^{\text{FEM}}$. Compare to the
analytic closed-form transformation for plane-stress orthotropic ply
\autocite{Jones1999MechanicsCompositeMaterials, DanielIshai2006}

$$
\frac{1}{E_x(\theta)} \;=\; \frac{c^4}{E_1} \;+\; \frac{s^4}{E_2}
\;+\; \left(\frac{1}{G_{12}} - \frac{2\nu_{12}}{E_1}\right) c^2 s^2,
\qquad c = \cos\theta, \; s = \sin\theta.
$$

At $\theta = 45^{\circ}$, $c^2 = s^2 = 1/2$ and $c^4 = s^4 = 1/4$, which
collapses to

$$
E_x\bigl(45^{\circ}\bigr) \;=\; \frac{4}{\dfrac{1}{E_1} + \dfrac{1}{E_2}
    + \dfrac{1}{G_{12}} - \dfrac{2\nu_{12}}{E_1}}.
$$

For the IM7/8552 card below this evaluates to $E_x(45^{\circ}) \approx
13.27$ GPa (numerical value computed by the runner, see section 9).
Pass if $|E_x^{\text{FEM}} - E_x^{\text{analytic}}| / E_x^{\text{analytic}}
\le 0.02$.

The two runs share mesh topology, material card, property card, and
solver configuration. Only the per-element material orientation differs.
That isolation is the point of this stage.

Verification only, not validation. No comparison to physical specimens
in this stage. Soden-Hinton-Kaddour 1998 supplies the material card
\autocite{SodenHintonKaddour1998}.

---

## 2. Material Card (IM7/8552 per Soden, Hinton, Kaddour 1998)

\autocite{SodenHintonKaddour1998}. SI base, units in the deck are mm,
ms, kg, N, MPa (the OpenRadioss "engineering" unit set; mass density
is in $10^{-3}\,\mathrm{g/mm^3}$ which equals $10^3\,\mathrm{kg/m^3}$;
LAW25 stiffness inputs are in MPa).

| Symbol | Quantity | Value | Units (deck) | Units (SI) |
|---|---|---|---|---|
| $\rho$ | density | $1.58 \times 10^{-3}$ | $\mathrm{g/mm^3}$ | $1580\,\mathrm{kg/m^3}$ |
| $E_1$ | longitudinal modulus | $171000$ | MPa | $171\,\mathrm{GPa}$ |
| $E_2 = E_3$ | transverse modulus | $9080$ | MPa | $9.08\,\mathrm{GPa}$ |
| $\nu_{12} = \nu_{13}$ | major Poisson ratio | $0.32$ | - | - |
| $\nu_{23}$ | through-thickness Poisson ratio | $0.50$ | - | - |
| $G_{12} = G_{13}$ | in-plane shear modulus | $5290$ | MPa | $5.29\,\mathrm{GPa}$ |
| $G_{23}$ | through-thickness shear modulus | $3025$ | MPa | $3.025\,\mathrm{GPa}$ |
| $X_T$ | longitudinal tensile strength | $2326$ | MPa | $2.326\,\mathrm{GPa}$ |
| $Y_T$ | transverse tensile strength | $64.7$ | MPa | $64.7\,\mathrm{MPa}$ |
| $S$ | in-plane shear strength | $92.3$ | MPa | $92.3\,\mathrm{MPa}$ |

The transverse $\nu_{23} = 0.50$ and $G_{23}$ values are not given
explicitly by Soden 1998 for IM7/8552 in the WWFE-I tabulation; we use
$\nu_{23} = 0.50$ (the WWFE-II convention for IM7/8552 in Kaddour and
Hinton 2013) and derive $G_{23} = E_2 / [2(1+\nu_{23})] = 9080/3 \approx
3025\,\mathrm{MPa}$ \autocite{KaddourHinton2013WWFEII}. These two
constants do not appear in the on-axis or 45 degree closed forms used
for the pass criterion; they affect only the through-thickness stress
state, which is uniformly small in a 1 mm thin coupon under in-plane
tension. Their precise values are not load bearing for this stage.

Strength inputs ($X_T$, $Y_T$, $S$) are present in the LAW25 deck but
are not exercised in this stage because the applied displacement is
kept in the linear regime (see section 5). They are included so the
deck is reusable for stage 6 (failure criteria) and stage 9 (laminate
first-ply-failure) without re-templating the material card.

LAW25 is the documented OpenRadioss orthotropic composite card on
solid elements \autocite{OpenRadioss_composite_intro}. Its solid-element
support is confirmed by the OpenRadioss audit row 7 and row G
(/PROP/TYPE6, /PROP/TYPE14, /PROP/TYPE20, /PROP/TYPE21, /PROP/TYPE22).
We use /PROP/TYPE14 (general orthotropic solid) per master plan
section 1 and audit row 9.

---

## 3. Geometry

ASTM D3039/D3039M-17 standard balanced symmetric coupons
\autocite{ASTM_D3039}.

Coupon 7A. On-axis 0 degree UD.
- Total length $L = 250\,\mathrm{mm}$.
- Width $w = 15\,\mathrm{mm}$.
- Thickness $t = 1\,\mathrm{mm}$ (single ply).
- Coordinate system: $x$ along the length (loading axis), $y$ along the
  width, $z$ through the thickness.

Coupon 7B. Off-axis 45 degree UD.
- Total length $L = 250\,\mathrm{mm}$.
- Width $w = 25\,\mathrm{mm}$ (wider per ASTM D3039 off-axis
  recommendation to mitigate the Pagano-Halpin oblique-end
  constraint, per \autocite{PaganoHalpin1968}; using a wider coupon
  reduces the parasitic shear-extension coupling at the grips).
- Thickness $t = 1\,\mathrm{mm}$ (single ply).
- Coordinate system identical to 7A.

Tabs and grips. ASTM D3039 specifies bonded composite or aluminum
tabs at each grip end. We deliberately replace the tab assembly with
a friction-only clamped-face simplification at each end of the coupon
(Dirichlet BCs on the end face nodes; see section 5). This is a known
simplification that biases the apparent modulus by introducing a small
Saint-Venant zone at each grip that is excluded from the gauge
extraction window. The bias is documented for off-axis tension in
\autocite{PaganoHalpin1968}. We mitigate it by extracting modulus from
a central gauge window of length $L_g = 100\,\mathrm{mm}$ (40 percent
of the total length, centered at $x = L/2$), which is well outside
the Saint-Venant decay zone of approximately one coupon width on
either side of each grip face for the on-axis case. For the 45 degree
case the decay length is longer (a few times $w$); the gauge window
is sized accordingly and the runner reports both a wide-window and
narrow-window extracted modulus to make the residual end-effect
explicit.

This simplification is acknowledged here so the spec is honest. A
follow-on stage with explicit tabs is straightforward but outside the
scope of this verification.

---

## 4. Mesh

Solid HEXA8 only, single ply through thickness. The mesh is structured
(swept extruded prism) so that nodes line up across the entire coupon
length, simplifying the gauge-window extraction.

Coupon 7A (on-axis, $L \times w \times t = 250 \times 15 \times 1$ mm).
- Length: 80 elements ($\Delta x = 3.125$ mm).
- Width: 8 elements ($\Delta y = 1.875$ mm).
- Thickness: 1 element ($\Delta z = 1.0$ mm).
- Element count: $80 \times 8 \times 1 = 640$ HEXA8.
- Aspect ratio: longest edge $\Delta x = 3.125$ mm, shortest $\Delta z = 1$ mm,
  ratio 3.125. Acceptable for linear elastic verification.

Coupon 7B (off-axis 45 degree, $L \times w \times t = 250 \times 25 \times 1$ mm).
- Length: 80 elements ($\Delta x = 3.125$ mm).
- Width: 12 elements ($\Delta y \approx 2.083$ mm).
- Thickness: 1 element ($\Delta z = 1.0$ mm).
- Element count: $80 \times 12 \times 1 = 960$ HEXA8.

Both meet the stage brief's minima ($\ge 6$ elements across width,
$\ge 40$ along length, single ply through thickness). The 7B width
count is increased to 12 to better resolve the through-width shear
gradient that off-axis tension introduces.

Mesh objectivity is not in scope for stage 7 (linear-elastic, no
softening). A single mesh per run is sufficient. The runner does
optionally rerun with $h/2$ width refinement to confirm
sub-1-percent convergence in $E_1$ and $E_x$ if requested.

Mesh generation. Pure structured grid, generated in Python by the
runner via direct node and element list construction, written
directly into the OpenRadioss starter deck (no GMSH dependency for
this stage; see runner section 7). This keeps the verification path
minimal and self-contained.

---

## 5. Boundary Conditions and Loading

Identical for 7A and 7B. The only difference between runs is the
material orientation (section 6).

Grip-A face (left end, $x = 0$): clamped. All nodes on this face have
$u_x = u_y = u_z = 0$ via /BCS on the grip-A node group.

Grip-B face (right end, $x = L = 250$ mm): displacement controlled.
All nodes on this face are constrained to $u_y = u_z = 0$ (no lateral
or through-thickness translation), and $u_x = u_x^{\max}$ is ramped
linearly from 0 to $u_x^{\max}$ over the implicit pseudo-time interval
$[0, 1]$.

Load magnitude. Target a maximum longitudinal strain
$\varepsilon_{xx}^{\max} = 2 \times 10^{-3}$ (0.2 percent). For
$L = 250$ mm this gives $u_x^{\max} = 0.5$ mm. For the on-axis run the
maximum stress is $E_1 \cdot \varepsilon = 171000 \cdot 0.002 = 342$ MPa,
which is 14.7 percent of $X_T = 2326$ MPa, deeply inside the linear
regime. For the off-axis run the apparent stress is
$E_x \cdot \varepsilon \approx 13270 \cdot 0.002 \approx 26.5$ MPa, with
the largest in-fiber-frame component being shear $\tau_{12}$ on the
order of $\sigma_{xx} \sin\theta \cos\theta = 13.3$ MPa (29 percent of
$S = 92.3$ MPa). Linear regime.

Solver. /IMPL static (linear elastic, geometric linear). We are in
the small-strain small-rotation regime; no need for /IMPL/NONLIN.
A single load step with sub-stepping disabled is sufficient; the
runner uses 5 implicit increments to give 5 sample points along
the linear $\sigma$-$\varepsilon$ curve so that the slope can be
extracted from a least-squares fit rather than a single secant.

Why displacement control rather than force control. Stage 7 is a
modulus-recovery test, not a strength test. Displacement control
gives a clean $\sigma_{xx} = F_{Bx} / (w t)$ readout where $F_{Bx}$
is the resultant reaction on grip-B (recovered from /TH/NODE on the
grip-B node set, summed). Force control would require a separate
applied-load card and would couple modulus extraction with
applied-traction-distribution choices.

---

## 6. Per-Element Material Orientation

This is the central pedagogy of stage 7. There are two documented
mechanisms in OpenRadioss for setting the material orientation on a
solid element governed by /PROP/TYPE14 plus /MAT/LAW25
\autocite{OpenRadioss_composite_intro, OpenRadioss_audit_row_8}.

Mechanism A. Per-element angle $\psi$ on /PROP/TYPE14.
The /PROP/TYPE14 card carries an in-plane orientation angle $\psi$
(degrees) measured between the element's first material direction
and a reference direction. With the orthotropy flag $\mathrm{Iorth} = 1$
(co-rotational orthotropic frame, the default for composite plies),
the material direction 1 is defined as the projection of the
reference vector $\vec{V}$ onto the element face, rotated by $\psi$
in the plane of the face. The reference vector $\vec{V}$ is supplied
on the property card; for a brick whose top face lies in the global
$xy$ plane, $\vec{V} = (1, 0, 0)$ aligns the unrotated material 1
direction with global $x$, and $\psi$ rotates it about the global $z$
axis (the through-thickness normal).

Mechanism B. /SKEW reference frame on /PROP/TYPE14.
A /SKEW card defines an orthonormal triad $(X', Y', Z')$ in the global
frame. /PROP/TYPE14 references this skew via the $\mathrm{Iskew\_id}$
field, and the material 1, 2, 3 axes are taken to align with $X', Y',
Z'$. To rotate the orthotropic frame about the global $z$ axis by
$\theta$, the skew is defined with $X' = (\cos\theta, \sin\theta, 0)$,
$Y' = (-\sin\theta, \cos\theta, 0)$, $Z' = (0, 0, 1)$. /SKEW is
documented for stage 8 verification work (audit row 8) as the
standard mechanism for per-element or per-group orientation.

Decision for stage 7. We use **Mechanism A** (per-element angle
$\psi$ on /PROP/TYPE14 with reference vector $\vec{V} = (1, 0, 0)$),
because it is the simpler and more localized form. The reference
vector is set on the property card; the per-property angle $\psi$
is set per run.

- Run 7A (on-axis 0 degree): $\psi = 0^{\circ}$. Material 1 aligns
  with global $x$ (loading axis). Fibers along loading axis.
- Run 7B (off-axis 45 degree): $\psi = 45^{\circ}$. Material 1 is
  rotated 45 degrees about global $z$ from the global $x$ axis.
  Fibers at 45 degrees to the loading axis.

If a future regression test wants per-element variation (for instance
a fiber-angle gradient), Mechanism B with one /SKEW per element group
is the documented path. We note this for stage 8, which verifies the
orientation transformation algebraically across $\theta \in
\{0, 15, 30, 45, 60, 75, 90\}^{\circ}$.

---

## 7. OpenRadioss Deck Skeleton

The deck is templated by the runner. The skeleton below uses Jinja2
placeholders in the form `{{ name }}`. Comments use `#` per the
OpenRadioss convention.

Starter deck `stage7_<run>_0000.rad`.

```
#RADIOSS STARTER
/BEGIN
stage7_{{ run }}_run
      2024         0
                  Mg                  mm                   s
                  Mg                  mm                   s
/UNIT/1
mass length time
Mg mm s
#---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|
/MAT/LAW25/1
IM7_8552_Soden1998
# rho_i        E11        E22        E33    nu12       nu13   nu23
{{ rho }}    {{ E1 }}   {{ E2 }}   {{ E2 }}  {{ nu12 }}  {{ nu12 }}  {{ nu23 }}
# G12         G23        G31        Iform   Ioff
{{ G12 }}   {{ G23 }}   {{ G12 }}    1       0
# Xt          Xc          Yt          Yc          S
{{ Xt }}     1700.0      {{ Yt }}    200.0       {{ S }}
# (additional LAW25 lines for damage and softening left at defaults)
0.0  0.0  0.0  0.0  0.0
0.0  0.0  0.0  0.0
#
/PROP/TYPE14/1
solid_UD_psi_{{ psi }}
# Iint  Iframe  Ihbe   Iplas  Itetra4  Itetra10
   1      1       0      0       0         0
# qa  qb  h  Vx  Vy  Vz  skew_id  Iorth  Icpre  Inpts  dn  Iss_rot  ...
0.0  0.0  0.0  1.0  0.0  0.0  0  1  0  0  0  0
# Phi (per-element / per-property orientation angle, degrees)
{{ psi }}
#
/SKEW/FIX/1
global
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
#
# nodes ----------------------------------------------------------
/NODE
{{ nodes_block }}
#
# elements -------------------------------------------------------
/BRICK/1
{{ bricks_block }}
#
# part / subset --------------------------------------------------
/PART/1
ud_coupon
1 1 0
#
# node groups ----------------------------------------------------
/GRNOD/NODE/1
gripA_face
{{ gripA_nodes }}
/GRNOD/NODE/2
gripB_face
{{ gripB_nodes }}
/GRNOD/NODE/3
gauge_window
{{ gauge_nodes }}
#
# kinematic BCs --------------------------------------------------
/BCS/1
gripA_clamped
111 111 1   # Tx Ty Tz / Rx Ry Rz / skew=1 (global)
#
/BCS/2
gripB_lateral
011 111 1
#
# imposed displacement on grip-B in x ----------------------------
/IMPDISP/1
gripB_pull
2 1 1
# fct_id  dir   skew_id  Tstart   Tstop   scale
   1      1       1       0.0     1.0     {{ u_max }}
#
/FUNCT/1
ramp
0.0 0.0
1.0 1.0
#
# implicit static -------------------------------------------------
/IMPL/PRINT/N
1
/IMPL/SOLVER/1
0  0  3  0  0
/IMPL/NONLIN/SMDISP
{{ n_steps }}  0.0  1.0  1.0e-3  1.0e-6  20  0
/IMPL/DT/STOP
1.0e-12   1.0
/IMPL/DTINI
{{ dt_ini }}
#
# requested time-history outputs ---------------------------------
/TH/NODE/1
gripB_reaction
2 0 0
DX FX
/TH/NODE/2
gauge_window_th
3 0 0
DX
#
/ANIM/BRICK/TENS/STRESS
/ANIM/BRICK/TENS/STRAIN
/ANIM/DT
0.2
#
/END
```

Engine deck `stage7_<run>_0001.rad`.

```
#RADIOSS ENGINE
/RUN/stage7_{{ run }}_run/1
1.0
#
/PRINT/-100
/HIS/DT
0.05
/STOP
0.0  1.0e-3  0.0  0.0
#
/END
```

Notes on the deck.

- The /MAT/LAW25 line ordering uses the audit-confirmed solid form:
  $\rho$, $E_{11}$, $E_{22}$, $E_{33}$, $\nu_{12}$, $\nu_{13}$, $\nu_{23}$,
  then $G_{12}$, $G_{23}$, $G_{31}$, $\mathrm{Iform}$, then strengths.
  Exact field positions and trailing-line damage/softening parameters
  vary slightly by OpenRadioss version; the runner consults the
  installed OpenRadioss reference to confirm before launching, and
  fails fast with a deck-validation error if the starter rejects the
  card.
- /PROP/TYPE14 uses $\mathrm{Iorth} = 1$, $\vec{V} = (1, 0, 0)$, and
  per-property $\Phi$ rotation angle (Mechanism A). For 7A,
  $\Phi = 0$. For 7B, $\Phi = 45$.
- /SKEW is defined as the global skew (identity) for use as the
  reference frame in /BCS and /IMPDISP. The optional Mechanism B
  alternative would define a second /SKEW with the rotation already
  baked in and would set $\mathrm{Iskew\_id}$ on /PROP/TYPE14
  instead of using $\Phi$.
- /IMPL/NONLIN/SMDISP requests a small-displacement implicit solve
  with $n$ load steps. For linear elastic this could equally be
  /IMPL/LINEAR; we use NONLIN with the geometric-nonlinear flag off
  so the runner can rerun the same coupon at large strain in a future
  stage without restructuring the deck.
- The OpenRadioss build must be MUMPS-linked. If only the explicit
  build is available the runner falls back to /DYREL dynamic
  relaxation (a documented quasi-static approximation, audit row 1).

---

## 8. Reference Solution

Closed forms for the linear-elastic homogeneous orthotropic UD ply
under uniaxial in-plane tension, from \autocite{Jones1999MechanicsCompositeMaterials}
chapter 2 and \autocite{DanielIshai2006} chapter 3.

On-axis (Run 7A). For $\theta = 0$ the loading axis coincides with
material direction 1, so

$$
E_x(0) \;=\; E_1 \;=\; 171\,000\,\mathrm{MPa} \;=\; 171\,\mathrm{GPa}.
$$

The transverse Poisson contraction is $\nu_{12} = 0.32$, so the
expected lateral strain is $\varepsilon_{yy} = -\nu_{12}\,\varepsilon_{xx}
= -6.4 \times 10^{-4}$ at the maximum applied $\varepsilon_{xx} = 2
\times 10^{-3}$.

Off-axis (Run 7B). For $\theta = 45^{\circ}$,

$$
E_x(45^{\circ}) \;=\; \frac{4}{\dfrac{1}{E_1}+\dfrac{1}{E_2}+\dfrac{1}{G_{12}}-\dfrac{2\nu_{12}}{E_1}}.
$$

Substituting $E_1 = 171\,000$, $E_2 = 9080$, $G_{12} = 5290$,
$\nu_{12} = 0.32$ (all MPa or dimensionless):

$$
\frac{1}{E_1} = 5.848 \times 10^{-6},\quad
\frac{1}{E_2} = 1.1013 \times 10^{-4},\quad
\frac{1}{G_{12}} = 1.8904 \times 10^{-4},\quad
\frac{2\nu_{12}}{E_1} = 3.7427 \times 10^{-6}.
$$

Sum: $5.848 \times 10^{-6} + 1.1013 \times 10^{-4} + 1.8904 \times
10^{-4} - 3.7427 \times 10^{-6} = 3.0153 \times 10^{-4}\,\mathrm{MPa}^{-1}$.

$$
E_x(45^{\circ}) \;=\; \frac{4}{3.0153 \times 10^{-4}\,\mathrm{MPa}^{-1}}
\;\approx\; 13\,265\,\mathrm{MPa} \;\approx\; 13.27\,\mathrm{GPa}.
$$

The runner re-evaluates this number at runtime from the material
constants in the LAW25 dictionary, so the closed form lives in code,
not in this prose. The number 13.27 GPa is shown here as a sanity
anchor only.

---

## 9. Pass Criterion (Per Master Plan Stage 7)

Both criteria must pass for stage 7 to be considered green.

C1 (on-axis). $\dfrac{|E_1^{\text{FEM}} - E_1|}{E_1} \le 0.02$.

C2 (off-axis). $\dfrac{|E_x^{\text{FEM}} - E_x^{\text{analytic}}(45^{\circ})|}{E_x^{\text{analytic}}(45^{\circ})} \le 0.02$.

Modulus extraction protocol.

1. From /TH/NODE on grip-B, recover the time history of the summed
   reaction $F_{Bx}(t)$.
2. Compute apparent stress $\sigma_{xx}(t) = F_{Bx}(t) / A$ with
   $A = w \cdot t$ (15 mm² for 7A, 25 mm² for 7B).
3. From /TH/NODE on the gauge window, recover the time history of
   the displacement $u_x$ at two cross-sections at $x_1 = (L-L_g)/2$
   and $x_2 = (L+L_g)/2$ (the gauge-window endpoints).
4. Compute apparent strain $\varepsilon_{xx}(t) =
   (u_x(x_2,t) - u_x(x_1,t)) / L_g$.
5. Linear-regress $\sigma_{xx}$ on $\varepsilon_{xx}$ over the
   five implicit increments. The slope is $E^{\text{FEM}}$.
6. Report the regression $R^2$. Require $R^2 \ge 0.9999$ as a
   sanity check that the response is in the linear regime.

The runner additionally reports the through-width-averaged
$\sigma_{xx}$ from /ANIM at midspan as an independent check on the
reaction-force-divided-by-area approach. The two should agree to
within 0.5 percent in the on-axis case; for the off-axis case they
may disagree by a couple of percent due to the parasitic
shear-extension coupling near the grips, in which case the
regression-based number from steps 1-5 is the canonical
$E_x^{\text{FEM}}$ for the pass test.

Stage failure modes (what could go wrong, what we will catch).

- LAW25 field ordering wrong for the installed OpenRadioss version:
  starter rejects the deck. Runner fails fast.
- /PROP/TYPE14 reference vector or $\mathrm{Iorth}$ flag wrong:
  apparent $E_1$ comes out as $E_2$ (or some intermediate). Pass test
  catches a $\sim 95$ percent error.
- /PROP/TYPE14 angle $\psi$ applied in the wrong plane (e.g. about
  $y$ instead of $z$): apparent $E_x$ at 45 degrees comes out as
  $E_2$ instead of $\sim 13.3$ GPa. Pass test catches this.
- Saint-Venant zone bias too large because $L_g$ is too generous:
  off-axis modulus shows 5-10 percent low. Mitigated by the dual
  wide-window/narrow-window report; if the narrow-window result
  passes but the wide-window does not, the spec is right and the
  end-effect is just bigger than expected, which is documented.

---

## 10. Run Sequence and Outputs

Both runs share the runner; see `runner.py` for the orchestration.

Sequence.

1. Build the on-axis deck (7A) by templating with $\psi = 0^{\circ}$,
   $w = 15$ mm, $L_g = 100$ mm.
2. Run starter, then engine, inside Lima (per master plan section 7).
3. Parse `T01` time-history into Python via `vortex-radioss` or the
   T01 ASCII export.
4. Extract $E_1^{\text{FEM}}$ per section 9.
5. Repeat for the 45 degree deck (7B) with $\psi = 45^{\circ}$,
   $w = 25$ mm.
6. Compute $E_x^{\text{analytic}}(45^{\circ})$ from the LAW25
   constants.
7. Print a pass/fail report with both relative errors and the
   regression $R^2$.
8. Optionally convert `.anim` to `.vtkhdf` via the Kitware tool
   and dump a midspan $\sigma_{xx}$ contour PNG via PyVista
   headless for documentation.
9. Write a CSV with columns
   `run, theta_deg, E_FEM_MPa, E_analytic_MPa, rel_err, R2, pass`
   to `tests/stage_07_UD_tow_D3039/results/stage7_summary.csv`
   for downstream Typst/CeTZ figure compilation.

Outputs.

- `runs/7A/stage7_7A_0000.rad`, `..._0001.rad`, `*.anim`, `*.T01`.
- `runs/7B/...`.
- `results/stage7_summary.csv` (canonical pass/fail table).
- `results/stage7_log.txt` (free-form console log).
- (Optional) `results/stage7_7A_midspan_sxx.png`,
  `results/stage7_7B_midspan_sxx.png`.

References.

\printbibliography[heading=none, keyword=stage7]

\bibliographystyle{plain}

% Per-stage citation keys used in this spec:
% ASTM_D3039
% SodenHintonKaddour1998
% HintonKaddourSoden2002
% KaddourHinton2013WWFEII
% Jones1999MechanicsCompositeMaterials
% DanielIshai2006
% PaganoHalpin1968
% (OpenRadioss audit / composite intro pages cited inline.)
