# Stage 08 -- Ply rotation transformation verification (off-axis stiffness sweep)

**Author.** J.C. Vaught
**Date.** 2026-04-29
**Stage row in `master_plan.md` Section 3.** Row 8.
**Verdict from end-to-end audit (`openradioss_endtoend_audit.md`, row 8).** PASS.
**Pass criterion.** Apparent uniaxial modulus $E_x(\theta)$ recovered from a uniaxial tensile FEM run on a single-ply solid coupon, sampled at $\theta \in \{0,15,30,45,60,75,90\}^\circ$, must agree with the closed-form Jones / Daniel-Ishai transformation
\begin{equation}
\frac{1}{E_x(\theta)} \;=\; \frac{c^4}{E_1} \;+\; \frac{s^4}{E_2} \;+\; \left(\frac{1}{G_{12}}-\frac{2\nu_{12}}{E_1}\right) c^2 s^2,
\label{eq:Extheta}
\end{equation}
with $c=\cos\theta$ and $s=\sin\theta$, to within **1 percent relative error at every angle**.

---

## 1. Goal and scope

Stage 7 verified that OpenRadioss reproduces the on-axis engineering constants $(E_1,E_2,\nu_{12},G_{12})$ of the Soden-Hinton-Kaddour 1998 IM7/8552 unidirectional carbon/epoxy ply when the fibers are aligned with the load. Stage 8 does the immediately next test. With the **same geometry, the same mesh, and the same boundary conditions**, only the per-element material orientation is changed. Seven runs are produced, one per fiber-to-load angle. The recovered apparent modulus $E_x(\theta)$ is compared against the analytic transformation in equation \eqref{eq:Extheta}.

This is the most direct test of OpenRadioss's per-element material orientation machinery. Stage 9 (cross-ply and quasi-isotropic laminates) and stage 11 (pseudo-woven mesoscale) both depend on the rotation pipeline working ply-by-ply. If stage 8 fails, those later stages cannot work, because their plies are simply this stage at multiple angles in a single deck.

The test is pure verification, not validation. The reference is a closed-form algebraic identity, not an experiment. Tolerance is therefore tight (1 percent).

---

## 2. Geometry

Identical to stage 7's $0^\circ$ tow coupon, in SI units throughout.

| Quantity | Value | Note |
|---|---|---|
| Length $L$ | $0.250$ m | gauge plus grip in one piece |
| Width $b$ | $0.025$ m | ASTM D3039-style rectangular coupon |
| Thickness $t$ | $0.001$ m | single ply, solid representation |
| Coordinate frame | $x$ along load, $y$ along width, $z$ through-thickness | same for all seven runs |

The geometry is a brick. There is no fillet, no grip tab, no notch. The load axis is $x$. The fiber direction in the local material frame is denoted $1$. The fiber-to-load angle $\theta$ is the rotation about $z$ that takes the global $x$ axis onto the local $1$ axis.

GMSH script `geom.geo` (deterministic, scriptable, no GUI):

```text
// stage_08 / geom.geo
SetFactory("OpenCASCADE");
L  = 0.250;
b  = 0.025;
t  = 0.001;
Box(1) = {0, -b/2, -t/2, L, b, t};

// physical groups for BCs and material assignment
Physical Volume("ply", 1) = {1};
Physical Surface("face_x0", 11) = {1};   // x = 0,    clamp / symmetric BC
Physical Surface("face_xL", 12) = {2};   // x = L,    pull face
Physical Surface("face_y0", 13) = {3};   // y = -b/2, symmetry candidate
Physical Surface("face_yb", 14) = {4};   // y = +b/2
Physical Surface("face_z0", 15) = {5};   // z = -t/2
Physical Surface("face_zt", 16) = {6};   // z = +t/2
```

The face IDs above are placeholders; the runner introspects the GMSH model to recover the actual surface tags and writes a deterministic mapping into the deck. The point is that one geometry is built once and reused for all seven angles; nothing about the mesh, the topology, or the BC surface IDs changes between runs.

---

## 3. Mesh

Same mesh as stage 7. Solid HEXA8 only, no shells, no tetrahedra. Structured grid through the brick.

| Direction | Element count | Element size |
|---|---|---|
| $x$ (length) | 100 | $0.0025$ m |
| $y$ (width) | 10 | $0.0025$ m |
| $z$ (thickness) | 2 | $0.0005$ m |

Element formulation: `Isolid = 24` (HEXA8 with full integration, hourglass-controlled), via `/PROP/TYPE14`. This matches stage 7 so that any orientation-dependent error is *not* an element-formulation artifact. The mesh has $100 \times 10 \times 2 = 2000$ HEXA8 elements and roughly $2{,}300$ nodes.

The off-axis runs ($\theta = 15,30,45,60,75^\circ$) develop in-plane shear stress $\sigma_{xy}$ in the global frame because the orthotropic compliance couples extension and shear when the material axes are not aligned with the load axis. The sample is long enough ($L/b = 10$) that Saint-Venant decay leaves a uniform stress region in the gauge centre regardless of $\theta$. The runner samples the gauge centre, not the grip vicinity, when extracting $E_x$.

---

## 4. Material card

Same Soden-Hinton-Kaddour 1998 IM7/8552 UD CFRP card as stage 7. SI units.

| Symbol | Value | Unit | Note |
|---|---|---|---|
| $\rho$ | $1{,}580$ | kg/m$^3$ | density |
| $E_1$ | $165{,}000\times 10^6$ | Pa | longitudinal (fiber) modulus |
| $E_2$ | $9{,}000\times 10^6$ | Pa | transverse modulus |
| $E_3$ | $9{,}000\times 10^6$ | Pa | through-thickness, taken equal to $E_2$ for transversely isotropic UD |
| $\nu_{12}$ | $0.34$ | --- | major Poisson |
| $\nu_{13}$ | $0.34$ | --- | $=\nu_{12}$ for transverse isotropy |
| $\nu_{23}$ | $0.5$ | --- | per Soden 1998 Table 1 |
| $G_{12}$ | $5{,}600\times 10^6$ | Pa | in-plane shear |
| $G_{13}$ | $5{,}600\times 10^6$ | Pa | $=G_{12}$ for transverse isotropy |
| $G_{23}$ | $3{,}000\times 10^6$ | Pa | per Soden 1998 |

Strength values are not needed for stage 8. Loading stays inside the linear-elastic regime by construction. Failure cards (`/FAIL/HASHIN`, `/FAIL/PUCK`) are not attached.

OpenRadioss card. `/MAT/LAW25` (COMPSH / CRASURV), with the elastic-only form of the card (no plastic flow on, no failure card). LAW25 is documented as compatible with solid properties `/PROP/TYPE6, TYPE14, TYPE20, TYPE21, TYPE22` per `AltairRadiossLAW25` and `AltairRadiossCompositeIntro`. We use `/PROP/TYPE14` (general solid) with the orthotropic flag `Iorth` enabled and the local 1-axis pointed by a `/SKEW` system.

---

## 5. Boundary conditions and loading

Identical for all seven runs. Only the orientation changes.

**Clamped end** ($x=0$). All nodes on `face_x0` constrained: $u_x=u_y=u_z=0$. (Three rotational DOF do not exist on solid bricks.)

**Pull end** ($x=L$). All nodes on `face_xL` driven by an imposed displacement $u_x = u_x^{\text{imp}}(t)$, with $u_y$ and $u_z$ left free so the Poisson lateral contraction and the off-axis shear are not over-constrained. The runner imposes $u_x^{\text{imp}}$ via `/IMPDISP` ramped linearly from $0$ to $\bar\varepsilon\,L$ over a fictitious unit pseudo-time. The strain target is $\bar\varepsilon = 0.002$ (0.2 percent), well inside the linear regime for IM7/8552 even at the most compliant orientation ($90^\circ$, $\bar\sigma_x \approx 18$ MPa).

**Side faces** (`face_y0`, `face_yb`, `face_z0`, `face_zt`). Free, no constraint. We deliberately do *not* impose a transverse displacement constraint, because that would over-constrain off-axis Poisson and induce a spurious $\sigma_{xy}$ that contaminates the apparent $E_x$ extraction.

**Solver.** `/IMPL/LINEAR` for static linear analysis (stage 1 already verified the implicit pipeline; this stage is also linear so the same activation applies). The MUMPS-linked build is required, see master plan section 10 risk 2.

---

## 6. Material orientation: the only knob that changes

For each $\theta_k \in \{0,15,30,45,60,75,90\}^\circ$ a separate deck is written that differs from the stage-7 $0^\circ$ deck only in one card block, the `/SKEW` definition and its `/PROP/TYPE14` reference.

Two equivalent OpenRadioss mechanisms are supported:

1. **`/SKEW/FIX` block referenced by `/PROP/TYPE14`.** A right-handed Cartesian skew is defined by the global coordinates of its origin plus two vectors. For an in-plane rotation about the global $z$ axis,
\begin{align}
\mathbf{e}_1^{\text{loc}} &= (\cos\theta,\;\sin\theta,\;0),\\
\mathbf{e}_2^{\text{loc}} &= (-\sin\theta,\;\cos\theta,\;0),\\
\mathbf{e}_3^{\text{loc}} &= (0,0,1).
\end{align}
The skew ID is referenced in `/PROP/TYPE14` via the orthotropic-frame flag `Iorth=1` and the `Iskew` field (per `AltairRadiossPropType14`).

2. **In-card angle.** `/PROP/TYPE14` accepts an in-plane angle field $\Phi$ (degrees) that rotates the orthotropic frame about the element's local $z$ axis once a skew is referenced. Setting the skew to global and $\Phi=\theta$ produces the same result as mechanism 1.

The runner uses mechanism 1 (`/SKEW/FIX`) because it makes the rotation explicit in the deck text and is robust against future Radioss changes to the in-card angle field. Mechanism 2 is a documented one-line alternative if a reviewer wants to confirm equivalence.

For all seven runs the geometry, mesh, BC surfaces, IMPDISP magnitude, and all LAW25 numerical parameters are bit-identical. The only diff between two decks is the three numbers $(\cos\theta,\sin\theta,-\sin\theta)$ inside the `/SKEW/FIX` block. This is the cleanest possible test of per-element material orientation.

---

## 7. Reference solution: derivation of $E_x(\theta)$

The closed form in equation \eqref{eq:Extheta} is the standard transverse-isotropic plane-stress lamina result. Both Jones (1999, *Mechanics of Composite Materials*, Chapter 2, equations 2.85--2.88) and Daniel-Ishai (2006, *Engineering Mechanics of Composite Materials*, Chapter 5) derive it in essentially the same way. The derivation is recorded in full here so the tolerance check in section 9 has an unambiguous reference.

### 7.1 Lamina compliance in material axes

For an in-plane orthotropic lamina under plane stress with material 1-axis along the fibers, the compliance relation between strains $\boldsymbol{\varepsilon}_{12} = (\varepsilon_1,\varepsilon_2,\gamma_{12})^\top$ and stresses $\boldsymbol{\sigma}_{12} = (\sigma_1,\sigma_2,\tau_{12})^\top$ is
\begin{equation}
\boldsymbol{\varepsilon}_{12} \;=\; \mathbf{S}\,\boldsymbol{\sigma}_{12},
\quad
\mathbf{S} \;=\;
\begin{pmatrix}
1/E_1 & -\nu_{12}/E_1 & 0 \\
-\nu_{12}/E_1 & 1/E_2 & 0 \\
0 & 0 & 1/G_{12}
\end{pmatrix}.
\label{eq:S12}
\end{equation}
This is symmetric because $\nu_{21}/E_2 = \nu_{12}/E_1$ (Maxwell-Betti).

### 7.2 Rotation matrices

Let $\theta$ be the angle from the global $x$ axis to the material 1-axis, measured counterclockwise about $z$, with $c\equiv\cos\theta$ and $s\equiv\sin\theta$. The strain transformation $\boldsymbol{\varepsilon}_{12} = T_\varepsilon\,\boldsymbol{\varepsilon}_{xy}$ and the stress transformation $\boldsymbol{\sigma}_{12} = T_\sigma\,\boldsymbol{\sigma}_{xy}$ use the engineering-shear conventions in Jones 1999 equations 2.81--2.84,
\begin{equation}
T_\sigma \;=\;
\begin{pmatrix}
c^2 & s^2 & 2cs \\
s^2 & c^2 & -2cs \\
-cs & cs & c^2-s^2
\end{pmatrix},
\quad
T_\varepsilon \;=\;
\begin{pmatrix}
c^2 & s^2 & cs \\
s^2 & c^2 & -cs \\
-2cs & 2cs & c^2-s^2
\end{pmatrix}.
\end{equation}
Note the factor-of-two difference between $T_\sigma$ and $T_\varepsilon$ that arises from using $\gamma_{12}$ rather than $\varepsilon_{12}$ (engineering shear strain rather than tensorial shear strain).

### 7.3 Off-axis compliance

The compliance in global axes is $\bar{\mathbf{S}} = T_\varepsilon^{-1}\,\mathbf{S}\,T_\sigma$. Algebra gives the standard result for the apparent uniaxial modulus along $x$,
\begin{equation}
\frac{1}{E_x(\theta)} \;=\; \bar{S}_{11}(\theta) \;=\; \frac{c^4}{E_1} \;+\; \frac{s^4}{E_2} \;+\; \left(\frac{1}{G_{12}}-\frac{2\nu_{12}}{E_1}\right)c^2 s^2.
\label{eq:Extheta_repeat}
\end{equation}
This is the criterion line of the spec.

### 7.4 Side checks recovered for free

The same algebra produces companion expressions that the runner reports for diagnostic value but does not gate the pass criterion on:
\begin{align}
\frac{1}{G_{xy}(\theta)} &\;=\; 2\!\left(\frac{2}{E_1}+\frac{2}{E_2}+\frac{4\nu_{12}}{E_1}-\frac{1}{G_{12}}\right) c^2 s^2 \;+\; \frac{1}{G_{12}}(c^2-s^2)^2,\\
\nu_{xy}(\theta) &\;=\; E_x(\theta)\!\left[\frac{\nu_{12}}{E_1}(c^4+s^4) - \left(\frac{1}{E_1}+\frac{1}{E_2}-\frac{1}{G_{12}}\right) c^2 s^2 \right].
\end{align}
Numerical values for the IM7/8552 card at the seven sweep angles, for sanity:

| $\theta$ (deg) | $E_x$ (GPa) | $\nu_{xy}$ | $G_{xy}$ (GPa) |
|---:|---:|---:|---:|
| 0  | 165.00 | 0.340 | 5.60 |
| 15 |  59.96 | 0.338 | 7.05 |
| 30 |  23.22 | 0.297 | 14.57 |
| 45 |  13.72 | 0.225 | 31.24 |
| 60 |  10.46 | 0.134 | 14.57 |
| 75 |   9.29 | 0.052 | 7.05 |
| 90 |   9.00 | 0.019 | 5.60 |

Numerical values produced by `runner.py --analytic-only` to four significant figures. The order-of-magnitude drop from $E_x(0^\circ)=165$ GPa to $E_x(45^\circ)=13.7$ GPa to $E_x(90^\circ)=9.0$ GPa is the single largest dynamic range of any verification stage in the progression; it is exactly why this stage is the canonical test of per-element material orientation. Note also the non-monotone behavior of $G_{xy}(\theta)$ -- it peaks at $45^\circ$ where the off-axis test is in pure shear in the rotated material frame -- which is a useful sanity gate that the rotation algebra is correct.

The runner recomputes the analytic table in code at runtime so the spec table is never relied on numerically.

---

## 8. OpenRadioss deck skeleton

Single starter deck `coupon_<theta>_0000.rad` and engine deck `coupon_<theta>_0001.rad`, templated per angle. The skeleton (cards listed in order, only the salient ones).

```text
#-- starter (coupon_<theta>_0000.rad)
/BEGIN
TITLE   stage_08_off_axis_theta_<theta>
#       run_id  version
        12345    2026
#-- units (SI)
                 kg                    m                    s
                 kg                    m                    s
#-- material LAW25 elastic-only IM7/8552 (Soden 1998)
/MAT/LAW25/1
IM7_8552_UD
#         rho_i
        1.580E+03
#       E_1            E_2           nu_12             G_12
       1.6500E+11    9.0000E+09    3.4000E-01      5.6000E+09
#       E_3           nu_13           nu_23            G_13            G_23
       9.0000E+09    3.4000E-01      5.0000E-01      5.6000E+09     3.0000E+09
# (no plastic / failure parameters: linear elastic only, all strengths set to large dummies)
#
#-- skew (in-plane rotation by theta about z)
/SKEW/FIX/9
ply_skew_theta_<theta>
#       Ox     Oy     Oz   Ax     Ay     Az    Bx     By     Bz
       0.000  0.000  0.000  <cos>  <sin>  0.000  <-sin> <cos>  0.000
#
#-- solid property: TYPE14 with orthotropic frame, references skew 9
/PROP/TYPE14/1
ply_solid_prop
#  Ihex Iframe Iorth Iplas    Phi  Iskew
       24      1     1     0   0.0       9
#
#-- part: associates property + material with the element group
/PART/1
ply_part 1 1
#
#-- mesh data is included from a separate file produced by gmsh+inp2rad
#include coupon_mesh.rad
#
#-- boundary conditions (clamp at x=0, free elsewhere)
/BCS/1
clamp_x0
                  X Y Z 0 0 0      <grnod_face_x0>
#
#-- imposed displacement at x=L: ramp to ubar = strain_bar * L
/IMPDISP/1
pull_xL_x
#  Ifunc  Dir  Iskew  Iframe  Sens
        7    1      0       0     0
#       grnod          scale (ubar)
        <grnod_face_xL>      5.0000E-04
#
#-- linear static implicit
/IMPL/LINEAR
#
/END
```

```text
#-- engine (coupon_<theta>_0001.rad)
/RUN/<runname>/1
                  1.0
/IMPL/LINEAR
/PRINT/-1
/H3D/DT
              0.0     1.0
/H3D/SOLID/STRESS/ALL
/H3D/SOLID/STRAIN/ALL
/STOP
```

The function ID `7` referenced by `/IMPDISP` is a unit-ramp tabulated function $f(t) = t$ for $t\in[0,1]$, written into the deck verbatim by the runner.

The deck is intentionally minimal. Stages 9 and 11 will reuse this skeleton with multiple `/SKEW` blocks and multiple parts.

---

## 9. Validation success criterion (gating logic)

Per master plan section 3 row 8, **engineering constants vs analytic within 1 percent**. Concretely, for every sweep angle $\theta_k$ the runner extracts a single number $E_x^{\text{FEM}}(\theta_k)$ from the linear regime and computes
\begin{equation}
r(\theta_k) \;\equiv\; \frac{|E_x^{\text{FEM}}(\theta_k) - E_x^{\text{analytic}}(\theta_k)|}{E_x^{\text{analytic}}(\theta_k)} \,.
\end{equation}
The stage **passes** if $r(\theta_k) \le 0.01$ for all $k$. Otherwise the stage **fails** and reports which angles were out of tolerance.

Extraction protocol (deterministic, no manual inspection):

1. Read the H3D engine output via the Kitware `openradioss-to-vtkhdf` converter (master plan section 5).
2. Mask elements whose centroid lies in the gauge sub-volume $|x-L/2|\le L/4$, $|y|\le b/2$, $|z|\le t/2$. (i.e. drop the grip influence zones near $x=0$ and $x=L$ inside one width.)
3. Compute the volume-averaged Cauchy stress $\langle\sigma_{xx}\rangle$ over the gauge sub-volume.
4. Compute the engineering strain $\bar\varepsilon$ as $u_x^{\text{imp}}/L$. Because the run is linear-elastic and $\bar\varepsilon=0.002$ is below the elastic limit at all angles (worst case $90^\circ$, $\bar\sigma_x\approx18$ MPa is one-sixth of $Y_T$), the slope $\langle\sigma_{xx}\rangle/\bar\varepsilon$ is the apparent modulus.
5. Report $E_x^{\text{FEM}}(\theta_k) = \langle\sigma_{xx}\rangle/\bar\varepsilon$.

Diagnostic outputs that do not gate the stage but are recorded:

- $\nu_{xy}(\theta_k)$ from $-\langle\varepsilon_{yy}\rangle/\bar\varepsilon$. Compared against the closed form in section 7.4 for sanity, with a 5 percent tolerance.
- A Typst+CeTZ figure `Ex_vs_theta.pdf` showing FEM points overlaid on the analytic curve, drawn from the CSV exported by the runner. Plot uses the brand colors per global preferences. Garnet (#73000A) for the analytic curve, Black for the FEM points.
- Volume-averaged $\sigma_{xx}$, $\sigma_{yy}$, $\tau_{xy}$, $\varepsilon_{xx}$, $\varepsilon_{yy}$, $\gamma_{xy}$ in the gauge for every angle, dumped to `summary.csv`.

The Typst+CeTZ figure is authored separately. The runner produces the CSV; the Typst document imports it and draws the figure (per global preferences on figure generation).

---

## 10. Cross-cutting notes and dependencies

**Cites.** Jones, R. M., *Mechanics of Composite Materials*, 2nd ed., 1999, Chapter 2 (transformation matrices, off-axis compliance) -- bib key `Jones1999MechanicsCompositeMaterials`. Daniel, I. M. and Ishai, O., *Engineering Mechanics of Composite Materials*, 2nd ed., 2006, Chapter 5 (worked off-axis examples with the same closed form) -- bib key `DanielIshai2006`. Lamina property card from Soden, Hinton, Kaddour 1998, *Composites Science and Technology* 58(7), 1011--1022 -- bib key `SodenHintonKaddour1998`.

**OpenRadioss documentation.** `/MAT/LAW25` -- bib key `AltairRadiossLAW25`; `/PROP/TYPE14` -- bib key `AltairRadiossPropType14`; `/IMPL/LINEAR` -- bib key `AltairRadiossImplLinear`; composite intro / failure-card compatibility -- bib key `AltairRadiossCompositeIntro`. Per master plan, MUMPS-linked OpenRadioss build is required for `/IMPL/LINEAR`; if the prebuilt binaries are explicit-only (master plan risk 2), fall back to explicit quasi-static loading (`/DYREL` plus mass scaling) -- the recovered $E_x$ is identical for a small enough loading rate, and the 1 percent criterion is still met because the regime is purely elastic.

**Hard requirements.** Solid elements only (no shells), SI units, no GUI, deck templated from Python -- all four are met by the skeleton in section 8 and the runner in `runner.py`.

**Upstream dependency.** Stage 7 has produced its $0^\circ$ deck and verified $E_1$ within 1 percent. This stage's $\theta=0^\circ$ deck is identical to that one, modulo trivially redundant `/SKEW/FIX` block defining the global frame. Recovering the same $E_x(0^\circ)=E_1$ in stage 8 is a sanity gate that the rotation infrastructure has not broken the on-axis case.

**Downstream dependency.** Stage 9 (cross-ply / quasi-iso laminate) uses one solid element per ply with a different `/SKEW/FIX` per ply. If stage 8 fails for a single $\theta_k$, stage 9 is blocked until the failure is diagnosed (likely either the `Iorth` flag interpretation or the skew-axis convention). Stage 11 (pseudo-woven mesoscale) has the same dependency at the tow level.

**Risk register.** Two specific failure modes to watch for, both documented in the audit. (i) `/PROP/TYPE14` `Iorth` flag interpretation: Radioss historically had two conventions for whether the orthotropic frame co-rotates with the element or stays fixed in space; for a small-strain linear-elastic test the two conventions are numerically indistinguishable, but the runner pins the convention by setting `Iframe=1` (co-rotational) explicitly. (ii) `/SKEW/FIX` definition by two vectors: Radioss documentation requires the two vectors $(\mathbf{A},\mathbf{B})$ where $\mathbf{B}$ is in the local 1-2 plane but need not be orthogonal to $\mathbf{A}$; the runner supplies orthonormal vectors anyway to remove ambiguity.

**What success looks like.** Seven completed implicit-static runs, one CSV row per run, $r(\theta_k)\le 0.01$ for $k=0\dots 6$, one Typst+CeTZ figure overlaying the seven FEM points on the smooth analytic curve, total wall-clock under 10 minutes on a laptop-class machine inside the Lima VM. With this, stage 9 starts.
