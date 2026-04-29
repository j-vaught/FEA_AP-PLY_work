# Stage 06 - Composite Failure Criteria Side-by-Side (Tsai-Wu / Hashin / Puck)

Author: J.C. Vaught
Date: 2026-04-29
Stage: 06 of 16
Pass criterion (master plan, Stage 6): "failure envelope matches WWFE-II"
within 5 percent of the strength values along the principal axes.

This stage is a *constitutive verification*. It is the simplest possible test
of OpenRadioss `/MAT/LAW25` plus `/FAIL/TSAIWU`, `/FAIL/HASHIN`, and
`/FAIL/PUCK` on solid bricks. Three runs are made on the same uniform UD ply
coupon under a sweep of in-plane biaxial stress paths in the
`(sigma_1, sigma_2)` plane plus pure shear `tau_12`. Each run uses one
failure card. The predicted failure envelope is plotted and compared to the
WWFE-II Part B (Kaddour and Hinton 2013) reference envelope for the same
criterion. The three criteria *should* give different envelopes; the goal is
to verify that OpenRadioss correctly implements all three, not that they
agree with each other.

The structure below follows the 10-section format used across all stages.

---

## 1. Goal and scope

Verify that OpenRadioss `/MAT/LAW25` + each of `/FAIL/TSAIWU`,
`/FAIL/HASHIN`, `/FAIL/PUCK`, run on a solid HEXA8 brick UD ply coupon,
reproduces the WWFE-II Part B reference failure envelope for that criterion
within 5 percent of the listed strength values along the principal stress
axes (uniaxial fiber tension, uniaxial fiber compression, uniaxial transverse
tension, uniaxial transverse compression, pure in-plane shear).

The six checked stress paths are.

1. Uniaxial longitudinal tension along the fiber axis ($\sigma_1 > 0$, all others zero) - must recover $X_T$ within 5 percent on each criterion.
2. Uniaxial longitudinal compression ($\sigma_1 < 0$) - must recover $X_C$ within 5 percent.
3. Uniaxial transverse tension ($\sigma_2 > 0$) - must recover $Y_T$ within 5 percent.
4. Uniaxial transverse compression ($\sigma_2 < 0$) - must recover $Y_C$ within 5 percent.
5. Pure in-plane shear ($\tau_{12}$, all normal stresses zero) - must recover $S_{12}$ within 5 percent.
6. Combined biaxial transverse-plus-shear path used to differentiate the three criteria, where Tsai-Wu, Hashin matrix, and Puck IFF mode A/B/C give *different* failure stresses; this is the "discriminator" point and is reported but not used as a pass gate (the three criteria are not expected to agree).

Out of scope. Through-thickness stresses ($\sigma_3, \tau_{13}, \tau_{23}$),
free-edge effects, ply rotation (Stage 8), laminate stacking (Stage 9),
post-failure damage softening (Stage 13), high-strain-rate effects (Stage
15). This is single-ply, single-element-thickness, quasi-static, in-plane.

Limit of the test. The implementation in `/FAIL/TSAIWU` and the LAW25
in-built Tsai-Wu surface use the standard 1971 Tsai-Wu form with
$F_{12}^*$ taken as $-0.5$ by default unless the user supplies it (LAW25
input field `F12`). Hashin in OpenRadioss is the 1980 form per Hashin's
*Journal of Applied Mechanics* paper. Puck in OpenRadioss follows the
1998 / 2002 Puck-Schurmann form for plane-stress IFF modes A, B, C with
the master fracture-plane angle approximated as $\theta_{fp}=53^\circ$
for compression-dominated transverse loading (the Puck recommended value).
Any deviation between OpenRadioss and the WWFE-II reference curve that
exceeds 5 percent on a principal axis is recorded as a finding and
referred back to the failure-card reference manual page.

## 2. Material selection - WWFE-II reference card

Selected material. **IM7/8552 carbon/epoxy unidirectional prepreg**, as
listed in the WWFE-II material table reproduced by Kaddour and Hinton 2013
*Journal of Composite Materials* 47(6-7), 925-966, "Maturity of 3D failure
criteria for fibre-reinforced composites: comparison between theories and
experiments: Part B of WWFE-II". The same material is also used by UofSC
for HVI testing (master plan section 2, table row "UofSC HVI material") so
this stage shares a material card with Stages 13-16 downstream.

Why IM7/8552 over T800/3900-2. Both are listed in WWFE-II. IM7/8552 is the
more frequently quoted card across WWFE-II Part B because it is the
benchmark used in the 3D failure-envelope figures (Kaddour-Hinton 2013,
Figures 5-12). It is also the material in the AP-PLY HVI panel (Vakili Rad
2020). T800/3900-2 is the LVI/CAI material (Kodagali 2023) and would be
substituted in Stage 14, not here.

WWFE-II Part B IM7/8552 card values (UD lamina). All values are the
WWFE-II tabulated entries; SI units.

| Property | Symbol | Value | Unit |
|---|---|---|---|
| Longitudinal Young's modulus | $E_1$ | 165.0 | GPa |
| Transverse Young's modulus | $E_2 = E_3$ | 9.0 | GPa |
| In-plane shear modulus | $G_{12} = G_{13}$ | 5.6 | GPa |
| Transverse shear modulus | $G_{23}$ | 2.8 | GPa |
| In-plane Poisson's ratio | $\nu_{12} = \nu_{13}$ | 0.34 | - |
| Transverse Poisson's ratio | $\nu_{23}$ | 0.50 | - |
| Longitudinal tensile strength | $X_T$ | 2560.0 | MPa |
| Longitudinal compressive strength | $X_C$ | 1590.0 | MPa |
| Transverse tensile strength | $Y_T$ | 73.0 | MPa |
| Transverse compressive strength | $Y_C$ | 185.0 | MPa |
| In-plane shear strength | $S_{12}$ | 90.0 | MPa |
| Density | $\rho$ | 1580.0 | kg/m^3 |

Citation. Kaddour, A. S. and Hinton, M. J. (2013), "Maturity of 3D failure
criteria for fibre-reinforced composites: comparison between theories and
experiments: Part B of WWFE-II", *Journal of Composite Materials*,
47(6-7), 925-966, Special Issue. Bibtex key
`KaddourHinton2013WWFEII` (in `references/test_progression_refs.bib`).

The same property set appears in the WWFE-II Part A (Kaddour, Hinton, Smith,
Li 2013) and is consistent with Soden, Hinton, and Kaddour 1998.

## 3. Geometry

Single UD ply coupon. The geometry is deliberately tiny because this is a
constitutive-card test, not a free-edge or structural test.

| Dimension | Value | Unit |
|---|---|---|
| In-plane length $L_1$ (fiber direction) | 20.0 | mm |
| In-plane width $L_2$ (transverse direction) | 20.0 | mm |
| Thickness $t$ (cured ply, master plan section 2) | 0.18 | mm |

Material orientation. Fiber direction along the global x-axis. Transverse
direction along y. Through-thickness along z. The `/PROP/TYPE14` brick
property's orthotropic frame is set by a `/SKEW/FIX` aligned with the global
axes (no rotation in this stage; rotation is Stage 8).

A hand drawing of the coupon for reviewer reference.

```
                     y (transverse, sigma_2)
                     ^
                     |
              +------+------+
             /|     /|     /|
            +------+------+ |
            | |    | |    | |
            | +----| +----| +
            |/     |/     |/
            +------+------+ ----> x (fiber, sigma_1)
                  /
                 /
                v z (thickness, 0.18 mm)
```

## 4. Mesh

Uniform 10 by 10 by 1 HEXA8 grid. Deliberately coarse because the goal is
to recover a uniform stress state across the coupon when uniform strain is
imposed at the edges, not to resolve any gradient.

| Mesh parameter | Value |
|---|---|
| Element type | HEXA8 (linear 8-node brick) |
| In-plane element count | 10 by 10 |
| Through-thickness elements | 1 |
| Total elements | 100 |
| Total nodes | 11 by 11 by 2 = 242 |
| Element size in-plane | 2.0 mm |
| Element size through-thickness | 0.18 mm |
| Aspect ratio | 11:1 in-plane vs through-thickness, acceptable because the test is in-plane only |

GMSH meshing is straightforward (transfinite quad of 10 by 10 extruded by
one element through thickness). The meshing script is generated by
`runner.py` section 3.1 below.

## 5. Material card and failure-card inputs

OpenRadioss material card. `/MAT/LAW25` (COMPSH) on `/PROP/TYPE14` solid
brick. LAW25 is selected (over LAW12 / LAW14) because LAW25 is the only
solid-composite card that accepts all three of `/FAIL/TSAIWU`,
`/FAIL/HASHIN`, `/FAIL/PUCK` (Altair LAW25 reference, AltairRadiossLAW25;
composite intro page, AltairRadiossCompositeIntro). LAW25 also computes
the Tsai-Wu index as its native yield surface, but for stage 6 we use the
explicit `/FAIL/TSAIWU` card so that the three criteria are applied in
exactly equivalent ways (each is a separate `/FAIL/...` block on the same
LAW25 base).

LAW25 input fields, IM7/8552 values (one card, shared across runs A/B/C).

```
/MAT/LAW25/1
IM7-8552 UD ply
# RHO_I
1580.0e-9                                     # kg/mm^3 (SI on mm-N-MPa-ms unit set)
# E11    E22    nu12    Iform                Eps_f1 Eps_f2
165000.0 9000.0 0.34    1                    0.0    0.0
# G12    G23    G31     Eps_t1 Eps_m1 d_max
5600.0   2800.0 5600.0  0.0    0.0    1.0
# Wp_max  Wp_ref Iflawp  ... (CRASURV plasticity off; pure elastic until /FAIL)
0.0       0.0    0
# sigma_1y_t  b_t  n_t  sigma_1_max_t  c_t  EPS_1t_max     # tensile fiber
0.0           0.0  0.0  0.0            0.0  1.0
# sigma_1y_c  b_c  n_c  sigma_1_max_c  c_c  EPS_1c_max     # compressive fiber
0.0           0.0  0.0  0.0            0.0  1.0
# Tsai-Wu on the LAW25 line (we set them but use /FAIL/TSAIWU for run A,
# and /FAIL/HASHIN for run B, /FAIL/PUCK for run C; LAW25's internal Tsai-Wu
# is left inactive by setting Iflawt=0)
# F11... not used here; keep zero
```

(The exact field order follows `help.altair.com/hwsolvers/rad/topics/solvers/rad/law25_composite_material_r.htm`. The `runner.py` produces the literal deck text; the table above lists the input numbers, not the deck syntax.)

Per-criterion failure-card inputs.

### 5.1 Run A - `/FAIL/TSAIWU`

The Tsai-Wu surface in `/FAIL/TSAIWU` is the 1971 quadratic form

$$
F_1\sigma_1 + F_2\sigma_2 + F_{11}\sigma_1^2 + F_{22}\sigma_2^2 + F_{66}\tau_{12}^2 + 2 F_{12} \sigma_1 \sigma_2 = 1
$$

with the standard strength definitions

$$
F_1 = \frac{1}{X_T} - \frac{1}{X_C},\quad F_2 = \frac{1}{Y_T} - \frac{1}{Y_C},\quad F_{11} = \frac{1}{X_T X_C},\quad F_{22} = \frac{1}{Y_T Y_C},\quad F_{66} = \frac{1}{S_{12}^2},\quad F_{12} = F_{12}^* \sqrt{F_{11} F_{22}}
$$

with $F_{12}^* \in [-1, 1]$. The default value used in WWFE-II reference
plots is $F_{12}^* = -0.5$ (Tsai-Wu's own recommendation for typical CFRP).

OpenRadioss `/FAIL/TSAIWU` inputs.

| Field | Value | Unit |
|---|---|---|
| `Sigma_1t` ($X_T$) | 2560.0 | MPa |
| `Sigma_1c` ($X_C$) | 1590.0 | MPa |
| `Sigma_2t` ($Y_T$) | 73.0 | MPa |
| `Sigma_2c` ($Y_C$) | 185.0 | MPa |
| `Sigma_12` ($S_{12}$) | 90.0 | MPa |
| `F12*` (interaction) | -0.5 | - |
| `Ifail` (erosion flag) | 1 (delete element on failure) | - |

### 5.2 Run B - `/FAIL/HASHIN`

The Hashin 1980 plane-stress criterion separates into four sub-criteria.
Each gives an index $f_*$; failure occurs when $f_* \ge 1$.

Fiber tensile failure ($\sigma_1 > 0$).

$$
f_{f}^{T} = \left(\frac{\sigma_1}{X_T}\right)^2 + \left(\frac{\tau_{12}}{S_{12}}\right)^2
$$

Fiber compressive failure ($\sigma_1 < 0$).

$$
f_{f}^{C} = \left(\frac{\sigma_1}{X_C}\right)^2
$$

Matrix tensile failure ($\sigma_2 > 0$).

$$
f_{m}^{T} = \left(\frac{\sigma_2}{Y_T}\right)^2 + \left(\frac{\tau_{12}}{S_{12}}\right)^2
$$

Matrix compressive failure ($\sigma_2 < 0$).

$$
f_{m}^{C} = \left(\frac{\sigma_2}{2 S_{23}}\right)^2 + \left[\left(\frac{Y_C}{2 S_{23}}\right)^2 - 1\right]\frac{\sigma_2}{Y_C} + \left(\frac{\tau_{12}}{S_{12}}\right)^2
$$

with $S_{23}$ the transverse (out-of-plane) shear strength. For IM7/8552
(WWFE-II Part B) we use $S_{23} = 50.0$ MPa (typical CFRP, also used by
Pinho et al. 2006 in their Hashin verification).

OpenRadioss `/FAIL/HASHIN` inputs.

| Field | Value | Unit |
|---|---|---|
| `Sigma_1t` ($X_T$) | 2560.0 | MPa |
| `Sigma_1c` ($X_C$) | 1590.0 | MPa |
| `Sigma_2t` ($Y_T$) | 73.0 | MPa |
| `Sigma_2c` ($Y_C$) | 185.0 | MPa |
| `Sigma_12` ($S_{12}$) | 90.0 | MPa |
| `Sigma_23` ($S_{23}$) | 50.0 | MPa |
| `Ifail` | 1 | - |

### 5.3 Run C - `/FAIL/PUCK`

Puck 1998 / 2002 separates fiber failure (FF) from inter-fiber failure (IFF)
with three IFF modes A, B, C.

Fiber failure (Puck 1998 magnification-corrected form, which for the
in-plane sweep at no through-thickness load reduces to maximum-stress).

$$
f_{FF}^{T} = \frac{\sigma_1}{X_T}\quad (\sigma_1 > 0),\qquad f_{FF}^{C} = \frac{|\sigma_1|}{X_C}\quad (\sigma_1 < 0).
$$

Inter-fiber failure modes. Define the transverse normal stress
$\sigma_n = \sigma_2$ and the in-plane shear $\tau_{nt} = \tau_{12}$ (we
fix the fracture plane to the master plane $\theta_{fp}=0$ for matrix
tension and $\theta_{fp}=53^\circ$ for matrix compression, the
Puck-recommended values for plane-stress).

Mode A: $\sigma_n \ge 0$ (transverse tension - matrix cracks on fracture plane normal to $\sigma_2$).

$$
f_{IFF,A} = \sqrt{\left[\left(\frac{1}{Y_T}-\frac{p_{\perp\parallel}^{(+)}}{S_{12}}\right)\sigma_n\right]^2 + \left(\frac{\tau_{nt}}{S_{12}}\right)^2} + \frac{p_{\perp\parallel}^{(+)}}{S_{12}}\sigma_n
$$

Mode B: $\sigma_n < 0$ and $|\sigma_n / \tau_{nt}| \le R_{\perp\perp}^A / |\tau_{nt,c}|$ (small transverse compression with shear; failure still on the master plane).

$$
f_{IFF,B} = \frac{1}{S_{12}}\left[\sqrt{\tau_{nt}^2 + (p_{\perp\parallel}^{(-)} \sigma_n)^2} + p_{\perp\parallel}^{(-)} \sigma_n \right]
$$

Mode C: $\sigma_n < 0$ and large negative $\sigma_n$ (transverse compression dominant, fracture plane rotates to about $\theta_{fp}=53^\circ$).

$$
f_{IFF,C} = \left[\left(\frac{\tau_{nt}}{2(1+p_{\perp\perp}^{(-)})S_{12}}\right)^2 + \left(\frac{\sigma_n}{Y_C}\right)^2\right]\cdot \frac{Y_C}{-\sigma_n}
$$

with the standard Puck inclination parameters

$$
p_{\perp\parallel}^{(+)} = 0.30,\quad p_{\perp\parallel}^{(-)} = 0.25,\quad p_{\perp\perp}^{(+)} = 0.20,\quad p_{\perp\perp}^{(-)} = 0.25,
$$

these being the WWFE-II / Puck-2002 recommended CFRP defaults.

OpenRadioss `/FAIL/PUCK` inputs.

| Field | Value | Unit |
|---|---|---|
| `Sigma_1t` ($X_T$) | 2560.0 | MPa |
| `Sigma_1c` ($X_C$) | 1590.0 | MPa |
| `Sigma_2t` ($Y_T$) | 73.0 | MPa |
| `Sigma_2c` ($Y_C$) | 185.0 | MPa |
| `Sigma_12` ($S_{12}$) | 90.0 | MPa |
| `p_pos_perp_par` | 0.30 | - |
| `p_neg_perp_par` | 0.25 | - |
| `p_pos_perp_perp` | 0.20 | - |
| `p_neg_perp_perp` | 0.25 | - |
| `Ifail` | 1 | - |

## 6. Boundary conditions and load sweep

Implicit static, `/IMPL/LINEAR` followed by `/IMPL/NONLIN` for the path
that crosses the failure surface (load proportional from 0 to a target
high enough to fail the coupon for every direction in the sweep).

Each "biaxial path" $(c_1, c_2, c_6)$ is a unit vector in
$(\sigma_1, \sigma_2, \tau_{12})$ space. For each path, the deck applies
the corresponding *strain* boundary condition

$$
\bar\varepsilon_1 = c_1 \cdot \alpha,\quad \bar\varepsilon_2 = c_2 \cdot \alpha,\quad \bar\gamma_{12} = c_6 \cdot \alpha
$$

via prescribed displacements on the four edge faces, with $\alpha$ ramped
linearly from 0 to a target high enough that every criterion fails. The
recorded "failure stress" at the path is the homogenized
$\bar\sigma = \mathbf{C}\bar\varepsilon$ (from `/MAT/LAW25` elastic
constants) at the load step where the failure index first reaches 1.

The discrete sweep used here is 24 paths around the
$(\sigma_1, \sigma_2)$ envelope plus 8 paths in the
$(\sigma_2, \tau_{12})$ envelope.

| Path | $c_1$ | $c_2$ | $c_6$ | Description |
|---|---|---|---|---|
| 1 | 1 | 0 | 0 | Pure longitudinal tension - recovers $X_T$ |
| 2 | -1 | 0 | 0 | Pure longitudinal compression - recovers $X_C$ |
| 3 | 0 | 1 | 0 | Pure transverse tension - recovers $Y_T$ |
| 4 | 0 | -1 | 0 | Pure transverse compression - recovers $Y_C$ |
| 5 | 0 | 0 | 1 | Pure in-plane shear - recovers $S_{12}$ |
| 6 | $\cos(15^\circ)$ | $\sin(15^\circ)$ | 0 | (LT) biaxial 15 deg in $(\sigma_1, \sigma_2)$ |
| 7 | $\cos(30^\circ)$ | $\sin(30^\circ)$ | 0 | 30 deg |
| 8-23 | sweep $\theta \in \{45, 60, 75, 105, 120, ..., 345\}^\circ$ in 15-deg steps | - | 0 | Full $(\sigma_1, \sigma_2)$ envelope |
| 24 | 0 | $\cos(\phi)$ | $\sin(\phi)$ | $\phi \in \{0, 30, 60, 90, 120, 150, 180, 210\}^\circ$ in $(\sigma_2, \tau_{12})$ envelope (eight paths) |

Boundary conditions (uniform-strain corner-driven). The four in-plane edges
of the brick are kinematically constrained to translate as rigid edges
(this is achieved with `/RBE2` rigid-edge constraints on each face,
slaved to a single master node per face, the master-node displacement set
by `/IMPDISP`). The through-thickness top and bottom faces have only z
constrained at three points to suppress rigid-body motion (the SP-three
rule: one node fully fixed, one node z+y fixed, one node z fixed).

This gives a uniform stress state inside the brick to within 0.5 percent
once the coupon is at least 5x5 elements (it is 10x10), which is well
inside the budget for a constitutive-card test. The stress is read from a
`/TH/PART` time-history block on the brick part.

Termination of each run. The starter / engine pair runs until either
(i) the failure index of the criterion reaches 1.0 in any element (the
`/FAIL/...` `Ifail=1` flag deletes the element, which `runner.py` detects
through the `/H3D/ELEM/SCAL` damage variable), or
(ii) the proportional load reaches a hard cap (e.g., 5 times $X_T$ along
the fiber path) - this protects against criteria that never trigger
because of an input-mistake in the failure card.

## 7. Reference solution - analytic envelopes for comparison

The pass criterion is "envelope matches WWFE-II Part B reference within 5
percent on each principal axis". The principal-axis points are
$(\sigma_1=X_T, 0, 0)$, $(\sigma_1=-X_C, 0, 0)$, $(0, Y_T, 0)$,
$(0, -Y_C, 0)$, $(0, 0, S_{12})$. Off-axis points (combined biaxial) are
plotted but not gated, because the three criteria are not expected to
agree off-axis - the WWFE-II found 30-200 percent disagreement there.

Closed-form analytic envelopes used by `runner.py` for the "WWFE-II
reference curve" overlay are computed below from primary-source equations,
then evaluated on the 32-path sweep of section 6.

### 7.1 Tsai-Wu (Tsai and Wu 1971)

The failure surface is a quadric in $(\sigma_1, \sigma_2, \tau_{12})$,

$$
F_1\sigma_1 + F_2\sigma_2 + F_{11}\sigma_1^2 + F_{22}\sigma_2^2 + F_{66}\tau_{12}^2 + 2 F_{12}\sigma_1 \sigma_2 = 1.
$$

For a path direction $(c_1, c_2, c_6)$ with $\sigma_i = c_i R$ where $R$
is the load magnitude at failure, this is a quadratic in $R$,

$$
[F_{11} c_1^2 + F_{22} c_2^2 + F_{66} c_6^2 + 2 F_{12} c_1 c_2] R^2 + [F_1 c_1 + F_2 c_2] R - 1 = 0,
$$

solved with the standard quadratic formula taking the positive root.
$F_{12}^* = -0.5$ unless overridden.

### 7.2 Hashin (Hashin 1980)

Each path is scaled by the smallest $R$ at which any of the four
sub-criteria reaches 1. For a path $(c_1, c_2, c_6)$, define
$R_{ff,T}$, $R_{ff,C}$, $R_{mf,T}$, $R_{mf,C}$ as the positive root of
each sub-criterion (with the natural sign convention applied to
$\sigma_1$, $\sigma_2$).

$$
R_{ff,T} = \left[\left(\frac{c_1}{X_T}\right)^2 + \left(\frac{c_6}{S_{12}}\right)^2\right]^{-1/2}\quad\text{(applies if } c_1 R > 0\text{)},
$$

$$
R_{ff,C} = \frac{X_C}{|c_1|}\quad\text{(applies if } c_1 R < 0\text{)},
$$

$$
R_{mf,T} = \left[\left(\frac{c_2}{Y_T}\right)^2 + \left(\frac{c_6}{S_{12}}\right)^2\right]^{-1/2}\quad\text{(applies if } c_2 R > 0\text{)},
$$

$$
R_{mf,C}\text{: positive root of }\left(\frac{c_2}{2 S_{23}}\right)^2 R^2 + \left[\left(\frac{Y_C}{2 S_{23}}\right)^2 - 1\right]\frac{c_2 R}{Y_C} + \left(\frac{c_6 R}{S_{12}}\right)^2 = 1\quad\text{(applies if } c_2 R < 0\text{)}.
$$

The path failure load is $R = \min(R_{ff,T}, R_{ff,C}, R_{mf,T}, R_{mf,C})$
over the active sub-criteria for that path.

### 7.3 Puck (Puck-Schurmann 1998, revised 2002)

Puck FF reduces to maximum-stress for the in-plane sweep with no
through-thickness load.

$$
R_{FF,T} = \frac{X_T}{c_1}\quad(c_1 > 0),\qquad R_{FF,C} = \frac{X_C}{|c_1|}\quad(c_1 < 0).
$$

Puck IFF mode A is active when $c_2 R \ge 0$. Substituting
$\sigma_n = c_2 R$, $\tau_{nt} = c_6 R$ in the mode-A expression and
setting $f_{IFF,A} = 1$,

$$
\sqrt{\left[\left(\frac{1}{Y_T} - \frac{p_{\perp\parallel}^{(+)}}{S_{12}}\right) c_2 R\right]^2 + \left(\frac{c_6 R}{S_{12}}\right)^2} + \frac{p_{\perp\parallel}^{(+)} c_2 R}{S_{12}} = 1,
$$

a quadratic in $R$. The positive root is $R_{IFF,A}$.

Puck IFF mode B is active when $c_2 R < 0$ and the
shear-to-compression ratio satisfies
$|\tau_{nt}/\sigma_n| \ge R_{\perp\perp}^A/(|\tau_{nt,c}|)$. The mode-B
failure stress is

$$
R_{IFF,B} = \frac{S_{12}}{\sqrt{c_6^2 + (p_{\perp\parallel}^{(-)} c_2)^2} + p_{\perp\parallel}^{(-)} c_2}.
$$

Puck IFF mode C is active when $c_2 R < 0$ and the ratio is below the
mode-B threshold. The mode-C failure stress comes from setting
$f_{IFF,C} = 1$ and solving the resulting quadratic for $R$.

The path failure load is $R = \min(R_{FF,T}, R_{FF,C}, R_{IFF,A},
R_{IFF,B}, R_{IFF,C})$ over the active modes.

The mode-switch threshold between modes A, B, C is the Puck recommendation
based on the slope of the fracture plane: mode A for $\sigma_n \ge 0$,
mode B for small negative $\sigma_n$ with shear, mode C for large negative
$\sigma_n$. The exact angular threshold depends on the inclination
parameters; for $p_{\perp\parallel}^{(-)} = 0.25$ the
mode-B-to-mode-C transition lies at
$|\tau_{nt}/\sigma_n| \approx 0.4$.

### 7.4 WWFE-II Part B reference points

The WWFE-II Part B paper (Kaddour-Hinton 2013) reproduces the predicted
envelopes for each criterion in their Figures 5-7 for the IM7/8552 system.
The principal-axis pass-gate values are exactly the strengths in section 2.
The reference for off-axis (combined-loading) points is the numerical
envelope produced by analytic equations 7.1-7.3 of this document, since
WWFE-II Figures 5-7 do not tabulate every point.

## 8. OpenRadioss deck skeleton (single-stage)

A single shared starter deck plus three engine decks (one per criterion).
Filenames `stage06_TSAIWU_0000.rad`, `stage06_HASHIN_0000.rad`,
`stage06_PUCK_0000.rad`.

```
#--- Stage 06 - LAW25 + /FAIL/<criterion> on a 10x10x1 HEXA8 UD ply
#--- Units: mm, ms, kg, MPa  (mass density 1.580e-6 kg/mm^3)
#---
/BEGIN
Stage06 LAW25 + /FAIL/{{CRITERION}}
2025
mm                          ms                            kg                            MPa
#---
/UNIT/1                     0
/MAT/LAW25/1
IM7-8552 UD ply (WWFE-II Part B)
1.580e-6
165000.0  9000.0  0.34  1   0.0   0.0
5600.0    2800.0  5600.0 0.0 0.0  1.0
0.0       0.0     0
0.0  0.0  0.0  0.0  0.0  1.0
0.0  0.0  0.0  0.0  0.0  1.0
#---
/PROP/TYPE14/1
1                                 1                                 0                                 0                                 0                                 0                                 0
#---
/SKEW/FIX/1
SKEW_GLOBAL                                                                                                                                                                                                            
0.0  0.0  0.0
1.0  0.0  0.0
0.0  1.0  0.0
#---
{% if CRITERION == 'TSAIWU' %}
/FAIL/TSAIWU/1
2560.0   1590.0   73.0   185.0   90.0   -0.5
1
{% elif CRITERION == 'HASHIN' %}
/FAIL/HASHIN/1
2560.0   1590.0   73.0   185.0   90.0   50.0
1
{% elif CRITERION == 'PUCK' %}
/FAIL/PUCK/1
2560.0   1590.0   73.0   185.0   90.0
0.30   0.25   0.20   0.25
1
{% endif %}
#---
/PART/1
COUPON_UD                                                                                                                                                                                                              
1   1   0   0   0   0   0   0
#---
*INCLUDE 'stage06_mesh.inc'      # 11x11x2 nodes, 100 HEXA8 elements
#--- BCs: rigid-edge MPC + master node IMPDISP per the path-vector
*INCLUDE 'stage06_bcs_path_{{PATH_ID}}.inc'
#---
/IMPL/LINEAR/1
/IMPL/NONLIN/1
/IMPL/DT/1
1.0e-3   1.0e-5   1.0e-2   100   1
/IMPL/SOLVER/1
1   1   0   0
#---
/TH/PART/1
COUPON_UD_TH
1
SX SY SZ SXY SXZ SYZ
#---
/H3D/ELEM/SCAL/DAMA           # damage scalar (failure-index proxy)
#---
/END
```

The `runner.py` script generates one such deck per (criterion, path)
combination, runs each in a Lima / Apptainer OpenRadioss container, and
records the load magnitude $R$ at which the H3D damage scalar crosses 1.0
for the first time (i.e., the failure stress along that path).

## 9. Validation success criterion - exact pass gate

For each criterion C in {TSAIWU, HASHIN, PUCK}.

1. Principal-axis test (HARD GATE).

   The five principal-axis strength values produced by OpenRadioss for criterion C must each match the WWFE-II Part B value from section 2 within 5 percent.

   $$
   \frac{|R_C^{\text{OpenRadioss, axis}} - R_C^{\text{WWFE-II, axis}}|}{R_C^{\text{WWFE-II, axis}}} \le 0.05
   $$

   for axis in {LT, LC, TT, TC, S}, with values $X_T = 2560$, $X_C = 1590$, $Y_T = 73$, $Y_C = 185$, $S_{12} = 90$ MPa.

2. Self-consistency of the criterion (HARD GATE).

   Tsai-Wu: along $\sigma_1 = X_T$, $\sigma_2 = 0$, $\tau_{12} = 0$ the index must equal 1.0 within 1 percent. (Internal sanity check on the Tsai-Wu surface.)

   Hashin: in the LT path, fiber-tensile sub-criterion $f_f^T$ must reach 1.0 first; in the TC path, matrix-compression sub-criterion $f_m^C$ must reach 1.0 first; etc. The mode flag from OpenRadioss must match the expected sub-criterion.

   Puck: along the TT path (transverse tension), mode A must trigger first; along the TC path (transverse compression), mode C must trigger first; along the LT path, FF tensile must trigger first. The Puck-mode flag must match.

3. Off-axis envelope shape (SOFT GATE - reported, not pass-gated).

   The full 32-path failure envelope for each criterion is plotted alongside the analytic envelope from section 7. Visually the OpenRadioss curve and the analytic curve should overlap to within plotting linewidth. RMS error over all 32 paths reported in the runner output table.

4. Cross-criterion comparison (REPORTED, not pass-gated).

   The three envelopes are plotted on a single figure. The well-known result that Tsai-Wu predicts a smooth ellipse, Hashin produces a piecewise-elliptic envelope with corners at sub-criterion transitions, and Puck shows a compressive bulge in the TC quadrant due to mode C, is verified visually. This is *not* a pass gate; it is the expected behavior that the three criteria differ.

If any of the principal-axis tests (gate 1) fails by more than 5 percent
the stage fails and the failure-card input fields are re-checked against
the `help.altair.com` reference page. A 1-2 percent discrepancy is
acceptable and is attributed to the sampling of strain-driven proportional
loading hitting the failure surface at a non-zero load step.

## 10. Run procedure and output artifacts

1. From the project root run `python tests/stage_06_composite_failure_criteria/runner.py --build`. This generates `stage06_mesh.inc`, the 32 BC include files, and the three starter decks.

2. `python tests/stage_06_composite_failure_criteria/runner.py --solve` runs all 96 (= 3 criteria * 32 paths) decks in series through the Lima Apptainer OpenRadioss container. Wall-clock time is tens of seconds per deck (single HEXA8 brick, implicit-static linear at most 100 steps), totalling ~30-60 minutes on an M-series Mac.

3. `python tests/stage_06_composite_failure_criteria/runner.py --analyze` parses the 96 H3D output files for the failure-load $R$ along each path, and writes one CSV per criterion to `tests/stage_06_composite_failure_criteria/out/`.

4. `python tests/stage_06_composite_failure_criteria/runner.py --plot` writes Typst data files for CeTZ to `figures/stage_06/` and emits the principal-axis pass-gate table to stdout.

Output artifacts (final).

| Artifact | Location | Format |
|---|---|---|
| Failure load per (criterion, path) | `out/envelope_TSAIWU.csv`, `out/envelope_HASHIN.csv`, `out/envelope_PUCK.csv` | CSV |
| Pass-gate table | `out/pass_gate.csv` and stdout | CSV / table |
| Cross-criterion plot data | `out/sigma1_sigma2_envelope.csv`, `out/sigma2_tau12_envelope.csv` | CSV |
| Analytic-vs-OpenRadioss overlay plot | `figures/stage_06/envelope_sigma1_sigma2.pdf` (Typst+CeTZ from `figures/stage_06/envelope.typ`) | PDF |
| Run log | `out/run_log.txt` | text |

Citations (consolidated at end of run log).

- Tsai, S. W. and Wu, E. M. (1971), bibtex `TsaiWu1971`.
- Hashin, Z. (1980), bibtex `Hashin1980`.
- Puck, A. and Schurmann, H. (1998), bibtex `PuckSchurmann1998`.
- Puck, A. and Schurmann, H. (2002), bibtex `PuckSchurmann2002`.
- Hinton, M. J., Kaddour, A. S., Soden, P. D. (2002), bibtex `HintonKaddourSoden2002`.
- Kaddour, A. S. and Hinton, M. J. (2013), bibtex `KaddourHinton2013WWFEII`.
- Soden, P. D., Hinton, M. J., Kaddour, A. S. (1998), bibtex `SodenHintonKaddour1998`.
- Pinho, S. T., Iannucci, L., Robinson, P. (2006), bibtex `PinhoPartI2006`, `PinhoPartII2006` (for Hashin / Puck FE-implementation reference).
- Altair Engineering (2025), `/MAT/LAW25` reference, bibtex `AltairRadiossLAW25`.
- Altair Engineering (2025), `/FAIL/HASHIN` reference, bibtex `AltairRadiossFailHashin`.
- Altair Engineering (2025), composite material introduction, bibtex `AltairRadiossCompositeIntro`.

End of spec.
