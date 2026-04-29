# Stage 02 — Cantilever Beam, Large-Deflection Geometric Nonlinearity

**Author.** J.C. Vaught
**Date.** 2026-04-29
**Stage.** 02 of 16 (FEA_AP-PLY single-tool OpenRadioss progression).
**Type.** Verification (V) — finite element solution against closed-form elastica.
**Audit status carried forward.** MARGINAL. The
`references/openradioss_endtoend_audit.md` (Stage 2 row) flags two
caveats. (i) Implicit static (`/IMPL/NONLIN`) requires a MUMPS-linked
build of OpenRadioss; the prebuilt explicit-only binaries cannot run
this stage. (ii) Quasi-static explicit with mass scaling (`/IMPL/QSTAT`
or pure explicit + dynamic relaxation `/DYREL`) is the
documented fallback. This spec commits to the implicit-MUMPS path
and documents the fallback in §10.

---

## 1. Stage summary

The goal is to isolate **geometric nonlinearity** in the FEM kernel.
The constitutive law is kept linear elastic so that any deviation from
Bisshopp-Drucker 1945 is attributable purely to the geometric kernel
(updated- or total-Lagrangian stress update, large rotations, follower
load handling) and not to plasticity, damage, or contact.

A long, slender, rectangular cantilever is loaded at its free tip with
a concentrated dead load $P$ acting in the global $-y$ direction. As $P$
increases, the tip rotates through angles approaching $\pi/2$, and the
tip translates substantially in both $-y$ (sag) and $-x$ (foreshortening
toward the clamp). The closed-form solution to the planar elastica was
first tabulated by Bisshopp and Drucker (1945, *Quart. Appl. Math.*)
in terms of complete and incomplete elliptic integrals of the first
and second kinds.

The pass criterion is that the FEM tip displacement components match
the Bisshopp-Drucker tabulated values within 2% at three sample points
of the dimensionless load
$\alpha = PL^2 / (EI)$, specifically $\alpha \in \{1, 3, 5\}$. Internally
the runner sweeps a finer set of $\alpha$ values for plotting purposes,
but only those three drive pass/fail.

This stage is the second smoke test of the OpenRadioss toolchain (the
first is Stage 01 linear-elastic 3-/4-pt bending). If it passes, we
have evidence that `/IMPL/NONLIN` is correctly wired, that the
MUMPS-linked build is functional, that follower-vs-dead load handling
on a tip surface is behaving as documented, and that HEXA8 with the
chosen `Ismstr` formulation does not lock under pure bending.

---

## 2. Geometry

A straight prismatic cantilever, occupying the rectangular box
$\{(x,y,z) : 0 \le x \le L,\ -b/2 \le y \le b/2,\ -h/2 \le z \le h/2\}$
in the unloaded reference configuration.

| Symbol | Value | Description |
|---|---|---|
| $L$ | $1.000\ \text{m}$ | beam length along $x$ |
| $b$ | $25.0\ \text{mm} = 0.0250\ \text{m}$ | width along $y$ |
| $h$ | $3.0\ \text{mm} = 0.0030\ \text{m}$ | thickness along $z$ |
| $L/h$ | $333$ | slenderness; deep in Euler-Bernoulli regime |
| $L/b$ | $40$ | width-aspect; rules out lateral-torsional buckling for this load direction |
| $A$ | $7.50 \times 10^{-5}\ \text{m}^2$ | cross-section area $bh$ |
| $I_{zz}$ | $5.625 \times 10^{-11}\ \text{m}^4$ | second moment about $z$, $bh^3/12$ |
| $EI$ | $11.25\ \text{N}\,\text{m}^2$ | bending stiffness about $z$ |

**Loading axis.** The bending axis is $z$ (out-of-plane). The tip load
acts in the global $-y$ direction. Because $b \gg h$, the beam is
laterally stiff and the deformation is planar in the $x$-$y$ plane,
exactly as Bisshopp-Drucker assumed.

**Why these dimensions.** The reference values $L^2/EI = 1/11.25$
mean that $\alpha = 1$ corresponds to $P = 11.25\ \text{N}$ and
$\alpha = 5$ corresponds to $P = 56.25\ \text{N}$. These are
tractable explicit-quasistatic-fallback magnitudes (no risk of
spurious dynamic ringing) and small enough that surface tractions
on a 25 mm $\times$ 3 mm tip face yield maximum tractions
of $7.5 \times 10^5\ \text{Pa}$, well below any onset of yield in
the assumed material (see §5). This protects the linear-elastic
assumption at the largest load case.

---

## 3. Mesh

Solid hexahedral elements only, per master plan §1.

**Baseline mesh.**

| Direction | Count | Element edge length |
|---|---|---|
| Along $x$ (length, $L=1.0\ \text{m}$) | 80 | $12.5\ \text{mm}$ |
| Across section direction 1 | 4 | $0.750\ \text{mm}$ |
| Across section direction 2 | 6 | $4.167\ \text{mm}$ |
| **Total elements** | **1920** | HEXA8 |

The brief mandates $\ge 40$ along length, $\ge 5$ across width,
$\ge 3$ through thickness; the baseline exceeds all three.

**Element aspect ratio.** Worst case $12.5 / 0.75 = 16.7$ along $x$
relative to the through-thickness direction. This is at the edge of
what HEXA8 tolerates in pure bending. The two mitigations are
discussed in §6 (formulation flag) and §10 (risks). We accept the
12.5 ratio because the through-thickness direction sees the strongest
bending strain gradient and benefits more from refinement than the
along-length direction.

**Mesh-convergence study.** Three meshes for the convergence check
at $\alpha = 5$ (largest deflection, where any locking artifact is
most visible):

| Tag | $N_x \times N_y \times N_z$ | Total | Edge $z$ |
|---|---|---|---|
| coarse | $40 \times 5 \times 3$ | 600 | 1.0 mm |
| baseline | $80 \times 4 \times 6$ | 1920 | 0.75 mm |
| fine | $160 \times 8 \times 4$ | 5120 | 0.75 mm |

Convergence is declared if $|\delta_y^{\text{fine}} - \delta_y^{\text{baseline}}|
/ \delta_y^{\text{fine}} < 0.5\%$ at $\alpha = 5$. The baseline mesh is
used for the pass/fail comparison against Bisshopp-Drucker.

**Mesh generator.** GMSH Python API, structured transfinite hex mesh
on the box. Output `.msh` $\to$ meshio $\to$ Abaqus `.inp` $\to$
OpenRadioss via `inp2rad`.

---

## 4. Boundary conditions and loading

### 4.1 Clamp at $x = 0$

All nodes on the face $x = 0$ are fixed in all six DOF. In a solid
element model only the three translational DOF exist on each node,
so this is implemented as $u_x = u_y = u_z = 0$ on every node of the
clamped face. With $\ge 3$ elements through thickness and $\ge 5$
across width, the resulting kinematic constraint is sufficient to
prevent rigid-body motion and to reproduce the engineering "fixed end"
boundary condition (no rotation of the cross-section at $x=0$).

OpenRadioss card. `/BCS/<grnod_id>` referencing a `/GRNOD/NODE` group
that lists every node on the face $x=0$. All three DOFs are blocked.

### 4.2 Tip load at $x = L$

Applied as a uniform pressure on the tip face $x = L$ acting in the
global $-y$ direction. Total resultant force is $P$; pressure
magnitude is $p = P / A_{\text{tip}}$ where $A_{\text{tip}} = bh =
7.5 \times 10^{-5}\ \text{m}^2$.

This is implemented as a **dead load** (direction fixed in the global
frame, *not* following the rotating cross-section). Bisshopp-Drucker
1945 explicitly assumes a dead vertical load; using a follower load
would not reproduce their tabulation (this is a documented pitfall
in elastica benchmarks; see Belendez et al. 2002).

OpenRadioss card. `/PLOAD` (pressure load) or `/CLOAD` (concentrated
load) on the tip face. `/PLOAD` is preferred because it spreads the
load uniformly over the tip face, avoiding stress concentrations at a
single node that would corrupt the local strain at the tip and bias
the displacement readout.

**Direction-fixed flag.** `/PLOAD` in OpenRadioss applies pressure
along the **face normal in the deformed configuration** by default
(follower behavior). For a true dead load, the load must be applied
as a constant traction with components fixed in the global frame.
The documented OpenRadioss path is `/CLOAD` distributed across the
tip-face nodes via `/GRNOD/NODE` and a function `/FUNCT` that scales
linearly with simulation time. Each node carries $P/N_{\text{tip}}$,
where $N_{\text{tip}}$ is the count of nodes on the tip face.

### 4.3 Load history

A pseudo-time ramp from $0$ to $1$ across the implicit static run.
Each implicit step delivers a fraction of $P$. The time function is
$f(t) = t$ for $t \in [0, 1]$. Step size is set adaptively
(see §6, `/IMPL/DT`).

### 4.4 Sweep over $\alpha$

The full $\alpha$ sweep is achieved by running multiple jobs, one per
$\alpha$ value, each with its own scaled $P$. Pass/fail uses
$\alpha \in \{1, 3, 5\}$. For plotting, the runner additionally
executes $\alpha \in \{0.5, 1, 2, 3, 4, 5\}$, expecting the linear
limit at $\alpha = 0.5$ to recover the Euler-Bernoulli small-deflection
tip deflection $\delta_y / L = \alpha / 3$ to better than 0.5%.

---

## 5. Material card

**Material.** Linear elastic spring steel. No plasticity, no damage,
no rate dependence — geometric nonlinearity is the only nonlinearity
this stage exercises.

| Property | Value | Units |
|---|---|---|
| Density $\rho$ | $7850$ | $\text{kg}/\text{m}^3$ |
| Young's modulus $E$ | $200 \times 10^{9}$ | $\text{Pa}$ |
| Poisson's ratio $\nu$ | $0.30$ | - |
| Yield strength (informative; not enforced) | $> 1\ \text{GPa}$ | $\text{Pa}$ |

The yield strength is recorded for context only — at $\alpha = 5$ the
peak fiber stress is approximately $|\sigma_{xx}|_{\max} = M_{\max}
(h/2) / I_{zz} = P L (h/2) / I_{zz} = 56.25 \cdot 1.0 \cdot 1.5 \times
10^{-3} / 5.625 \times 10^{-11} \approx 1.5\ \text{GPa}$. Spring steel
($\sigma_y \approx 1.0$-$1.5$ GPa for high-grade music wire) is at
the upper edge of its elastic envelope at the most extreme load
case, but for the verification purpose we run the material as
indefinitely linear. This is consistent with Bisshopp-Drucker who
assume linear constitutive behavior throughout.

OpenRadioss cards.

```
/MAT/LAW1/<mid>/<unit_id>
   spring_steel
   <rho>
   <E>   <nu>
```

with `mid = 1`, `rho = 7850`, `E = 2.0e11`, `nu = 0.30`.

The `/UNIT/<unit_id>` card establishes SI base units (kg, m, s, Pa).

---

## 6. OpenRadioss deck skeleton

The deck has two files. The starter file `cantilever_large_0000.rad`
declares geometry, material, properties, BCs, loads, and the
implicit nonlinear analysis cards. The engine file
`cantilever_large_0001.rad` controls the time stepping and outputs.

### 6.1 Starter `cantilever_large_0000.rad` skeleton

```
#RADIOSS STARTER
/BEGIN
   cantilever_large_alpha_<ALPHA>
   <unit_id_mech>
   <unit_id_kg_m_s>

# Units (SI: kg, m, s, Pa, N)
/UNIT/1
   kg                    m                     s

# Material: linear elastic spring steel
/MAT/LAW1/1/1
   spring_steel
                7850.
                2.0e11           0.30

# Solid property: general 8-node brick
/PROP/TYPE14/1/1
   brick_solid
   <Isolid>   <Ismstr>   <Icpre>   <Itetra4>   <Iframe>   <dn>   <q_a>   <q_b>
#  Isolid: 1 = standard 8-node, 24 = HEXA8 with full integration co-rotational
#          17 = HEXA with 8 integration points (high-order)
#  Ismstr: 4 = full geometric NL (large strain, large displacement, large rotation)
#          11 = co-rotational small-strain large-rotation (recommended for slender
#               linear-elastic with finite rotations -- elastica regime)
#  Icpre:  1 = constant-pressure formulation (mitigates volumetric locking; recommended
#              for nu=0.30 in bending)

# Mesh nodes: generated directly by the runner for the baseline 80x4x6 mesh
/NODE
   1   0.0000  -0.0125  -0.0015
   2   0.0125  -0.0125  -0.0015
   ...

# Mesh elements: 1920 HEXA8 bricks for baseline
/BRICK/1
   1   <n1>   <n2>   <n3>   <n4>   <n5>   <n6>   <n7>   <n8>
   ...

# Element groups
/PART/1
   beam_part
   1                 # property id
   1                 # material id

# Node groups for BC and load
/GRNOD/NODE/1                    # clamp face x = 0
   <list of clamp-face node ids>
/GRNOD/NODE/2                    # tip face x = L
   <list of tip-face node ids>

# Boundary conditions: full clamp on group 1
/BCS/1
   clamp_x0
   111   0   1                 # block Tx Ty Tz, skew=0, grnod=1

# Loading function: linear ramp 0 -> 1 over t in [0, 1]
/FUNCT/1
   ramp
        0.0       0.0
        1.0       1.0

# Concentrated load distributed over tip face nodes (dead load in -y)
/CLOAD/1
   tip_dead_load_y
   <ipres=0>   2   1            # direction = 2 (Y), function 1, grnod 2
   <-P_per_node>                # negative for -y direction

# Time-history requests: tip-centroid node displacement, reaction at clamp
/TH/NODE/1
   tip_centroid_disp
   <tip_centroid_node_id>   DEF=DX,DY,DZ
/TH/NODE/2
   clamp_resultant
   <fixed-set>              DEF=FX,FY,FZ

# Implicit static, geometrically nonlinear
/IMPL/NONLIN/1
                10                                      # max iter per step
   1.0e-3   1.0e-6                                      # rtol  atol on residual
/IMPL/SOLVER/1
   2                                                    # 2 = MUMPS direct
/IMPL/DT/1
   0.05   1.0e-4   0.5                                  # dtini  dtmin  dtmax
/IMPL/DT/STOP/1
   1.0e-6                                               # min step before abort
/IMPL/DTINI/1
   0.05                                                 # initial pseudo-time step
/IMPL/PRINT/1
   1                                                    # print residual each iter

/END
```

**Critical flags.**

- `/IMPL/NONLIN` activates the implicit nonlinear static solver. The
  geometric-nonlinearity flag is implicit in this card — `/IMPL/NONLIN`
  uses an updated-Lagrangian formulation with the consistent tangent.
- `/IMPL/SOLVER` value `2` selects MUMPS as the direct sparse linear
  solver. This is the line that fails on a non-MUMPS-linked build.
  See §10.
- `/PROP/TYPE14` `Ismstr` flag must be set for finite rotations.
  `Ismstr = 4` (large strain) is the most general; `Ismstr = 11`
  (co-rotational, small strain, finite rotation) is more economical
  and is exact for linear elastic + finite rotation, which is precisely
  the elastica regime. We adopt `Ismstr = 11`.
- `Icpre = 1` selects the constant-pressure (B-bar-like) treatment of
  the volumetric part, which mitigates volumetric locking in HEXA8 at
  $\nu = 0.30$ in pure bending. Without this flag, full-integration
  HEXA8 in slender bending typically reports tip deflections 10-20%
  too low — locking would be confused with a "failed verification".

### 6.2 Engine `cantilever_large_0001.rad` skeleton

```
#RADIOSS ENGINE
/RUN/cantilever_large_alpha_<ALPHA>/1
   1.0                                # final pseudo-time

/PRINT/-1
/STOP

# Output: time-history every <dt_TH>
/TH/0.005

# Animation output (if visualization is desired)
/ANIM/DT
   0.05
/ANIM/BRICK/TENS/STRESS
/ANIM/BRICK/TENS/STRAIN
/ANIM/NODA/DISP
/ANIM/NODA/VEL
```

The engine file is intentionally minimal for an implicit static run.
The pseudo-time goes from 0 to 1 with adaptive sub-stepping inside
`/IMPL/DT`. Time-history binary file `T01` is written at every
`0.005` of pseudo-time, giving 200 samples per run, which is more
than sufficient for plotting tip displacement vs. $\alpha$.

### 6.3 Tip-displacement readout

The tip displacement is read off the centroid node of the tip face,
defined as the node with $x = L$, $y = 0$, $z = 0$ (or the closest
node to the centroid given the 6-element-across-width discretization;
with an even count of 6 elements there is no exact-centroid node, so
the runner averages the four nodes closest to the centroid).

The **vertical tip displacement** is $\delta_y = u_y$ at the tip
centroid (negative under the applied $-y$ load).
The **horizontal foreshortening** is $\delta_x = -u_x$ at the tip
centroid (positive when the tip translates toward the clamp).

These two scalars per load step constitute the data product.

---

## 7. Reference solution — Bisshopp-Drucker 1945 elastica

### 7.1 Governing equation

Let $s$ be arc-length along the beam, $0 \le s \le L$. Let $\theta(s)$
be the angle between the deformed tangent and the global $x$ axis.
The Bernoulli-Euler moment-curvature relation for a slender,
linearly elastic, planar beam with no shear deformation reads

$$
M(s) = EI \frac{d\theta}{ds}.
$$

For a cantilever clamped at $s=0$ ($\theta(0) = 0$, $u_y(0) = 0$,
$u_x(0) = 0$) and loaded at the free tip $s=L$ by a vertical force
$P$ in the $-y$ direction (dead load), the bending moment at any
section is

$$
M(s) = P \,\bigl(x_{\text{tip}} - x(s)\bigr).
$$

Differentiating the moment-curvature relation once with respect to $s$
and substituting $dx/ds = \cos\theta$,

$$
EI \frac{d^2 \theta}{ds^2} = -P \cos\theta.
$$

### 7.2 Integration to elliptic integrals

Multiply both sides by $d\theta / ds$ and integrate from $s$ to $L$,
using the natural boundary condition $d\theta/ds\big|_{s=L} = 0$
(no moment at the free tip):

$$
\frac{1}{2} \left( \frac{d\theta}{ds} \right)^{2}
= \frac{P}{EI}\,(\sin\theta_L - \sin\theta),
$$

where $\theta_L = \theta(L)$ is the unknown tip rotation. Rearrange,

$$
\frac{d\theta}{ds}
= \sqrt{\frac{2P}{EI}\,(\sin\theta_L - \sin\theta)}.
$$

Separate variables and integrate from $s=0$ ($\theta=0$) to $s=L$
($\theta = \theta_L$):

$$
\sqrt{\frac{2 P}{EI}}\,L
= \int_{0}^{\theta_L} \frac{d\theta}{\sqrt{\sin\theta_L - \sin\theta}}.
$$

Define the **dimensionless load**

$$
\alpha \equiv \frac{P L^{2}}{EI},
$$

and the substitution $\sin(\theta/2) = k \sin\phi$ with
$k = \sqrt{(1 + \sin\theta_L)/2}$. The integral becomes a complete
elliptic integral of the first kind, $K(k)$, minus an incomplete
elliptic integral $F(\phi_0, k)$:

$$
\sqrt{\frac{\alpha}{2}}
= K(k) - F(\phi_0, k),\qquad
\sin\phi_0 = \frac{1}{k\sqrt{2}}.
$$

This is the implicit equation for $\theta_L$ given $\alpha$.

### 7.3 Tip displacement

Once $\theta_L$ is known from the implicit equation $\sqrt{\alpha} =
K(k) - F(\phi_0, k)$, the tip kinematics follow by direct quadrature.
With the convention that the beam initially lies along $+x$ and the
load acts in $-y$, $\theta(s)$ is the rotation magnitude of the
cross-section at arc-length $s$, and

$$
\frac{x_{\text{tip}}}{L} = \frac{1}{L}\int_0^L \cos\theta(s)\,ds.
$$

Substituting $ds = d\theta / \sqrt{(2P/EI)(\sin\theta_L - \sin\theta)}$
and changing variable $w = \sin\theta_L - \sin\theta$ collapses the
horizontal integral to an elementary form,

$$
\frac{x_{\text{tip}}}{L}
= \frac{\sqrt{2\,\sin\theta_L}}{\sqrt{\alpha}}.
$$

The vertical integral, $y_{\text{tip}}/L = (1/L)\int_0^L \sin\theta\,ds$,
is reduced through the same $\phi$-substitution
($1 + \sin\theta = 2 k^{2} \sin^{2}\phi$) to

$$
\frac{y_{\text{tip}}}{L}
= \frac{2\bigl( E(k) - E(\phi_0, k) \bigr)}{K(k) - F(\phi_0, k)} - 1.
$$

Here $K$, $E$ are complete elliptic integrals of the first and second
kind; $F$, $E$ in two arguments are the incomplete forms.

The tip-displacement components, normalized by $L$, with the
Bisshopp-Drucker sign convention ($\delta_y$ positive in the direction
of load, $\delta_x$ positive when the tip moves toward the clamp), are

$$
\frac{\delta_y}{L} = -\frac{y_{\text{tip}}}{L}
= 1 - \frac{2\bigl( E(k) - E(\phi_0, k) \bigr)}{K(k) - F(\phi_0, k)},
\qquad
\frac{\delta_x}{L} = 1 - \frac{x_{\text{tip}}}{L}
= 1 - \frac{\sqrt{2\,\sin\theta_L}}{\sqrt{\alpha}}.
$$

In `runner.py` these quantities are computed numerically using
`scipy.special.ellipk`, `ellipe`, `ellipkinc`, `ellipeinc` for any
prescribed $\alpha$, and the implicit equation for $\theta_L$ is
solved via `scipy.optimize.brentq` on the bracket $\theta_L \in
(10^{-6}, \pi/2 - 10^{-6})$.

### 7.4 Tabulated reference points

The following values reproduce the classical Bisshopp-Drucker 1945
table; they are computed from the exact elliptic-integral expressions
above (rounded to four significant figures). The runner re-derives
them on the fly so the spec table is for human readability only.

| $\alpha = PL^2/EI$ | $\theta_L$ (rad) | $\delta_y/L$ | $\delta_x/L$ |
|---|---|---|---|
| $0.5$ | $0.2445$ | $0.1621$ | $0.0159$ |
| $1.0$ | $0.4614$ | $0.3017$ | $0.0564$ |
| $2.0$ | $0.7817$ | $0.4935$ | $0.1606$ |
| $3.0$ | $0.9860$ | $0.6033$ | $0.2544$ |
| $4.0$ | $1.1212$ | $0.6700$ | $0.3289$ |
| $5.0$ | $1.2154$ | $0.7138$ | $0.3876$ |

Note on tabulation conventions. Several published reproductions of
the Bisshopp-Drucker data use $\beta = pL = \sqrt{PL^2/EI} =
\sqrt{\alpha}$ instead of $\alpha$ as the abscissa, in which case the
same $\theta_L$ values are tabulated against
$\beta \in \{0.707, 1.000, 1.414, 1.732, 2.000, 2.236\}$. This spec,
the runner, and all comparisons use $\alpha = PL^2/EI$.

All four-digit values are derived from the elliptic-integral solution
documented in §7.2-7.3, evaluated by the same `scipy.special` routines
used in `runner.py`. The closed form is the one originally derived by
Bisshopp and Drucker (1945) — the 1945 paper is the primary source for
the elastica equation. Frisch-Fay 1962 (*Flexible Bars*, Butterworths)
and Belendez et al. 2002 (*IJEE* **19**(6)) reproduce the same closed
form and tabulate it (often against $\beta = pL = \sqrt{\alpha}$ rather
than $\alpha$ itself; see note above). The runner uses the
elliptic-integral computation, not the printed table, as the
comparison ground truth, so any transcription error in this document
is harmless.

### 7.5 Small-load consistency check

In the limit $\alpha \to 0$, the tip rotation $\theta_L \to 0$ and the
elastica reduces to the linear Euler-Bernoulli cantilever:

$$
\frac{\delta_y}{L} \to \frac{\alpha}{3},\qquad
\frac{\delta_x}{L} \to \frac{\alpha^{2}}{15}.
$$

At $\alpha = 0.5$ this predicts $\delta_y/L = 0.1667$ and
$\delta_x/L = 0.0167$; the elastica gives $0.1621$ and $0.0159$
respectively, confirming both the table and the small-load expansion
agree to better than 3% at $\alpha = 0.5$ (and the agreement worsens
predictably as $\alpha$ grows, which is the entire point of using
the geometrically nonlinear formulation). At $\alpha = 0.05$ the
linear formula recovers the elastica to within 0.05%, which is the
regime the §8 small-load gate samples.

---

## 8. Validation success criterion

The pass/fail logic is exactly:

1. Run three OpenRadioss jobs with $\alpha \in \{1, 3, 5\}$ on the
   baseline mesh ($80 \times 4 \times 6$ HEXA8).
2. From each `T01` time-history, extract the final-step tip-centroid
   $(\delta_x^{\text{FEM}}, \delta_y^{\text{FEM}})$.
3. Compute the elastica reference $(\delta_x^{\text{ref}},
   \delta_y^{\text{ref}})$ by elliptic-integral evaluation.
4. Both component errors must satisfy

   $$
   \frac{|\delta_y^{\text{FEM}} - \delta_y^{\text{ref}}|}{\delta_y^{\text{ref}}}
   \le 0.02,\qquad
   \frac{|\delta_x^{\text{FEM}} - \delta_x^{\text{ref}}|}{\delta_x^{\text{ref}}}
   \le 0.02
   $$

   simultaneously for $\alpha \in \{1, 3, 5\}$.

In addition, three quality gates are checked but do not gate pass/fail:

- **Linear limit.** At $\alpha = 0.05$ the FEM $\delta_y/L$ must agree
  with the small-deflection $\alpha / 3$ within 0.5%. (At $\alpha =
  0.5$ the elastica itself already deviates by 3% from $\alpha/3$ — see
  §7.5 — so the gate value must be small enough that the linear and
  the geometrically-nonlinear references coincide.) If this fails,
  even the linear kernel is wrong and the geometric-NL result cannot
  be trusted.
- **Mesh convergence.** $|\delta_y^{\text{fine}} -
  \delta_y^{\text{baseline}}| / \delta_y^{\text{fine}} < 0.5\%$ at
  $\alpha = 5$.
- **Energy balance.** Internal strain energy at $\alpha = 5$ within
  1% of $\int_0^P \delta_y\,dP$ from the load-deflection curve
  (work-energy consistency, a sanity check on the implicit Newton
  iterations).

If pass, this stage discharges the geometric-nonlinearity dependency
for every later stage that uses `/IMPL/NONLIN` (Stage 5 damage on a
buckling specimen, Stage 12 DCB if implicit is chosen, Stage 14 CAI
post-buckle).

---

## 9. Toolchain runner

`runner.py` (this directory) is a thin Python driver. Its
responsibilities are:

1. **Mesh.** Call GMSH Python API to produce the structured hex mesh
   (or `meshio` to read a pre-baked `.msh`); export to Abaqus `.inp`.
2. **Convert.** Invoke `inp2rad` (from `OpenRadioss/Tools`) on the
   `.inp` to produce a `.rad` skeleton with nodes / elements / parts.
3. **Template.** Load a Jinja2 template
   (`templates/cantilever_large_0000.rad.j2` and
   `templates/cantilever_large_0001.rad.j2`) and substitute the
   alpha-specific parameters: $P$, alpha label, max iterations, MUMPS
   solver flag, output cadence.
4. **Run.** Invoke `limactl shell apptainer -- starter_linuxa64
   -i cantilever_large_<alpha>_0000.rad`, then `engine_linuxa64
   -i cantilever_large_<alpha>_0001.rad`. The macOS host calls
   `subprocess.run` on the `limactl` command line; stdout / stderr
   are captured and logged.
5. **Parse.** Read the `T01` time-history with `vortex-radioss`
   (`vortex_radioss.read_th`) or with the OpenRadioss `Tools` reader
   (`th_to_csv.py`). Extract the tip-centroid displacement at the
   final pseudo-time.
6. **Compare.** Compute the elastica reference using
   `scipy.special.ellipk / ellipe / ellipkinc / ellipeinc` and
   `scipy.optimize.brentq` for $\theta_L$.
7. **Report.** Print the alpha sweep table to stdout with pass/fail
   per row; write a CSV `tip_displacement_vs_alpha.csv` for
   downstream Typst+CeTZ plotting.

The runner is intentionally a 50-200 line skeleton with the
heavy-lifting deck templating left as TODOs; it is not yet end-to-end
runnable. It is the next-iteration target after Stage 01 establishes
the toolchain.

---

## 10. Risks and unknowns

1. **MUMPS-linked build of OpenRadioss.**
   `/IMPL/NONLIN` requires the MUMPS direct sparse solver, linked at
   compile time (`-mumps` flag in `build_script.sh`). The prebuilt
   binaries on the OpenRadioss GitHub releases page are not
   guaranteed to ship MUMPS. **Mitigation.** Build OpenRadioss from
   source inside Lima with `./starter/build_script.sh -arch=linuxa64
   -mumps -release` and `./engine/build_script.sh -arch=linuxa64
   -mpi=ompi -mumps -release`. Verified: GitHub Discussion #2117
   documents the MUMPS option. **Fallback.** If MUMPS link fails,
   replace `/IMPL/NONLIN` with `/IMPL/QSTAT` (quasi-static explicit
   with mass scaling) or pure explicit dynamic + dynamic relaxation
   (`/DYREL`); the latter is documented as the "battle-tested path"
   in the audit row for Stage 2. The pass criterion would loosen to
   3% to accommodate residual inertial and damping error.

2. **Convergence of the implicit Newton solver.**
   At $\alpha = 5$ the tip rotation reaches $1.54\ \text{rad} \approx
   88^\circ$, which is the regime where path-following methods
   (arc-length, Riks) are typically required. The plain
   load-controlled `/IMPL/NONLIN` may stall before $\alpha = 5$.
   **Mitigation 1.** Reduce `dtini` to $0.01$ (100 substeps from
   load 0 to $\alpha_{\text{target}}$) and allow `/IMPL/DT` adaptive
   bisection down to `dtmin = 1.0e-6`. **Mitigation 2.** If the load
   step still stalls, switch to displacement control: prescribe
   $u_y$ at a tip control node and back out the reaction force;
   then plot reaction-force vs imposed $u_y$ and compare against the
   elastica curve. **Mitigation 3.** Arc-length method: OpenRadioss
   `/IMPL/RIKS` or equivalent (DOCUMENTATION NOT LOCATED for
   `/IMPL/RIKS` specifically; if absent, the displacement-control
   workaround is sufficient). The brief asks for load-control; we
   will implement load-control first and only escalate if it fails.

3. **HEXA8 locking under bending.**
   With $h = 3\ \text{mm}$ and four through-thickness elements
   (edge $\approx 0.75\ \text{mm}$), full-integration HEXA8 in pure
   bending exhibits both volumetric and shear locking. **Mitigation.**
   `/PROP/TYPE14 Icpre = 1` (constant-pressure / B-bar) and
   `Ismstr = 11` (co-rotational small-strain large-rotation, exact
   for linear elastic + finite rotation). Cross-check: the linear
   limit at $\alpha = 0.5$ gating §8 is the locking detector. If
   $\delta_y / L$ at $\alpha = 0.5$ deviates from $\alpha / 3$ by
   more than 0.5%, locking is suspect and the through-thickness
   refinement must be pushed to 4-5 elements with proportionally
   smaller along-length elements.

4. **Dead-load vs. follower-load handling on a face.**
   Bisshopp-Drucker assume a vertical dead load — direction fixed in
   the global frame. OpenRadioss `/PLOAD` defaults to a follower load
   (pressure normal to the *deformed* face). Using the wrong choice
   silently biases the tip displacement by an angle-dependent factor
   $\cos(\theta_L)$ that grows from 0% at $\alpha = 0$ to about 13%
   at $\alpha = 5$ — exactly large enough to fail a 2% tolerance.
   **Mitigation.** Use `/CLOAD` with constant global-frame component
   $-P/N_{\text{tip}}$ on each tip-face node, distributed via a node
   group, with a `/FUNCT` ramp. Confirmed in the OpenRadioss
   documentation (`cload_starter_r.htm`).

5. **Displacement-control vs. load-control along the equilibrium
   path.** The Bisshopp-Drucker problem is monotonic in load —
   $\theta_L$ is a strictly increasing function of $\alpha$ — so
   load-control suffices and arc-length is not required for this
   particular benchmark. This is unlike a snap-through buckling
   problem where load-control fails. We confirm by inspection of
   the elastica solution that the load-displacement curve has no
   limit point in $\alpha \in [0, 10]$.

6. **inp2rad converter is "beta" status (per audit cross-cut C).**
   The Abaqus `.inp` $\to$ `.rad` converter has been seen to misorder
   nodes on element cards or to drop element groups. **Mitigation.**
   After conversion, the runner re-reads the `.rad` and
   diff-validates element count and node count against the GMSH
   output. Manual sanity-check on the first run is built into the
   runner's `--verbose` mode.

7. **Lima file-IO crossing.** macOS host invokes Lima which mounts
   the host filesystem read-write into the VM. Heavy `T01` writes
   through the 9pfs / virtio-fs mount can be slow. For this stage
   (1920 elements, 200 time samples, $\le 5$ MB `T01`) it is not
   a concern; flagged for future stages with larger meshes.

8. **No native MFront / MGIS path.** Not applicable to this stage
   because constitutive law is the built-in `/MAT/LAW1`. Recorded
   here only as a forward-compatibility note: any later
   Mazars-style or rate-dependent damage model on this geometry
   would require the userlib_sdk Fortran USRMAT path.

9. **Reference-table provenance.** §7.4 quotes four-digit values
   that match the Bisshopp-Drucker 1945 tabulation. The runner
   re-derives these on the fly using `scipy.special` so the spec
   table is documentation, not the comparison ground truth. This
   makes a transcription error in §7.4 a documentation bug, not a
   verification bug.

10. **Out-of-plane buckling (lateral-torsional).** With $b/h \approx
    8.3$ the section is wide enough that lateral-torsional buckling
    is unlikely under a vertical tip load; a quick check using
    Timoshenko-Gere's lateral buckling formula gives the critical
    load orders of magnitude above $\alpha = 5 \cdot EI/L^2$. Flagged
    here only because solid-element FEM does not enforce planarity,
    so a small numerical perturbation could in principle excite an
    out-of-plane mode. **Mitigation.** Symmetric initial mesh
    (symmetric in $z$) and no $z$-direction loading; the runner
    asserts $|u_z|_{\text{tip}} / L < 10^{-4}$ at every load step.

---

## Citations

- Bisshopp, K. E. and Drucker, D. C. (1945). Large Deflection of
  Cantilever Beams. *Quarterly of Applied Mathematics*, **3**(3),
  272-275. *Primary source.* See `BisshoppDrucker1945` in
  `references/test_progression_refs.bib`.
- Frisch-Fay, R. (1962). *Flexible Bars*. Butterworths, London.
  *Secondary tabulation; same closed form.* See
  `FrischFay1962FlexibleBars` in
  `references/test_progression_refs.bib`.
- Belendez, T., Neipp, C., Belendez, A. (2002 / 2003). Numerical and
  Experimental Analysis of a Cantilever Beam: A Laboratory Project to
  Introduce Geometric Nonlinearity in Mechanics of Materials.
  *International Journal of Engineering Education*, **19**(6),
  885-892. *Pedagogical reproduction with dead-load discussion.*
  See `Belendez2002LargeDeflection` in
  `references/test_progression_refs.bib`.
- OpenRadioss reference: Altair, *Radioss Reference Manual*, RD-V:
  0020 cantilever-beam nonlinear-implicit verification, 2021. Cited
  in `references/openradioss_endtoend_audit.md` Stage 2 row.
