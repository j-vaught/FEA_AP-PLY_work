# Stage 10 - Unidirectional Tow-wise Direct Mesoscale (Reframed)

Author. J.C. Vaught
Date. 2026-04-29
Solver. OpenRadioss (single-tool plan, see master_plan.md §1).
Element family. Solid only (HEXA8 preferred, TETRA10 fallback for the curved fiber-matrix interface).
Pass criterion. Effective $E_1$, $E_2$, $G_{12}$ each within 5 percent of Halpin-Tsai prediction at $V_f = 0.60$ with $\zeta=2$ (transverse) and $\zeta=1$ (shear), and inside the Hashin-Shtrikman 1963 bounds.

---

## 1. Reframe statement: direct mesoscale, not periodic homogenization

The original Stage 10 in `references/test_progression_literature.md` calls for periodic boundary conditions on a representative volume element (RVE), reading off the homogenized orthotropic stiffness from corner-DOF macro-strain drivers. The end-to-end OpenRadioss audit (`references/openradioss_endtoend_audit.md`, cross-cut F and stage rows 10 and 11) established that OpenRadioss exposes no general three-DOF linear-equation MPC keyword. The available constraint cards (`/BCS/CYCLIC` for cylindrical symmetry, `/RBE2` for kinematic rigid bodies, `/RBE3` for force interpolation, `/MPC` for three pre-canned joint types, `/RBODY` for rigid bodies) cannot express the generic constraint $u_i^{A} - u_i^{B} = \varepsilon_{ij} L_j$ that periodic homogenization requires. This is the single hardest gap in the audit.

The user-approved workaround is to reframe Stage 10 as a *direct mesoscale simulation* under uniform displacement boundary conditions, on a specimen that is significantly larger than the fiber diameter so that the volume-averaged response converges to the effective property. This is the kinematic uniform boundary condition (KUBC) of Hill-Mandel, also called *displacement-uniform* or *Dirichlet* homogenization. Under KUBC the apparent stiffness $\bar{C}^{\text{KUBC}}$ is a known *upper bound* on the effective stiffness $C^{\text{eff}}$, with the gap closing as the specimen size grows relative to the heterogeneity length scale. The complementary lower bound is given by static uniform boundary conditions (SUBC, traction-uniform), and periodic BCs land between the two and converge fastest. For a 100 micron specimen containing approximately 156 fibers of diameter 7 micron at $V_f = 0.60$, the KUBC overestimate of $E_2$ and $G_{12}$ has been reported in the literature to fall under 5 percent for similar fiber-matrix contrast (Sun and Vaidya 1996; Kanit, Forest, Galliet, Mounoury, Jeulin 2003 *Int. J. Solids Struct.* on the determination of the size of the representative volume element). This is acceptable for a verification stage whose tolerance is 5 percent.

What this reframe gives up. It does *not* compute the true periodic-homogenized moduli. The reported numbers are KUBC apparent moduli, and they will be biased high relative to true periodic homogenization. The tolerance against Halpin-Tsai (5 percent) is what makes this acceptable - Halpin-Tsai itself is a closed-form engineering approximation with similar error bands.

What it preserves. The verification spirit of the scaffold. The physical question - "does a discrete fiber-matrix solid mesh recover the right effective UD-ply stiffness?" - is answered. Stage 7 (UD tow tension on-axis and off-axis) has already verified the orthotropic stiffness card on a homogenized solid; Stage 10 verifies that the mesoscale architecture *underneath* that homogenized card is consistent. The user's ultimate goal is tow-wise direct simulation of the panel (Stage 16), which never invokes periodic BCs in the first place. Periodic homogenization was always a verification sidecar, not on the critical path.

Citations underpinning the reframe.

- Hill, R. (1963), "Elastic properties of reinforced solids: some theoretical principles", *J. Mech. Phys. Solids*, **11**, 357-372.
- Suquet, P. (1987), "Elements of homogenization for inelastic solid mechanics", in *Homogenization Techniques for Composite Media*, Springer Lecture Notes in Physics 272.
- Kanit, T., Forest, S., Galliet, I., Mounoury, V., Jeulin, D. (2003), "Determination of the size of the representative volume element for random composites", *Int. J. Solids Struct.*, **40**, 3647-3679 - quantitative KUBC vs SUBC vs PBC bias study.
- Hashin, Z. and Shtrikman, S. (1963) - HS bounds, used as the second-order check (`@HashinShtrikman1963`).
- Halpin, J.C. and Tsai, S.W. (1969), AFML-TR-67-423 (`@HalpinTsai1969`).
- Sun, C.T. and Vaidya, R.S. (1996), *Compos. Sci. Technol.*, **56**, 171-179 (`@SunVaidya1996`) - canonical RVE setup that we deviate from.
- Soden, P.D., Hinton, M.J., Kaddour, A.S. (1998), *Compos. Sci. Technol.*, **58**, 1011-1022 (`@SodenHintonKaddour1998`) - constituent properties source.

---

## 2. Geometry and GMSH script outline

### 2.1 Specimen definition

| Quantity | Symbol | Value | Source |
|---|---|---|---|
| Specimen length (fiber direction, $z$) | $L_z$ | 50 microns | reframe choice; $\sim 7 \cdot d_f$ |
| Specimen width (transverse, $x$) | $L_x$ | 100 microns | reframe choice; $\sim 14 \cdot d_f$ |
| Specimen depth (transverse, $y$) | $L_y$ | 100 microns | reframe choice; $\sim 14 \cdot d_f$ |
| Fiber diameter | $d_f$ | 7 microns | IM7-class, Soden 1998 |
| Fiber radius | $r_f$ | 3.5 microns | $d_f / 2$ |
| Volume fraction | $V_f$ | 0.60 | Soden 1998 reference systems |
| Packing | hex | hexagonal | deterministic, repeatable for verification |

The specimen is oriented with the fiber axis along $z$. Choose a right-handed coordinate frame with origin at the lower-corner $(x,y,z)=(0,0,0)$ and the body occupying $[0, L_x] \times [0, L_y] \times [0, L_z]$.

### 2.2 Hex-pack pitch derivation

For a hexagonal close-packed array of parallel fibers in the transverse plane, the relationship between volume fraction $V_f$, fiber diameter $d_f$, and hex pitch $p$ (the inter-fiber spacing along a primitive lattice vector) is

$$
V_f = \frac{\pi d_f^{2}}{2\sqrt{3}\, p^{2}}
\quad \Longrightarrow \quad
p = d_f \sqrt{\frac{\pi}{2\sqrt{3}\, V_f}}.
$$

Substituting $d_f = 7$ microns and $V_f = 0.60$,

$$
p = 7 \cdot \sqrt{\frac{\pi}{2 \cdot 1.7321 \cdot 0.60}} = 7 \cdot 1.2294 = 8.606\;\text{microns},
$$

with row vertical spacing $e_v = p\,\sqrt{3}/2 = 7.453$ microns. The hex array is laid out in $xy$ as rows: row $j$ has fibers at $x_i = (i + 0.5\,(j \bmod 2))\, p$ and $y_j = j\, e_v$. Counting fibers fully or partially inside a $100 \times 100$ micron square cross section gives approximately 156 fibers (also equal to $V_f \cdot A / A_{\text{fiber}} = 0.60 \cdot 10^4 / (\pi \cdot 3.5^2) = 155.9$). Documenting this count is part of the runner output. Periodicity of the lattice is *not* required by the stage because the stage is reframed away from periodic BCs.

### 2.3 GMSH script outline (Python API)

```python
# tests/stage_10_UD_mesoscale_direct/mesh_build.py
import gmsh
import math

L_x, L_y, L_z = 100e-6, 100e-6, 50e-6   # SI metres
d_f = 7e-6
r_f = d_f / 2.0
V_f = 0.60
p = d_f * math.sqrt(math.pi / (2.0 * math.sqrt(3.0) * V_f))   # 8.606 microns
ev = math.sqrt(3.0) / 2.0 * p                                  # 7.453 microns row vertical spacing

gmsh.initialize()
gmsh.model.add("UD_mesoscale_hex")

# 1. Build the rectangular block of matrix.
matrix_box = gmsh.model.occ.addBox(0, 0, 0, L_x, L_y, L_z)

# 2. Build the array of fiber cylinders extruded along z.
fibers = []
n_rows = int(L_y / ev) + 2
n_cols = int(L_x / p) + 2
for j in range(-1, n_rows):
    y_c = j * ev
    if y_c < -r_f or y_c > L_y + r_f:
        continue
    x_offset = 0.5 * p if (j % 2) else 0.0
    for i in range(-1, n_cols):
        x_c = i * p + x_offset
        if x_c < -r_f or x_c > L_x + r_f:
            continue
        cyl = gmsh.model.occ.addCylinder(x_c, y_c, 0, 0, 0, L_z, r_f)
        fibers.append((3, cyl))

# 3. Boolean intersect every fiber with the box (clip protruding parts),
#    boolean fragment to get conformal interfaces between fibers and matrix.
out, _ = gmsh.model.occ.fragment([(3, matrix_box)], fibers)
gmsh.model.occ.synchronize()

# 4. Tag physical groups: matrix (one group), fibers (one group).
all_volumes = gmsh.model.getEntities(dim=3)
fiber_tags, matrix_tags = [], []
for dim, tag in all_volumes:
    com = gmsh.model.occ.getCenterOfMass(dim, tag)
    # crude classifier: centerline of any fiber is within r_f of the COM in xy
    is_fiber = False
    for j in range(-1, n_rows):
        y_c = j * ev
        x_offset = 0.5 * p if (j % 2) else 0.0
        for i in range(-1, n_cols):
            x_c = i * p + x_offset
            if (com[0]-x_c)**2 + (com[1]-y_c)**2 < (0.9*r_f)**2:
                is_fiber = True; break
        if is_fiber: break
    (fiber_tags if is_fiber else matrix_tags).append(tag)
gmsh.model.addPhysicalGroup(3, fiber_tags, name="FIBER")
gmsh.model.addPhysicalGroup(3, matrix_tags, name="MATRIX")

# 5. Tag each of the six outer faces as physical surfaces for BC application.
#    Use bounding-box queries on surface centroids: x=0, x=L_x, y=0, y=L_y, z=0, z=L_z.
# (See runner for face-tag constants FACE_XMIN ... FACE_ZMAX.)

# 6. Mesh size targets: 4 elements across fiber diameter -> ~1.75 micron in fiber,
#    and ~3 micron in bulk matrix far from fibers.
gmsh.option.setNumber("Mesh.MeshSizeMin", 1.5e-6)
gmsh.option.setNumber("Mesh.MeshSizeMax", 3.0e-6)
gmsh.option.setNumber("Mesh.Algorithm3D", 10)   # HXT for tet meshing speed
gmsh.option.setNumber("Mesh.ElementOrder", 2)   # quadratic -> TETRA10
gmsh.model.mesh.generate(3)

# 7. Write Abaqus .inp for the inp2rad path into OpenRadioss.
gmsh.write("UD_mesoscale_hex.inp")
gmsh.finalize()
```

HEXA8 versus TETRA10. A hexahedral mesh with conformal cylindrical interfaces over 156 fibers is impractical to script in GMSH without per-fiber transfinite construction; the typical industrial path for this geometry is a TETRA10 mesh from the OCC kernel followed by `inp2rad`. We default to TETRA10. If the user later needs HEXA8 (for example, to remove the higher-order mass coupling), the GMSH script must be replaced with an extruded-quad cross-sectional mesh (sweep in $z$); this is a separate work item and is documented as a future option, not required for verification of $E_1$, $E_2$, $G_{12}$ at the 5 percent tolerance.

### 2.4 Volume fraction sanity check

After meshing, the runner computes the actual $V_f$ realized by the tessellation as $V_f^{\text{actual}} = \sum_{e \in \text{FIBER}} V_e / V_{\text{total}}$ and asserts $|V_f^{\text{actual}} - 0.60| < 0.005$. A mismatch larger than this indicates that fibers near the boundary have been clipped asymmetrically; in that case, increase $L_x$ and $L_y$ to integer multiples of the hex unit cell.

---

## 3. Mesh and convergence study on element count per fiber

The pass criterion needs at least four elements across the fiber diameter; this is the canonical micromechanics rule of thumb (Sun-Vaidya 1996). The convergence study sweeps three meshes:

| Mesh tag | Elements per fiber diameter | Approximate element edge in fiber | Approximate global node count |
|---|---|---|---|
| coarse | 3 | 2.33 microns | ~80k |
| medium | 4 | 1.75 microns | ~180k |
| fine | 6 | 1.17 microns | ~520k |

Quantity tracked across the sweep. The three effective constants $E_1^{\text{KUBC}}$, $E_2^{\text{KUBC}}$, $G_{12}^{\text{KUBC}}$.

Convergence assertion. The medium-to-fine relative change in each constant is less than 1 percent. If it is not, the medium mesh is rejected and the fine mesh is reported as the verification result.

Mesh objectivity. Linear elasticity has no length scale, so all three meshes should produce the same answer to within numerical noise. The convergence study is a sanity check on (i) the tessellation of the fiber-matrix interface and (ii) the TETRA10 quadrature.

---

## 4. Materials

Both phases are isotropic linear elastic, modeled with `/MAT/LAW1` on `/PROP/SOLID` (TYPE14) generic solid bricks. The carbon fiber is in reality transversely isotropic (axial modulus much greater than transverse), but the verification target Halpin-Tsai uses isotropic-fiber inputs, so an isotropic fiber here is the consistent choice. The limitation is documented; an extension to transversely isotropic IM7 (axial $E_{1f} = 263$ GPa, transverse $E_{2f} = 19$ GPa, $\nu_{12f}=0.20$, $G_{12f}=27.6$ GPa per Soden 1998) is a future work item once `/MAT/LAW12` solid orthotropic on a per-fiber-cylinder local frame is wired through `inp2rad`.

| Phase | Card | $E$ | $\nu$ | $\rho$ | Source |
|---|---|---|---|---|---|
| Fiber (T800/IM7-class carbon) | `/MAT/LAW1` | $230\times 10^{9}$ Pa | $0.20$ | $1780$ kg/m^3 | Soden 1998 Table for IM7 |
| Matrix (Hexcel 8552 epoxy) | `/MAT/LAW1` | $4.08 \times 10^{9}$ Pa | $0.39$ | $1300$ kg/m^3 | Soden 1998 Table for 8552 |

Density is not used in implicit static but is required by the OpenRadioss starter. SI units throughout (Pa, m, kg).

---

## 5. Boundary conditions and loading

Three independent linear static analyses are run, each with a different KUBC macro strain.

### 5.1 Common BC machinery

For a body $\Omega$ with bounding box $[0,L_x] \times [0,L_y] \times [0,L_z]$, the KUBC for a macro-strain tensor $\bar{\boldsymbol{\varepsilon}}$ prescribes on every node $\boldsymbol{x}$ of the boundary $\partial\Omega$,

$$
u_i(\boldsymbol{x}) = \bar{\varepsilon}_{ij}\, x_j \quad \text{for all } \boldsymbol{x} \in \partial\Omega.
$$

In OpenRadioss this is realized through `/BCS` cards on each of the six face node-sets. Per-node displacement values are written into a `/FUNCT` table and applied via an `/IMPDISP` card pointing at that table; alternatively, one `/IMPDISP` per spatial component with a linear ramp from 0 to the prescribed value.

### 5.2 Case A - Axial stretch $\bar{\varepsilon}_{zz} = \varepsilon_0$, all other components zero

Imposed displacements on the boundary (with $\varepsilon_0 = 10^{-3}$, well within the elastic regime for both phases),

$$
u_x = 0, \quad u_y = 0, \quad u_z = \varepsilon_0\, z, \qquad \forall \boldsymbol{x} \in \partial\Omega.
$$

Effective axial modulus extraction,

$$
E_1^{\text{KUBC}} = \frac{\langle \sigma_{zz} \rangle}{\bar{\varepsilon}_{zz}} = \frac{1}{V \varepsilon_0} \int_{\Omega} \sigma_{zz}\, dV.
$$

Effective Poisson contraction (recovered for completeness, not in pass criterion),

$$
\nu_{12}^{\text{KUBC}} = -\frac{\langle \sigma_{xx} \rangle + \langle \sigma_{yy} \rangle}{2\, \langle \sigma_{zz} \rangle}\Big|_{\bar{\varepsilon}_{xx}=\bar{\varepsilon}_{yy}=0}.
$$

(Note: under KUBC with the lateral faces clamped to zero displacement, $\langle \sigma_{xx} \rangle$ and $\langle \sigma_{yy} \rangle$ are reaction stresses, not free Poisson contraction; recovering $\nu_{12}$ from KUBC requires the dual SUBC run and is not part of the pass criterion.)

### 5.3 Case B - Transverse stretch $\bar{\varepsilon}_{xx} = \varepsilon_0$, all other components zero

$$
u_x = \varepsilon_0\, x, \quad u_y = 0, \quad u_z = 0, \qquad \forall \boldsymbol{x} \in \partial\Omega.
$$

$$
E_2^{\text{KUBC}} = \frac{\langle \sigma_{xx} \rangle}{\bar{\varepsilon}_{xx}}.
$$

### 5.4 Case C - Longitudinal shear $\bar{\varepsilon}_{xz} = \gamma_0/2$, all other components zero

The macro-strain tensor for engineering shear $\gamma_{12}$ is $\bar{\varepsilon}_{xz} = \bar{\varepsilon}_{zx} = \gamma_0/2$. The KUBC displacement field is

$$
u_x = \tfrac{\gamma_0}{2}\, z, \quad u_y = 0, \quad u_z = \tfrac{\gamma_0}{2}\, x, \qquad \forall \boldsymbol{x} \in \partial\Omega.
$$

$$
G_{12}^{\text{KUBC}} = \frac{\langle \sigma_{xz} \rangle}{\gamma_0}.
$$

(The factor of 2 disappears because $\langle \sigma_{xz} \rangle = G_{12} \cdot 2\bar{\varepsilon}_{xz} = G_{12}\, \gamma_0$.)

### 5.5 Stress averaging from OpenRadioss output

Output cards used. `/TH/BRICK` for time-history of element stress (overkill for a single linear step, but standardized) and `/ANIM/BRICK/TENS/STRESS` for the full stress tensor field on every element. The runner reads the animation file via the Kitware `openradioss-to-vtkhdf` converter, integrates the stress tensor weighted by element volume, and divides by total volume. For TETRA10 the integration weight is the cell volume read from the `.vtkhdf`.

---

## 6. OpenRadioss deck skeleton

Two files per loading case: `<job>_0000.rad` (starter) and `<job>_0001.rad` (engine). Generated from a Jinja2 template by the runner. Three cases (A, B, C) yield three job pairs.

```
#include "UD_mesoscale_hex.inc"      ; nodes + elements from inp2rad
/MAT/LAW1/1                           ; FIBER
fiber                                  ; title
1780.0                                 ; rho
230.0E9 0.20                           ; E nu
/MAT/LAW1/2                           ; MATRIX
matrix
1300.0
4.08E9  0.39
/PROP/SOLID/1                         ; TYPE14 generic solid for fiber
fiber_solid
0   0   0   0   0   0
/PROP/SOLID/2
matrix_solid
0   0   0   0   0   0
/PART/1
fiber_part
1   1
/PART/2
matrix_part
2   2
/IMPL/LINEAR                          ; static linear analysis
/IMPL/SOLVER/2                        ; MUMPS direct solver
/IMPL/PRINT/N-1                       ; one print per nonlinear iteration (none here)
/BCS/1                                ; case A: u_x=0 on all six face nodesets
1  0  0  111  GRNOD/ALLBOUND          ; fix DOF 1,2,3 = u_x,u_y; (DOF 4,5,6 left as native)
/IMPDISP/1                            ; case A: u_z = eps0 * z, function-driven
3  1  GRNOD/ALLBOUND  FUNCT/1
/FUNCT/1
ramp_uz
0.0   0.0
1.0   1.0
/TH/BRICK/1
all_bricks
GRBRIC/ALL  EPSXX  EPSYY  EPSZZ  EPSXY  EPSYZ  EPSZX  SIGXX  SIGYY  SIGZZ  SIGXY  SIGYZ  SIGZX
/ANIM/DT
0.0  1.0
/ANIM/BRICK/TENS/STRESS
/ANIM/BRICK/TENS/STRAIN
/ANIM/BRICK/VOL
/STOP
```

The above is the skeleton; the actual per-node `u_z = eps0 * z` is realized by *one /IMPDISP per face* with a uniform value (since $z$ is constant on the $z$-min and $z$-max faces) for the four faces orthogonal to $z$, plus an /IMPDISP that sweeps $z$ on the two lateral faces using a node-level variable function. OpenRadioss does not support a closed-form $u(\boldsymbol{x})$ on a face; the runner therefore emits *per-node* /IMPDISP entries for the lateral faces, generated programmatically by reading the face nodeset and computing $u_i = \bar{\varepsilon}_{ij} x_j$ at each node. This is a hundreds-of-nodes table, written by `runner.py` directly into the `.rad` file. The pattern is the same for cases B and C with the appropriate $\bar{\boldsymbol{\varepsilon}}$.

Implicit static is required (`/IMPL/LINEAR`); the build must link MUMPS. If the prebuilt OpenRadioss binaries are explicit-only, fall back to dynamic relaxation (`/DYREL`) with mass scaling, which gives the same linear-static answer at convergence (master plan §10 risk 2).

---

## 7. Reference solutions (Halpin-Tsai and Hashin-Shtrikman)

### 7.1 Constituent-derived numbers (SI)

| Quantity | Formula | Value |
|---|---|---|
| Fiber shear modulus | $G_f = E_f / [2(1+\nu_f)]$ | $95.83$ GPa |
| Matrix shear modulus | $G_m = E_m / [2(1+\nu_m)]$ | $1.468$ GPa |
| Fiber bulk modulus | $K_f = E_f / [3(1-2\nu_f)]$ | $127.78$ GPa |
| Matrix bulk modulus | $K_m = E_m / [3(1-2\nu_m)]$ | $6.182$ GPa |

### 7.2 Halpin-Tsai targets at $V_f = 0.60$

Rule of mixtures for $E_1$ (fiber-dominated, no $\zeta$ needed),
$$
E_1^{\text{ROM}} = V_f E_f + (1-V_f) E_m = 0.60 \cdot 230 + 0.40 \cdot 4.08 = 139.63\;\text{GPa}.
$$

Halpin-Tsai for $E_2$ with $\zeta = 2$,
$$
\eta_{E2} = \frac{E_f/E_m - 1}{E_f/E_m + \zeta} = \frac{56.37 - 1}{56.37 + 2} = 0.9486,
\qquad
E_2^{\text{HT}} = E_m\, \frac{1 + \zeta\, \eta_{E2}\, V_f}{1 - \eta_{E2}\, V_f} = 20.25\;\text{GPa}.
$$

Halpin-Tsai for $G_{12}$ with $\zeta = 1$,
$$
\eta_{G12} = \frac{G_f/G_m - 1}{G_f/G_m + 1} = \frac{65.30 - 1}{65.30 + 1} = 0.9698,
\qquad
G_{12}^{\text{HT}} = G_m\, \frac{1 + \eta_{G12}\, V_f}{1 - \eta_{G12}\, V_f} = 5.55\;\text{GPa}.
$$

### 7.3 Hashin-Shtrikman bounds

The HS bounds for the bulk and shear moduli of a two-phase isotropic mixture (fiber index $f$, matrix index $m$) are, with $V_f + V_m = 1$,

$$
K_{\text{HS}}^{\pm} = K_m + \frac{V_f}{(K_f - K_m)^{-1} + 3 V_m / (3 K_m + 4 G_m)}
\quad \text{(lower, by swap, upper)},
$$

$$
G_{\text{HS}}^{\pm} = G_m + \frac{V_f}{(G_f - G_m)^{-1} + 6 V_m (K_m + 2 G_m) / [5 G_m (3 K_m + 4 G_m)]}
\quad \text{(lower, by swap, upper)}.
$$

The lower bound is computed with $(K_m, G_m)$ as the comparison medium; the upper bound by swapping $f \leftrightarrow m$. These bounds apply to *isotropic* effective moduli; the UD ply is transversely isotropic, so the bounds are used as a *bracketing sanity check* on $E_2$ (transverse) and $G_{12}$ (in-plane shear), not as a sharp pass criterion.

The runner computes both bounds at $V_f = 0.60$ and asserts that the recovered $E_2^{\text{KUBC}}$ and $G_{12}^{\text{KUBC}}$ lie inside the bracket (allowing a small numerical tolerance of 1 percent on the bound values).

### 7.4 Optional Sun-Vaidya 1996 cross-check

Sun and Vaidya 1996 Table 1 reports for AS4/3501-6 ($V_f = 0.60$, $E_f = 235$ GPa, $E_m = 4.8$ GPa, $\nu_f = 0.20$, $\nu_m = 0.34$): $E_1 = 142.6$ GPa, $E_2 = 9.6$ GPa, $G_{12} = 6.0$ GPa, $\nu_{12} = 0.252$. We do not target Sun-Vaidya's $E_2$ specifically because their constituents differ from our 8552/T800-IM7 set; we report Sun-Vaidya as historical context only.

---

## 8. Pass criterion (formalized)

Let $\hat{X}$ denote the value extracted from the FEM run and $X^{\text{HT}}$ the Halpin-Tsai prediction. The stage passes if all four conditions hold.

1. $|\hat{E}_1 - E_1^{\text{ROM}}| / E_1^{\text{ROM}} \le 0.05$.
2. $|\hat{E}_2 - E_2^{\text{HT}}| / E_2^{\text{HT}} \le 0.05$.
3. $|\hat{G}_{12} - G_{12}^{\text{HT}}| / G_{12}^{\text{HT}} \le 0.05$.
4. $\hat{E}_2$ inside the HS bracket and $\hat{G}_{12}$ inside the HS bracket (each with a 1 percent expansion tolerance to absorb numerical noise on the bound).

A failed condition triggers the convergence-study escalation in §3 (refine to the fine mesh and re-test). A persistent failure on the fine mesh is reported as a stage failure with the residuals tabulated.

---

## 9. Honest limitations and caveats

1. KUBC bias. The reported moduli are the upper-bound apparent moduli, not the true effective moduli. Expected bias on $E_2$ and $G_{12}$ at our specimen size is small (order of 1-3 percent based on Kanit et al. 2003 for similar contrast ratios), and the 5 percent pass criterion is set with this bias in mind.
2. Isotropic-fiber assumption. Real carbon fibers are transversely isotropic. The Halpin-Tsai target also uses isotropic fibers, so the comparison is internally consistent, but the $E_2$ and $G_{12}$ recovered here are *not* the values one would get from a transversely isotropic fiber model.
3. No interface phase. Fiber-matrix interface adhesion is assumed perfect (mesh nodes shared at the interface). No cohesive zone, no interfacial slip. This is the conventional micromechanics assumption for stiffness prediction.
4. No fiber waviness or distribution randomness. Hex packing was chosen for verification repeatability; real composites show random packing with local volume-fraction fluctuations that lower $E_2$ by 5-10 percent. This is a known modeling choice, not a bug.
5. Implicit-build dependence. `/IMPL/LINEAR` requires a MUMPS-linked OpenRadioss build. Master plan §10 risk 2 acknowledges this.
6. inp2rad beta. Per the audit, the Abaqus-to-Radioss converter is in beta; spot-check material assignments and node ordering after every conversion.

---

## 10. Per-stage runner script

The companion file `runner.py` automates the pipeline. Sections 5 and 6 above pin the BCs and deck skeleton; the runner emits the three per-case decks, runs the OpenRadioss starter and engine inside the Lima Apptainer VM, converts the animation files to VTKHDF via the Kitware tool, computes volume-averaged stress and strain, derives the three effective moduli, evaluates the four pass conditions of §8, and writes a CSV row per case for the Typst-CeTZ report.

Inputs.
- `UD_mesoscale_hex.inp` from `mesh_build.py` (§2.3).
- Constituent properties in `materials.json` (§4 numbers, SI).
- `or_paths.json` with the Lima Apptainer command and OpenRadioss binary paths.

Outputs.
- `outputs/effective_moduli.csv` with columns `case, mesh, E_or_G_ref_GPa, E_or_G_FEM_GPa, rel_error, HS_lower_GPa, HS_upper_GPa, in_HS_bracket, pass`.
- `outputs/convergence.csv` with columns `mesh, E1_GPa, E2_GPa, G12_GPa, runtime_s, n_elements`.
- `figures/data/stage_10_summary.csv` for the Typst-CeTZ summary plot.
- A short `outputs/stage_10_log.txt` recording the per-case starter and engine returncodes.

Pass logic.

```
overall_pass = all([
    rel_err(E1)  <= 0.05,
    rel_err(E2)  <= 0.05,
    rel_err(G12) <= 0.05,
    E2_HS_lower * 0.99 <= E2_FEM <= E2_HS_upper * 1.01,
    G12_HS_lower * 0.99 <= G12_FEM <= G12_HS_upper * 1.01,
])
```

The runner is described in detail in `runner.py`; the next file in this directory.

---

## Bibliographic keys (from `references/test_progression_refs.bib`)

`@HalpinTsai1969`, `@HashinShtrikman1963`, `@SodenHintonKaddour1998`, `@SunVaidya1996`, `@XiaZhangEllyin2003`. Master plan §6 owns the file; no new keys are created in this stage.
