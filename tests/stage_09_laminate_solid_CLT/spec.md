# Stage 9 — Solid Laminate, Classical Lamination Theory Verification

**Author.** J.C. Vaught
**Date.** 2026-04-29
**Stage.** 9 of 16 (per `plan/master_plan.md` §3, V/E = V).
**Predecessors.** Stages 1-8 (linear elastic solid, geometric nonlinearity, isotropic dogbone, open-hole Kirsch, isotropic damage, Tsai-Wu/Hashin/Puck side-by-side, UD tow on-axis and off-axis, ply-rotation transformation).
**Successors.** Stage 10 (UD tow-wise direct mesoscale) reuses the per-element orientation machinery proved out here.

---

## 1. Goal and scope

Verify that a stack of solid HEXA8 plies, each with its own per-ply /MAT/LAW25 orthotropic card and per-ply /PROP/TYPE14 fiber orientation, reproduces the Classical Lamination Theory (CLT) in-plane laminate stiffness matrix $\mathbf{A}$ within 2 percent component-wise. This is the moment in the test progression where solid composite stacking is exercised in OpenRadioss for the first time.

The audit (`references/openradioss_endtoend_audit.md` row 9) gives this stage a clean PASS verdict. /MAT/LAW25 on a solid /PROP/TYPE14 brick is the documented OpenRadioss pathway for the "one solid element per ply" workflow. The user's hard requirement (solid elements only, no shell laminate stacks) is served directly here, bypassing the conventional shell stack approach (TYPE10, TYPE11, TYPE17, TYPE51, /PCOMPP) which is rejected by the master plan.

**Explicit deviation from "shell laminate" thinking.** Conventional CLT verification uses laminate shell elements with through-thickness Gauss integration (one integration point per ply on a single element layer). Here every ply is its own solid layer, so the through-thickness count is at least the number of plies (4 for cross-ply, 8 for quasi-iso) and is doubled or more after refinement. The cost per element is higher by one-to-two orders of magnitude than a shell stack, but this is the explicit user constraint. We document the per-element through-thickness count and the convergence as the through-thickness mesh is refined.

---

## 2. Geometry

Rectangular laminate plate, in-plane $L_x = L_y = 100$ mm, total thickness $t$ set by ply count and ply thickness. SI units throughout (m, kg, s, Pa).

**Layup A — cross-ply $[0/90]_s$.**
- 4 plies, ply thickness $t_p = 0.18$ mm (Kok `ap_ply_model_creation` default, `references/delft_apply.md`).
- Stacking from bottom (z = 0) to top (z = t): $0^\circ$, $90^\circ$, $90^\circ$, $0^\circ$.
- Total thickness $t = 4 \times 0.18 = 0.72$ mm.

**Layup B — quasi-isotropic $[0/+45/-45/90]_s$.**
- 8 plies, $t_p = 0.18$ mm.
- Stacking from bottom to top: $0$, $+45$, $-45$, $90$, $90$, $-45$, $+45$, $0$ (degrees).
- Total thickness $t = 8 \times 0.18 = 1.44$ mm.

The geometry is meshed in GMSH (Python API) as a structured hex grid. The mesh is partitioned into one named volume per ply so that each volume can be assigned its own /PROP/TYPE14 card with the correct fiber direction.

```
z = t  ┌─────────────────────────┐
       │ ply k = N (top, theta_N)│
       ├─────────────────────────┤
       │ ply k = N-1             │
       ├─────────────────────────┤
       │           ...           │
       ├─────────────────────────┤
       │ ply k = 1 (bottom)      │
z = 0  └─────────────────────────┘
       0 ────── L_x ────── L_x
```

---

## 3. Mesh

Solid HEXA8 brick elements, structured grid.

**Baseline mesh** (one solid element per ply, in-plane $20 \times 20$).
- In-plane element size $h_{xy} = 5$ mm, giving $20 \times 20 = 400$ elements per ply layer.
- Through-thickness: ONE solid element per ply (per the brief's first refinement level). Element thickness equals $t_p = 0.18$ mm. In-plane aspect ratio $h_{xy}/t_p \approx 27.8$ — high but acceptable for in-plane-stretching dominated CLT verification (no bending, $\mathbf{A}$ matrix only).
- Total elements: layup A $= 400 \times 4 = 1600$; layup B $= 400 \times 8 = 3200$.

**Refined mesh** (two solid elements per ply, $20 \times 20$ in-plane).
- Through-thickness: TWO solid elements per ply, each $t_p / 2 = 0.09$ mm thick.
- Total through-thickness elements: layup A $= 8 \ge 4$ (passes brief's "$\ge 4$ elements after first refinement"); layup B $= 16 \ge 4$.
- Total elements: layup A $= 3200$; layup B $= 6400$.

**Optional second refinement** (in-plane $40 \times 40$, two elements per ply).
- Used only if either of the first two meshes fails the 2 percent tolerance.

**Per-ply material orientation.** Each ply layer is a separate volume in the GMSH `.geo`. The OpenRadioss deck assigns one /PROP/TYPE14 + one /MAT/LAW25 per layer, with the orthotropic frame rotated by $\theta_k$ about the $z$-axis. The rotation is implemented via the property card's `Vx`, `Vy`, `Vz` first-direction vector (or via /SKEW), per `help.altair.com/hwsolvers/rad/topics/solvers/rad/prop_type14_solid_starter_r.htm`. `Iorth = 1` so the orthotropic frame co-rotates with the element.

**Convergence study.** The runner computes the FEM-recovered $\mathbf{A}$ matrix on both mesh levels and reports the L2 norm change between them. If the change between baseline (one per ply) and refined (two per ply) is greater than 2 percent in any component, the run is escalated to the second refinement.

---

## 4. Boundary conditions and loading

The laminate is loaded with three independent in-plane uniform-strain cases per layup. Each case applies a prescribed displacement on one pair of opposite edges and constrains rigid-body modes only on the other edges. The bottom face ($z = 0$) and the top face ($z = t$) are traction-free.

**Case 1 — uniaxial extension in $x$, $\bar{\varepsilon}_{xx} = \varepsilon_0$, $\bar{N}_{yy}$ free, $\bar{N}_{xy}$ free.**
- Edge at $x = 0$: $u_x = 0$ on all nodes.
- Edge at $x = L_x$: $u_x = \varepsilon_0 L_x$ on all nodes.
- Edge at $y = 0$ and $y = L_y$: $u_y$ unconstrained (free-edge), $u_x$ unconstrained.
- Two corner nodes pinned in $u_y$ and one corner in $u_z$ to remove rigid-body translation and the in-plane rotation.
- Top and bottom $z$-faces traction-free.

**Case 2 — uniaxial extension in $y$, $\bar{\varepsilon}_{yy} = \varepsilon_0$, $\bar{N}_{xx}$ free, $\bar{N}_{xy}$ free.**
- Edge at $y = 0$: $u_y = 0$.
- Edge at $y = L_y$: $u_y = \varepsilon_0 L_y$.
- Symmetric to case 1.

**Case 3 — in-plane shear, $\bar{\gamma}_{xy} = \gamma_0$.**
- Edge at $y = L_y$: $u_x = \gamma_0 L_y$, $u_y = 0$.
- Edge at $y = 0$: $u_x = 0$, $u_y = 0$.
- Edge at $x = 0$ and $x = L_x$: $u_y = 0$, $u_x$ free (matched displacement to satisfy the simple-shear kinematic).
- Top and bottom $z$-faces traction-free.

The applied strain is small enough to stay linear, $\varepsilon_0 = 10^{-4}$, $\gamma_0 = 10^{-4}$ (no geometric nonlinearity needed; /IMPL/LINEAR is sufficient).

**Why not pure-strain periodic BCs.** OpenRadioss has no general 3-DOF periodic BC keyword (audit cross-cut F, GAP). For the laminate $\mathbf{A}$-matrix recovery this is not a problem because the reference solution is a Kirchhoff thin-plate $\mathbf{A}$-matrix that assumes uniform mid-plane strain. The free-edge displacement BCs above produce uniform mid-plane strain in the interior of the plate (Saint-Venant decay zone is $\sim t$ from each loaded edge, so for $L_x / t \approx 70$ to $140$ the interior is uniform to better than 1 percent). The runner extracts $\mathbf{A}$ from the interior region only (cropping 5 mm from each edge), which removes the free-edge boundary-layer bias.

The biaxial load case mentioned in the brief is optional and would correspond to applying both case 1 and case 2 simultaneously; we implement cases 1, 2, 3 separately because three independent in-plane strain states are exactly what is needed to invert for the full $3 \times 3$ in-plane stiffness $\mathbf{A}$. The biaxial case is therefore redundant for $\mathbf{A}$-matrix recovery and not used.

---

## 5. Material card

**IM7/8552 unidirectional carbon/epoxy.** SI units. Values from Soden, Hinton, and Kaddour (1998), CMH-17 Volume 3 (2012), and Camanho-Maimi-Davila (2007). The Soden 1998 paper itself tabulates the closely related IM7/8551-7 system; the IM7/8552 refinement is in CMH-17 and Camanho et al. 2007. We use the CMH-17 / Camanho 2007 numbers as the operative card and call it "the Soden card for IM7/8552" per the brief.

| Property | Symbol | Value | Unit |
|---|---|---|---|
| Density | $\rho$ | $1580$ | kg/m$^3$ |
| Longitudinal modulus | $E_1$ | $171.42 \times 10^9$ | Pa |
| Transverse modulus | $E_2 = E_3$ | $9.08 \times 10^9$ | Pa |
| Major Poisson | $\nu_{12} = \nu_{13}$ | $0.32$ | — |
| Through-thickness Poisson | $\nu_{23}$ | $0.50$ | — |
| In-plane shear | $G_{12} = G_{13}$ | $5.29 \times 10^9$ | Pa |
| Through-thickness shear | $G_{23}$ | $E_2 / [2 (1 + \nu_{23})] = 3.027 \times 10^9$ | Pa |
| Long. tensile strength | $X_T$ | $2323.5 \times 10^6$ | Pa |
| Long. compressive strength | $X_C$ | $1200.1 \times 10^6$ | Pa |
| Trans. tensile strength | $Y_T$ | $62.3 \times 10^6$ | Pa |
| Trans. compressive strength | $Y_C$ | $199.8 \times 10^6$ | Pa |
| In-plane shear strength | $S_L$ | $92.3 \times 10^6$ | Pa |

The strengths are not exercised in this stage (loads are linear-elastic, $\varepsilon_0 = 10^{-4}$); they are listed because /MAT/LAW25 requires them on the card, and we keep the card consistent with stages 7 and 13. Same card as stage 7.

**OpenRadioss card.** /MAT/LAW25 (CRASURV formulation, Tsai-Wu base) on /PROP/TYPE14 (general orthotropic solid). One material and one property per ply layer, with the rotation built into /PROP/TYPE14 via the orthotropic frame vectors. `Iorth = 1` (frame co-rotates).

---

## 6. OpenRadioss deck skeleton

Two-file two-phase workflow per `references/openradioss_endtoend_audit.md` cross-cut A.

**Starter deck (`stage_09_<layup>_<refinement>_0000.rad`).**

```
#RADIOSS STARTER
/BEGIN
stage_09_<layup>_<refinement>
      2024         0
                 g                 Pa
                 m                 s
/UNIT/1
SI units
kg                  m                  s
/IMPL/LINEAR
/IMPL/PRINT/N-LIN/-1
# --- Nodes (from GMSH, written by mesh writer) ---
/NODE
... (one row per node) ...
# --- Elements (HEXA8) ---
/BRICK/1
... ply 1 elements ...
/BRICK/2
... ply 2 elements ...
... up to ply N ...
# --- Per-ply property + material (one block per ply, with theta_k baked into Vx,Vy,Vz) ---
/PROP/TYPE14/1
ply_1_theta_0deg
   Ihbe Ismstr Itetra4 Iframe   Iorth   Ip Icpre  ...
       Vx           Vy           Vz   skew_id   Iorth=1
/MAT/LAW25/1
ply_1_IM7_8552
       rho
       E_1     E_2     nu_12  ...
       G_12    G_23    G_13
       X_T     X_C     Y_T    Y_C    S_L
       ... (full LAW25 card with strengths zero or large to suppress damage)
# --- Repeat /PROP and /MAT for each ply with rotated (Vx,Vy,Vz) per theta_k ---
# --- Boundary conditions for case <c> in {1,2,3} ---
/BCS/1
... fixed-displacement BCs on the loaded edges ...
/IMPDISP/1
... prescribed displacement amplitude time-function ...
/FUNCT/1
linear ramp 0 -> eps0*L over t in [0, 1]
/END
```

**Engine deck (`stage_09_<layup>_<refinement>_0001.rad`).**

```
#RADIOSS ENGINE
/RUN/stage_09_<layup>_<refinement>/1
1.0
/IMPL/LINEAR
/PRINT/-1
/TFILE
0.05
/ANIM/DT
0.0  1.0
/ANIM/BRICK/STRAIN
/ANIM/BRICK/STRESS
/ANIM/BRICK/TENS/STRESS/ALL
/ANIM/BRICK/TENS/STRAIN/ALL
/H3D/DT
0.0  1.0
/STOP
```

**Per-ply orientation vectors.** For ply $k$ at angle $\theta_k$, the in-plane $V$ vector is $(\cos\theta_k, \sin\theta_k, 0)$ and the through-thickness reference is the global $z$-axis. The deck templater sets `Vx = cos(theta_k)`, `Vy = sin(theta_k)`, `Vz = 0` on /PROP/TYPE14 and lets OpenRadioss build the second in-plane axis automatically (per the documented Iframe = 0 default).

**Total decks generated.** For each (layup, refinement, case) triple = 2 layups $\times$ 2 refinement levels $\times$ 3 cases = 12 starter+engine pairs (plus optional second refinement, if needed: 6 more = 18 total). Each pair is its own subdirectory `runs/<layup>_<refinement>_case<c>/`.

---

## 7. Reference solution — analytic CLT $\mathbf{A}$-matrix

Closed-form, computed in Python and saved as CSV at run time. The algorithm follows Jones (1999, *Mechanics of Composite Materials*, 2nd ed., Chapter 4) and is cross-checked against Daniel and Ishai (2006, Chapter 5). Equations are Jones eq. (2.38) for $Q_{ij}$, Jones eq. (2.84-2.85) for $\bar{Q}_{ij}(\theta)$, and Jones eq. (4.24a) for $A_{ij}$.

**Step 1 — reduced stiffness in ply axes.** Assuming plane stress in the ply (the standard CLT assumption, valid for thin plies even when stacked),
$$
Q_{11} = \frac{E_1}{1 - \nu_{12} \nu_{21}}, \quad
Q_{22} = \frac{E_2}{1 - \nu_{12} \nu_{21}}, \quad
Q_{12} = \frac{\nu_{12} E_2}{1 - \nu_{12} \nu_{21}}, \quad
Q_{66} = G_{12},
$$
with $\nu_{21} = \nu_{12} E_2 / E_1$. For the IM7/8552 card, $\nu_{21} = 0.32 \times 9.08 / 171.42 = 0.01695$, $1 - \nu_{12}\nu_{21} = 0.99458$, giving
$Q_{11} = 172.35$ GPa, $Q_{22} = 9.130$ GPa, $Q_{12} = 2.922$ GPa, $Q_{66} = 5.29$ GPa, $Q_{16} = Q_{26} = 0$.

**Step 2 — rotated reduced stiffness.** For ply $k$ at angle $\theta_k$ measured counterclockwise from the laminate $x$-axis to the ply $1$-axis, with $c = \cos\theta_k$, $s = \sin\theta_k$,
$$
\bar{Q}_{11}(\theta) = Q_{11} c^4 + 2 (Q_{12} + 2 Q_{66}) s^2 c^2 + Q_{22} s^4,
$$
$$
\bar{Q}_{22}(\theta) = Q_{11} s^4 + 2 (Q_{12} + 2 Q_{66}) s^2 c^2 + Q_{22} c^4,
$$
$$
\bar{Q}_{12}(\theta) = (Q_{11} + Q_{22} - 4 Q_{66}) s^2 c^2 + Q_{12} (s^4 + c^4),
$$
$$
\bar{Q}_{66}(\theta) = (Q_{11} + Q_{22} - 2 Q_{12} - 2 Q_{66}) s^2 c^2 + Q_{66} (s^4 + c^4),
$$
$$
\bar{Q}_{16}(\theta) = (Q_{11} - Q_{12} - 2 Q_{66}) s c^3 + (Q_{12} - Q_{22} + 2 Q_{66}) s^3 c,
$$
$$
\bar{Q}_{26}(\theta) = (Q_{11} - Q_{12} - 2 Q_{66}) s^3 c + (Q_{12} - Q_{22} + 2 Q_{66}) s c^3.
$$

**Step 3 — laminate $\mathbf{A}$ matrix.** Sum over plies,
$$
A_{ij} = \sum_{k=1}^{N} \bar{Q}_{ij}(\theta_k) \, t_k,
$$
where $t_k = z_k - z_{k-1}$ is the thickness of ply $k$. For all plies of equal thickness $t_p$, $A_{ij} = t_p \sum_k \bar{Q}_{ij}(\theta_k)$.

**Reference $\mathbf{A}$ for layup A, $[0/90]_s$, $t_p = 0.18$ mm.** Symmetry gives $\bar{Q}(0) + \bar{Q}(90) = \bar{Q}(0) + \bar{Q}(0)^T_{1\leftrightarrow 2}$, which produces $A_{11} = A_{22} = 2 t_p (Q_{11} + Q_{22}) = 0.72 \text{ mm} \times (172.35 + 9.13) \text{ GPa} = 130.67 \times 10^6$ N/m, $A_{12} = 4 t_p Q_{12} = 4 \times 0.00018 \text{ m} \times 2.922 \times 10^9 \text{ Pa} = 2.104 \times 10^6$ N/m, $A_{66} = 4 t_p Q_{66} = 3.809 \times 10^6$ N/m, $A_{16} = A_{26} = 0$. Numerical values written to `cltA_layup_A_reference.csv` by `runner.py`.

**Reference $\mathbf{A}$ for layup B, $[0/+45/-45/90]_s$.** Quasi-isotropy gives $A_{11} = A_{22}$, $A_{16} = A_{26} = 0$, $A_{66} = (A_{11} - A_{12}) / 2$. Numerical values written to `cltA_layup_B_reference.csv` by `runner.py`.

**Cross-check.** For layup B, the apparent in-plane modulus $E_x = (A_{11} A_{22} - A_{12}^2) / (t A_{22})$ and apparent Poisson $\nu_{xy} = A_{12}/A_{22}$ are reported alongside the $\mathbf{A}$ entries; quasi-isotropy demands $E_x = E_y$. This cross-check matches the success criterion in `references/test_progression_literature.md` Stage 9 ("$E_x = (A_{11} A_{22} - A_{12}^2)/(t A_{22})$ within 1 percent").

---

## 8. Recovery of FEM $\mathbf{A}$-matrix and validation criterion

**FEM $\mathbf{A}$ recovery.** The laminate $\mathbf{A}$-matrix relates resultants $\mathbf{N} = (N_{xx}, N_{yy}, N_{xy})^T$ to mid-plane strains $\bar{\boldsymbol{\varepsilon}}^0 = (\bar{\varepsilon}_{xx}^0, \bar{\varepsilon}_{yy}^0, \bar{\gamma}_{xy}^0)^T$ via $\mathbf{N} = \mathbf{A} \bar{\boldsymbol{\varepsilon}}^0$. Three independent strain states give three independent columns of $\mathbf{N}$, from which $\mathbf{A}$ is read directly.

For each load case $c \in \{1, 2, 3\}$ the runner does the following.
1. Read the OpenRadioss `.anim` (or VTKHDF) at the final time step.
2. For each in-plane Cartesian element in the cropped interior region ($x \in [5, 95]$ mm, $y \in [5, 95]$ mm), pull the global Cartesian stress tensor components $\sigma_{xx}^e$, $\sigma_{yy}^e$, $\sigma_{xy}^e$ and the through-thickness coordinate of the element centroid.
3. Compute the laminate force resultants by through-thickness integration:
$$
N_{xx}^{(c)}(x, y) = \int_0^t \sigma_{xx}(x, y, z) \, dz \approx \sum_{e \in \text{column}(x,y)} \sigma_{xx}^e \, h_z^e,
$$
where $h_z^e$ is the through-thickness extent of element $e$, and the sum is over the elements stacked at the in-plane location $(x, y)$. Same for $N_{yy}^{(c)}$ and $N_{xy}^{(c)}$.
4. Average over the cropped interior to obtain the uniform $\bar{\mathbf{N}}^{(c)} = \langle \mathbf{N}^{(c)}(x, y)\rangle$, with the spatial standard deviation reported as a sanity check (the BC-induced edge boundary layer should already be cropped out, so the standard deviation should be below 1 percent of the mean for a converged mesh).
5. Read the applied mean strain $\bar{\boldsymbol{\varepsilon}}^{0,(c)}$ directly from the prescribed displacement (case 1: $(\varepsilon_0, 0, 0)$; case 2: $(0, \varepsilon_0, 0)$; case 3: $(0, 0, \gamma_0)$). Verify against the FEM-computed average mid-plane strain (the average in-plane Cauchy strain across all elements in the cropped interior); discrepancy must be below 1 percent.
6. Stack the three load cases column-wise to form $\mathbf{A}^{\text{FEM}} \in \mathbb{R}^{3\times 3}$:
$$
[\bar{\mathbf{N}}^{(1)}, \bar{\mathbf{N}}^{(2)}, \bar{\mathbf{N}}^{(3)}] = \mathbf{A}^{\text{FEM}} [\bar{\boldsymbol{\varepsilon}}^{0,(1)}, \bar{\boldsymbol{\varepsilon}}^{0,(2)}, \bar{\boldsymbol{\varepsilon}}^{0,(3)}].
$$
Inverting the right-hand strain matrix gives $\mathbf{A}^{\text{FEM}} = [\bar{\mathbf{N}}^{(1)}, \bar{\mathbf{N}}^{(2)}, \bar{\mathbf{N}}^{(3)}] \cdot \text{diag}(\varepsilon_0, \varepsilon_0, \gamma_0)^{-1}$.

**Pass criterion (per-component, 2 percent).** For every $(i, j) \in \{(1,1), (2,2), (1,2), (6,6)\}$:
$$
\frac{|A_{ij}^{\text{FEM}} - A_{ij}^{\text{CLT}}|}{|A_{ij}^{\text{CLT}}|} \le 0.02.
$$
The off-diagonal $A_{16}, A_{26}$ for both layups are CLT-zero; the FEM result must be smaller in magnitude than $0.02 \times \max(|A_{11}|, |A_{22}|, |A_{66}|)$ (equivalent two-percent tolerance, normalized to the largest in-plane component since dividing by zero is meaningless).

**Convergence criterion.** Recompute $\mathbf{A}^{\text{FEM}}$ on the refined mesh (two elements per ply); the L2 norm change between baseline and refined is reported. If the change is below 1 percent the result is mesh-converged at the baseline. If the baseline already passes the 2 percent CLT tolerance and the refined-vs-baseline change is below 1 percent, the stage is PASS.

**Cross-check on apparent moduli.** The runner also reports apparent $E_x, E_y, \nu_{xy}, G_{xy}$ from $\mathbf{A}^{\text{FEM}}$ (formulas: $E_x = (A_{11} A_{22} - A_{12}^2)/(t A_{22})$, $\nu_{xy} = A_{12}/A_{22}$, $G_{xy} = A_{66}/t$). These are compared to the textbook closed forms (Jones 1999 Section 4.3, Daniel-Ishai 2006 Chapter 5).

---

## 9. Citations

Primary (closed-form CLT and rotation matrices).
- Jones, R.M. (1999). *Mechanics of Composite Materials*, 2nd ed. Taylor and Francis. Chapters 2 and 4. Bib key `Jones1999MechanicsCompositeMaterials`. **Primary** for the $Q_{ij}$, $\bar{Q}_{ij}(\theta)$, and $A_{ij}$ definitions.

Secondary (cross-checks and worked examples).
- Daniel, I.M. and Ishai, O. (2006). *Engineering Mechanics of Composite Materials*, 2nd ed. Oxford University Press. Chapter 5. Bib key `DanielIshai2006`. Cross-check on $\bar{Q}(\theta)$ and apparent-modulus formulas.
- Reddy, J.N. (2003). *Mechanics of Laminated Composite Plates and Shells: Theory and Analysis*, 2nd ed. CRC Press. Bib key `Reddy2003LaminatedPlatesShells`. Reference for full $\mathbf{ABD}$ machinery; we exercise only $\mathbf{A}$ here.
- CMH-17 Volume 3 (2012). Composite Materials Handbook — Polymer Matrix Composites. Bib key `CMH17V3`. IM7/8552 material properties.

Material card.
- Soden, P.D., Hinton, M.J., Kaddour, A.S. (1998). "Lamina Properties, Lay-up Configurations and Loading Conditions for a Range of Fibre-reinforced Composite Laminates", *Composites Science and Technology* **58**, 1011-1022. Bib key `SodenHintonKaddour1998`. The IM7/8551-7 entry; IM7/8552 numbers we use are the CMH-17 / Camanho 2007 refinement.
- Camanho, P.P., Maimi, P., Davila, C.G. (2007). "Prediction of size effects in notched laminates using continuum damage mechanics", *Composites Science and Technology* **67**(13), 2715-2727. IM7/8552 baseline.

OpenRadioss documentation.
- Altair, /MAT/LAW25 reference. `help.altair.com/hwsolvers/rad/topics/solvers/rad/composite_material_intro_c.htm`.
- Altair, /PROP/TYPE14 reference. `help.altair.com/hwsolvers/rad/topics/solvers/rad/prop_type14_solid_starter_r.htm`.
- Altair, /IMPL/LINEAR. `help.altair.com/hwsolvers/rad/topics/solvers/rad/implicit_analysis_activation_r.htm`.

Test progression context.
- This repo, `references/test_progression_literature.md` Stage 9.
- This repo, `references/openradioss_endtoend_audit.md` row 9 (verdict PASS).
- This repo, `plan/master_plan.md` §3 row 9.

---

## 10. Outputs and artifacts

Outputs land under `tests/stage_09_laminate_solid_CLT/runs/` and `tests/stage_09_laminate_solid_CLT/results/`.

Per (layup, refinement, case) run.
- `runs/<layup>_<refinement>_case<c>/stage_09_*_0000.rad` — starter deck (templated, version-controllable).
- `runs/<layup>_<refinement>_case<c>/stage_09_*_0001.rad` — engine deck.
- `runs/<layup>_<refinement>_case<c>/stage_09_*A001.anim` — final-state animation file.
- `runs/<layup>_<refinement>_case<c>/stage_09_*T01` — time-history file.
- `runs/<layup>_<refinement>_case<c>/stage_09_*.h3d` — HyperWorks 3D output.
- `runs/<layup>_<refinement>_case<c>/stage_09_*.vtkhdf` — converted via Kitware `openradioss-to-vtkhdf`.

Per (layup, refinement) aggregate.
- `results/cltA_<layup>_reference.csv` — analytic CLT $\mathbf{A}$ matrix (3 rows, 3 cols, units N/m).
- `results/cltA_<layup>_<refinement>_FEM.csv` — FEM-recovered $\mathbf{A}$ matrix.
- `results/cltA_<layup>_<refinement>_compare.csv` — entry-wise relative error and pass/fail.
- `results/apparent_moduli_<layup>_<refinement>.csv` — $E_x, E_y, G_{xy}, \nu_{xy}$, FEM and CLT side by side.

Stage-level summary.
- `results/stage_09_summary.csv` — one row per (layup, refinement) with PASS/FAIL flag and per-component max relative error.
- `results/stage_09_summary.md` — narrative summary.

Plots (Typst + CeTZ per global preferences, no matplotlib).
- `figures/figA_layup_A_Aij_compare.typ` and `.pdf` — bar chart of $A_{ij}^{\text{FEM}}$ vs $A_{ij}^{\text{CLT}}$ for layup A; brand colors (Garnet for FEM, Atlantic for CLT).
- `figures/figB_layup_B_Aij_compare.typ` and `.pdf` — same for layup B.
- `figures/figC_convergence.typ` and `.pdf` — relative error vs. through-thickness elements per ply (1 vs 2 vs optional 4); brand color Garnet line, no rounded corners.
- `figures/figD_strain_field.typ` and `.pdf` — interior cropped region strain map for case 1, case 2, case 3 (one panel each).

CSV files are the data source; figures reference the CSVs via Typst's `csv()` import per the global figure-generation rule.

---

## Sanity checks before running

Before launching OpenRadioss the runner must verify the following.
1. Each ply layer has its own /PROP/TYPE14 with the orthotropic-frame $V$ vector matching its $\theta_k$.
2. The CLT reference $\mathbf{A}$ for layup B satisfies quasi-isotropy ($A_{11} = A_{22}$, $A_{16} = A_{26} = 0$, $A_{66} = (A_{11} - A_{12})/2$) to within $10^{-10}$.
3. The CLT reference $\mathbf{A}$ for layup A satisfies cross-ply symmetry ($A_{11} = A_{22}$, $A_{16} = A_{26} = 0$) to within $10^{-10}$.
4. The mesh has the expected through-thickness count: 4 (resp. 8) on baseline; 8 (resp. 16) on refined.
5. The applied displacement boundary conditions on case $c$ produce, at the final time step, a uniform mid-plane strain in the cropped interior; spatial standard deviation < 1 percent of mean.

If any check fails, the stage stops with a non-zero return code and a diagnostic on stdout. Sanity check 1 is the key one — historical bugs in deck templating put the wrong $\theta_k$ on a layer and produce a "FEM passes but for the wrong reason" failure mode. The runner deliberately re-derives $\theta_k$ from the property card $V$ vector and asserts it matches the layup definition.
