# Canonical FEM Verification Test Progression for AP-PLY Composite Ballistic Modeling

Author: J.C. Vaught
Date: 2026-04-29

This document defines a graduated, scientifically rigorous progression of finite element verification tests, building from the simplest linear-elastic beam through to a full pseudo-woven composite (AP-PLY) panel under ballistic projectile impact. Each stage isolates a single capability so that confidence is established in that piece of the FEM machinery before the next physics is introduced.

For every stage we record (i) the ASTM / ISO / MIL standard if applicable, (ii) the canonical reference paper(s) providing closed-form or experimental benchmark data, (iii) representative specimen dimensions, (iv) the failure or success criterion of the physical experiment, (v) any FEniCSx tutorial that already implements the problem, and (vi) the validation success criterion, i.e. the quantitative tolerance that the FEM result must meet against the reference solution to count as "verified".

The progression is deliberately additive. Every later stage assumes that all earlier stages have been verified.

---

## Stage 1 - Linear-Elastic Beam (Euler-Bernoulli vs Timoshenko)

Goal. Verify that the FEM linear-elastic kernel reproduces a closed-form solution for a slender (Euler-Bernoulli) and a moderately thick (Timoshenko) beam in three-point and four-point bending. This isolates the assembly of the stiffness matrix, integration, and Dirichlet/Neumann boundary conditions.

Standard. ASTM E290 (bend testing of metallic materials, ductility); ASTM D790 (flexural properties of unreinforced and reinforced plastics, three-point bending). For composite three-point flexure, ASTM D7264.

Reference. Timoshenko, S.P. and Gere, J.M., *Mechanics of Materials* (1972), provides the canonical closed-form deflection $\delta_{\max} = PL^3/(48 EI)$ for a simply supported beam under midspan load (Euler-Bernoulli). Timoshenko's correction adds the shear-deflection term $PL/(4 \kappa GA)$. For a higher-order benchmark, Reddy's "Refined Theories of Plates and Shells" and the closed-form Levinson/Reddy beam solutions are standard. The Timoshenko-Ehrenfest exact solutions in Elishakoff (2020) tabulate analytic deflections and rotations. For four-point bending of a constant-moment region, see Gere and Goodno (2012, *Mechanics of Materials*).

Specimen dimensions. Slender beam: $L=200$ mm, $b=10$ mm, $h=10$ mm ($L/h=20$, Euler-Bernoulli regime). Thick beam: $L=50$ mm, $b=10$ mm, $h=10$ mm ($L/h=5$, Timoshenko regime). Midspan load $P=100$ N, isotropic steel ($E=210$ GPa, $\nu=0.3$).

Failure / success criterion of the physical experiment. Onset of yield at the bottom-fiber midspan (3-pt) or in the constant-moment region (4-pt). Linear regime only at this stage, no failure modeled.

FEniCSx tutorial. (a) Bleyer, *Computational Mechanics Numerical Tours with FEniCSx*, "Linear elasticity" tour (https://bleyerj.github.io/comet-fenicsx/intro/linear_elasticity/linear_elasticity.html). (b) Dokken, *FEniCSx tutorial* - "The equations of linear elasticity" (https://jsdokken.com/dolfinx-tutorial/chapter2/linearelasticity.html). (c) Gupta, FEniCS Course Day 3 - "Beam bending" (https://abhigupta.io/fenics-workshop/2_notebooks/day-3/tutorials/1_beam_bending/).

Validation success criterion.
- 3-pt midspan deflection: $|\delta_{\text{FEM}} - \delta_{\text{analytic}}|/\delta_{\text{analytic}} \le 1\%$ on a converged mesh.
- 4-pt constant-moment region: bending stress $\sigma_{xx}$ within 1% of $My/I$ at the surface.
- For the thick-beam (Timoshenko) case, compare to Timoshenko-Gere closed form, RMS error on $w(x)$ along the beam axis $\le 2\%$.
- Mesh convergence study must show monotonic $h$-refinement convergence at the theoretical rate ($O(h^2)$ for linear Lagrange, $O(h^3)$ for quadratic).

---

## Stage 2 - Cantilever Beam under Tip Load (Geometric Nonlinearity)

Goal. Introduce geometric nonlinearity. Verify that the total Lagrangian / updated Lagrangian kernel reproduces the elastica solution at large tip rotations, while reducing to Euler-Bernoulli in the small-displacement limit.

Standard. No dedicated standard; benchmark is academic.

Reference. Bisshopp, K.E. and Drucker, D.C., "Large Deflection of Cantilever Beams", *Quarterly of Applied Mathematics*, **3**, 272-275 (1945). This is the canonical closed-form elastica benchmark; a tip-loaded cantilever's deflection is integrated in terms of complete and incomplete elliptic integrals of the first and second kind. Also Frisch-Fay, *Flexible Bars* (1962). Modern tabulations: Belendez et al. (2002, "Numerical and experimental analysis of a cantilever beam: a laboratory project to introduce geometric nonlinearity"), and the SIAM analytical solution paper (https://www.siam.org/media/yhvm4yyz/s135734pdf-copy.pdf).

Specimen dimensions. Slender beam: $L=1000$ mm, rectangular cross section $b=20$ mm, $h=2$ mm; isotropic steel. Tip load swept so that $PL^2/(EI) \in [0, 10]$, taking the maximum tip rotation $\theta_{\max}$ from $0$ to $\sim 80^{\circ}$.

Failure / success criterion of the physical experiment. None - elastic regime, no yield.

FEniCSx tutorial. (a) Dokken, "Hyperelasticity" (https://jsdokken.com/dolfinx-tutorial/chapter2/hyperelasticity.html). (b) NewFrac FEniCSx Training, "Finite Elasticity Part I" - clamped beam sagging under self-weight (https://newfrac.gitlab.io/newfrac-fenicsx-training/02-finite-elasticity/finite-elasticity-I.html).

Validation success criterion.
- Tip vertical deflection $\delta_v / L$ within 2% of Bisshopp-Drucker's tabulated values across the full sweep of dimensionless load $\alpha = PL^2/(EI)$.
- Tip horizontal deflection $\delta_h / L$ within 2% (this checks foreshortening - the small-displacement model gets it identically wrong, confirming the value of geometric nonlinearity).
- In the small-load limit $\alpha \to 0$, FEM result must collapse to the linear $\delta_v = PL^3/(3EI)$ within 0.5%.

---

## Stage 3 - Tensile Coupon, Isotropic (ASTM E8 Dogbone)

Goal. Verify uniaxial stress / strain / engineering stress-strain extraction on a dogbone, including elastic-plastic constitutive integration if a J2 plasticity model is added (optional substage 3b).

Standard. ASTM E8 / E8M (Standard Test Methods for Tension Testing of Metallic Materials).

Reference. ASTM E8/E8M-22 (https://store.astm.org/e0008_e0008m-22.html). For elastic-plastic verification, the canonical FEM benchmark is the "necking of a circular bar" problem of Simo and Hughes, *Computational Inelasticity* (1998), Chapter 3, which has well-established reference solutions.

Specimen dimensions (ASTM E8 standard rectangular sub-size). Reduced section width $w=12.5$ mm, gauge length $L_0=50$ mm, thickness $t=6$ mm, total length $\sim 200$ mm, grip width $25$ mm with a smooth fillet radius of $12.5$ mm. Round specimens: $D_0=12.5$ mm, $L_0=50$ mm.

Failure / success criterion of the physical experiment. Yield strength (0.2% offset), ultimate tensile strength (UTS), uniform elongation, and reduction-of-area at fracture. For elastic verification, no failure - just a constant-strain region in the gauge.

FEniCSx tutorial. (a) Bleyer, "Linear elasticity" (https://bleyerj.github.io/comet-fenicsx/intro/linear_elasticity/linear_elasticity.html). (b) Bleyer, "Nonlinear elasto-plasticity" tour (https://bleyerj.github.io/comet-fenicsx/) for the J2 plasticity substage.

Validation success criterion.
- Engineering stress in the gauge section uniform within $\pm 0.5\%$ of $F/A_0$ (confirms Saint-Venant has decayed past the fillets).
- Apparent Young's modulus extracted from $\sigma$-$\epsilon$ slope within 1% of the input $E$.
- Lateral strain ratio reproduces input Poisson's ratio within 1%.
- For the J2 substage, yield onset, hardening curve, and load-displacement to UTS within 3% RMS of an Abaqus / reference-code solution on the same mesh.

---

## Stage 4 - Tensile Coupon with Stress Concentration (Open-Hole Tension)

Goal. Verify that the FEM resolves a stress concentration with sufficient mesh refinement. The Kirsch closed-form solution for an infinite plate gives $K_t = 3$ for a circular hole under uniaxial tension; for a finite-width plate, Howland's series solution and Heywood's empirical formula give $K_t > 3$.

Standard. ASTM D5766 / D5766M (Open-Hole Tensile Strength of Polymer Matrix Composite Laminates). For metallic open-hole, ASTM E338 (sharp-notch tension) or general practice from Peterson's *Stress Concentration Factors*.

Reference. Kirsch, E.G. (1898), "Die Theorie der Elastizitat und die Bedurfnisse der Festigkeitslehre", *Zeitschrift des Vereines deutscher Ingenieure*, **42**, 797-807, gives $\sigma_{\theta\theta}/\sigma_{\infty} = 3$ at the hole edge for an infinite plate. For finite-width corrections: Howland, R.C.J. (1929-1930), *Phil. Trans. R. Soc. London* A 229, 49-86, and Peterson, R.E., *Stress Concentration Factors* (1974). For composite OHT: Camanho, P.P., Maimi, P., Davila, C.G., "Prediction of size effects in notched laminates using continuum damage mechanics", *Composites Science and Technology*, **67**(13), 2715-2727 (2007); Lavoie, J.A., Soutis, C., Morton, J., "Apparent strength scaling in continuous fiber composite laminates", *Composites Science and Technology*, **60**, 283-299 (2000); Green, B.G., Wisnom, M.R., Hallett, S.R., "An experimental investigation into the tensile strength scaling of notched composites", *Composites Part A*, **38**, 867-878 (2007).

Specimen dimensions. Composite OHT: rectangular plate $L=200-300$ mm, $w=36$ mm, hole diameter $d=6$ mm ($w/d=6$), thickness 2-4 mm. The metallic Kirsch verification uses a thin plate of height/width = 5 (approximately infinite) with a single 6 mm hole.

Failure / success criterion of the physical experiment. Net-section ultimate tensile strength (open-hole tensile strength, OHT) reported in MPa; load divided by gross cross section disregarding the hole, per ASTM D5766. Failure typically initiates at the hole edge at the $\theta=\pm 90^{\circ}$ position.

FEniCSx tutorial. (a) Bleyer, "Isotropic and orthotropic plane stress elasticity" tour (https://bleyerj.github.io/comet-fenicsx/tours/linear_problems/isotropic_orthotropic_elasticity/isotropic_orthotropic_elasticity.html) - this directly uses a perforated plate. (b) Bleyer, "2D linear elasticity" demo (https://comet-fenics.readthedocs.io/en/latest/demo/elasticity/2D_elasticity.py.html).

Validation success criterion.
- Hoop stress $\sigma_{\theta\theta}$ at the hole edge converges to $K_t \cdot \sigma_{\infty}$ within 2%, where $K_t$ comes from Howland's finite-width correction (i.e. $K_t \approx 3.0$ for the infinite-plate limit, with explicit width corrections from Heywood for $w/d=6$).
- Stress decay along the net-section ligament matches Kirsch's polynomial decay $\sigma_{\theta\theta}(r) = (\sigma_{\infty}/2)[2 + (a/r)^2 + 3(a/r)^4]$ within 3%.
- Mesh refinement study at the hole edge with at least four element layers within $0.1 a$.

---

## Stage 5 - Tensile Dogbone with Damage (Mazars / Lemaitre / Isotropic Damage)

Goal. Introduce the softening branch and continuum damage. Verify that an isotropic damage model produces strain softening that is mesh-objective when fracture-energy regularization is used.

Standard. No dedicated test standard; verification is done against published damage benchmarks. ASTM E399 / E1820 cover fracture toughness.

Reference. Mazars, J. (1986), "A description of micro- and macroscale damage of concrete structures", *Engineering Fracture Mechanics*, **25**(5/6), 729-737. Pijaudier-Cabot, G. and Mazars, J., "Damage Models for Concrete", *Handbook of Materials Behavior Models* (2001). For Lemaitre: Lemaitre, J. (1985), "A continuous damage mechanics model for ductile fracture", *Journal of Engineering Materials and Technology*, **107**, 83-89. de Souza Neto, E.A., Peric, D., Owen, D.R.J., *Computational Methods for Plasticity* (2008), Chapters 12-13 give the canonical FEM benchmark of a Lemaitre 1D bar with strain localization, with closed-form softening response. Modern benchmarking with energy regularization: Vassaux, M. et al. (2022), "A modified Mazars damage model with energy regularization", *Engineering Fracture Mechanics*.

Specimen dimensions. Notched dogbone, gauge $L_0 = 50$ mm, reduced width $w=12.5$ mm, with a small geometric imperfection (1% width reduction) at midspan to localize damage.

Failure / success criterion of the physical experiment. Peak load and post-peak softening branch; total dissipated energy $G_f$ matches input fracture energy.

FEniCSx tutorial. (a) Bleyer, "Gradient damage models" tour series (https://bleyerj.github.io/comet-fenicsx/) which includes Mazars-type models with non-local regularization. (b) Through MFront/MGIS-FEniCS integration: Helfer et al., "Integration of MFront constitutive laws into FEniCS - phase-field damage" (https://thelfer.github.io/mgis/web/mgis_fenics_phase_field.html).

Validation success criterion.
- Peak load within 3% of the analytical or reference Mazars / Lemaitre solution (1D bar).
- Post-peak softening: total dissipated energy $\int \sigma \, d\epsilon \cdot V_{\text{loc}} = G_f \cdot A_{\text{frac}}$ within 5%, on at least three meshes (h, h/2, h/4).
- Mesh objectivity: load-displacement curves on h, h/2, h/4 must coincide within 5% RMS once the characteristic length has been calibrated.

---

## Stage 6 - Tensile Dogbone with Multiple Failure Criteria (Tsai-Wu / Hashin / Puck)

Goal. On a single composite coupon under combined loading, demonstrate that the implementation distinguishes the failure criteria. Tsai-Wu gives a single quadratic failure index; Hashin separates fiber-tension, fiber-compression, matrix-tension, matrix-compression; Puck distinguishes fiber failure (FF) from inter-fiber failure (IFF) on a fracture plane.

Standard. ASTM D3039 (tensile properties), used for coupon-level data; ASTM D3518 for $\pm 45^{\circ}$ in-plane shear; ASTM D6641 for compressive strength.

Reference. Tsai, S.W. and Wu, E.M. (1971), "A general theory of strength for anisotropic materials", *Journal of Composite Materials*, **5**, 58-80. Hashin, Z. (1980), "Failure criteria for unidirectional fiber composites", *Journal of Applied Mechanics*, **47**, 329-334. Puck, A. and Schurmann, H. (1998), "Failure analysis of FRP laminates by means of physically based phenomenological models", *Composites Science and Technology*, **58**, 1045-1067; revised in *CST* **62**, 1633-1662 (2002). Pinho, S.T., Iannucci, L., Robinson, P. (2006), "Physically based failure models and criteria for laminated fibre-reinforced composites with emphasis on fibre kinking. Part I: Development; Part II: FE implementation", *Composites Part A*, **37**, 63-73 and 766-777. Comparative benchmark: the World-Wide Failure Exercise (WWFE-I, II, III) of Hinton, Kaddour, and Soden, *Composites Science and Technology* special issues (1998, 2002, 2013). The WWFE provides tabulated experimental envelopes against which all three criteria can be plotted.

Specimen dimensions. Single $\pm 45^{\circ}$ tensile coupon per ASTM D3518: $L=250$ mm, $w=25$ mm, $t \approx 2$ mm, eight plies, T800/Cycom CFRP or similar, with the $\pm 45^{\circ}$ stacking ensuring all three criteria predict different stress-state interactions.

Failure / success criterion of the physical experiment. The three criteria are scored on (a) predicted first-ply-failure load, (b) predicted failure mode, (c) fracture-plane angle. Experimental match is documented in WWFE-I/II Test Cases.

FEniCSx tutorial. No published FEniCSx tutorial for combined Tsai-Wu / Hashin / Puck. Implementation is via post-processing of the per-ply stress field from a CLT/FEM solution; the per-ply stress comes directly from the orthotropic tour above (Stage 8 below).

Validation success criterion.
- First-ply-failure load from each criterion reproduces the original authors' tabulated values from WWFE Test Case 4 (typically a $\pm 45^{\circ}$ E-glass/MY750 epoxy laminate) within 5%.
- Failure mode flag (FF / IFF / matrix tension / matrix compression) matches the criterion's published prediction on the same load case.
- Tsai-Wu, Hashin, and Puck must each pass their own self-consistency test (e.g. uniaxial fiber-direction tension recovers $X_T$ exactly).

---

## Stage 7 - Anisotropic Single-Tow Tension (UD CFRP, Fiber-Aligned and Off-Axis)

Goal. Verify orthotropic stiffness and strength prediction in fiber, transverse, and off-axis directions. This isolates the assembly of the orthotropic stiffness tensor and the rotation matrix.

Standard. ASTM D3039 / D3039M (Tensile Properties of Polymer Matrix Composite Materials).

Reference. ASTM D3039/D3039M-17 (https://store.astm.org/d3039_d3039m-17.html). Soden, P.D., Hinton, M.J., Kaddour, A.S. (1998), "Lamina properties, lay-up configurations and loading conditions for a range of fibre-reinforced composite laminates", *Composites Science and Technology*, **58**, 1011-1022 - this is the canonical materials-data paper for the WWFE and the standard set of UD-CFRP / UD-GFRP properties used worldwide. Daniel, I.M. and Ishai, O., *Engineering Mechanics of Composite Materials* (2006) - the standard textbook with worked closed-form examples for $0^{\circ}$, $90^{\circ}$, and off-axis tension. Pagano, N.J. and Halpin, J.C. (1968), "Influence of end constraint in the testing of anisotropic bodies", *Journal of Composite Materials*, **2**, 18-31.

Specimen dimensions. Per ASTM D3039: $0^{\circ}$ UD: $L=250$ mm, $w=15$ mm, $t=1$ mm; $90^{\circ}$ UD: $L=175$ mm, $w=25$ mm, $t=2$ mm. Off-axis (e.g. $10^{\circ}$, $30^{\circ}$, $45^{\circ}$, $60^{\circ}$): per ASTM D3039 + Pagano-Halpin oblique-grip correction.

Failure / success criterion of the physical experiment. Longitudinal modulus $E_1$, transverse modulus $E_2$, in-plane Poisson's ratio $\nu_{12}$, and longitudinal/transverse tensile strengths $X_T$, $Y_T$.

FEniCSx tutorial. (a) Bleyer, "Isotropic and orthotropic plane stress elasticity" (https://bleyerj.github.io/comet-fenicsx/tours/linear_problems/isotropic_orthotropic_elasticity/isotropic_orthotropic_elasticity.html). (b) "Orthotropic linear elasticity" tour (https://comet-fenics.readthedocs.io/en/latest/demo/elasticity/orthotropic_elasticity.py.html).

Validation success criterion.
- $0^{\circ}$ test: extracted $E_1$ within 1% of input.
- $90^{\circ}$ test: extracted $E_2$ within 2% (transverse-direction stress concentration around grips makes this slightly looser).
- Off-axis: apparent modulus $E_x(\theta)$ reproduces the closed-form transformation $1/E_x = \cos^4\theta/E_1 + (1/G_{12} - 2\nu_{12}/E_1)\cos^2\theta\sin^2\theta + \sin^4\theta/E_2$ within 2% RMS for $\theta \in \{0, 10, 30, 45, 60, 90\}^{\circ}$.

---

## Stage 8 - Anisotropic Ply with Rotation (Stiffness Transformation Verification)

Goal. Verify the implementation of the in-plane rotation matrix on the $6\times 6$ stiffness $C$ (or $3\times 3$ plane-stress $Q$) for an arbitrary ply orientation $\theta$.

Standard. None directly; algebraic verification.

Reference. Daniel and Ishai, *Engineering Mechanics of Composite Materials* (2006), Chapter 5 - Closed-form $\bar{Q}_{ij}(\theta)$ from $Q_{ij}$. Jones, R.M., *Mechanics of Composite Materials*, 2nd ed. (1999), gives the canonical $T_{\sigma}$ and $T_{\epsilon}$ rotation matrices and their use in $\bar{Q} = T_{\epsilon}^{-1} Q T_{\epsilon}$ (or equivalent forms).

Specimen dimensions. Single rectangular ply, $L=100$ mm, $w=100$ mm, $t=0.125$ mm, single element / very coarse mesh. Test cases at $\theta = \{0, 15, 30, 45, 60, 75, 90\}^{\circ}$.

Failure / success criterion. None - pure algebraic verification.

FEniCSx tutorial. Same as Stage 7. The orthotropic-elasticity tour exposes the rotation step explicitly.

Validation success criterion.
- All six (or three) entries of $\bar{Q}(\theta)$ from FEM match the textbook formula entry-wise to machine precision (relative error $< 10^{-10}$, since this is purely algebraic).
- Round-trip verification: rotate by $\theta$, then by $-\theta$, recover the original $Q$ to machine precision.
- Apparent uniaxial modulus $E_x(\theta)$ from a single-element tensile FEM matches the closed-form formula to machine precision.

---

## Stage 9 - Multilayer Laminate (Cross-Ply, then Quasi-Isotropic) - Classical Lamination Theory

Goal. Verify the laminate-level $A$, $B$, $D$ matrices computed from per-ply contributions match Classical Lamination Theory (CLT). Verify the through-thickness ply-stress recovery from laminate strains.

Standard. ASTM D3039 (tensile), ASTM D3518 (in-plane shear), ASTM D6856 (sampling and statistics for composites).

Reference. Reddy, J.N., *Mechanics of Laminated Composite Plates and Shells: Theory and Analysis*, 2nd ed., CRC Press (2003). Jones, R.M. (1999), Chapters 4-6. Tsai, S.W., *Theory of Composites Design* (1992). MIL-HDBK-17-3 / CMH-17 Volume 3, Chapter 4 gives worked closed-form examples for cross-ply and quasi-isotropic laminates.

Specimen dimensions. (a) Cross-ply $[0/90]_S$, four plies, $t_{\text{ply}}=0.125$ mm, total $t=0.5$ mm, $L=250$ mm, $w=25$ mm. (b) Quasi-isotropic $[0/45/-45/90]_S$, eight plies, total $t=1$ mm. Material: T800/Cycom or Soden-Hinton-Kaddour reference UD properties.

Failure / success criterion of the physical experiment. Effective laminate moduli $E_x$, $E_y$, $G_{xy}$, $\nu_{xy}$; first-ply-failure load and last-ply-failure load.

FEniCSx tutorial. No dedicated CLT tour in FEniCSx is published. Implementation is either (a) layered shell elements with through-thickness Gauss integration, or (b) full 3D solid stack with one element per ply. The Bleyer "Linear shell model" tour (https://bleyerj.github.io/comet-fenicsx/tours/shells/linear_shell/linear_shell.html) is the closest starting point.

Validation success criterion.
- $A$, $B$, $D$ matrices match CLT closed form to within $10^{-8}$ relative error (it is just an integration of a piecewise-constant $\bar{Q}(z)$).
- Cross-ply: $B \equiv 0$ to machine precision, and the in-plane response is uncoupled from bending.
- Quasi-iso: $A_{16} = A_{26} \equiv 0$, $A_{11} = A_{22}$, in-plane behavior is isotropic. $E_x = (A_{11} A_{22} - A_{12}^2)/(t A_{22})$ within 1% of CLT.
- Per-ply stress recovery matches Reddy's Chapter 7 worked example to within 2% RMS.

---

## Stage 10 - Mesoscale RVE of UD Composite (Fiber + Matrix, Periodic BCs)

Goal. From a fiber-and-matrix microstructure, recover the homogenized UD-ply orthotropic stiffness via periodic boundary conditions on a representative volume element (RVE). Verify against analytical Halpin-Tsai and Mori-Tanaka micromechanics.

Standard. None directly; ASTM D3878 (definitions of terms) is the closest. Verification is against published RVE / micromechanics benchmarks.

Reference. Sun, C.T. and Vaidya, R.S. (1996), "Prediction of composite properties from a representative volume element", *Composites Science and Technology*, **56**, 171-179 - the canonical reference for fiber-matrix RVEs with periodic boundary conditions, gives the engineering constants for IM7/8552 and AS4/3501-6 type systems. Xia, Z., Zhang, Y., Ellyin, F. (2003), "A unified periodical boundary conditions for representative volume elements of composites and applications", *International Journal of Solids and Structures*, **40**, 1907-1921. Hashin, Z. and Shtrikman, S. (1963), "A variational approach to the theory of the elastic behaviour of multiphase materials", *J. Mech. Phys. Solids*, **11**, 127-140 - the classical bounds. Halpin, J.C. and Tsai, S.W. (1969), "Effects of environmental factors on composite materials", AFML-TR-67-423 - the engineering closed-form micromechanics. Bohm, H.J. (2021), "A short introduction to basic aspects of continuum micromechanics", ILSB Report.

Specimen / RVE dimensions. 2D plane-strain or generalized plane-strain square cell with one or more fibers. Single-fiber square-packing RVE of side $a = 2 R_f / \sqrt{V_f \pi}$, with $R_f = 3.5\,\mu$m, $V_f = 0.6$. Or a multi-fiber cell with random hexagonal packing, side length $5-10 R_f$.

Failure / success criterion of the physical experiment. Six engineering constants of the homogenized UD ply: $E_1$, $E_2$, $\nu_{12}$, $\nu_{23}$, $G_{12}$, $G_{23}$.

FEniCSx tutorial. (a) Bleyer, "Periodic homogenization of linear elasticity" (https://bleyerj.github.io/comet-fenicsx/tours/homogenization/periodic_elasticity/periodic_elasticity.html) - implements exactly this stage with `dolfinx_mpc`. (b) Older FEniCS version: https://comet-fenics.readthedocs.io/en/latest/demo/periodic_homog_elas/periodic_homog_elas.html.

Validation success criterion.
- $E_1$ (rule of mixtures, fiber-dominated): within 2% of $V_f E_f + (1-V_f) E_m$.
- $E_2$, $G_{12}$ within 5% of Halpin-Tsai's prediction with $\xi=2$ (transverse) and $\xi=1$ (shear).
- All six constants must lie inside the Hashin-Shtrikman bounds.
- Reproduce Sun-Vaidya 1996 Table 1 / Table 2 (their AS4/3501-6 RVE) within 5% RMS.

---

## Stage 11 - Mesoscale RVE of Woven / Pseudo-Woven Composite

Goal. Replace the single-fiber RVE with a tow-and-matrix (mesoscale) RVE that captures interlacing of warp and weft tows (woven), or the AP-PLY tow-undulation pattern (pseudo-woven). Recover the laminate-level ABD matrix (Kirchhoff thin-plate assumption) or the full 3D homogenized stiffness.

Standard. None directly; verification is against published woven-RVE benchmarks.

Reference. Lomov, S.V. et al. (2007), "Meso-FE modelling of textile composites: road map, data flow and algorithms", *Composites Science and Technology*, **67**, 1870-1891 - the canonical paper defining the meso-FE workflow with the TexGen / WiseTex preprocessor. Sherburn, M. (PhD thesis, Nottingham, 2007) - TexGen reference. Espadas-Escalante, J.J. et al. (2017), "On the application of periodic boundary conditions for multiscale analysis of textile composites", *Composite Structures*. Catalanotti, G. et al. (2024), "A mesoscale computational approach to predict ABD matrix of thin woven composites", *Composite Structures* (https://www.sciencedirect.com/science/article/abs/pii/S0263822324001594). For AP-PLY specifically: Falco, O., Lopes, C.S. et al. (2017), "Prediction of delamination onset and growth for AP-PLY composite laminates using the finite element method", *Composites Part A*, **103**, 230-243; Falco, O., Mayugo, J.A., Lopes, C.S. et al. (2022), "Tensile response of AP-PLY composites: a multiscale experimental and numerical study", *Composites Part A*, **161**, 107101.

Specimen / RVE dimensions. Plain-weave RVE: $L_x = L_y = 4$ mm (one weave repeat), tow width 2 mm, tow thickness 0.2 mm, undulation amplitude 0.1 mm, $V_f^{\text{tow}}=0.7$ inside the tow. AP-PLY RVE: per Falco-Lopes 2017, a unit cell of $\sim 8 \times 8$ mm with through-thickness undulating tows.

Failure / success criterion of the physical experiment. Stiffness coefficients $A_{ij}$, $B_{ij}$, $D_{ij}$ of the woven / AP-PLY mesoscale unit cell, plus damage initiation locations.

FEniCSx tutorial. The Bleyer "Periodic homogenization of linear elasticity" tour generalizes; the new ingredient is mesh generation, typically via TexGen → mesh import. No published FEniCSx tutorial does woven RVEs end-to-end.

Validation success criterion.
- ABD matrix from FEM matches Catalanotti et al. 2024 Table 2 (their plain-weave benchmark) within 5%.
- For AP-PLY, predicted laminate $E_x$, $E_y$ in tension within 5% of Falco-Lopes 2022 Table 4 measured values.
- Macro stress-strain in pure tension matches Falco-Lopes 2022 Figure 7 within 5% RMS in the elastic regime.

---

## Stage 12 - Cohesive Zone Delamination (DCB Mode I, ENF Mode II)

Goal. Verify cohesive-zone modeling (CZM) by reproducing the analytic load-displacement curve of a Double Cantilever Beam (DCB, mode I) and an End-Notched Flexure (ENF, mode II) specimen with known $G_{Ic}$ and $G_{IIc}$.

Standard. ASTM D5528 / D5528M (DCB, mode I); ASTM D7905 / D7905M (ENF, mode II); ASTM D6671 (Mixed-Mode Bending, MMB).

Reference. ASTM D5528-21 (https://store.astm.org/d5528_d5528m-21.html). ASTM D7905-19 (https://store.astm.org/d7905_d7905m-19e01.html). Camanho, P.P., Davila, C.G., de Moura, M.F. (2003), "Numerical simulation of mixed-mode progressive delamination in composite materials", *Journal of Composite Materials*, **37**, 1415-1438 - the canonical FEM benchmark for CZM with DCB / ENF / MMB. Turon, A., Davila, C.G., Camanho, P.P., Costa, J. (2007), "An engineering solution for mesh size effects in the simulation of delamination using cohesive zone models", *Engineering Fracture Mechanics*, **74**, 1665-1682. Krueger, R. (2004), "Virtual crack closure technique: history, approach, and applications", *Applied Mechanics Reviews*, **57**, 109-143 - the canonical VCCT reference for verifying CZM against energy-based fracture mechanics. Krueger, R. and Carvalho, N.V. (2018), benchmark NASA TM on quasi-static delamination prediction.

Specimen dimensions. DCB per ASTM D5528: $L=125$ mm, $w=20-25$ mm, $h=3-5$ mm (per arm), pre-crack $a_0=50$ mm, T300/977-2 or AS4/PEEK. ENF per ASTM D7905: $L=140$ mm, span $2L_{\text{span}}=100$ mm, $w=20-25$ mm, $h=3-5$ mm, pre-crack $a_0=20-40$ mm.

Failure / success criterion of the physical experiment. Mode I and mode II interlaminar fracture toughness $G_{Ic}$ and $G_{IIc}$, extracted via Modified Beam Theory, Compliance Calibration, or Modified Compliance Calibration (MBT/CC/MCC) per ASTM D5528 / D7905.

FEniCSx tutorial. (a) Bleyer, "Cohesive zone modeling of debonding and bulk fracture" (https://bleyerj.github.io/comet-fenicsx/tours/interfaces/intrinsic_czm/intrinsic_czm.html). (b) Bleyer, "Cohesive zone modeling restricted to an interface" (https://bleyerj.github.io/comet-fenicsx/tours/interfaces/czm_interface_only/czm_interface_only.html).

Validation success criterion.
- DCB load-displacement curve: peak load within 5%, propagation plateau within 3% of the analytic Modified Beam Theory closed form $P(a) = (b/a)\sqrt{(G_{Ic} E_1 h^3 b)/(12)}$ (or its corrected form per ASTM D5528).
- ENF: peak load and propagation behavior within 5% of the closed-form $P(a) = (2 b h^2/3) \sqrt{(2 G_{IIc} E_1)/(3 a^2 + L_{\text{span}}^2 / 4)}$ or equivalent.
- $G_{Ic}$ extracted from FEM via energy release rate matches the input value within 2%.
- Mesh objectivity: at least three meshes show converged load-displacement curves once the cohesive zone is resolved by 3-5 elements (Camanho-Davila-Turon mesh size criterion).

---

## Stage 13 - Compression After Impact (CAI)

Goal. Combine impact damage and compressive residual strength. Verify that an externally introduced delamination / damage field in a quasi-isotropic plate produces a residual compressive strength matching ASTM D7137 measurements.

Standard. ASTM D7137 / D7137M (Compressive Residual Strength of Damaged Polymer Matrix Composite Plates). Used in tandem with ASTM D7136 (low-velocity impact, see Stage 14) to define the input damage state. NASA RP-1142 also documents the CAI test.

Reference. ASTM D7137-17 (https://store.astm.org/d7137_d7137m-17.html). Soutis, C., Curtis, P.T. (1996), "Prediction of the post-impact compressive strength of CFRP laminated composites", *Composites Science and Technology*, **56**, 677-684 - canonical analytic / semi-empirical CAI strength model. Davies, G.A.O., Hitchings, D., Ankersen, J. (2006), "Predicting delamination and debonding in modern aerospace composite structures", *Composites Science and Technology*, **66**, 846-854. Camanho, P.P. et al. (2008), CMH-17-derived CAI benchmark.

Specimen dimensions. Per ASTM D7137: $L=152.4$ mm (6 in), $w=101.6$ mm (4 in), $t \approx 5$ mm, quasi-isotropic $[45/0/-45/90]_{4S}$ T800/3900-2 or IM7/8552.

Failure / success criterion of the physical experiment. Compressive residual strength after a specified impact energy (typically 6.7 J/mm). Failure is local sublaminate buckling at the impact site followed by global compression failure.

FEniCSx tutorial. None directly. Built up from the orthotropic elasticity (Stage 7), CLT (Stage 9), and CZM (Stage 12) tours.

Validation success criterion.
- Residual compressive strength within 10% of the published Soutis-Curtis 1996 sublaminate-buckling closed-form for the same damage area.
- Buckling mode shape (local sublaminate vs. global) matches the published experimental observation for the same damage size.
- For a baseline (undamaged) specimen, compressive strength matches ASTM D6641 reference data within 5%.

---

## Stage 14 - Low-Velocity Impact (Drop Tower, ASTM D7136)

Goal. Couple explicit transient dynamics, contact, and progressive damage. Verify a drop-tower simulation against ASTM D7136 measurements: peak contact force, contact duration, absorbed energy, and projected delamination area (typically by C-scan).

Standard. ASTM D7136 / D7136M (Damage Resistance of a Fiber-Reinforced Polymer Matrix Composite to a Drop-Weight Impact Event). Companion ISO 18352.

Reference. ASTM D7136-20 (https://store.astm.org/d7136_d7136m-20.html). Olsson, R. (2001), "Analytical prediction of large mass impact damage in composite laminates", *Composites Part A*, **32**, 1207-1215. Olsson, R. (2010), "Closed form prediction of peak load and delamination onset under small mass impact", *Composite Structures*, **94**, 100-107. Davies, G.A.O., Olsson, R. (2004), "Impact on composite structures", *Aeronautical Journal*, **108**, 541-563 - canonical analytic model. Lopes, C.S., Camanho, P.P., Gurdal, Z., Maimi, P., Gonzalez, E.V. (2009), "Low-velocity impact damage on dispersed stacking sequence laminates. Part II: Numerical simulations", *Composites Science and Technology*, **69**, 937-947 - the canonical FEM benchmark with full intra-laminar damage and inter-laminar CZM.

Specimen dimensions. Per ASTM D7136: $L=150$ mm, $w=100$ mm, $t \approx 5$ mm, quasi-isotropic. Hemispherical impactor diameter 16 mm (5/8 in), impact mass 5.5 kg, drop heights yielding 6.7 J/mm of laminate thickness (typical aerospace requirement).

Failure / success criterion of the physical experiment. Peak force, force-time history, contact duration, absorbed energy, projected delamination area (C-scan), back-face indentation.

FEniCSx tutorial. None - explicit dynamics / contact is not in the standard FEniCSx tutorial set. Bleyer, "Explicit dynamics with mass lumping" (https://bleyerj.github.io/comet-fenicsx/) is the closest starting point. For contact, custom implementation is required.

Validation success criterion.
- Peak contact force within 10% of measured / Olsson-Davies analytic prediction.
- Contact duration within 10%.
- Absorbed energy within 15%.
- Projected delamination area within 20% (this is the loosest tolerance because it depends sensitively on stacking sequence and through-thickness CZM resolution).
- Back-face displacement history matches Lopes-Camanho 2009 Figure 9-12 within 10% RMS.

---

## Stage 15 - High-Velocity / Ballistic Impact (V50, Residual Velocity)

Goal. Validate the high-strain-rate constitutive law, projectile-target contact, and full penetration physics. Predict the V50 ballistic limit and the residual velocity at over-match.

Standard. MIL-STD-662F (V50 Ballistic Test for Armor); NIJ Standard-0101.06 / 0101.07 (Ballistic Resistance of Body Armor); STANAG 2920 (NATO Allied Engineering Publication, fragment-simulating projectile testing); ASTM F1233 (Security Glazing). For armor against fragments: MIL-DTL-46593 (FSP geometry).

Reference. MIL-STD-662F (https://everyspec.com/MIL-STD/MIL-STD-0500-0699/MIL-STD-662F_6718/). NIJ Standard-0101.06 (https://www.nist.gov/system/files/documents/oles/ballistic.pdf). The Recht-Ipson model: Recht, R.F. and Ipson, T.W. (1963), "Ballistic perforation dynamics", *Journal of Applied Mechanics*, **30**, 384-390 - canonical analytic model, $V_r^2 = (V_i^2 - V_{50}^2)$ for over-match velocities. Lambert, J.P. and Jonas, G.H. (1976), Lambert-Jonas approximation, *BRL Report* No. 1852. For composites: Cunniff, P.M. (1992), "An analysis of the system effects in woven fabrics under ballistic impact", *Textile Research Journal*, **62**, 495-509. Naik, N.K., Shrirao, P. (2004), "Composite structures under ballistic impact", *Composite Structures*, **66**, 579-590. Karthikeyan, K., Russell, B.P. (2014), *International Journal of Impact Engineering*, on Dyneema panels.

Specimen dimensions. Composite armor panel $L \times w = 152 \times 152$ mm or $300 \times 300$ mm, thickness 4-12 mm. Projectile: 0.30 cal FSP (Fragment Simulating Projectile, $m=2.85$ g) per MIL-DTL-46593, or 9 mm FMJ per NIJ Level II/IIIA.

Failure / success criterion of the physical experiment. V50 (the velocity for 50% probability of penetration) bracketed within $\pm 38$ m/s ($\pm 125$ ft/s) per MIL-STD-662F protocol; residual velocity $V_r$ for over-match shots; backface deformation (BFD) for non-penetrating shots per NIJ.

FEniCSx tutorial. None - this is at the edge of explicit-dynamics capability. For inspiration: Bleyer's explicit dynamics tour, plus custom strain-rate-dependent constitutive law, contact, and erosion.

Validation success criterion.
- V50 within $\pm 5\%$ of measured value (typical aerospace acceptance).
- Residual velocity at $V_i = 1.5 \cdot V_{50}$ within $\pm 10\%$ of Recht-Ipson prediction.
- Backface deformation within 15% for sub-V50 shots.
- Number of plies penetrated, in qualitative agreement with experimental cross-section (within $\pm 1$ ply).

---

## Stage 16 - Pseudo-Woven (AP-PLY) Panel under Projectile Impact (Final Goal)

Goal. Combine all preceding capabilities (orthotropic elasticity, CLT, mesoscale undulation, CZM, progressive damage, explicit dynamics, contact) to predict the ballistic response of an AP-PLY panel.

Standard. NIJ Standard-0101.06 / 0101.07; MIL-STD-662F; STANAG 2920.

Reference. Falco, O. et al. (2017, 2022) (above). Lopes-Camanho-AP-PLY family of papers; Costa et al. (2024), "Low velocity impact response of Automated Fiber Placement Advanced Placed Ply composites", *Composites Part B*, **283**, 111585 (https://www.sciencedirect.com/science/article/pii/S0266353824002069). Tavares, R.P., Bouvet, C., Lopes, C.S., et al. - AP-PLY ballistic and dynamic response work, including the *Composites Science and Technology* 2022 paper "Dynamic response of Advanced Placed Ply composites" (https://www.sciencedirect.com/science/article/pii/S135983682200720X).

Specimen dimensions. AP-PLY panel $300 \times 300$ mm, $t = 4-8$ mm, with the AP-PLY tow-placement architecture per Lopes-Falco. Impactor: 0.30 cal FSP at velocities sweeping V50.

Failure / success criterion of the physical experiment. V50, residual velocity, projected delamination area (C-scan), penetrated-ply count, BFD; modal failure pattern showing crack-arrest by tow undulations (a known AP-PLY feature, per Falco-Lopes).

FEniCSx tutorial. None.

Validation success criterion.
- V50 within $\pm 7\%$ of measured (looser than Stage 15 because AP-PLY architecture introduces additional uncertainty).
- Crack-arrest by undulations qualitatively reproduced (delamination area should be smaller than baseline tape-laminate of equivalent layup, per Falco-Lopes 2024 finding).
- Damage pattern (visualized via per-ply scalar damage variable) matches the experimentally observed C-scan within 25% in projected area.
- Energy partition (kinetic energy absorbed by intra-laminar damage vs delamination) matches the experimentally inferred values within 20%.

---

## Cross-Cutting Notes

1. Mesh refinement protocol. Every stage from 5 onwards must include a mesh-convergence study with at least three refinement levels and a tabulated convergence rate.

2. Material property provenance. Where possible, use the WWFE / Soden-Hinton-Kaddour 1998 material data so that all stages share a common baseline material card, ensuring inter-stage consistency.

3. Verification vs validation. Stages 1, 2, 4 (Kirsch part), 7-9 are *verification* (FEM vs closed form). Stages 3, 5, 6, 10-16 are *validation* (FEM vs experiment).

4. Ascending complexity. Each stage adds at most one new physics capability: linearity → geometric nonlinearity → material plasticity → damage softening → multiple criteria → anisotropy → multilayer → multiscale → interfaces → impact → ballistic.

5. AP-PLY-specific reference set. The Lopes-Falco group at IMDEA / Edinburgh / Oxford / Porto is the canonical source for AP-PLY data. Their papers (Falco 2017 delamination; Falco 2022 tensile multiscale; Costa 2024 LVI; Tavares 2022 dynamic) form the validation backbone for Stages 11 and 16.
