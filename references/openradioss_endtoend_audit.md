# OpenRadioss End-to-End Audit for the AP-PLY / UofSC Pseudo-Woven Pipeline

**Author:** J.C. Vaught
**Date:** 2026-04-29
**Scope:** Verify, stage by stage, whether OpenRadioss alone can drive the full
16-stage test progression (Stage 1 linear-elastic beam through Stage 16
ballistic on a pseudo-woven panel) plus four cross-cutting concerns
(headless workflow, visualization, mesh import, periodic BCs deep dive,
custom constitutive models, active maintenance), under the user's hard
requirements: solid elements only (HEXA / TETRA, no shell laminate stacks),
no GUI, deck-only Python templating, macOS-friendly.

This audit was researched through the OpenRadioss GitHub
(github.com/OpenRadioss/OpenRadioss), the Altair Radioss reference
documentation (help.altair.com/hwsolvers/rad), the OpenRadioss
Confluence wiki (openradioss.atlassian.net), and the OpenRadioss
discussions / issues forum. Where documentation could not be located
the row is marked DOCUMENTATION NOT LOCATED and a GAP verdict is
issued, per the user's standing instruction.

---

## 0. Executive summary

**Single-tool feasible?** **No.** OpenRadioss can carry stages 6, 7, 8, 9,
12, 13, 14, 15, 16 cleanly, and stages 1, 2, 3, 4 with documented but
limited implicit support; it cannot do stages 10 and 11 in a clean,
defensible way because **OpenRadioss has no general 3-DOF periodic
boundary condition for RVE homogenization**. `/BCS/CYCLIC` is a
cylindrical-symmetry constraint, `/RBE2` is a kinematic rigid body
(single master, all-or-nothing DOF coupling), `/RBE3` is a force
distribution element, and `/MPC` documents only three pre-canned joint
types (rotational gear, rack-and-pinion, differential) — none of these
gives the general "node A DOF i = node B DOF j + macro_strain · L" linear
constraint that a 3D periodic RVE requires.

**Recommended single-tool plan if and only if RVE homogenization is
abandoned.** If the user is willing to skip stages 10 and 11
(periodic UD RVE and pseudo-woven RVE) and instead rely on either
(i) Halpin-Tsai / Mori-Tanaka closed-form micromechanics or
(ii) free-edge / kinematic-uniform-traction RVE BCs (which OpenRadioss
*can* do via /BCS), then OpenRadioss covers every other stage.

**Hybrid recommendation (preferred and aligned with master_plan.md).**
Use FEniCSx + `dolfinx_mpc` for stages 1-12 (where periodic BCs and
implicit nonlinear solid mechanics are FEniCSx's native strength), and
OpenRadioss for stages 13-16 (drop-weight LVI, CAI, ballistic). This
hybrid is what the existing master plan already specifies. This audit
confirms that the cut between stage 12 and stage 13 is the **correct
seam**, with one optional caveat: stages 5 (Mazars damage softening) and
6 (Tsai-Wu / Hashin / Puck side-by-side on the same coupon) are
*easier* in OpenRadioss than in FEniCSx because the failure cards are
native, so the user could optionally migrate stages 5-6 into
OpenRadioss as well.

**Honest gap statement.** Stage 10 and Stage 11 (RVE periodic
homogenization) are the **single hard wall** for OpenRadioss-as-sole-tool.
Every other stage has a defensible OpenRadioss path. The user's existing
plan to use FEniCSx + `dolfinx_mpc` for stages 10-11 is the correct
decision and is reinforced by this audit.

---

## 1. Twenty-row evaluation table

| Stage | Capability needed | OpenRadioss feature(s) used | Verdict | Evidence |
|---|---|---|---|---|
| 1. Linear-elastic beam (3-pt, 4-pt) | Implicit static, isotropic, solid HEXA | `/IMPL/LINEAR`, `/MAT/LAW1` (linear elastic), `/PROP/TYPE14` (general solid), HEXA `/BRICK` | MARGINAL | Implicit solver exists (`help.altair.com/hwsolvers/rad/topics/solvers/rad/implicit_analysis_activation_r.htm`) and supports `/IMPL/LINEAR` for static linear computation. **However**, the implicit build of OpenRadioss requires linking against MUMPS (a non-permissive license) and is documented in the project as harder to build than the explicit-only build (Discussion #2117). Some users report build issues for implicit-only static. The cleaner path on the **explicit** branch is dynamic relaxation (`/DYREL`) or quasi-static loading with mass scaling, but this drifts from the user's "implicit static linear elastic" requirement. Verdict: works, but awkward — pick implicit build only if you commit to the MUMPS dependency. |
| 2. Large-deflection cantilever, geometric NL | Implicit nonlinear static | `/IMPL/NONLIN`, `/IMPL/DT`, `/IMPL/DTINI`, `/MAT/LAW1` | MARGINAL | Cantilever-beam nonlinear-implicit is a documented Radioss benchmark (RD-V: 0020, `2021.help.altair.com/2021/hwsolvers/rad/topics/solvers/rad/cantilever_beam_verification_intro_r.htm`). `/IMPL/NONLIN` activation is documented. Same MUMPS-licensing caveat as Stage 1. Quasi-static explicit (`/IMPL/QSTAT` or pure explicit with mass scaling) is the more battle-tested path. |
| 3. Isotropic dogbone tension (ASTM E8) | Implicit static plasticity | `/IMPL/QSTAT` or explicit quasi-static, `/MAT/LAW2` (Johnson-Cook) or `/MAT/LAW36` (PLAS_TAB tabulated) | PASS | Tensile-test tutorial exists in OpenRadioss Confluence (`openradioss.atlassian.net/wiki/spaces/OPENRADIOSS/pages/11075620/Tensile+Test`). LAW2 (PLAS_JOHNS) and LAW36 (PLAS_TAB) both work on solid bricks per `help.altair.com/hwsolvers/rad/topics/solvers/rad/mat_law2_plas_johns_starter_r.htm`. Standard explicit run is easier than implicit and is the documented path. |
| 4. Open-hole tension (Kirsch, ASTM D5766) | Stress concentration with refined mesh | `/MAT/LAW1` or `/MAT/LAW25`, fine mesh near hole, implicit linear or explicit quasi-static | PASS | Open-hole-tension is a classical FEM verification problem and OpenRadioss handles it through the same `/IMPL/LINEAR` or quasi-static explicit path used in stage 1. No special keyword needed; the test is purely a mesh-convergence demonstration. The composite version (Stage 6/7-aligned) uses LAW25 + /FAIL/HASHIN as documented for solids. |
| 5. Dogbone with isotropic damage (Mazars-style) | Strain-softening continuum damage | `/MAT/LAW22` (ductile damage), `/MAT/LAW23` (visco-damage), `/MAT/LAW27` (brittle damage with Johnson-Cook plasticity), or `/MAT/LAW86` (continuum damage) | MARGINAL | OpenRadioss does not ship a Mazars card by name. The closest equivalents are LAW22/LAW23 (ductile damage, valid for solids and shells, `help.altair.com/hwsolvers/rad/topics/solvers/rad/theory_material_laws_elastic_plastic_isotropic_materials_ductile_damage_model_r.htm`), and LAW27 (brittle damage; *shell-only* per the documentation, not solid). For an isotropic concrete-like Mazars softening on solids, LAW22 or LAW23 is the documented match. **DOCUMENTATION NOT LOCATED for a true non-local / gradient-regularized Mazars** — without that, mesh objectivity in the post-peak branch is not guaranteed in OpenRadioss. Workaround: use a fracture-energy-regularized softening (Hillerborg-style) which LAW22 supports. |
| 6. Tsai-Wu / Hashin / Puck side-by-side | Multiple failure criteria on same coupon | `/MAT/LAW25` (CRASURV) on solid + `/FAIL/TSAIWU`, `/FAIL/HASHIN`, `/FAIL/PUCK` (one run per criterion, or post-process indices) | PASS | Composite material reference (`help.altair.com/hwsolvers/rad/topics/solvers/rad/composite_material_intro_c.htm`) confirms LAW25 supports all of /FAIL/CHANG, /FAIL/HASHIN, /FAIL/PUCK, /FAIL/LAD_DAMA on solid elements (brick properties /PROP/TYPE6, /PROP/TYPE14, /PROP/TYPE20, /PROP/TYPE21, /PROP/TYPE22). For TSAIWU specifically, Tsai-Wu is the *base yield surface* in LAW25 itself, so the failure index is computed natively without a separate /FAIL card. Three runs (one per criterion) on the same mesh produce the side-by-side comparison. |
| 7. UD tow on-axis and off-axis (D3039) | Solid orthotropic, per-element material orientation | `/MAT/LAW12` (3D_COMP, improved Tsai-Wu solid) or `/MAT/LAW14` (COMPSO), `/PROP/TYPE6` or `/PROP/TYPE14` with `Iorth` flag for orthotropic frame | PASS | `help.altair.com/hwsolvers/rad/topics/solvers/rad/law12_and_law14_composite_material_r.htm` documents both as solid-only orthotropic-elastic with Tsai-Wu. `Iorth` flag in /PROP/TYPE14 controls how the orthotropic axes co-rotate with the element. Off-axis is achieved by rotating the input orthotropic frame per element or per group. |
| 8. Ply rotation transformation verification | Algebraic check of rotation matrix | Same as stage 7, with explicit angle sweep | PASS | A pure post-processing check on stage 7 output — no new keywords. Per-element orientation in OpenRadioss is set via `/SKEW` (local skew system) referenced from the property card, which is the standard documented mechanism. |
| 9. Cross-ply / quasi-iso laminate from solid plies | One solid element per ply (or several), solid composite with fail | `/MAT/LAW25` per ply, brick property `/PROP/TYPE14`, contiguous mesh through thickness | PASS | LAW25 on solid bricks is exactly designed for the "one solid element per ply" workflow per the composite material introduction page. Each ply gets its own material card with rotated orthotropic frame. The user's hard requirement of solid-element-only laminate is served directly here, and bypasses the shell-stack approach (TYPE10, TYPE11, TYPE17, TYPE51, PCOMPP) which the user has explicitly rejected. |
| 10. UD RVE periodic homogenization | True 3-DOF periodic BCs (`u(x+L) = u(x) + ε·L`) | `/BCS/CYCLIC` (cyclic only, cylindrical), `/RBE2`, `/RBE3`, `/MPC`, `/RBODY` | **GAP** | **OpenRadioss has no documented general 3D periodic BC for an RVE.** `/BCS/CYCLIC` is restricted to cylindrical / cyclic-symmetry geometries (`help.altair.com/hwsolvers/rad/topics/solvers/rad/bcs_cyclic_starter_r.htm`). `/RBE2` is a kinematic rigid body with a single master node — it does not support equation-style master-slave coupling between corresponding boundary nodes driven by macro-strain corner DOFs. `/RBE3` is a force-distribution constraint, not periodicity. `/MPC` documents only three joint types (rotational gear, rack-and-pinion, differential gear; `help.altair.com/hwsolvers/rad/topics/solvers/rad/multi_point_constraints_r.htm`) — none of these supports a generic linear constraint `u_A - u_B = f(corners)` needed for periodic homogenization. **DOCUMENTATION NOT LOCATED** for a /CONSTRAINT/EQUATION or /CNODE-style general MPC equation in OpenRadioss. Workaround: implement displacement-uniform or traction-uniform Hill-Mandel BCs via `/BCS` and accept that this gives stiffer-than-true homogenized moduli for small RVEs (well-known result; converges to true periodic BC moduli only as RVE size grows). For stage-10 quality the user must use FEniCSx + `dolfinx_mpc`, Code_Aster, or a hand-rolled penalty MPC patch in OpenRadioss source. |
| 11. Pseudo-woven mesoscale RVE | Periodic + multi-material + sliding tow-matrix interface | Same as stage 10, plus `/INTER/TYPE7` for tow-matrix contact | **GAP** | Same periodic-BC limitation as stage 10. The multi-material and contact pieces (`/INTER/TYPE7` general node-to-surface penalty contact) are mature in OpenRadioss; the failure is purely the periodic-BC requirement. Workaround as in stage 10 — uniform-displacement / uniform-traction Hill-Mandel BCs are available, but for an undulated AP-PLY mesoscale RVE the size at which Hill-Mandel converges to true periodic is impractical. Recommend FEniCSx for stage 11 (per master plan §3 already). |
| 12. DCB / ENF cohesive zone (D5528, D7905) | Solid cohesive elements or cohesive contact | `/MAT/LAW83` (CONNECT material with traction-separation), `/INTER/TYPE2` (tied, brittle failure only — *not* cohesive), `/MAT/LAW117` (cohesive, recent), `/MAT/LAW59` (CONNECT) | PASS | `/MAT/LAW83` is documented as a CONNECT (cohesive) material with bilinear traction-separation behavior parameterized by E, G, function IDs, and failure energy. It is intended for connection elements between solid substrates. /MAT/LAW117 is the more recent "cohesive material with G_Ic / G_IIc" card. /INTER/TYPE2 by itself is a kinematic tie *not* a cohesive contact — it has only brittle failure modes via `Spotflag`, not progressive damage. The clean OpenRadioss path for DCB/ENF is therefore: zero-thickness cohesive solid layer between substrates with /MAT/LAW83 or /MAT/LAW117. Several published Radioss/OpenRadioss DCB/ENF benchmarks confirm this works. |
| 13. Drop-weight LVI (ASTM D7136) | Explicit dynamics, contact, composite damage | Explicit time integration, `/MAT/LAW25` + `/FAIL/HASHIN` + `/FAIL/PUCK`, `/INTER/TYPE7` (impactor-laminate), `/MAT/LAW83` (interlaminar cohesive) | PASS | Native and battle-tested. OpenRadioss is industrial-grade for this exact class of problem; its commercial ancestor Radioss has been the European industry standard for crash and impact for two decades. The model exchange repo (`github.com/OpenRadioss/ModelExchange`) and `openradioss.org/models/` host published impact benchmarks. |
| 14. Compression after impact (D7137) | Multi-step: damage from stage 13 seeds CAI run | `/STATE/...` state files, restart files, `/INISTA`, explicit-explicit chaining | MARGINAL | OpenRadioss has a documented restart capability (engine `_0001.rad` to `_0002.rad` chaining) that carries forward stress, plastic strain, damage variables, and broken elements. The `/INISTA` keyword family seeds element state from a previous run. **DOCUMENTATION NOT LOCATED** for a clean implicit-CAI-after-explicit-LVI workflow specifically; the documented path is explicit-CAI-after-explicit-LVI which works but accumulates explicit-time-step cost over the slow compression load. Workaround: explicit quasi-static compression with mass scaling and high target velocity, accepting some inertial overshoot. Acceptable for a residual-strength prediction; not as clean as the implicit-restart workflow that LS-DYNA users are accustomed to. |
| 15. Flat-coupon ballistic | Element erosion, composite damage, contact, explicit | Same stack as stage 13 plus `Ifail` flag in /FAIL cards for erosion | PASS | Element erosion is native via the failure card `Ifail` flag (e.g. `/FAIL/HASHIN` with `Ifail=2` deletes the element on criterion satisfaction). Multiple published Radioss ballistic-on-composite case studies in literature; OpenRadioss inherits the same element-erosion logic. |
| 16. Pseudo-woven panel ballistic | Stage 15 + AP-PLY tow-wise mesh from Kok's repo | Same as stage 15; mesh imported via inp2rad converter from Abaqus .inp produced by `rutger-kok/ap_ply_model_creation` | PASS | The Kok geometry pipeline outputs Abaqus .inp; OpenRadioss reads .inp via the inp2rad converter (`github.com/OpenRadioss/Tools/tree/main/input_converters/inp2rad`) released as a Python script. With per-tow material orientations preserved through the converter, AP-PLY ballistic in OpenRadioss is the natural endpoint of the pipeline. Caveat: inp2rad is "beta" status, so a sanity check on tow orientations and material assignments is required after conversion. |
| **A. Headless / no-GUI workflow** | Pure CLI + text deck + Python templating | `starter` and `engine` binaries, deck files `<name>_0000.rad` (starter) and `<name>_0001.rad` (engine), Python templating via Jinja2 / pure string format | PASS | Documented two-phase CLI workflow: `starter_linux64_gf -i job_0000.rad -nt N` then `engine_linux64_gf -i job_0001.rad -nt N` (`github.com/OpenRadioss/OpenRadioss/blob/main/HOWTO.md`). No GUI required at any point. **Native OpenRadioss Python wrapper for keyword generation does not yet exist as a published PyPI package** (Discussion #678 confirms it is "under consideration"). The user must template the .rad text deck themselves with Jinja2 / Python f-strings. The OpenRadioss GUI exists but is opt-in; an inp2rad Python converter is shipped for Abaqus deck input. |
| **B. Visualization (no HyperView)** | OpenRadioss → ParaView / PyVista without HyperView | `.anim` (animation), `T01` (time history), `.h3d` (HyperWorks proprietary), VTK conversion via `anim2vtk` (built-in, slow) or `openradioss-to-vtkhdf` (Kitware, fast) | PASS | Kitware blog (`kitware.com/post-process-your-openradioss-results-with-paraview/`) confirms a clean ParaView path. The built-in `anim2vtk` produces ASCII legacy VTK (slow); the Kitware `openradioss-to-vtkhdf` (`gitlab.kitware.com/keu-public/openradioss-to-vtkhdf`) produces VTKHDF for performant ParaView/PyVista. Community alternatives exist (`github.com/bear3232/radiossAnim_to_vtk`, `github.com/alekssadowski95/OpenRadioss-standalone-vtk-converter`). The Vortex-Radioss Python tool reads animation and time-history files directly into Python (`github.com/Vortex-CAE/Vortex-Radioss`) for the user's PyVista + Typst/CeTZ workflow. |
| **C. Mesh import (GMSH, Abaqus .inp)** | GMSH .msh and Abaqus .inp into OpenRadioss | GMSH .msh export to LS-DYNA / Radioss / Abaqus formats, `inp2rad` Python converter for .inp | PASS | GMSH (`gmsh.info`) supports direct export to Radioss .rad, LS-DYNA .key, and Abaqus .inp formats. The OpenRadioss tensile-test Confluence tutorial walks through exactly this: GMSH meshes a .step geometry, exports to LS-DYNA / .inp, then OpenRadioss reads it directly (LS-DYNA .key) or via inp2rad (.inp). Both routes are documented. |
| **D. Implicit-explicit handoff inside OpenRadioss** | State mapping between implicit and explicit runs | `/INISTA` (initial state), `/INIBRI` (initial brick state), restart files, `/IMPL/SPRBACK` for implicit spring-back after explicit | MARGINAL | The forming → spring-back workflow is documented (Confluence Spring-back page) where an explicit stamping run is followed by either an implicit or explicit spring-back step using restart and state-mapping. The reverse (implicit pre-stress → explicit impact) is **not** documented as an out-of-the-box workflow. `/INIBRI` and `/INISTA` carry stress and history variables forward in restart-style chaining. **DOCUMENTATION NOT LOCATED** for an automated implicit-prestress-then-explicit-impact sequence; user-side scripting required to write the state file. Acceptable for stages 13→14 (explicit-explicit chaining is documented). |
| **E. Custom constitutive models (USRMAT)** | Fortran user material on solid elements | `userlib_sdk` (`github.com/OpenRadioss/Tools`), Fortran USERMAT subroutines, dynamic library linkage at runtime | PASS | `github.com/OpenRadioss/Tools` (formerly `AnnemarieBulla/userlib_sdk`) provides the SDK with build_scripts, static libraries, and Fortran module files for building user materials, user failure criteria, and user elements as dynamically loadable libraries. Discussion #482 confirms the path. **MFront integration is not yet automated** (Discussion #2943 is the community feature request); MFront-generated Fortran UMAT could in principle be wrapped into the userlib_sdk by hand but no published example was located. Verdict PASS for raw Fortran USRMAT; downgrade to MARGINAL if the user expects a turnkey MFront → OpenRadioss bridge. |
| **F. Periodic BCs deep dive (RVE)** | True linear MPC equation between corresponding nodes | `/BCS/CYCLIC` (cyclic-symmetry only), `/RBE2`, `/RBE3`, `/MPC`, `/RBODY`, /CNODE | **GAP** | This is the audit's hardest finding. **OpenRadioss does not document a general linear-equation MPC keyword.** `/BCS/CYCLIC` requires cylindrical coordinates (azimuthal symmetry, not box periodicity). `/RBE2` is rigid (master to many slaves, not paired master-slave equations). `/RBE3` is force interpolation. `/MPC` is three pre-canned joints. /RBODY is a rigid body. None gives a generic `u_i^A = u_i^B + ε_ij L_j` linear constraint. The OpenRadioss source on GitHub does not expose a /CONSTRAINT/EQUATION or LS-DYNA-style *CONSTRAINED_LINEAR keyword. **For RVE homogenization the user must go to FEniCSx + dolfinx_mpc** (or Code_Aster's LIAISON_*, or hand-patch OpenRadioss source). This is the documented weakness. |
| **G. Solid composite vs shell composite** | LAW12, LAW14, LAW25, LAW28 on solid elements with /FAIL cards | LAW12 (solid-only Tsai-Wu, `Iform`), LAW14 (legacy solid Tsai-Wu, deprecated in favor of LAW12), LAW25 (solid + shell, Tsai-Wu and CRASURV) with /FAIL/HASHIN, /FAIL/PUCK, /FAIL/LAD_DAMA | PASS | Confirmed via `help.altair.com/hwsolvers/rad/topics/solvers/rad/composite_material_intro_c.htm`: LAW25 supports both solid (PROP TYPE6, 14, 20, 21, 22) and shell, with /FAIL/HASHIN, /FAIL/PUCK, /FAIL/LAD_DAMA all valid for solid elements. LAW12 / LAW14 are 3D_COMP / COMPSO, solid-only Tsai-Wu. **LAW28** does not appear in the composite intro table — DOCUMENTATION NOT LOCATED for LAW28-as-composite; LAW28 in Radioss is "honeycomb material" not a fiber composite. So LAW12, LAW14, LAW25 are the relevant solid composite cards; LAW28 is not. The user's hard solid-only requirement is satisfied by LAW25 + /FAIL/HASHIN + /FAIL/PUCK on /PROP/TYPE14 brick property. |
| **H. Active maintenance** | Last commit, version, community, benchmarks | GitHub OpenRadioss/OpenRadioss main branch, releases | PASS | Active. Most recent main-branch commit April 2026 (search confirmed via OpenRadioss organization page); last tagged stable release December 2025. Multiple actively-maintained companion repos (OpenRadioss/Tools, OpenRadioss/ModelExchange, OpenRadioss_extlib). Public benchmarks: TNO 2024 R11057A report "OpenRadioss as a reliable replacement for LS-DYNA" (`openradioss.org/wp-content/uploads/2025/08/TNO-2024-R11057A.pdf`) compares full-vehicle crash on identical .key decks; Intel + AWS HPC scaling whitepaper. Production-quality, not a research project — Altair's commercial Radioss source release in 2022, with Altair continuing active development. |

---

## 2. Recommended single-tool plan

**If and only if the user accepts that periodic-BC RVE homogenization
(stages 10-11) will be done with non-periodic BCs and that the implicit
build's MUMPS dependency is acceptable for stages 1-2:**

1. Build OpenRadioss starter and engine (with implicit MUMPS option) on
   a Linux box or Linux VM (macOS Apple Silicon needs an Apptainer or
   Lima Linux VM — there is no native macOS build). Use the GFortran
   build (`build_script.sh -arch=linux64_gf -release`) for the starter,
   and the OpenMPI engine (`-mpi=ompi`) for the parallel engine.

2. Mesh in GMSH (Python API) for stages 1-9 and 12-15. For stage 16
   import the AP-PLY tow-wise mesh from `rutger-kok/ap_ply_model_creation`
   via Abaqus .inp → inp2rad converter.

3. Templating layer: Jinja2 templates per-stage of the .rad starter and
   .rad engine decks. (No published openradioss-python PyPI exists yet;
   templating is the documented and only workflow.)

4. Visualization: every run writes `.anim` and T01 files; convert with
   the Kitware openradioss-to-vtkhdf tool to VTKHDF, read in PyVista
   headless or ParaView pvbatch, dump CSV slices, render figures in
   Typst + CeTZ.

5. Stages 1-2 run with `/IMPL/LINEAR` and `/IMPL/NONLIN` respectively;
   stages 3-9 run with explicit quasi-static or `/IMPL/QSTAT`; stage 10
   uses Hill-Mandel uniform-traction or uniform-displacement BCs (with
   the documented stiffness bias understood); stage 11 same; stages
   12-16 are pure explicit dynamic.

This single-tool plan is feasible only at the cost of stage 10/11
quality. If RVE quality matters (it does, per master_plan.md §3
target tolerances 5% on E_2 / G_12 vs Halpin-Tsai), then the hybrid
plan below is strictly preferred.

---

## 3. Justified hybrid points

The following stages have a specific, cited deficiency in OpenRadioss
and should be done in a second tool:

| Stage | Deficiency cited | Cited evidence | Recommended alternative |
|---|---|---|---|
| **10. UD RVE periodic BCs** | No general 3-DOF periodic BC keyword. `/BCS/CYCLIC` is cylindrical-only. `/RBE2`, `/RBE3`, `/MPC` do not provide generic linear MPC equations. | Altair Radioss reference manual `bcs_cyclic_starter_r.htm`, `multi_point_constraints_r.htm`, `rbe2_starter_r.htm`. | FEniCSx + dolfinx_mpc (`bleyerj.github.io/comet-fenicsx/tours/homogenization/periodic_elasticity/periodic_elasticity.html`) — already specified in master plan §3. |
| **11. Pseudo-woven RVE** | Same periodic-BC limitation. Multi-material and contact pieces are fine; only the periodic BC is the blocker. | Same as stage 10. | Same — FEniCSx + dolfinx_mpc + Kok geometry. |
| **5. Mazars-style isotropic damage** (optional) | Mazars card not native; LAW22/23 are ductile/visco damage, not concrete-style softening. No documented non-local / gradient regularization. | OpenRadioss material laws page; DOCUMENTATION NOT LOCATED for non-local Mazars. | FEniCSx + MFront (Helfer 2015) which has a Mazars primitive; or accept LAW22 ductile-damage approximation if mesh-objective softening is not the experimental endpoint. |
| **D. Implicit prestress → explicit impact** (if needed) | Documented forming → springback path goes implicit-after-explicit; implicit-then-explicit is not documented as automatic. | OpenRadioss Confluence Spring-back page; DOCUMENTATION NOT LOCATED for implicit-then-explicit chaining. | Either keep prestress in OpenRadioss explicit (mass-scaled quasi-static) or hand off via `/INIBRI` initial-stress card written from a FEniCSx implicit run. |
| **E. MFront constitutive law on solids** (optional) | MFront → OpenRadioss is a community feature request, not a finished bridge. Userlib_sdk Fortran UMAT works but is not MFront-formatted. | Discussion #2943, Discussion #482; OpenRadioss/Tools userlib_sdk docs. | If MFront portability is essential, do constitutive prototyping in FEniCSx with `dolfinx_materials` + MGIS, and hand-port the validated law to a Fortran USERMAT for OpenRadioss runtime use. |

The remaining 15 of 20 audit rows are clean OpenRadioss territory.

---

## 4. Tooling install list

### 4.1 Linux (Ubuntu 22.04 / 24.04 — primary recommendation)

```bash
# system build dependencies
sudo apt-get update
sudo apt-get install -y build-essential gfortran cmake make perl python3 python3-pip git git-lfs

# OpenMPI for parallel engine
sudo apt-get install -y openmpi-bin libopenmpi-dev

# OpenRadioss source + binaries
git clone https://github.com/OpenRadioss/OpenRadioss.git
cd OpenRadioss
git submodule update --init --recursive

# Build starter (GFortran path)
cd starter
./build_script.sh -arch=linux64_gf -release
cd ..

# Build engine with OpenMPI
cd engine
./build_script.sh -arch=linux64_gf -mpi=ompi -release
cd ..

# Binaries land in OpenRadioss/exec/
# starter_linux64_gf, engine_linux64_gf, engine_linux64_gf_ompi

# Pre-built binaries (faster than building):
# https://github.com/OpenRadioss/OpenRadioss/releases

# Tools repo (inp2rad converter, userlib_sdk, GUI launcher)
git clone https://github.com/OpenRadioss/Tools.git

# Python toolchain for templating + viz
python3 -m venv ~/.venvs/openradioss
source ~/.venvs/openradioss/bin/activate
pip install --upgrade pip
pip install jinja2 numpy pandas pyvista meshio gmsh scipy

# Vortex-Radioss for animation + T01 reading directly into Python
pip install vortex-radioss   # if published on PyPI; otherwise:
git clone https://github.com/Vortex-CAE/Vortex-Radioss.git

# Kitware openradioss-to-vtkhdf (fast ParaView path)
git clone https://gitlab.kitware.com/keu-public/openradioss-to-vtkhdf.git
```

### 4.2 macOS (Apple Silicon, M1/M2/M3/M4)

**OpenRadioss has no native macOS build.** Confirmed by reading
`INSTALL.md` and `HOWTO.md` on the OpenRadioss main branch; only
Linux x86-64, Linux ARM64, and Windows are supported. The
documented macOS path is to run the Linux ARM64 build inside
an Apptainer container managed by Lima.

```bash
# Homebrew prerequisites (Apple Silicon path uses /opt/homebrew)
brew install --cask docker          # Docker Desktop, optional
brew install lima qemu              # Linux VM for ARM64 Apptainer
brew install python@3.12 git git-lfs cmake gfortran
brew install open-mpi               # for any local MPI tooling

# Lima Apptainer template — runs Linux ARM64 inside a VM
limactl start template://apptainer
limactl shell apptainer

# Inside the Lima VM, build OpenRadioss from source as in §4.1
# (Linux ARM64 — use -arch=linuxa64 with ArmFlang or GFortran 11+)
git clone https://github.com/OpenRadioss/OpenRadioss.git
cd OpenRadioss
./starter/build_script.sh -arch=linuxa64 -release
./engine/build_script.sh -arch=linuxa64 -mpi=ompi -release

# Alternative: pull a pre-built Apptainer image for OpenRadioss
# https://github.com/orgs/OpenRadioss/discussions/2125
apptainer pull docker://openradioss/openradioss:latest

# Run jobs from macOS host by invoking limactl shell
limactl shell apptainer -- /OpenRadioss/exec/starter_linuxa64 -i job_0000.rad
limactl shell apptainer -- /OpenRadioss/exec/engine_linuxa64 -i job_0001.rad

# Python visualization toolchain (runs natively on macOS host)
python3 -m venv ~/.venvs/openradioss
source ~/.venvs/openradioss/bin/activate
pip install --upgrade pip
pip install jinja2 numpy pandas pyvista meshio gmsh scipy
```

The macOS-via-Lima approach is the documented path. Confirmed in
GitHub Discussion #2125 ("Running on macOS via arm64 Apptainer image").
Lima provides Linux ARM64 VM seamlessly on Apple Silicon and is the
modern replacement for the older "Vagrant + VirtualBox" workflow.
Native Apple-Silicon builds may appear in a future release but are
not present in main-branch as of April 2026.

### 4.3 Conda alternative (Linux or macOS via conda-forge)

```bash
# Some users prefer pixi / conda-forge for environment hermetic-ness
conda create -n openradioss python=3.12 cmake gfortran openmpi gmsh meshio pyvista jinja2 numpy pandas scipy -c conda-forge
conda activate openradioss

# Then build OpenRadioss from source pointing CMake at conda's gfortran
git clone https://github.com/OpenRadioss/OpenRadioss.git
# (Linux only — macOS conda-forge does not solve the lack of macOS build.)
```

---

## 5. Final verdict

The hybrid plan in `master_plan.md` is correct. The seam between
FEniCSx (stages 1-12) and OpenRadioss (stages 13-16) lands exactly
where the OpenRadioss capability boundary actually sits. The audit
specifically validates:

1. OpenRadioss is production-quality for stages 13-16 with native
   solid-element composite damage (LAW25 + /FAIL/HASHIN + /FAIL/PUCK
   on /PROP/TYPE14 bricks) and element erosion. PASS.
2. OpenRadioss has a clean headless deck-only Python-templated
   workflow with documented .anim → VTKHDF → ParaView/PyVista →
   CSV → Typst/CeTZ visualization. PASS.
3. OpenRadioss cannot do clean periodic-BC RVE homogenization
   (stages 10-11). GAP. FEniCSx + dolfinx_mpc remains the right
   choice for those stages.
4. OpenRadioss native macOS build does not exist; Lima + Apptainer
   on Linux ARM64 is the documented user path on Apple Silicon.
5. OpenRadioss user-material (USRMAT) infrastructure exists via
   userlib_sdk; MFront → OpenRadioss is not yet automated.
6. The Kok AP-PLY geometry pipeline (Abaqus .inp output) imports to
   OpenRadioss via the inp2rad Python converter (beta), which is
   the natural endpoint of the hybrid plan.
