# Open-Source Composites Pipeline Ecosystem Survey

**Author:** J.C. Vaught
**Date:** 2026-04-29
**Status:** Wide ecosystem audit. Companion to `fea_alternatives.md`. This file
documents what was missed by the FEM-mainstream survey, with emphasis on
particle / meshfree methods, draping, DEM, GPU Python, cohesive / fracture
codes, visualisation, composite material libraries, and a brutal license audit.

All commercial, restricted-academic, and government-cleared tools are
**rejected**. AGPL is flagged as a *conditional accept* because it forbids
SaaS-style hosting of unmodified-source forks.

---

## 0. Scope and Method

The previous survey (`fea_alternatives.md`) catalogued fifteen mainstream
FEM packages: CalculiX, Code_Aster, Kratos, Akantu, MoFEM, deal.II, FEAP,
Project Chrono, Karamelo, Uintah, Anura3D, GetFEM, NGSolve, OOFEM, FEniCSx.
This file goes wider, covering:

1. Meshfree / particle methods (peridynamics, SPH, MPM)
2. Preform / draping / forming open-source codes
3. Discrete-element method codes and DEM-FEM bridges
4. Modern Python / GPU FEM and physics simulators
5. Cohesive / fracture-specific codes
6. Visualisation beyond ParaView / PyVista
7. Composite material model libraries
8. License audit of every entry
9. Updated end-to-end stack recommendation

For every package listed below:

- License is verified against a primary source where possible (GitHub LICENSE
  file, project homepage). Where the primary source could not be reached,
  the entry is marked **UNVERIFIED**.
- "Commercial use OK?" answers a hypothetical: can a small consulting firm
  use the tool to deliver paid work to an external client?
- "Network-use clause?" flags AGPL-style copyleft that triggers on SaaS
  deployment.
- "Composite ballistic readiness" is rated 0-5 against the AP-PLY pipeline
  specifically, not the package's general capability.

---

## 1. Meshfree and Particle Methods

### 1.1 Peridynamics

#### Peridigm (Sandia National Laboratories)

- **Repository:** github.com/peridigm/peridigm
- **License:** **3-clause BSD** (verified via project description, JOSS-style
  paper Littlewood-Parks-Mitchell-Silling 2023). Permissive, commercial use
  OK, no network clause.
- **Language / build:** C++ on top of Trilinos. Cubit mesh input, ParaView
  output. No native Python API. Driven by ASCII text input file.
- **Capabilities:** Bond-based, ordinary state-based, and non-ordinary
  state-based peridynamic models. Pre-cracks, contact between bodies,
  thermomechanics. Massively parallel via MPI through Trilinos.
- **Composite support:** None natively in core. The literature (Hu 2011 and
  followers, Madenci-Oterkus 2014, Hu 2023) extends bond-based PD to
  fiber-reinforced lamina with separate fiber and matrix bonds — these
  models are *not* in the upstream Peridigm distribution and require
  source-level addition.
- **Ballistic / impact:** Strong for dynamic fragmentation in homogeneous
  brittle solids; impact case studies in published literature exist
  (Parks et al. 2008, Silling-Askari 2005). Composite ballistic with
  Peridigm requires the fiber-matrix bond model to be coded into a custom
  material module.
- **Maintenance:** Active commits as of late 2025; managed by Sandia CCR.
  No formal versioned releases ("No releases published" on GitHub master).
- **Composite ballistic readiness:** **3 / 5** (good fragmentation physics,
  but composite constitutive law must be added).

#### PDLAMMPS (LAMMPS PERI package)

- **Repository:** lammps.sandia.gov; PERI package shipped with LAMMPS proper.
- **License:** **GPL v2** (LAMMPS license).
- **Capabilities:** Peridynamic atom style with `peri/pmb`, `peri/lps`,
  `peri/eps`, `peri/ves` pair styles. Originally Mike Parks (Sandia), with
  EPS / VES extensions by Rahman and Foster (UTSA). Documented in
  Parks-Lehoucq-Plimpton-Silling 2008 user guide.
- **Composite support:** No fiber-matrix peridynamic model natively. Authors
  note this is a research extension (LaRC05-PD by Hu 2023 is paper-only).
- **Ballistic / impact:** Built-in dynamic-fracture demonstrations
  (Kalthoff-Winkler etc.). LAMMPS scripting is text-based, no native Python
  API but `mpi4py`-style coupling is feasible.
- **Composite ballistic readiness:** **2 / 5**.

#### PeriPy (Alan Turing Institute / Exeter / Cambridge)

- **Repository:** github.com/alan-turing-institute/PeriPy
- **License:** **MIT** (Boys et al. 2021, JOSS / CMAME). Permissive, no
  restrictions.
- **Capabilities:** Bond-based peridynamics with OpenCL backend (CPU + GPU
  agnostic). Pure Python user interface with Cython + OpenCL kernels.
  Documented support for **composite and interface material models** via
  multiple bond types per simulation — exactly the structure needed for
  fiber + matrix bonds in a CFRP lamina.
- **Composite support:** *Best-in-class for an open-source Python-driven
  peridynamic composite code.* The `nbond_types` parameter directly
  encodes per-bond elastic and damage parameters; users supply fiber and
  matrix bond moduli at construction.
- **Ballistic / impact:** Bond-based, so explicit-dynamic by construction.
  Authors target sensitivity / UQ workflows, so single ballistic-impact
  cases are not the publication focus, but the kernel is suited.
- **Maintenance:** Active through 2024 commits; small group.
- **Composite ballistic readiness:** **4 / 5** (Python + OpenCL +
  multi-bond-type makes this the lightest-weight viable peridynamic path).

#### CabanaPD (ORNL)

- **Repository:** github.com/ORNL/CabanaPD
- **License:** **3-clause BSD**.
- **Capabilities:** Peridynamics built on the Cabana / Kokkos performance-
  portable particle library. Bond-based PMB, state-based LPS,
  thermomechanics, plasticity, DEM-PD contact, multi-material systems
  (added in v0.4, 2025). MPI parallelism and GPU portability via Kokkos
  (CUDA / HIP / SYCL).
- **Composite support:** Multi-material systems exist; explicit composite
  ply law not advertised but the architecture supports per-bond
  parameters.
- **Ballistic / impact:** Designed for HPC at scale. ORNL group publishes
  on dynamic fracture; specific composite ballistic case studies were not
  located.
- **Maintenance:** Very active; 2025 v0.4 released.
- **Composite ballistic readiness:** **4 / 5** (HPC-class, BSD, multi-material).

### 1.2 Smoothed Particle Hydrodynamics (SPH) and related

#### DualSPHysics

- **Repository:** github.com/DualSPHysics/DualSPHysics
- **License:** **LGPL** (verified on repo). Commercial use OK; modifications
  to DualSPHysics itself must be redistributed under LGPL.
- **Capabilities:** GPU SPH (CUDA + OpenMP). Free-surface flows, fluid-
  structure interaction, dynamic boundaries, frictional contacts with mesh
  bodies. Primarily fluid-mechanics focused.
- **Composite ballistic readiness:** **1 / 5** (not a structural code;
  could host a projectile-vs-armour fluid-impact case but not composite
  ply mechanics).

#### PySPH (Indian Institute of Technology Bombay)

- **Repository:** github.com/pypr/pysph
- **License:** **3-clause BSD**.
- **Capabilities:** Pure Python SPH framework, Cython + PyOpenCL backends.
  WCSPH, IISPH, DFSPH, GTVF, EDAC, transport-velocity. Includes solid-
  mechanics SPH (Total-Lagrangian SPH for elastic / elastoplastic solids)
  per the documentation.
- **Composite ballistic readiness:** **2 / 5** (Total-Lagrangian SPH is
  a candidate path; composite anisotropic constitutive law is not in the
  shipped library).

#### SPlisHSPlasH (Bender et al., RWTH Aachen)

- **Repository:** github.com/InteractiveComputerGraphics/SPlisHSPlasH
- **License:** **MIT**.
- **Capabilities:** WCSPH, PCISPH, PBF, IISPH, DFSPH, PF; strong viscosity,
  surface-tension, vorticity. Computer-graphics fluid-simulation focus.
- **Composite ballistic readiness:** **0 / 5** (graphics fluids).

### 1.3 Material Point Method (MPM)

#### CB-Geo MPM (University of Cambridge)

- **Repository:** github.com/cb-geo/mpm
- **License:** **MIT** (verified — `license.md` is permissive plus a
  Developer Certificate of Origin for contributions).
- **Capabilities:** C++ MPM with MPI parallelism and HPC focus. Single-phase
  and two-phase soil models, B-bar and CPDI, contact algorithms.
  Geomechanics-oriented.
- **Composite support:** None natively.
- **Composite ballistic readiness:** **1 / 5**.

#### Taichi Elements (Hu et al., MIT / former Taichi Graphics)

- **Repository:** github.com/taichi-dev/taichi_elements
- **License:** **MIT** (verified).
- **Capabilities:** MLS-MPM continuum physics in Taichi. Multi-material:
  fluid, sand, snow, elastic / elastoplastic solids. GPU on CUDA / Metal /
  Vulkan via Taichi runtime.
- **Composite support:** None natively. Anisotropic damage is *not* in
  Taichi Elements; it lives in `ziran2020` (AnisoMPM) which is GPL-3 and
  research-only.
- **Composite ballistic readiness:** **2 / 5** (extensible but research-
  level; would need an AP-PLY-specific anisotropic damage law).

#### AnisoMPM / ziran2020 (UPenn Graphics)

- **Repository:** github.com/penn-graphics-research/ziran2020
- **License:** **Apache 2.0** with patent grant (verified: ziran2020 LICENSE).
  Commercial use OK.
- **Capabilities:** SIGGRAPH 2020 anisotropic damage MPM (Wolper et al.).
  Geometric CDM with structural tensors for transversely-isotropic and
  orthotropic materials. Fracture demonstrated on cheese, oranges, meat,
  composite-like materials.
- **Composite support:** Yes, *anisotropic damage in MPM is exactly the
  AP-PLY problem statement* — but the code is research-grade, no
  documentation for ballistic at large scale, no validation against
  ASTM-coupon data.
- **Composite ballistic readiness:** **3 / 5** (right physics, immature
  packaging, single-paper-codebase).

#### NVIDIA Warp

- **Repository:** github.com/NVIDIA/warp
- **License:** **Apache 2.0**. Commercial use OK; patent grant included.
- **Capabilities:** Python framework that JIT-compiles Python kernels to
  CUDA / CPU. Built-in `warp.fem` Galerkin FEM (diffusion, elasticity,
  fluid, level-sets, contact-via-nonconforming-mesh). MPM / particle-in-
  cell quadratures supported as a primitive. Recent v1.12 (2026).
- **Composite support:** None native; Warp is a *primitive library*, not
  a turnkey composites code. AP-PLY laws would be coded in Warp Python
  kernels.
- **Composite ballistic readiness:** **3 / 5** (excellent foundation,
  significant DIY required).

#### GeoTaichi (Yihao Shi et al., Hong Kong Polytechnic / Wuhan Univ.)

- **Repository:** github.com/Yihao-Shi/GeoTaichi
- **License:** **GPL v3.0**. Commercial use OK only if downstream is also
  GPL — *flag for the user's preference*.
- **Capabilities:** Taichi-based DEM, MPM, MPDEM (coupled MP-DEM), and FEM.
  Geophysical / geotechnical focus. v0.4 (Aug 2025) added more constitutive
  models. GPU + CPU portable.
- **Composite support:** None for fiber composites; soil and granular focus.
- **Composite ballistic readiness:** **1 / 5** (geomechanics).

#### Matter (Blatny & Gaume, EPFL / Swiss Slope Stability)

- **Repository:** unspecified in search; published EGUsphere 2025.
- **License:** **UNVERIFIED**.
- **Capabilities:** Open-source MPM solver for granular matter — cohesive
  and compressible granular media. Geomechanics-focused.
- **Composite ballistic readiness:** **1 / 5**.

#### DiffTaichi

- **Repository:** github.com/taichi-dev/difftaichi
- **License:** Inherits Taichi (**MIT**). Commercial use OK.
- **Capabilities:** Differentiable MPM and ten other simulators, ICLR 2020.
  Compatible-Particle-In-Cell + MLS-MPM with differentiation through
  source-code transformation.
- **Composite support:** None native.
- **Composite ballistic readiness:** **2 / 5** (relevant if differentiable
  parameter calibration of the AP-PLY model is the goal).

### 1.4 Particle / FEM hybrids

#### Kratos `MPMApplication` and `DEMApplication`

(Already in primary survey; reaffirmed here as the strongest open-source
single-stack candidate.) BSD license, Python-driven, MPM + FEM + DEM in
the same framework. v10.4 (Dec 2025).

#### KratosMultiphysics MPM-FEM and MPM-DEM coupling

Documented in CIMNE publications (Singer et al. 2022 KratosMPM; Coppe 2023
MPM-FEM coupling for protective structures). Confirmed under the same BSD
license. The MPM-DEM partitioned coupling (Coppe et al.) targets
mass-movement hazards but the architecture transfers to projectile
fragmentation against composite plates.

---

## 2. Preform / Draping / Forming Simulation

This is a genuinely thin slice of open-source software. Most production
draping code is commercial (PAM-FORM, AniForm, FiberSIM, Digimat-FE).

### 2.1 TexGen (University of Nottingham)

- **Repository:** github.com/louisepb/TexGen + SourceForge mirror.
- **License:** **GPL v2** (verified). Commercial use OK only if downstream
  is GPL — *flag for the user*.
- **Capabilities:** Geometric textile modelling. Defines yarn paths,
  cross-sections, periodic unit cells for plain weaves, twills, satins,
  3D woven, 3D orthogonal, knit. Outputs `.vtu` (VTK), `.stl`, ABAQUS `.inp`,
  ANSYS `.cdb`, OpenFOAM-ready meshes. **TetGen + Triangle bundled** for
  volumetric and surface meshing.
- **Status:** v3.13.1 released Aug 2023; Python 3 port complete; modest
  recent commit volume.
- **AP-PLY relevance:** TexGen does *not* model AP-PLY's pseudo-woven
  tow-overlap-and-tuck topology directly, but its yarn-path API is general
  enough to generate a custom AP-PLY tow centerline + cross-section, then
  emit a periodic mesh suitable for FEniCSx / Kratos / OOFEM.
- **License classification:** **ACCEPT** for academic and internal use;
  *caveat* that any modified TexGen redistributed must remain GPL.

### 2.2 KinDrape (Krogh, Bak, Lindgaard, Aalborg University)

- **Repository:** github.com/chrkrogh/KinDrape
- **License:** **UNVERIFIED** (LICENSE file present but content not
  inspected); paper requests Zenodo DOI citation.
- **Capabilities:** Kinematic (fishnet) draping of a fabric onto a mould.
  MATLAB primary, Python and Octave ports. Educational / minimal-code
  implementation.
- **AP-PLY relevance:** Kinematic draping ignores fabric shear stiffness
  and yields a first-pass tow-orientation map only. For AP-PLY this is a
  *useful preliminary* but does not capture the tucking that defines the
  pseudo-woven architecture.
- **Composite ballistic readiness:** N/A (not a stress code); useful as
  geometry generator.
- **License classification:** **CONDITIONAL ACCEPT** — verify LICENSE file
  before commercial use.

### 2.3 fabrics-drape-data (virtualtextiles)

- **Repository:** github.com/virtualtextiles/fabrics-drape-data
- **License:** **UNVERIFIED**; declared as an open exchange format for
  textile mechanical parameters.
- **Capabilities:** Data-format only, with example Python utilities. Not a
  solver.
- **AP-PLY relevance:** Could be the calibrated-parameter container for
  AP-PLY tow properties.

### 2.4 WiseTex (KU Leuven, Verpoest / Lomov)

- **Repository:** Not on GitHub. Distributed by KU Leuven MTM via per-user
  email request.
- **License:** **REJECTED** (no clear public license; "WiseTex - ABAQUS
  convertor" page indicates restricted distribution; KU Leuven LRD
  technology transfer office handles licensing). Per the user's
  no-restricted-academic-license rule, **REJECTED**.

### 2.5 TexMind (commercial)

- **License:** **REJECTED** (commercial product including TexMind Braider
  and Loop3D).

### 2.6 Honest gap statement on AFP path-planning

Despite repeated searches, **no production-quality OSI-licensed AFP
path-planning code was located**. Academic codes (Halbritter dissertation,
Stewart 2019 cad-journal) are described in publications but not released
as open source. AddPath (commercial) is Python-scripted but not OSI-
licensed. The user's "AP-PLY tow-by-tow path generation" must therefore
either be (i) custom-coded on top of TexGen primitives, (ii) custom-coded
from a NURBS surface and a Lau 2023 / Nagaraj 2017 algorithm, or
(iii) borrowed from a closed-source tool with permission.

This is a real open-source gap. It does *not* go away by switching to
FEniCSx or Kratos — neither has a draping module either. It is a separate
problem.

### 2.7 Ad-hoc draping via cloth simulation (Blender, MPM cloth)

Blender's cloth solver (`bpy.ops.ptcache.bake`) is **GPL v2/v3 dual**.
Output is FBX / Alembic, which can be converted to a triangle mesh for
FEM import. This is feasible for *visualisation* of AP-PLY fabric drape
but is not a validated mechanical-draping code. **NOT RECOMMENDED** as a
load path for stress analysis.

---

## 3. Discrete Element Method

### 3.1 Yade (Université Grenoble Alpes / national consortium)

- **Repository:** gitlab.com/yade-dev/trunk
- **License:** **GPL v2**. Commercial use OK only if downstream is GPL —
  flag.
- **Capabilities:** Open-source DEM with Python scripting (the *primary*
  user interface). Particle shapes: spheres, polyhedral, clumps, deformable
  polyhedra. Cohesion, breakage, contact-law plugins.
- **DEM-FEM bridges:** `dem-fem-coupling` (Stránský 2013, repo
  github.com/stranskyjan/dem-fem-coupling) couples Yade with **OOFEM**.
  Both C++ + Python; documented FEMxDEM hierarchical multiscale tutorial
  in Yade docs. Yade-OpenFOAM (CFD-DEM) also exists (`dpkn31/Yade-OpenFOAM-coupling`).
- **AP-PLY relevance:** DEM is not the natural method for ply mechanics,
  but Yade-OOFEM coupling is a credible path for tow-particle interaction
  with continuum-FEM matrix.

### 3.2 LIGGGHTS-PUBLIC (DCS Computing, originally TU Graz)

- **Repository:** github.com/CFDEMproject/LIGGGHTS-PUBLIC
- **License:** **GPL v2 or later**.
- **Capabilities:** DEM extending LAMMPS. Built-in cohesion, mesh-particle
  contact, heat transfer. Commercial successor *Aspherix* is closed; the
  PUBLIC fork is maintained and free. Status: maintenance-only since DCS
  pivoted to Aspherix.
- **AP-PLY relevance:** Same as Yade — not the natural fit for ply
  composite mechanics.

### 3.3 MercuryDPM (University of Twente / Manchester)

- **Repository:** github.com/MercuryDPM (multiple mirrors)
- **License:** **3-clause BSD** (verified via official site).
- **Capabilities:** Object-oriented C++ DEM. Strong granular-flow focus.
  Active.

### 3.4 MUSEN (Technische Universität Hamburg)

- **Repository:** github.com/msolids/musen
- **License:** **3-clause BSD** (verified).
- **Capabilities:** GPU-accelerated DEM with CUDA. Cohesive bonds for
  agglomerate breakage. Granular-process focus.

### 3.5 CFDEMcoupling (DCS Computing)

- **License:** **GPL v3** for the open-source version.
- **Capabilities:** OpenFOAM + LIGGGHTS coupling for CFD-DEM. Not directly
  relevant to AP-PLY ballistic but useful if a fluid-loaded preform
  step were needed.

### 3.6 PhasicFlow (open-source GPU-DEM, Norouzi et al.)

- **License:** **UNVERIFIED**, declared open-source in JOSS-style 2023
  paper. Multi-architecture (CPU / GPU / heterogeneous).
- **Status:** Younger than Yade / LIGGGHTS but actively developed.

### 3.7 Project Chrono DEM module

(Already in primary survey.) BSD-3, supports flexible-body coupling, but
not designed for fine-mesh continuum failure.

---

## 4. Modern Python and GPU FEM / Physics

### 4.1 JAX-FEM (Tianju Xue, HKUST)

- **Repository:** github.com/deepmodeling/jax-fem (and the older
  tianjuxue/jax-am)
- **License:** **GNU GPL v3** (verified). For commercial use, contact
  Tianju Xue. **CAUTION**: this is a *strong copyleft* license — any
  composite simulation code statically or dynamically linked to JAX-FEM
  must also be GPL'd if redistributed.
- **Capabilities:** Differentiable, GPU-accelerated 3D FEM in pure Python
  on JAX. Mesh and finite-element abstractions; automatic differentiation
  for inverse problems; ~10x faster than commercial reference on a 7.7 M
  DOF tensile case. Quasi-static implicit; explicit dynamics not the
  primary use case.
- **Composite support:** No native composite damage; user-implemented.
- **Composite ballistic readiness:** **2 / 5** (great for inverse design /
  parameter ID; not a ballistic explicit code).

### 4.2 Firedrake

- **Repository:** github.com/firedrakeproject/firedrake
- **License:** **LGPL v3 or later** (verified).
- **Capabilities:** Code-generation FEM with high-order DG, extruded
  meshes, mature solver composability. Sister project to FEniCSx, shares
  UFL. Active 2026 release stream.
- **Composite support:** No off-the-shelf ply law; same DIY status as
  FEniCSx.
- **Composite ballistic readiness:** **2 / 5** (mathematician-grade, no
  contact / erosion).

### 4.3 scikit-fem (Tom Gustafsson, Aalto)

- **Repository:** github.com/kinnala/scikit-fem
- **License:** **3-clause BSD** (verified). Permissive.
- **Capabilities:** Pure Python finite-element assembly library. Excellent
  for academic prototyping; the user supplies the time integrator,
  contact, damage. No GPU.
- **Composite ballistic readiness:** **1 / 5** (a teaching / research
  toolbox, not a production ballistic code).

### 4.4 FElupe (Andreas Dutzler, TU Graz)

- **Repository:** github.com/adtzlr/felupe
- **License:** **GPL v3**.
- **Capabilities:** Finite-element nonlinear solid mechanics in pure Python
  (NumPy / SciPy). Hyperelasticity, large deformation, AceGen-style
  symbolic differentiation. **Has been used for cord-rubber composite
  air-spring damage analysis** (companion repo `adtzlr/fiberreinforcedrubber`).
- **Composite support:** Hyperelastic fiber-reinforced rubber demo; not
  CFRP-grade.
- **Composite ballistic readiness:** **1 / 5** (quasi-static, GPL).

### 4.5 SfePy

- **Repository:** github.com/sfepy/sfepy
- **License:** **New BSD (3-clause)**. Permissive.
- **Capabilities:** Python FEM driver on NumPy / SciPy / PETSc / mpi4py.
  Multi-scale and homogenization examples. Implicit primary; explicit
  feasible.
- **Composite / fracture support:** No native cohesive / XFEM as far as
  the search returned. PyFEM (`jjcremmers/PyFEM`) ships cohesive-zone
  examples and is a separate Python FEM educational code.
- **Composite ballistic readiness:** **1 / 5**.

### 4.6 Brax / MJX / Genesis (differentiable robotics simulators)

- **Brax:** Apache 2.0; rigid-body only. Not a continuum-mechanics code.
- **MJX:** Apache 2.0; MuJoCo XLA in JAX; rigid body, articulated.
- **Genesis:** Apache 2.0; differentiable MPM tool solver — *but no
  composite damage mechanics*; SIGGRAPH-style cloth, soft-body, and
  rigid-body. Not a ballistic-impact code.
- **Composite ballistic readiness:** **0 / 5** for all three.

### 4.7 PolyFEM (Jérémie Dumas et al.)

- **Repository:** github.com/polyfem/polyfem (+ polyfem-python bindings).
- **License:** **MIT** (project itself; some linked third-party libs vary).
- **Capabilities:** Polyvalent C++ FEM. **IPC contact** built in via the
  `ipc-toolkit` (separately MIT-licensed). Strong implicit nonlinear
  elastodynamics. Curved-mesh support.
- **AP-PLY relevance:** IPC contact is *intersection- and inversion-free*,
  attractive for projectile-vs-laminate interaction at low-to-moderate
  velocities. PolyFEM is implicit, so for explicit ballistic the IPC math
  must be rewritten — there is no production explicit-IPC code yet.
- **Composite ballistic readiness:** **2 / 5** (good for low-velocity
  drop-weight contact; not for explicit ballistic).

### 4.8 OpenSees / OpenSeesPy

- **Repository:** github.com/OpenSees/OpenSees and zhuminjie/OpenSeesPy
- **License:** **REJECTED for commercial redistribution.** OpenSees is
  free for research, education, and internal use, but **commercial
  redistribution requires a license** per the project FAQ. This violates
  the user's "no restricted-academic-license" rule.
- **Capabilities:** Earthquake-engineering-focused FEM with strong
  Python (OpenSeesPy) interface. Not designed for composite ballistic.

### 4.9 OpenRadioss (Altair, AGPL fork of Radioss)

- **Repository:** github.com/OpenRadioss/OpenRadioss
- **License:** **GNU AGPL v3.0** (verified). **CONDITIONAL ACCEPT.**
  - Commercial use is allowed.
  - **Network use clause: any modified version made available over a
    network must publish source.** This is the "ASP loophole closer."
  - Linking your own pre/post-processing scripts to AGPL OpenRadioss does
    *not* in itself trigger AGPL on your scripts (they communicate via
    text input deck), but a custom user material library compiled against
    the OpenRadioss userlib_sdk is at risk of derivative-work classification.
  - **Verdict:** Acceptable for internal research and consulting *if*
    no SaaS / web-API hosting of OpenRadioss is planned.
- **Capabilities:** Industrial-grade explicit dynamics solver, open-sourced
  by Altair Sept 2022. Same `MAT_LAW_25` (CRASURV composite Tsai-Wu),
  `MAT_LAW_15` (Chang-Chang), and `LAW_27`-and-friends as commercial
  Radioss. Element erosion via failure cards (`/FAIL/HASHIN`,
  `/FAIL/TSAIWU`, `/FAIL/PUCK`, `/FAIL/JOHNSON`, etc.) is **native**.
  Reads its own `.rad` plus LS-DYNA `.k/.key` and Abaqus `.inp`. Python
  interface added 2024 for keyword scripting (`-python` mode for the
  starter). Latest stable 20260319 build.
- **Composite ballistic readiness:** **5 / 5.** This is the only OSI-
  licensed code on this list with production-quality composite ballistic
  material cards and element erosion already wired in. **The license
  caveat is the only reason this is not the obvious one-stop answer.**
- **Bibliography:** OpenRadioss does not have a single canonical paper
  (it is a software release rather than a publication); see Altair user
  documentation and Wikipedia entry for Radioss for the closed-source
  ancestor.

---

## 5. Cohesive / Fracture-Specific Codes

### 5.1 Akantu (already in primary survey, reaffirmed)

LGPL v3, gold standard for dynamic-cohesive-element insertion.
Vocialta-Richart-Molinari 2017 IJNME for tempered-glass fragmentation.
Kaliske 2024 / Richart 2024 JOSS for the v5.0 release.

### 5.2 MoFEM XFEM and cohesive

(Already in primary survey.) MIT core, LGPL extensions. Cohesive interface
elements documented; XFEM crack growth in nuclear-graphite case studies.
Implicit-quasi-static, not ballistic.

### 5.3 dolfinx_materials (Jeremy Bleyer, ENPC)

- **Repository:** github.com/bleyerj/dolfinx_materials
- **License:** **Creative Commons BY-SA 4.0** (verified — *content* license).
  *This is a documentation / examples license, not a code license.* The
  underlying dolfinx (LGPL v3) governs the linked solver.
- **Capabilities:** dolfinx Python add-on for non-UFL constitutive laws.
  Wraps MFront / MGIS to call user materials inside FEniCSx. Existing
  examples include J2-plasticity, gradient damage, phase-field.

### 5.4 MFront / TFEL + MGIS (Cyril Helfer, CEA)

- **Repository:** github.com/thelfer/tfel and github.com/thelfer/MFrontGenericInterfaceSupport
- **License:** **LGPL v3** for TFEL (the MFront code generator);
  **LGPL v3 / GPL v3 dual** for MGIS (the C++ interface library).
  Permissive enough that user materials can be embedded in proprietary
  codes.
- **Capabilities:** MFront generates user-material subroutines for Code_
  Aster, Cast3M, Abaqus / Calculix, Ansys, Europlexus, MFEM, FEniCS /
  FEniCSx (via MGIS), and JuliaFEM. Multiphase fiber-reinforced models
  exist in the MGIS-FEniCS examples (`mgis_fenics_multiphase_model`).
  No off-the-shelf Hashin / Puck / LaRC03 — these must be authored in
  MFront's `.mfront` DSL, but the framework is *the* canonical OSI-
  licensed bridge for constitutive-law portability.
- **Composite ballistic readiness:** **n/a** (it is a constitutive-law
  generator, not a solver) but **essential infrastructure** for any
  multi-solver pipeline.

### 5.5 newfrac/fenicsx-fracture (NewFrac ITN)

- **Repository:** github.com/newfrac/fenicsx-fracture
- **License:** **MIT** (verified).
- **Capabilities:** Phase-field, cohesive, gradient-damage examples in
  dolfinx. Educational set, not a production solver.

### 5.6 Cardamom (DLR aerospace)

**REJECTED — UNVERIFIED.** Repeated searches (DLR open-source portfolio,
Github, OSTI) did not return a public OSI-licensed package called
"Cardamom" relevant to composites or ballistic impact. There is a
*commercial* Cardamom (Bell Helicopter) and a research-grade `CARDAMOM`
high-fidelity gas-phase chemistry code. **Mark as not located** — either
the user's recall is of a closed-source tool or the project does not
exist as named under an OSI license. Excluded from final stack.

### 5.7 CalculiX cohesive zone (`*COHESIVE`)

CalculiX 2.20+ supports zero-thickness cohesive sections (built-in `CIP`
elements). GPL v2. Mature for implicit; explicit cohesive is documented
but rarely exercised in published composite ballistic.

---

## 6. Visualisation Beyond ParaView / PyVista

### 6.1 Mayavi (Enthought)

- **License:** **3-clause BSD** (Enthought-licensed).
- **Status:** v4.8.3 (2024). Tightly integrated with VTK; older API style;
  Jupyter integration via X3D.

### 6.2 K3D-Jupyter

- **Repository:** github.com/K3D-tools/K3D-jupyter
- **License:** **MIT** (verified).
- **Capabilities:** WebGL 3D plots in Jupyter. Cloud-points, isosurfaces,
  voxels, mesh, VTK, volume rendering.

### 6.3 Trame (Kitware)

- **Repository:** github.com/Kitware/trame
- **License:** **Apache 2.0** for trame core; **3-clause BSD** for trame-vtk.
- **Capabilities:** Pure-Python web-app framework for VTK / ParaView / vtk.js
  visualisations. Recent IEEE CG&A 2025 paper.

### 6.4 F3D

- **Repository:** github.com/f3d-app/f3d
- **License:** **3-clause BSD** (Kitware SAS + Migliore + Westphal +
  F3D-APP Foundation).
- **Capabilities:** Cross-platform desktop / CLI mesh-and-volume viewer.
  Supports a wide variety of formats (STL, OBJ, GLTF, VTK, Exodus, OpenVDB).

### 6.5 In-situ visualisation: Catalyst and Ascent

- **ParaView Catalyst:** **3-clause BSD**. Conduit-based API.
- **Ascent (LLNL):** **3-clause BSD**. Built on VTK-h / VTK-m. In-situ
  rendering for exascale simulations.
- **VTK-h:** **3-clause BSD**.

### 6.6 SENSEI

- **Repository:** github.com/SENSEI-insitu/SENSEI
- **License:** **3-clause BSD**.
- **Capabilities:** In-situ analysis abstraction layer; bridges
  Catalyst / Ascent / Libsim / ADIOS.

### 6.7 VTK-to-Typst-CeTZ workflow

For the user's Typst figure standard, the recommended chain is:

1. Run solver → write VTK / XDMF.
2. Read with `meshio` (BSD-3) or `pyvista` (MIT).
3. Slice / probe / sample on a regular grid using PyVista's `sample` and
   `clip` filters.
4. Export sampled fields to **CSV** or **JSON** with `numpy.savetxt` or
   `json.dumps`.
5. In Typst, import CSV via `csv("path.csv")` and draw with **CeTZ**
   per the user's preferences.

This pipeline is fully OSI-licensed end to end.

---

## 7. Composite-Specific Open Material-Model Libraries

| Model | Open implementations |
|---|---|
| Hashin 2D | CalculiX UMAT examples; Code_Aster forum scripts; community Abaqus UMATs (Khorashad / ammarkh95 — VUMAT, Abaqus-only); MFront `.mfront` recipe (community); **OpenRadioss `/FAIL/HASHIN`** (native) |
| Hashin 3D | Khorashad/3DHashin_VUMAT (Abaqus-only); MDPI 2024 enhanced PDM (Abaqus VUMAT); **OpenRadioss /FAIL/HASHIN** (3D mode) |
| Puck | Community Abaqus UMATs (closed-formulation); **OpenRadioss `/FAIL/PUCK`** (native); not in CalculiX, Code_Aster, Kratos out of the box |
| LaRC03 / LaRC04 / LaRC05 | Pinho et al. originals are paper-only; Abaqus UMAT community implementations exist; **PeriPy bond-based LaRC05-PD** (Hu 2023, paper-only at this writing); OpenRadioss does not ship LaRC05 natively |
| Tsai-Wu | CalculiX *failure index*; **OpenRadioss `/FAIL/TSAIWU`** (native); Code_Aster post-process |
| Christensen / Sun-Vaidya | No turnkey OSI implementation located |
| CDM (general orthotropic damage) | OOFEM (multi-model); Akantu (Mazars-orthotropic); MFront (user-coded); Kratos `SerialParallelRuleOfMixtures` |
| Cohesive (mode-I/II/III, bilinear / exponential) | CalculiX, Akantu, MoFEM, OpenRadioss `/MAT/COHESIVE`, FEniCSx (Bleyer comet tour), OOFEM |

**Bottom line:** if the user wants a Hashin / Puck / Tsai-Wu damage card
that *already exists, runs explicit, supports element erosion, and is
OSI-licensed*, the only candidate is **OpenRadioss (AGPL)**. Every other
OSI-licensed FEM code requires the user to author the failure card
themselves, typically via MFront / dolfinx_materials.

---

## 8. License Audit Summary

| Tool | License | Commercial OK? | Network clause? | Verdict |
|---|---|---|---|---|
| Peridigm | BSD-3 | Y | N | ACCEPT |
| PDLAMMPS / LAMMPS PERI | GPL-2 | Y (downstream GPL) | N | ACCEPT (copyleft caveat) |
| PeriPy | MIT | Y | N | ACCEPT |
| CabanaPD | BSD-3 | Y | N | ACCEPT |
| DualSPHysics | LGPL | Y | N | ACCEPT |
| PySPH | BSD-3 | Y | N | ACCEPT |
| SPlisHSPlasH | MIT | Y | N | ACCEPT (irrelevant) |
| CB-Geo MPM | MIT | Y | N | ACCEPT (irrelevant for composites) |
| Taichi Elements | MIT | Y | N | ACCEPT |
| AnisoMPM (ziran2020) | Apache-2 | Y | N | ACCEPT |
| NVIDIA Warp | Apache-2 | Y | N | ACCEPT |
| GeoTaichi | GPL-3 | Y (downstream GPL) | N | ACCEPT (copyleft caveat) |
| DiffTaichi | MIT | Y | N | ACCEPT |
| Matter (MPM) | UNVERIFIED | ? | ? | DEFER |
| TexGen | GPL-2 | Y (downstream GPL) | N | ACCEPT (copyleft caveat) |
| KinDrape | UNVERIFIED | ? | ? | DEFER (verify file) |
| fabrics-drape-data | UNVERIFIED | ? | ? | DEFER |
| WiseTex | restricted academic | N | N | **REJECT** |
| TexMind | commercial | N | N | **REJECT** |
| Yade | GPL-2 | Y (downstream GPL) | N | ACCEPT (copyleft caveat) |
| LIGGGHTS-PUBLIC | GPL-2+ | Y (downstream GPL) | N | ACCEPT |
| MercuryDPM | BSD-3 | Y | N | ACCEPT |
| MUSEN | BSD-3 | Y | N | ACCEPT |
| CFDEMcoupling | GPL-3 | Y (downstream GPL) | N | ACCEPT |
| PhasicFlow | UNVERIFIED | ? | ? | DEFER |
| Project Chrono | BSD-3 | Y | N | ACCEPT |
| JAX-FEM | GPL-3 | Y (commercial license needed for redist; contact author) | N | **CAUTION** |
| Firedrake | LGPL-3+ | Y | N | ACCEPT |
| scikit-fem | BSD-3 | Y | N | ACCEPT |
| FElupe | GPL-3 | Y (downstream GPL) | N | ACCEPT (copyleft caveat) |
| SfePy | BSD-3 | Y | N | ACCEPT |
| Brax | Apache-2 | Y | N | ACCEPT (irrelevant) |
| MJX | Apache-2 | Y | N | ACCEPT (irrelevant) |
| Genesis | Apache-2 | Y | N | ACCEPT (irrelevant for ballistic) |
| PolyFEM | MIT | Y | N | ACCEPT |
| IPC Toolkit | MIT | Y | N | ACCEPT |
| OpenSees / OpenSeesPy | restricted | **N** | N | **REJECT** (commercial redist needs license) |
| OpenRadioss | **AGPL-3** | Y | **YES** | **CONDITIONAL ACCEPT** (no SaaS hosting) |
| MFront / TFEL | LGPL-3 | Y | N | ACCEPT |
| MGIS | LGPL-3 / GPL-3 | Y | N | ACCEPT |
| dolfinx_materials | CC-BY-SA-4 (docs) + LGPL-3 (linked dolfinx) | Y | N | ACCEPT |
| newfrac/fenicsx-fracture | MIT | Y | N | ACCEPT |
| MFEM | BSD-3 | Y | N | ACCEPT |
| MOOSE | LGPL-2.1 | Y | N | ACCEPT |
| Tahoe | UNVERIFIED (Sandia open-source) | ? | ? | DEFER |
| Albany | BSD-3 (Sandia open-source) | Y | N | ACCEPT (heavyweight) |
| ExaConstit | BSD-3 (LLNL) | Y | N | ACCEPT (crystal plasticity, not composites) |
| FreeFEM | LGPL | Y | N | ACCEPT |
| Salome / Salome-Meca | LGPL-2.1+ | Y | N | ACCEPT |
| FreeCAD | LGPL-2 | Y | N | ACCEPT |
| Mayavi | BSD-3 | Y | N | ACCEPT |
| K3D-Jupyter | MIT | Y | N | ACCEPT |
| Trame core | Apache-2 | Y | N | ACCEPT |
| Trame-vtk | BSD-3 | Y | N | ACCEPT |
| F3D | BSD-3 | Y | N | ACCEPT |
| ParaView Catalyst | BSD-3 | Y | N | ACCEPT |
| Ascent | BSD-3 | Y | N | ACCEPT |
| SENSEI | BSD-3 | Y | N | ACCEPT |
| meshio | MIT | Y | N | ACCEPT |
| pyvista | MIT | Y | N | ACCEPT |
| Gmsh | GPL-2 | Y (downstream GPL) | N | ACCEPT (copyleft caveat) |
| pygmsh | GPL-3 | Y (downstream GPL) | N | ACCEPT |
| PyMesh | MPL-2 | Y | N | ACCEPT |

### Export-control notes

None of the *open-source* packages in this list are themselves under U.S.
ITAR or EAR restrictions, because publication on a public website is by
definition public-domain technical data under the EAR §734.7 / 734.8
"published" exclusion. **However**, applying any of these tools to model
defense articles (e.g., NIJ Level III+ ballistic body armour, hard armour
plates fielded by U.S. forces) can pull the *output data* (drawings, CAD,
material cards) into ITAR USML Cat XII — *the simulation tool is fine, the
results may not be*. This is on the user, not on the open-source code.

The MAT_162 model itself, for the LS-DYNA usage, was developed under U.S.
Army contract and the parameter sets for specific armour-grade composites
have been distributed under DoD distribution-statement controls (Distribution
B / C / E). **MAT_162 *parameter cards* may be export-controlled even when
the solver is not.** This is irrelevant if the user calibrates their own
ply card from publicly published coupon data.

---

## 9. Updated Single-Stack Recommendation

The previous survey recommended **FEniCSx for coupon + Kratos MPM for
ballistic**, with Akantu as a fragmentation-physics alternative. After the
wider survey, that recommendation **stands with three additions**.

### 9.1 The single most important addition: OpenRadioss

OpenRadioss did not appear in the original FEM-mainstream survey because
the survey was scoped to FEM libraries written for the academic research
ecosystem. OpenRadioss is an *industrial* solver — Altair's commercial
Radioss, dropped to AGPL v3 in 2022. Because the AGPL only triggers
on network-served distribution, and because the user's pipeline is
local-Python + cluster-batch, **OpenRadioss is the closest open-source
match to the LS-DYNA + MAT_162 workflow that has historically dominated
composite ballistic publications.** It ships:

- Native composite shell elements with arbitrary stack-ups.
- `MAT_LAW_25` (CRASURV: Tsai-Wu with damage), `MAT_LAW_15` (Chang-Chang),
  shell composite damage, solid composite damage.
- `/FAIL/HASHIN`, `/FAIL/TSAIWU`, `/FAIL/PUCK`, `/FAIL/JOHNSON` failure
  cards with element erosion *as native features*.
- Explicit dynamics, contact, impactor / projectile macros.
- Python keyword scripting since 2024 (`-python` mode).
- Cross-import of LS-DYNA `.k` and Abaqus `.inp` decks.

The cost is the AGPL clause: **the user must not host an OpenRadioss-as-a-
service web app** without publishing source. For an internal research
pipeline this is a non-issue.

### 9.2 The PeriPy / CabanaPD addition

Where the deliverable is *fragmentation physics* and not *bulk penetration
depth*, an open peridynamic code is the right answer:

- **PeriPy** if Python ergonomics and OpenCL-portable GPU dominate, and
  the simulation fits in single-node memory.
- **CabanaPD** if the simulation is HPC-cluster-class, multi-material,
  and the team is comfortable with C++ + Kokkos.

Both ship the multi-bond-type machinery to encode the AP-PLY fiber +
matrix bond structure that Hu 2023 and earlier composite-PD papers
require.

### 9.3 The MFront / MGIS infrastructure addition

For *every* solver in this stack, MFront + MGIS is the recommended way to
write a ply constitutive law once and run it inside Code_Aster, FEniCSx,
MFEM, Calculix, Cast3M, JuliaFEM, and (with effort) Akantu. **This is the
correct portability layer.** A Hashin / Puck / LaRC05 `.mfront` file is
the proper open-source artefact — not an Abaqus VUMAT.

### 9.4 The complete single-stack pipeline (revised)

| Stage | Primary tool | License | Backup |
|---|---|---|---|
| 1. Coupon tensile (D3039 / D3518) | FEniCSx + dolfinx_materials + MFront | LGPL + LGPL | OOFEM (LGPL); Code_Aster (GPL) |
| 2. Lamina / RVE homogenisation | FEniCSx periodic homogenisation; TexGen for textile RVE geometry | LGPL + GPL-2 | Kratos `SerialParallelRuleOfMixtures` (BSD); MFront periodic homogenisation |
| 3. AP-PLY tow geometry / draping | TexGen for woven primitives + custom Python tow-stacker; KinDrape for kinematic first-pass; MANUAL geometry from Lau 2023 / Nagaraj 2017 algorithms | GPL-2 + UNVERIFIED + custom | **Genuine open-source gap** for AFP-style path-planning |
| 4. Low-velocity drop-weight | OOFEM (Newmark + cohesive); FEniCSx-explicit (custom UL central difference); PolyFEM (IPC implicit) | LGPL / LGPL / MIT | OpenRadioss (AGPL) |
| 5. High-velocity ballistic + erosion | **OpenRadioss `MAT_LAW_25` + `/FAIL/HASHIN`** | AGPL-3 | Kratos MPM (BSD); Akantu cohesive (LGPL); Uintah MPM (MIT); CabanaPD (BSD); PeriPy (MIT) |
| 6. Fragmentation post-processing | meshio + pyvista → CSV/JSON → Typst+CeTZ | MIT + MIT + Apache-2 | Mayavi; F3D; Trame |
| 7. Constitutive-law portability | MFront + MGIS | LGPL-3 | dolfinx_materials |
| 8. In-situ viz (HPC runs) | ParaView Catalyst or Ascent | BSD-3 | SENSEI |

### 9.5 If "purely BSD/MIT, no AGPL, no GPL" is required

The user's no-restricted-license stance does *not* automatically forbid
AGPL or GPL — those are still OSI-approved free software. But if a
strict permissive-only stack is desired:

| Stage | Tool | Comment |
|---|---|---|
| 1-2 | FEniCSx (LGPL — borderline) | Or scikit-fem (BSD), MFEM (BSD), Albany (BSD), Kratos (BSD) |
| 3 | **GAP** | TexGen is GPL-2; no BSD draping code located |
| 4 | PolyFEM (MIT) + IPC Toolkit (MIT) | Implicit only |
| 5 | **Kratos MPM (BSD)** | Replaces OpenRadioss; performance and material-card maturity sacrificed |
| 6 | pyvista (MIT) + meshio (MIT) | OK |
| 7 | dolfinx_materials (CC-BY-SA + LGPL) | LGPL via dolfinx |
| 8 | Catalyst / Ascent (BSD-3) | OK |

In the strict-permissive case, **OpenRadioss is replaced by Kratos MPM**
at the cost of (i) losing pre-built `/FAIL/HASHIN` and (ii) longer
calibration of `SerialParallelRuleOfMixtures`. **TexGen has no BSD/MIT
substitute** — this is a real gap.

### 9.6 Honest gap declarations

The wider survey did not close the following gaps:

1. **No BSD/MIT-licensed AP-PLY-aware draping or AFP path-planning code
   exists.** TexGen is the only open-source textile-geometry tool of any
   maturity, and it is GPL-2. KinDrape is academic-minimal.
2. **No turnkey OSI-licensed composite ballistic solver matches LS-DYNA
   + MAT_162 in feature parity.** OpenRadioss is the closest, but it is
   AGPL. Kratos MPM is BSD but requires substantial constitutive-law
   work. Akantu is LGPL but requires recasting failure as cohesive.
3. **No production OSI implementation of LaRC03/04/05 was located.**
   Pinho-Camanho-Davila models are paper-only in the OSI ecosystem; the
   only public LaRC05 implementations are Abaqus VUMATs (closed-tool).
4. **DLR Cardamom (as named) was not located as an OSI-licensed package.**
   Either misremembered or non-public; flagged for follow-up by the user.
5. **WiseTex is restricted-academic and rejected by the user's policy.**
   No OSI substitute for KU Leuven's WiseTex meso-FE textile workflow.

These gaps are real, not papered-over. They are also opportunities — a
BSD-licensed AP-PLY tow generator and an MFront-authored LaRC05 are
reasonable open-source contributions in the user's project arc.

---

## Comparison Matrix (extended, sorted by composite-ballistic readiness)

Legend: `Y` = native, production; `~` = via user code or extension; `N` = absent.
"BR" = composite ballistic readiness 0-5.

| Tool | License | Coupon | Layup | Damage | Contact | Explicit | Erosion | Composite Lit. | BR |
|---|---|---|---|---|---|---|---|---|---|
| OpenRadioss | AGPL-3 | Y | Y | Y (Hashin/Tsai/Puck) | Y | Y | Y | strong (industrial) | **5** |
| Kratos MPM | BSD-3 | Y | Y (RoM) | Y | Y | Y | Y (MPM) | moderate | 4 |
| Akantu | LGPL-3 | Y | ~ | Y (CZM/PF) | Y | Y | Y (cohesive) | fragmentation | 4 |
| Uintah MPM | MIT | Y | ~ | Y (JC, CDM) | Y | Y | Y | metals strong | 4 |
| CabanaPD | BSD-3 | Y | ~ | Y (PD) | Y (DEM) | Y | Y (PD) | dynamic frac | 4 |
| PeriPy | MIT | Y | ~ (multi-bond) | Y (PD) | ~ | Y | Y (PD) | composite-PD lit | 4 |
| OOFEM | LGPL-2.1 | Y | Y | Y | Y | Y | ~ | weak | 3 |
| Peridigm | BSD-3 | Y | DIY | Y | Y | Y | Y (PD) | brittle frac | 3 |
| Karamelo | GPL-2 | ~ | DIY | DIY | Y (MPM) | Y | Y | weak | 3 |
| AnisoMPM | Apache-2 | Y | ~ | Y (aniso) | ~ | Y | Y | research | 3 |
| Code_Aster | GPL-3 | Y | Y | ~ (post) | Y | Y | N | weak | 3 |
| CalculiX | GPL-2 | Y | Y | ~ (UMAT) | Y | Y | ~ patch | weak | 3 |
| Kratos FEM | BSD-3 | Y | Y | Y | Y | Y | ~ | moderate | 3 |
| MoFEM | MIT | Y | ~ | ~ | Y | ~ | N | none | 2 |
| FEniCSx | LGPL-3 | Y | ~ | DIY | DIY | DIY | DIY | none for impact | 2 |
| Firedrake | LGPL-3 | Y | DIY | DIY | DIY | DIY | DIY | none | 2 |
| PolyFEM | MIT | Y | DIY | DIY | Y (IPC) | DIY | N | low-vel | 2 |
| JAX-FEM | GPL-3 | Y | DIY | DIY | DIY | DIY | DIY | inverse design | 2 |
| MFEM | BSD-3 | Y | DIY | DIY | DIY | Y | DIY | none | 2 |
| MOOSE | LGPL-2.1 | Y | DIY | Y (gen.) | Y | Y | DIY | nuclear | 2 |
| GetFEM | LGPL-3 | Y | ~ | DIY | Y | ~ | N | none | 2 |
| FElupe | GPL-3 | Y | ~ | DIY | DIY | N | N | rubber | 1 |
| scikit-fem | BSD-3 | Y | DIY | DIY | DIY | DIY | N | none | 1 |
| SfePy | BSD-3 | Y | DIY | DIY | DIY | ~ | N | none | 1 |
| Project Chrono | BSD-3 | Y | ~ | N | Y (MBD) | Y | ~ MBD | MBD only | 1 |
| Genesis | Apache-2 | N | N | N | Y | Y | N | robotics | 0 |
| Brax / MJX | Apache-2 | N | N | N | Y | Y | N | robotics | 0 |
| Yade | GPL-2 | N | N | N | Y (DEM) | Y | Y (DEM) | granular | 0 |
| LIGGGHTS-PUBLIC | GPL-2+ | N | N | N | Y (DEM) | Y | Y (DEM) | granular | 0 |
| OpenSees | restricted | Y | ~ | ~ | Y | Y | DIY | earthquake | **REJECTED** |
| WiseTex | restricted | textile | textile | N | N | N | N | textile | **REJECTED** |
| TexMind | commercial | textile | textile | N | N | N | N | textile | **REJECTED** |

---

## 10. References (added to `ecosystem_refs.bib`)

See `ecosystem_refs.bib` in this directory for new entries; existing keys
in `fenicsx_refs.bib`, `impact_refs.bib`, `test_progression_refs.bib`,
and `uofsc_refs.bib` are *not duplicated*.

Key additions cover Peridigm, PeriPy, CabanaPD, OpenRadioss, MFront /
MGIS, TexGen, KinDrape, AnisoMPM, NVIDIA Warp, Trame, and the IPC
Toolkit, among others.

---

## 11. Bottom Line (refresh)

- **One-code answer (composite ballistic, no GUI, ITAR-clean):
  OpenRadioss.** AGPL is the only friction. The user must accept the
  AGPL terms (no SaaS hosting) and the loss of the GPL-incompatible
  coupling that LS-DYNA users are familiar with.

- **Two-code answer (preferred): FEniCSx + MFront + MGIS for stages 1-3,
  OpenRadioss for stages 4-5, pyvista + Typst/CeTZ for stage 6.**
  This minimises license complexity, keeps the user inside the FEniCSx
  ecosystem they prefer, and delegates the explicit-impact problem to a
  code that was built for it.

- **Three-code permissive-only answer:** FEniCSx + MFront + Kratos MPM,
  with TexGen for textile geometry (accepting GPL-2 downstream
  contamination of any modified TexGen). Sacrifices feature parity vs
  OpenRadioss for license cleanliness.

- **Fragmentation-physics answer:** **PeriPy or CabanaPD** for the
  ballistic stage, replacing both OpenRadioss and Kratos MPM. Most
  rigorous fragmentation, smallest community, longest constitutive-law
  development effort.

- **Genuine open-source gaps:** AP-PLY-aware draping; production LaRC05
  implementation; OSI-licensed AFP path-planning. The user's project is
  positioned to *contribute* code into these gaps rather than consume it.
