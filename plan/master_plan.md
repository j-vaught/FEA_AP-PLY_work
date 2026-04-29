# FEA_AP-PLY — master plan (committed)

**Status.** Direction committed. Single open-source FEA package, solid elements throughout, no GUI. Pre-implementation; nothing has been built yet.

**Goal.** A pure open-source, no-GUI, code-only pipeline that simulates, from coupon up to ballistic,

1. tensile of a single composite tow,
2. tensile of a full layup with epoxy,
3. impact of composites,
4. preform-draping deformation of a pseudo-woven onto a tool (deferred side track),
5. and ultimately a high-velocity projectile penetrating a fully cured UofSC pseudo-woven panel.

---

## 1. Direction — single tool, all solid elements

**OpenRadioss** (AGPL-3, Altair, dropped to OSI in Sept 2022) for every stage. Driven entirely from Python text-templating and the command-line `starter` / `engine` runner. Verified end-to-end by `references/openradioss_endtoend_audit.md`.

| Concern | Choice |
|---|---|
| Solver | OpenRadioss (implicit static + explicit dynamic in one binary) |
| Element family | Solid only — HEXA8 / TETRA10 / pentahedra; no shell layups, no laminate shell elements |
| Composite damage | LAW25 + /FAIL/HASHIN, /FAIL/PUCK, /FAIL/TSAIWU on /PROP/TYPE14 solid bricks |
| Cohesive zones | /INTER/TYPE2 surface cohesive + solid cohesive elements (DCB / ENF) |
| Element erosion | Native, on solid composite cards |
| Implicit static | Native, requires MUMPS-linked build |
| Multi-step (LVI → CAI) | Native explicit-explicit state-mapping (state file restart) |
| Custom material | userlib_sdk (Fortran) — path to port Kok's VUMAT if native cards prove insufficient |
| Geometry | `rutger-kok/ap_ply_model_creation` (LGPL-2.1) → GMSH `.msh` → Abaqus `.inp` → OpenRadioss via `inp2rad` |
| Visualization | OpenRadioss native output → Kitware `openradioss-to-vtkhdf` → ParaView `pvbatch` headless or PyVista headless |
| Plots | CSV export → Typst + CeTZ |
| Install (macOS) | Lima + Apptainer Linux ARM64 image (no native macOS build); Linux is bare-metal supported |

**Why single-tool OpenRadioss won out.** The end-to-end audit confirmed OpenRadioss covers every stage of the test progression except generic 3D periodic BCs for RVE homogenization. We replaced those two stages with direct mesoscale simulation (see §3), preserving the verification spirit of the scaffold without a second tool.

**Rejected and why.**
- FEniCSx — no production contact, no element erosion, no off-the-shelf Hashin/Puck. Forced a hybrid that the user rejected.
- Code_Aster — weak explicit dynamics; no element erosion; no native Hashin/Puck.
- CalculiX — limited cohesive zones, weaker composite damage cards.
- Akantu, Kratos, OOFEM — composite ply law would need to be user-coded for ballistic.
- Abaqus / LS-DYNA / EUROPLEXUS / Ansys / PAM-FORM — closed source.

**Single honest gap kept on the record.** OpenRadioss has no generic 3D periodic-BC keyword (`/BCS/CYCLIC` is cylindrical, `/RBE2/3` and `/MPC` are insufficient for periodic homogenization). Stages 10-11 are reframed accordingly.

---

## 2. Architecture identity — one architecture, three names

UofSC pseudo-woven (PW) / meso-architectured composite (MAC) is the same architecture as TU Delft AP-PLY. UofSC also called it "clutch laminate" (Van Tooren and Sockalingam 2017). Z. Gurdal is the human bridge. Confirmed against Nagelsmit 2013 PhD, Vakili Rad 2020 MSc, Kodagali 2023 PhD, Zheng 2016 PhD, Kok 2022/2023/2024.

**Locked geometry parameters** (Kok `ap_ply_model_creation` defaults plus UofSC theses).

| Parameter | Value | Source |
|---|---|---|
| Slit-tape width | 6.35 / 10 / 12.7 mm | Kok defaults |
| Cured ply thickness | 0.18 / 0.18125 / 0.213 mm | Kok defaults |
| Undulation ratio | 0.09 - 0.104 | Kok defaults |
| Tape spacing s | 1, 2, 3 (one-over-one to three-over-three) | Kok defaults |
| Reference material | SHD VTC401 carbon/epoxy | Kok defaults |
| UofSC HVI material | Hexcel IM7G / 8552-1 | Vakili Rad 2020 |
| UofSC LVI/CAI material | Toray T800-SC-24K / P2362W | Kodagali 2023 |
| AFP machine (UofSC) | "LYNX" gantry, 8 tow channels | Vakili Rad 2020 |
| Cure cycle | 177 °C, 0.62 MPa, 6 h | Kodagali 2023 |
| Reference panel | 1143 × 838 mm, 24 ply, 4.55 mm | Kodagali 2023 |
| Naming convention | `[fiber angles][placement seq.][angle shift][tow width]` | Kodagali 2023 |

---

## 3. Test progression (16 stages, single-tool OpenRadioss, solid elements)

Each stage isolates one capability and has a quantitative pass criterion. V/E marks Verification (closed form) or Validation (experiment). The only stages not in the original sketch are 10 and 11, which are reframed from RVE periodic homogenization to direct mesoscale simulation (justification at the bottom of this section).

| # | Stage | Standard | OpenRadioss feature(s) | Pass criterion | V/E |
|---|---|---|---|---|---|
| 1 | Linear-elastic 3-pt and 4-pt bending of a solid beam | — | /MAT/LAW1, /PROP/TYPE14 solid HEXA, /IMPL static | tip deflection within 1% of Euler-Bernoulli | V |
| 2 | Cantilever, large-deflection geometric nonlinearity | — | /IMPL static, geom. nonlin. on | within 2% of Bisshopp-Drucker 1945 | V |
| 3 | Isotropic dogbone tension | ASTM E8 | /MAT/LAW2 or /MAT/LAW36, /IMPL static | uniform σ within 1%, ε within 1% | V |
| 4 | Open-hole tension (Kirsch) | ASTM D5766 | /MAT/LAW1, mesh refinement at hole | $K_t$ within 2% at hole edge | V |
| 5 | Dogbone with isotropic ductile damage | — | /MAT/LAW22 or /MAT/LAW23 (Mazars substitute) | mesh-objective load-deflection within 5% RMS | E |
| 6 | Dogbone with Tsai-Wu, Hashin, Puck side-by-side | — | /MAT/LAW25 + /FAIL/TSAIWU, /FAIL/HASHIN, /FAIL/PUCK | failure envelope matches WWFE-II | E |
| 7 | UD tow tension on-axis and off-axis | ASTM D3039 | /MAT/LAW25 orthotropic on solid HEXA | $E_1, E_2, \nu_{12}, G_{12}$ within 2% of Soden card | V |
| 8 | Ply-rotation transformation | — | per-element material orientation | engineering constants vs. analytic within 1% | V |
| 9 | Cross-ply and quasi-iso laminate (solid plies) | — | stack of solid HEXA plies, one or more elements per ply | stiffness within 2% of CLT reference | V |
| 10 | UD tow-wise direct mesoscale (reframed) | — | discrete fiber / matrix solid mesh, displacement BCs | effective $E_1, E_2, G_{12}$ within 5% of Halpin-Tsai | E |
| 11 | Pseudo-woven mesoscale direct (tow + matrix, reframed) | — | Kok `ap_ply_model_creation` → inp2rad, displacement BCs | homogenized stiffness within 10% of Kok 2022 reported values | E |
| 12 | DCB and ENF cohesive zone | ASTM D5528, D7905 | /INTER/TYPE2 + solid cohesive elements | DCB peak load within 5% | E |
| 13 | Drop-weight LVI | ASTM D7136 | explicit dynamic, contact, /FAIL/HASHIN with erosion | peak force within 10%, delam area within 20% | E |
| 14 | Compression after impact (CAI) | ASTM D7137 | explicit-explicit chain via state restart from stage 13 | residual strength within 10% | E |
| 15 | Flat-coupon ballistic | NIJ 0101.06 / MIL-STD-662F | explicit dynamic + composite erosion | V50 within 5%, residual velocity within 10% | E |
| 16 | Pseudo-woven panel ballistic | ASTM D8101 | Kok geometry + /FAIL/HASHIN + erosion + projectile contact | V50 within 7% of Vakili Rad 2020 | E |

**Why stages 10-11 are reframed.** Periodic boundary conditions for RVE homogenization are an OpenRadioss documentation gap (audit cross-cut F). The physical question — "does the mesoscale architecture recover the right macroscale stiffness?" — is answered equivalently by direct simulation under uniform displacement BCs on a sufficiently large mesoscale specimen, comparing effective stiffness to Halpin-Tsai (stage 10) and Kok 2022 reported values (stage 11). This preserves the scaffold's verification function while staying inside one tool. The user's ultimate goal is tow-wise direct simulation of the panel anyway (stage 16); periodic homogenization was always a verification sidecar, not on the critical path.

**Forming/draping side track** (preform fitting a tool) is deferred. No BSD/MIT-licensed AP-PLY draping code exists; it would be a separate research track in TexGen / KinDrape (both copyleft / academic-minimal).

---

## 4. Implementation pipeline (one direction, no handoff)

```
                                           ┌────────────────────────────────────┐
[Constituent + ply + cure data]            │  Python deck templating            │
        │                                  │  (jinja2 / f-strings on .rad text) │
        ▼                                  └────────────────┬───────────────────┘
[Geometry: Kok ap_ply_model_creation]                       │
        │                                                   │
        ▼                                                   ▼
[GMSH mesh (.msh)] ─── meshio ───► [Abaqus .inp] ── inp2rad ──► [OpenRadioss .rad deck]
                                                                    │
                                                                    ▼
                                                 OpenRadioss `starter` + `engine`
                                                  (implicit + explicit, one binary,
                                                   driven from CLI / Python)
                                                                    │
                                                                    ▼
                                            OpenRadioss `.anim`, `T01`, `.h3d`
                                                                    │
                                                                    ▼
                                       Kitware `openradioss-to-vtkhdf` converter
                                                                    │
                                                                    ▼
                                           ParaView `pvbatch` headless OR PyVista
                                                                    │
                                                                    ▼
                                              PNG / MP4 + CSV (numerical export)
                                                                    │
                                                                    ▼
                                                 Typst + CeTZ vector plots
```

No GUI anywhere. Each box is text-scriptable.

---

## 5. Visualization

In-script diagnostic frames: PyVista headless, written from `.vtkhdf` straight to PNG / MP4.

Publication frames: ParaView `pvbatch` Python scripts, also reading `.vtkhdf`.

Quantitative plots: every numerical result (load-displacement, V50 vs. layup, residual velocity vs. impact velocity) exported as CSV; plots authored in Typst + CeTZ per global preferences.

---

## 6. Reference inventory

All under `references/`:

| File | Contents |
|---|---|
| `openradioss_endtoend_audit.md` + `openradioss_refs.bib` | End-to-end stage-by-stage capability proof; 33 entries |
| `fenicsx_capability_audit.md` + `fenicsx_refs.bib` | What FEniCSx cannot do — kept for the rejection record |
| `uofsc_pseudowoven.md` + `uofsc_refs.bib` | UofSC PW lineage; 17 entries |
| `delft_apply.md` + `delft_refs.bib` | Delft AP-PLY origin; 28 entries |
| `fea_alternatives.md` + `impact_refs.bib` | Mainstream OS FEA alternatives — kept for the rejection record; 22 entries |
| `test_progression_literature.md` + `test_progression_refs.bib` | 16-stage canonical benchmarks |
| `open_source_ecosystem_survey.md` + `ecosystem_refs.bib` | Wider OS ecosystem, peridynamics fallbacks; 35+ entries |

No duplicate keys across the seven `.bib` files.

**Open-data gap (still real)**: nothing on 4TU.ResearchData or Zenodo; no public micro-CT; no public AFP G-code. Ground truth must be digitized from Nagelsmit 2013, Vakili Rad 2020, Kodagali 2023, Zheng 2016/2024, and Kok 2022/2023/2024.

---

## 7. Toolchain install (macOS user)

```bash
# 1. Lima + Apptainer (OpenRadioss has no native macOS build).
brew install lima
limactl start --name=or template://default
limactl shell or sudo apt update && sudo apt install -y apptainer build-essential gfortran

# 2. OpenRadioss prebuilt binaries (or compile from source with MUMPS link).
# https://openradioss.atlassian.net/wiki/spaces/OPRD/pages/47546369/Build+OpenRadioss

# 3. Geometry pipeline (LGPL-2.1).
git clone https://github.com/rutger-kok/ap_ply_model_creation.git
git clone https://github.com/rutger-kok/composite_cdm_ap_ply.git   # reference VUMAT
pip install gmsh meshio numpy scipy

# 4. Mesh bridge.
# Use meshio + Abaqus .inp + OpenRadioss inp2rad converter (ships with OpenRadioss).

# 5. Visualization (headless).
pip install pyvista[all]
# ParaView (headless): brew install --cask paraview     # macOS native pvbatch is fine; renders the .vtkhdf produced inside Lima
# Kitware openradioss-to-vtkhdf: clone from kitware repo (cited in audit).

# 6. Plots.
# Typst + CeTZ already part of user's standard toolchain.
```

The Lima + Apptainer step is the only real friction; once running, OpenRadioss inside Lima reads / writes files on the macOS filesystem via the shared mount.

---

## 8. Phase 2 — per-test deep-dive agents (next move)

Once this plan is approved, spawn 16 opus agents in parallel, one per stage. Each produces a spec at `tests/stage_NN/spec.md` covering:

- Geometry (numerical dimensions + GMSH script).
- Mesh (element type, target size, refinement zones).
- Boundary conditions and loading.
- Material card (numerical values, citation source).
- OpenRadioss deck skeleton (specific keyword cards listed).
- Reference solution (closed-form, paper figure, or experimental dataset).
- Validation success criterion (the line in the table above).
- Per-stage runner script (Python that templates the .rad deck and invokes starter + engine inside Lima).

---

## 9. Phase 3 — first runnable

Stage 1 (linear-elastic solid beam, OpenRadioss implicit static, closed-form check). End-to-end smoke test of the toolchain: GMSH → meshio → inp2rad → OpenRadioss (Lima) → vtkhdf → PyVista PNG → CSV → Typst CeTZ figure. If this fails the toolchain is broken; if it passes, stages 2-9 are mostly variations.

---

## 10. Risks and unknowns (revised)

1. **macOS install friction** — Lima + Apptainer is well-documented but adds one layer between IDE and solver. File I/O performance crosses the VM boundary.
2. **MUMPS-linked OpenRadioss build** — implicit static (stages 1-2, 6, 12) requires linking MUMPS at build time; the prebuilt binaries may be explicit-only. Either compile from source or accept dynamic relaxation as a workaround.
3. **No public AP-PLY raw data** — V50 ground truth is one number per Vakili Rad 2020 panel; no per-shot residual-velocity table. Stage 15 / 16 calibration is bounded by what the theses publish.
4. **MFront ↔ OpenRadioss bridge does not exist** — userlib_sdk (Fortran) is the only custom-material path. Reasonable, but means any custom CDM is Fortran-77, not Python or UFL.
5. **Periodic-BC gap is closed by reframing** — but if a reviewer challenges the equivalence between direct mesoscale and periodic homogenization, the response is documented in §3 (uniform-displacement BCs on a sufficiently large mesoscale specimen recover the same effective stiffness for the ply-scale property; this is a documented technique in the homogenization literature).
6. **Kitware vtkhdf converter is recent** — confirm the converter version against the OpenRadioss output format we will be using; mismatch is plausible.
7. **AGPL-3 contamination** — only matters if we host OpenRadioss as a service. Default plan: distribute only `.rad` decks and Python templating scripts, not OpenRadioss derivatives.

---

## 11. Phase 2.6 — Kok geometry port (committed)

The Kok preprocessor (`rutger-kok/ap_ply_model_creation`, LGPL-2.1) is an Abaqus CAE Python plugin and requires Abaqus to run. This violates the open-source-only constraint at the geometry-generation step for stages 11 and 16. The committed resolution is to port the AP-PLY laydown logic to a standalone Python module backed by `gmsh.model.occ.*` (OpenCASCADE via GMSH).

**Approach.** Clean-room reimplementation from the public algorithmic descriptions in Nagelsmit 2013 PhD Chapter 2 (§2.2 Concept, §2.3 Pattern Details) and Kok 2022 Composites Part A §3.1 Microscale model. Section pointers verified by direct-PDF read in the port-plan agent run; the earlier "§3.4" / "§2" pointers in this document were wrong and have been corrected. Kok's GitHub repo is acknowledged as prior art in the port's README but not read or transcribed during implementation, so the port can ship under MIT/BSD rather than inheriting LGPL.

**Output**: a standalone CLI `kok_geom --config kok_config.json --out panel.{msh,inp}` that consumes the same JSON schema stages 11 and 16 expect, produces a meshable solid model with per-tow physical groups, and writes a JSON sidecar with per-element-set fiber orientation vectors.

**Module location**: `/Volumes/MacShare/Code/FEA_AP-PLY/src/kok_geom/`.

**Validation**: stage 11 must reproduce Kok 2022 reported homogenized stiffness within the 10% tolerance already in the master plan. No direct mesh-vs-Kok comparison is needed (and would require Abaqus); the homogenized stiffness check is the validation criterion.

**Scope**: clean-room AP-PLY tow-laydown for the master-plan's locked geometry parameters (slit-tape 6.35 / 10 / 12.7 mm, ply thickness 0.18 / 0.18125 / 0.213 mm, undulation 0.09-0.104, tape spacing 1/2/3, the Kok 2022 baseline stacking sequences). Out-of-scope: AFP machine path optimization, defect modeling, multi-axial preform draping (these are separate research tracks).

A planning agent will produce `plan/kok_port_plan.md` next, with module breakdown, milestone gates, and validation strategy. Implementation is scheduled after Phase 3 (stage 1 toolchain smoke test) so we have a verified OpenRadioss + Lima pipeline before the port lands.

---

## 12. Sign-off requested

§1 (single tool, all solid), §3 (16 stages including reframed 10-11), §4 (single direction, no handoff), §11 (Kok port committed, clean-room from papers, Python + gmsh-OCC, MIT/BSD).
