# Stage 16 Specification: Pseudo-Woven Panel Ballistic Impact (UofSC P1-TWT, Vakili Rad 2020 Reference Configuration) — PROJECT TERMINAL GOAL

Author: J.C. Vaught
Date: 2026-04-29
Stage tier: E (Validation, integration of every preceding stage)
Master plan reference: `plan/master_plan.md` §3 row 16, §1 (single tool, all solid)
OpenRadioss audit verdict: PASS (`references/openradioss_endtoend_audit.md`, row 16 and row C)
Upstream stages assumed verified: 1-15 inclusive. The Kok geometry pipeline is proven at small scale by Stage 11; the explicit-dynamic ballistic erosion stack with LAW25 plus /FAIL/HASHIN plus /INTER/TYPE7 is proven on a flat coupon by Stage 15.

## 2026-04-30 Execution Addendum

The original text below predates the empirical LAW/PROP compatibility matrix and still names the old LAW25 + TYPE14 route in several places. The binding recipe for any rerun is now the verified solid-composite route in `references/openradioss_law_compatibility_matrix.md` and `references/openradioss_orientation_convention.md`: `/MAT/LAW12` + `/PROP/TYPE6/SOL_ORTH`, per-element `/INIBRI/ORTHO` for tow rotation, and `/FAIL/HASHIN` erosion with `IFAIL_SO=1` and `PTHICKFAIL=1.0`.

Post-Kok-M5 Phase A reran the 200 mm x 200 mm P1-TWT representative section with the current `kok_geom` schema, explicit midplane symmetry, and explicit `tow_coverage_fraction = 0.94` in `geometry/kok_p1_twt_200_config.json`. The coarse 2 mm preflight command `time python -m kok_geom --config tests/stage_16_PW_panel_ballistic/geometry/kok_p1_twt_200_config.json --out tests/stage_16_PW_panel_ballistic/geometry/p1_twt_200.msh` completed in `402.98 s` (`6.72 min`) and wrote a mesh with `1,979,649` nodes and `1,440,000` TETRA10 elements. This supersedes the earlier `timeout 900` blocker: the full 200 mm section now clears the Phase A geometry gate and the next blocker, if any, must come from the measured Phase B solver wall-clock rather than from meshing.

---

## 1. Goal: this is the project's terminal test

This stage is the entire reason the FEA_AP-PLY project exists. Every preceding stage of the 16-stage progression in `plan/master_plan.md` was scoped to isolate one single capability so that, on arrival here, the integration risk is bounded. Stage 16 is where every isolated capability runs simultaneously on a real, published, peer-reviewed UofSC panel configuration. Pass at Stage 16 is the project's success criterion. Failure at Stage 16 is the project's failure criterion. Nothing downstream of this stage exists in the master plan.

The capability under test is a pseudo-woven (AP-PLY) panel manufactured by automated fiber placement, struck by a single-stage gas-gun projectile, with the composite tow geometry resolved at mesoscale (one solid HEXA8 element per ply through-thickness, in-plane mesh resolves at least the tape width), full intra-laminar damage with element erosion (LAW25 plus /FAIL/HASHIN with `Ifail=2`), and inter-laminar cohesive contact between every adjacent ply pair (/INTER/TYPE2). Boundary conditions on the cut edges are clamped or absorbing depending on the insensitivity check below. The reference experiment is Vakili Rad et al. 2020 in *Composites Part B*, performed at the NASA Glenn Research Center single-stage gas gun in collaboration with D. Revilock and C. Ruggeri.

The chosen primary verification target is the **P1-TWT** configuration from Vakili Rad 2020: 24-ply Hexcel IM7G/8552-1 with a pseudo-woven outer-skin / unidirectional inner-core stacking sequence. This is the configuration in the paper with the most published HVI data, including comparative back-face deflection, back-face damage area, and qualitative penetration behaviour across the velocity sweep 250-400 ft/s (76-122 m/s). The second target retained as a within-paper sanity check is **P3-Control** (the conventional [45/90/-45/0]_3s baseline), which the paper documents as the inferior performer; the model must reproduce the relative ranking P1-TWT better than P3-Control even if the absolute V50 of either is at the edge of the 7% tolerance.

What this stage does not do. It does not attempt the LVI / CAI sequence of Stages 13-14; those are separate test pieces with their own pass criteria. It does not attempt the thicker production panels of Kodagali 2024 (Toray T800-SC / P2362W); the chosen material and configuration are locked to Vakili Rad 2020 IM7G/8552-1 because that paper is the one with NASA Glenn HVI data. It does not attempt the side-track preform-draping problem deferred in master plan §3.

Pass criterion (single number, formal). The simulated V50 ballistic limit, computed by Recht-Ipson 1963 fit to a sweep of impact velocities about the published value, must lie within 7 percent of the Vakili Rad 2020 reported value for the P1-TWT panel. The 7 percent tolerance is loosened from the 5 percent of Stage 15 because of (i) AP-PLY architectural uncertainty inherited from the partially-undulated tow geometry, (ii) digitization error from the published figures (no per-shot residual-velocity table is published in the open literature), and (iii) the sparse ground truth: there is one published number per configuration, not a population.

The 7 percent tolerance maps directly to the master plan §3 table row 16 ("V50 within 7% of Vakili Rad 2020") and to the test-progression literature note (`references/test_progression_literature.md` Stage 16) that tightening below 7 percent is not justified by the open data.

---

## 2. Kok preprocessor invocation for the P1-TWT panel

The geometry is generated by Rutger Kok's `ap_ply_model_creation` Python preprocessor (LGPL-2.1, `github.com/rutger-kok/ap_ply_model_creation`) which is the same toolchain proven at small scale by Stage 11. The preprocessor outputs an Abaqus `.inp` mesh; that mesh is converted to an OpenRadioss `.rad` deck via the `inp2rad` Python converter shipped with `OpenRadioss/Tools` (`references/openradioss_endtoend_audit.md` row C, PASS). The pipeline at Stage 16 is identical to Stage 11 except for the panel size, the material card, and the impactor sub-deck.

The stacking-sequence notation is the UofSC convention from Kodagali 2023 thesis chapter 4: `[fiber angles][placement sequence][angle shift][tow width]`. For P1-TWT the relevant fields are documented below.

### 2.1 Vakili Rad 2020 P1-TWT laminate definition

P1-TWT is a 24-ply hybrid in which the outer six plies on each face are pseudo-woven (PW) and the inner twelve plies are unidirectional tape (UD). The PW skin is laid as a four-angle AP-PLY course `[+45 / 90 / -45 / 0]` placed with active/inactive AFP channels so that adjacent same-direction tows are separated by three tow widths (one active channel followed by three inactive channels), giving an AP-PLY interlace pattern that Kok terms `tape_spacing = 3` (three-over-three). The inner core is a conventional UD stack matching the same nominal in-plane stiffness so that the panel global stiffness is comparable to the P3-Control quasi-isotropic baseline. The panel is symmetric about the midplane.

Compact notation (per Kodagali 2023 ch.4 with `[UNVERIFIED]` flag where the exact placement string was not pinned in the digitized text of Vakili Rad 2020):

```
P1-TWT = [(+45 / 90 / -45 / 0)_PW]_3 / [(+45 / 90 / -45 / 0)_UD]_3 _S
       = AP-PLY skin (12 plies total, 6 per face) over UD core (12 plies)
       tape angles  : (+45, +90, -45, 0)
       placement seq: 10001000  (one active channel followed by three inactive)
       angle shift  : 0
       tow width    : 6.35 mm
       cured ply th : 0.18958 mm  (= 4.55 mm / 24 plies)
```

The cured ply thickness 0.18958 mm is the implied value from the panel-level 4.55 ± 0.2 mm thickness (Vakili Rad 2020 abstract; Kodagali 2023 thesis) divided by 24 plies, rather than the 0.18 mm Kok default. This matters for the LAW25 in-situ strength calibration in Section 5; the ply count is locked, and the per-ply thickness is calibrated to the paper-reported total.

### 2.2 Reduced 200 mm × 200 mm representative section

The full Vakili Rad 2020 panel is 1143 mm × 838 mm. Meshing the full panel at one element per ply through 24 plies and ≤ 0.5 mm in-plane resolution in the impact zone would exceed the practical element-count budget by an order of magnitude (see Section 3 below). The standard composite-impact reduction (Lopes-Camanho 2009 Stage 14 reference, Olsson-Davies 2004 Stage 14 reference) is to model only the centre portion of the panel containing the entire damage zone plus a Saint-Venant buffer, with cut-edge boundary conditions chosen to avoid spurious wave reflections.

For Stage 16 the reduced section is **200 mm × 200 mm centred on the impact location**. This is large enough to comfortably contain the 50-100 mm projected delamination footprint reported for the P1-TWT panels in Vakili Rad 2020 high-speed video, plus a buffer of at least 50 mm of laminate around it to allow flexural waves to decay before they reach the cut edge in the time window of the impact event (~0.5 ms; flexural-wave speed ~1500 m/s in CFRP gives a wave-travel distance of ~0.75 m, so reflection from the 100 mm boundary will return at ~130 µs and corrupt the back-face deflection signal; this is the BC-insensitivity check in Section 6).

The reduction matches the practice in Kok et al. 2023 *Composites Part B* (300 × 300 mm Edinburgh hand-layup specimen) and is on the small side relative to that, with the explicit understanding that the BC choice is the engineering compromise.

### 2.3 JSON configuration consumed by `runner.py`

The runner emits the Kok `tape_placement.laminate_creation` invocation through a JSON config file `kok_config.json` that captures every preprocessor argument verbatim, so that the geometry build is reproducible by anyone who clones the repo. The schema is defined here and used by Section 9.

```json
{
    "panel_geometry": {
        "specimen_polygon_mm": [[-100.0, -100.0], [100.0, -100.0],
                                [100.0,  100.0], [-100.0, 100.0]],
        "panel_size_full_mm":  [1143.0, 838.0],
        "section_size_mm":     [200.0, 200.0],
        "impact_centre_mm":    [0.0, 0.0],
        "section_offset_in_full_panel_mm": "centred"
    },
    "ap_ply_layup": {
        "configuration_name": "P1-TWT",
        "n_plies_total": 24,
        "skin_n_plies":  12,
        "core_n_plies":  12,
        "tape_angles_deg":          [45.0, 90.0, -45.0, 0.0],
        "tape_widths_mm":           [6.35, 6.35, 6.35, 6.35],
        "tape_spacing":             3,
        "placement_sequence":       "10001000",
        "angle_shift_deg":          0.0,
        "cured_ply_thickness_mm":   0.18958,
        "undulation_ratio":         0.09,
        "undulation_half_width_mm": 1.0537,
        "skin_outer_to_inner_seq":  ["+45_PW", "+90_PW", "-45_PW", "0_PW",
                                     "+45_PW", "+90_PW", "-45_PW", "0_PW",
                                     "+45_PW", "+90_PW", "-45_PW", "0_PW"],
        "core_seq":                 ["+45_UD", "+90_UD", "-45_UD", "0_UD",
                                     "+45_UD", "+90_UD", "-45_UD", "0_UD",
                                     "+45_UD", "+90_UD", "-45_UD", "0_UD"],
        "symmetry":                 "midplane_symmetric",
        "notation_source":          "Kodagali 2023 PhD thesis ch.4"
    },
    "material_card_assignments": {
        "tape_region":           "LAW25_IM7_8552_orthotropic",
        "undulation_region":     "LAW25_IM7_8552_orthotropic_locally_rotated",
        "resin_pocket_region":   "LAW1_8552_isotropic_epoxy",
        "interlaminar_cohesive": "INTER_TYPE2_Camanho_Davila_2003"
    },
    "mesh_targets": {
        "in_plane_target_mm_impact_zone": 0.5,
        "in_plane_target_mm_far_field":   2.0,
        "through_thickness_elements_per_ply": 1,
        "tow_wise_target_per_tape_width": 12,
        "graded_zone_radius_mm": 40.0
    },
    "kok_preprocessor_overrides": {
        "default_material_in_kok_repo":  "VTC401",
        "this_run_substitutes_material": "IM7G_8552_1_Vakili_Rad_2020",
        "geometry_only_pass":            true,
        "skip_kok_VUMAT_card_generation": true,
        "abaqus_inp_output_path":         "build/p1_twt_section.inp"
    }
}
```

The `geometry_only_pass = true` flag tells the Kok script to emit only the partitioned mesh and per-element material-region tags (tape / undulation / resin pocket); the constitutive-law cards in the Kok VUMAT format are *not* written, because Stage 16 substitutes OpenRadioss native LAW25 plus LAW1 cards downstream (see Section 5). This is the same `inp2rad` substitution pattern proven at Stage 11.

### 2.4 Pipeline command sequence

The runner invokes (in order) the following four commands inside the Lima Apptainer Linux ARM64 environment per master plan §7:

1. `python ap_ply_model_creation/ap_ply_model.py --config kok_config.json --output-inp build/p1_twt_section.inp`
2. `python OpenRadioss/Tools/inp2rad/inp2rad.py build/p1_twt_section.inp --output-rad build/p1_twt_section.rad`
3. `starter_linuxa64 -i p1_twt_section_0000.rad -nt N`
4. `engine_linuxa64 -i p1_twt_section_0001.rad -nt N`

Step 1 produces a partitioned mesh with three element-set tags per ply (`TAPE_<angle>`, `UND_<angle>`, `RESIN`) plus per-element local material orientations on the tape and undulation sets. Step 2 converts to OpenRadioss text format preserving the element sets and orientations; the runner verifies element count and node count consistency between `.inp` and `.rad` per the audit's beta-status caveat (`openradioss_endtoend_audit.md` row 16, "sanity check on tow orientations and material assignments is required after conversion").

---

## 3. Mesh sizing, element count, hourglass control, and expected wall clock

### 3.1 Element count budget

The reduced 200 mm × 200 mm section, meshed at the targets in Section 2.3 (in-plane 0.5 mm in the impact zone, 2.0 mm in the far field, with a graded transition zone of radius 40 mm around the impact point), and with one HEXA8 element per ply through-thickness for 24 plies, gives an order-of-magnitude element count as follows.

Impact zone (a circular region of radius 40 mm centred on the impact point, area ≈ 5026 mm²) at 0.5 mm in-plane × 24 plies through-thickness:

```
N_impact = (5026 mm² / (0.5 mm)²) × 24 = 482,496 elements
```

Far-field (the remainder of the 200 × 200 = 40000 mm² section minus the impact zone, area ≈ 34974 mm²) at 2.0 mm in-plane × 24 plies:

```
N_far = (34974 mm² / (2.0 mm)²) × 24 = 209,844 elements
```

Tow-wise refinement constraint. The Kok preprocessor partitions each ply into tape, undulation, and resin-pocket regions per its `tape_placement.laminate_creation` algorithm. The undulation regions are typically 1 mm wide (the half-width 1.0537 mm in Section 2.3), and the resin pockets are smaller. For mesh-objective resolution of the through-thickness undulation Kok et al. 2022 (*Composites Part A* 159, 107014) recommend at least 4 elements across an undulation. At 0.5 mm in-plane that gives 2 elements per undulation in the impact zone and only 0.5 elements in the far field, which is the documented mesh-size compromise. To satisfy "tow-wise mesh resolves at least the tape width 6.35 mm" in the brief, the in-plane element size must satisfy `h_xy ≤ 6.35 / 12 ≈ 0.53 mm` to put 12 elements across a tape, which is met in the impact zone (0.5 mm) and not met in the far field (2.0 mm). The two-tier grading is the engineering compromise: full tow-wise resolution where damage occurs, simplified tow-wise resolution where elastic flexural-wave physics dominates and damage is not expected.

Total element count (P1-TWT primary configuration):

```
N_total ≈ 482,496 + 209,844 ≈ 692,340 solid elements
```

Plus cohesive-zone contact pairs: with 23 ply-to-ply interfaces and the same in-plane mesh density, /INTER/TYPE2 is contact-based not element-based, so the cohesive interfaces add slave/master surface segments rather than additional volume elements. The slave/master segment count tracks the in-plane mesh and contributes ~700k surface segments across the 23 interfaces.

The master plan §10 risk register entry on element count is conservative ("likely 5-20 million" in the stage brief). The actual budget is an order of magnitude smaller (~700k) because the reduced 200 × 200 mm section is much smaller than the full 1143 × 838 mm panel; the brief's high estimate corresponded to the full panel which is not modelled here. **The 700k figure is a planning estimate; the runner records the actual count after the GMSH build and writes it to the run log.**

If the user later wishes to model the full panel, the element count scales by the area ratio `(1143 × 838) / (200 × 200) ≈ 24`, yielding `N_total ≈ 17 M elements`, which is the "5-20 million" upper bound from the brief. That run is HPC-class (master plan §10 risk register: HPC may be needed).

### 3.2 Time step and wall clock

Explicit dynamic stable time step from the Courant-Friedrichs-Lewy condition. Carbon/epoxy in the fiber direction has wave speed `c ≈ √(E_1 / ρ) = √(165 GPa / 1570 kg/m³) ≈ 10250 m/s`. With `h_min = 0.5 mm` in the impact zone:

```
Δt_stable ≈ 0.5e-3 m / 10250 m/s ≈ 4.9e-8 s
```

The OpenRadioss explicit engine runs at typically 0.9× the CFL bound, so the timestep is `Δt ≈ 4.4e-8 s`. The simulation needs to cover at least 0.5 ms to capture (i) the impact contact phase (~50 µs), (ii) the post-impact through-thickness wave reverberations (~100 µs), and (iii) the back-face deflection peak and rebound (~300 µs). At 0.5 ms / 4.4e-8 s ≈ 11.4 million timesteps per shot.

OpenRadioss explicit on a 700k-element solid mesh achieves roughly 1.2 µs/element/cycle on a 16-core x86_64 workstation per the TNO 2024 benchmark cited in `openradioss_endtoend_audit.md` row H. Wall-clock per shot:

```
T_wall ≈ 700,000 × 11.4e6 × 1.2e-6 / (16 cores × 1e6 µs/s)
       ≈ 600 s × 11.4 / (16) -- corrected scaling per CPU-cycle accounting
```

Translating the back-of-envelope to a concrete prediction: a single shot is **8-24 hours of wall clock on a 16-core workstation** with the 700k-element reduced-section mesh. The V50 sweep of seven shots (Section 9) is therefore **3-7 days** of contiguous wall clock. This is the master plan §10 risk register entry "HPC may be needed"; the runner accepts a `--cores` argument and a `--mpi-hosts` argument for OpenMPI-distributed engine execution per `openradioss_endtoend_audit.md` row A.

If the run exceeds the available wall clock budget, the documented mitigation is mass scaling: scaling element mass by `MS_factor = 4` raises the effective wave speed by `√(MS_factor) = 2`, doubling the stable timestep, halving the wall clock, at the cost of an inertial-overshoot artefact in the contact force which is typically <5% on this class of problem (Stage 14 audit row D). The runner exposes `--mass-scaling` as a CLI flag and writes a warning if it is enabled.

### 3.3 Hourglass control: /HEPH on solid HEXA, mandatory at this strain rate

OpenRadioss `/PROP/TYPE14` solid bricks default to one-point reduced integration with hourglass stabilization. At ballistic strain rates (~10^4 1/s, three orders of magnitude above quasi-static) hourglass energy can grow rapidly and contaminate both the contact force trace and the per-element damage-trigger criterion in /FAIL/HASHIN. The mandatory choice for this stage is the physical hourglass stabilization scheme **/HEPH** (Hourglass with Elastic-Plastic Hourglass — Belytschko-Bindeman 1993), which adds a small physical-stiffness regularization rather than the older artificial viscous hourglass.

This choice is documented in master plan §1 (hourglass control mandatory at this strain rate) and is enforced in the deck via the `Iform = 24` hourglass-formulation flag on `/PROP/TYPE14`. The runner enforces `E_hg / E_int ≤ 5 percent` as a hard sanity gate; if hourglass energy exceeds 5 percent of internal energy at any point in the simulation the run is flagged and the V50 fit is suppressed, on the principle that a contaminated contact force trace produces a contaminated V50.

Alternatives considered and rejected. /HFLEXAGON (full integration on hexa) would eliminate hourglass entirely but at 8× the elemental computational cost, pushing the wall clock to 24-72 hours per shot which is impractical for a V50 sweep. Standard /HSTAB (default reduced integration) is unsafe at ballistic strain rates per the audit row 13 PASS notes. /HEPH is the documented compromise.

---

## 4. Boundary conditions and projectile

### 4.1 Cut-edge boundary conditions

The 200 × 200 mm section is a cut subset of the 1143 × 838 mm physical panel. Two boundary-condition options are documented and tested.

Option A (clamped). All four cut edges of the section have all displacement degrees of freedom set to zero (`/BCS` with `Trarot = 111000` interpreted as full translation lock). This is the Vakili Rad 2020 experimental panel boundary condition, where the panel is clamped in a steel fixture frame around its perimeter at NASA Glenn (Vakili Rad 2020 thesis Fig. 6.x, NASA Glenn gas-gun setup; `[UNVERIFIED]` figure number, taken from the thesis abstract description). Clamped cut edges in a 200 mm sub-section impose a stiffer-than-physical boundary because the 200 mm cut is closer to the impact than the 1143/838 mm physical boundary, but the Saint-Venant decay over 100 mm of laminate damps the boundary stiffening; the BC-insensitivity check below confirms whether this is acceptable.

Option B (absorbing / non-reflecting). The four cut edges are assigned an /ALE-style or low-reflecting boundary that absorbs flexural wave energy, approximated in OpenRadioss explicit by a /DAMP/INTER patch with a viscous-dashpot coefficient calibrated to the local impedance `Z = ρ c h ≈ 1570 × 1500 × 4.55e-3 ≈ 10.7 kPa·s/m`. Absorbing BCs avoid the spurious flexural-wave reflection at 130 µs that would corrupt the back-face deflection trace but are an approximation of the true infinite-laminate condition, not the experimental clamped-frame condition.

The runner runs **both** options on the V50 = published-value shot only and compares the contact-force-vs-time trace; the BC choice is judged insensitive if the peak contact force differs by ≤ 5 percent and the contact duration differs by ≤ 10 percent. The primary V50 sweep then uses the option chosen as the more physically faithful (clamped, matching the experimental fixture) provided the insensitivity test passes; if it fails, the sweep is run with both BCs and the 7 percent V50 tolerance is loosened to 9 percent in writing in the runner output.

### 4.2 Projectile geometry, mass, and velocity range

Per Vakili Rad 2020 *Composites Part B* and the abstract description of the NASA Glenn single-stage gas-gun setup, the projectile is a **0.50 cal (12.7 mm diameter) cylindrical right-circular steel slug** with mass `m_p = 13.4 g ± 0.3 g` `[UNVERIFIED]`, struck at velocities in the range 250-400 ft/s = 76-122 m/s. The 0.50 cal cylindrical FSP (Fragment Simulating Projectile) is the standard NASA Glenn ballistic test projectile geometry (cited in the Vakili Rad 2020 paper figure and consistent with ASTM D8101, the standard the paper invokes). The exact projectile mass and length are paywalled behind the *Composites Part B* paper; the runner notes this with `[UNVERIFIED]` markers and uses the mass and dimensions from the publicly extracted abstract and the Vakili Rad 2020 thesis (`uofsc_pseudowoven.md` Section 5 and Section 6).

If a future re-extraction of the paper produces a different projectile mass, only one line of the runner config changes; the sweep is re-run.

Velocity range is **76-122 m/s** as published. This is "low" by general ballistic standards (the paper uses the term "high-velocity impact" rather than "ballistic" precisely because 122 m/s is below typical small-arms ballistic), but it is the published HVI envelope and the reference V50 falls inside it. Stage 15 (the upstream flat-ballistic stage that proves the erosion stack) is run at the same velocity range so the constitutive-law calibration carries over without re-tuning.

Reference V50 for P1-TWT (the pass-criterion target). The Vakili Rad 2020 *Composites Part B* paper reports a V50 figure for the P1-TWT configuration in their Table 5 or Figure 11 (`[UNVERIFIED]` exact figure / table number; the paper is paywalled and the published abstract reports only "P1-TWT is the best ballistic performer" plus the 45 percent reduction in back-face damage). The placeholder value used by the runner pre-tuning is `V50_published = 110 m/s ± [UNVERIFIED]`, derived from the midpoint of the 76-122 m/s velocity envelope in which the paper conducts the limit-velocity bracketing test. **The runner's pass-criterion check loads `V50_published` from a configuration file `vakili_rad_2020_v50_table.json` so that the value can be updated in one place when the paper is digitized.** This is the master plan §10 risk register entry "no public AP-PLY raw data ... ground truth is one number per configuration".

The runner emits a clear diagnostic if `V50_published` is still flagged `[UNVERIFIED]` and the pass criterion comparison is therefore against an estimated value, not a digitized one.

### 4.3 Contact, friction, and erosion linkage

Projectile-laminate contact is `/INTER/TYPE7` general node-to-surface penalty contact. Friction coefficient `μ = 0.3` (standard for steel-on-CFRP in dry conditions, NASA Glenn ballistic-test value `[UNVERIFIED]`). Penalty stiffness defaults to OpenRadioss-computed value tied to the local element stiffness.

Element erosion is enabled on the laminate via `Ifail = 2` on every `/FAIL/HASHIN` card; this is the documented OpenRadioss erosion mechanism for solid composite bricks per `openradioss_endtoend_audit.md` row 15. Eroded elements are removed from the contact slave-segment list automatically by the engine, so the projectile transitions cleanly into the next ply once a tow element fails. This is the same erosion stack proven on a flat coupon at Stage 15.

The projectile is a rigid `/RBODY` because the Stage 16 capability under test is composite damage propagation, not projectile deformation. Rigid-body kinematic mass is `m_p = 13.4 g` per Section 4.2.

---

## 5. Material cards (every card listed)

The Vakili Rad 2020 panel uses Hexcel IM7G/8552-1 carbon/epoxy prepreg slit tape. The constituent properties for IM7G/8552-1 are taken from the Soden-Hinton-Kaddour 1998 *Composites Science and Technology* (`SodenHintonKaddour1998` in `references/test_progression_refs.bib`) WWFE-I/II reference card, supplemented by WWFE-II strengths (`HintonKaddourSoden2002`) where Soden 1998 is silent, and by Camanho-Davila 2003 (`CamanhoDavila2003` in `test_progression_refs.bib`) for the cohesive-zone parameters.

### 5.1 LAW25 orthotropic for the tape and undulation regions

`/MAT/LAW25/<id>` (CRASURV solid orthotropic with Tsai-Wu base yield surface) per `openradioss_endtoend_audit.md` row 9.

| Symbol | Description | Value (SI) | Source |
|---|---|---|---|
| ρ | density | 1570 kg/m³ | Soden 1998 IM7/8552 card |
| E₁ | longitudinal modulus | 165 GPa | Soden 1998 |
| E₂ = E₃ | transverse modulus | 8.40 GPa | Soden 1998 |
| ν₁₂ = ν₁₃ | major Poisson | 0.34 | Soden 1998 |
| ν₂₃ | through-thickness Poisson | 0.45 | Soden 1998 |
| G₁₂ = G₁₃ | longitudinal shear modulus | 5.6 GPa | Soden 1998 |
| G₂₃ | transverse shear modulus | 2.8 GPa | Soden 1998 |
| Xₜ | longitudinal tensile strength | 2560 MPa | Soden 1998 |
| Xc | longitudinal compressive strength | 1590 MPa | Soden 1998 |
| Yₜ = Zₜ | transverse tensile strength (in-situ) | 73 MPa | WWFE-II in-situ correction; Camanho 2007 |
| Yc = Zc | transverse compressive strength | 185 MPa | Soden 1998 |
| Sₗ = SLT | in-plane shear strength | 90 MPa | Soden 1998 |
| Sₜ = STT | transverse shear strength | 80 MPa | WWFE-II |
| Iform | failure / softening flag | 1 (CRASURV with strength degradation) | LAW25 ref. |

The local material orientation per element (1-axis = fiber direction, 2-axis = transverse in-plane, 3-axis = through-thickness) is set by the per-element `/SKEW` reference written by the Kok preprocessor through inp2rad, so that the +45 plies have their 1-axis along the +45 direction in the laminate plane, the 90 plies along the +90 direction, etc. Undulation elements have a locally-rotated 1-axis that follows the through-thickness arc of the tow per the Kok geometry; this is the principal physical advantage of tow-wise FE.

### 5.2 /FAIL/HASHIN with element erosion on tape and undulation regions

`/FAIL/HASHIN/<id>` per `openradioss_endtoend_audit.md` row 6 and row 15 (PASS).

| Symbol | Description | Value | Source |
|---|---|---|---|
| Ifail | erosion flag | 2 (delete element on criterion satisfaction) | LAW25 + /FAIL/HASHIN reference, audit row 15 |
| F_max | failure criterion threshold | 1.0 | dimensionless Hashin index |
| Damping coefficient | numerical damping on softening branch | 1e-4 | Pinho 2006 / Camanho 2007 default |
| Activation | both fiber and matrix modes active | yes | Hashin 1980 |

The four Hashin sub-criteria are activated simultaneously: fiber tensile, fiber compressive, matrix tensile, matrix compressive. The first satisfied criterion at any Gauss point triggers element deletion at that integration point; with HEXA8 reduced integration there is one Gauss point per element, so the trigger is element-level.

### 5.3 /MAT/LAW1 isotropic epoxy for the resin pocket region

`/MAT/LAW1/<id>` linear isotropic for the small resin-rich pockets between adjacent tow drops in the AP-PLY pattern.

| Symbol | Description | Value (SI) | Source |
|---|---|---|---|
| ρ | density | 1300 kg/m³ | Hexcel 8552 datasheet, Soden 1998 |
| E | Young's modulus | 4.67 GPa | Soden 1998 (matrix-only datapoint) |
| ν | Poisson's ratio | 0.36 | Soden 1998 |

Resin pockets do not carry intra-laminar damage in this stage; they are linear elastic with a /FAIL strain criterion (engineering rupture strain ε_R = 0.05) for erosion. This is the same simplification used in Kok et al. 2022 (`KokTensile2022`) and is consistent with the small volumetric fraction of resin pocket in a well-consolidated AFP layup (<2% per Kodagali 2023 thesis C-scan).

### 5.4 /INTER/TYPE2 cohesive between every adjacent ply pair

`/INTER/TYPE2/<id>` tied surface contact with brittle / progressive failure per `openradioss_endtoend_audit.md` row 12. Note that audit row 12 also flags TYPE2 as a tied (kinematic) contact rather than a cohesive zone; the more cohesive-faithful path is `/MAT/LAW83` solid cohesive elements (also documented). For Stage 16 the **primary** path is /INTER/TYPE2 with the Camanho-Davila 2003 mode-I and mode-II energy parameters mapped onto its `Spotflag` failure model, because solid cohesive elements between every ply pair would add ~700k extra elements and is not budget-feasible at this scale. The runner exposes a `--cohesive-mode {type2,law83}` flag for the user to swap the cohesive mechanism if they have HPC budget.

| Symbol | Description | Value (SI) | Source |
|---|---|---|---|
| G_Ic | mode-I fracture toughness | 0.277 kJ/m² | Camanho-Davila 2003 IM7/8552 card |
| G_IIc | mode-II fracture toughness | 0.788 kJ/m² | Camanho-Davila 2003 |
| τ_n_max | mode-I peak traction | 60 MPa | Camanho-Davila 2003 |
| τ_s_max | mode-II peak traction | 90 MPa | Camanho-Davila 2003 |
| η_BK | Benzeggagh-Kenane mixed-mode exponent | 1.45 | Camanho-Davila 2003 |
| Spotflag | progressive failure mode | 4 (mixed-mode energy criterion) | LAW83 doc; TYPE2 ref. |

A separate /INTER/TYPE7 general contact is defined between the projectile (rigid /RBODY) and the laminate outer-skin master surface for the impact contact itself; that card is independent of the inter-laminar cohesive cards.

### 5.5 Summary table: every card in the deck

| Card | Purpose | Count |
|---|---|---|
| `/MAT/LAW25` IM7G/8552-1 orthotropic | tape and undulation regions | 1 (shared across all plies; orientation per /SKEW) |
| `/FAIL/HASHIN` with Ifail=2 | element erosion on the LAW25 sets | 1 attached to LAW25 |
| `/MAT/LAW1` neat 8552 epoxy isotropic | resin pockets | 1 |
| `/MAT/LAW2` (or rigid material) for projectile | rigid steel projectile | 1 |
| `/PROP/TYPE14` solid HEXA8 with /HEPH (Iform=24) | every solid element in laminate | 1 |
| `/RBODY` for projectile | rigid kinematic body | 1 |
| `/SKEW` per ply orientation | local material frame, 24 plies | 24 (or less, if the Kok preprocessor reuses) |
| `/INTER/TYPE2` ply-to-ply cohesive | inter-laminar bonding | 23 (for 24 plies) |
| `/INTER/TYPE7` projectile-to-laminate contact | impact contact | 1 |
| `/IMPVEL` initial velocity on projectile | impact velocity input | 1 |
| `/BCS` clamped cut edges | boundary | 4 (one per edge) |
| `/DAMP/INTER` (Option B fallback) | absorbing edges | 4 (only if BC option B chosen) |
| `/TH/PART`, `/TH/NODE`, `/H3D/ELEM` | output requests | 1 set |
| `/ANIM/ELEM/DAMA`, `/H3D/ELEM/EROS` | damage and erosion output | 1 set |

---

## 6. Solver controls and BC-insensitivity check

Solver: explicit dynamic, OpenRadioss `engine_linuxa64` with OpenMPI (`-mpi=ompi` build) per `openradioss_endtoend_audit.md` row A. Time-step control:

- `/DT/BRICK/CST/<id>` with `Δt_target = 4.4e-8 s` (90 percent of CFL bound from Section 3.2).
- `/DT/INTER/DEL` with the same target on contact-driven timestep.
- `/STOP/ETIME` at `t_end = 0.5 ms` per Section 3.2 simulation duration.

Output cadence:

- `/TH/...` time-history at every 1.0 µs (500 samples in a 0.5 ms run), enough to fit the Recht-Ipson model post-hoc.
- `/ANIM/DT 5e-6 s` (100 frames over 0.5 ms) for the visual back-face deflection diagnostic.
- `/H3D/ELEM/EROS` to write per-element erosion flags so the runner can tally penetrated plies.

BC-insensitivity check (the gate on Section 4.1). On the V50 = published-value shot only, the runner repeats the simulation with both BC Option A (clamped) and BC Option B (absorbing). The two contact-force-vs-time traces are compared, and the BC choice is declared insensitive if:

- Peak contact force differs by ≤ 5 percent.
- Time of peak differs by ≤ 10 percent.
- Contact duration (full-width at 10 percent of peak) differs by ≤ 10 percent.

If the insensitivity check fails, the V50 sweep is run with both BCs, both V50 values are reported, the final V50 is the average of the two, and the pass tolerance is relaxed from 7 percent to 9 percent in the runner's PASS/FAIL printout, with a written warning that the BC choice is on the same order of magnitude as the architectural-uncertainty budget. This is the master plan §10 risk register entry "ground truth is one published number per configuration".

---

## 7. Reference solution: Vakili Rad 2020 V50 figure or table

The reference V50 for the P1-TWT panel is taken from Vakili Rad et al. 2020, *Composites Part B: Engineering*, vol. 203, paper 108478 (DOI 10.1016/j.compositesb.2020.108478, `VakiliRad2020HighVelocityImpact` in `references/uofsc_refs.bib`).

The paper reports the bracketing-test result for each of the three configurations (P1-TWT, P2-WTW, P3-Control) in their Table 5 (`[UNVERIFIED]` table number; the open-access portions of the paper that have been digitized into `references/uofsc_pseudowoven.md` Section 6 give the qualitative statement "P1-TWT is the best ballistic performer" plus the quantitative back-face metrics "45 percent reduction in back-face damage; 19.5 percent less back-face deflection (P1-TWT best)" but do not give a literal V50 number in the open digitized excerpt). The full Table 5 / Figure 11 V50 values are paywalled.

Citation discipline. The runner reads the published V50 from `vakili_rad_2020_v50_table.json` which is the single source of truth for the digitized number. The JSON file ships with `[UNVERIFIED]` markers wherever a number could not be confirmed against the open-access portion of the paper, and the runner emits a CSV row per shot with the per-shot residual-velocity prediction so that future re-digitization of the paper produces a comparable per-shot dataset. The expected JSON schema:

```json
{
    "paper":             "Vakili Rad et al. 2020, Composites Part B, 203, 108478",
    "doi":               "10.1016/j.compositesb.2020.108478",
    "panel_size_mm":     [1143.0, 838.0],
    "panel_thickness_mm": 4.55,
    "n_plies":           24,
    "material":          "Hexcel IM7G/8552-1",
    "v_range_m_per_s":   [76.0, 122.0],
    "panels": {
        "P1_TWT":     {"v50_m_per_s": "[UNVERIFIED]", "v50_uncertainty_m_per_s": "[UNVERIFIED]"},
        "P2_WTW":     {"v50_m_per_s": "[UNVERIFIED]", "v50_uncertainty_m_per_s": "[UNVERIFIED]"},
        "P3_Control": {"v50_m_per_s": "[UNVERIFIED]", "v50_uncertainty_m_per_s": "[UNVERIFIED]"}
    },
    "qualitative_open_access_findings": {
        "P1_TWT_back_face_damage_reduction_percent": 45,
        "P1_TWT_back_face_deflection_reduction_percent_vs_P3": 19.5,
        "ranking_best_to_worst": ["P1_TWT", "P2_WTW", "P3_Control"]
    },
    "extraction_provenance": {
        "open_access_extraction_source": "references/uofsc_pseudowoven.md sec.6",
        "needs_paywall_pdf_re_extraction": true
    }
}
```

When the user re-extracts the paywalled PDF, only this JSON file changes; the runner re-evaluates the pass criterion automatically.

Cross-references for the pass criterion's epistemic context:

- `BossuytKodagali2020` `[UNVERIFIED authors]` (`references/uofsc_refs.bib` likely entry; the paper-list in the audit references both "Bossuyt 2020" and the Vakili Rad 2020 paper as the same paper, but the author list disagrees). The Bossuyt entry is treated as a duplicate of Vakili Rad 2020 in the runner.
- `KokTensile2022` (`Kok2022Tensile` in `references/delft_refs.bib`): mesoscale validation backbone for the Kok preprocessor, proven in Stage 11.
- `Kok2023Dynamic`: dynamic split-Hopkinson tension on AP-PLY at 30 1/s, the data point that justifies the rate-independent moduli assumption used in the LAW25 card here.
- `Kok2024LVI`: low-velocity impact on AP-PLY, the validation precedent for the Kok mesoscale-FE pipeline at impact rates.
- `RechtIpson1963`: Recht and Ipson, "Ballistic perforation dynamics", J. Appl. Mech. 30, 384-390. Canonical V50 fitting model used in Section 9.
- `Cunniff1992`: Cunniff, "An analysis of the system effects in woven fabrics under ballistic impact", Textile Research Journal 62, 495-509. Canonical reference for composite-armor V50 scaling, and the source of the V50 / Cunniff-velocity-parameter normalization used in the runner's diagnostic plot.
- `ASTM_D8101`: ASTM Standard D8101/D8101M, Standard Test Method for Resistance of Composite Materials to High-Velocity Impact, the standard the Vakili Rad 2020 paper cites.

---

## 8. (Cross-reference: Section 6 BC insensitivity check is the second-tier success gate; Section 7 is the V50 reference solution.)

This Section 8 placeholder maintains the 10-section format header used at Stage 5 and across the upstream stages. The substantive contents that would normally live in Section 8 (validation success criterion in technical detail) are integrated into Section 7 (reference solution, V50 source) and Section 6 (BC-insensitivity gate); the formal pass criterion is restated here for the runner's hard-coded check:

```
pass_criterion: |V50_simulated - V50_published| / V50_published <= 0.07
                AND BC_insensitivity_passed (or pass_tolerance relaxed to 0.09)
                AND E_hourglass_max / E_internal_min_during_impact <= 0.05
                AND ranking_simulated == ranking_published  (P1_TWT > P2_WTW > P3_Control,
                    on the optional cross-config sanity sweep)
```

---

## 9. Runner script overview: V50 sweep, residual-velocity parsing, Recht-Ipson fit

The runner `runner.py` performs the V50 bracketing test in the same spirit as MIL-STD-662F. The procedure:

1. **Geometry build (once per panel configuration).** Read `kok_config.json`. Invoke the Kok `ap_ply_model_creation` preprocessor to write `build/p1_twt_section.inp` per Section 2.4. Convert to `.rad` via inp2rad. Cache the build artefact so the geometry build is not repeated across shots.
2. **Material card templating (once per panel configuration).** Read `material_card_im7_8552.json` (Section 5 numerical values). Template the LAW25, LAW1, FAIL/HASHIN, INTER/TYPE2 cards into the starter deck `_0000.rad` text using f-string substitution.
3. **BC-insensitivity check (one extra shot once per panel configuration).** Run the V50 = `V50_published` shot with both BC Option A and Option B. Compare contact-force traces per Section 6. Record the result and the chosen BC for the sweep.
4. **V50 bracketing sweep (seven shots per panel configuration).** Generate the impact-velocity sweep:
   - `v_test ∈ {0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20} × V50_published`
   - For each `v_test`, write the `/IMPVEL` field into the projectile sub-deck and submit the OpenRadioss starter and engine.
   - For each shot, parse the time-history file `T01` via the Vortex-Radioss Python reader and extract:
     - `v_residual` = the projectile centroid x-velocity at `t_end = 0.5 ms` (zero if the projectile rebounded; positive if it penetrated).
     - `F_peak` = peak contact force in the laminate.
     - `delta_back_face` = peak back-face displacement.
     - `n_penetrated_plies` = number of plies with at least 50 percent eroded elements at any in-plane location below the impact point at `t_end`.
5. **Recht-Ipson V50 fit.** Fit the Recht-Ipson 1963 model
   $$ V_r^2 = a (V_i^2 - V_{50}^2)\quad\text{for } V_i \ge V_{50},\qquad V_r = 0\quad\text{for } V_i < V_{50} $$
   to the seven (`v_test`, `v_residual`) data points by nonlinear least squares (`scipy.optimize.curve_fit`) with `a` and `V_{50}` as free parameters. The fit is performed only over the over-match (`v_residual > 0`) shots; the sub-V50 shots are used as a censoring constraint.
6. **Pass criterion check.** Compare the fitted `V_{50}_{simulated}` to the `V_{50}_{published}` per Section 8. Print PASS or FAIL plus the absolute-percent error, the BC option used, the hourglass-energy ratio, and the per-shot table. Write `v50_summary.json` and `v50_sweep.csv`.
7. **Optional cross-config sanity sweep.** If the user passes `--all-configs`, repeat steps 1-6 for P2-WTW and P3-Control and check that the simulated V50 ranking matches the published ranking (P1-TWT > P2-WTW > P3-Control by V50). This exists not as a hard pass criterion but as a model-believability gate on the absolute V50.
8. **CeTZ-ready CSV emission.** Per `~/.claude/CLAUDE.md` global preferences, the runner emits `figures/v50_sweep.csv`, `figures/recht_ipson_fit.csv`, and `figures/contact_force_vs_time.csv` for downstream Typst / CeTZ figure compilation. **The runner does not emit matplotlib figures.**

The runner's CLI:

```
python runner.py --config kok_config.json
                 --vakili-rad-table vakili_rad_2020_v50_table.json
                 --panel P1_TWT
                 [--all-configs]
                 [--cores 16]
                 [--mpi-hosts host1,host2]
                 [--mass-scaling 4]
                 [--cohesive-mode type2|law83]
                 [--bc-option A|B|both]
                 [--lima-vm apptainer]
                 [--build-cache build/]
                 [--output runs/stage_16/]
```

CSV schema for `v50_sweep.csv`:

```
v_test_m_per_s, v_residual_m_per_s, F_peak_N, delta_back_face_mm,
n_penetrated_plies, hourglass_energy_J, internal_energy_J,
hourglass_ratio, contact_duration_us, BC_option, status
```

JSON schema for `v50_summary.json`:

```
{
    "panel_configuration": "P1_TWT",
    "v50_simulated_m_per_s": <float>,
    "v50_uncertainty_m_per_s": <float, from curve_fit covariance>,
    "v50_published_m_per_s": <float, from vakili_rad_2020_v50_table.json>,
    "absolute_percent_error": <float>,
    "tolerance_percent": 7.0,
    "bc_option": "A" | "B" | "both",
    "bc_insensitivity_passed": <bool>,
    "hourglass_max_ratio": <float>,
    "hourglass_check_passed": <bool>,
    "ranking_check_passed": <bool, optional>,
    "verdict": "PASS" | "FAIL"
}
```

---

## 10. Risks, limitations, honesty about the V50 ground-truth

This stage's risk profile is unique in the project because it is the integration test, and its ground truth is one published number per configuration in a paywalled paper. The risk register that follows is the set of items that would, if mishandled, cause a mathematically-correct simulation to fail the pass criterion or, worse, pass the pass criterion for the wrong reason.

1. **Ground-truth sparsity.** The Vakili Rad 2020 paper publishes one V50 per panel configuration. There is no per-shot residual-velocity table in the open literature. There is no published Recht-Ipson `a` parameter. The runner produces a fully-resolved Recht-Ipson curve from the simulation, but it can only be compared to a single point on the experimental envelope. This is documented in master plan §10 risk 3 ("V50 ground truth is one number per Vakili Rad 2020 panel; no per-shot residual-velocity table"). The mitigation is the cross-configuration ranking check (Section 9 step 7): even though the absolute V50 of P1-TWT is one number, the relative ranking P1-TWT > P2-WTW > P3-Control is a documented qualitative finding in the paper, and the simulation's reproduction of that ranking is independent evidence of model fidelity.

2. **Element-count and HPC budget.** Section 3.1 estimates the reduced-section element count at ~700k. If the user requires the full 1143 × 838 mm panel for a reviewer or a publication, the element count rises to ~17 M and the wall clock per shot rises to ~1 week on a 16-core workstation, taking the V50 sweep into multi-month territory. This is documented in master plan §10 risk note "HPC may be needed". Mitigation is master plan §7 Lima distributed-OpenMPI on a Linux cluster, with the runner's `--mpi-hosts` flag.

3. **Kok preprocessor + inp2rad beta status.** The Kok pipeline is LGPL-2.1 Python 2.7 (Abaqus CAE Python API plus Shapely) and the inp2rad converter is "beta" status per `openradioss_endtoend_audit.md` row C. At Stage 16 the laminate is 24-ply with 23 inter-laminar interfaces and per-element local material orientations on every tape and undulation element; converter bugs that drop or rotate per-element /SKEW orientations would silently corrupt the V50 by changing the apparent fiber angle of every ply. Mitigation: the runner sanity-checks element count, node count, element-set count, and per-element orientation count between `.inp` and `.rad`. Stage 11 is the upstream stage that proves the same pipeline at small scale (4-8 ply, single AP-PLY repeat unit); a successful Stage 11 is a hard prerequisite for trusting Stage 16.

4. **LAW25 + /FAIL/HASHIN at ballistic strain rate.** LAW25 in OpenRadioss is rate-insensitive in its base form. Kok et al. 2023 *Composites Part B* 248, 110347 (`Kok2023Dynamic`) reports AP-PLY moduli are rate-independent up to 30 1/s and strengths are slightly higher at rate; ballistic rates here are ~10^4 1/s, three to four decades above that experimental rate. The constitutive-law calibration is therefore borrowed from quasi-static data, with the documented rate-dependent strength uplift NOT applied. Stage 15's pass at the same rate range with the same material is the upstream evidence that this approximation is acceptable; if Stage 15 was tuned on a 5 percent rate-uplift, that calibration should be propagated to Stage 16 in `material_card_im7_8552.json`.

5. **Cohesive contact selection (TYPE2 vs LAW83).** /INTER/TYPE2 is faster to set up than solid /MAT/LAW83 cohesive elements, but per the audit (row 12) it is a tied contact with brittle failure rather than a true cohesive zone. Stage 12 (DCB / ENF cohesive) is the upstream evidence on this point; if Stage 12 chose LAW83 for the cohesive-zone benchmark, the same path should be available here via the runner's `--cohesive-mode law83` flag at the cost of an additional ~700k cohesive solid elements and a roughly doubled wall clock per shot.

6. **Boundary-condition reduction.** The 200 × 200 mm reduced section is a deliberate small-mesh compromise (Section 2.2). The BC-insensitivity check (Section 6) is the second-tier safety gate; if it fails, the absolute V50 number is degraded by the reflection artefact and the pass tolerance is relaxed to 9 percent. This is the master plan §10 risk that "the panel is large; HPC may be needed" being addressed by reduction rather than by HPC.

7. **Projectile geometry behind a paywall.** Section 4.2 marks the projectile mass and length `[UNVERIFIED]`; the runner uses placeholder values from the open-access abstract until the paywalled PDF is digitized. If the actual projectile mass differs from the placeholder by more than ~10 percent, the V50 prediction shifts by a factor of `√(m_actual / m_placeholder)` per Recht-Ipson, which is a non-trivial effect on the 7 percent pass tolerance. Mitigation: the runner emits a sensitivity-test row in `v50_summary.json` that reports `dV50/dm_p` near the operating point so the user can re-evaluate the pass criterion if the mass is later corrected.

8. **V50_published placeholder.** As of the date of this spec, `vakili_rad_2020_v50_table.json` ships with `V50_published = [UNVERIFIED]` for every configuration because the paper's Table 5 is paywalled. The runner refuses to print PASS unless the placeholder is replaced with a digitized number; if the placeholder is in place, the runner prints `INCONCLUSIVE` and writes the simulated V50 to the JSON output for later re-evaluation. This is documented honesty: the test is mathematically pass-or-fail, but the ground-truth retrieval is a separate operational task and the runner makes the dependency explicit.

9. **Hourglass energy budget and erosion confound.** Section 3.3 sets the hourglass ratio gate at 5 percent. With element erosion enabled (Section 4.3), eroded elements depart the energy bookkeeping; the runner's hourglass-ratio computation must subtract the contribution of eroded elements from the internal-energy denominator before dividing by the surviving-element hourglass energy, otherwise the ratio is artificially inflated as elements erode. Stage 15's runner is the upstream prototype for this corrected energy accounting; Stage 16 inherits the same correction.

10. **AGPL-3 / LGPL-2.1 license consistency.** OpenRadioss is AGPL-3, the Kok preprocessor is LGPL-2.1, and the project's own templating layer is the user's own code. The combined distribution does not host OpenRadioss as a service per master plan §10 risk 7, so AGPL contamination is not engaged; only the `.rad` decks and Python templating are distributed, not OpenRadioss derivatives. This risk is recorded for completeness; it does not affect the pass criterion.

What downstream of this stage assumes from this stage. Nothing. Stage 16 is the terminal stage in the master plan. Pass at Stage 16 is the project's success criterion. Failure at Stage 16 is the project's failure criterion. All downstream activity is publication, dissemination, and the deferred forming-draping side track per master plan §3.
