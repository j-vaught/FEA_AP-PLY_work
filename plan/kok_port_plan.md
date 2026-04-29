# Kok geometry preprocessor — clean-room Python port plan

**Author.** J.C. Vaught
**Date.** 2026-04-29
**Status.** Planning document. No implementation in scope. No source code from the upstream LGPL-2.1 project has been read or transcribed during the writing of this plan.
**Master plan reference.** `plan/master_plan.md` §11 (Phase 2.6 Kok geometry port).
**Consolidation cross-reference.** `plan/consolidation_review.md` item E (the Abaqus dependency this port resolves).
**License posture.** Clean-room from open-access papers. Target license MIT. Authoritative algorithmic sources: Nagelsmit 2013 PhD thesis (`Nagelsmit2013Thesis`) and Kok 2022 *Composites Part A* (`Kok2022Tensile`). Both fetched and verified in this session; section and figure references below are cross-checked against the PDF text.

---

## Source-citation honesty up front

The brief commissioning this plan named the algorithmic basis as "Nagelsmit 2013 §3.4" and "Kok 2022 §2". I directly fetched both PDFs in this session and the actual section numbering is different. I am preserving the original section numbers in `[CORRECTED]` notes wherever the brief's section pointer does not match the source.

| Brief said | Verified location in the source | Action |
|---|---|---|
| Nagelsmit 2013 §3.4 documents the laydown algorithm | The laydown concept and the `[(θ1/θ2)n×w/(θ3/θ4)n×w]_{3S}` notation are defined in Chapter 2 (Sections 2.2 "Concept" and 2.3 "Pattern Details") of Nagelsmit 2013, NOT §3.4. Chapter 3 of Nagelsmit 2013 is "Impact Behaviour of AP-PLY Laminates" and contains no laydown algorithm. | Cite Nagelsmit 2013 Chapter 2 throughout this plan. |
| Kok 2022 §2 documents the tow-wise geometry definition with undulation profiles | Kok 2022 §2 ("Materials and methods") covers the as-manufactured laydown (10 mm tape width, 3-tape-width same-set gap, laydown sequence in Figure 1) but the *idealized geometric model* (straight-tow / undulation / resin-rich region split, undulation profile equation, undulation ratio χ = t/Lu, average out-of-plane angle φ_avg integral) is in §3.1 ("Microscale model") and **Figure 4**, NOT §2. | Cite §2 for laydown parameters, §3.1 + Figure 4 for the idealized geometry equations. |

The rest of this plan uses the corrected pointers without further annotation.

---

## 1. Goal and scope

### 1.1 Restated objective

Replace the `rutger-kok/ap_ply_model_creation` Abaqus-CAE Python plugin (LGPL-2.1, requires a paid Abaqus license to run) with a clean-room Python module that consumes the same JSON configuration, performs the AP-PLY tow-laydown algorithm, and emits a meshable solid model in GMSH `.msh` format with per-tow physical groups and a sidecar JSON file recording per-physical-group fiber orientation. The port is the geometry-generation step inside the FEA pipeline shown in `plan/master_plan.md` §4. Stages 1-10 do not depend on it. Stages 11 and 16 do.

### 1.2 In scope

The port reproduces the AP-PLY laydown for the locked geometry parameter ranges from `plan/master_plan.md` §2.

| Parameter | Locked value(s) | Authoritative source |
|---|---|---|
| Slit-tape width $w$ | $\{6.35, 10, 12.7\}$ mm | Kok 2022 §2 (10 mm hand-slit); Nagelsmit 2013 §2.5 (1/4 in. = 6.35 mm); Vakili Rad 2020 thesis (1/4 in. AFP); Kok preprocessor README parameter ranges |
| Cured ply thickness $t_p$ | $\{0.18, 0.18125, 0.213\}$ mm | Nagelsmit 2013 §2.4 (AS4/8552 thermoset 0.18 mm; AS4/PEEK thermoplastic 0.125 mm cited but 0.18125 retained as the Kok README default per `references/delft_apply.md` §3); Zheng & Kassapoglou 2024 (0.18 mm AS4/epoxy); the 0.213 mm value is the alternate Kok default |
| Undulation ratio $\chi = t/L_u$ | $0.09 \le \chi \le 0.104$ | Kok 2022 Figure 4 (definition); Kok 2022 §3.4 reports $\chi = 0.0683$ for press-cured VTC401, but the master plan parameter range bounds 0.09-0.104 from the Kok preprocessor defaults documented in `references/delft_apply.md` §3 |
| Tape spacing $s$ (interlace level) | $\{1, 2, 3\}$ | Nagelsmit 2013 §2.3 ("skip factor" $n = 1, 2, 3$); the Kok preprocessor exposes the same parameter under the name `tape_spacing` per the Kok README |
| Stacking sequences | Cross-ply $[0/90]_{2S}$, quasi-iso $[0/45/90/-45]_S$ (Kok 2022 §2 baselines), the four-angle UofSC variants $[+45/+90/-45/0]$ with placement string `10001000` (Kodagali 2023 §3-4), Vakili Rad 2020 P1-TWT, P2-WTW, P3-Control 24-ply layups | Kok 2022 §2; Vakili Rad 2020 thesis ch. 6-7; Kodagali 2023 thesis ch. 4 |

The port emits, for each configuration, a single GMSH mesh file plus an orientations sidecar that downstream stages consume verbatim. The mesh contains: solid tow bodies (per-ply), undulation transition regions where tow paths cross, resin-rich pockets between adjacent tows on the same pass, and a free outer-boundary box that bounds the laminate.

### 1.3 Out of scope

The port is deliberately narrow.

1. **AFP machine-path optimization.** No G-code is parsed, no head kinematics is simulated, no air-time or compaction-roller model is included. The laydown is described by the JSON configuration and the per-ply angle plus placement string; the physical machine path that would realize that laydown is treated as the operator's responsibility, identical to how Nagelsmit 2013 §2.5 ("the programming of the fibre placement machine is different ... involves a lot of hand labour of the operator") and Kok 2022 §2 ("laid up by hand in a process emulating automated fiber placement") describe the as-manufactured panels. Any future AFP path-planning module is a separate research track.
2. **Defect modeling.** Tow gaps caused by tape-width tolerance, tow drops, tow overlaps from steering errors, voids, fiber misalignment, and resin-starved regions are not in scope. The port produces an as-designed (perfect) geometry. Defect-aware tow-wise modeling (Li 2021 *Composites Part A* tow-wise framework) is acknowledged in `references/uofsc_pseudowoven.md` §7 as a separate body of work.
3. **Preform draping.** The port assumes the laminate is flat. Curved-tool draping (the side-track explicitly deferred in `plan/master_plan.md` §3) is not in scope. A future draping stage would add a kinematic mapping from the 2-D laminate plane to a 3-D mold surface; nothing in the port should foreclose such an extension, but no code is written for it.
4. **Cure-cycle warpage.** The Vakili Rad 2019 *Composite Structures* "clutch laminate" paper documents up to 58 percent warpage reduction from the AP-PLY architecture, but the cause (residual cure stress relief) requires a coupled thermo-mechanical solve which is outside the geometry-generation step. The port produces nominal flat laminate geometry; any cure-warpage prediction is an OpenRadioss postprocess.
5. **Damage / cohesive interface modeling.** The Kok upstream repository ships both a geometry preprocessor (`ap_ply_model_creation`, LGPL-2.1) and a constitutive-law VUMAT (`composite_cdm_ap_ply`, LGPL-2.1, Fortran). The port reproduces only the *first*. The user's solver-side damage cards (LAW25 + /FAIL/HASHIN per `plan/master_plan.md` §1) are independent of both upstream artifacts.

### 1.4 Why clean-room

The Kok upstream project is LGPL-2.1. Reading and transcribing its source would propagate LGPL-2.1 onto every downstream artifact (`LGPLv21`). The user requires open-source-only with a permissive license. The clean-room posture — read the published papers (Nagelsmit 2013 thesis, open access; Kok 2022 *Composites Part A*, CC-BY) but never the upstream source code — keeps the port out of the LGPL derivative-work boundary and lets it ship under MIT (`MITLicense`). The Kok upstream README and top-level file list (`Kok2022PreprocRepo`) are read for orientation only; that is not a derivative-work threshold under LGPL-2.1 §0 and §2 because no expressive content from the source is incorporated.

---

## 2. Algorithmic decomposition (pseudocode)

This section names every primitive the port must implement, in pseudocode rather than Python. Each primitive is followed by the section/figure reference it derives from. Open questions — places where the published papers are silent or ambiguous — are listed at the end.

### 2.1 Configuration intake

The port consumes a JSON configuration file (schema in §3 below) describing: the laminate's in-plane dimensions, the per-ply nominal fiber angle, the per-ply placement string (active/inactive AFP channels), the tape spacing $s$, the tape width(s) $w$, the cured ply thickness $t_p$, the undulation ratio $\chi$, and any per-ply angle shift. The placement string follows the UofSC convention `[fiber angles][placement seq.][angle shift][tow width]` from Kodagali 2023 thesis ch. 4 (cited in `references/uofsc_pseudowoven.md` §5.1, table row "Placement-sequence notation"). The Nagelsmit 2013 ch. 2.3 notation `[(θ1/θ2)n×w/(θ3/θ4)n×w]_{3S}` is supported as an aliased form (§3.4 below).

### 2.2 Coordinate system and ply layering

Define the global laminate frame $(X, Y, Z)$ with $Z$ normal to the laminate. The laydown is built in $Z$ from $z = 0$ at the bottom face up to $z = n_p \cdot t_p$ at the top face, where $n_p$ is the ply count. Ply $i$ occupies $z \in [(i-1) \cdot t_p, \, i \cdot t_p]$ in the nominal (no-undulation) idealization. Nagelsmit 2013 §2.4.2 reports that AP-PLY thicknesses match unidirectional baselines to within measurement scatter, so we use a constant nominal $t_p$ across plies; deviations from nominal thickness are in scope only as the undulation perturbation defined in §2.5 below.

### 2.3 Per-ply tow path generation

For each ply $i$ with nominal angle $\theta_i$ and placement string $P_i$, generate a list of tow centerline polylines in the in-plane frame.

```
for each ply i in 1 .. n_p:
    place tow centerlines parallel to direction theta_i across the panel
    spacing between consecutive centerlines is (s + 1) * w  for tape spacing s
        (one active channel followed by s inactive channels;
         next active channel begins (s+1)*w later)
    apply per-ply origin offset to realize the placement string P_i
        (the offset shifts the active-channel pattern across the bandwidth so
         that subsequent passes at the same angle fill the gaps left by the
         previous pass)
    each centerline is a straight line in the (X, Y) plane
        (no curved AFP steering in scope per §1.3)
    clip each centerline to the panel polygon
```

Source: Nagelsmit 2013 §2.2 ("Concept") describes this as "instead of placing multiple fibre bands directly next to each other, bands are placed in two directions leaving space in between adjacent bands. These spaces are filled up with alternating bands in both directions, such that on every location two bands are positioned on top of each other." Nagelsmit 2013 §2.3 names the parameters: bandwidth $w$, skip factor $n$ (= $s$ in our notation), angle $\theta$, number of layers to interweave. Kok 2022 §2 ("In both laminates a gap of three tow widths was left between tows placed in the same set") confirms the $s = 3$ case used by Edinburgh.

### 2.4 Per-tow cross-section profile (straight-tow regions)

For each tow centerline, build a straight-extruded cross-section. The cross-section in the plane orthogonal to the centerline is a rectangle of width $w$ and height $t_p$.

```
for each tow centerline:
    define a rectangle of width w and height t_p in the plane
       perpendicular to the centerline direction
    sweep the rectangle along the centerline to produce a swept solid
    the resulting solid is the "straight-tow" body of this tow
```

Source: Kok 2022 §3.1 ("Microscale model") and Figure 4 explicitly identify the straight-fiber-tow region as a rectangular cross-section. The squared-edge approximation matches Nagelsmit 2013 §2.3 micrograph observation that "the edge of the tape is smeared out by the undulation, smoothening the angle of the undulating tow". The smoothing is captured by the undulation primitive (§2.5), not the straight-tow primitive.

### 2.5 Undulation primitive at tow-tow crossings

When two tows from different plies cross — i.e., an active channel of ply $j$ at angle $\theta_j$ overlaps in $(X, Y)$ with an active channel of ply $i$ at angle $\theta_i$ for $j \neq i$ and $\theta_i \neq \theta_j$ — the *physical* manufactured laminate forces one of the two tows to deviate locally in $Z$ to make room for the other. This is the "undulation" or "crimp" region.

The Kok 2022 idealization (Figure 4) defines the undulation by:

- An **undulation length** $L_u$ in the in-plane direction along the undulating tow's centerline.
- An **undulation amplitude** equal to the ply thickness $t_p$, because the undulating tow climbs over a parallel tow of nominal thickness $t_p$.
- An **undulation ratio** $\chi = t_p / L_u$, defined explicitly in Kok 2022 Figure 4.
- An **out-of-plane angle profile** $\varphi(x)$ along the undulation, parameterized as
  $\varphi(x) = \tan^{-1}\bigl(\,\partial z(x) / \partial x\,\bigr)$
  where $z(x)$ is the local through-thickness centerline shape. Kok 2022 Figure 4 schematic indicates this is computed from the tow centerline traversing one ply thickness over half the undulation length, i.e., $z(x)$ rises by $t_p$ over $L_u / 2$ and the *average* angle is
  $\varphi_{\text{avg}} = (1 / (L_u/2)) \int_0^{L_u/2} \varphi(x)\, dx$
  Kok 2022 §3.1 uses $\varphi_{\text{avg}}$ as the tow-frame rotation about the in-plane perpendicular axis, applied uniformly inside the undulation region (each undulation "unit cell" gets a single constant out-of-plane rotation). The exact closed-form $z(x)$ shape is **not specified** in Kok 2022 §3.1 (see open question OQ-1 below). What Kok 2022 commits to is: (i) the ratio $\chi = t / L_u$, (ii) the climb amplitude $t_p$, (iii) the in-plane location of the undulation between two tow boundaries, and (iv) the average-angle integral that produces the tow-frame rotation used by the constitutive law.

Pseudocode for the undulation primitive:

```
for each pair (ply j active tow at angle theta_j, ply i active tow at angle theta_i)
    where theta_i != theta_j and the two tows overlap in (X,Y):

    locate the in-plane intersection region where one tow crosses over the other
    set L_u from the configured undulation ratio chi = t_p / L_u
        (so L_u = t_p / chi)
    half of the undulation (length L_u/2) sits above the crossover line on the
        undulating tow's centerline; the other half (length L_u/2) sits below
    inside the undulation length, the tow centerline climbs in z by t_p
    the in-plane footprint of the undulation region is a rectangle of width w
        (tape width) and length L_u centered on the crossover
    compute phi_avg as above and store with the undulation region
    the resulting solid is the "undulation" body of this crossover
```

The Nagelsmit 2013 §2.3 "Undulations" subsection observes the same primitive empirically: "the undulating area of one tow is bandwidth × thickness". Nagelsmit 2013 §2.4.1 reports a 5-degree (thermoplastic AS4/PEEK) and 7-degree (thermoset AS4/8552) undulation angle in the as-manufactured laminate, which are consistent with $\arctan(2\chi)$ for $\chi \approx 0.04 - 0.06$ on those material systems. For our locked $\chi$ range 0.09-0.104, the corresponding $\varphi_{\text{avg}}$ is approximately 10-12 degrees.

### 2.6 Resin-rich pocket primitive

Where two adjacent tows on the same pass are separated by less than one full active-channel spacing (the placement-string boundary between an active and inactive channel), or where the laydown leaves a triangular cusp at a tow boundary, the manufactured laminate fills the space with neat resin. The Kok 2022 §3.1 idealization treats these as a third region type (the "resin-rich" unit cell shown in Figure 4) that is a *mixture* of unidirectional tow constituent and pure resin, not a pure-resin pocket. The port follows that idealization: the resin-rich region is the residual volume in each ply that is not occupied by either a straight-tow body (§2.4) or an undulation body (§2.5). Whether this volume is treated as pure resin (LAW1) or as a mixed constituent at the constitutive-law level is a downstream decision (`tests/stage_11_PW_mesoscale_direct/spec.md` §2.2.2 chooses pure resin LAW1 for stage 11; the port emits the region as a separate physical group either way).

```
for each ply i:
    let ply_volume = panel_polygon × [z_i_bottom, z_i_top]
    let occupied = union of all straight-tow bodies and undulation bodies in ply i
    let resin_pocket = ply_volume - occupied
    tag resin_pocket as a separate physical group "RESIN_PLY_<i>"
```

The Kok 2022 §3.1 description of the resin-rich region as a mixed-constituent unit cell (Figure 4 "Resin Tow" and "Resin + Tow" unit cells) is acknowledged but not enforced at the geometry layer; the constitutive choice is made downstream.

### 2.7 Layer-to-layer angle shift (the UofSC `angle_shift` parameter)

The Kodagali 2023 thesis ch. 4 notation `[fiber angles][placement seq.][angle shift][tow width]` includes an explicit `angle shift` parameter $a_s$ that rotates the placement-string pattern of ply $i+m$ (where $m$ is the number of fiber-angle directions) relative to the placement-string pattern of ply $i$. For $a_s = 0$ the second period of the layup repeats the first period verbatim (the Kok-default behavior used at stage 11). For $a_s \neq 0$ the second period is offset by $a_s$ tape-widths in the laydown plane, shifting where the active and inactive channels land across the panel. The port supports $a_s = 0$ for compatibility with Kok defaults and supports nonzero $a_s$ for compatibility with the UofSC notation.

```
for each ply i:
    period_index = floor((i - 1) / m)   # m is the number of fiber-angle directions
    angle_index_within_period = (i - 1) mod m
    nominal_angle = config.angles[angle_index_within_period]
    placement_origin_offset = period_index * a_s * w
    apply origin offset and placement string to generate ply i tow centerlines
```

Source: Kodagali 2023 thesis ch. 4 (cited in `references/uofsc_pseudowoven.md` §5.1, "placement-sequence notation" row). Nagelsmit 2013 §2.3 does not name this parameter explicitly but its "totally interwoven" pattern (§2.2 last subsection) requires it implicitly: "after step 3 of the original series pattern, a third orientation is added".

### 2.8 Laminate stacking enumeration

The full laydown is the ordered concatenation of plies 1 through $n_p$, each defined by (nominal angle, placement string, tape width, tape spacing, ply thickness, undulation ratio, period-shifted origin from §2.7). The port iterates plies in order and accumulates the solid bodies described in §2.4-2.6. Mid-plane symmetry (the `_S` subscript in Nagelsmit 2013 §2.3.1 notation) is supported as a generator that mirrors the first half of the stack to produce the second half.

### 2.9 Open questions to resolve before implementation

The published papers do not fully specify every geometric primitive. The port must commit to a specific choice for each open question and document that choice in the runner output for reproducibility.

- **OQ-1 (undulation centerline shape).** Kok 2022 Figure 4 schematically depicts a smooth undulation between two straight-tow segments and gives $\varphi(x) = \tan^{-1}(\partial z / \partial x)$ as the local out-of-plane angle. The exact $z(x)$ shape — whether sinusoid, cubic Hermite spline, quintic, parabola, or piecewise linear — is **not** given. The downstream constitutive consequence is bounded because Kok 2022 §3.1 averages the angle over half the undulation length to define $\varphi_{\text{avg}}$, which is a single scalar per undulation region; thus any centerline shape that integrates to the same $\varphi_{\text{avg}}$ produces the same constitutive behavior. The port commits to a **half-period sinusoid** $z(x) = (t_p / 2) \cdot (1 - \cos(\pi x / L_u))$ for $x \in [0, L_u]$ as the canonical choice (it is the simplest $C^1$ bump that integrates analytically to $\varphi_{\text{avg}} = \arctan(2 t_p / L_u) / 2 \approx \arctan(\chi)$ for small $\chi$). The choice is documented in the runner log; alternative shapes can be added later under a config flag without breaking downstream consumers because $\varphi_{\text{avg}}$ is the only thing the constitutive law sees.
- **OQ-2 (in-plane tape edge geometry).** Nagelsmit 2013 §2.3 ("Undulations") notes that the as-manufactured tow edge is "smeared out by the undulation, smoothening the angle of the undulating tow", so the cured laminate has a slightly filleted tow edge rather than a perfectly square slit-tape edge. Kok 2022 §3.1 idealizes the tow as a rectangular cross-section throughout, including in the undulation region. The port follows the Kok 2022 idealization (square cross-section) because (i) the smearing is a sub-mesh-resolution effect at the locked tape widths 6.35-12.7 mm and target mesh sizes 0.06-0.5 mm in the impact zone (`tests/stage_11_PW_mesoscale_direct/spec.md` §3.2), and (ii) introducing a fillet adds a free parameter (fillet radius) that no source pins down. The choice is documented in the port's README.
- **OQ-3 (undulation footprint shape in plan view).** Kok 2022 Figure 4 shows the undulation as a single rectangular footprint (width $w$, length $L_u$) centered on the crossover. For non-orthogonal angles ($\theta_i - \theta_j \neq 90$ deg) the actual physical crossover is a parallelogram, not a rectangle. The port treats the undulation footprint as a parallelogram of width $w$ in the perpendicular direction to the undulating tow's centerline, with length $L_u$ along that centerline; this is the natural extension of Kok 2022 Figure 4 to non-orthogonal layups and is consistent with Nagelsmit 2013 §2.4.1 reporting that "the angle between the interwoven plies does not change this unit cell significantly except for the visual appearance".
- **OQ-4 (resin-pocket constituent treatment).** §2.6 above notes Kok 2022 Figure 4 idealizes the resin-rich region as a mixed unit cell, while `tests/stage_11_PW_mesoscale_direct/spec.md` §2.2.2 treats it as pure-resin LAW1. The port emits the region as its own physical group `RESIN_PLY_<i>`; the choice of constitutive law is downstream.
- **OQ-5 (tape-spacing edge cases).** For tape spacing $s = 1$ (one-over-one), every active channel has exactly one inactive channel to its side and the placement string `1010...` describes one ply orientation. For $s = 2$ the string is `100100...` and for $s = 3$ the string is `10001000...` (the UofSC-Vakili Rad string). Higher values are documented but not exercised by the master plan; the port supports them by formula but tests them only at $s \in \{1, 2, 3\}$.
- **OQ-6 (panel boundary clipping).** Tow centerlines extend across the entire panel polygon. Tows that clip the panel boundary are truncated. The port's truncation is hard-edged (the tow body is cut by the panel boundary plane). Realistic AFP edge effects (head retraction, partial tow drops) are out of scope per §1.3.

Each open question will be reflected in the configuration schema (§3) so that future refinement of any one question does not require a config-format change.

---

## 3. JSON config schema

### 3.1 Schema overview

The configuration is a single JSON object with three top-level keys: `panel`, `laydown`, and `mesh`. Optional `output` and `validation` keys carry reproducibility metadata.

```jsonc
{
  "panel": {
    "size_x_mm": 25.0,
    "size_y_mm": 25.0,
    "n_plies": 4,
    "polygon_mm": [[ -12.5, -12.5 ], [ 12.5, -12.5 ], [ 12.5, 12.5 ], [ -12.5, 12.5 ]],
    "symmetry": "none"
  },
  "laydown": {
    "fiber_angles_deg": [0.0, 45.0, -45.0, 90.0],
    "placement_sequence": "1010",
    "angle_shift_deg": 0.0,
    "tape_width_mm": 6.35,
    "cured_ply_thickness_mm": 0.18,
    "undulation_ratio": 0.09,
    "tape_spacing": 1
  },
  "mesh": {
    "in_plane_target_mm_impact_zone": 0.06,
    "in_plane_target_mm_far_field": 0.06,
    "through_thickness_target_mm": 0.045,
    "graded_zone_radius_mm": 0.0
  },
  "output": {
    "msh_path": "build/panel.msh",
    "orientations_json_path": "build/orientations.json",
    "msh_format": "msh4_ascii"
  }
}
```

Every field is required unless explicitly optional. `panel.polygon_mm` defaults to a centered axis-aligned rectangle of `panel.size_x_mm` by `panel.size_y_mm`; arbitrary convex polygons are accepted but not exercised in this master plan.

### 3.2 Field definitions and constraints

| Field | Type | Constraint | Source / consumer |
|---|---|---|---|
| `panel.size_x_mm`, `panel.size_y_mm` | float | > 0, $\le 1200$ for stage 16's 1143 mm panel | stages 11 (25 mm), 16 (200 mm reduced section, 1143 × 838 mm full) |
| `panel.n_plies` | int | $\ge 1$, even when `symmetry == "midplane"` | stages 11 (4), 16 (24) |
| `panel.polygon_mm` | list of (x,y) | convex, closed implicitly (last point auto-connects to first) | the panel-polygon clipping primitive in §2.3 |
| `panel.symmetry` | "none" or "midplane" | mid-plane mirror generator per §2.8 | stage 16's 24-ply symmetric layup |
| `laydown.fiber_angles_deg` | list of floats | length is the number of fiber-angle directions $m$ | Kodagali 2023 ch. 4 notation; Kok 2022 §2 (cross-ply $m=2$, quasi-iso $m=4$) |
| `laydown.placement_sequence` | string of `0` and `1` | length must be a multiple of $m \cdot (s+1)$; `1` = active channel, `0` = inactive | UofSC notation ; example: `"10001000"` for 4-angle `s=3` (Vakili Rad 2020 P1-TWT) |
| `laydown.angle_shift_deg` | float | per-period origin offset in tape-widths, see §2.7 | UofSC notation |
| `laydown.tape_width_mm` | float OR list | single value or per-ply list, length must equal `n_plies` if list | locked range $\{6.35, 10, 12.7\}$ from `master_plan.md` §2 |
| `laydown.cured_ply_thickness_mm` | float | locked range $\{0.18, 0.18125, 0.213\}$ | `master_plan.md` §2 |
| `laydown.undulation_ratio` | float | locked range $0.09 \le \chi \le 0.104$ | `master_plan.md` §2; Kok 2022 Figure 4 |
| `laydown.tape_spacing` | int | $\in \{1, 2, 3\}$ | `master_plan.md` §2; Nagelsmit 2013 §2.3 (skip factor) |
| `mesh.in_plane_target_mm_impact_zone` | float | > 0 | stage 11 ≤ 0.06 mm; stage 16 = 0.5 mm |
| `mesh.through_thickness_target_mm` | float | $\le t_p / 3$ (3 elements through tape) | stage 11 §3.2 |
| `mesh.graded_zone_radius_mm` | float | 0 disables grading | stage 16 = 40 mm |
| `output.msh_format` | string | `"msh4_ascii"`, `"msh4_binary"`, or `"msh2_ascii"` | meshio bridge prefers msh4 |

### 3.3 Backwards-compatible Nagelsmit / UofSC notation

The Kodagali 2023 notation `[fiber angles][placement sequence][angle shift][tow width]` and the Nagelsmit 2013 notation `[(θ1/θ2)n×w/(θ3/θ4)n×w]_{kS}` are both supported as *aliased* string forms via an optional `laydown.shorthand` field. When present, `shorthand` is parsed by a small grammar in `parser.py` (§4 below) and overrides the explicit fields. Examples:

- `"shorthand": "[0,45,-45,90][1010][0][6.35]"` is equivalent to the explicit form in §3.1.
- `"shorthand": "[(45/-45)1×0.5/(90/0)1×0.5]_3S"` (Nagelsmit 2013 §2.3.1 example, units: tape-widths in inches expressed as fractions) is equivalent to a 4-angle 6-ply repeat with $s=1$ on 0.5-inch (12.7 mm) tape, repeated 3 times symmetrically.

The Nagelsmit form is offered as a literature-compatibility convenience; the canonical form is the explicit JSON in §3.1.

### 3.4 Worked example configurations

#### 3.4.1 Stage 11 — Kok 2022 quasi-iso 4-ply mesoscale block (small)

```json
{
  "panel": {
    "size_x_mm": 25.0,
    "size_y_mm": 25.0,
    "n_plies": 4,
    "symmetry": "none"
  },
  "laydown": {
    "fiber_angles_deg": [0.0, 45.0, -45.0, 90.0],
    "placement_sequence": "1010",
    "angle_shift_deg": 0.0,
    "tape_width_mm": 6.35,
    "cured_ply_thickness_mm": 0.18,
    "undulation_ratio": 0.09,
    "tape_spacing": 1
  },
  "mesh": {
    "in_plane_target_mm_impact_zone": 0.06,
    "in_plane_target_mm_far_field": 0.06,
    "through_thickness_target_mm": 0.045,
    "graded_zone_radius_mm": 0.0
  },
  "output": {
    "msh_path": "build/stage11_block.msh",
    "orientations_json_path": "build/stage11_orientations.json",
    "msh_format": "msh4_ascii"
  }
}
```

Matches `tests/stage_11_PW_mesoscale_direct/spec.md` §2.2.1 exactly (four-angle quasi-iso, 6.35 mm tape, 25 mm × 25 mm × 0.72 mm block, $s=1$, $\chi=0.09$).

#### 3.4.2 Stage 16 — UofSC Vakili Rad 2020 P1-TWT 24-ply 200 mm reduced section

```json
{
  "panel": {
    "size_x_mm": 200.0,
    "size_y_mm": 200.0,
    "n_plies": 24,
    "symmetry": "midplane"
  },
  "laydown": {
    "fiber_angles_deg": [45.0, 90.0, -45.0, 0.0],
    "placement_sequence": "10001000",
    "angle_shift_deg": 0.0,
    "tape_width_mm": 6.35,
    "cured_ply_thickness_mm": 0.18958,
    "undulation_ratio": 0.09,
    "tape_spacing": 3
  },
  "mesh": {
    "in_plane_target_mm_impact_zone": 0.5,
    "in_plane_target_mm_far_field": 2.0,
    "through_thickness_target_mm": 0.18958,
    "graded_zone_radius_mm": 40.0
  },
  "output": {
    "msh_path": "build/p1_twt_section.msh",
    "orientations_json_path": "build/p1_twt_orientations.json",
    "msh_format": "msh4_binary"
  }
}
```

Matches `tests/stage_16_PW_panel_ballistic/spec.md` §2.3 except for the Vakili Rad cured ply thickness 0.18958 mm (= 4.55 mm / 24) and the placement string `10001000` for $s=3$ (the documented Kodagali 2023 ch. 4 P1-TWT laydown).

#### 3.4.3 Edinburgh Kok 2022 cross-ply 32-ply panel (literature cross-check)

```json
{
  "panel": {
    "size_x_mm": 40.0,
    "size_y_mm": 40.0,
    "n_plies": 32,
    "symmetry": "midplane"
  },
  "laydown": {
    "fiber_angles_deg": [0.0, 90.0],
    "placement_sequence": "1010",
    "angle_shift_deg": 0.0,
    "tape_width_mm": 10.0,
    "cured_ply_thickness_mm": 0.213,
    "undulation_ratio": 0.104,
    "tape_spacing": 3
  },
  "mesh": {
    "in_plane_target_mm_impact_zone": 1.5,
    "in_plane_target_mm_far_field": 1.5,
    "through_thickness_target_mm": 0.07,
    "graded_zone_radius_mm": 0.0
  },
  "output": {
    "msh_path": "build/edinburgh_xpap.msh",
    "orientations_json_path": "build/edinburgh_xpap_orientations.json",
    "msh_format": "msh4_ascii"
  }
}
```

This config reproduces the Edinburgh Kok 2022 cross-ply ($[0/90]_{2S}$) 40 mm × 40 mm "approximate RVE" of Kok 2022 §2 with the 10 mm hand-slit tape and 3-tape-width gap also from Kok 2022 §2. It is used in the validation strategy (§7c) as a literature-cross-check configuration.

---

## 4. Module structure

The Python package lives under `src/kok_geom/` and is laid out as small, single-responsibility modules. No module exceeds ~400 lines of code (the largest is the OCC primitive layer, §4.7). The package targets Python 3.11+, depends on `gmsh` (LGPL but explicitly carved out as a system dependency the way `numpy` is — the geometry is GMSH-driven, not GMSH-derived), `numpy`, `meshio`, `pydantic` (for JSON schema validation), and `pytest` (for tests).

```
src/kok_geom/
├── __init__.py                  # public API surface; exposes generate_mesh(config_path) -> (msh_path, orientations_json_path)
├── cli.py                       # argparse-based CLI entry point; `python -m kok_geom --config <path> --out <dir>`
├── parser.py                    # JSON schema validation; Nagelsmit / UofSC shorthand grammar; produces typed config object
├── tow_path.py                  # per-ply tow-centerline generation (§2.3); pure-Python, no GMSH dependency
├── tow_section.py               # straight-tow rectangular cross-section primitive (§2.4); pure Python, returns coordinate arrays
├── undulation.py                # undulation primitive: footprint + half-period-sinusoid centerline + phi_avg integral (§2.5, §2.9 OQ-1)
├── ply_assembler.py             # per-ply union of tow bodies, undulation bodies, and resin pockets (§2.6)
├── occ_kernel.py                # GMSH-OCC adapter layer; wraps gmsh.model.occ.* calls behind a small typed interface
├── fragment.py                  # boolean-cascade strategy: tile, fragment, fuse, deduplicate (§5.4)
├── mesher.py                    # mesh-size field, Mesh.OptimizeNetgen, output writer
├── orientations.py              # builds the per-physical-group orientation map and writes orientations.json (§6)
├── config.py                    # pydantic models for the JSON schema in §3
└── __main__.py                  # delegates to cli.main()

tests/
├── unit/
│   ├── test_parser.py           # parser round-trip + shorthand grammar (§10)
│   ├── test_tow_path.py         # tow-path geometry: counts and angles for s in {1,2,3} (§10)
│   ├── test_undulation.py       # undulation centerline equation, phi_avg integral
│   ├── test_orientations.py     # orientation-vector unit-norm and per-group cardinality
│   └── test_msh_loadable.py     # exported .msh round-trips through meshio (§10)
└── integration/
    ├── test_kok2022_xpap.py     # 40 × 40 mm cross-ply config (§3.4.3) builds and meshes
    └── test_stage11_smoke.py    # the §3.4.1 config builds, meshes, and orientations.json validates
```

| Module | Responsibility (one line) |
|---|---|
| `parser.py` | Validate JSON config; parse Nagelsmit / UofSC shorthand; produce typed config |
| `tow_path.py` | Generate per-ply straight tow centerlines from angle, placement string, and tape spacing |
| `tow_section.py` | Build rectangular tow cross-section and the swept-solid coordinate description |
| `undulation.py` | Locate tow-tow crossings; build undulation footprints; compute φ_avg per region |
| `ply_assembler.py` | Per-ply union of straight-tow + undulation + resin-pocket; emits typed body list |
| `occ_kernel.py` | Thin wrapper over `gmsh.model.occ` so tests can mock the OCC layer |
| `fragment.py` | Tile + fragment + fuse strategy to keep boolean-cascade time bounded |
| `mesher.py` | Apply mesh-size fields and write the GMSH .msh output |
| `orientations.py` | Build {physical group → fiber direction unit vector}; write sidecar JSON |
| `cli.py` | argparse entry; calls `parser → tow_path → undulation → ply_assembler → fragment → mesher → orientations` |
| `config.py` | pydantic models for the schema of §3 |

---

## 5. GMSH-OCC implementation strategy

### 5.1 Primitive-to-API mapping

Each primitive named in §2 maps to a small set of `gmsh.model.occ.*` calls. The mapping is the explicit subject of `occ_kernel.py` (§4.7).

| Primitive (§2 reference) | GMSH-OCC API call(s) | Notes |
|---|---|---|
| Panel bounding box (the laydown domain) | `gmsh.model.occ.addBox(x0, y0, z0, dx, dy, dz)` | One per laminate; sets the outer clip envelope used by the resin-pocket subtraction in §2.6. |
| Straight tow body (§2.4) | `gmsh.model.occ.addBox(...)` then `gmsh.model.occ.rotate(...)` about the panel-Z axis by $\theta_i$ | A box of size $w \times L_{\text{tow}} \times t_p$ rotated to the tow direction is the cleanest OCC primitive for a straight rectangular tow. Cheaper than `addPipe` on a polyline. |
| Tow centerline polyline (intermediate) | `gmsh.model.occ.addLine` and `gmsh.model.occ.addWire` | Used only as a debug/check primitive in `tow_path.py`. |
| Undulation footprint (§2.5) in plan view | `gmsh.model.occ.addRectangle` (parallelogram via four points + `addPlaneSurface`) | The parallelogram footprint (§2.9 OQ-3) is built from four explicit points. |
| Undulation through-thickness body | Either `gmsh.model.occ.addPrism` (sweeping the parallelogram footprint along $Z$ by $t_p$) or `gmsh.model.occ.addThruSections` (lofting between bottom and top parallelograms with a $Z$ offset). The half-period sinusoid centerline is approximated as a polyline with at least 8 segments per half-undulation per stage 11 §3.2 mesh resolution. | The lofted surface is what gives the undulation its tilted-axis character. The constitutive layer reads the average tilt φ_avg, which is computed analytically by `undulation.py` and stored in the orientation sidecar — it does not have to be inferred from the mesh. |
| Resin pocket (§2.6) | `gmsh.model.occ.cut(panel_volume, union_of_tow_and_undulation_volumes)` | The remaining ply volume after cutting out occupied regions. |
| Tow-tow boolean fragmentation (resolving overlaps at undulation crossings) | `gmsh.model.occ.fragment(all_tow_bodies, all_undulation_bodies)` | The fragment operator is the OCC primitive that splits and re-glues bodies so their interfaces are conformal. This is the heaviest call in the pipeline. |
| Per-tow physical group | `gmsh.model.addPhysicalGroup(dim=3, tags=[volume_tag], name="TOW_PLY_<i>_TAG_<j>")` | One physical group per tow, named by ply index and tow index within that ply. |
| Per-undulation physical group | `gmsh.model.addPhysicalGroup(dim=3, tags=[volume_tag], name="UNDUL_PLY_<i>_PLY_<j>_TAG_<k>")` | One per undulation crossover, named by the two crossing ply indices and a within-ply crossover index. |
| Per-resin-pocket physical group | `gmsh.model.addPhysicalGroup(dim=3, tags=[volume_tag], name="RESIN_PLY_<i>")` | One per ply; the pocket is unioned within a ply to keep the group count tractable. |
| Mesh-size field (in-plane) | `gmsh.model.mesh.field.add("Distance", id)` then `gmsh.model.mesh.field.add("Threshold", id2)` referencing the impact-point as the distance source; `gmsh.model.mesh.setSize(...)` for fine control at tow boundaries. | Implements the impact-zone-vs-far-field grading from `tests/stage_16_PW_panel_ballistic/spec.md` §3.1. |
| Final synchronization | `gmsh.model.occ.synchronize()` then `gmsh.model.mesh.generate(3)` | Synchronize is required before any meshing or physical-group assignment can see OCC entities. |
| Cleanup of duplicate entities | `gmsh.model.occ.removeAllDuplicates()` | Called after every fragment cascade tile to keep the entity table small. |
| Export | `gmsh.write(output_msh_path)` | GMSH chooses the format from the file extension; `msh_format` config field maps to `"Mesh.MshFileVersion"` and `"Mesh.Binary"` options before `gmsh.write`. |
| Per-element orientation reference vector | NOT a GMSH primitive. The port writes the orientation as a sidecar JSON file (§6.2) keyed on physical-group name. The downstream `inp2rad` step then translates to OpenRadioss `/SKEW/FIX` per `tests/stage_11_PW_mesoscale_direct/spec.md` §2.2.2. | GMSH does not natively carry per-element material orientation; this is the responsibility of the constitutive layer downstream. |

### 5.2 Why OCC and not the GMSH built-in geo kernel

The GMSH built-in `geo` kernel performs only constructive solid geometry (CSG) on box-like primitives without true boolean operators. AP-PLY laydown is fundamentally a boolean problem: tows from different plies overlap in the in-plane domain at every crossing, and the *fragment* primitive that splits and re-glues those overlaps is OCC's `BOPAlgo_Builder` underneath the GMSH wrapper (`OpenCASCADE2026`). Using the geo kernel would force a manual partition of every crossover, which is exactly what the Kok upstream does in Shapely 2-D and projects through Abaqus's CSG. The OCC kernel performs the same operation directly in 3-D inside GMSH, which is why the master plan §4 names "GMSH-OCC" specifically.

### 5.3 The boolean-cascade scaling concern

The Vakili Rad 2020 P1-TWT 1143 × 838 mm panel at $s = 3$ tape spacing has approximately
$$
N_{\text{tow}} \approx \frac{1143 \cdot 838}{(s+1) \cdot w \cdot w \cdot \tan(\pi / m_{\text{eff}})} \cdot n_p \approx \mathcal{O}(10^4)
$$
tows over 24 plies. Each tow generates 4 to 12 undulations with crossing tows in adjacent plies. The total number of OCC bodies before fragmentation is therefore $\mathcal{O}(10^5)$, and the post-fragment body count rises further as fragments split bodies along intersection curves. OCC's `BOPAlgo_Builder` time complexity is documented (`OpenCASCADE2026`) as super-linear in the number of input shapes, with practical observations of $\mathcal{O}(N \log N)$ to $\mathcal{O}(N^{1.5})$ depending on the spatial distribution of intersections. A single boolean call on $10^5$ bodies is impractical (multi-hour wall clock and likely to trigger OCC's internal precision-ε failure modes).

The reduced 200 × 200 mm section used by stage 16 cuts the tow count by a factor of $24$ to $\mathcal{O}(10^3)$ tows, which is at the upper bound of what a single-cascade boolean can handle in a few minutes. Stage 11's 25 × 25 mm block has $\mathcal{O}(10^2)$ tows — a single boolean call is fine.

### 5.4 Strategy: panel tiling with intra-tile fragmentation

The port adopts a deliberate tiling strategy.

1. **Generate a unit cell.** Find the smallest rectangular tile that, when array-replicated in $X$ and $Y$, reproduces the full panel laydown. The unit cell size depends on $w$, $s$, and the angle set; for a quasi-iso four-angle laydown with $s = 1$ on $w = 6.35$ mm tape it is $4w \times 4w \approx 25 \times 25$ mm — exactly the stage 11 block size. For the Vakili Rad P1-TWT $s = 3$ it is $8w \times 8w \approx 50 \times 50$ mm.
2. **Build the unit cell as a single OCC fragment cascade.** All tows, undulations, and resin pockets inside the unit cell are unioned into a single fragmented body. Boolean cost is $\mathcal{O}(10^2)$ shapes regardless of panel size.
3. **Replicate the unit cell across the panel via OCC `copy` + `translate`.** Each replica is a separate fragmented body that does not need re-fragmentation internally.
4. **Stitch adjacent unit cells.** At the four shared edges of each unit cell, a single fragment-and-fuse along the edge plane resolves the conformal mesh between adjacent tiles. This is $\mathcal{O}(N_{\text{tile}})$ boolean calls each on $\mathcal{O}(10)$ shapes (the boundary tow stubs from each of two adjacent tiles).
5. **Apply the panel-boundary clip.** A final boolean cut against the panel polygon trims any tow stubs that extend beyond the panel edge.
6. **Synchronize and mesh.** A single `gmsh.model.occ.synchronize()` then a single mesh-generate call.

The total OCC effort for the 200 × 200 mm stage 16 reduced section under this strategy is approximately $4 \times 4 = 16$ unit cells × O($10^2$) shapes per fragment, plus $\sim 24$ tile-edge fuses × O($10$) shapes each. This is comfortably under a 10-minute wall clock on a laptop. The full 1143 × 838 mm panel is $\sim 23 \times 17 = 391$ unit cells, which is still tractable (maybe an hour of OCC time) and is the correct strategy if the user later opts into the full panel.

### 5.5 Mesh-size field across the impact-zone-vs-far-field grading

The mesh-size field is not part of the OCC step but is the natural follow-on. The port implements the grading from `tests/stage_16_PW_panel_ballistic/spec.md` §3.1 by:

- Defining a `Distance` field whose source is the impact-centre point.
- Defining a `Threshold` field that maps distance ≤ `mesh.graded_zone_radius_mm` to `mesh.in_plane_target_mm_impact_zone` and distance > `mesh.graded_zone_radius_mm` to `mesh.in_plane_target_mm_far_field`, with a smooth ramp between (default GMSH `Threshold` field uses a linear ramp on a `LcMin`–`LcMax` interval).
- Setting that threshold field as the background mesh size with `gmsh.model.mesh.setBackgroundField`.

For stage 11 there is no grading: in-plane and far-field targets are identical (0.06 mm) and `graded_zone_radius_mm` is zero.

### 5.6 Element type

GMSH's default 3-D mesher produces TETRA4 (linear tets). The downstream consumer (`tests/stage_11_PW_mesoscale_direct/spec.md` §3.1) prefers HEXA8 on `/PROP/TYPE14` and falls back to TETRA10. The port emits whatever GMSH chooses by default; the `mesher.py` module exposes a config flag `mesh.recombine_3d` that, when true, attempts `gmsh.model.mesh.recombine()` to convert tet meshes into hex-dominant. The flag defaults to false because all-HEXA recombination on the AP-PLY partitioned geometry is documented to fail frequently for non-orthogonal undulation regions; the stage-11 spec §3.3 acknowledges TETRA10 as the documented fallback.

---

## 6. Output format

### 6.1 GMSH `.msh` file

The port emits GMSH MSH version 4 (ASCII or binary, configurable). MSH4 is the GMSH default since 2019 and is fully supported by `meshio` (`Schloemer2024Meshio`). Per-physical-group naming follows §5.1. The mesh contains:

- 3-D elements (TETRA4 by default, optionally hex-dominant).
- 2-D physical surfaces only on the panel boundary (for downstream BC group definition).
- 1-D and 0-D entities only as required by the OCC dimension hierarchy (no explicit 1-D physical groups are emitted).

The `.msh → .inp` conversion that stages 11 and 16 already expect (`tests/stage_11_PW_mesoscale_direct/spec.md` §2.3 path B) is performed by `meshio convert panel.msh panel.inp` with no port-specific handling required. `meshio` preserves physical-group names as Abaqus `*ELSET` and `*NSET` blocks, which is exactly what `inp2rad` parses to produce `/PROP` references downstream.

### 6.2 Orientation sidecar `orientations.json`

Per-element fiber orientation is not natively supported by either the GMSH MSH4 format or the Abaqus INP `*ORIENTATION` round-trip through `meshio`. The port emits a sidecar JSON file mapping physical-group name to a unit vector (in the global $X$, $Y$, $Z$ frame) that is the local 1-axis (fiber direction) for that group's elements. The downstream OpenRadioss `/SKEW/FIX` block is built by the runner from this sidecar; the `inp2rad` step does not need to round-trip orientation through the INP file.

```jsonc
{
  "schema_version": 1,
  "panel_config_hash": "<sha256 of the input config JSON>",
  "groups": [
    {
      "name": "TOW_PLY_1_TAG_0",
      "kind": "straight_tow",
      "ply_index": 1,
      "nominal_angle_deg": 0.0,
      "fiber_direction_unit_vector": [1.0, 0.0, 0.0],
      "transverse_in_plane_unit_vector": [0.0, 1.0, 0.0],
      "through_thickness_unit_vector": [0.0, 0.0, 1.0]
    },
    {
      "name": "UNDUL_PLY_1_PLY_2_TAG_0",
      "kind": "undulation",
      "ply_indices": [1, 2],
      "nominal_angles_deg": [0.0, 45.0],
      "phi_avg_deg": 10.2,
      "fiber_direction_unit_vector": [0.9842, 0.0, 0.1771],
      "transverse_in_plane_unit_vector": [0.0, 1.0, 0.0],
      "through_thickness_unit_vector": [-0.1771, 0.0, 0.9842]
    },
    {
      "name": "RESIN_PLY_1",
      "kind": "resin_pocket",
      "ply_index": 1,
      "isotropic": true
    }
  ]
}
```

The `transverse_in_plane_unit_vector` and `through_thickness_unit_vector` are derived from the fiber direction by the right-hand rule (perpendicular within the laminate plane, then $Z$-component as the third axis). For undulation groups, the φ_avg field carries the average out-of-plane angle from Kok 2022 Figure 4 directly, so downstream consumers can reconstruct the rotated frame without re-running the integral.

The runner script in `tests/stage_11_PW_mesoscale_direct/runner.py` already templates one `/SKEW/FIX` per orientation (§6.1 of that spec); the orientation sidecar is consumed in `runner.py`'s mesh-include generation step.

---

## 7. Validation strategy

Validation proceeds in three layers, each with a binary pass/fail decision and an explicit honesty caveat.

### 7.1 Layer A — geometry sanity checks

Cheap, deterministic, run on every build. All three checks must pass before the mesh is written.

1. **Total tow volume vs. nominal.** For each ply, compute the analytic straight-tow volume per the placement string,
   $V_{\text{tow,nominal}} = N_{\text{active}} \cdot w \cdot L_{\text{tow}} \cdot t_p$
   where $N_{\text{active}}$ is the count of `1`s in the placement string per period times the number of periods across the panel, and $L_{\text{tow}}$ is the tow length after panel-boundary clipping. Sum the OCC body volumes of every tow body in the ply (`gmsh.model.occ.getMass(dim=3, tag)` returns volume for unit-density solids). The two values must agree to within 1 percent.
2. **Resin pocket plus tow plus undulation = ply nominal.** For each ply, the sum of the OCC volumes of every body assigned to that ply (straight-tow plus undulation plus resin pocket) must equal the ply's nominal volume $V_{\text{ply}} = A_{\text{panel}} \cdot t_p$ to within 1 percent. Failure indicates a boolean leak.
3. **Tow count check.** The number of physical groups named `TOW_PLY_<i>_*` must match the analytic count from the placement string and panel size. Failure indicates a fragmentation error or a rotation off-by-one.

These three together are necessary but not sufficient. They catch arithmetic and boolean errors; they do not catch a wrong undulation amplitude or a wrong angle.

### 7.2 Layer B — visual cross-section overlay

For each canonical configuration (Kok 2022 cross-ply, Kok 2022 quasi-iso, Vakili Rad P1-TWT), render the unit cell with PyVista headless (`Sullivan2019PyVista`) at the same orientation and crop as Kok 2022 Figure 4 (or Nagelsmit 2013 Figure 2.5 for the AS4/8552 thermoset undulation micrograph). The render is saved to `figures/validation/`. A reviewer overlays the render against the published figure and confirms qualitative match.

This is **explicitly qualitative**. The published figures are schematics (Kok 2022 Figure 4) or single micrographs at one as-manufactured location (Nagelsmit 2013 Figure 2.5). The match quality is "right region types in roughly the right places" rather than pixel-level. The deliverable is a side-by-side image embedded in the port's README. If a future quantitative metric (Hausdorff distance to a digitized published outline) is desired it can be added; nothing in the master plan requires it.

### 7.3 Layer C — quantitative homogenized stiffness (the binding criterion)

The binding validation is the stage 11 acceptance criterion already in the master plan: stage 11 must reproduce Kok 2022 Table 4 numerical $E_x \approx 53.3$ GPa for the quasi-iso AP-PLY layup within a 10 percent tolerance. The port itself is not a solver — it produces geometry; the geometry feeds stage 11; stage 11 produces $E_x^{\text{FEM}}$; the comparison is `tests/stage_11_PW_mesoscale_direct/spec.md` §7.2:

$$
\frac{|E_x^{\text{FEM}} - E_x^{\text{Kok2022}}|}{E_x^{\text{Kok2022}}} \le 0.10
$$

If stage 11 fails this criterion, the port is wrong somewhere and the milestones (§8) regress to the failure point. The port owner is responsible for diagnosing whether the failure traces to (i) the geometry layer (this port), (ii) the material card, (iii) the BC choice, or (iv) the inp2rad orientation translation. The Layer A and Layer B checks above are designed to localize geometry-layer bugs before stage 11 runs.

The 10 percent window is documented in `tests/stage_11_PW_mesoscale_direct/spec.md` §7.2 to absorb (i) the Hill uniform-displacement bias for a 25 mm specimen with a 12.7 mm interlace period, (ii) the IM7/8552-vs-VTC401 architectural correction, and (iii) the mesh-discretization residual. None of those three is a port-side error if it fails. But a 30 percent error would be a port-side error — the geometry would be visibly wrong.

### 7.4 Optional: Edinburgh cross-ply cross-check

The §3.4.3 Edinburgh cross-ply config reproduces Kok 2022's other published configuration ($[0/90]_{2S}$). If the user wants a second binding numerical comparison, the port should produce a mesh that, when run through the same stage 11 pipeline at a different stacking sequence, recovers Kok 2022 Table 4's cross-ply numerical $E_x$ (which the table reports separately from the QI value). This is a stage-11-with-different-config experiment, not a port-side change. It is listed here for completeness; the master plan does not require it.

---

## 8. Milestones with explicit gates

Five milestones, each with a binary "done when" criterion. Total budget 4 weeks of focused work.

### 8.1 M1 — Single-tow path and cross-section (week 1)

**Scope.** Implement `parser.py`, `tow_path.py`, `tow_section.py`, the `occ_kernel.py` adapter, and a minimal `cli.py` that takes a single-ply, single-angle, $s = 1$ config and emits a GMSH `.msh` with a single rectangular tow swept along a straight centerline at the configured angle.

**Done when.** The §10.1 unit tests `test_parser.py`, `test_tow_path.py`, and a smoke test that builds a single tow and asserts (i) the OCC volume equals the analytic $w \cdot L_{\text{tow}} \cdot t_p$ to within 0.1 percent, (ii) `meshio` round-trips the `.msh` to `.inp` and back without error, all pass.

### 8.2 M2 — Single-ply array (week 2)

**Scope.** Implement `ply_assembler.py` for the straight-tow-and-resin-pocket subset (no undulations yet). The output is a single ply of straight tows with the placement string applied and the resin pockets correctly subtracted.

**Done when.** The Layer A geometry sanity checks (§7.1) all three pass on a single-ply $s = 1$ and $s = 3$ configuration; the unit test `test_orientations.py` confirms one orientation per tow physical group; a visual render of the ply matches the schematic in Kok 2022 Figure 1 ("Pass 1: 0°") qualitatively.

### 8.3 M3 — Multi-ply assembly with undulations and orientation tags (week 3)

**Scope.** Implement `undulation.py` (the §2.5 primitive plus §2.9 OQ-1 sinusoid centerline), extend `ply_assembler.py` to multi-ply with mid-plane symmetry, integrate `orientations.py` for the sidecar JSON.

**Done when.** A 4-ply $[0/45/-45/90]$ quasi-iso build (the §3.4.1 stage-11 config) completes Layer A sanity checks; the visual cross-section overlay (Layer B, §7.2) qualitatively matches Kok 2022 Figure 4; `orientations.json` validates against the schema in §6.2 and the unit-norm check in `test_orientations.py` passes.

### 8.4 M4 — Tile-and-fragment scaling and full mesh export (week 4 first half)

**Scope.** Implement `fragment.py` (the §5.4 tile-and-fuse strategy), wire it through the multi-ply assembler. Run on the §3.4.2 stage-16 reduced section to confirm the boolean cascade completes in under 10 minutes on a laptop.

**Done when.** The §3.4.2 P1-TWT 200 × 200 mm × 24-ply mesh builds in less than 10 minutes; the resulting `.msh` is loadable by `meshio` and `gmsh -open`; element count is within the master plan's stated 700k order-of-magnitude budget (`tests/stage_16_PW_panel_ballistic/spec.md` §3.1).

### 8.5 M5 — Stage 11 validation (week 4 second half)

**Scope.** Run the full §3.4.1 stage-11 pipeline end-to-end: port produces mesh, mesh feeds stage-11 runner, runner produces $E_x^{\text{FEM}}$. Compare to Kok 2022 Table 4 quasi-iso $E_x \approx 53.3$ GPa.

**Done when.** Stage 11's pass criterion (`tests/stage_11_PW_mesoscale_direct/spec.md` §7.2) is satisfied: $|E_x^{\text{FEM}} - 53.3| / 53.3 \le 0.10$. Same for $E_y$. $G_{xy}$ within 10 percent of CLT-corrected estimate.

The only gate that *blocks* the rest of the project is M5. M1-M4 are internal milestones. If M5 fails by a wide margin (>30 percent error), the cause is investigated layer by layer, beginning with the sanity checks of §7.1 — most likely a sign error, a units error, or a wrong placement-string rotation in the angle-shift logic of §2.7.

### 8.6 Schedule risk

The 4-week budget assumes: (i) one full-time owner, (ii) GMSH-OCC is installed and working before M1 begins, (iii) the user has already passed stage 1 (the toolchain smoke test from `plan/master_plan.md` §9, scheduled before this port lands), and (iv) no scope creep into the §1.3 out-of-scope items.

---

## 9. Risks and unknowns

### 9.1 Honest known unknowns

- **R1. Undulation profile equation choice.** Kok 2022 §3.1 commits to $\varphi_{\text{avg}}$ via an integral but does not pin the centerline shape (§2.9 OQ-1). The port's choice (half-period sinusoid) is the simplest credible one; it may turn out the actual Kok upstream uses a cubic spline or a piecewise linear with different φ_avg integral. The bound on error is small because the constitutive law sees only φ_avg, but this is the principal source of port-vs-upstream geometric divergence and a sceptic could fail to be convinced without an explicit citation. *Mitigation.* Document the choice in the port README; offer a config-flag override if a future paper pins the shape.
- **R2. OCC boolean stability at large tow counts.** §5.3-5.4 commits to a tile-and-fuse strategy. OCC's `BOPAlgo_Builder` has documented precision-ε failures (`OpenCASCADE2026`) when shapes intersect at small angles or with sub-millimeter tolerances. The undulation primitive at $\chi = 0.09$ has a ~5-degree out-of-plane angle; at ~0.18 mm ply thickness and 6.35 mm tape width the smallest geometric feature is the half-undulation-length (~1 mm) by ply-thickness intersection (~0.18 mm). Both are well above OCC's default tolerance (1e-7), so this is a soft risk; `removeAllDuplicates` is called after every cascade as the standard mitigation. *Fallback.* If OCC fails on a specific config, the failing tile is rebuilt with explicit `gmsh.model.occ.healShapes` before fragmentation.
- **R3. meshio `.msh → .inp` orientation preservation.** §6 commits to the sidecar-JSON path because per-element orientation does not survive the meshio round-trip (`Schloemer2024Meshio` does not handle Abaqus `*ORIENTATION` blocks at the per-element level). The downstream `inp2rad` step parses `*ELSET` and `*NSET` from the INP correctly but the orientation must be supplied separately. *Mitigation.* The runner reads `orientations.json` and writes `/SKEW/FIX` blocks directly; this is documented in `tests/stage_11_PW_mesoscale_direct/spec.md` §6.1 and §6.4.
- **R4. Tape-edge geometry (square vs. filleted).** §2.9 OQ-2 commits to square edges. Nagelsmit 2013 §2.3 documents that as-manufactured edges are smoothed by the undulation (a partial fillet effect). At the locked tape widths of 6.35-12.7 mm and target mesh sizes of 0.06-0.5 mm, the fillet is a sub-mesh feature; downstream constitutive results should be insensitive at the 1 percent level. The square-edge approximation matches the Kok 2022 §3.1 idealization. *Risk.* If a reviewer asks for a fillet, it is a non-trivial geometry change (fillet radius is a free parameter no source pins down). *Mitigation.* Document the choice; offer no fillet config knob; cite Kok 2022 §3.1 as the precedent.
- **R5. The Kok 2022 reported $E_x \approx 53.3$ GPa was obtained on VTC401, not IM7/8552.** Stage 11 substitutes IM7/8552 per the consolidation review canonical material card (item A). Stage 11 §7.1 documents the architectural-correction calculation that brings the IM7/8552 prediction onto the VTC401 reference. The port has no role in this correction — it is purely material — but the validation criterion in §7.3 of *this* plan inherits the correction. *Risk.* If the architectural correction is wrong (e.g., an error in the rule-of-mixtures step), stage 11 will fail by a margin that looks like a port bug. *Mitigation.* Run the optional Edinburgh cross-check (§7.4) where the material is VTC401 throughout; if that comparison passes, the architectural correction is the source of any IM7/8552 mismatch and not the port.
- **R6. Master plan scope drift.** If the user adds a curved-tool draping requirement (the §3 deferred side-track), the port's flat-laminate assumption breaks. *Mitigation.* The port's primitives (§2) are flat-only by design; a draping extension is a separate feature. The configuration schema (§3) does not have a "tool surface" field, so adding one is a backwards-compatible extension that does not invalidate existing configs.

### 9.2 Known knowns (not risks)

- The port is a one-shot build: each config produces one mesh, written to disk, consumed by stages 11/16 verbatim. There is no streaming, no incremental update, no checkpoint/resume. Wall-clock budget per build is dominated by the OCC fragment cascade (§5.3).
- The license is MIT. No upstream LGPL-2.1 source has been read or transcribed.
- The geometry parameters are locked by `master_plan.md` §2 to the small set listed in §1.2 above. The port supports only those values within the schema range; extending the range is straightforward but not in scope.
- The port does not need to call OpenRadioss directly; it ends at the `.msh` plus `orientations.json` boundary. The runner scripts in `tests/stage_11_*/runner.py` and `tests/stage_16_*/runner.py` consume the port's output.

---

## 10. Test fixtures

A small, focused pytest suite that should pass before stage 11 runs. Each test is independent; no test depends on a real OpenRadioss install.

### 10.1 Unit tests

**T1 — `test_parser.py::test_roundtrip`.** Load each of the §3.4 example configs from JSON, validate against the pydantic schema, serialize back to JSON, and assert the two JSON dictionaries are byte-identical (after deterministic key ordering). This catches schema drift.

**T2 — `test_parser.py::test_shorthand`.** Parse the Nagelsmit 2013 §2.3.1 shorthand `[(45/-45)1×0.5/(90/0)1×0.5]_3S` and assert it expands to the same internal config object as the equivalent explicit JSON. Likewise for the Kodagali 2023 ch. 4 shorthand `[0,45,-45,90][1010][0][6.35]`. Failure indicates a grammar error.

**T3 — `test_tow_path.py::test_single_tow_volume`.** Build a single-ply, single-tow config (1 mm × 25 mm × 0.18 mm tow at 0 deg in a 25 × 25 × 0.18 mm panel), call the OCC pipeline, and assert the OCC mass equals $1.0 \cdot 25.0 \cdot 0.18 = 4.5$ mm³ to within 0.1 percent. This catches sign errors in the tow-section sweep.

**T4 — `test_tow_path.py::test_interlace_pattern_s1_s2_s3`.** For each $s \in \{1, 2, 3\}$, build a single-ply, single-angle, full-panel config and assert the count of `TOW_PLY_1_*` physical groups equals the analytic prediction $\lfloor L_{\text{panel}} / ((s+1) \cdot w) \rfloor$. This catches off-by-one errors in the placement-string expansion.

**T5 — `test_tow_path.py::test_angle_shift`.** With `angle_shift_deg = 1.0` (one tape-width shift per period) and a 4-angle 8-ply layup, assert that ply 5's tow centerlines are translated by exactly $w$ relative to ply 1's tow centerlines along the in-plane perpendicular axis. This catches errors in §2.7's period-shift logic.

**T6 — `test_undulation.py::test_phi_avg_integral`.** Numerically integrate the chosen sinusoid $\varphi(x) = \tan^{-1}(\partial z / \partial x)$ for a unit case ($t_p = 0.18$, $L_u = 2.0$, $\chi = 0.09$) and assert the analytic $\varphi_{\text{avg}}$ matches a SymPy-derived closed-form to within 1e-6. This catches errors in the OQ-1 commitment.

**T7 — `test_undulation.py::test_undulation_count`.** For the §3.4.1 4-ply config, assert the number of `UNDUL_*` physical groups equals the analytic crossover count (every active tow on ply $j$ that crosses an active tow on ply $j-1$ where $j \ge 2$). For the simple $[0/45/-45/90]$ $s=1$ case the analytic count is computable by hand and is hard-coded into the test.

**T8 — `test_orientations.py::test_unit_norms`.** Build the §3.4.1 config end-to-end, parse the resulting `orientations.json`, and assert each `fiber_direction_unit_vector` has $|\vec{v}| = 1$ within 1e-9. Catches numerical drift in the angle-to-vector conversion.

**T9 — `test_orientations.py::test_one_orientation_per_group`.** Assert that the set of physical-group names in the `.msh` is exactly the set of `groups[*].name` in `orientations.json`. Catches set-mismatch errors.

**T10 — `test_msh_loadable.py::test_meshio_roundtrip`.** Build the §3.4.1 config, load the `.msh` with `meshio`, write to `.inp`, load the `.inp` back, and assert the element count and the set of cell-block names round-trip exactly. This catches MSH4-vs-MSH2 format incompatibilities.

### 10.2 Integration tests

**I1 — `test_kok2022_xpap.py::test_full_build`.** Build the §3.4.3 Edinburgh cross-ply 32-ply config end-to-end. Assert: (a) build completes in < 5 minutes on the test runner, (b) all Layer A sanity checks (§7.1) pass, (c) `.msh` is loadable by meshio. This is the pre-stage-11 integration smoke test on a configuration with abundant published data.

**I2 — `test_stage11_smoke.py::test_full_build`.** Build the §3.4.1 stage-11 config end-to-end. Assert all Layer A sanity checks pass and the `.msh` is loadable. This is the immediate prerequisite for the M5 milestone.

### 10.3 Test execution order

The tests are independent but a sensible execution order is:
1. T1, T2 (parser layer; cheapest, fail-fast on config errors).
2. T3, T4, T5 (tow-path layer).
3. T6, T7 (undulation layer).
4. T8, T9 (orientation layer).
5. T10 (output format).
6. I1, I2 (integration).

Total wall-clock budget for the full suite is < 15 minutes on a laptop. Stage 11's `runner.py` checks all unit tests pass before invoking the port.

---

## Appendix A — Citation cross-reference

This plan cites only items that resolve cleanly against the source PDFs verified in this session, plus items already in the project's bib files. The verified primary sources are:

- **Nagelsmit 2013 PhD thesis** (`Nagelsmit2013Thesis` in `references/delft_refs.bib`, 167 pages, 76 MB PDF). Algorithmic content for AP-PLY laydown is in **Chapter 2 (Sections 2.2 and 2.3)**, NOT §3.4 as the brief stated. Section pointers in this plan are corrected.
- **Kok 2022 *Composites Part A*** (`Kok2022Tensile` in `references/delft_refs.bib`, 29-page accepted manuscript, CC-BY). The geometry idealization (straight-tow / undulation / resin-rich, undulation ratio χ = t/Lu, average angle φ_avg integral) is in **§3.1 ("Microscale model") and Figure 4**, NOT §2 as the brief stated. Section 2 ("Materials and methods") covers the as-manufactured laydown parameters. Section pointers in this plan are corrected.
- **Kok 2023 *Composites Part B*** (`Kok2023Dynamic` in `references/delft_refs.bib`, 12-page open-access PDF). Used here as a secondary source for the three-region idealization (§4.1, Figure 7) which confirms the Kok 2022 geometric definitions.

The Kok upstream README and top-level file list (`Kok2022PreprocRepo` in `references/delft_refs.bib`) was read for orientation only — no source code from that repository has been read or transcribed. The port acknowledges this prior art in its README without inheriting its license.

New bib entries used by this plan that were not already in the project's bib files are catalogued in `references/kok_port_refs.bib`.

---

## Appendix B — Master plan and consolidation review cross-references

This plan satisfies the §11 commitment of `plan/master_plan.md` and resolves item E of `plan/consolidation_review.md`. Specifically:

- `master_plan.md` §11 promise: "A planning agent will produce `plan/kok_port_plan.md` next, with module breakdown, milestone gates, and validation strategy." → Delivered as the document you are reading.
- `master_plan.md` §11 scope: "clean-room AP-PLY tow-laydown for the master-plan's locked geometry parameters". → §1.2 above lists each locked parameter with its master-plan source row.
- `master_plan.md` §11 schedule: "Implementation is scheduled after Phase 3 (stage 1 toolchain smoke test)". → §8.6 above re-affirms this prerequisite.
- `consolidation_review.md` item E recommendation: "path (1) — Port Kok's tow-laydown logic to standalone GMSH + `gmsh.model.occ` Python. Estimated effort: 2-4 weeks of focused work." → §8 above commits to a 4-week budget across five milestones.
- Stage 11 dependence: `tests/stage_11_PW_mesoscale_direct/spec.md` §2.2 and §2.3 expect a JSON-driven Kok preprocessor invocation that produces a partitioned `.inp` plus orientation tags. → §6 above commits to `.msh` plus orientation sidecar; the `meshio` `.msh → .inp` step is a one-liner the runner already templates.
- Stage 16 dependence: `tests/stage_16_PW_panel_ballistic/spec.md` §2.3 and §2.4 expect the same JSON-driven invocation at panel scale. → §3.4.2 above gives the matching config; §5.4 commits to the tile-and-fragment strategy that keeps stage 16's 200 mm reduced section tractable; §9.1 R2 acknowledges the OCC scaling concern.

Both downstream consumers are satisfied by the port's `.msh` plus `orientations.json` output, with no further coupling between the port and the OpenRadioss layer. The port's success criterion is bound by stage 11's 10 percent stiffness-comparison window per §7.3 of this plan; the project's terminal success criterion (stage 16's 7 percent V50 window) inherits the port's correctness through stage 11.
