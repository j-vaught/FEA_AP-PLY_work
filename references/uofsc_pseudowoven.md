# UofSC Pseudo-Woven (Meso-Architectured) Composites: Research Inventory

Compiled 2026-04-29. Sources: Web searches, ScholarCommons (UofSC ETD), TU Delft Repository, University of Edinburgh Research Explorer, ResearchGate metadata, SAMPE digital library index. Primary text extracted from publicly hosted PDFs (Kodagali 2023 dissertation; Vakili Rad 2020 thesis; Zheng & Kassapoglou 2024; Kok et al. 2023).

> Style note: Where a fact is unverified or could not be confirmed against a primary source, the line is prefixed `UNVERIFIED:`. Do not propagate any unverified line into a publication without re-checking.

---

## 1. Executive summary

The architecture the user is calling "pseudo-woven" is a UofSC-coined branding of the broader concept that TU Delft introduced in 2010-2013 as **AP-PLY (Advanced Placed Ply)**. AP-PLY was originated by C. Kassapoglou (TU Delft, 2010 ISCM talk) and developed in M.H. Nagelsmit's TU Delft PhD (2013); the original SAMPE Journal paper "AP-PLY: A new fibre placement architecture for fabric replacement" appeared in 2011. At UofSC, the same idea was first published in 2017 by Michel J.L. van Tooren and Subramani Sockalingam under the term **"clutch laminates"** (ASC 2017; *Composite Structures* 2019), then re-branded by the Sockalingam group as **"pseudo-woven (PW)"** and **"meso-architectured composites (MAC)"** beginning with C. Vakili Rad's 2020 MSc thesis. Production has occurred on the McNair Aerospace Center's gantry-based **LYNX AFP machine** using 1/4 in. (6.35 mm) wide carbon/epoxy slit tape (T800-SC-24K/P2362W and Hexcel IM7G/8552-1), cured in autoclave at 350 deg F / 90 psi (177 deg C / 0.62 MPa). Headline experimental dataset includes warpage, uniaxial tension, low-velocity impact (ASTM D7136, 25-55 J), high-velocity impact (ASTM D8101, 250-400 ft/s, performed at NASA Glenn / Langley with NASA co-authors), open-hole tension/compression, and compression-after-impact. Independent FEM and dynamic-impact work on AP-PLY is being done at TU Delft (Kassapoglou group, cohesive-zone modeling) and the University of Edinburgh (Martinez-Hergueta group, multiscale tow-wise modeling). No public open-data repository (GitHub, Zenodo, Mendeley Data) for the UofSC layups, G-code, or scan data has been located as of the date above.

---

## 2. Definitive name(s)

- **At UofSC (Sockalingam group):** "pseudo-woven (PW)" laminates and "meso-architectured composites (MAC)". MAC is used in the more recent (2023+) work; PW is used in the 2019-2020 papers.
- **Earlier UofSC term (van Tooren / Sockalingam, 2017-2019):** "clutch laminates" (synonymous with what later became PW/MAC).
- **At TU Delft (originator, Kassapoglou and Nagelsmit, 2010-2013):** "AP-PLY" (Advanced Placed Ply). AP-PLY is the dominant term in the EU literature (Edinburgh, Delft, Oxford, JRC).
- **Other observed phrasings:** "non-conventional AFP composites" (Edinburgh 2021 short-beam-shear paper); "tow-wise composites" (when discussed in the modeling literature).

The UofSC and Delft variants share the same core idea (interlace tows in a tape-laid pattern to mimic a woven fabric) but differ in detail: Delft typically describes it as one staggered tow-set followed by an angled, gap-filling second set; UofSC formalizes it via an `[fiber angles][placement sequence][angle shift][tow width]` notation for arbitrary AFP courses with active/inactive channels.

---

## 3. Principal investigators, students, and collaborators

### UofSC core team

| Person | Role / Affiliation | Topic |
| --- | --- | --- |
| Subramani ("Mani") Sockalingam | Associate Professor, Mechanical Engineering, UofSC; Director, MSMOMS Group; affiliated with McNair Aerospace Center | Advisor on every UofSC PW/MAC paper since 2017 |
| Michel J.L. van Tooren | (then) Endowed Chair, SmartState Center for Multifunctional Materials and Structures, UofSC | Co-PI on the 2017-2019 "clutch laminate" papers |
| Cyrus Vakili Rad | MSc 2020, Mechanical Engineering, UofSC (advisor: Sockalingam); now at Wisk Aero (per LinkedIn) | First UofSC dissertation on PW; first author of *Composites Part B* 2020 paper |
| Karan Kodagali | PhD 2023, Mechanical Engineering, UofSC (advisor: Sockalingam) | Author of MAC framework, *Composites Part B* 2024 paper, follow-on Composites Communications 2025 |
| Frank D. Thomas | UofSC collaborator (2017-2019) | Co-author of clutch-laminate papers |
| Brandon Seay | UofSC, listed on 2019 *Composite Structures* paper | Manufacturing |
| Ramy Harik | Professor, Mechanical Engineering, UofSC; McNair Center; AFP process expert (now Director of Clemson Composites Center per one secondary source - UNVERIFIED, retains UofSC affiliation in 2020 *Composites Part B*) | Co-author on HVI 2020 paper, AFP path planning/process |
| Zafer Gurdal | Endowed Chair, McNair Center, UofSC | Co-author on the 2020 SAMPE and 2024 *Composites Part B* papers; legacy AP-PLY connection via his TU Delft tenure |
| Eric Miller | UofSC | Co-author on Kodagali et al. 2024 *Composites Part B* |
| Julie Roark, Robert Bass, Emily Kaufmann | UofSC research staff/students (acknowledged in Vakili Rad 2020 thesis) | Manufacturing support |

### NASA collaborators (2020 HVI paper)

| Person | Affiliation |
| --- | --- |
| Duane Revilock | NASA Glenn Research Center |
| Charles Ruggeri | NASA Glenn Research Center |
| (Single-stage gas-gun HVI testing at NASA; one P1-TWT panel was manufactured at NASA, the rest at UofSC McNair) | |

### TU Delft (parallel/originator line)

| Person | Affiliation |
| --- | --- |
| Christos Kassapoglou | Associate Professor, Aerospace Structures, TU Delft (originator of AP-PLY, 2010 ISCM) |
| Marius H. Nagelsmit | TU Delft PhD 2013, primary developer of AP-PLY thesis |
| Weiling Zheng | Northwestern Polytechnical University (China) / TU Delft (Kassapoglou 2024 co-author) |

### University of Edinburgh / Oxford / JRC (modeling)

| Person | Affiliation |
| --- | --- |
| Francisca Martinez-Hergueta | Reader, University of Edinburgh (multiscale modeling, tow-wise FE) |
| Rutger Kok | PhD, Edinburgh (lead author 2021-2024 Edinburgh papers) |
| Antonio Pellegrino | University of Oxford (dynamic experiments) |
| Marco Peroni | EU JRC, Ispra (Split-Hopkinson bar) |
| R. Cuvillo, J. Pernas-Sanchez, J. Artero Guerrero, V. Rodriguez-Garcia, R. Guzman de Villoria | UC3M / IMDEA / collaborators (2024 LVI paper) |

---

## 4. Annotated bibliography

> Format: full citation, then 3-5 sentence summary of contribution. DOIs verified against publisher landing pages. Where text was extracted directly from a downloaded PDF, the citation has a (`* primary text confirmed`) marker.

### 4.1 Foundational AP-PLY (TU Delft / SAMPE)

1. **Nagelsmit, M.H., Kassapoglou, C., & Gurdal, Z. (2011).** *AP-PLY: A new fibre placement architecture for fabric replacement.* SAMPE Journal, 47(2), 36-45.
   The originating paper. Defines AP-PLY as a fibre-placed laminate where parallel tows are deposited with gaps the width of one tow; a second course at an angle nests in those gaps; subsequent courses fill remaining gaps. Reports 15% gain in compression-after-impact (CAI) strength, smaller delamination footprints, and a lower BVID energy threshold relative to baseline angle-ply.
   *Note that Z. Gurdal was at TU Delft at the time; he later moved to UofSC McNair, providing a direct intellectual lineage between the Delft and UofSC programs.*

2. **Nagelsmit, M.H. (2013).** *Fibre Placement Architectures for Improved Damage Tolerance.* Doctoral thesis, TU Delft. https://doi.org/10.4233/uuid:f8c8ad4e-5d52-4bfd-8449-2748bad26673
   The full PhD development of AP-PLY. Includes geometry definition, a parametric description of the staggered/interspaced placement, manufacturing trials, and impact characterization. Foundational reference for every later AP-PLY paper.

### 4.2 UofSC line (clutch laminate -> pseudo-woven -> meso-architectured)

3. **Van Tooren, M.J.L., & Sockalingam, S. (2017).** *Non-conventional Composite Laminates - Clutch Laminates.* Proceedings of the American Society for Composites - 32nd Technical Conference. DOI: 10.12783/asc2017/15199.
   First UofSC publication of the architecture. Calls it "clutch laminates" because the through-thickness undulations create a mechanical clutch effect. Carbon/epoxy slit tape laid on an AFP machine via tow skips. Uses 2D shell-element FE with a sub-cell representation of the heterogeneous lay-up. Notes that prediction of deformation, damage, and impact response remains a significant open problem.

4. **Vakili Rad, C., Thomas, F.D., Seay, B., & Sockalingam, S. (2019).** *Manufacturing and characterization of novel clutch non-conventional fiber-reinforced composite laminates.* Composite Structures, 213, 173-185. https://doi.org/10.1016/j.compstruct.2018.12.080
   Journal-length expansion of the 2017 ASC paper. Manufacturing on the UofSC McNair Center AFP. Demonstrates that the clutch architecture can reduce post-cure warpage in asymmetric layups. 3D shell FE sub-cell model used for residual-stress / warpage prediction.

5. **Vakili Rad, C., Thomas, F., Sockalingam, S., & Gurdal, Z. (2019).** *Low Velocity Impact Response of Hybrid Pseudo-Woven Fiber-Reinforced Composite Laminates.* SAMPE 2019 (Charlotte, NC). Paper 189-TP19-1473.
   First UofSC paper to drop the "clutch" branding in favor of "pseudo-woven." Drop-tower LVI on hybrid PW/UD layups at 30 J and 50 J. Compares damage vs. the [45/90/-45/0]s baseline.

6. **Vakili Rad, C. (2020).** *The Manufacturing and Characterization of Pseudo-Woven Carbon Fiber Composite Architectures for Enhanced Damage Tolerance.* MSc Thesis, University of South Carolina. https://scholarcommons.sc.edu/etd/5708 (`* primary text confirmed`)
   Master's thesis; the most accessible single reference for the as-manufactured geometry and process. Reports 4-ply PW achieves 58% warpage reduction vs. asymmetric baseline. Documents three 24-ply hybrid configurations (P1-TWT, P2-WTW, P3-Control). Slit-tape and machine details (see Section 5 below) are pulled from this document. Acknowledges Boeing, SC Space Grant Consortium, and McNair Aerospace Center funding/support.

7. **Vakili Rad, C., Kodagali, K., Roark, J., Revilock, D., Ruggeri, C., Harik, R., & Sockalingam, S. (2020).** *High velocity impact response of hybridized pseudo-woven carbon fiber composite architectures.* Composites Part B: Engineering, 203, 108478. https://doi.org/10.1016/j.compositesb.2020.108478
   The first peer-reviewed journal paper using "pseudo-woven" terminology from UofSC. ASTM D8101 single-stage gas-gun, 250-400 ft/s. Uses three 24-ply Hexcel IM7G/8552-1 layups: P1-TWT (PW outside, UD inside), P2-WTW (UD outside, PW inside), P3-Control [45/90/-45/0]_3s. P1-TWT is the best ballistic performer. Includes DIC and high-speed video.

8. **Kodagali, K., Vakili Rad, C., Sockalingam, S., Gurdal, Z., & Miller, E. (2024).** *Low velocity impact and compression-after-impact response of hybrid pseudo-woven meso-architectured carbon/epoxy composite laminates manufactured via automated fiber placement.* Composites Part B: Engineering, 271, 111154. https://doi.org/10.1016/j.compositesb.2023.111154
   The flagship LVI/CAI paper. Switches material to Toray T800-SC-24K / P2362W epoxy on McNair LYNX AFP. Up to 37% higher critical delamination load and 26% higher residual compressive strength after 55 J impact compared to QI control. Three MAC sub-laminate variants studied: 4-ply [45/90/-45/0]_MAC, 2-ply [45/-45]_MAC, 2-ply [0/90]_MAC.

9. **Kodagali, K. (2023).** *Studies of Damage Tolerance in Automated Fiber Placement Based Heterogeneous Meso-Architectured Carbon/Epoxy Composite Laminates.* Doctoral dissertation, University of South Carolina. https://scholarcommons.sc.edu/etd/7650 (`* primary text confirmed`)
   The single most complete reference. Defines the AFP notation [angles][placement sequence][angle shift][tow width]. Describes the LYNX AFP machine at McNair, lists material specs, manufactures 1143 mm x 838.2 mm panels at 4.55 +/- 0.2 mm thickness with <2% porosity. Includes a Python/ABAQUS scripting workflow that generates a shell mesh and partitions it for the heterogeneous tow path - this is the closest thing to "open" geometry we have. LVI 15-55 J, HVI 250-400 ft/s, OHT, OHC.

10. **Kodagali, K., Vakili Rad, C., & Sockalingam, S. (2025).** *Low-velocity impact response of hybrid pseudo-woven carbon/epoxy thin composite laminates manufactured via automated fiber placement.* Composites Communications (early-2025 issue, Elsevier S2352492825030170).
    UNVERIFIED: full bibliographic record (volume/page/DOI suffix beyond the PII) was not directly confirmable through publisher pages from this run. Search snippets indicate it extends the 2024 work to thin-laminate (sub-2 mm) configurations.

11. **Chakraborty, D., Kodagali, K., Kattil, S.R., Miller, D., Sockalingam, S., & Sutton, M.A. (2024).** *On tack development during automated placement of uncured thermoset carbon/epoxy tows.* Manufacturing Letters, 40, 109-112.
    Adjacent process-mechanics paper from the Sockalingam group. Relevant if you want to model the cohesive layer between freshly placed tows in the as-manufactured (pre-cure) state, which determines whether undulations relax during cure.

### 4.3 Independent FEM / experimental work on AP-PLY (Edinburgh, Delft, Oxford)

12. **Li, X., Bai, R., Hou, S., Ma, P., Wang, J., & Lessard, L. (2021).** *Tow-wise modelling of non-conventional automated fibre placement composites: short beam shear study.* Composites Part A: Applied Science and Manufacturing, 151, 106660. https://doi.org/10.1016/j.compositesa.2021.106660
    First validation of the tow-wise modelling (TWM) framework on AP-PLY. Shows TWM intrinsically resolves gaps, overlaps, tow drops, thickness variations, and as-manufactured tow trajectories. Two AP-PLY configs with different crimp angles plus a UD baseline tested in short-beam-shear; DIC and CT validation.

13. **Li, X., et al. (2021).** *Automation of tow wise modelling for automated fibre placement and filament wound composites.* Composites Part A, 145, 106378. https://doi.org/10.1016/j.compositesa.2021.106378
    Companion / preceding paper detailing the TWM automation toolchain.

14. **Kok, R., Martinez-Hergueta, F., et al. (2022).** *Tensile response of AP-PLY composites: A multiscale experimental and numerical study.* Composites Part A: Applied Science and Manufacturing, 159, 107014. https://doi.org/10.1016/j.compositesa.2022.107014
    Multiscale FE framework that explicitly resolves through-thickness tow undulations. Validated on quasi-static tension. Compares AP-PLY against angle-ply baseline. Foundational reference for how to RVE-mesh a pseudo-woven geometry.

15. **Kok, R., Peroni, M., Martinez-Hergueta, F., & Pellegrino, A. (2023).** *Dynamic response of Advanced Placed Ply composites.* Composites Part B: Engineering, 248, 110347. https://doi.org/10.1016/j.compositesb.2022.110347 (`* primary text confirmed`)
    Split-Hopkinson tensile bar at 30 1/s. AP-PLY moduli are strain-rate independent; strengths slightly higher than baseline. Notes that hand layup of 10 mm slit tows (not AFP) was used for the Edinburgh specimens, with a "gap of three tow widths" between same-set tows. Cured 4 bar / 120 deg C / 120 min in a hot press. Specimens are 300 x 300 mm panels.

16. **Martinez-Hergueta, F., Kok, R., Cuvillo, R., Pernas-Sanchez, J., Rodriguez-Garcia, V., Artero Guerrero, J., & Guzman de Villoria, R. (2023).** *Impact Response of Automated Fibre Placement Advanced Placed Ply Composites.* Proc. 23rd International Conference on Composite Materials (ICCM23), Belfast, July 30 - Aug 4 2023, paper 339.
    Conference precursor to the 2024 *Composites Science and Technology* paper.

17. **Kok, R., Cuvillo, R., Rodriguez-Garcia, V., Pernas-Sanchez, J., Artero Guerrero, J., Guzman de Villoria, R., & Martinez-Hergueta, F. (2024).** *Low velocity impact response of Automated Fiber Placement Advanced Placed Ply composites.* Composites Science and Technology, 253, 110636. https://doi.org/10.1016/j.compscitech.2024.110636
    LVI + CAI on triaxial vs. quasi-isotropic vs. cross-ply AP-PLY. Cross-ply suffers a 49.1% strength loss at 50 J; triaxial and QI AP-PLY contain delamination much better. Identifies tow debonding, matrix cracking of impacted yarns, and bounded delamination as the failure mechanisms.

18. **Zheng, W., & Kassapoglou, C. (2024).** *Experimental tests and numerical simulation of delamination and fiber breakage in AP-PLY composite laminates.* Journal of Reinforced Plastics and Composites, 44(19-20), 1702-1714. https://doi.org/10.1177/07316844241245469 (`* primary text confirmed`)
    Three-point bending test of 24-layer AP-PLY laminate (AS4 carbon / epoxy). Tow width 6.35 mm, ply thickness 0.18 mm. Lay-up `[90_5 / 0_5 / (0_w / 90_f) / 90 / 0 / 0_5 / 90_5]`. Specimen 202 x 191.5 x 4.59 mm. Cohesive-zone delamination model in ABAQUS; Kc = 1.58 kJ/m^2. Predicted delamination initiation load is within 20% of test (~47-49 N/mm normalized to width). Cured per manufacturer spec in autoclave.

19. **Earlier delamination analysis paper (TU Delft):** *Delamination Analysis of A Class of AP-PLY Composite Laminates.* TU Delft repository file 9022491. UNVERIFIED authors and venue; the file did not parse cleanly via WebFetch and exact metadata could not be extracted.

---

## 5. Geometry / architecture details

### 5.1 UofSC LYNX AFP layups (Vakili Rad 2020 + Kodagali 2023)

| Parameter | Value | Source |
| --- | --- | --- |
| AFP machine | Gantry-based "LYNX" AFP, McNair Aerospace Center, UofSC | Kodagali 2023 thesis Fig 3.1; Vakili Rad 2020 thesis |
| Number of tow channels per course | 8 | Kodagali 2023 thesis ch.4 |
| Tow / slit-tape width (primary) | **6.35 mm (0.25 in / 1/4 in)** | Both theses; consistent across all UofSC PW/MAC papers |
| Tow / slit-tape width (alternate, parametric study) | 12.7 mm (0.5 in / 1/2 in) explored as a "doubled tow" check | Kodagali 2023 ch.3 Fig. 3.8 |
| Material - 2024 *Composites Part B* (Kodagali) | **Toray T800-SC-24K / P2362W-19** carbon/epoxy slit tape; FAW 191 gsm; 35.5% resin by weight | Kodagali 2023 thesis Table 4.2 |
| Material - 2020 *Composites Part B* (Vakili Rad) | **Hexcel IM7G/8552-1** carbon/epoxy slit tape, 1/4 in. | Vakili Rad 2020 thesis ch.6-7 |
| Material - earlier exploratory panels (Vakili Rad) | **AS4/8552** (hand-layup baseline) and **T-800SC-24K-10E** (early AFP trials) | Vakili Rad 2020 thesis |
| Cure cycle (autoclave) | **177 deg C / 350 deg F** at **0.62 MPa / 90 psi** for **6 h** (Vakili Rad early trials) or per manufacturer spec for IM7/8552 (later panels) | Vakili Rad 2020 thesis; Kodagali 2023 thesis |
| As-manufactured panel size | 1143 mm x 838.2 mm (45 in x 33 in) | Kodagali 2023 thesis ch.4 |
| Cured panel thickness (24-ply) | 4.55 +/- 0.2 mm | Kodagali 2023 thesis |
| Porosity (post-C-scan) | < 2% | Kodagali 2023 thesis |
| Placement-sequence notation | `[fiber angles][placement sequence (active/inactive channels)][angle shift a_s][tow width]` e.g. `[45,90,-45,0][10001000][0][6.35 mm]` | Kodagali 2023 thesis ch.3-4 |
| Gap between same-direction tows in MAC | 3 tow widths (i.e. one active channel followed by 3 inactive) for 4-ply MAC; 1 tow width (alternating 1010) for 2-ply MAC | Kodagali 2023 thesis ch.4 |
| Number of fiber angles m | 2 (e.g. [0/90], [45/-45]) or 4 ([45/90/-45/0]) | Kodagali 2023 thesis |
| Stacking sequences | 24-ply hybrid; P1-TWT (PW outer skins, UD inner core), P2-WTW (UD outer, PW inner), P3-Control [45/90/-45/0]_3s | Vakili Rad 2020 thesis; Kodagali 2024 *Composites Part B* |
| Symmetry strategy | Mirrored sub-laminates (e.g. [0,-45,90,45]) used to ensure overall layup symmetry | Kodagali 2023 thesis |
| Fiber volume fraction | All configurations contain identical material volume; nominal Vf is ply-level (UNVERIFIED specific number) | Kodagali 2023 thesis |
| RUC definition | Heuristic: "the unit where the stacking sequence (stiffness components) is approximately the same" - a Python/ABAQUS script generates the partitioned shell-element panel | Kodagali 2023 thesis ch.3 |

### 5.2 Edinburgh layups (Kok et al. 2023 dynamic paper)

| Parameter | Value |
| --- | --- |
| Tow width | 10 mm (slit by hand from a roll of prepreg) |
| Material | SHD Composites VTC401 (carbon/epoxy) |
| Manufacturing | Hand layup (NOT AFP), guided by 3D-printed alignment templates |
| Gap between same-set tows | 3 tow widths (= 30 mm) |
| Cure | Hot press, 4 bar, 120 deg C, 120 min |
| Panel size | 300 x 300 mm |
| Tow undulation length scale | ~1.5 mm (sub-DIC resolution of 3 mm) |

### 5.3 Delft (Zheng & Kassapoglou 2024 three-point-bending paper)

| Parameter | Value |
| --- | --- |
| Tow width | 6.35 mm |
| Ply thickness | 0.18 mm |
| Material | AS4 carbon / epoxy |
| Layup | `[90_5 / 0_5 / (0_w / 90_f) / 90 / 0 / 0_5 / 90_5]` 24-layer |
| Specimen | 202 x 191.5 x 4.59 mm (after machining) |
| Cure | Autoclave per material manufacturer spec |
| Cohesive Kc (mode-mixed, mode I dropped under shear-dominated 3-pt bend) | 1.58 kJ/m^2 |

---

## 6. Mechanical characterization data published

| Test | Standard / setup | Source(s) | Headline result |
| --- | --- | --- | --- |
| Warpage (asymmetric 4-ply) | Curvature measurement | Vakili Rad 2020 thesis; van Tooren & Sockalingam 2017 | 58% warpage reduction vs. asymmetric UD baseline |
| Uniaxial tension | ASTM D3039 (UNVERIFIED standard - thesis says "uniaxial tension") | Vakili Rad 2020; Vakili Rad et al. 2019 *Comp. Struct.* | Similar UTS to UD; higher strain-to-failure |
| Low-velocity impact | ASTM D7136, instrumented drop tower; 25-55 J | Vakili Rad 2020 (30, 50 J); Kodagali 2024 (15-55 J) | Up to 37% higher critical delamination load and 26% higher CAI strength after 55 J |
| Compression after impact (CAI) | ASTM D7137 (UNVERIFIED standard) | Kodagali 2024 *Comp. Part B* | 26% higher residual compressive strength vs. control after 55 J |
| High-velocity impact (ballistic) | ASTM D8101, single-stage gas gun, 250-400 ft/s, NASA Glenn | Vakili Rad et al. 2020 *Comp. Part B* | 45% reduction in back-face damage; 19.5% less back-face deflection (P1-TWT best) |
| Open-hole tension (OHT) | (Standard not pinned in extracted text - UNVERIFIED) | Kodagali 2023 thesis ch.5 | 7% increase in OHT strength; up to 16% reduction in strain near hole |
| Open-hole compression (OHC) | UNVERIFIED standard | Kodagali 2023 thesis ch.5 | Reported but specific numbers UNVERIFIED |
| Short-beam shear (ILSS) | EN 2563 / ASTM D2344 (UNVERIFIED specific standard) | Li et al. 2021 *Comp. Part A* (NOT UofSC, but on AP-PLY) | Used to validate tow-wise FE; full numbers in paper |
| High-strain-rate tension | Bespoke Split-Hopkinson bar at JRC Ispra; 30 1/s | Kok et al. 2023 *Comp. Part B* | Modulus rate-independent; strength slightly higher than baseline |

---

## 7. Modeling work

| Reference | Approach | Length scale | Code |
| --- | --- | --- | --- |
| Van Tooren & Sockalingam 2017; Vakili Rad et al. 2019 | 2D / 3D shell elements, sub-cell partitioning | Macro / laminate | ABAQUS (implied) |
| Kodagali 2023 thesis ch.3 | Python -> ABAQUS shell-element script that partitions a panel by tow path; composite-layup feature applies tow-by-tow stacking | Mesoscale shell | ABAQUS Python API |
| Li et al. 2021 (Edinburgh) | Tow-wise FE: each tow is meshed explicitly with as-manufactured trajectory; resin-rich pockets and undulations resolved | Mesoscale solid | UNVERIFIED (likely ABAQUS) |
| Kok et al. 2022, 2023 (Edinburgh) | Multiscale CDM with idealized "straight tow / resin-rich / undulation" unit-cell regions; isostrain to micro; CDM per region | Mesoscale to macro | ABAQUS C3D8R + custom CDM (UMAT/VUMAT - UNVERIFIED specific UMAT) |
| Zheng & Kassapoglou 2024 (TU Delft) | Cohesive-zone (mode II, shear-dominated 3pt-bend) for delamination + max-stress for fiber break | Macro/3D laminate | ABAQUS (cohesive elements) |
| Kok et al. 2024 (Edinburgh) | Multiscale + matrix cracking + delamination + tow debonding | Mesoscale | ABAQUS Explicit (UNVERIFIED) |

**No public RVE meshes have been located.** All modeling above appears to be performed in-house. The closest thing to a reusable artifact is the script described in Kodagali 2023 ch.3 / 4, which is described in narrative form but is not (as far as could be verified) released to a public repository.

---

## 8. Open-data inventory

| Asset | Location | Status |
| --- | --- | --- |
| Vakili Rad 2020 MSc thesis (PDF, 5.0 MB) | https://scholarcommons.sc.edu/etd/5708 (PDF: https://scholarcommons.sc.edu/cgi/viewcontent.cgi?article=6797&context=etd) | Open Access (downloaded, parsed) |
| Kodagali 2023 PhD dissertation (PDF, 10.2 MB) | https://scholarcommons.sc.edu/etd/7650 (PDF: https://scholarcommons.sc.edu/cgi/viewcontent.cgi?article=8555&context=etd) | Open Access (downloaded, parsed) |
| Nagelsmit 2013 PhD (TU Delft) | https://repository.tudelft.nl/islandora/object/uuid:f8c8ad4e-5d52-4bfd-8449-2748bad26673 | Open Access |
| Zheng & Kassapoglou 2024 paper PDF | https://pure.tudelft.nl/ws/portalfiles/portal/189816257/zheng-kassapoglou-2024-...pdf | Open Access (downloaded, parsed) |
| Kok et al. 2023 dynamic-response paper (Edinburgh ORA copy) | https://ora.ox.ac.uk/objects/uuid:6fc6b837-c2e7-46fd-9e72-596a6631d919 | Open Access |
| ICCM23 paper (Martinez-Hergueta et al. 2023) | https://www.research.ed.ac.uk/en/publications/impact-response-of-automated-fibre-placement-advanced-placed-ply- | Open metadata; PDF availability UNVERIFIED |
| GitHub / Zenodo / Mendeley Data with G-code, AFP `.cnc`, layup `.json`, or scan data | None located | No public artifact found in this run |
| Sockalingam group page | https://sc.edu/study/colleges_schools/engineering_and_computing/research/research_directory/sockalingam_group.php | Link; no datasets linked |
| neXt McNair / harik.org | https://nextusc.com (one URL redirected to an unrelated site at fetch time - may be stale); http://harik.org/tag/afp/ (connection refused at fetch time) | UNVERIFIED current state |

**Recommendation for FEM modeling.** Without an open machine-path file, the best path is to (a) use the placement-sequence notation from Kodagali 2023 ch.4 (explicitly given for `[45,90,-45,0][10001000][0][6.35 mm]`) plus the panel size 1143 x 838 mm, and (b) reproduce the Python/ABAQUS partition workflow described in ch.3-4 of the dissertation. Alternatively the Edinburgh hand-layup (10 mm tow, 30 mm gap, 300 x 300 mm) is fully specified in Kok et al. 2023 and is dimensionally simpler for an RVE.

---

## 9. Comparison to Delft AP-PLY

- **Lineage.** Z. Gurdal moved from TU Delft (where AP-PLY originated) to UofSC McNair, where he co-authors UofSC pseudo-woven papers (2020 *Comp. Part B*; 2024 *Comp. Part B*). This is the explicit human bridge between the two programs.
- **Terminology.** Delft and Edinburgh consistently use "AP-PLY"; UofSC uses "pseudo-woven" or "MAC". Kok et al. 2023 explicitly say AP-PLY "creates through-thickness reinforcements by interlacing fiber tows in a pseudo-woven architecture", confirming the equivalence.
- **Definition difference.** Vakili Rad 2020 thesis flags an explicit distinction: in true woven, every tow undulates; in **pseudo-woven** the first and last passes are still placed straight, leaving two non-undulating tows. This is the same as Delft's AP-PLY in practice but the UofSC convention names the asymmetry directly.
- **Manufacturing detail.** Delft (Nagelsmit) and UofSC (Vakili Rad / Kodagali) both use AFP. Edinburgh's experimental specimens used **hand layup of slit tows**, not AFP, despite being labelled "AP-PLY".
- **No paper directly publishes a head-to-head comparison** of the two architectures with the same materials and same test protocol; the literature treats them as the same architecture under two names.

---

## 10. Open questions / unknowns

1. **Exact slit-tape thickness (per ply).** Vakili Rad 2020 thesis quotes IM7/8552 with FAW; Kodagali 2023 quotes T800-SC FAW = 191 gsm and 35.5% resin. Neither directly reports the cured ply thickness in the extracted text, though the 24-ply / 4.55 mm panel implies ~0.190 mm/ply. (Zheng & Kassapoglou 2024 use 0.18 mm.) **Confirm from the source PDFs section-by-section before using in FEM.**
2. **Fiber volume fraction Vf.** Stated as "the same across all configurations" but a numerical value is not extracted (UNVERIFIED).
3. **Detailed cure cycle.** Vakili Rad 2020 mentions one early panel inadvertently held at dwell for a "prolonged time". Kodagali 2023 says "nominal" 177 deg C / 0.62 MPa / 6 h - not a full ramp/hold/ramp cycle.
4. **AFP machine specifics.** "LYNX AFP" in Kodagali 2023 - manufacturer (Electroimpact? Coriolis? Mikrosam?) is **not stated** in the extracted text. Worth a follow-up email to Sockalingam (sockalin@cec.sc.edu) or M. van Tooren.
5. **OHC, OHT exact specimen geometries and standards.** Reported but not pinned in the extracted text.
6. **Public datasets / G-code / scan data.** None located; recommend reaching out directly. McNair Center, Sockalingam group, and Harik group webpages do not currently link any public artifact.
7. **2025 thin-laminate paper bibliographic completeness.** PII S2352492825030170 in *Composites Communications* (early 2025) was indicated by search snippets; full citation not confirmed.
8. **Earlier "Delamination Analysis of A Class of AP-PLY Composite Laminates"** (TU Delft repository 9022491) - PDF was too large to re-parse cleanly; authors and venue could not be confirmed (UNVERIFIED).
9. **Whether the Sockalingam Python -> ABAQUS partitioning script has been or will be released.** Not located in any public repo as of this run.
10. **NIAR / Hexcel reference data** Vakili Rad 2020 cites "NIAR study on the properties of IM7/8552" for CAI - the specific NIAR report number is not pulled here.

---

## 11. Source URLs (consolidated)

- McNair Aerospace Center: https://sc.edu/about/centers_institutes/mcnair/
- Sockalingam (UofSC faculty): https://sc.edu/study/colleges_schools/engineering_and_computing/faculty-staff/subramani_sockalingam.php
- MSMOMS group: https://sc.edu/study/colleges_schools/engineering_and_computing/research/research_directory/sockalingam_group.php
- Vakili Rad 2020 thesis: https://scholarcommons.sc.edu/etd/5708
- Kodagali 2023 thesis: https://scholarcommons.sc.edu/etd/7650
- Nagelsmit 2013 thesis: https://repository.tudelft.nl/islandora/object/uuid:f8c8ad4e-5d52-4bfd-8449-2748bad26673
- Vakili Rad et al. 2020 *Composites Part B* (HVI): https://www.sciencedirect.com/science/article/abs/pii/S1359836820335265
- Kodagali et al. 2024 *Composites Part B* (LVI/CAI): https://www.sciencedirect.com/science/article/abs/pii/S1359836823006571
- Kok et al. 2023 *Composites Part B* (dynamic): https://www.sciencedirect.com/science/article/pii/S135983682200720X (and ORA mirror linked above)
- Kok et al. 2024 *Composites Science and Technology* (LVI): https://www.sciencedirect.com/science/article/pii/S0266353824002069
- Zheng & Kassapoglou 2024 JRPC: https://doi.org/10.1177/07316844241245469
- Li et al. 2021 *Composites Part A* (TWM short-beam-shear): https://doi.org/10.1016/j.compositesa.2021.106660
- Kok et al. 2022 *Composites Part A* (tensile multiscale): https://doi.org/10.1016/j.compositesa.2022.107014
- Van Tooren & Sockalingam 2017 ASC clutch laminate: https://dpi-proceedings.com/index.php/asc32/article/view/15199
- Vakili Rad et al. 2019 *Composite Structures* clutch: https://www.sciencedirect.com/science/article/abs/pii/S0263822318334354
- SAMPE 2019 paper Vakili Rad et al.: https://www.nasampe.org/store/viewproduct.aspx?id=13743102
