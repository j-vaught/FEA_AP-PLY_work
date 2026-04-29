# AP-PLY (Advanced Placed Ply) — Literature & Data Inventory

Curated for the FEA_AP-PLY modeling project. All entries below were located via direct web/database search; items I could not verify are flagged `[UNVERIFIED]`. ScienceDirect abstract pages frequently blocked direct fetches, so abstracts for those papers are summarized from search-result snippets and cross-referenced metadata (Google Scholar, repository preprints, ResearchGate listings) rather than full text.

---

## 1. Executive Summary

**What AP-PLY is.** AP-PLY (Advanced Placed Ply, sometimes "Advanced Placed Ply" or "AP-Ply") is an interlaced automated fibre placement (AFP) laminate architecture in which adjacent UD tows are placed with deliberate gaps; subsequent passes at different angles fill those gaps so that tows from different orientations physically interlace through the thickness, mimicking a woven fabric while still being made on a standard AFP head with no machine modifications. The interlacing creates through-thickness fibre connectivity that arrests delamination and improves CAI/impact response over baseline UD tape laminates.

**Origin and naming.** The architecture was conceived and named at TU Delft in collaboration with the Netherlands Aerospace Centre (NLR). The earliest located public artifact is a 2010 talk by C. Kassapoglou (TU Delft) at ISCM titled "AP-PLY: A new fibre placement architecture for fabric replacement", followed by the journal article M. Nagelsmit, C. Kassapoglou & Z. Gürdal, "AP-PLY: A New Fibre Placement Architecture for Fabric Replacement," *SAMPE Journal* 47(2), March/April 2011. The full architecture, manufacturing trials, and impact results are presented in M.H. Nagelsmit's PhD thesis, *Fibre Placement Architectures for Improved Damage Tolerance* (TU Delft, 2013, supervisor Z. Gürdal). NLR's public communications credit the joint TU Delft-NLR development.

**IP / patent.** No granted US/EP patent with the literal name "AP-PLY" was located through Google/Patents searches. NLR press releases describe AP-PLY as a jointly developed material concept but do not cite a patent number. Treat the IP situation as `[UNVERIFIED — no patent located]`.

**Principal investigators and lineages.** Three research lineages dominate the literature:

1. **TU Delft (origin).** Z. Gürdal (supervisor, then at TUD/USC), C. Kassapoglou (TUD), M.H. Nagelsmit (PhD 2013), W. Zheng (PhD 2016, supervisor Kassapoglou). NLR collaboration on manufacturing. Relevant chairs: Aerospace Structures & Computational Mechanics, Aerospace Manufacturing Technologies (Dransfeld is current chair, but no AP-PLY paper authored by him was located).
2. **University of Edinburgh (modeling + dynamic testing).** F. Martinez-Hergueta and F. Teixeira-Dias supervise; R. Kok (PhD c.2022) is the principal student, in collaboration with M. Peroni (JRC Ispra), A. Pellegrino (Oxford), and the Madrid IMDEA / UC3M group (R. Cuvillo, V. Rodriguez-Garcia, J. Pernas-Sánchez, J. Artero-Guerrero, R. Guzman de Villoria).
3. **University of South Carolina / McNair (independent pseudo-woven line).** M.J.L. van Tooren and S. Sockalingam coined "clutch laminate" / "pseudo-woven meso-architectured composite" (PW-MAC), with Z. Gürdal also affiliated. Students: C. Vakili Rad (MS 2020), K. Kodagali (PhD 2023). Their architecture is conceptually identical (tow-skip + interlace) but they generally avoid the TUD "AP-PLY" trademark-style name.
4. **UNSW / Deakin (Australia, tow-wise modeling).** X. Li, S.A. Brown, M. Joosten, G.M. Pearce; high-fidelity tow-wise FEM of AP-PLY/clutch laminates.

**Modeling state of the art.** Tow-level (a.k.a. tow-wise) FEM is the consensus mesoscale tool. Two open-source/openly-documented frameworks exist for AP-PLY:
- Rutger Kok's **Abaqus VUMAT + Python pre-processor** (`composite_cdm_ap_ply` and `ap_ply_model_creation` GitHub repos, LGPL-2.1, last updated 2022). Generates 3D solid-element AP-PLY geometry via Shapely, assigns tape/undulation/resin-rich material regions, and runs continuum damage mechanics with cohesive interactions tow-by-tow. Documented in Kok et al., *Composites Part A* 159 (2022) 106989.
- The Li / Brown / Joosten / Pearce **tow-wise modeling (TWM) algorithm** (UNSW/Deakin, *Composites Part A* 2021) — automated geometry from AFP tow paths, validated against AP-PLY short-beam-shear tests.

**What's missing.** No publicly archived raw test dataset, micro-CT scan, or G-code/AFP recipe was found on 4TU.ResearchData, Zenodo, or any institutional repository for AP-PLY. The two Kok GitHub repos contain code only — no mechanical raw data. The original Nagelsmit 2013 thesis PDF (72.8 MB at TU Delft repository) is open-access and is the most complete source for the original geometry definitions, but a full read was not possible within this session.

---

## 2. Annotated Bibliography

### 2.1 TU Delft / NLR — Origin & Foundational Work

**Nagelsmit, Kassapoglou & Gürdal (2010, ISCM talk).**
Kassapoglou, C. *AP-PLY: A new fibre placement architecture for fabric replacement*. Talk at International Symposium on Composites Manufacturing (ISCM), 1 May 2010. TU Delft Research Portal record (no DOI). This is the earliest public mention of the AP-PLY name located in a search; it predates the SAMPE Journal article and frames the concept as "fabric replacement" — i.e., a fibre placement layup that delivers fabric-like damage resistance without using fabric.

**Nagelsmit, Kassapoglou & Gürdal (2011, SAMPE Journal).**
Nagelsmit, M., Kassapoglou, C., Gürdal, Z. "AP-PLY: A New Fibre Placement Architecture for Fabric Replacement," *SAMPE Journal* 47(2), March/April 2011, pp. 36–45. ResearchGate publication 297311997. The seed publication. Describes the placement scheme: parallel UD tows are deposited with gaps equal to one or more tow widths; subsequent passes at +θ / -θ / 90 fill the gaps and cross over previously placed tows, producing a pseudo-woven structure. Reports first CAI improvements over UD baselines and notes compatibility with thermoset/thermoplastic prepreg slit tape and dry fibre with no AFP machine modification.

**NLR Executive Summary / report (2010, NLR-TP-XXX).**
A Dutch Aerospace Centre (NLR) executive-summary technical report by Nagelsmit & Kassapoglou attributed to ~2010 was located via core.ac.uk PDF (`files01.core.ac.uk/download/pdf/53034028.pdf`); the document repeats the SAMPE Journal results. The exact NLR-TP report number was not extractable in this session — `[UNVERIFIED report number]`. The same press communication says the joint TU Delft-NLR work shows ~10% structural-weight saving potential.

**Nagelsmit, M.H. (2013, PhD thesis).**
Nagelsmit, M.H. *Fibre Placement Architectures for Improved Damage Tolerance.* PhD dissertation, Delft University of Technology, 2013. Supervisor: Z. Gürdal. DOI: 10.4233/uuid:f8c8ad4e-5d52-4bfd-8449-2748bad26673. PDF "Thesis_Nagelsmit.pdf" (72.8 MB) is open access at the TU Delft repository. Defines AP-PLY architecture, manufacturing trials, and CAI test program; reports up to 15% increase in compression-after-impact strength for one configuration plus smaller delamination footprints and lower BVID thresholds. This is the canonical source for the original parameters.

**Zheng, W. (2016, PhD thesis).**
Zheng, W. *Delamination Analysis of A Class of AP-PLY Composite Laminates.* PhD dissertation, Delft University of Technology, 2016. DOI: 10.4233/uuid:26d1ec89-e6e6-4d5f-b4d0-2591c4aa406a. Supervisor C. Kassapoglou. Develops 2D plane-stress models with cohesive elements at woven interfaces and an analytical SERR / beam-theory method including a stiffness-discontinuity correction factor. Investigates effect of crimp angle, number of woven plies, number of straight filling plies, and through-thickness woven location.

**Zheng & Kassapoglou (2017, Composites Part A).**
Zheng, W., Kassapoglou, C. "Prediction of delamination onset and growth for AP-PLY composite laminates using the finite element method," *Composites Part A* 101 (2017) 381–393. DOI: 10.1016/j.compositesa.2017.06.032. Couples interlaminar-stress recovery with a max-stress criterion for delamination onset and 2D cohesive-element FEM for propagation; parametric study of woven angle, number of woven plies, and number of straight filler plies on delamination response.

**Zheng & Kassapoglou (2019, Journal of Composite Materials).**
Zheng, W., Kassapoglou, C. "Energy method for the calculation of the energy release rate of delamination in composite beams," *Journal of Composite Materials*, 2019. DOI: 10.1177/0021998318785952. Energy-based beam-theory SERR for cracks in laminates with variable through-thickness stiffness; correction factor for small-crack accuracy; bimaterial-tip singularity. Forms the analytical basis for AP-PLY tailoring.

**Zheng & Kassapoglou (2024, Journal of Reinforced Plastics and Composites).**
Zheng, W., Kassapoglou, C. "Experimental tests and numerical simulation of delamination and fiber breakage in AP-PLY composite laminates," *J. Reinforced Plastics and Composites* 44(19-20) (2024) 1702–1714. DOI: 10.1177/07316844241245469. Three-point-bending tests on AP-PLY beams plus cohesive-zone FEM; predicted damage-initiation load within ~20% of test. Open access PDF on TU Delft pure repository (path `pure.tudelft.nl/.../zheng-kassapoglou-2024-...`).

### 2.2 University of Edinburgh — Multiscale Modeling & Dynamic Response

**Kok, Martinez-Hergueta & Teixeira-Dias (2022, Composites Part A).**
Kok, R., Martinez-Hergueta, F., Teixeira-Dias, F. "Tensile response of AP-PLY composites: A multiscale experimental and numerical study," *Composites Part A* 159 (2022) 106989. DOI: 10.1016/j.compositesa.2022.106989. Open access (CC-BY). Funded by Royal Society UK grant RGS/R2/180091. The reference paper for the open-source modeling pipeline: cross-ply and quasi-isotropic AP-PLY tensile coupons, baseline angle-ply controls, 3D multiscale FEM with explicit through-thickness undulations, Tan-style failure criteria, Maimi CDM. Confirms stiffness ≈ baseline, strength is layup-dependent.

**Kok, Peroni, Martinez-Hergueta & Pellegrino (2023, Composites Part B).**
Kok, R., Peroni, M., Martinez-Hergueta, F., Pellegrino, A. "Dynamic response of Advanced Placed Ply composites," *Composites Part B: Engineering* 248 (2023) 110347. DOI: 10.1016/j.compositesb.2022.110347. Bespoke split-Hopkinson tension bar at JRC Ispra, sized for AP-PLY's large RVE. Tensile loading at ~30 s⁻¹ vs quasi-static. Moduli rate-independent; strengths marginally higher than baseline. Improved damage tolerance under dynamic loading from suppression of out-of-plane tow straightening.

**Peroni, Kok, Martinez-Hergueta & Pellegrino (2022, conference, LWAG Freiburg).**
Peroni, M., Kok, R., Martinez-Hergueta, F., Pellegrino, A. "Hopkinson tensile testing on large specimens: application to Advanced Placed Ply AFP composites," Light-Weight Armour for Defense & Security (LWAG) conference, Freiburg, 2022. Edinburgh Research Explorer record. Companion paper describing the experimental apparatus designed around AP-PLY's >10 mm characteristic length scale. `[UNVERIFIED proceedings DOI]`

**Martinez-Hergueta, Kok, Cuvillo, Rodriguez-Garcia, Pernas-Sánchez, Artero-Guerrero, Guzman de Villoria (2023, ICCM23).**
"Impact response of automated fibre placement Advanced Placed Ply composites," 23rd International Conference on Composite Materials, Belfast, 2023. Edinburgh Research Explorer + ResearchGate publication 373218580. 50 J low-velocity impact + CAI on triaxial, quasi-isotropic, and cross-ply AP-PLY panels. Cross-ply showed 49.1% residual-strength loss at 50 J; triaxial and QI AP-PLY confined the delamination footprint.

**Kok, Cuvillo, Rodriguez-Garcia, Pernas-Sánchez, Artero-Guerrero, Guzman de Villoria, Martinez-Hergueta (2024, Composites Science and Technology).**
"Low velocity impact response of Automated Fiber Placement Advanced Placed Ply composites," *Composites Science and Technology* 253 (2024) 110636. DOI: 10.1016/j.compscitech.2024.110636. Journal extension of the ICCM23 paper; matrix cracking, delamination, and tow debonding identified as the controlling failure mechanisms; numerical model captures all three.

**Kok, R. (c. 2022, PhD thesis, University of Edinburgh).**
Kok, R. *Advanced Placed Ply laminate architectures for improved impact tolerance: experimental characterization and simulation.* PhD thesis, University of Edinburgh, c. 2022/2023. Supervisors: F. Martinez-Hergueta and F. Teixeira-Dias. The exact ERA (era.ed.ac.uk) handle could not be located through search — `[UNVERIFIED ERA handle]`. The thesis underpins the four Kok-coauthored journal/conference papers above.

### 2.3 University of South Carolina / McNair — "Clutch" / "Pseudo-Woven" Lineage

The USC group describes architecturally identical laminates under different names ("clutch laminate", "pseudo-woven", "meso-architectured composite", PW-MAC); all use AFP tow-skip + interlace. They cite Nagelsmit and credit the architecture, but generally do not use the AP-PLY trademark.

**Van Tooren & Sockalingam (2017, ASC32).**
Van Tooren, M.J.L., Sockalingam, S. "Non-conventional Composite Laminates – Clutch Laminates," Proceedings of the American Society for Composites 32nd Technical Conference, 2017. DOI: 10.12783/asc2017/15199. Coins "clutch laminate"; demonstrates dramatic warpage reduction (4-ply asymmetric layups can be self-balanced via interlacing).

**Vakili Rad, Thomas, Seay, van Tooren & Sockalingam (2019, Composite Structures).**
"Manufacturing and characterization of novel clutch non-conventional fiber-reinforced composite laminates," *Composite Structures* 210 (2019) 391–401. DOI: 10.1016/j.compstruct.2018.11.057 (`S0263822318334354`). Carbon/epoxy slit tape via gantry AFP at McNair Aerospace Center, USC; sub-cell shell-FE modeling.

**Vakili Rad, C. (2020, MS thesis, USC).**
Vakili Rad, C. *The Manufacturing and Characterization of Pseudo-Woven Carbon Fiber Composite Architectures for Enhanced Damage Tolerance.* MS thesis, University of South Carolina, Spring 2020. Supervisor: S. Sockalingam. Open access at scholarcommons.sc.edu/etd/5708. Reports up to 58% warpage reduction in 4-ply pseudo-woven vs asymmetric baseline.

**Bossuyt et al. (2020, Composites Part B).**
"High velocity impact response of hybridized pseudo-woven carbon fiber composite architectures," *Composites Part B: Engineering* 207 (2021), DOI prefix `S1359836820335265`. IM7/8552 slit tape, 24-ply hybrid (PW-UD-PW), gas-gun ASTM D8101 impact at 250–400 ft/s. `[Author list unverified — searches consistently returned the title and journal but the author block remained unresolved]`

**Kodagali, Vakili Rad, Sockalingam, Gürdal & Miller (2024, Composites Part B).**
"Low velocity impact and compression-after-impact response of hybrid pseudo-woven meso-architectured carbon/epoxy composite laminates manufactured via automated fiber placement," *Composites Part B* 271 (2024) 111154. DOI: 10.1016/j.compositesb.2023.111154 (`S1359836823006571`). Toray T800-SC-24K / 3900 carbon/epoxy slit tape, gantry AFP at McNair. Hybrid PW-UD-PW config achieves +37% critical delamination load and +26% residual compressive strength after 55 J impact.

**Kodagali, K. (2023, PhD thesis, USC).**
Kodagali, K. *Studies of Damage Tolerance in Automated Fiber Placement Based Heterogeneous Meso-Architectured Carbon/Epoxy Composite Laminates.* PhD dissertation, University of South Carolina, Fall 2023. Supervisor: S. Sockalingam. Open access at scholarcommons.sc.edu/etd/7650.

**Kodagali, Vakili Rad, Sockalingam (2025, Materials Today Communications).**
"Low-velocity impact response of hybrid pseudo-woven carbon/epoxy thin composite laminates manufactured via automated fiber placement," DOI prefix `S235249282503017X`. Thin-laminate variant.

### 2.4 UNSW / Deakin — Tow-Wise Modeling

**Li, Brown, Joosten & Pearce (2021, Composites Part A — automation).**
Li, X., Brown, S.A., Joosten, M., Pearce, G.M. "Automation of tow wise modelling for automated fibre placement and filament wound composites," *Composites Part A* 146 (2021) 106382 (`S1359835X21001718`). DOI: 10.1016/j.compositesa.2021.106382 `[UNVERIFIED exact volume — DOI prefix confirmed]`. Algorithm to read AFP tow paths and build tow-wise FE meshes capturing gaps, overlaps, drops, waviness; AP-PLY as a worked example.

**Li, Brown, Joosten & Pearce (2021, Composites Part A — short beam shear).**
Li, X., Brown, S.A., Joosten, M., Pearce, G.M. "Tow wise modelling of non-conventional automated fibre placement composites: short beam shear study," *Composites Part A* (Dec 2021) (`S1359835X21004784`). Validates the TWM mesh against AP-PLY ASTM D2344 short-beam-shear specimens; correctly predicts crack deflection and mixed-mode delamination.

**Joosten, Li, Huang (2021, ASC36).**
"3D Printed Continuous Fibre Composite Research at Deakin University: Design and Analysis Methods for UD, Hybrid and Pseudo-Woven Ply Architectures," ASC 36th, 2021. DOI: 10.12783/asc36/35755. Extends pseudo-woven ply concepts to continuous-fibre 3D printing.

**Li, X. (2022, dissertation).**
Li, X. *Mesoscale Numerical Modelling and Failure Prediction of Automated Fibre Placement Composites.* Dissertation, c. 2022 — institution and exact date `[UNVERIFIED]`; cited via Google Scholar.

### 2.5 Open-Source Software (the two repositories most useful for this project)

**Kok, R. — `composite_cdm_ap_ply` (GitHub, 2022).**
URL: https://github.com/rutger-kok/composite_cdm_ap_ply. License: LGPL-2.1. Language: Fortran 77 fixed-format. 3D continuum damage mechanics VUMAT for Abaqus/Explicit, three source files: `composite_cdm_ap_ply.for` (main), `resin_damage.for`, `rotation_matrix.for`, plus `state_vars.xlsx`. Tan et al. (2015) failure criteria + Maimi CDM evolution + Shah et al. (2021) 3D woven extensions. Required material card includes the AP-PLY-specific parameters `und_angle1`, `und_angle2` (in-plane undulation angles for the two interlaced constituents), `und_length`, `cpt` (cured ply thickness).

**Kok, R. — `ap_ply_model_creation` (GitHub, 2022).**
URL: https://github.com/rutger-kok/ap_ply_model_creation. License: LGPL-2.1. Language: Python 2.7 (Abaqus/CAE Python API + Shapely). Files: `ap_ply_model.py` (top-level driver), `tape_placement.py` (geometry generation), `ap_ply_3d.py`, `ap_ply_elastic.py`, `ap_ply_materials.py`, `sigc.py`, `abaqus_model.py`, plus `impact/` and `tensile/` subdirectories. `tape_placement.laminate_creation()` is the geometric core: given `tape_angles`, `tape_widths`, `tape_spacing` (number of tape-width gaps between tows in a single pass), `cured_ply_thickness`, `undulation_ratio` (height/length of an undulation), and the specimen polygon, it produces tow-wise polygons that mark tape regions, undulation regions, and resin pockets, then assigns four material types (tape elastic, tape damage, undulation damage, resin-rich damage). The default material is **SHD Composites VTC401** UD carbon/epoxy with as-coded properties E11=124.35 GPa, E22=E33=7.231 GPa, ν12=ν13=0.339, ν23=0.374, G12=G13=3.268 GPa, G23=2.632 GPa, Xt=2.550 GPa, Xc=1.102 GPa, Yt=0.131 GPa (in-situ), Yc=0.184 GPa, SL=0.122 GPa, ST=0.08273 GPa, G1+=0.133 N/mm, G1-=0.095 N/mm, G2+=0.00038 N/mm, G6=0.00162 N/mm, ρ=1.59e-6 kg/mm³. A second material `HiTapeUD210` (Hexcel HiTape UD210) is also coded.

### 2.6 Other relevant work

**Bergsma, Koussios, Bersee, Lindstedt & Parlevliet (2006, ASC).**
"Trinity essence for composite design," Proceedings of the American Society for Composites 21st Conference, Lancaster PA, 17-20 Sep 2006, ISBN 1-932078-60-6. Predates AP-PLY but is the Bergsma group's design-philosophy paper that motivated subsequent fibre-architecture work at TU Delft.

**Lopes-style cohesive-element matrix-cracking models** (cited in Zheng's thesis but not specific to AP-PLY) and the Tan-Falzon-Chiu-Price (2015) and Shah-Megat-Yusoff et al. (2021) failure-criteria papers underpin the Edinburgh CDM. Maimi-Camanho-Mayugo-Davila (2007 Parts I & II, *Mech. of Materials* 39, DOIs `10.1016/j.mechmat.2007.03.005` and `.006`) provide the damage-evolution law.

---

## 3. Geometry / Architecture Table

All numerical values below come from the indicated source. Default values in the Kok pre-processor are the most explicit and usable for FEM input.

| Parameter | Symbol | Value(s) | Source |
|---|---|---|---|
| Tape (tow) width | w | 6.35 mm (1/4 in slit tape, common AFP) | Multiple papers; Edinburgh Composites Part A 2024 quotes 6.35 mm |
| Tape (tow) width | w | 10.0 mm | Kok `ap_ply_model.py` `__main__` example |
| Tape (tow) width | w | 12.7 mm (1/2 in slit tape) | Kok `tape_placement.py` `__main__` example; Nagelsmit thesis (1/8" slit tape, 4-up = 12.7 mm bandwidth) |
| Cured ply thickness | cpt | 0.18 mm | Kok `ap_ply_model.py` example (typical for 200 gsm UD) |
| Cured ply thickness | cpt | 0.18125 mm | Kok comment in source |
| Cured ply thickness | cpt | 0.213 mm | Kok `tape_placement.py` example |
| Cured ply thickness | cpt | ≈ 0.18 mm | Edinburgh papers (consistent with VTC401 UD) |
| Undulation ratio | u_ratio = h/L | 0.09 | Kok `ap_ply_model.py` example |
| Undulation ratio | u_ratio | 0.104 | Kok `tape_placement.py` example |
| Undulation half-width | u_width = (cpt / u_ratio) / 2 | 1.0 mm (default) ≈ cpt/(2·u_ratio) | Kok pre-processor formula |
| Tape spacing (gaps between tows in one pass) | s | 1, 2, 3 (integer multiples of w) | Kok pre-processor `tape_spacing` argument |
| Tape angles (tested) | θ | (0, 90), (0, 45, 90, -45), triaxial (0, +60, -60) | Kok pre-processor + Edinburgh Comp Sci Tech 2024 |
| Number of plies (RVE) | l_plies | 4 (default test), 8, 16, 24 | Kok + USC papers |
| Specimen RVE (Kok tensile) | — | typically 40 × 40 mm region modeled | Kok 2022 supplementary code |
| BVID-equivalent impact energy | — | 50 J (Edinburgh CAI), 55 J (USC PW-MAC) | Edinburgh/USC papers |
| AFP material (default) | — | SHD VTC401 UD carbon/epoxy | Kok Composites Part A 2022 |
| AFP material (alt) | — | Hexcel HiTape UD210 | Kok pre-processor |
| AFP material (USC) | — | Toray T800-SC-24K / 3900 prepreg slit tape | Kodagali 2024 |
| AFP material (USC alt) | — | IM7/8552 slit tape | Bossuyt 2020 / Kodagali earlier work |
| Interlacing pattern | — | One-over-one (`s=1`), two-over-two (`s=2`), three-over-three (`s=3`) — in Kok terminology this is the number of tape-width gaps left between adjacent tows in a single pass; subsequent passes (offset by `n_shift × w`) fill those gaps | Kok `tape_placement.laminate_creation` source |
| AFP machine recipe | — | No open G-code or open recipe found | search returned no public AFP program file |

`[NOT FOUND]` Open-data items: explicit lay-up sequences for each Edinburgh test, panel-level coupon dimensions, raw load-displacement curves, and micro-CT scan files were not located in any public dataset.

---

## 4. Open-Data Inventory

| Asset | Where | Status |
|---|---|---|
| Nagelsmit 2013 PhD thesis (full geometry, manufacturing, CAI data) | TU Delft repository, DOI 10.4233/uuid:f8c8ad4e-5d52-4bfd-8449-2748bad26673 | Open access PDF (72.8 MB) |
| Zheng 2016 PhD thesis | TU Delft repository, DOI 10.4233/uuid:26d1ec89-e6e6-4d5f-b4d0-2591c4aa406a | Open access |
| Zheng & Kassapoglou 2024 paper | TU Delft pure (`pure.tudelft.nl/.../zheng-kassapoglou-2024-...`) | Open access PDF |
| Kok et al. 2022 paper (Composites Part A) | Elsevier ScienceDirect | Open access (CC-BY) |
| Vakili Rad MS thesis 2020 | scholarcommons.sc.edu/etd/5708 | Open access |
| Kodagali PhD thesis 2023 | scholarcommons.sc.edu/etd/7650 | Open access |
| Kok GitHub `composite_cdm_ap_ply` (Abaqus VUMAT) | github.com/rutger-kok/composite_cdm_ap_ply | LGPL-2.1, code only |
| Kok GitHub `ap_ply_model_creation` (Python pre-processor) | github.com/rutger-kok/ap_ply_model_creation | LGPL-2.1, code only |
| 4TU.ResearchData AP-PLY dataset | searched with `AP-PLY`, `Bergsma`, `Nagelsmit` | **None found** |
| Zenodo AP-PLY dataset | searched | **None found** |
| Micro-CT scans (any of the above papers) | searched | **None located in public archives** |
| AFP G-code / machine-recipe files | searched | **None found in public sources** |

---

## 5. Open Questions / Things To Verify

1. **Patent.** Confirm whether TU Delft / NLR ever filed an EP, NL, or US patent on AP-PLY. Search Espacenet by inventor (Nagelsmit, Kassapoglou, Gürdal) and applicants (TUD, NLR/Stichting NLR). Currently `[UNVERIFIED]`.
2. **Original SAMPE Journal page numbers** for Nagelsmit/Kassapoglou/Gürdal 2011. The article exists in nxtbook v.47 issue 2; pp. 36–45 is my best estimate from the index entry but should be confirmed by reading the issue.
3. **NLR-TP report number** for the Nagelsmit/Kassapoglou 2010 NLR document on AP-PLY. The full identifier (typically NLR-TP-YYYY-NNN) was not extracted.
4. **Kok PhD thesis** ERA (era.ed.ac.uk) handle and exact deposit year. The thesis is referenced in multiple places but the canonical archive record was not located.
5. **Bossuyt et al. 2020 author block** for the *Composites Part B* paper "High velocity impact response of hybridized pseudo-woven carbon fiber composite architectures." Search results consistently returned the title but not the author list.
6. **Tow-wise modeling 2021 papers — verify authors.** Independent verification of Li, Brown, Joosten, Pearce as authors of both `S1359835X21001718` and `S1359835X21004784` is recommended (one search returned these names, others did not).
7. **A "Tailoring of AP-PLY composite laminates for improved performance in the presence of delaminations" (Zheng, Kassapoglou & Zheng, *Composite Structures* 2019)** appears in one Google Scholar listing (cited 21 times) but was not corroborated by direct fetch of an Elsevier or repository record. If real, this would be the third Zheng-Kassapoglou journal paper (in addition to 2017 and 2024). `[UNVERIFIED — flagged for follow-up]`
8. **Open mechanical raw data.** No raw stress-strain, CT, or DIC data appears to be archived. For an FEM validation task, contacting the corresponding authors (Kok, Martinez-Hergueta, Sockalingam, Kassapoglou) is the recommended next step.
9. **Drechsler / TUM connection.** No AP-PLY paper authored by Drechsler was located. Drechsler's AFP defect-modeling work (e.g., Böckl, Wedel, Misik, Drechsler 2023) is adjacent but does not study the AP-PLY architecture per se.

---

## 6. Recommended starting points for FEM modeling

For a tow-level finite element model of AP-PLY in an open-source FEA stack, the best literature-to-code pipeline is:

1. Download the Nagelsmit 2013 thesis to lock the canonical interlacing geometry definition.
2. Mirror the Kok `ap_ply_model_creation` Python (Shapely) algorithm for generating tow polygons. Replace the Abaqus calls with a meshing path that targets your open-source solver (e.g., Gmsh + CalculiX, Gmsh + code_aster, or FEniCS via a custom mesher).
3. Use the Kok 2022 *Composites Part A* paper as the validation benchmark (open-access, has tensile + multiscale FEM results) and the Li/Pearce 2021 short-beam-shear paper as the ILSS validation case.
4. For impact / CAI, use the Edinburgh 2024 *Composites Sci Tech* paper and the USC Kodagali 2024 *Composites Part B* paper as cross-validation targets at two energy levels (50 J and 55 J).
5. Borrow the Tan/Maimi/Shah failure-law constants directly from the Kok VUMAT until you replace them with your own material card.
