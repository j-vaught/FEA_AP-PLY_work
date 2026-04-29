# Phase 2 consolidation review

All 16 per-test specs are written. Total ≈ 7,900 lines of `spec.md` + 13,000 lines of skeleton `runner.py`. Each spec follows the 10-section format with primary citations and a quantitative pass criterion. This document is the cross-stage consistency pass — not new work, but reconciling choices that drifted across the 16 parallel agents.

## What every spec passed

- Solid HEXA8 (with TETRA10 fallback explicitly documented where mesh topology forced it — stage 10 only).
- SI units throughout.
- 10-section format honored.
- Primary references cited; no fabricated entries.
- Runner skeleton with Lima / Apptainer wrapping for macOS.
- Brand-color Typst/CeTZ figures, no matplotlib.

## Items requiring a consolidation decision

### A. Material card source for IM7/8552

Eight stages depend on a single canonical IM7/8552 elastic + strength card. Three different cited sources are in play.

| Stage | Source cited | Numerical drift |
|---|---|---|
| 6 | WWFE-II Part B, Kaddour-Hinton 2013 | E1=165, E2=9.0, G12=5.6 |
| 7, 8, 13, 14, 15, 16 | Soden-Hinton-Kaddour 1998 | E1=171.4, E2=9.08, G12=5.29 |
| 9 | Camanho-Maimi-Davila 2007 / CMH-17 V3 | E1=171.4, E2=9.08, G12=5.29 (same numbers, different cite) |
| 10 | Soden 1998 (fiber + matrix isotropic only) | fiber E=230, matrix E=4.08 |
| 11 | Soden 1998 + Kok 2022 architectural correction | E1=161, E2=11.4 (Kok's VTC401-aligned values) |

**Proposed canonical card** (single source of truth, save once at `references/material_cards/im7_8552.json`):
- Elastic + Poisson + densities → **Soden-Hinton-Kaddour 1998** (the WWFE benchmark card used by 7 of the 8 composite-bearing stages).
- Strengths → **Kaddour-Hinton 2013 (WWFE-II Part B)** (richer than Soden alone for Hashin/Puck/Tsai-Wu side-by-side at stage 6).
- Cohesive (G_Ic, G_IIc, T_n, T_s, η_BK) → **Camanho-Davila 2003 (NASA TM-2002-211737)**.

**Stages to revise**: 9 (re-cite Soden 1998 even though numerical values are unchanged); 11 (clarify that the Kok 2022 column is the *cross-check target*, not the *input card*).

### B. Per-element fiber-orientation mechanism

OpenRadioss offers two equivalent ways to set the orthotropic frame on `/PROP/TYPE14`:
- **Mechanism A**: `Phi` angle in the property card, with `Iorth=1` and a reference vector.
- **Mechanism B**: `/SKEW/FIX` referenced via `Iskew` slot.

| Stage | Mechanism |
|---|---|
| 7 | A |
| 8 | B |
| 9 | A |
| 11 | B |
| 15 | B |

**Proposed canonical mechanism**: **B (`/SKEW/FIX`)**. Used in 3 of 5 stages including the more complex (mesoscale, ballistic) ones. Stage 8's analysis chose B explicitly for "robust against future Radioss-card layout drift." Stage 7's analysis noted A was used for simplicity but B is documented as the alternative.

**Stages to revise**: 7 and 9 — switch to /SKEW/FIX, keep A documented as alternative.

### C. Cohesive zone card for delamination

| Stage | Card chosen |
|---|---|
| 12 | /INTER/TYPE2 primary + /MAT/LAW117 fallback |
| 13 | /INTER/TYPE2 |
| 14 | inherits stage 13 |
| 15 | /INTER/TYPE2 between every adjacent ply |
| 16 | /INTER/TYPE2 |

**Stage 12's audit finding** (the most important one across all 16 specs): `/INTER/TYPE2` was historically a **kinematic-tie-with-brittle-failure** card, NOT a traction-separation cohesive law. The cohesive-capable behavior requires `Spotflag=25` in modern OpenRadioss versions. If we run on an older starter, every inter-ply interface in stages 13-16 silently becomes a brittle tie — which would catastrophically over-predict V50.

**Proposed standard**: at toolchain bring-up time, runner.py must inspect the OpenRadioss starter version and, if pre-Spotflag-25, switch to **/MAT/LAW117 zero-thickness solid cohesive elements** automatically. Stage 12's runner already does this; stages 13, 14, 15, 16 do not.

**Stages to revise**: 13, 14, 15, 16 — add the same auto-fallback that stage 12 has, plus a startup-time version check that fails loudly if the starter is too old AND LAW117 isn't available either.

### D. Test-progression literature stage 13/14 swap

`references/test_progression_literature.md` has stages 13 and 14 in the opposite order from `master_plan.md` (literature §13=CAI, §14=LVI vs. master plan §13=LVI, §14=CAI). The stage 13 agent caught this and followed the master plan. **Action**: update the literature doc so the row order matches master plan.

### E. Kok preprocessor → Abaqus dependency

Stages 11 and 16 invoke `rutger-kok/ap_ply_model_creation` through `abaqus cae noGUI`. The Kok repo is LGPL-2.1 source code, but it is an Abaqus CAE Python plugin and *requires* Abaqus to run. This violates the open-source-only constraint at the geometry-generation step.

**Three resolution paths**:
1. **Port** Kok's tow-laydown logic to standalone GMSH + `gmsh.model.occ` Python. Estimated effort: 2-4 weeks of focused work; the algorithm is documented in Nagelsmit 2013 PhD §3.4 and Kok 2022 Composites Part A §2. Output is a meshable solid model with per-tow material orientation.
2. **Alternate open generator**: the wider ecosystem survey identified TexGen (GPL-2, Nottingham) for woven textile geometry. TexGen does not natively understand AP-PLY but its meta-modeller could be extended. Effort comparable to (1) and probably more total work.
3. **Accept Abaqus-only for geometry generation**, then run OpenRadioss on the resulting `.inp`. Cleanest engineering path but breaks the "no closed source anywhere" constraint.

**Recommendation**: path (1). The user has already committed to single-tool OpenRadioss for the *solver*; the geometry generator is a separable concern, and a clean Python port has compounding value for stage 16 plus any future architectural variants. **This is the largest hidden cost in the project.**

### F. Hourglass control conventions

Different stages use different hourglass-control choices because the strain-rate regimes are different.

| Stage | Hourglass choice | Rationale |
|---|---|---|
| 1, 3, 4 | none / default | implicit static, no concern |
| 2 | Ismstr=11, Icpre=1 | small-strain co-rotational + constant-pressure for HEXA8 bending lock |
| 9 | default | implicit static, multi-ply solid stack |
| 13 | implicit-then-explicit, default in implicit | LVI is short-duration but the impactor is hemispherical (smooth) |
| 15, 16 | Isolid=24 + Ihq=8 (HEPH Belytschko-Bindeman 1993) | ballistic — viscous Ihq=1 explicitly rejected because it absorbs ballistic KE |

These differences are **intentional and correct** — keep as documented.

### G. Sources for V50 ground truth

Stage 15 and stage 16 both flag the V50 reference value as `[UNVERIFIED]` because:
- Stage 15: a published IM7/8552 4 mm ballistic V50 in the open literature is needed.
- Stage 16: Vakili Rad 2020 V50 number is paywalled; Kodagali 2023 PhD has it but needs page-level digitization.

**Action item**: a small task to digitize the V50 numbers from the UofSC theses (Vakili Rad 2020 MSc, Kodagali 2023 PhD) — both are open-access on Scholar Commons. This is a 30-minute literature task before stages 15/16 can return a PASS/FAIL verdict instead of INCONCLUSIVE.

### H. Computational budget reality check (stage 16)

Stage 16's spec predicts **8-24 hours per shot on 16 cores, 3-7 days per V50 sweep** for the reduced 200×200 mm representative section (~692k elements). Full 1143×838 mm panel is ~17M elements and not feasible without HPC.

**Action item**: confirm whether HPC access (UofSC HPC, NSF ACCESS, etc.) is available before stage 16 is attempted. Stages 1-15 are all laptop-feasible; stage 16 is the only HPC-class task in the ladder.

---

## Punch list (in priority order)

1. **Standardize the IM7/8552 material card** (item A) — write `references/material_cards/im7_8552.json` and have all 8 dependent stages reference it.
2. **Standardize cohesive card with auto-fallback** (item C) — push stage 12's runner pattern to stages 13, 14, 15, 16.
3. **Fix test_progression_literature.md stage 13/14 swap** (item D) — small file edit.
4. **Standardize fiber-orientation mechanism on /SKEW/FIX** (item B) — update stages 7 and 9.
5. **Decide Kok geometry path** (item E) — this is the largest decision left and unblocks stages 11 and 16.
6. **Digitize V50 ground truth** (item G) — read the UofSC theses for the actual numbers.
7. **Confirm HPC availability** (item H) — gates stage 16 timing.

---

## Phase 3 entry point

Stage 1 (linear-elastic beam) is the toolchain smoke test. Independent of items A-H above (uses /MAT/LAW1 isotropic, no composite, no cohesive, no impact). Recommended starting point regardless of how items A-H resolve. Spec at `tests/stage_01_beam_bending/spec.md`, runner at `tests/stage_01_beam_bending/runner.py`, both ready to flesh out (Jinja deck template + GMSH mesh builder are the open `NotImplementedError` TODOs).

The Phase 3 first-runnable plan is unchanged: **Stage 1 implementation** — write `deck.rad.j2`, write the GMSH builder, install Lima + OpenRadioss, run end-to-end on the toolchain, capture closed-form pass.
