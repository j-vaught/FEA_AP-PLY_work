# Codex autonomous brief — FEA_AP-PLY implementation

You are running autonomously as `codex exec --dangerously-bypass-approvals-and-sandbox` on a Linux compute host. Implement and validate all 16 stages of the FEA test ladder defined in this repository, in order, committing and pushing progress after each stage. This brief is the only message you will receive — read it fully before starting.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. Author for written documentation: `J.C. Vaught`. Do not mention any AI provider in commits, comments, or docs.

## 1. Required reading (in this order)

1. `README.md` — repo overview.
2. `plan/master_plan.md` — single-tool OpenRadioss + all-solid commitment, 16-stage progression, toolchain.
3. `plan/consolidation_review.md` — known cross-stage consistency items A-H. **Apply items A, B, C, D before / during the relevant stages**:
   - A: standardize the IM7/8552 card on Soden 1998 + WWFE-II 2013 strengths + Camanho-Davila 2003 cohesive. Save canonical card to `references/material_cards/im7_8552.json` and have stages reference it.
   - B: standardize fiber-orientation mechanism on `/SKEW/FIX` referenced via `Iskew` in `/PROP/TYPE14`.
   - C: `/INTER/TYPE2` requires `Spotflag=25` on modern OpenRadioss for traction-separation cohesive. **Probe the starter version on first run**; if cohesive form unavailable, auto-fall-back to `/MAT/LAW117` zero-thickness solid cohesive elements throughout stages 12-16.
   - D: `references/test_progression_literature.md` has stages 13/14 swapped relative to `master_plan.md`. Fix the literature doc so §13=LVI, §14=CAI to match the master plan.
4. `plan/kok_port_plan.md` — Kok geometry port plan. **Implement only after stages 1-10 pass.** Implementation lives at `src/kok_geom/` and is required for stages 11 and 16.
5. `references/openradioss_endtoend_audit.md` — confirms what OpenRadioss can / cannot do per stage; honor the MARGINAL caveats.
6. `tests/stage_NN_<name>/spec.md` — the per-stage specification (geometry, mesh, BCs, material, deck skeleton, reference solution, pass criterion). Read at the start of each stage.

## 2. Compute environment

```bash
# conda env (already installed; activate before any python invocation)
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
# python 3.12, gmsh 4.15, meshio 5.3, pyvista 0.47, numpy, scipy, jinja2, pytest, openmpi, vtk

# OpenRadioss runtime env vars (export before EVERY starter/engine invocation)
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH

# Binaries
$OR/exec/starter_linux64_gf       # preprocessor
$OR/exec/engine_linux64_gf        # solver
$OR/exec/engine_linux64_gf_ompi   # MPI variant
$OR/exec/engine_linux64_gf_sp     # single-precision variant
$OR/exec/anim_to_vtk_linux64_gf   # animation -> ASCII VTK

# Known-good .rad format reference (use these as templates; fixed-width columns matter)
ls /mnt/storage/j-vaught/openradioss/OpenRadioss-src/qa-tests/
```

The host is `comech-2422`, Ubuntu 22.04, 112 cores, 502 GB RAM, 4× NVIDIA RTX 6000 Ada (not used by OpenRadioss). 10.8 TB free at `/mnt/storage/j-vaught/`.

## 3. Per-stage workflow

For each stage `NN ∈ {01..16}`, in order:

1. `git pull --rebase` (sync with any human edits).
2. Read `tests/stage_NN_<name>/spec.md` end to end.
3. Inspect the existing `tests/stage_NN_<name>/runner.py` skeleton — most have explicit `NotImplementedError` TODOs (Jinja deck templates, GMSH mesh builders, T01 parsers). Fill those in. Do **not** rewrite skeletons that already work.
4. Implement the deck template by **copying skeleton structure from `qa-tests/`** and modifying. `.rad` files use fixed-width columns (10I10 for ints, F20.x for floats); do NOT hand-roll format from documentation.
5. Implement the GMSH mesh builder using `gmsh.model.occ.*` for solid bodies. Solid HEXA8 only (TETRA10 is acceptable fallback only where HEXA cannot mesh, e.g. stage 10's circular fibers — document the choice).
6. Run end-to-end: GMSH → meshio (or direct `/NODE`+`/BRICK` block writer) → starter → engine → `.anim` (gunzip) → `anim_to_vtk` → `.vtk` → PyVista → numerical extraction → comparison vs reference.
7. Produce these artifacts under `tests/stage_NN_<name>/`:
   - `runs/<rootname>_*.rad` — generated decks (do not commit large ones; .gitignore covers it)
   - `runs/<rootname>A001.vtk` through last frame — VTK output (commit first frame only)
   - `results/results.json` — `{"stage": NN, "verdict": "PASS|FAIL|INCONCLUSIVE", "metrics": {...}, "reference": {...}, "tolerance": {...}, "git_sha": "..."}`. Always commit.
   - `results/timeseries.csv` — measured vs reference. Always commit.
   - `figures/*.typ` — Typst + CeTZ plot source (brand colors per `~/.claude/CLAUDE.md`; no matplotlib). Always commit.
   - `run.log` — full stdout+stderr from the run. Commit.
8. After each stage:
   - `git add` the new artifacts; `git commit -m "stage NN <name>: <verdict>, <one-sentence summary>"`; `git push`.
   - If `results.json["verdict"] == "PASS"`: continue to next stage.
   - If `verdict == "FAIL"` (numerics genuinely don't match reference): write `tests/stage_NN_<name>/blocker.md` documenting the failure mode. **Stop. Do not proceed.** Push the blocker.md and exit.
   - If `verdict == "INCONCLUSIVE"` (toolchain limitation, not numerical): write `blocker.md`, commit, but proceed to the next stage if the dependency chain allows.

## 4. Constraints (hard)

- **Solid elements only**: HEXA8 / TETRA10. No shells, no laminate shell elements, no `/PROP/TYPE1` shell.
- **Single tool**: OpenRadioss for the solver. No FEniCSx, no Code_Aster, no commercial codes.
- **No closed source**: MIT/BSD/Apache/LGPL/AGPL only. No Abaqus, no LS-DYNA, no Ansys.
- **No GUI**: everything from CLI / Python. No HyperView, no manual ParaView clicking.
- **No matplotlib**: all plots are Typst + CeTZ from CSV. Brand colors from `~/.claude/CLAUDE.md`. No rounded edges.
- **No fabricated numerics**: if a reference value cannot be located, mark `[UNVERIFIED]` and produce `INCONCLUSIVE`.
- **All units SI** (metres / kilograms / seconds / Newtons / Pascals) inside the deck. Honor `/UNIT` block.
- **All paths absolute** in scripts. Do not assume `cwd`.

## 5. Stage-specific notes

| # | Stage | Notes |
|---|---|---|
| 01 | beam_bending | Toolchain smoke test. Material 6061-T6 isotropic. Pass within 1% of Euler-Bernoulli for both 3-pt and 4-pt. Implicit static needs MUMPS-linked starter; if MUMPS unavailable in the prebuilt, use dynamic relaxation as documented in spec. |
| 02 | cantilever_large | Bisshopp-Drucker 1945 reference values verified in spec to 4 digits. /CLOAD distributed force, NOT /PLOAD (follower-load bias up to 13% at α=5). Co-rotational Ismstr=11 + constant-pressure Icpre=1 to fight HEXA bending lock. |
| 03 | iso_dogbone_E8 | LAW2 (Johnson-Cook with rate-effects off) chosen over LAW36 in the spec. |
| 04 | open_hole_kirsch | Primary target Howland 1929 K_tg=3.035 at 2a/W=0.10 (finite-width corrected), Kirsch K_t=3.0 secondary. Butterfly O-grid, Richardson sweep at N_θ ∈ {32, 64, 128}. |
| 05 | dogbone_damage | LAW22 ductile damage substitutes Mazars (per audit MARGINAL note). Split pass criterion: 5% RMSE pre-damage-onset (mandatory), post-peak divergence reported but NOT failed (local-CDM limitation). |
| 06 | composite_failure_criteria | LAW25 + /FAIL/TSAIWU + /FAIL/HASHIN + /FAIL/PUCK side-by-side. 5% gate on 5 principal-axis points only — off-axis envelopes are reported but not gated (criteria are not expected to agree off-axis). |
| 07 | UD_tow_D3039 | Per consolidation item B, use /SKEW/FIX referenced via Iskew (NOT in-property Phi angle). |
| 08 | ply_rotation | /SKEW/FIX per ply, sweep θ ∈ {0,15,30,45,60,75,90}. |
| 09 | laminate_solid_CLT | Per consolidation item A, re-cite Soden 1998 (numerics unchanged from CMH-17). Per item B, use /SKEW/FIX. |
| 10 | UD_mesoscale_direct | KUBC (kinematic uniform BCs), NOT periodic. TETRA10 acceptable for the curved fiber boundaries. |
| 11 | PW_mesoscale_direct | **Requires Kok geometry port to be implemented first** (see §6 below). Direct uniform-displacement BCs, target Kok 2022 E_x ≈ 53.3 GPa within 10%. |
| 12 | DCB_ENF_cohesive | Critical: probe starter version. If `/INTER/TYPE2` cohesive (`Spotflag=25`) unavailable, fall back to `/MAT/LAW117` zero-thickness solid cohesive. Process zone uses E_3=11.4 GPa (NOT E_1=161). |
| 13 | LVI_D7136 | Apply consolidation item C: same `/INTER/TYPE2` ↔ `/MAT/LAW117` fallback as stage 12. State-file output (`_S0001.sta` + `_S0001` binary + `_0000.rad` mesh) sha256-hashed and emitted as `state_handoff_manifest.json` for stage 14. |
| 14 | CAI_D7137 | Chained explicit-explicit from stage 13 final state. KE/IE ≤ 0.05 runtime check; auto-bisect velocity if it fails. |
| 15 | flat_ballistic | Hourglass control: Isolid=24 reduced 1-pt + Ihq=8 HEPH (Belytschko-Bindeman 1993). Viscous Ihq=1 explicitly forbidden — it absorbs ballistic KE. V50 reference is `[UNVERIFIED]` until digitized from open-access UofSC theses. |
| 16 | PW_panel_ballistic | Project terminal goal. Reduced 200×200 mm representative section (~692k elements, 8-24 h/shot on 16 cores, 3-7 d per V50 sweep). HPC budget aware — if a single shot would exceed 24 h on this host, document and downscale, do not abandon. V50 reference paywalled; mark `[UNVERIFIED]` and report INCONCLUSIVE if the published value cannot be located in the open theses (Vakili Rad 2020 Scholar Commons /etd/5708/, Kodagali 2023 /etd/7650/). |

## 6. Kok geometry port (after stage 10, before stage 11)

`plan/kok_port_plan.md` is the binding plan. **Clean-room from Nagelsmit 2013 Ch. 2 §§2.2-2.3 and Kok 2022 §3.1.** Do NOT read or transcribe code from `rutger-kok/ap_ply_model_creation`. Implement under `src/kok_geom/` per the planned 12-module package layout. CLI signature: `kok_geom --config kok_config.json --out panel.{msh,inp}`. Produces a `.msh` plus `orientations.json` sidecar.

Validation gate: stage 11 must reproduce Kok 2022 reported in-plane modulus E_x ≈ 53.3 GPa within 10%. If it does not, the port is wrong somewhere — debug it before stage 16.

## 7. Output expectations (what the human checks when you finish)

- `tests/stage_NN_*/results/results.json` for every stage that was attempted, with `verdict` + `metrics`.
- `tests/stage_NN_*/run.log` for every stage that was run.
- `tests/stage_NN_*/blocker.md` for any FAIL or INCONCLUSIVE.
- A top-level `RESULTS_SUMMARY.md` you write at the end (or when blocked) listing each stage's verdict in one line, total wall-clock per stage, and any open issues.
- Every artifact committed and pushed to `origin/main` continuously.

## 8. When to stop

- All 16 stages PASS — write `RESULTS_SUMMARY.md`, commit, push, and exit.
- Any stage genuinely FAILs (numerics) — stop immediately.
- A toolchain blocker that cannot be worked around in this session (e.g. requires a sudo apt install) — write `blocker.md`, commit, push, exit.
- HPC runtime would exceed 48 h for a single stage — downscale per spec.md guidance, mark INCONCLUSIVE if even downscaled cannot complete, push, continue.
- You hit your token / time budget — push everything and exit gracefully.

Begin with stage 01.
