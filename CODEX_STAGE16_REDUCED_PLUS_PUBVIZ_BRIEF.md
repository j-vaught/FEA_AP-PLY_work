# Codex stage 16 reduced-section + final pubviz brief — post-Kok-M5

You run autonomously on `comech-2422` in tmux session `codex-finish`. This is a combined two-part finish pass after the Kok M5 work landed:

- **Stage 11 PW_mesoscale_direct PASS** at commit `566c7be` (Ex 49.04 GPa, Ey 53.29 GPa, Gxy 22.53 GPa).
- **Stage 16 PW_panel_ballistic INCONCLUSIVE-compute-bound** at commit `584177d` — measured 35.57 h/shot at `-nt 32` on the 200 mm M5 deck.

Two things to finish.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. No AI provider mentions.

## Hard SMP scaling finding (binding — do not waste compute on this)

A `-nt 100` probe on the same Phase B deck was just measured. SMP scaling is **flat** between 32 and 100 cores:

```
NC=100   -nt 32:  ELAPSED 45.93 s   REMAINING 157101.49 s
NC=100   -nt 100: ELAPSED 45.30 s   REMAINING 154924.76 s

NC=300   -nt 32:  ELAPSED 116.52 s  REMAINING 132762.79 s   (~36.9 h)
NC=300   -nt 100: ELAPSED 115.52 s  REMAINING 131629.47 s   (~36.6 h)
```

Implication: this 1.44 M-element + 23-cohesive-interface deck is memory-bandwidth-bound, not core-bound. **Throwing more cores will not reduce wall-clock**. Adding cores beyond 32 is not a path to V50 sweep feasibility.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH
```

Use `-nt 32` for solver runs (no benefit beyond that per the SMP finding above).

---

## Part 1 — Stage 16 reduced-section retry

Reduce the *element count* (not the cores) to bring single-shot wall-clock under the 6 h/shot decision gate. Two independent levers, applied together:

### Lever A — Smaller representative section

Spec.md §2.2 documents the Saint-Venant rationale: minimum section side is **8× tape width = 51 mm**. The 200 mm choice was conservative buffer. Drop to **75 mm × 75 mm** (still 1.5× the 51 mm minimum and centered on the impact spot). Update:

- `tests/stage_16_PW_panel_ballistic/geometry/kok_p1_twt_75_config.json` — copy of `kok_p1_twt_200_config.json` with `size_x_mm=75`, `size_y_mm=75`, all other fields preserved (24 plies, midplane symmetry, tape spacing 3, tow_coverage_fraction 0.94).
- The clamp boundary moves with the section. The strike-zone definition stays at the geometric center.
- Section-size geometric scaling: 75/200 = 0.375 → element count drops ~`0.375²` ≈ 14%, so naive estimate is `1.44M × 0.14 ≈ 200k TETRA10`.

### Lever B — Coarser global mesh

The 2 mm in-plane mesh resolution was conservative. Drop to **3 mm in-plane** for the reduced-section deck. This is still finer than ply thickness (0.190 mm × 24 = 4.55 mm panel) and resolves the tape width (6.35 mm). Update:

- In `kok_p1_twt_75_config.json` set `mesh.in_plane_target_mm_impact_zone = 3.0` (matches `mesh.in_plane_target_mm` to keep uniform globally; we don't have graded meshing).
- Through-thickness should remain one element per ply minimum (per spec).

### Phase A1 — geometry meshing time-budget

```
time python -m kok_geom \
  --config tests/stage_16_PW_panel_ballistic/geometry/kok_p1_twt_75_config.json \
  --out tests/stage_16_PW_panel_ballistic/geometry/p1_twt_75.msh
```

Gate: ≤15 min wall-clock to mesh. Capture node/element count. Commit `stage 16 phase A1: 75 mm reduced + 3 mm mesh, <N> nodes / <M> TETRA10 in <T> min`.

### Phase B1 — single-shot wall-clock measurement

Reuse the existing `runner.py run_phase_b` plumbing — point it at the new mesh + orientations sidecar. Run starter + engine with `-nt 32`. Run the engine for **at least NC=2000** to get a clean projected `REMAINING TIME` value, then check the gate.

- If projected total wall-clock ≤ 6 h/shot → proceed to Phase C1.
- If 6 < T ≤ 12 h/shot → run the **5-shot reduced sweep** (multipliers 0.85, 0.95, 1.00, 1.05, 1.15).
- If T > 12 h/shot → reduce further: try **50 mm × 50 mm + 4 mm mesh** as Lever C, repeat Phase A1/B1.
- If 50 mm + 4 mm still fails → INCONCLUSIVE-compute-bound terminal.

Commit `stage 16 phase B1: <reduced-section spec>, single-shot measured <T> h on <N>-elem deck`.

### Phase C1 — V50 sweep (only if Phase B1 clears 12 h/shot)

V50 multipliers around `V50_published ≈ 100 m/s` (see spec.md §2). Run, parse residual velocity from T01, fit Recht-Ipson 1963, compute simulated V50, compare to published with 7% tolerance (9% if BC-insensitivity check fails).

Commit `stage 16 PW_panel_ballistic: <verdict>, post-Kok-M5 reduced section`.

---

## Part 2 — Publication-grade figures for stages where runs worked

Regenerate / upgrade figures for stages with PASS verdicts and post-M5 deltas. **Read** `~/.claude/CLAUDE.md` for brand colors (Garnet `#73000A`, Horseshoe `#65780B`, etc.); **no rounded corners**, **no matplotlib** (Typst+CeTZ for 2D; PyVista headless OSMesa for 3D; ffmpeg for MP4).

Stages 01-09 (PASS) already have publication-grade composite + HD motion figures (committed `213e6c3`, `5d88d1d`, `8a6f841`, `ebbfa41`, `8937061`, `331e12d`, `6519795`, `cad507c`, `039f92d`). **Do not re-render** unless the existing files are missing or broken.

### Stage 11 (post-M5 PASS) — full upgrade required

Stage 11 currently has only a wireframe blocker figure from when it FAILed (commit `5d88d1d`). Now that it PASSes with the new geometry, replace with publication-grade artifacts:

- `tests/stage_11_PW_mesoscale_direct/figures/stage_11_composite.png` — 2400×1600 four-panel:
  - Top-left: deformed mesh under axial-x BC, coloured by displacement magnitude (mm). Use brand garnet → atlantic gradient.
  - Top-right: same warp coloured by von Mises stress (MPa).
  - Bottom-left: same warp coloured by axial fiber-direction stress σ₁ (MPa) for the 0° tows.
  - Bottom-right: 2D Typst+CeTZ bar chart of `(measured, target, error%)` for `E_x`, `E_y`, `G_xy`. Header band: `Stage 11 PW_mesoscale_direct — PASS — post-Kok-M5`.
  - Each colourbar with units and brand colors.
- `tests/stage_11_PW_mesoscale_direct/figures/stage_11_motion_hd.mp4` — 1920×1080, 30 fps, ~10 s, H.264. Three-panel synchronized time evolution: (axial-x, transverse-y, shear-xy) deformed meshes coloured by displacement; moving time cursor at the bottom. Title card 1 s; end card 1 s with verdict and the three measured moduli.
- If the stage 11 VTK output is single-frame (because the runner is quasi-static), regenerate the deck with `/ANIM/DT` set to write 50+ frames into `runs/<rootname>_pubviz/` (do **not** touch the verification artifacts under `runs/stage11_*A001.vtk`).

### Stage 16 — Phase A geometry render

Add `tests/stage_16_PW_panel_ballistic/figures/stage_16_phase_a_geometry.png` — 2400×1600, two-panel:

- Left: full 200 mm M5 panel mesh, isometric view, colored by ply (24 plies → 4 angle classes → 4 brand colors: garnet 0°, atlantic +45°, horseshoe -45°, congaree 90°). Cohesive interfaces optionally shown as thin contrast lines.
- Right: zoom-in on the impact-center 25 mm × 25 mm patch with the same coloring, so the AP-PLY interlock pattern is visible.
- Footer band: `Stage 16 — Phase A — 200 mm M5 mesh, 1.44 M TETRA10, post-Kok-M5`.

If Phase B1/C1 land VTK animation output, also produce `stage_16_phase_b_motion_hd.mp4` for the single-shot impact (panel back-face deflection in time, projectile penetration depth as overlay).

### Cross-stage summary update

Update `figures/cross_stage_summary.typ` (existing) to flip stage 11 from FAIL (red) to PASS (green) and stage 16 from `INCONCLUSIVE-toolchain` (amber) to either `PASS`, `INCONCLUSIVE-compute-bound` (amber, with new measured 35.57 h/shot annotation), or whatever Phase C1 resolves to. Recompile the PDF.

Also update `figures/discovery_kok_port.typ` to flip the M1-M4 timeline → M1-M5, with the stage 11 PASS arithmetic on the right (`Ex 49.04 / Ey 53.29 / Gxy 22.53 GPa`).

### RESULTS_REPORT.md

Update the report to reference the new stage 11 figures and the stage 16 phase A geometry render. Replace the stage 11 wireframe blocker link with the new composite. Update the final tally table to show stage 11 PASS, stage 16 INCONCLUSIVE-compute-bound (with the measured wall-clock annotation if Part 1 doesn't land PASS).

Commit per artifact group: `pubviz stage 11: composite + HD motion (post-M5)`, `pubviz stage 16: phase A geometry render`, `pubviz: cross-stage summary + discovery figures (M5)`, `final: results report figure links (M5)`.

---

## Constraints (hard, unchanged)

- Solid HEXA8 / TETRA10 only.
- Pure Python; gmsh.model.occ.* for any new geometry; MIT-licensed.
- LAW12 + TYPE6/SOL_ORTH + INIBRI/ORTHO + FAIL/HASHIN + INTER/TYPE2.
- Brand colors strict; no rounded corners; no matplotlib.
- Typst + CeTZ for 2D plots; PyVista headless OSMesa for 3D renders; ffmpeg for MP4.
- No GUI, no closed-source.
- Author = `j-vaught <jvaught@sc.edu>`.
- Do not change verification verdicts (results.json) of already-PASS stages 01-09.

## When to stop

- Part 1 lands PASS or measured-INCONCLUSIVE on the reduced section → commit and proceed to Part 2.
- Part 2 lands all required figures → commit and exit.
- Token / time budget — push everything in flight, exit gracefully.

Begin with Part 1 Phase A1 (build the 75 mm reduced-section config + mesh).
