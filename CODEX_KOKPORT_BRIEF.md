# Codex Kok geometry port brief — FEA_AP-PLY

You run autonomously on `comech-2422` in tmux session `codex-kok`. Your job is to **implement the standalone Python AP-PLY geometry port** specified in `plan/kok_port_plan.md`. This unblocks stages 11 (PW mesoscale) and 16 (PW panel ballistic — the project's terminal goal).

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. No AI provider mentions.

## License posture (binding)

- **Clean-room reimplementation**. Algorithm comes from publicly published papers only:
  - Nagelsmit 2013 PhD thesis (TU Delft, open-access PDF) — Ch. 2 §§2.2 (Concept) and 2.3 (Pattern Details).
  - Kok 2022 *Composites Part A* (CC-BY accepted manuscript) — §3.1 (Microscale model).
  - Vakili Rad 2020 MSc thesis and Kodagali 2023 PhD thesis (open access at scholarcommons.sc.edu) — for naming convention and parameter ranges.
- **DO NOT read or copy any code from `rutger-kok/ap_ply_model_creation` (LGPL-2.1).** Inspect only the public README and top-level file list for orientation; do not transcribe code.
- Output is **MIT-licensed**. Cite Kok's repo in the README as prior art, not as a source.

## Required reading (in order)

1. `plan/kok_port_plan.md` — your binding implementation plan (713 lines). Read it end to end. Every milestone, module, and validation gate is specified there.
2. `references/uofsc_pseudowoven.md` and `references/delft_apply.md` — architecture context and confirmed geometry parameters.
3. `tests/stage_11_PW_mesoscale_direct/spec.md` — the first downstream consumer; defines the validation gate (E_x ≈ 53.3 GPa within 10%).
4. `tests/stage_16_PW_panel_ballistic/spec.md` — the terminal consumer; UofSC P1-TWT panel configuration.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply  # has gmsh 4.15, meshio 5.3, pyvista 0.47, numpy, scipy
export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
```

## Implementation milestones (from `plan/kok_port_plan.md`)

Work strictly in order. After each milestone, commit + push.

### M1 — Single-tow path + cross-section (target: 1 week)

Done when: `kok_geom.tow.Tow` class can produce a 3D solid model of one slit-tape with:
- Configurable tape width (default 6.35 mm)
- Cured ply thickness (default 0.18 mm)
- Undulation profile per Kok 2022 Figure 4 / equation form (the χ = t/Lu undulation ratio with Lu derived per the paper; default χ = 0.09)
- Per-tow material orientation as a unit vector `(cos θ, sin θ, 0)`
- Output: a `gmsh.model.occ` solid + the orientation vector

Unit tests T1, T2, T3, T6 from `kok_port_plan.md` §10 must pass.

Commit: `kok_geom M1: single-tow path + undulation primitive`.

### M2 — Single-ply array (target: 1 week)

Done when: `kok_geom.ply.Ply` can lay down a complete ply (all parallel slit-tapes covering a finite in-plane extent) with:
- Tape spacing parameter s ∈ {1, 2, 3} controlling interlace pattern position
- Resin pockets between adjacent tapes (Boolean fragment to fill gaps)
- All tapes share the same fiber direction θ_ply

Unit tests T4, T7 must pass.

Commit: `kok_geom M2: single-ply array with resin pockets`.

### M3 — Multi-ply assembly (target: 1 week)

Done when: `kok_geom.laminate.Laminate` can stack multiple plies with:
- Per-ply fiber angle from a stacking-sequence string (e.g. `[0/+45/-45/90]`)
- Per-ply angle shift to produce the AP-PLY interlacing pattern (per Nagelsmit 2013 §2.3)
- Per-tow physical groups in the GMSH model so each tow gets its own material orientation tag in the output

Unit test T5 must pass.

Commit: `kok_geom M3: multi-ply assembly with per-tow orientation tags`.

### M4 — Full mesh export and stage 11 validation (target: 1 week)

Done when:
- `kok_geom` exposes a CLI: `kok_geom --config kok_config.json --out panel.{msh,inp}`
- Output: GMSH MSH4 file + `orientations.json` sidecar mapping each per-tow physical group to its `(Vx, Vy, Vz)` unit vector
- meshio handles the .msh → .inp conversion (orientations passed via the sidecar JSON since meshio doesn't preserve per-element orientation)
- Integration test I1 (parser round-trip) and I2 (stage 11 mock-run) pass
- **The output is then plugged into `tests/stage_11_PW_mesoscale_direct/runner_richviz.py` (or runner.py) — re-run stage 11. Goal: E_x ≈ 53.3 GPa within 10% per Kok 2022.**

Unit tests T8, T9, T10 must pass.

Commit: `kok_geom M4: full pipeline + stage 11 PASS`.

If stage 11 PASSes, also flip its `results/results.json` from INCONCLUSIVE to PASS and commit `stage 11 PW_mesoscale_direct: PASS, post-Kok-port`. Push.

## Module layout

Per `plan/kok_port_plan.md` §4, under `src/kok_geom/`:

```
src/kok_geom/
  __init__.py
  config.py        # pydantic schema for kok_config.json
  geometry/
    tow.py         # M1: single-tow path + undulation
    ply.py         # M2: single-ply array
    laminate.py    # M3: multi-ply assembly
  occ_backend.py   # gmsh.model.occ wrapper layer
  mesh.py          # GMSH meshing + physical group emission
  io.py            # MSH4 export + orientations.json sidecar
  cli.py           # `kok_geom` entry point
tests/
  test_tow.py
  test_ply.py
  test_laminate.py
  test_io.py
  test_cli.py
```

`pyproject.toml` (new) declares the `kok_geom` package and `kok_geom` console_script.

## Constraints (hard)

- Pure Python. No Fortran, no C extensions, no compiled binaries.
- Solid HEXA8 / TETRA10 only in the output mesh.
- MIT license; do not transcribe Kok's LGPL code.
- gmsh.model.occ.* is the geometry kernel.
- All units SI in code; expose mm/m as user choice in the JSON config.
- pytest passes for every milestone before commit.
- Author = `j-vaught <jvaught@sc.edu>`.

## Tile-and-fragment scaling strategy (per `kok_port_plan.md` §5)

For panels with > 100 tows, use the documented tile-and-fragment approach:

1. Build one tape as a single OCC primitive.
2. Array-replicate to a unit cell of the interlace pattern.
3. Boolean fragment within the tile only (small cascade).
4. Translate-replicate the tile to fill the panel.
5. Boolean fragment between adjacent tiles only at their interface.

This avoids OCC choking on a single 10⁴-tow boolean cascade. Implement it as part of M3 / M4 if the M3 small-block test takes more than a few minutes per build.

## When to stop

- Reaches M4 done — stage 11 PASSes — write `RESULTS_KOKPORT_SUMMARY.md`, commit, push, exit.
- Hard blocker: gmsh.model.occ refuses an operation, or pydantic / meshio APIs need patching beyond reasonable scope. Write a `blocker.md` under `src/kok_geom/`, commit, push, exit.
- Token / time budget — push everything in flight (even partial milestones), exit gracefully.

Begin with M1.
