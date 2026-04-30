# Stage 16 Reduced-Section Status - Phase B1 gate clears on 75 mm deck

Author: J.C. Vaught

The old `200 mm x 200 mm` post-M5 blocker remains useful as the full-section
upper-bound measurement, but it is no longer the active gate for this pass.
Part 1 of the reduced-section retry was executed with both requested levers:

1. Section size reduced from `200 mm x 200 mm` to `75 mm x 75 mm`.
2. Uniform in-plane mesh relaxed from `2.0 mm` to `3.0 mm`.

## Phase A1 - reduced geometry build

Command run:

```bash
python -m kok_geom \
  --config tests/stage_16_PW_panel_ballistic/geometry/kok_p1_twt_75_config.json \
  --out tests/stage_16_PW_panel_ballistic/geometry/p1_twt_75.msh
```

Measured Phase A1 output:

- Mesh wall-clock: `58.36 s` = `0.97 min`.
- Mesh size: `290,521` nodes, `207,936` `TETRA10`.
- Orientation sidecar groups: `1,633`.
- Undulation metadata records: `638`.

The reduced geometry therefore clears the `<= 15 min` meshing gate easily.

## Phase B1 - reduced single-shot wall-clock probe

The Stage 16 runner was extended with a narrow `--stop-at-cycle` mode so the
engine can be measured at a clean checkpoint without wasting hours on a full
shot during the gate test.

Command run:

```bash
python tests/stage_16_PW_panel_ballistic/runner.py \
  --mesh tests/stage_16_PW_panel_ballistic/geometry/p1_twt_75.msh \
  --orientations tests/stage_16_PW_panel_ballistic/geometry/orientations.json \
  --run-dir tests/stage_16_PW_panel_ballistic/runs/stage_16_post_m5/phase_b_single_shot_75_c2000 \
  --threads 32 \
  --stop-at-cycle 2000
```

Measured reduced-deck counts:

- Total nodes after interface duplication: `427,509`.
- Composite `TETRA10`: `195,348`.
- Resin `TETRA10`: `12,588`.
- Projectile `TETRA4`: `2,496`.
- `/INTER/TYPE2` interfaces: `23`.
- Interface surface quads: `33,212`.

Measured solver timings:

- Starter return code: `0`.
- Starter wall-clock: `12.64 s`.
- Engine checkpoint: `NC = 2000`, `T = 3.0867e-06 s`, `DT = 1.5433e-09 s`.
- Engine elapsed at checkpoint: `78.55 s`.
- Engine remaining at checkpoint: `12,645.32 s` = `3.51 h`.
- Projected total shot wall-clock from the checkpoint: `3.53 h`.

Interpretation:

- The reduced `75 mm x 75 mm`, `3 mm` mesh is **not compute-bound** against the
  `6 h/shot` decision gate on the available `-nt 32` workstation path.
- Phase B1 therefore clears and the next required step is **Phase C1 V50
  sweep**, not further geometric reduction.
- Stage 16 remains `INCONCLUSIVE` in this pass only because the V50 sweep was
  not completed yet; the wall-clock blocker moved from "single shot infeasible"
  to "reduced shot feasible, calibration sweep still pending."
