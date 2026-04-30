# Stage 16 Phase A Diagnostics - 200 mm M5 mesh clears the geometry gate

Author: J.C. Vaught

The old pre-M5 Stage 16 blocker is superseded. Post-Kok-M5, the full reduced
`200 mm x 200 mm` P1-TWT section now meshes successfully within the one-hour
Phase A cap.

Configuration updates applied before the rerun:

- `panel.symmetry = "midplane"` to reflect the symmetric 24-ply stack.
- `laydown.tow_coverage_fraction = 0.94` explicitly, matching the working
  post-M5 Stage 11 geometry rather than relying on the older default.
- Mesh settings remain the coarse Phase A preflight values:
  `2.0 mm` in-plane and through-thickness targets, `element_order = 2`.

Measured command:

```bash
/usr/bin/time -p timeout 3600 python -m kok_geom \
  --config tests/stage_16_PW_panel_ballistic/geometry/kok_p1_twt_200_config.json \
  --out tests/stage_16_PW_panel_ballistic/geometry/p1_twt_200.msh
```

Measured result:

- Wall-clock: `402.98 s` = `6.72 min`.
- Nodes: `1,979,649`.
- TETRA10 elements: `1,440,000`.
- Orientation groups: `9,844`.
- Undulation records: `4,466`.
- Output sizes: `272 MiB` `.msh`, `261 MiB` `.inp`, `8.9 MiB` orientations JSON.

Interpretation:

- The 200 mm representative section now clears the Phase A geometry gate with
  large margin against both the `60 min` hard cap and the `30 min` target.
- The old `timeout 900` result is no longer relevant for Stage 16 budgeting.
- The next decision gate must come from a **measured** Phase B shot on this
  post-M5 200 mm mesh at `-nt 32`, not from any more geometry extrapolation.
