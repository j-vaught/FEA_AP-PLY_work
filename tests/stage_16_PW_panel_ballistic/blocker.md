# Stage 16 Phase B Diagnostics - 200 mm M5 shot is compute-bound at `-nt 32`

Author: J.C. Vaught

The old pre-M5 Stage 16 blocker is superseded twice over:

1. Phase A proved that the post-M5 `200 mm x 200 mm` P1-TWT section meshes in
   minutes, so geometry is no longer the gating problem.
2. Phase B proved that the actual OpenRadioss shot on that real post-M5 deck is
   still compute-bound at `-nt 32`, but now with a measured live solver
   trajectory instead of the stale pre-M5 extrapolation.

Phase A geometry facts retained from the accepted rerun:

- `panel.symmetry = "midplane"` to reflect the symmetric 24-ply stack.
- `laydown.tow_coverage_fraction = 0.94` explicitly, matching the working
  post-M5 Stage 11 geometry rather than relying on the older default.
- The accepted Phase A mesh itself remains:
  `1,979,649` nodes, `1,440,000` TETRA10, `402.98 s` (`6.72 min`) wall-clock.

Phase B deck path:

- Replaced the stale pre-M5 `runner.py` with a direct post-M5 deck writer that
  reads `geometry/p1_twt_200.msh` plus `geometry/orientations.json`.
- The writer duplicates only the `23` interior ply planes, remaps the upper-ply
  connectivity on each interface, and emits:
  `LAW12 + TYPE6 + /INIBRI/ORTHO + /FAIL/HASHIN`, `23` `/INTER/TYPE2` cards,
  a rigid tetra4 steel projectile, `/INTER/TYPE7` impact contact, and clamped
  cut-edge BCs.
- Resulting starter-deck size and counts:
  `605 MiB` starter deck, `2,909,493` total nodes after interface duplication,
  `1,343,052` composite TETRA10, `96,948` resin TETRA10, `2,496` projectile
  TETRA4, and `230,000` interface surface quads.

Measured Phase B live run:

```bash
python tests/stage_16_PW_panel_ballistic/runner.py
```

Measured starter result on the actual Phase B run:

- Starter return code: `0` (PASS).
- Starter wall-clock: `150.82 s` = `2.51 min`.

Measured engine checkpoints on the same `-nt 32` run before manual stop:

- `NC = 100`, `T = 1.4615e-07 s`, elapsed `45.93 s`, remaining
  `157101.49 s` = `43.64 h`.
- `NC = 200`, `T = 2.9229e-07 s`, elapsed `81.16 s`, remaining
  `138759.07 s` = `38.54 h`.
- `NC = 300`, `T = 4.3844e-07 s`, elapsed `116.52 s`, remaining
  `132762.79 s` = `36.88 h`.
- `NC = 400`, `T = 5.8458e-07 s`, elapsed `152.31 s`, remaining
  `130115.46 s` = `36.14 h`.
- `NC = 500`, `T = 7.3073e-07 s`, elapsed `187.42 s`, remaining
  `128051.43 s` = `35.57 h`.

Manual-stop bookkeeping for that run:

- Runner-reported starter wall-clock: `150.82 s`.
- Runner-reported engine wall-clock before stop: `211.89 s`.
- Runner total wall-clock before stop: `362.71 s` = `6.05 min`.
- Engine return code after manual kill: `3`.

Interpretation:

- The deck is no longer blocked on syntax or geometry. The full post-M5 Stage 16
  stack starts cleanly through starter and enters the real engine on `32`
  threads.
- The compute bottleneck is now measured on the real deck, not inferred from
  the old pre-M5 mesh. By `NC=500`, the live engine is still projecting
  `35.57 h` remaining, before even reaching `0.001 ms` of simulated time.
- Even if the remaining-time estimate improved materially after contact, it is
  already far outside the `6 h/shot` Phase B gate. Launching the Stage 16 sweep
  would be irresponsible use of compute.
- Stage 16 therefore stops here with:
  `INCONCLUSIVE: compute-bound, measured on the post-Kok-M5 200 mm deck at -nt 32`.
