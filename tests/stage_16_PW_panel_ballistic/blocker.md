# Stage 16 Blocker - Kok geometry and V50 sweep exceed workstation budget

Author: J.C. Vaught

The reduced 200 mm x 200 mm P1-TWT Kok geometry was attempted with the clean-room `kok_geom` CLI:

```bash
python -m kok_geom \
  --config tests/stage_16_PW_panel_ballistic/geometry/kok_p1_twt_200_config.json \
  --out tests/stage_16_PW_panel_ballistic/geometry/p1_twt_200.msh
```

The 200 mm build was capped with `timeout 900` and did not write a mesh within 15 minutes.

Per the brief, I downscaled to a 100 mm x 100 mm representative section with the same 24-ply `[+45/90/-45/0]`, `10001000`, 6.35 mm tape, 4.55 mm total-thickness P1-TWT definition:

```bash
python -m kok_geom \
  --config tests/stage_16_PW_panel_ballistic/geometry/kok_p1_twt_100_config.json \
  --out tests/stage_16_PW_panel_ballistic/geometry/p1_twt_100.msh
```

That downscaled preflight completed after roughly 12 minutes at a coarse 2 mm global TETRA10 mesh. The resulting mesh has:

- Nodes: `1,120,178`.
- TETRA10 elements: `603,481`.
- Physical groups: `120`.
- Orientation-sidecar groups: `120`.
- Undulation metadata records: `207`.

This is still outside the practical V50 sweep budget. Stage 11 measured about `258 s` for one `0.2 ms` explicit case on `25,918` TETRA10 elements. Scaling by element count and the `0.5 ms` Stage 16 duration gives a coarse lower-bound estimate:

```text
258 s * (603481 / 25918) * (0.5 ms / 0.2 ms) = 15074 s = 4.19 h per shot
```

That estimate excludes projectile contact, Hashin erosion, state/history output, and any mesh refinement beyond the coarse 2 mm preflight. The required spec mesh has 0.5 mm impact-zone sizing; because the current `kok_geom` mesher does not yet implement graded-zone sizing, enforcing 0.5 mm globally would increase element count by approximately `16x`, pushing a single shot well beyond the workstation budget.

`V50_published` remains `[UNVERIFIED]` in the stage spec. Per the brief, the rough center would be `100 m/s`, but no V50 sweep was launched because even the 100 mm downscaled preflight exceeds the per-shot wall-clock budget.

Verdict: `INCONCLUSIVE`. Recommendation: implement tile-and-fragment geometry scaling plus true graded meshing, then run Stage 16 on HPC or a distributed OpenRadioss build. The current 100 mm preflight mesh is useful for element-count planning only, not for validation.
