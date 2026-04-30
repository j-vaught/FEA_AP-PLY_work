# Stage 15 Blocker - flat ballistic deck generator invalid

Author: J.C. Vaught

Stage 15 was attempted with a single 300 m/s strike:

```bash
OPENRADIOSS_STARTER=/mnt/storage/j-vaught/openradioss/OpenRadioss/exec/starter_linux64_gf \
OPENRADIOSS_ENGINE=/mnt/storage/j-vaught/openradioss/OpenRadioss/exec/engine_linux64_gf \
OPENRADIOSS_NT=16 \
python tests/stage_15_flat_ballistic/runner.py --single 300
```

The runner rendered `tests/stage_15_flat_ballistic/runs/Vs300/flatballistic_Vs300_0000.rad` and invoked starter, but starter returned code `2` immediately:

```text
ERROR ID : 100201
File is not a valid deck, mandatory first card "#RADIOSS STARTER" missing.
```

The generated deck is not currently a valid native OpenRadioss deck. It uses string IDs such as `/MAT/LAW25/MAT_PLY01`, `key = value` field syntax, a missing generated mesh include (`mesh.rad`), and the stale `/MAT/LAW25` + `/PROP/TYPE14` composite path. The required post-matrix recipe for this stage is `/MAT/LAW12` + `/PROP/TYPE6/SOL_ORTH` + `/INIBRI/ORTHO` with `/FAIL/HASHIN` erosion (`IFAIL_SO=1`, `PTHICKFAIL=1.0`) and HEPH-compatible solid-hourglass settings on a valid HEXA8 mesh.

Verdict: `INCONCLUSIVE`. The next fix is a deck-generator rewrite before any V50 sweep can be run.
