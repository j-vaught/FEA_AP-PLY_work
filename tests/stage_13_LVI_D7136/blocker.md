# Stage 13 Blocker - LVI mesh/deck pipeline stale

Author: J.C. Vaught

The stage 13 runner was attempted in two steps:

```bash
python tests/stage_13_LVI_D7136/runner.py --render-only
python tests/stage_13_LVI_D7136/runner.py --execute --n-proc 16
```

Render reached GMSH, but the mesh conversion path failed with `mesh pipeline incomplete ('pyramid')`. The generated GMSH mesh contains element types that the current meshio/Abaqus/inp2rad bridge path does not handle for the OpenRadioss deck expected by the runner. The runner therefore rendered a placeholder mesh block.

The subsequent starter attempt failed immediately because the rendered starter deck does not begin with the mandatory `#RADIOSS STARTER` card. Even after that formatting issue, the runner is stale relative to the binding post-rotation recipe: it templates `/MAT/LAW25` plus `/PROP/TYPE14` and `/SKEW`, while the verified recipe for stages 13-16 is `/MAT/LAW12` plus `/PROP/TYPE6/SOL_ORTH` plus per-element `/INIBRI/ORTHO`, with `/FAIL/HASHIN` erosion using `IFAIL_SO=1` and `PTHICKFAIL=1.0`.

No valid OpenRadioss LVI solve or stage-14 state handoff was produced.

Verdict: `INCONCLUSIVE`. The next fix is a stage-13 runner rewrite in the post-matrix style: direct structured HEXA8 mesh, LAW12 + TYPE6 plies, `/INIBRI/ORTHO` per ply element, Hashin erosion on the verified matrix row, and a valid `#RADIOSS STARTER` deck.
