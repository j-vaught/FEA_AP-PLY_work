# Stage 12 Blocker - cohesive runner mesh/deck path not bounded

Author: J.C. Vaught

The existing stage 12 runner did not reach an OpenRadioss starter deck within a bounded attempt. Command used:

```bash
python tests/stage_12_DCB_ENF_cohesive/runner.py --threads 16 --approach A
```

The process remained CPU-bound in the GMSH full-resolution DCB mesh build for more than 2.5 minutes before writing any file under `_work_stage12`. I terminated that attempt to preserve the phase-B coverage budget.

The current runner is also stale relative to the verified composite recipe: its templated material path still describes `/MAT/LAW25` with `/PROP/TYPE14`, while the binding project recipe is `/MAT/LAW12` with `/PROP/TYPE6/SOL_ORTH`. The script further depends on an `inp2rad` geometry path and placeholder cohesive templating that does not yet build a validated `/INTER/TYPE2` or `/MAT/LAW117` OpenRadioss deck end-to-end.

Closed-form references were evaluated successfully:

- DCB simple beam theory peak load: `150.468 N`.
- DCB modified beam theory peak load: `145.239 N`.
- ENF simple beam theory peak load: `627.922 N`.
- ENF corrected beam theory peak load: `602.805 N`.
- Mode-I process-zone length: `2.787 mm`, Turon three-element limit `0.929 mm`.

Verdict: `INCONCLUSIVE`. The next fix should rewrite this runner directly in the post-matrix style used by stages 07-11: structured mesh generated in the runner, LAW12 + TYPE6 arms, and a starter-checked cohesive `/INTER/TYPE2` deck with LAW117 solid-cohesive fallback.
