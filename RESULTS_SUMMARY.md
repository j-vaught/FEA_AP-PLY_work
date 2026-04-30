# FEA_AP-PLY Results Summary

**Author.** J.C. Vaught
**Date.** 2026-04-30

| Stage | Verdict | Wall-clock / attempt | Notes |
|---|---:|---:|---|
| 01 beam_bending | PASS | 911.4 s | M1 HEXA8 solid runs passed Euler-Bernoulli gates: 3pt error 0.968%, 4pt error 0.968%. |
| 02 cantilever_large | INCONCLUSIVE | 84.4 s for passing alpha=1 probe | MUMPS-linked implicit engine is installed and reaches DMUMPS. Baseline alpha=1 passes; alpha=3/5 remain blocked by nonlinear timestep-limit convergence. |
| 03 iso_dogbone_E8 | PASS | 78.6 s | Medium HEXA8 solid dogbone load sweep passed uniform stress, gauge strain, yield-onset, and apparent-modulus gates. |
| 04 open_hole_kirsch | PASS | 7.5 s | Ntheta=64 implicit solid O-grid met the Howland Kt gate: Kt=3.0956 vs 3.035, error 1.996%; far-field stress error 0.585%. |
| 05 dogbone_damage | PASS | 2376.0 s | LAW22 notched-dogbone mesh ladder passed pre-onset RMSE gates: coarse-medium 4.843%, medium-fine 3.736%; medium ligament EPSP max 0.1229 exceeded the 0.05 damage-onset threshold. |
| 06 composite_failure_criteria | PASS | 5.4 s | LAW12 + TYPE6/SOL_ORTH solid composite decks start, run, convert to VTK, and register TSAIWU, HASHIN, and PUCK failure cards on the verified matrix row. |
| 07 UD_tow_D3039 | PASS | 9.6 s | LAW12 + TYPE6/SOL_ORTH decks use the verified `Ip=3`, `Iorth=0`, `Phi=theta` property convention. 0 deg error 0.403%; 45 deg off-axis error 0.114% using engineering strain from displacement. |
| 08 ply_rotation | PASS | 14.0 s | Seven-angle LAW12 + TYPE6/SOL_ORTH sweep passed the analytic Ex(theta) gate; max error 0.815%. VTK cell `Stra[0]` is not used for the gate. |
| 09 laminate_solid_CLT | PASS | 35.1 s | Cross-ply `[0/90]s` and quasi-isotropic `[0/+45/-45/90]s` solid-ply stacks passed the CLT A-matrix gate; max component error 1.085%. |
| 10 UD_mesoscale_direct | INCONCLUSIVE | 6.3 s | Self-contained 20x20x4 isotropic fiber/matrix HEXA8 voxel cell runs, but KUBC bias dominates the Halpin-Tsai gate; transverse error 11.189%, shear error 99.934%. |
| 11 PW_mesoscale_direct | FAIL | 792.1 s | Kok CLI stage-11 block runs through OpenRadioss with LAW12 + TYPE6 + `/INIBRI/ORTHO`, but stiffness misses Kok 2022 targets: Ex 16.706 GPa, Ey 16.968 GPa, Gxy 0.0229 GPa. |
| 12 DCB_ENF_cohesive | INCONCLUSIVE | 161 s bounded attempt | Existing runner did not reach starter; full-resolution GMSH/deck path and cohesive templating need rewrite against the verified LAW12 + TYPE6 recipe. |
| 13 LVI_D7136 | INCONCLUSIVE | bounded render + starter attempt | Render reached GMSH but mesh conversion failed on `pyramid`; starter then rejected the placeholder deck because `#RADIOSS STARTER` was missing. No state handoff. |
| 14 CAI_D7137 | INCONCLUSIVE | immediate dependency stop | Correctly stopped before deck rendering because Stage 13 did not produce a `lvi_d7136_*.sta` final damaged state. |
| 15 flat_ballistic | INCONCLUSIVE | single 300 m/s starter attempt | Runner rendered an invalid native OpenRadioss deck with stale LAW25 + TYPE14 syntax and missing mandatory first card. V50 sweep not reached. |
| 16 PW_panel_ballistic | INCONCLUSIVE | 900 s 200 mm timeout; 12 min 100 mm preflight | 200 x 200 mm Kok build timed out. 100 x 100 mm downscale wrote a coarse TETRA10 mesh with 1,120,178 nodes and 603,481 elements, giving a lower-bound estimate of 4.19 h per shot before contact/damage overhead. |

## Open Issues

Stage 02 no longer fails on missing MUMPS. The active blocker is robust completion of the `alpha=3` and `alpha=5` large-deflection implicit paths. See `tests/stage_02_cantilever_large/blocker.md`.

Stage 10 is classified as `INCONCLUSIVE`, not `FAIL`, because the current direct KUBC cell is being compared against Halpin-Tsai RVE targets. A periodic or mixed-boundary homogenization workflow is needed before this is a decisive material validation.

Stage 11 proves the solver-side Kok-geometry path can run, but the clean-room geometry port does not yet reproduce Kok-style physical undulation volumes or mixed tow/resin unit cells. That gap dominates the stiffness miss.

Stages 12, 13, and 15 need runner rewrites to the verified LAW12 + TYPE6/SOL_ORTH pattern. Stage 14 is blocked only by the missing Stage 13 restart state. Stage 16 needs graded meshing and HPC/distributed OpenRadioss before a defensible V50 sweep.
