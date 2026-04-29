# FEA_AP-PLY Results Summary

**Author.** J.C. Vaught
**Date.** 2026-04-29

| Stage | Verdict | Wall-clock | Notes |
|---|---:|---:|---|
| 01 beam_bending | PASS | 911.4 s | M1 HEXA8 solid runs passed Euler-Bernoulli gates: 3pt error 0.968%, 4pt error 0.968%. |
| 02 cantilever_large | INCONCLUSIVE | 84.4 s for passing alpha=1 probe | MUMPS-linked implicit engine is installed and reaches DMUMPS. Baseline alpha=1 passes; alpha=3/5 do not yet complete due nonlinear timestep-limit convergence. |
| 03 iso_dogbone_E8 | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 04 open_hole_kirsch | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 05 dogbone_damage | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 06 composite_failure_criteria | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 07 UD_tow_D3039 | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 08 ply_rotation | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 09 laminate_solid_CLT | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 10 UD_mesoscale_direct | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 11 PW_mesoscale_direct | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 12 DCB_ENF_cohesive | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 13 LVI_D7136 | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 14 CAI_D7137 | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 15 flat_ballistic | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |
| 16 PW_panel_ballistic | not attempted | - | Stopped after Stage 02 nonlinear convergence blocker. |

## Open Issues

Stage 02 no longer fails on missing MUMPS. The active blocker is robust completion of the `alpha=3` and `alpha=5` large-deflection implicit paths. See `tests/stage_02_cantilever_large/blocker.md` for attempted solver controls and failure modes.
