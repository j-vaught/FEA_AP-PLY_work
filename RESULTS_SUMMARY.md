# FEA_AP-PLY Results Summary

**Author.** J.C. Vaught
**Date.** 2026-04-29

| Stage | Verdict | Wall-clock | Notes |
|---|---:|---:|---|
| 01 beam_bending | PASS | 911.4 s | M1 HEXA8 solid runs passed Euler-Bernoulli gates: 3pt error 0.968%, 4pt error 0.968%. |
| 02 cantilever_large | INCONCLUSIVE | 3.5 s explicit probe | Nonlinear implicit engine reports `Fatal error: MUMPS required`; explicit fallback probe ran but was not a gated quasi-static verification. |
| 03 iso_dogbone_E8 | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 04 open_hole_kirsch | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 05 dogbone_damage | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 06 composite_failure_criteria | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 07 UD_tow_D3039 | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 08 ply_rotation | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 09 laminate_solid_CLT | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 10 UD_mesoscale_direct | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 11 PW_mesoscale_direct | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 12 DCB_ENF_cohesive | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 13 LVI_D7136 | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 14 CAI_D7137 | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 15 flat_ballistic | not attempted | - | Stopped after Stage 02 toolchain blocker. |
| 16 PW_panel_ballistic | not attempted | - | Stopped after Stage 02 toolchain blocker. |

## Open Issues

Stage 02 needs either a MUMPS-linked OpenRadioss build for `/IMPL/NONLIN` or a validated long-run explicit dynamic-relaxation plan with kinetic-energy checks. The current prebuilt engine cannot run the intended nonlinear implicit verification.
