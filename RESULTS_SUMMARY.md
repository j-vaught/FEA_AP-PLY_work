# FEA_AP-PLY Results Summary

**Author.** J.C. Vaught
**Date.** 2026-04-29

| Stage | Verdict | Wall-clock | Notes |
|---|---:|---:|---|
| 01 beam_bending | PASS | 911.4 s | M1 HEXA8 solid runs passed Euler-Bernoulli gates: 3pt error 0.968%, 4pt error 0.968%. |
| 02 cantilever_large | INCONCLUSIVE | 84.4 s for passing alpha=1 probe | MUMPS-linked implicit engine is installed and reaches DMUMPS. Baseline alpha=1 passes; alpha=3/5 do not yet complete due nonlinear timestep-limit convergence. |
| 03 iso_dogbone_E8 | PASS | 78.6 s | Medium HEXA8 solid dogbone load sweep passed uniform stress, gauge strain, yield-onset, and apparent-modulus gates. |
| 04 open_hole_kirsch | PASS | 7.5 s | Ntheta=64 implicit solid O-grid met the Howland Ktg gate: Kt=3.0956 vs 3.035, error 1.996%; far-field stress error 0.585%. Ligament and Richardson diagnostics are reported as non-gating finite-strip checks. |
| 05 dogbone_damage | PASS | 2376.0 s | LAW22 notched-dogbone mesh ladder passed pre-onset RMSE gates: coarse-medium 4.843%, medium-fine 3.736%; medium ligament EPSP max 0.1229 exceeded the 0.05 damage-onset threshold. |
| 06 composite_failure_criteria | PASS | 5.4 s | Post-matrix LAW12 + TYPE6/SOL_ORTH solid composite decks start, run, convert to VTK, and register TSAIWU, HASHIN, and PUCK failure cards on the verified matrix row. Principal-axis strengths match the WWFE-II card values used by the failure cards. |
| 07 UD_tow_D3039 | INCONCLUSIVE | 1.2 s | Canonical 0 deg and 45 deg D3039 LAW25 + TYPE14 decks fail at starter with ERROR 3047. TYPE6/SOL_ORTH proxy starters succeed, but no modulus gate is evaluated because the required TYPE14 path is blocked. |
| 08 ply_rotation | INCONCLUSIVE | 4.3 s | Seven-angle LAW25 + TYPE14 ply-rotation sweep fails at starter with ERROR 3047 for every angle. TYPE6/SOL_ORTH proxy starters succeed; analytic Ex(theta) table is reported without FEM modulus claims. |
| 09 laminate_solid_CLT | INCONCLUSIVE | 1.1 s | Cross-ply and quasi-isotropic per-ply LAW25 + TYPE14 laminate probes fail at starter with ERROR 3047. Analytic CLT A-matrices are reported; TYPE6/SOL_ORTH proxy starters succeed but are not used for the gate. |
| 10 UD_mesoscale_direct | not attempted | - | Pending. |
| 11 PW_mesoscale_direct | not attempted | - | Pending; Kok geometry port required after stage 10. |
| 12 DCB_ENF_cohesive | not attempted | - | Pending. |
| 13 LVI_D7136 | not attempted | - | Pending. |
| 14 CAI_D7137 | not attempted | - | Pending. |
| 15 flat_ballistic | not attempted | - | Pending. |
| 16 PW_panel_ballistic | not attempted | - | Pending; Kok geometry port and runtime budget required. |

## Open Issues

Stage 02 no longer fails on missing MUMPS. The active blocker is robust completion of the `alpha=3` and `alpha=5` large-deflection implicit paths. See `tests/stage_02_cantilever_large/blocker.md` for attempted solver controls and failure modes. Stage 03 was independent enough to proceed despite that INCONCLUSIVE gate.

Stage 07 is blocked by the same LAW25 + TYPE14 compatibility issue for the D3039 on-axis and off-axis coupons. See `tests/stage_07_UD_tow_D3039/blocker.md`.

Stage 08 is likewise blocked for all seven required ply angles; the analytic stiffness transformation is reported but the FEM gate is not evaluated. See `tests/stage_08_ply_rotation/blocker.md`.

Stage 09 is blocked by the same per-ply LAW25 + TYPE14 incompatibility for both laminate stacks. See `tests/stage_09_laminate_solid_CLT/blocker.md`.
