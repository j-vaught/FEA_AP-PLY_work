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
| 07 UD_tow_D3039 | PASS | 9.6 s | Post-rotation LAW12 + TYPE6/SOL_ORTH decks use the verified `Ip=3`, `Iorth=0`, `Phi=theta` property convention. 0 deg error 0.403%; 45 deg off-axis error 0.114% using engineering strain from displacement. |
| 08 ply_rotation | PASS | 14.0 s | Seven-angle LAW12 + TYPE6/SOL_ORTH sweep uses `Ip=3`, `Iorth=0`, `Phi=theta` on a compact 4x4x1 HEXA8 block. Max analytic Ex(theta) error 0.815%; VTK cell `Stra[0]` is not used for the gate. |
| 09 laminate_solid_CLT | PASS | 35.1 s | Cross-ply `[0/90]s` and quasi-isotropic `[0/+45/-45/90]s` solid-ply stacks use per-ply LAW12 + TYPE6 with `Ip=3`, `Iorth=0`, `Phi=theta_ply`. Six starter+engine+VTK jobs passed the CLT A-matrix gate; max component error 1.085%. |
| 10 UD_mesoscale_direct | FAIL | 6.3 s | Self-contained isotropic fiber/matrix HEXA8 voxel cell starts, runs, converts to VTK, and post-processes. Axial ROM gate passes at 0.026% error, but transverse Halpin-Tsai error is 11.189% and longitudinal shear error is 99.934%; see `tests/stage_10_UD_mesoscale_direct/blocker.md`. |
| 11 PW_mesoscale_direct | not attempted | - | Pending; Kok geometry port required after stage 10. |
| 12 DCB_ENF_cohesive | not attempted | - | Pending. |
| 13 LVI_D7136 | not attempted | - | Pending. |
| 14 CAI_D7137 | not attempted | - | Pending. |
| 15 flat_ballistic | not attempted | - | Pending. |
| 16 PW_panel_ballistic | not attempted | - | Pending; Kok geometry port and runtime budget required. |

## Open Issues

Stage 02 no longer fails on missing MUMPS. The active blocker is robust completion of the `alpha=3` and `alpha=5` large-deflection implicit paths. See `tests/stage_02_cantilever_large/blocker.md` for attempted solver controls and failure modes. Stage 03 was independent enough to proceed despite that INCONCLUSIVE gate.

Stage 10 now runs end-to-end on this host rather than relying on the stale Lima/VTKHDF wrapper, but the direct mesoscale stiffness gate fails numerically. The current blocker is the transverse/shear mismatch of the isotropic fiber/matrix voxel cell against Halpin-Tsai targets. Per workflow, stages 11-16 were not advanced after this genuine stage-10 gate miss.
