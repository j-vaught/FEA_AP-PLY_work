# OpenRadioss LAW12/LAW14 TYPE6 Orientation Convention

Author: J.C. Vaught

Date: 2026-04-29

## Verdict

The working convention for the current OpenRadioss build is property/initial-state
orientation on `/PROP/TYPE6`, not material-card orientation on `/MAT/LAW12` or
`/MAT/LAW14`.

Recommended choices:

1. Uniform off-axis coupon: `/PROP/TYPE6` with `Ip = 3`, `Iorth = 0`,
   `Phi = theta`. This matched all five angles for both LAW12 and LAW14 with
   maximum relative error 0.24%.
2. Per-element fiber orientation: `/INIBRI/ORTHO`, one element card per brick,
   with global material axis 1 set to `(cos(theta_e), sin(theta_e), 0)` and axis
   2 set to `(-sin(theta_e), cos(theta_e), 0)`. A neutral `/PROP/TYPE6` with
   `Ip = 1`, `Iorth = 0`, `Phi = 0` is sufficient. This also matched all five
   angles with maximum relative error 0.24%.
3. Equivalent property vector option: `/PROP/TYPE6` with
   `Vx,Vy,Vz = (cos(theta), sin(theta), 0)`, `Ip = 13`, `Iorth = 0`. This is
   useful when a vector is preferable to `Phi`, but it is still property-level.

Do not use the previous stage-07 convention `skew_ID = 1`, `Ip = 1`,
`Iorth = 1`, `Phi = 0` for off-axis solids. In the probe it held the response
near the 0-degree fiber modulus for off-axis angles because `Ip = 1` does not
consume the skew frame.

## Primary Sources

- `hm_cfg_files/config/CFG/radioss2020/MAT/3d_comp_12.cfg`: LAW12 grammar; no
  orientation fields.
- `hm_cfg_files/config/CFG/radioss2020/MAT/matl14_compso.cfg`: LAW14 grammar; no
  orientation fields.
- `hm_cfg_files/config/CFG/radioss2025/PROP/prop_p6_sol_orth.cfg`: TYPE6/SOL_ORTH
  grammar; owns `Vx,Vy,Vz`, `skew_ID`, `Ip`, `Iorth`, `Phi`, `Px,Py,Pz`.
- `hm_cfg_files/config/CFG/radioss110/TABLE/inibri_ortho*.cfg`: `/INIBRI/ORTHO`
  grammar for per-brick axes.
- `starter/source/properties/solid/hm_read_prop06.F` and
  `starter/source/elements/solid/solide/smorth3.F`: starter implementation of
  `Ip = 0, 1, 2, 3, 11, 12, 13, 20, 21, 23, 24`.
- qa-test examples cited in `references/openradioss_orientation_recon.md`.

## Probe Method

The driver is `references/orientation_probe_decks/run_orientation_probes.py`.
It generated 4 x 4 x 1 HEXA8 blocks, ran `starter -check`, then starter, engine,
and `anim_to_vtk`, and compared apparent axial modulus to:

```text
1/E_x(theta) = c^4/E1 + s^4/E2 + (1/G12 - 2*nu12/E1) c^2 s^2
```

The exact analytic values used from `references/material_cards/im7_8552.json`
are:

| theta | analytic E_x, GPa |
|---:|---:|
| 0 | 171.400 |
| 30 | 22.267 |
| 45 | 13.277 |
| 60 | 10.303 |
| 90 | 9.080 |

These 30-degree and 60-degree values are the exact result of the stated formula
and JSON constants; they differ slightly from the approximate values in the
rotation brief.

Primary comparison metric:

```text
E_global = mean(VTK cell sigma_xx over middle cells) /
           engineering strain from VTK displacement difference
```

Important measurement caveat: for rotated TYPE6 solids, the VTK cell tensor
strain component `Stra[0]` is not a reliable global axial engineering strain.
Using `mean(sigma_xx) / mean(Stra[0])` recreates the apparent 45-degree failure.
For example, LAW12 + `Ip = 3`, `Iorth = 0`, `Phi = 45` gives
`E_global = 13.271 GPa` but the cell-tensor ratio gives `35.326 GPa`.

Full row-level results are saved in:

- `references/orientation_probe_decks/orientation_probe_results.csv`
- `references/orientation_probe_decks/orientation_probe_results.json`

## Recommended Snippets

### Uniform property angle

Use this when every brick in the property has the same fiber angle. For
theta = 45 degrees:

```radioss
/MAT/LAW12/1
IM7_8552_LAW12_orientation_probe
#              RHO_I
                1580
#                E11                 E22                 E33
        171400000000          9080000000          9080000000
#               NU12                NU23                NU31
                0.32                 0.5     0.0169521586931
#                G12                 G23                 G31
          5290000000       3026666666.67          5290000000
#           SIGMA_T1            SIGMA_T2            SIGMA_T3               DELTA
          2560000000            73000000            73000000                0.05
#                  B                   n                fmax               Wpref
                   1                   1                   1                   1
#          sigma_1yt           sigma_2yt           sigma_1yc           sigma_2yc
          2560000000            73000000          1590000000           185000000
#         sigma_12yt          sigma_12yc          sigma_23yt          sigma_23yc
            90000000            90000000            90000000            90000000
#          sigma_3yt           sigma_3yc          sigma_13yt          sigma_13yc
            73000000           185000000            90000000            90000000
#              alpha                  Ef                   c          EPS_RATE_0   STRFLAG
                   0                   0                   0                   0         1
/PROP/TYPE6/1
type6_sol_orth_phi_ip3
#   Isolid    Ismstr               Icpre  Itetra10     Inpts   Itetra4    Iframe                  Dn
        24         4                   1         0         0         0         2                   0
#                 qa                  qb                   h
                   0                   0                   0
#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth
                   1                   0                   0         0         3         0
#                Phi                 Px                  Py                  Pz
                  45                   0                   0                   0
#             dt_min   istrain      IHKT
                   0         0         0
```

Change only `Phi` for other uniform angles. No `/SKEW` block is needed; keep
`skew_ID = 0`.

### Per-element `/INIBRI/ORTHO`

Use this when each brick needs its own fiber direction. The property may remain
neutral; `/INIBRI/ORTHO` supplies the actual axes. For theta = 45 degrees and
brick 1:

```radioss
/PROP/TYPE6/1
type6_sol_orth_inibri
#   Isolid    Ismstr               Icpre  Itetra10     Inpts   Itetra4    Iframe                  Dn
        24         4                   1         0         0         0         2                   0
#                 qa                  qb                   h
                   0                   0                   0
#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth
                   1                   0                   0         0         1         0
#                Phi                 Px                  Py                  Pz
                   0                   0                   0                   0
#             dt_min   istrain      IHKT
                   0         0         0
/INIBRI/ORTHO
#  brick_ID  Nb_layer  Isolnod Prop_type    Isolid
#                  X1                  Y1                  Z1                  X2                  Y2
#                  Z2
         1         1         8         6        24
      0.707106781187      0.707106781187                   0     -0.707106781187      0.707106781187
                   0
```

Repeat the `/INIBRI/ORTHO` element card for each brick, replacing axis 1 and
axis 2 with the element's own `(cos(theta_e), sin(theta_e), 0)` and
`(-sin(theta_e), cos(theta_e), 0)`.

## Recipe Evidence

LAW12, `Ip = 3`, `Iorth = 0`, `Phi = theta`:

| theta | E_global GPa | analytic GPa | rel err % | cell tensor ratio GPa |
|---:|---:|---:|---:|---:|
| 0 | 171.239 | 171.400 | 0.09 | 171.398 |
| 30 | 22.213 | 22.267 | 0.24 | 104.892 |
| 45 | 13.271 | 13.277 | 0.04 | 35.326 |
| 60 | 10.301 | 10.303 | 0.02 | 16.149 |
| 90 | 9.065 | 9.080 | 0.16 | 9.081 |

LAW12, `/INIBRI/ORTHO`, `Ip = 1`, `Iorth = 0`:

| theta | E_global GPa | analytic GPa | rel err % | cell tensor ratio GPa |
|---:|---:|---:|---:|---:|
| 0 | 171.239 | 171.400 | 0.09 | 171.398 |
| 30 | 22.213 | 22.267 | 0.24 | 104.892 |
| 45 | 13.271 | 13.277 | 0.04 | 35.326 |
| 60 | 10.301 | 10.303 | 0.02 | 16.149 |
| 90 | 9.065 | 9.080 | 0.16 | 9.081 |

LAW14 produced the same values for these recipes.

## Comparison Table

The table below is the full combination-level result matrix. The CSV/JSON files
contain every angle row with starter and engine return codes.

| LAW | Mechanism | Iorth | Ip | Verdict | angle passes | max err % |
|---|---|---:|---:|---|---:|---:|
| LAW12 | `M1_phi_ip1` | 0 | 1 | FAIL | 1/5 | 94.71 |
| LAW12 | `M1_phi_ip1` | 1 | 1 | FAIL | 2/5 | 336.43 |
| LAW12 | `M1_phi_ip3` | 0 | 3 | PASS | 5/5 | 0.24 |
| LAW12 | `M1_phi_ip3` | 1 | 3 | FAIL | 1/5 | 1786.57 |
| LAW12 | `M2_skew_ip0` | 0 | 0 | FAIL | 1/5 | 94.67 |
| LAW12 | `M2_skew_ip0` | 1 | 0 | FAIL | 1/5 | 94.74 |
| LAW12 | `M2_skew_requested_ip1` | 0 | 1 | FAIL | 1/5 | 94.71 |
| LAW12 | `M2_skew_requested_ip1` | 1 | 1 | FAIL | 1/5 | 1786.57 |
| LAW12 | `M2_skew_requested_ip3` | 0 | 3 | FAIL | 1/5 | 1785.90 |
| LAW12 | `M2_skew_requested_ip3` | 1 | 3 | FAIL | 1/5 | 94.73 |
| LAW12 | `M3_vector_ip11` | 0 | 11 | FAIL | 1/5 | 94.71 |
| LAW12 | `M3_vector_ip11` | 1 | 11 | FAIL | 1/5 | 1786.57 |
| LAW12 | `M3_vector_ip13` | 0 | 13 | PASS | 5/5 | 0.24 |
| LAW12 | `M3_vector_ip13` | 1 | 13 | FAIL | 1/5 | 1786.57 |
| LAW12 | `M3_vector_requested_ip1` | 0 | 1 | FAIL | 1/5 | 94.71 |
| LAW12 | `M3_vector_requested_ip1` | 1 | 1 | FAIL | 1/5 | 1786.57 |
| LAW12 | `M3_vector_requested_ip3` | 0 | 3 | FAIL | 1/5 | 1785.90 |
| LAW12 | `M3_vector_requested_ip3` | 1 | 3 | FAIL | 1/5 | 94.73 |
| LAW12 | `M4_vector_point_ip24` | 0 | 24 | PASS | 5/5 | 1.59 |
| LAW12 | `M4_vector_point_ip24` | 1 | 24 | FAIL | 1/5 | 1786.53 |
| LAW12 | `M5_inibri_ortho_ip1` | 0 | 1 | PASS | 5/5 | 0.24 |
| LAW12 | `M5_inibri_ortho_ip1` | 1 | 1 | PASS | 5/5 | 0.76 |
| LAW12 | `M5_inibri_ortho_ip20` | 0 | 20 | PASS | 5/5 | 0.24 |
| LAW12 | `M5_inibri_ortho_ip20` | 1 | 20 | PASS | 5/5 | 0.76 |
| LAW12 | `M6_point_ip21` | 0 | 21 | PASS | 4/5 | 2.40 |
| LAW12 | `M6_point_ip21` | 1 | 21 | FAIL | 1/5 | 1788.62 |
| LAW12 | `M6_vphi_ip23` | 0 | 23 | BLOCKED | 0/5 | blocked |
| LAW12 | `M6_vphi_ip23` | 1 | 23 | BLOCKED | 0/5 | blocked |
| LAW14 | `M1_phi_ip1` | 0 | 1 | FAIL | 1/5 | 94.71 |
| LAW14 | `M1_phi_ip1` | 1 | 1 | FAIL | 2/5 | 336.43 |
| LAW14 | `M1_phi_ip3` | 0 | 3 | PASS | 5/5 | 0.24 |
| LAW14 | `M1_phi_ip3` | 1 | 3 | FAIL | 1/5 | 1786.57 |
| LAW14 | `M2_skew_ip0` | 0 | 0 | FAIL | 1/5 | 94.67 |
| LAW14 | `M2_skew_ip0` | 1 | 0 | FAIL | 1/5 | 94.74 |
| LAW14 | `M2_skew_requested_ip1` | 0 | 1 | FAIL | 1/5 | 94.71 |
| LAW14 | `M2_skew_requested_ip1` | 1 | 1 | FAIL | 1/5 | 1786.57 |
| LAW14 | `M2_skew_requested_ip3` | 0 | 3 | FAIL | 1/5 | 1785.90 |
| LAW14 | `M2_skew_requested_ip3` | 1 | 3 | FAIL | 1/5 | 94.73 |
| LAW14 | `M3_vector_ip11` | 0 | 11 | FAIL | 1/5 | 94.71 |
| LAW14 | `M3_vector_ip11` | 1 | 11 | FAIL | 1/5 | 1786.57 |
| LAW14 | `M3_vector_ip13` | 0 | 13 | PASS | 5/5 | 0.24 |
| LAW14 | `M3_vector_ip13` | 1 | 13 | FAIL | 1/5 | 1786.57 |
| LAW14 | `M3_vector_requested_ip1` | 0 | 1 | FAIL | 1/5 | 94.71 |
| LAW14 | `M3_vector_requested_ip1` | 1 | 1 | FAIL | 1/5 | 1786.57 |
| LAW14 | `M3_vector_requested_ip3` | 0 | 3 | FAIL | 1/5 | 1785.90 |
| LAW14 | `M3_vector_requested_ip3` | 1 | 3 | FAIL | 1/5 | 94.73 |
| LAW14 | `M4_vector_point_ip24` | 0 | 24 | PASS | 5/5 | 1.59 |
| LAW14 | `M4_vector_point_ip24` | 1 | 24 | FAIL | 1/5 | 1786.53 |
| LAW14 | `M5_inibri_ortho_ip1` | 0 | 1 | PASS | 5/5 | 0.24 |
| LAW14 | `M5_inibri_ortho_ip1` | 1 | 1 | PASS | 5/5 | 0.76 |
| LAW14 | `M5_inibri_ortho_ip20` | 0 | 20 | PASS | 5/5 | 0.24 |
| LAW14 | `M5_inibri_ortho_ip20` | 1 | 20 | PASS | 5/5 | 0.76 |
| LAW14 | `M6_point_ip21` | 0 | 21 | PASS | 4/5 | 2.40 |
| LAW14 | `M6_point_ip21` | 1 | 21 | FAIL | 1/5 | 1788.62 |
| LAW14 | `M6_vphi_ip23` | 0 | 23 | BLOCKED | 0/5 | blocked |
| LAW14 | `M6_vphi_ip23` | 1 | 23 | BLOCKED | 0/5 | blocked |

## Mechanism Labels

- `M1_phi_ip1`, `M1_phi_ip3`: `Phi` field on `/PROP/TYPE6`.
- `M2_skew_requested_*`: `skew_ID` populated while `Ip` remains 1 or 3.
- `M2_skew_ip0`: documented `skew_ID` mode with `Ip = 0`.
- `M3_vector_requested_*`: vector populated while `Ip` remains 1 or 3.
- `M3_vector_ip11`, `M3_vector_ip13`: documented vector projection modes.
- `M4_vector_point_ip24`: source-revealed `V + P` mode.
- `M5_inibri_ortho_*`: `/INIBRI/ORTHO` per-brick axes.
- `M6_point_ip21`, `M6_vphi_ip23`: source-revealed point/vector modes.
