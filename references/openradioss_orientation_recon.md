# OpenRadioss LAW12/LAW14 and TYPE6 Orientation Recon

Author: J.C. Vaught

Date: 2026-04-29

This note records the local OpenRadioss input grammar and qa-test examples used for the rotation-convention probe. Sources are the OpenRadioss CFG files and qa-tests in `/mnt/storage/j-vaught/openradioss/OpenRadioss-src/`.

## Source Files Read

The brief's unversioned `CFG/RADIOSS/...` path maps in this build to the versioned CFG hierarchy under `hm_cfg_files/config/CFG/radioss20xx/`. The 2026 CFG directory only contains deltas for these cards, so the latest full material/property definitions are:

- `hm_cfg_files/config/CFG/radioss2020/MAT/3d_comp_12.cfg`
- `hm_cfg_files/config/CFG/radioss2020/MAT/matl14_compso.cfg`
- `hm_cfg_files/config/CFG/radioss2025/PROP/prop_p6_sol_orth.cfg`
- `hm_cfg_files/config/CFG/radioss110/TABLE/inibri_ortho.cfg`
- `hm_cfg_files/config/CFG/radioss110/TABLE/inibri_ortho_sub.cfg`

## Documented Input Grammar

### `/MAT/LAW12` (`3D_COMP`)

CFG evidence: `3d_comp_12.cfg` defines the Radioss 2020 format at lines 293-350. It accepts `/MAT/LAW12/%d` or `/MAT/3D_COMP/%d`, then the following cards.

| Card | Fields | Category | Notes |
|---|---|---|---|
| Header/title | `/MAT/LAW12/{mat_ID}`, `TITLE` | identifier | Alias `/MAT/3D_COMP/{mat_ID}` is also exported. |
| Density | `MAT_RHO`, optional `Refer_Rho` | material constant | `Refer_Rho` is gated by the reference-density option. |
| Elastic moduli | `MAT_EA`, `MAT_EB`, `MAT_EC` | material constant | Young's modulus directions 1, 2, 3. |
| Poisson ratios | `MAT_PRAB`, `MAT_PRBC`, `MAT_PRCA` | material constant | 12, 23, 31 notation. |
| Shear moduli | `MAT_GAB`, `MAT_GBC`, `MAT_GCA` | material constant | 12, 23, 31 notation. |
| Tensile strengths/damage | `sigma_t1`, `sigma_t2`, `sigma_t3`, `delta` | strength/damage | Composite tensile strengths and tensile damage parameter. |
| Plastic hardening | `MAT_BETA`, `n`, `fmax`, `Wpref` | numerical/material flag | `BETA` is not an orientation angle here; it is the plasticity hardening parameter. |
| Yield strengths 1/2 | `sigma_1yt`, `sigma_2yt`, `sigma_1yc`, `sigma_2yc` | strength | Tension/compression yield strengths. |
| Shear yield strengths | `sigma_12yt`, `sigma_12yc`, `sigma_23yt`, `sigma_23yc` | strength | Shear traction/compression values. |
| Yield strengths 3/13 | `sigma_3yt`, `sigma_3yc`, `sigma_13yt`, `sigma_13yc` | strength | Present in LAW12. |
| Fiber/rate | `alpha`, `Ef`, `c`, `EPS_RATE_0`, `STRFLAG` | numerical/material flag | Fiber volume fraction, fiber modulus, strain-rate coefficient/reference, and rate flag. |

Orientation fields on LAW12: none. No `skew_ID`, `Phi`, `Vx/Vy/Vz`, `Px/Py/Pz`, `Iorth`, or `Ip` fields exist on the material card in the CFG.

### `/MAT/LAW14` (`COMPSO`)

CFG evidence: `matl14_compso.cfg` defines the Radioss 2020 format at lines 276-330. It accepts `/MAT/LAW14/%d` or `/MAT/COMPSO/%d`, then the following cards.

| Card | Fields | Category | Notes |
|---|---|---|---|
| Header/title | `/MAT/LAW14/{mat_ID}`, `TITLE` | identifier | Alias `/MAT/COMPSO/{mat_ID}` is also exported. |
| Density | `MAT_RHO`, optional `Refer_Rho` | material constant | `Refer_Rho` is gated by the reference-density option. |
| Elastic moduli | `MAT_EA`, `MAT_EB`, `MAT_EC` | material constant | Young's modulus directions 1, 2, 3. |
| Poisson ratios | `MAT_PRAB`, `MAT_PRBC`, `MAT_PRCA` | material constant | 12, 23, 31 notation. |
| Shear moduli | `MAT_GAB`, `MAT_GBC`, `MAT_GCA` | material constant | 12, 23, 31 notation. |
| Tensile strengths/damage | `SIGMA_T1`, `SIGMA_T2`, `SIGMA_T3`, `DELTA` | strength/damage | Composite tensile strengths and damage parameter. |
| Plastic hardening | `B`, `n`, `fmax`, `Wpref` | numerical/material flag | `B` is hardening, not orientation. |
| Yield strengths 1/2 | `sigma_1yt`, `sigma_2yt`, `sigma_1yc`, `sigma_2yc` | strength | Tension/compression yield strengths. |
| Shear yield strengths | `sigma_12yt`, `sigma_12yc`, `sigma_23yt`, `sigma_23yc` | strength | Shear traction/compression values. |
| Fiber/rate | `ALPHA`, `E_f`, `c`, `EPS_RATE_0`, `ICC` | numerical/material flag | `ICC` maps to `STRFLAG`. |

Orientation fields on LAW14: none. No `skew_ID`, `Phi`, `Vx/Vy/Vz`, `Px/Py/Pz`, `Iorth`, or `Ip` fields exist on the material card in the CFG.

### `/PROP/TYPE6` / `/PROP/SOL_ORTH`

CFG evidence: `prop_p6_sol_orth.cfg` defines the Radioss 2025 format at lines 348-413. It accepts `/PROP/TYPE6/%d` or `/PROP/SOL_ORTH/%d`, then the following cards.

| Card | Fields | Category | Notes |
|---|---|---|---|
| Header/title | `/PROP/TYPE6/{prop_ID}`, `TITLE` | identifier | Alias `/PROP/SOL_ORTH/{prop_ID}` is also valid. |
| Solid formulation | `Isolid`, `Ismstr`, `Icpre`, `Itetra10`, `Inpts`, `Itetra4`, `Iframe`, `Dn` | numerical flags | `Inpts` encodes R/S/T integration counts as hundreds/tens/ones when displayed. |
| Viscosity/hourglass | `qa`, `qb`, `h` | numerical flags | Bulk and hourglass viscosity. |
| Orientation vector and flags | `Vx`, `Vy`, `Vz`, `skew_ID`, `Ip`, `Iorth` | orientation | `Vx/Vy/Vz` is the reference vector; `skew_ID` references a skew frame; `Ip` is the reference plane/mode; `Iorth` is the orthotropic system formulation. |
| Orientation angle/point | `Phi`, `Px`, `Py`, `Pz` | orientation | `Phi` is the orthotropic angle with the first reference-plane direction. `Px/Py/Pz` exist in the 2022+ grammar and are exposed when `Ip > 20` in the GUI. |
| Time/distortion | `deltaT_min`, `vdef_min`, `vdef_max`, `ASP_max`, `COL_min` | numerical flags | Radioss 2025 format. Earlier 2022 format used `deltaT_min`, `Istrain`, `Ihkt` here. |
| Sol2SPH/distortion | `Ndir`, `sphpartID`, `Icontrol` | numerical flags | Radioss 2025 format. |

`Ip` options documented in the GUI block:

| `Ip` | Meaning |
|---:|---|
| 0 | Use `skew_ID`. |
| 1 | Plane `(r,s)` plus angle `Phi`. |
| 2 | Plane `(s,t)` plus angle `Phi`. |
| 3 | Plane `(t,r)` plus angle `Phi`. |
| 11 | Plane `(r,s)` plus orthogonal projection of `Vx,Vy,Vz` on plane `(r,s)`. |
| 12 | Plane `(s,t)` plus orthogonal projection of `Vx,Vy,Vz` on plane `(s,t)`. |
| 13 | Plane `(t,r)` plus orthogonal projection of `Vx,Vy,Vz` on plane `(t,r)`. |

`Iorth` options documented in the GUI block:

| `Iorth` | Meaning |
|---:|---|
| 0 | First orthotropy axis maintained at constant angle with respect to the orthonormal co-rotational element coordinate system. |
| 1 | First orthotropy direction is constant with respect to non-orthonormal isoparametric coordinates. |

Orientation is therefore controlled by `/PROP/TYPE6` or initial-state cards, not by `/MAT/LAW12` or `/MAT/LAW14`.

Starter source evidence adds modes that are not fully exposed by the GUI enum. In
`starter/source/properties/solid/hm_read_prop06.F`, the property reader accepts
`Ip = 20, 21, 23, 24` and rejects `Ip = 22` and values above 24. In
`starter/source/elements/solid/solide/smorth3.F`, those modes are implemented as:

| `Ip` | Source implementation |
|---:|---|
| 20 | Orthotropic directions from element connectivity, using brick edge 1-2 for axis 1 and edge 1-4 for axis 2. |
| 21 | Axis 1 from point `P` toward the element centroid; axis 2 is the in-plane perpendicular. |
| 23 | `Vx,Vy,Vz` plus `Phi`; starter raises `ERROR ID 1920` if the vector degenerates against the source frame. |
| 24 | `Vx,Vy,Vz` plus point `P`; source builds a frame from the vector and the element-centroid-to-point direction. |

The same source file converts `Ip = 0` plus a valid `skew_ID` into a negative
internal skew pointer before orthotropic directions are initialized.

### `/INIBRI/ORTHO`

CFG evidence: `inibri_ortho.cfg` defines `/INIBRI/ORTHO`; `inibri_ortho_sub.cfg` defines the subobject card at lines 108-127.

| Card | Fields | Category | Notes |
|---|---|---|---|
| Header | `/INIBRI/ORTHO` | identifier | Initial brick orthotropy state. |
| Element card | `brick_ID`, `Nb_layer`, `Isolnod`, `Prop_type`, `Isolid` | element/state metadata | For `/PROP/TYPE6`, `Prop_type = 6`. |
| Direction 1 | `X1`, `Y1`, `Z1` | orientation | First axis of orthotropy reference system in global coordinates. |
| Direction 2 | `X2`, `Y2`, `Z3` | orientation | Second axis in global coordinates. The CFG label says `Z3`, but GUI text calls it `Z2`; engine output comments write `X2,Y2` then `Z2`. |

This card is a plausible per-element orientation mechanism and is included in the empirical probe.

## qa-test Reconnaissance

Search scope: `/mnt/storage/j-vaught/openradioss/OpenRadioss-src/qa-tests/`.

### 1. LAW12 orthotropic direction check

Source: `qa-tests/miniqa/LOIS/LOI12/test1/data/readme.txt` describes the deck as `/MAT/LAW12 material orthotropic direction check`. The starter deck uses `LAW12` at lines 61-82 and three `/PROP/TYPE6` cards at lines 202-239.

Observed orientation syntax:

```radioss
/PROP/TYPE6/2
#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth
                   0                   0                   0         0         1         0
#                phi
                   0

/PROP/TYPE6/4
#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth
                   0                   0                   0         0         2         0
#                phi
                   0

/PROP/TYPE6/6
#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth
                   0                   0                   0         0         3         0
#                phi
                   0
```

This is not an off-axis `Phi` test; it varies `Ip` while the three brick connectivities differ.

### 2. LAW14 different material orientations

Source: `qa-tests/miniqa/LOIS/LOI14/data/readme.txt` describes the deck as `Test /MAT/LAW14 with different material orientations`. The starter deck uses `LAW14` at lines 59-78 and three `/PROP/TYPE6` cards at lines 192-229.

Observed orientation syntax:

```radioss
/PROP/TYPE6/2
#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth
                   0                   0                   0         0         1         0
#                phi
                   0

/PROP/TYPE6/4
#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth
                   0                   0                   0         0         2         0
#                phi
                   0

/PROP/TYPE6/6
#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth
                   0                   0                   0         0         3         0
#                phi
                   0
```

Like the LAW12 qa test, this exercises `Ip` and element connectivity, not a nonzero `Phi`.

### 3. LAW12 + SOL_ORTH property `Phi = 90`

Source: `qa-tests/miniqa/PROP/PROPSOL_ORTH_Istrain/data/Test_traction_0000.rad`. The deck uses `LAW12` at lines 28-39 and `/PROP/SOL_ORTH/1` at lines 77-83.

Observed orientation syntax:

```radioss
/PROP/SOL_ORTH/1
Fil
        24         0                   0
                 1.1                0.05                 0.1
                 1.0                 0.0                 0.0         0         0
                90.0
           1000000.0         1
```

This is the only qa-test hit found with `LAW12`/`LAW14` and a visibly nonzero TYPE6/SOL_ORTH `Phi` value. It is old-format SOL_ORTH syntax, so the line after the vector/flags card contains only `Phi`, not `Phi Px Py Pz`.

## Immediate Probe Implications

- `LAW12` and `LAW14` do not own fiber-frame rotation slots; all orientation mechanisms to test are property-level or initial-state-level.
- The canonical current `/PROP/TYPE6` card supports exactly the requested `Vx/Vy/Vz`, `skew_ID`, `Ip`, `Iorth`, `Phi`, and `Px/Py/Pz` slots.
- qa-tests provide confidence that `Ip = 1, 2, 3` and old-format `Phi = 90` are accepted, but they do not prove the analytic off-axis transformed modulus for a uniaxial coupon.
- `/INIBRI/ORTHO` is present in CFG grammar and supports per-brick axes for `Prop_type = 6`; it must be checked empirically because no qa-test deck using `/INIBRI/ORTHO` was found.
- Source-only `Ip = 20, 21, 23, 24` modes are included in the empirical sweep as M6/source-driven candidates.
