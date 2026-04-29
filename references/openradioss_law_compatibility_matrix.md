# OpenRadioss LAW x /PROP x /FAIL compatibility matrix

This matrix is empirical for the OpenRadioss build on `comech-2422`. It supersedes the earlier audit assumption that `/MAT/LAW25` or other solid composite laws are acceptable on `/PROP/TYPE14` in this build.

The decisive result is:

- `/PROP/TYPE14` is blocked for the relevant composite/orthotropic solid laws: `LAW12`, `LAW14`, `LAW25`, `LAW28`, `LAW53`, and `LAW128` all fail with `ERROR ID : 3047` or `3046`.
- `/PROP/TYPE6` (`/PROP/SOL_ORTH`) is the compatible solid-brick property for the practical composite candidates: `LAW12`, `LAW14`, `LAW25`, `LAW28`, `LAW53`, and `LAW128` all parse on the same HEXA8 probe.
- `/FAIL/HASHIN`, `/FAIL/PUCK`, and `/FAIL/TSAIWU` parse with the `TYPE6` solid-brick composite rows above. They do not rescue `TYPE14`; `TYPE14` remains blocked before the failure model matters.

## Verified On

- OpenRadioss starter: `/mnt/storage/j-vaught/openradioss/OpenRadioss/exec/starter_linux64_gf`
- Starter version:
  - `OpenRadioss Starter`
  - `Platform release : linux64`
  - `Platform info : Linux 64 bits, GNU compiler`
  - `Time of build : 14:07:11`
  - `Date of build : Apr 29 2026`
  - `CommitID : User Build`
  - `Reader : 20260429_8d9da7aa`
- MUMPS-linked: yes, using the `$OR/exec` binaries from the requested OpenRadioss tree.
- Date: 2026-04-29
- Environment:
  - `OR=/mnt/storage/j-vaught/openradioss/OpenRadioss`
  - `RAD_CFG_PATH=$OR/hm_cfg_files`
  - `RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64`
  - `LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH`

## Probe Method

Probe decks are in `references/lawmatrix_probe_decks/`. Each run used:

- one HEXA8 solid brick for `TYPE14` and `TYPE6`;
- one 4-node shell for `TYPE1`;
- one material id (`/MAT/LAW*`);
- optional `/FAIL/*` card on the same material id;
- a minimal `/BCS` clamp and `/CLOAD`;
- `starter_linux64_gf -i deck_*_0000.rad -check`.

OpenRadioss requires starter input files to end in `_0000.rad`; the probe decks keep that suffix. Full raw results are in `references/lawmatrix_probe_decks/results.csv`; the corresponding `*.out` listings are retained beside each deck.

`OK` means the starter returned zero errors. Some accepted 2026-format minimal decks carry warning IDs `100214` and/or `100217`; these are retained in `results.csv` and the listings. They did not change the material/property/failure compatibility verdicts.

The brief's `/FAIL/JOHNSON_COOK` request maps to the native OpenRadioss keyword `/FAIL/JOHNSON`. `/FAIL/USER1` was tested literally as `/FAIL/USER1/1`.

## Candidate Enumeration

Candidates were selected from the native Radioss CFG tree under `/mnt/storage/j-vaught/openradioss/OpenRadioss-src/hm_cfg_files/config/CFG/`, searching `/MAT/LAW*` definitions and headers for composite, orthotropic, anisotropic, transversely isotropic, Tsai, Hashin, Puck, ply, fiber, and fabric terms.

Included native LAW candidates:

- `LAW1` (`ELAST`) and `LAW2` (`PLAS_JOHNS`) as baselines.
- `LAW12` (`3D_COMP`), `LAW14` (`COMPSO`), `LAW25` (`COMPSH`), `LAW28` (`HONEYCOMB`), `LAW53` (`TSAI_TAB`), `LAW58` (`FABR_A`), `LAW128` (`HILL_VISC_PLAST`) as practical anisotropic/composite-adjacent candidates.
- `LAW42` (`OGDEN`) and `LAW62` (`VISC_HYP`) as non-composite controls.
- `LAW123`, `LAW125`, `LAW127`, `LAW132`, and `LAW158` because current CFGs expose laminated/fabric material names matching the search terms.

Excluded after CFG search:

- `LAW17`: no native `/MAT/LAW17` CFG was found in this build's Radioss native material hierarchy.
- `LAW55`: no native `/MAT/LAW55` CFG was found. LS-DYNA `*MAT_054/055` compatibility files exist under Keyword971 trees, but those are not native OpenRadioss `/MAT/LAW55` cards.

## Matrix

Legend:

- `B####` means `BLOCKED:<error_id>`.
- In `/FAIL` columns, the pair is `TYPE14/TYPE6`.
- `TYPE6` is `/PROP/TYPE6` / `/PROP/SOL_ORTH`, an orthotropic solid property, not a shell property.
- `TYPE1` is a plain shell property probe for contrast.
- Non-composite baselines may syntactically accept composite `/FAIL` cards; that is a parser result, not a physical recommendation.

| LAW | Name | TYPE14 solid | TYPE6/SOL_ORTH solid | TYPE1 shell | /FAIL/HASHIN TYPE14/TYPE6 | /FAIL/PUCK TYPE14/TYPE6 | /FAIL/TSAIWU TYPE14/TYPE6 | Other /FAIL TYPE14/TYPE6 | Notes |
|---|---|---:|---:|---:|---:|---:|---:|---|---|
| LAW1 | ELAST | OK | OK | OK | OK/OK | OK/OK | OK/OK | CHANG OK/OK; JOHNSON OK/OK; USER1 B1130/B1130 | Isotropic baseline; composite failures are parser-only n/a physically. |
| LAW2 | PLAS_JOHNS | OK | OK | OK | OK/OK | OK/OK | OK/OK | CHANG OK/OK; JOHNSON OK/OK; USER1 B1130/B1130 | Johnson-Cook plastic baseline; composite failures are parser-only n/a physically. |
| LAW12 | 3D_COMP | B3047 | OK | B3046 | B3047/OK | B3047/OK | B3047/OK | CHANG B3047/OK; JOHNSON B3047/OK; USER1 B1130/B1130 | Best verified solid composite row. Use `TYPE6`, not `TYPE14`. |
| LAW14 | COMPSO | B3047 | OK | B3046 | B3047/OK | B3047/OK | B3047/OK | CHANG B3047/OK; JOHNSON B3047/OK; USER1 B1130/B1130 | Legacy composite solid orthotropic row. Use `TYPE6`, not `TYPE14`. |
| LAW15 | CHANG | B3046 | B3046 | B1052 | B3046/B3046 | B3046/B3046 | B3046/B3046 | CHANG B3046/B3046; JOHNSON B3046/B3046; USER1 B1130/B1130 | Did not parse with these solid or plain shell property probes. |
| LAW25 | COMPSH | B3047 | OK | B1052 | B3047/OK | B3047/OK | B3047/OK | CHANG B3047/OK; JOHNSON B3047/OK; USER1 B1130/B1130 | `TYPE14` blocker reproduced. `TYPE6` solid proxy is accepted. |
| LAW28 | HONEYCOMB | B3047 | OK | B3046 | B3047/OK | B3047/OK | B3047/OK | CHANG B3047/OK; JOHNSON B3047/OK; USER1 B1130/B1130 | Orthotropic honeycomb, not a fiber-ply composite recommendation. |
| LAW42 | OGDEN | OK | OK | OK | OK/OK | OK/OK | OK/OK | CHANG OK/OK; JOHNSON OK/OK; USER1 B1130/B1130 | Hyperelastic control; composite failures are parser-only n/a physically. |
| LAW53 | TSAI_TAB | B3047 | OK | B3046 | B3047/OK | B3047/OK | B3047/OK | CHANG B3047/OK; JOHNSON B3047/OK; USER1 B1130/B1130 | Tsai-tabulated anisotropic row; `TYPE6` only. |
| LAW58 | FABR_A | B3046 | B3046 | B3047 | B3046/B3046 | B3046/B3046 | B3046/B3046 | CHANG B3046/B3046; JOHNSON B3046/B3046; USER1 B1130/B1130 | Did not parse with these solid or plain shell property probes. |
| LAW62 | VISC_HYP | OK | OK | OK | OK/OK | OK/OK | OK/OK | CHANG OK/OK; JOHNSON OK/OK; USER1 B1130/B1130 | Visco-hyperelastic control; composite failures are parser-only n/a physically. |
| LAW123 | LAMINATED_FRACTURE_DAIMLER_PINHO | B3033 | B3033 | B3033 | B3033/B3033 | B3033/B3033 | B3033/B3033 | CHANG B3033/B3033; JOHNSON B3033/B3033; USER1 B3033/B3033 | Representative card blocked in material-law parsing before property compatibility is useful. |
| LAW125 | LAMINATED_COMPOSITE | B3068 | B3068 | B3068 | B3068/B3068 | B3068/B3068 | B3068/B3068 | CHANG B3068/B3068; JOHNSON B3068/B3068; USER1 B3068/B3068 | Representative card blocked in material-law parsing before property compatibility is useful. |
| LAW127 | ENHANCED_COMPOSITE | B3068 | B3068 | B3068 | B3068/B3068 | B3068/B3068 | B3068/B3068 | CHANG B3068/B3068; JOHNSON B3068/B3068; USER1 B3068/B3068 | Representative card blocked in material-law parsing before property compatibility is useful. |
| LAW128 | HILL_VISC_PLAST | B3047 | OK | B3047 | B3047/OK | B3047/OK | B3047/OK | CHANG B3047/OK; JOHNSON B3047/OK; USER1 B1130/B1130 | Anisotropic viscoplastic metal, not a fiber-ply composite recommendation. |
| LAW132 | LAMINATED_FRACTURE_DAIMLER_CAMANHO | B3033 | B3033 | B3033 | B3033/B3033 | B3033/B3033 | B3033/B3033 | CHANG B3033/B3033; JOHNSON B3033/B3033; USER1 B3033/B3033 | Representative card blocked in material-law parsing before property compatibility is useful. |
| LAW158 | FABR_NL | B1581 | B1581 | B1581 | B1581/B1581 | B1581/B1581 | B1581/B1581 | CHANG B1581/B1581; JOHNSON B1581/B1581; USER1 B1130/B1130 | Representative nonlinear fabric card blocked in material-law parsing. |

## Error ID Legend

- `3047`: material/property compatibility error. This is the decisive `TYPE14` blocker for `LAW12`, `LAW14`, `LAW25`, `LAW28`, `LAW53`, and `LAW128`.
- `3046`: material/element compatibility error.
- `1052`: property-set error; in the retained listings this is often accompanied by an additional material/property compatibility error.
- `3033`, `3068`, `1581`: material-law input errors for the representative cards used for those advanced laws.
- `1130`: user-interface error for literal `/FAIL/USER1/1`.
- Some blocked listings also include secondary `ERROR ID : 307`; the table records the first blocking error id emitted by the classifier.

## Recommended Substitutions

These recommendations cite only rows that are directly verified above.

| Downstream need | Verified substitution | Matrix evidence | Notes |
|---|---|---|---|
| Linear orthotropic solid, no damage | `/MAT/LAW12` + `/PROP/TYPE6` (`/PROP/SOL_ORTH`), no `/FAIL` | `LAW12` row: `TYPE6/SOL_ORTH solid = OK`; `TYPE14 solid = B3047` | Keeps HEXA8 solid elements. Do not keep the old `TYPE14` requirement. `LAW14 + TYPE6` is also verified, but `LAW12` is the preferred 3D composite row. |
| Solid composite with Hashin failure | `/MAT/LAW12` + `/PROP/TYPE6` + `/FAIL/HASHIN` | `LAW12` row: `/FAIL/HASHIN TYPE14/TYPE6 = B3047/OK` | Listing `deck_LAW12_TYPE6_HASHIN_0000.out` prints `FAILURE MODEL TYPE = HASHIN` and `IFAIL_SO = 1`. |
| Solid composite with Tsai-Wu failure | `/MAT/LAW12` + `/PROP/TYPE6` + `/FAIL/TSAIWU` | `LAW12` row: `/FAIL/TSAIWU TYPE14/TYPE6 = B3047/OK` | The separate `/FAIL/TSAIWU` card parses on the verified solid `TYPE6` path. |
| Solid composite with Puck failure | `/MAT/LAW12` + `/PROP/TYPE6` + `/FAIL/PUCK` | `LAW12` row: `/FAIL/PUCK TYPE14/TYPE6 = B3047/OK` | Listing `deck_LAW12_TYPE6_PUCK_0000.out` prints `FAILURE MODEL TYPE = PUCK` and `IFAIL_SO = 1`. |
| Solid composite with element erosion for ballistic | `/MAT/LAW12` + `/PROP/TYPE6` + `/FAIL/HASHIN` with `IFAIL_SH=0`, `PTHICKFAIL=1.0`, `IFAIL_SO=1` | `LAW12` row: `/FAIL/HASHIN TYPE14/TYPE6 = B3047/OK`; listing confirms solid deletion flag | This is only a starter compatibility recommendation, not a calibrated ballistic validation. It preserves the all-solid HEXA path through `TYPE6/SOL_ORTH`. |

## Practical Stop Rule For Future Decks

For this OpenRadioss build, do not spend time trying to make `LAW12`, `LAW14`, `LAW25`, `LAW28`, `LAW53`, or `LAW128` work with `/PROP/TYPE14` for composite solid bricks. The starter rejects those combinations. Move the solid composite decks to `/PROP/TYPE6` (`/PROP/SOL_ORTH`) or document a blocker if a downstream stage still requires `/PROP/TYPE14` specifically.
