# Codex build-and-resume brief — FEA_AP-PLY

You are running autonomously on `comech-2422` in tmux session `codex-fea`. The previous run stopped at stage 02 with `Fatal error: MUMPS required` — the prebuilt OpenRadioss linux64_gf binaries do not include MUMPS, so `/IMPL/NONLIN` cannot run. Your job is to **build a MUMPS-linked OpenRadioss from source, then resume stages 02-16**.

**Author identity for git commits:** `j-vaught <jvaught@sc.edu>`. Do not mention any AI provider.

## Required reading (in order)

1. `README.md`
2. `plan/master_plan.md`
3. `CODEX_BRIEF.md` — the original autonomous brief; everything in it still applies after the build.
4. `RESULTS_SUMMARY.md` — current state: stage 01 PASS, stage 02 INCONCLUSIVE.
5. `tests/stage_02_cantilever_large/blocker.md` — exact MUMPS error.
6. `/mnt/storage/j-vaught/openradioss/OpenRadioss-src/` — has the OpenRadioss source already sparse-checked-out for `qa-tests` only. You will need to expand the sparse-checkout to include the full source tree (`common_source`, `starter`, `engine`, `extlib`, `CMake_Compilers`, `Compiling_tools`, `hm_cfg_files`, `reader`, `scripts`, `tools`). Use `git sparse-checkout disable` to get the full tree.

## Compute environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate feaapply
# python 3.12, gmsh, meshio, pyvista, numpy, scipy, jinja2, openmpi, cmake, make
# system gcc 13.4, gfortran 11.4

export OR=/mnt/storage/j-vaught/openradioss/OpenRadioss
export OR_SRC=/mnt/storage/j-vaught/openradioss/OpenRadioss-src
export RAD_CFG_PATH=$OR/hm_cfg_files
export RAD_H3D_PATH=$OR/extlib/h3d/lib/linux64
export LD_LIBRARY_PATH=$OR/extlib/hm_reader/linux64:$LD_LIBRARY_PATH
```

User account `j-vaught`, no sudo. 502 GB RAM, 112 cores, 10.8 TB free at `/mnt/storage/j-vaught/`.

## Phase A — Build MUMPS-linked OpenRadioss

A1. Install MUMPS into the `feaapply` conda env:
```bash
mamba install -y -c conda-forge "mumps>=5.7" "metis" "scotch" "scalapack"
# pick the variant that matches openmpi (already installed). If conda-forge has multiple
# mumps builds (mpi vs openmp vs serial), prefer the openmpi-linked one to match the
# existing openmpi in the env.
```
After install, verify: `mamba list | grep -E "mumps|metis|scotch|scalapack"` and `find ~/miniforge3/envs/feaapply -name "libdmumps*" -o -name "libmumps*"`.

A2. Expand OpenRadioss source:
```bash
cd $OR_SRC
git sparse-checkout disable
git pull
ls -d common_source starter engine extlib CMake_Compilers Compiling_tools
```

A3. Read the build instructions in `$OR_SRC/HOWTO.md` (or `INSTALL.md`, or `README.md`) for the implicit-with-MUMPS path. The standard sequence is:
```bash
cd $OR_SRC/starter
./build_script.sh -arch=linux64_gf -nt 16
cd $OR_SRC/engine
./build_script.sh -arch=linux64_gf -mpi=ompi -mumps_root=<MUMPS_PREFIX> -nt 16
```
where `<MUMPS_PREFIX>` is the conda env prefix containing the installed mumps headers and libs (typically `~/miniforge3/envs/feaapply`). Verify the build flag by reading the `engine/build_script.sh` help output first.

A4. The new binaries land at `$OR_SRC/exec/`. Confirm:
```bash
$OR_SRC/exec/starter_linux64_gf -version
$OR_SRC/exec/engine_linux64_gf -help | head -20
```

A5. Place the MUMPS-linked binaries so the runners can pick them up. Strategy (your choice; document it):
- **Replace** the prebuilt `$OR/exec/{starter,engine}_linux64_gf` with the new ones (back up the prebuilt as `_prebuilt.bak`), keeping the same paths so all existing `runner.py` scripts work unchanged. **Recommended.**
- Or **sidecar**: keep both, switch via env var `OR_EXEC_DIR=$OR_SRC/exec`. Update each runner's `OR_DIR` constant.

A6. Smoke-test the new MUMPS-linked engine by running the stage 02 implicit case from the existing runner:
```bash
python tests/stage_02_cantilever_large/runner.py --case implicit_coarse_a1p00
```
The engine must NOT emit `Fatal error: MUMPS required`. If it does, debug the build/link.

A7. Commit + push:
```bash
git add CODEX_BUILD_BRIEF.md   # already authored; just include if not yet committed
git commit -m "build: MUMPS-linked OpenRadioss for /IMPL/NONLIN"
git push
```
(Do NOT commit the OpenRadioss-src tree itself — it is outside the repo. Just commit any FEA_AP-PLY changes that document the new binary location, env vars, or runner adjustments.)

## Phase B — Resume stages 02-16

B1. `git pull --rebase` to sync any human edits.

B2. Re-run stage 02 end-to-end with the MUMPS-linked engine. Expected pass criterion (from `tests/stage_02_cantilever_large/spec.md`): tip displacement within 2% of Bisshopp-Drucker 1945 tabulated values at α∈{1, 3, 5}.

B3. **Update `tests/stage_02_cantilever_large/results/results.json`** with the new run's `verdict` (target PASS), commit `stage 02 cantilever_large: PASS, ...`, push.

B4. **Remove the stage 02 blocker.md** if it now passes (or update it to record that the toolchain issue was resolved).

B5. **Update `RESULTS_SUMMARY.md`** to reflect stage 02 = PASS and that the rest of the ladder is in progress.

B6. Continue stage 03 → 16 per the rules in `CODEX_BRIEF.md` (which still applies in full):
- One stage at a time, in order.
- Commit + push after each stage.
- FAIL → write `blocker.md`, stop.
- INCONCLUSIVE → write `blocker.md`, proceed if dependency chain allows.
- Solid elements only, no GUI, no closed source, no matplotlib, all SI.
- Honour cross-stage consolidation items A-D from `plan/consolidation_review.md`.
- The Kok geometry port (Phase 2.6) lands AFTER stage 10, before stage 11.

## When to stop

- All 16 stages PASS — write final `RESULTS_SUMMARY.md`, commit, push, exit.
- Stage genuinely FAILs (numerics) — write blocker, push, exit.
- Toolchain blocker that is NOT MUMPS (already resolved) — write blocker, push, exit.
- HPC runtime would exceed 48 h on a single shot — downscale per spec, mark INCONCLUSIVE, push, continue.
- Token / time budget — push everything and exit gracefully.

Begin with Phase A.
