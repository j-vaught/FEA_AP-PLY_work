"""
Stage 05 runner: notched dogbone with isotropic ductile damage (LAW22).

Author: J.C. Vaught
Date: 2026-04-29

Runs three mesh-density variants of the notched dogbone deck (coarse / medium /
fine), invokes OpenRadioss starter and engine for each, reads back the
load-displacement history and per-element damage variable, and computes the
windowed RMS difference between the curves up to damage onset. The mesh-
objectivity verdict (PASS or FAIL) is printed at the end.

The pre-damage-onset window is the principal pass criterion: the LAW22 card is
a local (non-regularized) CDM model, so post-peak softening is mesh-dependent
by construction. The runner reports the post-peak RMSE for transparency, but
the verdict is decided exclusively on the [0, u_D] window. See spec.md sections
5 and 8 for the rationale tying this to the OpenRadioss audit MARGINAL verdict.

Toolchain. GMSH (Python API) builds the geometry and writes an Abaqus .inp.
OpenRadioss inp2rad converts to .rad. OpenRadioss starter and engine run the
job inside Lima on macOS, or natively on Linux. Vortex-Radioss reads the T01
time history and the .anim damage field. Numpy and pandas do the postprocessing.

This script does not draw figures; per project preference, plots are authored
in Typst + CeTZ from the CSV the runner writes.

Usage examples.

    python runner.py                         # default: all three meshes, A36
    python runner.py --metal dp780           # alternate parameter block
    python runner.py --only medium           # run a single mesh
    python runner.py --skip-solve            # postprocess existing T01 files
    python runner.py --engine-cores 8        # parallel engine
    python runner.py --solver explicit       # force the dynamic-relaxation path

The runner is deliberately a single file; per stage 1 / 2 conventions in the
test suite, each stage owns its own runner with no cross-stage imports beyond
the standard library plus numpy / pandas / vortex-radioss / gmsh.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

STAGE_DIR = Path(__file__).resolve().parent
WORK_DIR = STAGE_DIR / "work"
RESULTS_DIR = STAGE_DIR / "results"

# Lima VM name from master_plan.md section 7.
LIMA_VM = os.environ.get("OR_LIMA_VM", "or")

# OpenRadioss binary names (Linux ARM64 inside Lima per audit section 4.2).
STARTER_BIN = os.environ.get("OR_STARTER_BIN", "starter_linuxa64")
ENGINE_BIN = os.environ.get("OR_ENGINE_BIN", "engine_linuxa64")
INP2RAD_BIN = os.environ.get("OR_INP2RAD_BIN", "inp2rad")

# Mesh resolutions from spec.md section 4.
MESH_SIZES_MM = {
    "coarse": 1.00,
    "medium": 0.50,
    "fine": 0.25,
}

# Damage-onset detection threshold per spec.md section 8.
D_ONSET_THRESHOLD = 0.01

# Mesh-objectivity tolerance per spec.md section 8.
RMSE_PASS_TOLERANCE = 0.05  # 5% of peak load


# ---------------------------------------------------------------------------
# Material parameter blocks (spec.md section 5)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MaterialCard:
    """LAW22 material parameters in consistent SI (kg, m, s, Pa)."""

    name: str
    rho: float           # density, kg/m^3
    E: float             # Young's modulus, Pa
    nu: float            # Poisson ratio
    sigma_y: float       # yield stress, Pa
    hardening_pts: list  # list of (eps_p, sigma) tuples in Pa
    eps_p_d: float       # damage onset plastic strain
    eps_p_r: float       # full-damage / element-deletion plastic strain
    d_max: float = 0.999


def material_a36() -> MaterialCard:
    return MaterialCard(
        name="A36",
        rho=7850.0,
        E=2.00e11,
        nu=0.30,
        sigma_y=2.50e8,
        hardening_pts=[
            (0.000, 2.50e8),
            (0.020, 3.00e8),
            (0.050, 3.50e8),
            (0.100, 4.00e8),
            (0.200, 4.50e8),
            (0.500, 4.50e8),
        ],
        eps_p_d=0.05,
        eps_p_r=0.50,
    )


def material_dp780() -> MaterialCard:
    return MaterialCard(
        name="DP780",
        rho=7850.0,
        E=2.00e11,
        nu=0.30,
        sigma_y=5.00e8,
        hardening_pts=[
            (0.000, 5.00e8),
            (0.020, 6.20e8),
            (0.050, 7.00e8),
            (0.100, 7.80e8),
            (0.200, 7.80e8),
        ],
        eps_p_d=0.04,
        eps_p_r=0.20,
    )


# ---------------------------------------------------------------------------
# Geometry parameters (spec.md section 3) -- millimetres
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Geometry:
    L_tot_mm: float = 200.0
    L_g_mm: float = 60.0
    w_g_mm: float = 12.5
    w_n_mm: float = 8.0
    r_n_mm: float = 2.25
    w_grip_mm: float = 25.0
    r_f_mm: float = 12.5
    t_mm: float = 6.0


GEOM = Geometry()


# ---------------------------------------------------------------------------
# GMSH meshing
# ---------------------------------------------------------------------------

def build_mesh(h_mm: float, out_inp: Path) -> dict:
    """Build the notched-dogbone HEXA8 mesh and export an Abaqus .inp.

    Returns a small dict with element count, node count, and a few derived
    quantities used by inp2rad sanity checks.

    The geometry is built as a 2D outline with the pair of semicircular notches
    at midspan, then extruded through-thickness to HEXA8 with the number of
    layers chosen to keep aspect ratio near unity in the ligament.
    """
    import gmsh

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(f"stage05_h{h_mm:.2f}")

        g = GEOM
        # 2D outline points (mm). Half the dogbone is mirrored about y=0; we
        # build the full outline directly so the mesh is symmetric without
        # any explicit symmetry plane (see spec.md section 3 rationale).
        L = g.L_tot_mm
        Lg = g.L_g_mm
        wg = g.w_g_mm
        wn = g.w_n_mm
        rn = g.r_n_mm
        wgrip = g.w_grip_mm
        rf = g.r_f_mm

        # x-coordinates of key transitions
        x_grip_R = (L - Lg) / 2.0 - rf      # end of straight grip (left)
        x_gauge_L = (L - Lg) / 2.0          # start of gauge straight
        x_notch_C = L / 2.0                 # notch center
        x_gauge_R = (L + Lg) / 2.0          # end of gauge straight
        x_grip_L2 = (L + Lg) / 2.0 + rf     # start of straight grip (right)

        # Build half-outline (top, y>0) and mirror by adding y<0 points later.
        # Use OpenCASCADE for boolean cuts of the semicircular notches.
        occ = gmsh.model.occ

        # Base rectangle at grip width
        grip_L = occ.addRectangle(0.0, -wgrip / 2.0, 0.0, x_grip_R, wgrip)
        grip_R = occ.addRectangle(x_grip_L2, -wgrip / 2.0, 0.0, L - x_grip_L2, wgrip)

        # Gauge rectangle (reduced width)
        gauge = occ.addRectangle(x_gauge_L, -wg / 2.0, 0.0, Lg, wg)

        # Fillet trapezoids: build each shoulder as a rectangle then subtract
        # the fillet quarter-disks.  This is a stable OCC pattern that does
        # not produce sliver faces.
        shL_x = x_grip_R
        shL_w = x_gauge_L - x_grip_R
        shoulder_L = occ.addRectangle(shL_x, -wgrip / 2.0, 0.0, shL_w, wgrip)
        shoulder_R = occ.addRectangle(x_gauge_R, -wgrip / 2.0, 0.0, shL_w, wgrip)

        # Subtract fillet quarter-disks on the four shoulder corners.
        # Top-left fillet: disk centered at (x_grip_R + rf, +wg/2 + rf)
        d1 = occ.addDisk(x_grip_R + rf, +wg / 2.0 + rf, 0.0, rf, rf)
        d2 = occ.addDisk(x_grip_R + rf, -wg / 2.0 - rf, 0.0, rf, rf)
        d3 = occ.addDisk(x_gauge_R - rf, +wg / 2.0 + rf, 0.0, rf, rf)
        d4 = occ.addDisk(x_gauge_R - rf, -wg / 2.0 - rf, 0.0, rf, rf)

        # Fuse the base shapes
        all_outline, _ = occ.fuse(
            [(2, grip_L), (2, gauge), (2, grip_R), (2, shoulder_L), (2, shoulder_R)],
            [],
            removeObject=True, removeTool=True,
        )
        # Subtract the four corner disks to form fillets. Note: only the
        # quadrants that cut into the grip rectangle should be inverted; this
        # is an approximation acceptable at the verification-mesh resolution.
        outline_after_fillet, _ = occ.cut(
            all_outline,
            [(2, d1), (2, d2), (2, d3), (2, d4)],
            removeObject=True, removeTool=True,
        )

        # Subtract the symmetric notches at midspan.
        notch_top = occ.addDisk(x_notch_C, +wg / 2.0, 0.0, rn, rn)
        notch_bot = occ.addDisk(x_notch_C, -wg / 2.0, 0.0, rn, rn)
        outline_final, _ = occ.cut(
            outline_after_fillet,
            [(2, notch_top), (2, notch_bot)],
            removeObject=True, removeTool=True,
        )

        occ.synchronize()

        # Extrude through-thickness to solids.
        n_through = max(int(round(g.t_mm / h_mm)), 4)
        # Mesh-size policy: refine in the gauge ligament, coarsen in grips.
        # We do this via a Box field around the notch ligament.
        gmsh.model.mesh.field.add("Box", 1)
        gmsh.model.mesh.field.setNumber(1, "VIn", h_mm)
        gmsh.model.mesh.field.setNumber(1, "VOut", max(2.0 * h_mm, 1.5))
        gmsh.model.mesh.field.setNumber(1, "XMin", x_notch_C - 1.5 * Lg / 4.0)
        gmsh.model.mesh.field.setNumber(1, "XMax", x_notch_C + 1.5 * Lg / 4.0)
        gmsh.model.mesh.field.setNumber(1, "YMin", -wg / 2.0)
        gmsh.model.mesh.field.setNumber(1, "YMax", +wg / 2.0)
        gmsh.model.mesh.field.setAsBackgroundMesh(1)
        gmsh.option.setNumber("Mesh.Algorithm", 8)        # Frontal-Delaunay for quads
        gmsh.option.setNumber("Mesh.RecombineAll", 1)     # quads then hexes on extrude
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)

        # Extrude: returns list of (dim, tag) of generated entities.
        extruded = occ.extrude(
            outline_final, 0.0, 0.0, g.t_mm,
            numElements=[n_through], recombine=True,
        )
        occ.synchronize()

        # Physical groups for boundary condition tagging.
        # Left grip face: x = 0 plane
        # Right grip face: x = L plane
        eps_tol = 1e-3
        bnd = gmsh.model.getEntities(2)
        left_faces, right_faces = [], []
        for dim, tag in bnd:
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
            if abs(xmin - 0.0) < eps_tol and abs(xmax - 0.0) < eps_tol:
                left_faces.append(tag)
            elif abs(xmin - L) < eps_tol and abs(xmax - L) < eps_tol:
                right_faces.append(tag)
        if not left_faces or not right_faces:
            raise RuntimeError(
                "GMSH could not identify left/right grip faces by bounding-box "
                "filter. Tighten eps_tol or inspect the OCC model."
            )
        gmsh.model.addPhysicalGroup(2, left_faces, tag=101, name="LEFT_GRIP")
        gmsh.model.addPhysicalGroup(2, right_faces, tag=102, name="RIGHT_GRIP")

        # Volume physical group (the part).
        vols = [t for d, t in gmsh.model.getEntities(3)]
        gmsh.model.addPhysicalGroup(3, vols, tag=1, name="DOGBONE")

        gmsh.model.mesh.generate(3)
        # Convert tetrahedra (if any sneaked in) to hexes via subdivide:
        # GMSH 4.x subdivision of HEX-EXTRUSION should already be all-hex;
        # we leave a guard in case the OCC topology produced wedges.
        try:
            gmsh.model.mesh.recombine()
        except Exception:
            pass

        # Counts
        node_tags, _, _ = gmsh.model.mesh.getNodes()
        elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
        n_nodes = int(len(node_tags))
        n_elems = int(sum(len(t) for t in elem_tags))

        out_inp.parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(out_inp))

    finally:
        gmsh.finalize()

    return {"n_nodes": n_nodes, "n_elems": n_elems}


# ---------------------------------------------------------------------------
# Deck templating
# ---------------------------------------------------------------------------

STARTER_TEMPLATE = """\
#RADIOSS STARTER
/BEGIN
Stage 05 dogbone with isotropic ductile damage (LAW22) -- mesh {mesh_label}
      0       0
/UNIT/1
kg                  m                   s
/MAT/LAW22/1
{mat_name}_LAW22
{rho:>20.6e}{E:>20.6e}{nu:>20.6f}
{sigma_y:>20.6e}{eps_p_d:>20.6f}{eps_p_r:>20.6f}{d_max:>20.6f}
# (hardening function 1 referenced below)
/FUNCT/1
HARDENING_{mat_name}
{hardening_pairs}
/PROP/TYPE14/1
SOLID_GENERAL
{ihkt:>10d}{isolid:>10d}{ismstr:>10d}
/PART/1
DOGBONE_PART     1     1
# nodes and bricks below come from inp2rad conversion
#include "mesh_block.rad"
/GRNOD/SURF/1
LEFT_GRIP_NODES   101
/GRNOD/SURF/2
RIGHT_GRIP_NODES  102
/BCS/1
LEFT_ENCASTRE
1 1 1 1 1 1                 1
/BCS/2
RIGHT_TRANSV_FIX
0 1 1 0 0 0                 2
/IMPDISP/1
RIGHT_X_PULL
1 0 0                       2       3
# function 3 is the displacement ramp, end value u_end (m)
/FUNCT/3
DISP_RAMP
              0.0                 0.0
              1.0     {u_end_m:.6e}
{solver_block}
/TH/PART/1
PART_OUTPUT       1
DEF DISP FORC ENER
/TH/NODE/1
GRIP_NODE
{control_node:>10d}
DEF DISP FORC
/ANIM/ELEM/DAMA
/H3D/ELEM/DAMA
/END
"""

ENGINE_TEMPLATE = """\
#RADIOSS ENGINE
/RUN/STAGE05_{mesh_label}/1
              1.0
/TFILE
            0.001
/ANIM/DT
             0.02
/H3D/DT
             0.02
/STOP
/END
"""

SOLVER_BLOCK_IMPL = """\
/IMPL/QSTAT
/IMPL/SOLVER/MUMPS
/IMPL/NONLIN/SMDISP
       1                                    50      1.0e-3
/IMPL/DT/STOP
       1.0e-4              1.0
"""

SOLVER_BLOCK_EXPL = """\
/DT/BRICK/1
              0.9          1.0e-7
/DAMP/1
GLOBAL_RAYLEIGH
              0.05            0.0
/MASS/SCAL
              1.0e-7
"""


def render_starter(
    mat: MaterialCard,
    mesh_label: str,
    u_end_m: float,
    solver: str,
    control_node: int,
) -> str:
    """Render the starter deck text. Field widths follow Radioss fixed format."""
    hardening_pairs = "\n".join(
        f"{ep:>20.6e}{sig:>20.6e}" for ep, sig in mat.hardening_pts
    )
    if solver == "implicit":
        solver_block = SOLVER_BLOCK_IMPL
    elif solver == "explicit":
        solver_block = SOLVER_BLOCK_EXPL
    else:
        raise ValueError(f"unknown solver mode {solver!r}")
    return STARTER_TEMPLATE.format(
        mesh_label=mesh_label,
        mat_name=mat.name,
        rho=mat.rho,
        E=mat.E,
        nu=mat.nu,
        sigma_y=mat.sigma_y,
        eps_p_d=mat.eps_p_d,
        eps_p_r=mat.eps_p_r,
        d_max=mat.d_max,
        hardening_pairs=hardening_pairs,
        ihkt=2,             # tabulated isotropic hardening
        isolid=14,          # /PROP/TYPE14 solid
        ismstr=10,          # large strain solid formulation
        u_end_m=u_end_m,
        solver_block=solver_block,
        control_node=control_node,
    )


def render_engine(mesh_label: str) -> str:
    return ENGINE_TEMPLATE.format(mesh_label=mesh_label)


# ---------------------------------------------------------------------------
# Toolchain invocation (Lima + OpenRadioss)
# ---------------------------------------------------------------------------

def _have_lima() -> bool:
    return shutil.which("limactl") is not None


def _have_native_or() -> bool:
    return shutil.which(STARTER_BIN) is not None and shutil.which(ENGINE_BIN) is not None


def _or_command(argv: list[str], cwd: Path) -> subprocess.CompletedProcess:
    """Run an OpenRadioss command natively or inside Lima."""
    if _have_native_or():
        return subprocess.run(argv, cwd=str(cwd), capture_output=True, text=True)
    if _have_lima():
        full = ["limactl", "shell", LIMA_VM, "--workdir", str(cwd), "--"] + argv
        return subprocess.run(full, capture_output=True, text=True)
    raise RuntimeError(
        "Neither native OpenRadioss binaries nor Lima detected. Install per "
        "master_plan.md section 7."
    )


def detect_implicit_capability() -> bool:
    """Probe the engine binary for MUMPS / implicit support.

    Heuristic: invoke `engine -h` and grep for IMPL keywords. If unavailable,
    return False so the runner falls back to explicit dynamic relaxation per
    spec.md section 6 Option B.
    """
    try:
        proc = _or_command([ENGINE_BIN, "-h"], cwd=STAGE_DIR)
    except Exception:
        return False
    text = (proc.stdout or "") + (proc.stderr or "")
    return "IMPL" in text.upper() and "MUMPS" in text.upper()


def run_inp2rad(inp_path: Path, rad_out: Path) -> None:
    """Convert Abaqus .inp to OpenRadioss .rad via OpenRadioss inp2rad."""
    cwd = inp_path.parent
    proc = _or_command([INP2RAD_BIN, str(inp_path.name)], cwd=cwd)
    if proc.returncode != 0:
        raise RuntimeError(
            f"inp2rad failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    # inp2rad writes <name>_0000.rad next to the .inp; we standardize the name.
    converted = cwd / (inp_path.stem + "_0000.rad")
    if not converted.exists():
        # Some inp2rad versions name it differently; fall back to glob.
        candidates = sorted(cwd.glob(inp_path.stem + "*.rad"))
        if not candidates:
            raise RuntimeError("inp2rad produced no .rad output")
        converted = candidates[0]
    rad_out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(converted, rad_out)


def run_starter(deck_0000: Path, n_threads: int = 1) -> None:
    proc = _or_command(
        [STARTER_BIN, "-i", deck_0000.name, "-nt", str(n_threads)],
        cwd=deck_0000.parent,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"OpenRadioss starter failed on {deck_0000}:\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )


def run_engine(deck_0001: Path, n_threads: int = 4) -> None:
    proc = _or_command(
        [ENGINE_BIN, "-i", deck_0001.name, "-nt", str(n_threads)],
        cwd=deck_0001.parent,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"OpenRadioss engine failed on {deck_0001}:\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )


# ---------------------------------------------------------------------------
# Time-history extraction (Vortex-Radioss reader)
# ---------------------------------------------------------------------------

def read_time_history(t01_path: Path) -> pd.DataFrame:
    """Return DataFrame with columns ('time_s', 'displacement_m', 'force_N').

    The runner pulls the right-grip control-node displacement and the
    right-grip resultant reaction force.
    """
    try:
        from vortex_radioss import T01Reader  # type: ignore
    except ImportError:
        # Fallback: parse the ASCII T01 dump if Vortex-Radioss is unavailable.
        return _read_time_history_ascii(t01_path)

    reader = T01Reader(str(t01_path))
    df = reader.as_dataframe()
    # Expected columns include Node DX and Node FX for the control node;
    # exact column naming depends on Vortex-Radioss version.
    disp_col = next(c for c in df.columns if "DX" in c.upper() and "GRIP" in c.upper())
    force_col = next(c for c in df.columns if "FX" in c.upper() and "GRIP" in c.upper())
    out = pd.DataFrame({
        "time_s": df["time"].to_numpy(),
        "displacement_m": df[disp_col].to_numpy(),
        "force_N": df[force_col].to_numpy(),
    })
    return out


def _read_time_history_ascii(t01_path: Path) -> pd.DataFrame:
    """Fallback ASCII T01 parser; conservative regex on whitespace columns."""
    rows = []
    with open(t01_path, "r", errors="ignore") as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                t = float(parts[0]); u = float(parts[1]); f = float(parts[2])
            except ValueError:
                continue
            rows.append((t, u, f))
    if not rows:
        raise RuntimeError(f"No parseable rows in {t01_path}")
    return pd.DataFrame(rows, columns=["time_s", "displacement_m", "force_N"])


def read_max_damage(anim_path: Path) -> np.ndarray:
    """Return array of max(D) over the gauge ligament at each animation frame.

    Tries Vortex-Radioss; falls back to returning a synthetic zero array (so
    the runner can still report load-displacement objectivity if the damage
    output is unreadable for any reason).
    """
    try:
        from vortex_radioss import AnimReader  # type: ignore
    except ImportError:
        return np.zeros(0)
    reader = AnimReader(str(anim_path))
    frames = reader.frames()
    out = np.zeros(len(frames))
    for k, fr in enumerate(frames):
        d = fr.scalar("DAMA")
        if d is None or len(d) == 0:
            out[k] = 0.0
        else:
            out[k] = float(np.max(d))
    return out


# ---------------------------------------------------------------------------
# Mesh-objectivity computation (spec.md section 8)
# ---------------------------------------------------------------------------

def windowed_rmse(
    u_a: np.ndarray, f_a: np.ndarray,
    u_b: np.ndarray, f_b: np.ndarray,
    u_lo: float, u_hi: float,
    f_norm: float,
) -> float:
    """Compute RMSE of f_a(u) vs f_b(u) over [u_lo, u_hi], normalized by f_norm.

    Both curves are resampled onto a shared monotone displacement grid via
    linear interpolation (curves are monotone-non-decreasing in displacement
    by construction of the imposed-displacement BC).
    """
    if u_hi <= u_lo:
        return float("nan")
    u_grid = np.linspace(u_lo, u_hi, 256)
    fa = np.interp(u_grid, u_a, f_a)
    fb = np.interp(u_grid, u_b, f_b)
    diff = (fa - fb) / max(abs(f_norm), 1e-30)
    return float(np.sqrt(np.mean(diff * diff)))


def detect_damage_onset_disp(
    u_medium: np.ndarray, dmax_medium_per_t: np.ndarray, t_medium: np.ndarray,
) -> float:
    """Return the displacement (m) at which max(D) first exceeds threshold,
    interpolated against the medium-mesh time history. If damage data are
    unavailable (empty array), returns the displacement at peak load on the
    medium mesh as a documented fallback.
    """
    if dmax_medium_per_t.size == 0:
        # No damage data; use peak-load displacement as the proxy onset.
        return float(u_medium[int(np.argmax(np.abs(np.gradient(u_medium))))])
    # Find first index where damage crosses threshold.
    idx = np.argmax(dmax_medium_per_t > D_ONSET_THRESHOLD)
    if dmax_medium_per_t[idx] <= D_ONSET_THRESHOLD:
        # Threshold never crossed; return last displacement.
        return float(u_medium[-1])
    # Interpolate between idx-1 and idx in time, then map to displacement.
    if idx == 0:
        return float(u_medium[0])
    t_cross = np.interp(
        D_ONSET_THRESHOLD,
        [dmax_medium_per_t[idx - 1], dmax_medium_per_t[idx]],
        [t_medium[idx - 1], t_medium[idx]],
    )
    return float(np.interp(t_cross, t_medium, u_medium))


# ---------------------------------------------------------------------------
# Per-mesh orchestration
# ---------------------------------------------------------------------------

@dataclass
class MeshRun:
    label: str
    h_mm: float
    work_dir: Path
    inp_path: Path
    rad_starter: Path
    rad_engine: Path
    t01_path: Path
    anim_path: Path
    n_nodes: int = 0
    n_elems: int = 0
    solver: str = "explicit"
    fd_curve: Optional[pd.DataFrame] = None
    dmax_per_t: np.ndarray = field(default_factory=lambda: np.zeros(0))


def run_one_mesh(
    label: str, h_mm: float, mat: MaterialCard, solver: str,
    engine_cores: int, skip_solve: bool,
) -> MeshRun:
    work = WORK_DIR / label
    work.mkdir(parents=True, exist_ok=True)
    inp = work / f"dogbone_{label}.inp"
    starter = work / f"job_{label}_0000.rad"
    engine_deck = work / f"job_{label}_0001.rad"
    t01 = work / f"job_{label}T01"
    anim = work / f"job_{label}A001"

    run = MeshRun(
        label=label, h_mm=h_mm, work_dir=work, inp_path=inp,
        rad_starter=starter, rad_engine=engine_deck,
        t01_path=t01, anim_path=anim, solver=solver,
    )

    if not skip_solve:
        # 1. Build mesh.
        meta = build_mesh(h_mm, inp)
        run.n_nodes = meta["n_nodes"]
        run.n_elems = meta["n_elems"]
        # 2. Convert to OpenRadioss .rad.
        mesh_block = work / "mesh_block.rad"
        run_inp2rad(inp, mesh_block)
        # 3. Render starter and engine decks.
        starter_text = render_starter(
            mat, label, u_end_m=6.0e-3, solver=solver, control_node=1,
        )
        starter.write_text(starter_text)
        engine_deck.write_text(render_engine(label))
        # 4. Run starter then engine.
        run_starter(starter, n_threads=1)
        run_engine(engine_deck, n_threads=engine_cores)

    # 5. Read time history.
    if t01.exists():
        run.fd_curve = read_time_history(t01)
    else:
        # Locate any T01-suffixed file in the workdir
        candidates = sorted(work.glob("*T01*"))
        if candidates:
            run.fd_curve = read_time_history(candidates[0])
        else:
            print(f"  [warn] {label}: no T01 file found", file=sys.stderr)
            run.fd_curve = pd.DataFrame(columns=["time_s", "displacement_m", "force_N"])

    if anim.exists():
        run.dmax_per_t = read_max_damage(anim)
    else:
        candidates = sorted(work.glob("*A0*"))
        if candidates:
            run.dmax_per_t = read_max_damage(candidates[0])

    return run


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--metal", choices=["a36", "dp780"], default="a36")
    p.add_argument("--only", choices=list(MESH_SIZES_MM.keys()), default=None,
                   help="run only one mesh (debug aid)")
    p.add_argument("--skip-solve", action="store_true",
                   help="skip GMSH/inp2rad/starter/engine; postprocess only")
    p.add_argument("--engine-cores", type=int, default=4)
    p.add_argument("--solver", choices=["auto", "implicit", "explicit"],
                   default="auto",
                   help="auto = probe MUMPS, fall back to explicit if absent")
    args = p.parse_args(argv)

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    mat = material_a36() if args.metal == "a36" else material_dp780()

    if args.solver == "auto":
        if not args.skip_solve and detect_implicit_capability():
            solver = "implicit"
        else:
            solver = "explicit"
    else:
        solver = args.solver
    print(f"[stage 05] material = {mat.name}")
    print(f"[stage 05] solver path = {solver} "
          f"({'/IMPL/QSTAT' if solver == 'implicit' else 'explicit dynamic relaxation'})")

    labels = [args.only] if args.only else list(MESH_SIZES_MM.keys())
    runs: dict[str, MeshRun] = {}
    for label in labels:
        h = MESH_SIZES_MM[label]
        print(f"[stage 05] running mesh '{label}' (h = {h:.2f} mm) ...")
        t0 = time.time()
        runs[label] = run_one_mesh(
            label, h, mat, solver,
            engine_cores=args.engine_cores, skip_solve=args.skip_solve,
        )
        dt = time.time() - t0
        print(f"  done in {dt:.1f} s "
              f"(nodes = {runs[label].n_nodes}, elems = {runs[label].n_elems})")

    # ---- Mesh-objectivity analysis (only meaningful with all three) ----
    if set(runs) != set(MESH_SIZES_MM):
        print("[stage 05] partial mesh set; skipping objectivity verdict.")
        return 0

    coarse = runs["coarse"].fd_curve
    medium = runs["medium"].fd_curve
    fine = runs["fine"].fd_curve
    if any(df is None or df.empty for df in (coarse, medium, fine)):
        print("[stage 05] missing time-history data; cannot compute RMSE.",
              file=sys.stderr)
        return 2

    # Damage-onset displacement from the medium mesh.
    u_D = detect_damage_onset_disp(
        medium["displacement_m"].to_numpy(),
        runs["medium"].dmax_per_t,
        medium["time_s"].to_numpy(),
    )
    # Peak load on the medium mesh used as normalization.
    f_peak_medium = float(np.max(np.abs(medium["force_N"].to_numpy())))
    u_max = min(
        float(coarse["displacement_m"].max()),
        float(medium["displacement_m"].max()),
        float(fine["displacement_m"].max()),
    )

    rmse_pre_cm = windowed_rmse(
        coarse["displacement_m"].to_numpy(), coarse["force_N"].to_numpy(),
        medium["displacement_m"].to_numpy(), medium["force_N"].to_numpy(),
        u_lo=0.0, u_hi=u_D, f_norm=f_peak_medium,
    )
    rmse_pre_mf = windowed_rmse(
        medium["displacement_m"].to_numpy(), medium["force_N"].to_numpy(),
        fine["displacement_m"].to_numpy(), fine["force_N"].to_numpy(),
        u_lo=0.0, u_hi=u_D, f_norm=f_peak_medium,
    )
    rmse_post_cm = windowed_rmse(
        coarse["displacement_m"].to_numpy(), coarse["force_N"].to_numpy(),
        medium["displacement_m"].to_numpy(), medium["force_N"].to_numpy(),
        u_lo=u_D, u_hi=u_max, f_norm=f_peak_medium,
    )
    rmse_post_mf = windowed_rmse(
        medium["displacement_m"].to_numpy(), medium["force_N"].to_numpy(),
        fine["displacement_m"].to_numpy(), fine["force_N"].to_numpy(),
        u_lo=u_D, u_hi=u_max, f_norm=f_peak_medium,
    )

    pass_pre = (rmse_pre_cm <= RMSE_PASS_TOLERANCE) and \
               (rmse_pre_mf <= RMSE_PASS_TOLERANCE)
    verdict = "PASS" if pass_pre else "FAIL"

    # Write CSV (combined load-displacement) and summary JSON for Typst+CeTZ.
    u_grid = np.linspace(0.0, u_max, 1024)
    fc = np.interp(u_grid, coarse["displacement_m"], coarse["force_N"])
    fm = np.interp(u_grid, medium["displacement_m"], medium["force_N"])
    ff = np.interp(u_grid, fine["displacement_m"],   fine["force_N"])
    dc = np.interp(u_grid,
                   np.linspace(0, u_max, runs["coarse"].dmax_per_t.size or 1),
                   runs["coarse"].dmax_per_t if runs["coarse"].dmax_per_t.size
                   else np.zeros(1))
    dm = np.interp(u_grid,
                   np.linspace(0, u_max, runs["medium"].dmax_per_t.size or 1),
                   runs["medium"].dmax_per_t if runs["medium"].dmax_per_t.size
                   else np.zeros(1))
    df_ = np.interp(u_grid,
                    np.linspace(0, u_max, runs["fine"].dmax_per_t.size or 1),
                    runs["fine"].dmax_per_t if runs["fine"].dmax_per_t.size
                    else np.zeros(1))

    out_csv = RESULTS_DIR / "mesh_objectivity.csv"
    pd.DataFrame({
        "displacement_m": u_grid,
        "F_coarse_N": fc, "F_medium_N": fm, "F_fine_N": ff,
        "D_coarse": dc, "D_medium": dm, "D_fine": df_,
    }).to_csv(out_csv, index=False)

    summary = {
        "stage": 5,
        "material": mat.name,
        "solver": solver,
        "u_damage_onset_m": u_D,
        "F_peak_medium_N": f_peak_medium,
        "RMSE_pre_onset_coarse_vs_medium": rmse_pre_cm,
        "RMSE_pre_onset_medium_vs_fine": rmse_pre_mf,
        "RMSE_post_onset_coarse_vs_medium": rmse_post_cm,
        "RMSE_post_onset_medium_vs_fine": rmse_post_mf,
        "tolerance_pre_onset": RMSE_PASS_TOLERANCE,
        "verdict_pre_onset": verdict,
        "post_onset_divergence_expected": True,
        "post_onset_divergence_note": (
            "LAW22 is a local CDM model; mesh-dependent post-peak softening is "
            "expected and is not part of the pass criterion. Nonlocal "
            "regularization per Pijaudier-Cabot and Mazars 1989 is not "
            "available in OpenRadioss (DOCUMENTATION NOT LOCATED, audit row 5)."
        ),
    }
    out_json = RESULTS_DIR / "mesh_objectivity_summary.json"
    out_json.write_text(json.dumps(summary, indent=2))

    print()
    print("=" * 72)
    print(f"[stage 05] mesh-objectivity verdict (pre-damage-onset window)")
    print("=" * 72)
    print(f"  damage-onset displacement u_D = {u_D*1e3:.4f} mm")
    print(f"  peak load (medium)             = {f_peak_medium:.3f} N")
    print(f"  RMSE pre  (coarse vs medium)   = {rmse_pre_cm*100:.3f} % of peak")
    print(f"  RMSE pre  (medium vs fine)     = {rmse_pre_mf*100:.3f} % of peak")
    print(f"  tolerance                      = {RMSE_PASS_TOLERANCE*100:.1f} % of peak")
    print(f"  VERDICT                        = {verdict}")
    print()
    print(f"  RMSE post (coarse vs medium)   = {rmse_post_cm*100:.3f} % of peak  "
          f"[expected divergence; documented]")
    print(f"  RMSE post (medium vs fine)     = {rmse_post_mf*100:.3f} % of peak  "
          f"[expected divergence; documented]")
    print()
    print(f"  CSV:  {out_csv}")
    print(f"  JSON: {out_json}")

    return 0 if pass_pre else 1


if __name__ == "__main__":
    sys.exit(main())
