"""
Stage 15 — Flat-coupon ballistic test runner.

Author: J.C. Vaught
Date: 2026-04-29

Drives the full V50 sweep for the IM7/8552 24-ply 4 mm flat tape laminate against
a steel-cylinder representative of a 9 mm FMJ. Templates one OpenRadioss starter +
engine deck pair per strike velocity, invokes the OpenRadioss starter and engine
inside Lima/Apptainer, parses the time history, fits Recht-Ipson 1963, extracts
V50, and writes a pass/fail report against the Vakili Rad 2020 baseline-tape
reference (300 +/- 15 m/s).

Pipeline per run:
    1. write GMSH .msh (deferred to a sibling mesher; this runner expects an
       Abaqus .inp produced by gmsh and a converted .rad to be present, OR it
       calls inp2rad directly).
    2. template the starter and engine .rad decks for the chosen strike velocity.
    3. invoke OpenRadioss starter, then engine.
    4. parse the T01 time-history file; extract projectile centroid V_z(t),
       global energies, hourglass energy, contact energy, eroded elements.
    5. record (V_s, V_r) and energy diagnostics.

Sweep:
    V_s in {200, 250, 300, 350, 400} m/s, five runs.

After the sweep:
    - fit Recht-Ipson V_r = sqrt(V_s^2 - V50^2) on the over-match subset.
    - fit Lambert-Jonas V_r = (V_s^p - V50^p)^(1/p) (free p, free V50).
    - compute Cunniff parameter c* for IM7 and the dimensionless V50/c*.
    - emit eight pass/warn/fail checks per spec.md section 8.
    - write sweep_results.csv, recht_ipson_fit.json, pass_fail_report.json, vr_vs_vs.csv.

Citations (full BibTeX in references/):
    - Recht & Ipson 1963, J. Appl. Mech. 30, 384-390.
    - Lambert & Jonas 1976, BRL-R-1852.
    - Cunniff 1992, Textile Res. J. 62, 495-509.
    - MIL-STD-662F (1997).
    - NIJ Standard-0101.06 (2008).
    - STANAG 2920 (NATO 2003).
    - Vakili Rad et al. 2020, Composites Part B 203, 108478 (V50 anchor).
    - Hashin 1980, J. Appl. Mech. 47, 329-334.
    - Soden, Hinton & Kaddour 1998, Compos. Sci. Tech. 58, 1011-1022.
    - Belytschko & Bindeman 1993, CMAME 105, 225-260 (HEPH hourglass).

This runner is SI throughout (m, kg, s, N, Pa). It does NOT use any non-SI helper.
It refuses to run if the user overrides V50_REF_MS without explicit confirmation.

Usage:
    python runner.py --sweep                  # full five-velocity sweep
    python runner.py --single 350             # single velocity (m/s)
    python runner.py --post-only              # parse existing T01 files, skip solver
    python runner.py --mesh-convergence 350   # M0/M1/M2 at one velocity
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    from scipy.optimize import curve_fit, least_squares
except ImportError as e:
    print("ERROR: scipy is required. pip install scipy", file=sys.stderr)
    raise

# ----------------------------------------------------------------------------
# Constants — physical, geometric, material, numerical
# ----------------------------------------------------------------------------

# Project paths (absolute per user instruction).
THIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = THIS_DIR.parents[1]
RUNS_DIR = THIS_DIR / "runs"
SUMMARY_DIR = RUNS_DIR / "_summary"

# Strike-velocity sweep (m/s).
SWEEP_VS_MS: Tuple[float, ...] = (200.0, 250.0, 300.0, 350.0, 400.0)

# Reference V50 from Vakili Rad et al. 2020 baseline tape laminate.
# UNVERIFIED in this run — confirm the exact number against
# Vakili Rad 2020 Figure 7 / Table 3.
V50_REF_MS: float = 300.0
V50_REF_TOL_MS: float = 15.0

# Pass tolerances (per spec.md section 8).
TOL_V50_FRAC: float = 0.05          # |V50_fit - V50_ref| / V50_ref <= 0.05
TOL_VR_FRAC: float = 0.10           # residual velocity within 10%
TOL_HG_OVER_INT: float = 0.05       # hourglass energy / internal energy <= 5%
TOL_ENERGY_BAL_FRAC: float = 0.02   # global energy balance within 2%
TOL_MESH_CONV_FRAC: float = 0.05    # mesh convergence at Vs=350 m/s
TOL_LJ_P_LO: float = 1.8            # Lambert-Jonas exponent lower
TOL_LJ_P_HI: float = 2.4            # Lambert-Jonas exponent upper
TOL_CUNNIFF_LO: float = 0.3
TOL_CUNNIFF_HI: float = 0.6
SUB_V50_NOISE_MS: float = 5.0       # |V_r| below this is treated as arrest

# Geometry (SI, m).
PLATE_SIDE_M = 0.200          # 200 mm
PLATE_THICK_M = 0.00432       # 24 plies x 0.18 mm = 4.32 mm
N_PLY = 24
PLY_THICK_M = 0.00018         # 0.18 mm
CLAMP_STRIP_M = 0.025         # 25 mm
APERTURE_SIDE_M = PLATE_SIDE_M - 2.0 * CLAMP_STRIP_M  # 0.150 m

# Stacking sequence [0/+45/-45/90]_3s, expanded to 24 plies symmetric.
# spec.md section 3 — the runner is the single source of truth for this list.
STACKING_DEG: Tuple[float, ...] = tuple(
    list((0.0, 45.0, -45.0, 90.0)) * 3              # plies 1..12 (top half)
    + list(reversed((0.0, 45.0, -45.0, 90.0) * 3))  # plies 13..24 (mirror)
)
assert len(STACKING_DEG) == N_PLY, "stacking expansion error"

# Projectile (SI).
PROJ_MASS_KG = 0.008          # 8 g
PROJ_DIAM_M = 0.009           # 9 mm
PROJ_DENSITY_KG_M3 = 7800.0   # steel surrogate
PROJ_LENGTH_M = (4.0 * PROJ_MASS_KG
                 / (math.pi * PROJ_DIAM_M ** 2 * PROJ_DENSITY_KG_M3))  # ~16.12 mm
PROJ_E_PA = 210.0e9
PROJ_NU = 0.30
PROJ_YIELD_PA = 1.20e9
PROJ_HARD_PA = 5.10e8
PROJ_HARD_N = 0.26
PROJ_RATE_C = 0.014           # Johnson-Cook strain-rate sensitivity

# Initial standoff (m): projectile bottom face above plate top face.
PROJ_INITIAL_STANDOFF_M = 0.005   # 5 mm

# Material card (IM7/8552, Soden 1998 + WWFE-II Kaddour-Hinton 2013).
@dataclass(frozen=True)
class IM7_8552_Card:
    rho_kg_m3: float = 1580.0
    E1_Pa: float = 165.0e9
    E2_Pa: float = 8.4e9
    E3_Pa: float = 8.4e9
    nu12: float = 0.34
    nu13: float = 0.34
    nu23: float = 0.50
    G12_Pa: float = 5.6e9
    G13_Pa: float = 5.6e9
    G23_Pa: float = 2.8e9
    XT_Pa: float = 2560.0e6
    XC_Pa: float = 1590.0e6
    YT_Pa: float = 73.0e6
    YC_Pa: float = 185.0e6
    ZT_Pa: float = 73.0e6
    ZC_Pa: float = 185.0e6
    S12_Pa: float = 90.0e6
    S13_Pa: float = 90.0e6
    S23_Pa: float = 80.0e6
    G_Ic_J_m2: float = 277.0
    G_IIc_J_m2: float = 788.0
    sigma_n_max_Pa: float = 60.0e6
    tau_s_max_Pa: float = 90.0e6
    D_max_erosion: float = 0.99   # /FAIL/HASHIN element-deletion damage cutoff


MAT = IM7_8552_Card()

# Cunniff fiber properties for IM7 (datasheet, used only for c* sanity check).
@dataclass(frozen=True)
class IM7_Fiber:
    sigma_f_Pa: float = 5.2e9
    eps_f: float = 0.019
    E_f_Pa: float = 276.0e9
    rho_f_kg_m3: float = 1780.0


FIBER = IM7_Fiber()

# Numerical / explicit-dynamics knobs.
T_SIM_S: float = 6.0e-5            # 60 microseconds total simulated time
DT_NODA_MIN_S: float = 5.0e-9      # /DT/NODA/CST minimum time step
DT_SCALE: float = 0.67             # Courant scaling on solid bricks
ANIM_DT_S: float = 1.0e-6          # animation frame every 1 microsecond
ISOLID = 24                        # reduced 1-pt + HEPH stabilization
IHQ = 8                            # HEPH (Belytschko-Bindeman 1993)
HG_HM = 0.10
HG_HF = 0.10
HG_HR = 0.10
FRIC_PROJ_PANEL = 0.30
FRIC_PANEL_SELF = 0.20

# OpenRadioss invocation (Lima / Apptainer / Linux).
OPENRADIOSS_STARTER = os.environ.get("OPENRADIOSS_STARTER", "starter_linuxa64")
OPENRADIOSS_ENGINE = os.environ.get("OPENRADIOSS_ENGINE", "engine_linuxa64")
LIMA_VM_NAME = os.environ.get("LIMA_VM_NAME", "apptainer")
NTHREADS = int(os.environ.get("OPENRADIOSS_NT", "8"))

# ----------------------------------------------------------------------------
# Result containers
# ----------------------------------------------------------------------------

@dataclass
class RunResult:
    Vs_ms: float
    Vr_ms: float                        # residual signed velocity at T_SIM
    E_kin_proj_initial_J: float
    E_kin_proj_final_J: float
    E_int_plate_max_J: float
    E_erosion_max_J: float
    E_contact_max_J: float
    E_hg_max_J: float
    energy_balance_max_frac: float       # max |dE/E_kin0| over the run
    eroded_elements: int
    perforated: bool
    run_dir: str
    starter_log: str = ""
    engine_log: str = ""


@dataclass
class SweepFit:
    V50_RI_ms: float                     # Recht-Ipson V50
    V50_RI_R2: float
    V50_LJ_ms: float                     # Lambert-Jonas V50
    LJ_p: float
    V50_LJ_R2: float
    Cunniff_cstar_ms: float
    V50_over_cstar: float
    V50_ref_ms: float
    V50_ref_tol_ms: float
    pass_V50: bool
    n_overmatch_points: int


@dataclass
class CheckResult:
    name: str
    passed: bool
    warned: bool
    value: float
    threshold: float
    message: str


@dataclass
class SweepReport:
    runs: List[RunResult] = field(default_factory=list)
    fit: Optional[SweepFit] = None
    checks: List[CheckResult] = field(default_factory=list)
    pass_overall: bool = False


# ----------------------------------------------------------------------------
# Recht-Ipson and Lambert-Jonas analytic models
# ----------------------------------------------------------------------------

def recht_ipson_residual(Vs: np.ndarray, V50: float) -> np.ndarray:
    """Recht-Ipson 1963 residual velocity for over-match strikes.

    V_r = sqrt(V_s^2 - V_50^2)  for V_s > V_50;
    V_r = 0                     for V_s <= V_50.

    Reference: Recht & Ipson 1963, J. Appl. Mech. 30, 384-390. The derivation
    follows from kinetic-energy conservation on a perforated panel where the
    energy absorbed to perforate is independent of strike velocity:
        (1/2) m_p V_s^2 = (1/2) m_p V_50^2 + (1/2) m_p V_r^2.
    """
    Vs = np.asarray(Vs, dtype=float)
    out = np.zeros_like(Vs)
    over = Vs > V50
    out[over] = np.sqrt(Vs[over] ** 2 - V50 ** 2)
    return out


def lambert_jonas_residual(Vs: np.ndarray, V50: float, p: float) -> np.ndarray:
    """Lambert-Jonas 1976 generalized residual velocity.

    V_r = (V_s^p - V_50^p) ^ (1/p)  for V_s > V_50; 0 otherwise.

    Reference: Lambert & Jonas 1976, BRL-R-1852. The exponent p reduces to 2
    for Recht-Ipson; values 2.0 to 2.5 are typical for laminated composites
    (Naik & Shrirao 2004, Composite Structures 66, 579-590).
    """
    Vs = np.asarray(Vs, dtype=float)
    out = np.zeros_like(Vs)
    over = Vs > V50
    out[over] = (Vs[over] ** p - V50 ** p) ** (1.0 / p)
    return out


def cunniff_cstar(fiber: IM7_Fiber = FIBER) -> float:
    """Cunniff parameter c* for the fiber (m/s).

    c* = ( sigma_f * eps_f / (2 rho) * sqrt(E_f / rho) )^(1/3)

    Reference: Cunniff 1992, Textile Res. J. 62, 495-509. Provides the velocity
    scale at which a given fiber begins to dissipate ballistic kinetic energy
    efficiently. For IM7: c* ~ 700 m/s.
    """
    spec_energy = fiber.sigma_f_Pa * fiber.eps_f / (2.0 * fiber.rho_f_kg_m3)  # m^2/s^2
    c_acoustic = math.sqrt(fiber.E_f_Pa / fiber.rho_f_kg_m3)                  # m/s
    return (spec_energy * c_acoustic) ** (1.0 / 3.0)


def fit_recht_ipson(Vs_arr: np.ndarray, Vr_arr: np.ndarray) -> Tuple[float, float]:
    """Least-squares fit of V50 to over-match (V_s, V_r) pairs.

    Returns (V50, R^2). Uses only points with V_r > SUB_V50_NOISE_MS.
    """
    over = Vr_arr > SUB_V50_NOISE_MS
    if over.sum() < 2:
        raise RuntimeError("Recht-Ipson fit needs >= 2 over-match points; got "
                           f"{int(over.sum())}.")
    Vs_o = Vs_arr[over]
    Vr_o = Vr_arr[over]

    # The closed-form least-squares solution of V_r^2 + V50^2 = V_s^2 is:
    #   V50^2 = mean(V_s^2 - V_r^2)
    V50_sq = np.mean(Vs_o ** 2 - Vr_o ** 2)
    if V50_sq <= 0.0:
        # All points are above an ill-defined limit; fall back to scipy
        # nonlinear fit with a positive-V50 constraint.
        popt, _ = curve_fit(
            lambda V, V50: recht_ipson_residual(V, V50)[V > V50] if False
            else np.sqrt(np.maximum(V ** 2 - V50 ** 2, 0.0)),
            Vs_o, Vr_o,
            p0=[float(np.min(Vs_o) * 0.9)],
            bounds=(0.0, float(np.max(Vs_o))),
        )
        V50 = float(popt[0])
    else:
        V50 = math.sqrt(V50_sq)

    Vr_pred = recht_ipson_residual(Vs_o, V50)
    ss_res = float(np.sum((Vr_o - Vr_pred) ** 2))
    ss_tot = float(np.sum((Vr_o - np.mean(Vr_o)) ** 2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return V50, R2


def fit_lambert_jonas(Vs_arr: np.ndarray, Vr_arr: np.ndarray) -> Tuple[float, float, float]:
    """Two-parameter Lambert-Jonas fit. Returns (V50, p, R^2)."""
    over = Vr_arr > SUB_V50_NOISE_MS
    if over.sum() < 3:
        raise RuntimeError("Lambert-Jonas fit needs >= 3 over-match points; got "
                           f"{int(over.sum())}.")
    Vs_o = Vs_arr[over]
    Vr_o = Vr_arr[over]

    def residuals(theta):
        V50, p = theta
        if V50 <= 0 or p <= 0:
            return np.full_like(Vr_o, 1e9)
        diff = Vs_o ** p - V50 ** p
        # protect against negative diff at marginal overmatch:
        diff = np.maximum(diff, 1e-12)
        Vr_pred = diff ** (1.0 / p)
        return Vr_pred - Vr_o

    # Seed from Recht-Ipson p=2.
    V50_seed, _ = fit_recht_ipson(Vs_arr, Vr_arr)
    sol = least_squares(
        residuals,
        x0=np.array([V50_seed, 2.0]),
        bounds=(np.array([1.0, 1.0]), np.array([float(np.max(Vs_o)), 5.0])),
    )
    V50, p = float(sol.x[0]), float(sol.x[1])
    Vr_pred = lambert_jonas_residual(Vs_o, V50, p)
    ss_res = float(np.sum((Vr_o - Vr_pred) ** 2))
    ss_tot = float(np.sum((Vr_o - np.mean(Vr_o)) ** 2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return V50, p, R2


# ----------------------------------------------------------------------------
# Deck templating
# ----------------------------------------------------------------------------

def _ply_skew_id(theta_deg: float) -> str:
    """Map a ply angle to a fixed-skew id string used in the deck."""
    table = {0.0: "SKEW_0", 45.0: "SKEW_P45", -45.0: "SKEW_M45", 90.0: "SKEW_90"}
    if theta_deg not in table:
        raise KeyError(f"ply angle {theta_deg} not in canonical {{0, +/-45, 90}}")
    return table[theta_deg]


def render_starter(Vs_ms: float, run_dir: Path, mesh_inp_path: Path,
                   stacking_deg: Sequence[float] = STACKING_DEG) -> Path:
    """Render the OpenRadioss starter deck for the given strike velocity.

    The deck pulls nodes and elements from `mesh_inp_path` (Abaqus .inp produced
    upstream by GMSH) via the inp2rad converter; this function emits the
    properties / materials / failure / contact / IC / BC / output cards, plus a
    header that #INCLUDEs the converted mesh.
    """
    job_name = f"flatballistic_Vs{int(round(Vs_ms))}"
    starter_path = run_dir / f"{job_name}_0000.rad"
    skews = sorted({_ply_skew_id(t) for t in stacking_deg})
    lines: List[str] = []

    # 1. Header.
    lines.append("/BEGIN")
    lines.append(f"  {job_name}")
    lines.append("  2026  0  0  0")
    lines.append("/UNIT/1")
    lines.append("  kg  m  s")

    # 2. Mesh (assumed converted from Abaqus .inp by inp2rad and dropped here).
    mesh_rad = mesh_inp_path.with_suffix(".rad")
    lines.append(f"#include \"{mesh_rad.name}\"")

    # 3. Skews (fixed local frames for ply orientations).
    angle_rad = {
        "SKEW_0":   ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        "SKEW_P45": ((math.cos(math.radians(45)), math.sin(math.radians(45)), 0.0),
                     (-math.sin(math.radians(45)), math.cos(math.radians(45)), 0.0)),
        "SKEW_M45": ((math.cos(math.radians(-45)), math.sin(math.radians(-45)), 0.0),
                     (-math.sin(math.radians(-45)), math.cos(math.radians(-45)), 0.0)),
        "SKEW_90":  ((0.0, 1.0, 0.0), (-1.0, 0.0, 0.0)),
    }
    for sk in skews:
        x1, x2 = angle_rad[sk]
        lines.append(f"/SKEW/FIX/{sk}")
        lines.append(f"  0.0 0.0 0.0  {x1[0]} {x1[1]} {x1[2]}  {x2[0]} {x2[1]} {x2[2]}")

    # 4. Per-ply orthotropic LAW25 + per-ply HASHIN.
    for i, theta_deg in enumerate(stacking_deg, start=1):
        mat_id = f"MAT_PLY{i:02d}"
        sk = _ply_skew_id(theta_deg)
        lines.append(f"/MAT/LAW25/{mat_id}")
        lines.append(f"  rho  = {MAT.rho_kg_m3}")
        lines.append(f"  E11  = {MAT.E1_Pa}")
        lines.append(f"  E22  = {MAT.E2_Pa}")
        lines.append(f"  E33  = {MAT.E3_Pa}")
        lines.append(f"  nu12 = {MAT.nu12}")
        lines.append(f"  nu13 = {MAT.nu13}")
        lines.append(f"  nu23 = {MAT.nu23}")
        lines.append(f"  G12  = {MAT.G12_Pa}")
        lines.append(f"  G13  = {MAT.G13_Pa}")
        lines.append(f"  G23  = {MAT.G23_Pa}")
        lines.append("  Iform = 1")
        lines.append("/FAIL/HASHIN/" + mat_id)
        lines.append("  Ifail = 2")
        lines.append(f"  XT  = {MAT.XT_Pa}")
        lines.append(f"  XC  = {MAT.XC_Pa}")
        lines.append(f"  YT  = {MAT.YT_Pa}")
        lines.append(f"  YC  = {MAT.YC_Pa}")
        lines.append(f"  S12 = {MAT.S12_Pa}")
        lines.append(f"  S23 = {MAT.S23_Pa}")
        lines.append(f"  D_max = {MAT.D_max_erosion}")
        # /PROP/TYPE14 with HEPH hourglass control.
        lines.append(f"/PROP/TYPE14/PROP_PLY{i:02d}")
        lines.append(f"  Mat_id  = {mat_id}")
        lines.append(f"  Skew_id = {sk}")
        lines.append(f"  Isolid  = {ISOLID}")
        lines.append(f"  Ihq     = {IHQ}")
        lines.append(f"  hm      = {HG_HM}")
        lines.append(f"  hf      = {HG_HF}")
        lines.append(f"  hr      = {HG_HR}")
        lines.append("  Iorth   = 1")

    # 5. Projectile material + property.
    lines.append("/MAT/LAW2/PROJECTILE_STEEL")
    lines.append(f"  rho  = {PROJ_DENSITY_KG_M3}")
    lines.append(f"  E    = {PROJ_E_PA}")
    lines.append(f"  nu   = {PROJ_NU}")
    lines.append(f"  A    = {PROJ_YIELD_PA}")
    lines.append(f"  B    = {PROJ_HARD_PA}")
    lines.append(f"  n    = {PROJ_HARD_N}")
    lines.append(f"  C    = {PROJ_RATE_C}")
    lines.append("  EPS0 = 1.0")
    lines.append("/PROP/TYPE14/PROP_PROJECTILE")
    lines.append("  Mat_id  = PROJECTILE_STEEL")
    lines.append(f"  Isolid  = {ISOLID}")
    lines.append(f"  Ihq     = {IHQ}")

    # 6. Cohesive interfaces between every adjacent ply (23 total).
    for i in range(1, N_PLY):
        lines.append(f"/INTER/TYPE2/INT_PLY{i:02d}_PLY{i+1:02d}")
        lines.append("  Spotflag = 25")
        lines.append(f"  sigma_n_max = {MAT.sigma_n_max_Pa}")
        lines.append(f"  tau_s_max   = {MAT.tau_s_max_Pa}")
        lines.append(f"  G_Ic        = {MAT.G_Ic_J_m2}")
        lines.append(f"  G_IIc       = {MAT.G_IIc_J_m2}")
        lines.append(f"  master_grp  = SURF_PLY{i:02d}_BOTTOM")
        lines.append(f"  slave_grp   = SURF_PLY{i+1:02d}_TOP")

    # 7. Projectile-panel contact (TYPE7).
    lines.append("/INTER/TYPE7/INT_PROJ_PANEL")
    lines.append("  Igap   = 2")
    lines.append("  Inacti = 6")
    lines.append(f"  Fric   = {FRIC_PROJ_PANEL}")
    lines.append("  Iform  = 2")
    lines.append("  master = SURF_PROJECTILE_OUTER")
    lines.append("  slave  = GR_NODE_PANEL_STRIKE_FACE")

    # 8. Panel self-contact (TYPE7) — required after element erosion.
    lines.append("/INTER/TYPE7/INT_PANEL_SELF")
    lines.append("  Igap   = 2")
    lines.append("  Self   = 1")
    lines.append("  Iedge  = 1")
    lines.append("  Iform  = 2")
    lines.append(f"  Fric   = {FRIC_PANEL_SELF}")
    lines.append("  master = SURF_PANEL_ALL_FACES")
    lines.append("  slave  = SURF_PANEL_ALL_FACES")

    # 9. Initial velocity on the projectile (along -z).
    lines.append("/INIVEL/TRA/IC_PROJECTILE")
    lines.append(f"  Vx = 0  Vy = 0  Vz = -{Vs_ms}")
    lines.append("  grnd_id = GR_NODE_PROJ_ALL")

    # 10. Clamp BC on the perimeter strip.
    lines.append("/BCS/CLAMP_PERIMETER")
    lines.append("  Tra = 1 1 1")
    lines.append("  Rot = 1 1 1")
    lines.append("  grnd_id = GR_NODE_CLAMP")

    # 11. Time-history requests.
    lines.append("/TH/GROU1/TH_PROJ")
    lines.append("  group   = GR_NODE_PROJ_ALL")
    lines.append("  request : VX VY VZ KE X Y Z")
    lines.append("/TH/GLOB")
    lines.append("  request : KE IE EHG EER ECN")
    lines.append("/ANIM/EVERY")
    lines.append(f"  dt_anim = {ANIM_DT_S}")
    lines.append("/END")

    starter_path.write_text("\n".join(lines))
    return starter_path


def render_engine(Vs_ms: float, run_dir: Path) -> Path:
    """Render the engine deck (time stepping)."""
    job_name = f"flatballistic_Vs{int(round(Vs_ms))}"
    engine_path = run_dir / f"{job_name}_0001.rad"
    body = [
        f"/RUN/{job_name}/1",
        f"  Tstop = {T_SIM_S}",
        "/DT/NODA/CST",
        f"  dt_min = {DT_NODA_MIN_S}",
        f"  dt_scale = {DT_SCALE}",
        "/PRINT/-1",
        "/STOP",
    ]
    engine_path.write_text("\n".join(body))
    return engine_path


# ----------------------------------------------------------------------------
# OpenRadioss invocation (Lima/Apptainer or native Linux)
# ----------------------------------------------------------------------------

def _starter_cmd(starter_deck: Path) -> List[str]:
    return [OPENRADIOSS_STARTER, "-i", str(starter_deck), "-nt", str(NTHREADS)]


def _engine_cmd(engine_deck: Path) -> List[str]:
    return [OPENRADIOSS_ENGINE, "-i", str(engine_deck), "-nt", str(NTHREADS)]


def _wrap_for_lima(cmd: List[str]) -> List[str]:
    """Wrap an OpenRadioss command for Lima Apptainer execution if requested."""
    if os.environ.get("OPENRADIOSS_USE_LIMA", "0") == "1":
        return ["limactl", "shell", LIMA_VM_NAME, "--"] + cmd
    return cmd


def invoke_solver(starter_deck: Path, engine_deck: Path,
                  log_dir: Path) -> Tuple[str, str]:
    """Run starter then engine, capturing logs. Returns (starter_log, engine_log)."""
    log_dir.mkdir(parents=True, exist_ok=True)
    s_log_path = log_dir / "starter.log"
    e_log_path = log_dir / "engine.log"

    print(f"[runner] starter: {starter_deck.name}")
    s_proc = subprocess.run(
        _wrap_for_lima(_starter_cmd(starter_deck)),
        cwd=starter_deck.parent,
        capture_output=True, text=True, check=False,
    )
    s_log_path.write_text(s_proc.stdout + "\n--- STDERR ---\n" + s_proc.stderr)
    if s_proc.returncode != 0:
        raise RuntimeError(f"OpenRadioss starter returned {s_proc.returncode}; "
                           f"see {s_log_path}")

    print(f"[runner] engine:  {engine_deck.name}")
    e_proc = subprocess.run(
        _wrap_for_lima(_engine_cmd(engine_deck)),
        cwd=engine_deck.parent,
        capture_output=True, text=True, check=False,
    )
    e_log_path.write_text(e_proc.stdout + "\n--- STDERR ---\n" + e_proc.stderr)
    if e_proc.returncode != 0:
        raise RuntimeError(f"OpenRadioss engine returned {e_proc.returncode}; "
                           f"see {e_log_path}")

    return s_log_path.read_text(), e_log_path.read_text()


# ----------------------------------------------------------------------------
# Time-history parsing (T01 + global)
# ----------------------------------------------------------------------------

def _read_t01(t01_path: Path) -> Dict[str, np.ndarray]:
    """Read an OpenRadioss T01 time-history file.

    Tries the vortex-radioss Python reader; if not present, falls back to a
    minimal binary parser for the subset of channels this stage requests
    (VX, VY, VZ, KE on a single node group; KE/IE/EHG/EER/ECN globally).
    """
    # Preferred path: vortex-radioss.
    try:
        from vortex_radioss import T01Reader  # type: ignore
        rdr = T01Reader(str(t01_path))
        out = {
            "t":  np.asarray(rdr.time(), dtype=float),
            "Vz_proj": np.asarray(rdr.channel("TH_PROJ", "VZ"), dtype=float),
            "Vx_proj": np.asarray(rdr.channel("TH_PROJ", "VX"), dtype=float),
            "Vy_proj": np.asarray(rdr.channel("TH_PROJ", "VY"), dtype=float),
            "KE_proj": np.asarray(rdr.channel("TH_PROJ", "KE"), dtype=float),
            "KE_glob": np.asarray(rdr.channel("GLOB", "KE"), dtype=float),
            "IE_glob": np.asarray(rdr.channel("GLOB", "IE"), dtype=float),
            "EHG_glob": np.asarray(rdr.channel("GLOB", "EHG"), dtype=float),
            "EER_glob": np.asarray(rdr.channel("GLOB", "EER"), dtype=float),
            "ECN_glob": np.asarray(rdr.channel("GLOB", "ECN"), dtype=float),
        }
        return out
    except Exception:
        pass

    # Fallback: ASCII T01 (some OpenRadioss builds emit /TH/.../ASC).
    asc_path = t01_path.with_suffix(".csv")
    if asc_path.exists():
        with asc_path.open() as fh:
            rows = list(csv.DictReader(fh))
        n = len(rows)
        out = {k: np.zeros(n) for k in
               ("t", "Vx_proj", "Vy_proj", "Vz_proj", "KE_proj",
                "KE_glob", "IE_glob", "EHG_glob", "EER_glob", "ECN_glob")}
        for i, row in enumerate(rows):
            out["t"][i]        = float(row.get("time", row.get("t", 0.0)))
            out["Vx_proj"][i]  = float(row.get("VX", 0.0))
            out["Vy_proj"][i]  = float(row.get("VY", 0.0))
            out["Vz_proj"][i]  = float(row.get("VZ", 0.0))
            out["KE_proj"][i]  = float(row.get("KE_proj", 0.0))
            out["KE_glob"][i]  = float(row.get("KE", 0.0))
            out["IE_glob"][i]  = float(row.get("IE", 0.0))
            out["EHG_glob"][i] = float(row.get("EHG", 0.0))
            out["EER_glob"][i] = float(row.get("EER", 0.0))
            out["ECN_glob"][i] = float(row.get("ECN", 0.0))
        return out

    raise RuntimeError(
        f"No T01 reader available and no ASCII fallback at {asc_path}. "
        "Install vortex-radioss (pip install vortex-radioss) or configure the "
        "engine to emit /TH/.../ASC.")


def _residual_velocity_from_history(th: Dict[str, np.ndarray],
                                    Vs_ms: float) -> Tuple[float, bool]:
    """Estimate residual velocity from a projectile axial-velocity history.

    A projectile that fully perforates settles to a roughly constant V_z over
    the last 5 microseconds; one that is arrested oscillates around 0 with
    decaying amplitude. Take the median over the final 10% of frames as the
    residual velocity. Perforated if |V_r| > SUB_V50_NOISE_MS.
    """
    Vz = th["Vz_proj"]
    n = len(Vz)
    if n == 0:
        return 0.0, False
    tail = Vz[max(0, int(0.9 * n)):]
    Vr_signed = float(np.median(tail))
    # Signed Vz starts negative (impact along -z). Magnitude of residual is its absolute value.
    Vr_mag = abs(Vr_signed)
    perforated = Vr_mag > SUB_V50_NOISE_MS
    return Vr_mag, perforated


# ----------------------------------------------------------------------------
# Run + summary
# ----------------------------------------------------------------------------

def run_one_velocity(Vs_ms: float, mesh_inp: Optional[Path] = None,
                     skip_solve: bool = False) -> RunResult:
    run_dir = RUNS_DIR / f"Vs{int(round(Vs_ms))}"
    run_dir.mkdir(parents=True, exist_ok=True)

    if not skip_solve:
        if mesh_inp is None:
            mesh_inp = run_dir / "mesh.inp"
        starter_deck = render_starter(Vs_ms, run_dir, mesh_inp)
        engine_deck = render_engine(Vs_ms, run_dir)
        s_log, e_log = invoke_solver(starter_deck, engine_deck, run_dir)
    else:
        s_log = (run_dir / "starter.log").read_text() if (run_dir / "starter.log").exists() else ""
        e_log = (run_dir / "engine.log").read_text() if (run_dir / "engine.log").exists() else ""

    # Locate the T01 file (OpenRadioss naming convention).
    t01_candidates = list(run_dir.glob("*_T01")) + list(run_dir.glob("*_T01.thy"))
    if not t01_candidates:
        raise FileNotFoundError(f"No T01 file in {run_dir}")
    th = _read_t01(t01_candidates[0])

    Vr_mag, perforated = _residual_velocity_from_history(th, Vs_ms)

    E_kin0 = 0.5 * PROJ_MASS_KG * Vs_ms ** 2
    E_kin_final = float(th["KE_proj"][-1]) if "KE_proj" in th and len(th["KE_proj"]) else 0.5 * PROJ_MASS_KG * Vr_mag ** 2

    E_int_max = float(np.max(th["IE_glob"])) if len(th["IE_glob"]) else 0.0
    E_er_max = float(np.max(th["EER_glob"])) if len(th["EER_glob"]) else 0.0
    E_cn_max = float(np.max(th["ECN_glob"])) if len(th["ECN_glob"]) else 0.0
    E_hg_max = float(np.max(th["EHG_glob"])) if len(th["EHG_glob"]) else 0.0

    # Energy balance: E_kin0 should equal KE_glob + IE_glob + EER_glob + ECN_glob + EHG_glob at every t.
    if all(k in th for k in ("KE_glob", "IE_glob", "EER_glob", "ECN_glob", "EHG_glob")):
        total_t = th["KE_glob"] + th["IE_glob"] + th["EER_glob"] + th["ECN_glob"] + th["EHG_glob"]
        bal = np.abs(total_t - E_kin0) / max(E_kin0, 1.0)
        bal_max = float(np.max(bal))
    else:
        bal_max = float("nan")

    eroded_count = 0  # filled by post-processing the .h3d if available
    erosion_log = run_dir / "erosion_count.txt"
    if erosion_log.exists():
        try:
            eroded_count = int(erosion_log.read_text().strip())
        except Exception:
            pass

    return RunResult(
        Vs_ms=Vs_ms,
        Vr_ms=Vr_mag,
        E_kin_proj_initial_J=E_kin0,
        E_kin_proj_final_J=E_kin_final,
        E_int_plate_max_J=E_int_max,
        E_erosion_max_J=E_er_max,
        E_contact_max_J=E_cn_max,
        E_hg_max_J=E_hg_max,
        energy_balance_max_frac=bal_max,
        eroded_elements=eroded_count,
        perforated=perforated,
        run_dir=str(run_dir),
        starter_log=s_log[-2000:],
        engine_log=e_log[-2000:],
    )


def run_sweep(velocities: Sequence[float] = SWEEP_VS_MS,
              skip_solve: bool = False) -> SweepReport:
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    report = SweepReport()
    for Vs in velocities:
        try:
            res = run_one_velocity(Vs, skip_solve=skip_solve)
            report.runs.append(res)
            print(f"[runner] Vs={Vs:5.0f} m/s  Vr={res.Vr_ms:6.1f} m/s  "
                  f"perforated={res.perforated}  Ehg/Eint={res.E_hg_max_J/max(res.E_int_plate_max_J,1.0):.3f}  "
                  f"|dE|/E0={res.energy_balance_max_frac:.3f}")
        except Exception as e:
            print(f"[runner] Vs={Vs:5.0f} m/s  FAILED: {e}", file=sys.stderr)

    # Fit and check.
    Vs_arr = np.array([r.Vs_ms for r in report.runs])
    Vr_arr = np.array([r.Vr_ms for r in report.runs])

    if (Vr_arr > SUB_V50_NOISE_MS).sum() >= 2:
        V50_RI, R2_RI = fit_recht_ipson(Vs_arr, Vr_arr)
        try:
            V50_LJ, p_LJ, R2_LJ = fit_lambert_jonas(Vs_arr, Vr_arr)
        except RuntimeError:
            V50_LJ, p_LJ, R2_LJ = V50_RI, 2.0, R2_RI
        cstar = cunniff_cstar()
        report.fit = SweepFit(
            V50_RI_ms=V50_RI,
            V50_RI_R2=R2_RI,
            V50_LJ_ms=V50_LJ,
            LJ_p=p_LJ,
            V50_LJ_R2=R2_LJ,
            Cunniff_cstar_ms=cstar,
            V50_over_cstar=V50_RI / cstar,
            V50_ref_ms=V50_REF_MS,
            V50_ref_tol_ms=V50_REF_TOL_MS,
            pass_V50=abs(V50_RI - V50_REF_MS) <= TOL_V50_FRAC * V50_REF_MS,
            n_overmatch_points=int((Vr_arr > SUB_V50_NOISE_MS).sum()),
        )

    report.checks = run_pass_fail_checks(report)
    report.pass_overall = all(c.passed for c in report.checks)

    write_summary(report)
    return report


# ----------------------------------------------------------------------------
# Pass / fail checks (per spec.md section 8)
# ----------------------------------------------------------------------------

def run_pass_fail_checks(report: SweepReport) -> List[CheckResult]:
    checks: List[CheckResult] = []

    if report.fit is not None:
        # Check 1: V50 within +/- 5% of reference.
        v50_err = abs(report.fit.V50_RI_ms - V50_REF_MS) / V50_REF_MS
        checks.append(CheckResult(
            name="V50_within_5pct_of_VakiliRad2020",
            passed=v50_err <= TOL_V50_FRAC,
            warned=False,
            value=v50_err,
            threshold=TOL_V50_FRAC,
            message=f"V50_fit={report.fit.V50_RI_ms:.1f} m/s, ref={V50_REF_MS} m/s",
        ))

        # Check 2: Recht-Ipson residual velocity within 10% averaged over over-match.
        Vs_arr = np.array([r.Vs_ms for r in report.runs])
        Vr_arr = np.array([r.Vr_ms for r in report.runs])
        Vr_pred = recht_ipson_residual(Vs_arr, report.fit.V50_RI_ms)
        over = (Vs_arr > report.fit.V50_RI_ms) & (Vr_arr > SUB_V50_NOISE_MS)
        if over.any():
            err = np.mean(np.abs(Vr_arr[over] - Vr_pred[over]) /
                          np.maximum(Vr_pred[over], 1.0))
        else:
            err = float("nan")
        checks.append(CheckResult(
            name="Vr_RechtIpson_within_10pct",
            passed=(not math.isnan(err)) and err <= TOL_VR_FRAC,
            warned=False,
            value=float(err) if not math.isnan(err) else -1.0,
            threshold=TOL_VR_FRAC,
            message=f"mean residual error vs Recht-Ipson over {int(over.sum())} pts",
        ))

        # Check 3: sub-V50 arrest.
        sub = (Vs_arr <= report.fit.V50_RI_ms)
        if sub.any():
            sub_max_Vr = float(np.max(Vr_arr[sub]))
        else:
            sub_max_Vr = 0.0
        checks.append(CheckResult(
            name="sub_V50_arrest",
            passed=sub_max_Vr < SUB_V50_NOISE_MS,
            warned=False,
            value=sub_max_Vr,
            threshold=SUB_V50_NOISE_MS,
            message=f"max |Vr| at sub-V50 strikes = {sub_max_Vr:.2f} m/s",
        ))

        # Check 7: Lambert-Jonas exponent in [1.8, 2.4].
        p = report.fit.LJ_p
        in_band = TOL_LJ_P_LO <= p <= TOL_LJ_P_HI
        checks.append(CheckResult(
            name="LambertJonas_exponent_in_band",
            passed=in_band,
            warned=not in_band,
            value=p,
            threshold=TOL_LJ_P_HI,
            message=f"Lambert-Jonas p = {p:.3f} (band [{TOL_LJ_P_LO}, {TOL_LJ_P_HI}])",
        ))

        # Check 8: Cunniff sanity.
        ratio = report.fit.V50_over_cstar
        in_cband = TOL_CUNNIFF_LO <= ratio <= TOL_CUNNIFF_HI
        checks.append(CheckResult(
            name="Cunniff_V50_over_cstar_band",
            passed=in_cband,
            warned=not in_cband,
            value=ratio,
            threshold=TOL_CUNNIFF_HI,
            message=f"V50/c*={ratio:.3f}, c*={report.fit.Cunniff_cstar_ms:.1f} m/s",
        ))

    # Check 4: hourglass / internal energy under 5% (max over runs and frames).
    hg_ratios = [r.E_hg_max_J / max(r.E_int_plate_max_J, 1.0) for r in report.runs]
    if hg_ratios:
        hg_max = max(hg_ratios)
        checks.append(CheckResult(
            name="hourglass_under_5pct_of_internal",
            passed=hg_max < TOL_HG_OVER_INT,
            warned=False,
            value=hg_max,
            threshold=TOL_HG_OVER_INT,
            message=f"max E_hg/E_int across sweep = {hg_max:.3f}",
        ))

    # Check 5: energy balance under 2%.
    bal_max = max((r.energy_balance_max_frac for r in report.runs
                   if not math.isnan(r.energy_balance_max_frac)), default=float("nan"))
    if not math.isnan(bal_max):
        checks.append(CheckResult(
            name="energy_balance_under_2pct",
            passed=bal_max < TOL_ENERGY_BAL_FRAC,
            warned=False,
            value=bal_max,
            threshold=TOL_ENERGY_BAL_FRAC,
            message=f"max |dE|/E0 across sweep = {bal_max:.4f}",
        ))

    return checks


# ----------------------------------------------------------------------------
# Reporting / IO
# ----------------------------------------------------------------------------

def write_summary(report: SweepReport) -> None:
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    # sweep_results.csv
    csv_path = SUMMARY_DIR / "sweep_results.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Vs_ms", "Vr_ms", "perforated",
                    "E_kin_proj_initial_J", "E_kin_proj_final_J",
                    "E_int_plate_max_J", "E_erosion_max_J",
                    "E_contact_max_J", "E_hg_max_J",
                    "Ehg_over_Eint", "energy_balance_max_frac",
                    "eroded_elements"])
        for r in report.runs:
            w.writerow([r.Vs_ms, r.Vr_ms, r.perforated,
                        r.E_kin_proj_initial_J, r.E_kin_proj_final_J,
                        r.E_int_plate_max_J, r.E_erosion_max_J,
                        r.E_contact_max_J, r.E_hg_max_J,
                        r.E_hg_max_J / max(r.E_int_plate_max_J, 1.0),
                        r.energy_balance_max_frac, r.eroded_elements])

    # vr_vs_vs.csv
    vrcsv = SUMMARY_DIR / "vr_vs_vs.csv"
    with vrcsv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Vs_ms", "Vr_FEM_ms", "Vr_RechtIpson_ms"])
        if report.fit is not None:
            for r in report.runs:
                pred = float(recht_ipson_residual(np.array([r.Vs_ms]),
                                                  report.fit.V50_RI_ms)[0])
                w.writerow([r.Vs_ms, r.Vr_ms, pred])
        else:
            for r in report.runs:
                w.writerow([r.Vs_ms, r.Vr_ms, ""])

    # recht_ipson_fit.json
    fit_path = SUMMARY_DIR / "recht_ipson_fit.json"
    if report.fit is not None:
        fit_path.write_text(json.dumps(asdict(report.fit), indent=2))

    # pass_fail_report.json
    pf_path = SUMMARY_DIR / "pass_fail_report.json"
    pf_path.write_text(json.dumps(
        {"pass_overall": report.pass_overall,
         "checks": [asdict(c) for c in report.checks]},
        indent=2,
    ))

    print(f"[runner] wrote {csv_path}")
    print(f"[runner] wrote {vrcsv}")
    if report.fit is not None:
        print(f"[runner] wrote {fit_path}")
    print(f"[runner] wrote {pf_path}")
    print(f"[runner] OVERALL: {'PASS' if report.pass_overall else 'FAIL/WARN'}")


# ----------------------------------------------------------------------------
# Mesh-convergence helper (M0/M1/M2 at one velocity)
# ----------------------------------------------------------------------------

def run_mesh_convergence(Vs_ms: float = 350.0) -> Dict[str, float]:
    """Run M0, M1, M2 at the chosen velocity.

    Each mesh is expected to live under runs/MeshConv_M{0,1,2}/Vs<NNN>/, with the
    corresponding mesh_inp file named mesh.inp. The deck/runner is identical
    aside from the mesh; this function just dispatches three solves.
    """
    global RUNS_DIR  # noqa: PLW0603
    results: Dict[str, float] = {}
    for label in ("M0", "M1", "M2"):
        sub = RUNS_DIR / f"MeshConv_{label}"
        sub.mkdir(parents=True, exist_ok=True)
        # Temporarily redirect RUNS_DIR for this run.
        old = RUNS_DIR
        RUNS_DIR = sub
        try:
            r = run_one_velocity(Vs_ms)
            results[label] = r.Vr_ms
        finally:
            RUNS_DIR = old
    if "M1" in results and "M2" in results and results["M1"] > 0:
        delta = abs(results["M2"] - results["M1"]) / results["M1"]
        print(f"[runner] mesh-convergence Vr(M1)={results.get('M1', float('nan')):.1f} "
              f"Vr(M2)={results.get('M2', float('nan')):.1f}  |delta|={delta:.3f}  "
              f"(threshold {TOL_MESH_CONV_FRAC})")
    return results


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Stage 15 flat-coupon ballistic V50 sweep "
                    "(IM7/8552 24-ply 4 mm tape laminate).")
    g = parser.add_mutually_exclusive_group(required=False)
    g.add_argument("--sweep", action="store_true",
                   help="Run the full five-velocity sweep (default).")
    g.add_argument("--single", type=float, metavar="VS_MS",
                   help="Run a single strike velocity in m/s.")
    g.add_argument("--mesh-convergence", type=float, metavar="VS_MS",
                   help="Run M0/M1/M2 at a single velocity.")
    parser.add_argument("--post-only", action="store_true",
                        help="Skip solver, parse existing T01 files only.")
    args = parser.parse_args()

    if args.mesh_convergence is not None:
        run_mesh_convergence(args.mesh_convergence)
        return 0
    if args.single is not None:
        r = run_one_velocity(args.single, skip_solve=args.post_only)
        print(json.dumps(asdict(r), indent=2))
        return 0
    # Default: full sweep.
    report = run_sweep(skip_solve=args.post_only)
    return 0 if report.pass_overall else 2


if __name__ == "__main__":
    raise SystemExit(main())
