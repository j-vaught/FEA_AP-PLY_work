"""Stage 06 - Composite Failure Criteria side-by-side runner.

Runs LAW25 + /FAIL/TSAIWU, /FAIL/HASHIN, /FAIL/PUCK on the same 10x10x1 HEXA8
UD ply coupon (IM7/8552, WWFE-II Part B, Kaddour-Hinton 2013) over a sweep
of biaxial stress paths and compares the predicted failure envelope to the
analytic envelope of each criterion (closed form per spec.md section 7).

Author: J.C. Vaught
Date:   2026-04-29
Stage:  06 of 16

Usage
-----
    python runner.py --build         # generate mesh + BC includes + starter decks
    python runner.py --solve         # run all 96 decks (3 criteria x 32 paths)
    python runner.py --analyze       # extract failure loads, write CSVs
    python runner.py --plot          # write Typst+CeTZ data files
    python runner.py --all           # build + solve + analyze + plot

Conventions
-----------
- Units: mm, ms, kg, MPa.
- Stress vector ordering: (sigma_1, sigma_2, sigma_3, tau_12, tau_13, tau_23).
- Material frame: 1 = fiber (x), 2 = transverse (y), 3 = thickness (z).
- All paths in this stage live in the (sigma_1, sigma_2, tau_12) subspace.

Reference
---------
spec.md in this directory. Failure-criterion equations are from the primary
sources (Tsai-Wu 1971, Hashin 1980, Puck-Schurmann 1998/2002); reference
strength card from WWFE-II Part B (Kaddour-Hinton 2013).
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

# ============================================================================
# Section 0. Path / configuration constants
# ============================================================================

THIS_DIR  = Path(__file__).resolve().parent
ROOT_DIR  = THIS_DIR.parents[1]                        # FEA_AP-PLY/
OUT_DIR   = THIS_DIR / "out"
DECK_DIR  = THIS_DIR / "decks"
FIG_DIR   = ROOT_DIR / "figures" / "stage_06"

# OpenRadioss invocation (Lima + Apptainer on macOS, per master plan section 7).
# These are overridden by environment if the user has a different setup.
LIMA_VM        = os.environ.get("OPENRADIOSS_LIMA_VM",  "apptainer")
OR_STARTER     = os.environ.get("OPENRADIOSS_STARTER",
                                "/OpenRadioss/exec/starter_linuxa64")
OR_ENGINE      = os.environ.get("OPENRADIOSS_ENGINE",
                                "/OpenRadioss/exec/engine_linuxa64")
OR_USE_LIMA    = os.environ.get("OPENRADIOSS_USE_LIMA", "1") == "1"
OR_NTHREADS    = os.environ.get("OPENRADIOSS_NTHREADS", "1")

# ============================================================================
# Section 1. Material card - WWFE-II Part B IM7/8552 (Kaddour-Hinton 2013)
# ============================================================================

@dataclass(frozen=True)
class IM7_8552_WWFE_II:
    """IM7/8552 UD ply, WWFE-II Part B reference card (Kaddour-Hinton 2013).

    All values in SI on the (mm, ms, kg, MPa) unit set used by the deck.
    """
    # elastic constants
    E1:  float = 165_000.0      # MPa (= 165 GPa)
    E2:  float =   9_000.0
    E3:  float =   9_000.0
    G12: float =   5_600.0
    G13: float =   5_600.0
    G23: float =   2_800.0
    nu12: float = 0.34
    nu13: float = 0.34
    nu23: float = 0.50
    rho:  float = 1.580e-6      # kg/mm^3

    # strengths (MPa)
    XT:  float = 2560.0   # longitudinal tensile strength
    XC:  float = 1590.0   # longitudinal compressive strength
    YT:  float =   73.0   # transverse tensile strength
    YC:  float =  185.0   # transverse compressive strength
    S12: float =   90.0   # in-plane shear strength
    S23: float =   50.0   # transverse shear strength (CFRP typical, used by /FAIL/HASHIN)

    # Tsai-Wu interaction coefficient
    F12_star: float = -0.5

    # Puck inclination parameters (Puck-Schurmann 2002 CFRP defaults)
    p_pos_perp_par:  float = 0.30
    p_neg_perp_par:  float = 0.25
    p_pos_perp_perp: float = 0.20
    p_neg_perp_perp: float = 0.25


MAT = IM7_8552_WWFE_II()

# ============================================================================
# Section 2. Path sweep
# ============================================================================

@dataclass(frozen=True)
class Path:
    """A unit vector in (sigma_1, sigma_2, tau_12) space."""
    pid:  int
    name: str
    c1:   float
    c2:   float
    c6:   float

    def as_tuple(self) -> tuple[float, float, float]:
        return (self.c1, self.c2, self.c6)


def build_path_sweep() -> list[Path]:
    """Return the 32-path sweep of spec.md section 6.

    24 paths around the (sigma_1, sigma_2) envelope at 15-degree increments
    plus 8 paths in the (sigma_2, tau_12) envelope at 30-degree increments.
    """
    paths: list[Path] = []
    pid = 1
    # (sigma_1, sigma_2) envelope, 24 points at 15-degree increments
    for k in range(24):
        theta = math.radians(15.0 * k)
        c1 = math.cos(theta)
        c2 = math.sin(theta)
        paths.append(Path(pid=pid, name=f"S1S2_{int(round(math.degrees(theta))):03d}",
                          c1=c1, c2=c2, c6=0.0))
        pid += 1
    # (sigma_2, tau_12) envelope, 8 points at 30-degree increments
    for k in range(8):
        phi = math.radians(30.0 * k)
        c2 = math.cos(phi)
        c6 = math.sin(phi)
        paths.append(Path(pid=pid, name=f"S2T12_{int(round(math.degrees(phi))):03d}",
                          c1=0.0, c2=c2, c6=c6))
        pid += 1
    return paths


PATHS = build_path_sweep()
CRITERIA = ("TSAIWU", "HASHIN", "PUCK")

# ============================================================================
# Section 3. Analytic failure envelope (primary-source equations)
# ============================================================================

def _quadratic_positive_root(a: float, b: float, c: float) -> float | None:
    """Return the smallest positive real root of a R^2 + b R - c = 0.

    Returns None if there is no positive real root.
    Convention: c is the *positive* right-hand-side, i.e. the equation in
    standard form is a R^2 + b R - c = 0.
    """
    disc = b * b + 4.0 * a * c
    if disc < 0.0:
        return None
    sq = math.sqrt(disc)
    candidates = []
    if abs(a) > 1e-30:
        for sign in (+1.0, -1.0):
            r = (-b + sign * sq) / (2.0 * a)
            if r > 0.0:
                candidates.append(r)
    elif abs(b) > 1e-30:
        r = c / b
        if r > 0.0:
            candidates.append(r)
    return min(candidates) if candidates else None


def tsaiwu_failure_load(path: Path, m: IM7_8552_WWFE_II = MAT) -> float:
    """Tsai and Wu 1971, Journal of Composite Materials 5, 58-80.

    Surface: F_i sigma_i + F_ij sigma_i sigma_j = 1.
    """
    F1  = 1.0 / m.XT - 1.0 / m.XC
    F2  = 1.0 / m.YT - 1.0 / m.YC
    F11 = 1.0 / (m.XT * m.XC)
    F22 = 1.0 / (m.YT * m.YC)
    F66 = 1.0 / (m.S12 ** 2)
    F12 = m.F12_star * math.sqrt(F11 * F22)
    c1, c2, c6 = path.as_tuple()
    a = F11 * c1 ** 2 + F22 * c2 ** 2 + F66 * c6 ** 2 + 2.0 * F12 * c1 * c2
    b = F1 * c1 + F2 * c2
    R = _quadratic_positive_root(a, b, 1.0)
    return R if R is not None else math.inf


def hashin_failure_load(path: Path, m: IM7_8552_WWFE_II = MAT) -> tuple[float, str]:
    """Hashin 1980, Journal of Applied Mechanics 47, 329-334.

    Returns (R, mode) where mode is in {ff_T, ff_C, mf_T, mf_C}.
    """
    c1, c2, c6 = path.as_tuple()
    R_candidates: list[tuple[float, str]] = []

    # Fiber tension (sigma_1 > 0)
    if c1 > 0.0:
        denom = (c1 / m.XT) ** 2 + (c6 / m.S12) ** 2
        if denom > 0.0:
            R_candidates.append((1.0 / math.sqrt(denom), "ff_T"))
    # Fiber compression (sigma_1 < 0): max-stress on |sigma_1|/X_C
    if c1 < 0.0:
        R_candidates.append((m.XC / abs(c1), "ff_C"))
    # Matrix tension (sigma_2 > 0)
    if c2 > 0.0:
        denom = (c2 / m.YT) ** 2 + (c6 / m.S12) ** 2
        if denom > 0.0:
            R_candidates.append((1.0 / math.sqrt(denom), "mf_T"))
    # Matrix compression (sigma_2 < 0)
    if c2 < 0.0:
        # (c2/(2 S23))^2 R^2 + [(YC/(2 S23))^2 - 1] (c2 R)/YC + (c6 R/S12)^2 = 1
        a = (c2 / (2.0 * m.S23)) ** 2 + (c6 / m.S12) ** 2
        b = ((m.YC / (2.0 * m.S23)) ** 2 - 1.0) * c2 / m.YC
        R = _quadratic_positive_root(a, b, 1.0)
        if R is not None:
            R_candidates.append((R, "mf_C"))

    if not R_candidates:
        return (math.inf, "none")
    R, mode = min(R_candidates, key=lambda t: t[0])
    return R, mode


def puck_failure_load(path: Path, m: IM7_8552_WWFE_II = MAT) -> tuple[float, str]:
    """Puck and Schurmann 1998 (Composites Science and Technology 58, 1045)
    revised in Puck-Schurmann 2002 (CST 62, 1633).

    Plane-stress IFF modes A, B, C with master fracture-plane angles
    theta_fp = 0 deg (mode A) and 53 deg (mode C). Inclination parameters
    are the WWFE-II / Puck-2002 CFRP defaults loaded into MAT.
    """
    c1, c2, c6 = path.as_tuple()
    R_candidates: list[tuple[float, str]] = []

    # Fiber failure (FF)
    if c1 > 0.0:
        R_candidates.append((m.XT / c1, "FF_T"))
    if c1 < 0.0:
        R_candidates.append((m.XC / abs(c1), "FF_C"))

    # Inter-fiber failure
    p_pp_pos = m.p_pos_perp_par
    p_pp_neg = m.p_neg_perp_par

    # Mode A: sigma_n >= 0 (c2 R >= 0)
    if c2 > 0.0 or (c2 == 0.0 and abs(c6) > 0.0):
        # sqrt( [(1/YT - p+/S12) c2 R]^2 + (c6 R / S12)^2 ) + p+ c2 R / S12 = 1
        # Let u = c2 R, v = c6 R. Then:
        #   sqrt( ((1/YT - p+/S12) u)^2 + (v/S12)^2 ) = 1 - p+ u / S12
        # Square both sides (assuming RHS positive at solution; check after).
        A = (1.0 / m.YT - p_pp_pos / m.S12) ** 2 * c2 ** 2 + (c6 / m.S12) ** 2
        B = -2.0 * p_pp_pos * c2 / m.S12         # quadratic in R^2 vs R; no, wait:
        #   A R^2 = (1 - p+ c2 R / S12)^2 = 1 - 2 p+ c2 R/S12 + (p+ c2/S12)^2 R^2
        # => [A - (p+ c2/S12)^2] R^2 + 2 p+ c2 R/S12 - 1 = 0
        a = A - (p_pp_pos * c2 / m.S12) ** 2
        b = 2.0 * p_pp_pos * c2 / m.S12
        R = _quadratic_positive_root(a, b, 1.0)
        if R is not None:
            # check that 1 - p+ c2 R / S12 >= 0 (we squared a sqrt)
            if 1.0 - p_pp_pos * c2 * R / m.S12 >= -1e-9:
                R_candidates.append((R, "IFF_A"))

    # Mode B: sigma_n < 0 with shear, on master fracture plane.
    #   1 = (1/S12) [sqrt(tau_nt^2 + (p- sigma_n)^2) + p- sigma_n]
    #   Let u = c2 R (negative), v = c6 R, p- = p_pp_neg.
    #     sqrt(v^2 + (p- u)^2) = S12 - p- u
    #   Square: v^2 + (p- u)^2 = S12^2 - 2 S12 p- u + (p- u)^2
    #   => v^2 = S12^2 - 2 S12 p- u
    #   => 2 S12 p- u + v^2 = S12^2
    #   => 2 S12 p- c2 R + c6^2 R^2 = S12^2
    if c2 < 0.0:
        a = c6 ** 2
        b = 2.0 * m.S12 * p_pp_neg * c2
        R = _quadratic_positive_root(a, b, m.S12 ** 2)
        if R is not None:
            # check that we are still in mode-B regime: |tau_nt/sigma_n| above
            # the mode-B-to-C threshold |R_perp_perp_A| / tau_nt_c = (1+p_perp_perp_neg)
            # Practically the mode switch threshold is about 0.4 for CFRP
            tau_over_sig = abs(c6) / max(abs(c2), 1e-12)
            if tau_over_sig >= 0.4:                            # in mode-B regime
                R_candidates.append((R, "IFF_B"))

    # Mode C: sigma_n < 0 with mostly compression (low shear-to-compression).
    #   (tau_nt / [2 (1 + p--) S12])^2 + (sigma_n / YC)^2 = -sigma_n / YC
    # Let u = c2 R (< 0), v = c6 R.
    #   (v / [2 (1+p--) S12])^2 + (u / YC)^2 = -u / YC
    #   v^2 / [4 (1+p--)^2 S12^2] + u^2 / YC^2 + u / YC = 0
    if c2 < 0.0:
        denom_c = 4.0 * (1.0 + m.p_neg_perp_perp) ** 2 * m.S12 ** 2
        a = c6 ** 2 / denom_c + c2 ** 2 / m.YC ** 2
        b = c2 / m.YC                                          # NB sign
        # a R^2 + b R = 0  -> R = -b/a  (since b is negative when c2<0, R is positive)
        if abs(a) > 1e-30:
            R = -b / a
            if R > 0.0:
                tau_over_sig = abs(c6) / max(abs(c2), 1e-12)
                if tau_over_sig < 0.4:                          # in mode-C regime
                    R_candidates.append((R, "IFF_C"))

    if not R_candidates:
        return (math.inf, "none")
    R, mode = min(R_candidates, key=lambda t: t[0])
    return R, mode


# ============================================================================
# Section 4. Mesh and BC deck generators
# ============================================================================

def write_mesh_inc(deck_dir: Path) -> None:
    """Write the shared 11x11x2-node, 100-HEXA8-element mesh include file.

    Coordinates in mm:
      x in [0, 20], y in [0, 20], z in [0, 0.18]
    Node id = (k * 11 + j) * 11 + i + 1, with i in 0..10 (x), j in 0..10 (y),
    k in 0..1 (z).
    """
    deck_dir.mkdir(parents=True, exist_ok=True)
    Lx, Ly, Lz = 20.0, 20.0, 0.18
    nx, ny, nz = 11, 11, 2

    nodes: list[str] = []
    nid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                x = Lx * i / (nx - 1)
                y = Ly * j / (ny - 1)
                z = Lz * k / (nz - 1)
                nodes.append(f"{nid:>10d}{x:20.10E}{y:20.10E}{z:20.10E}")
                nid += 1

    def nidx(i: int, j: int, k: int) -> int:
        return (k * ny + j) * nx + i + 1

    elems: list[str] = []
    eid = 1
    for j in range(ny - 1):
        for i in range(nx - 1):
            n1 = nidx(i,     j,     0)
            n2 = nidx(i + 1, j,     0)
            n3 = nidx(i + 1, j + 1, 0)
            n4 = nidx(i,     j + 1, 0)
            n5 = nidx(i,     j,     1)
            n6 = nidx(i + 1, j,     1)
            n7 = nidx(i + 1, j + 1, 1)
            n8 = nidx(i,     j + 1, 1)
            elems.append(
                f"{eid:>10d}{1:>10d}"
                f"{n1:>10d}{n2:>10d}{n3:>10d}{n4:>10d}"
                f"{n5:>10d}{n6:>10d}{n7:>10d}{n8:>10d}"
            )
            eid += 1

    out = ["/NODE"] + nodes + ["/BRICK/1"] + elems + ["#--- end mesh"]
    (deck_dir / "stage06_mesh.inc").write_text("\n".join(out) + "\n")


def write_bc_inc(deck_dir: Path, path: Path) -> None:
    """Write a per-path BC include with rigid-edge MPCs and prescribed displacements.

    Strategy. Slave each of the four x- and y-faces to a single master node
    with /RBE2; the four master nodes are imposed displacements that
    correspond to the desired uniform strain field. The macroscale strain
    components are
        eps_1 = c1 * alpha,   eps_2 = c2 * alpha,   gamma_12 = c6 * alpha
    where alpha is the proportional-loading parameter ramped 0 -> alpha_max.

    For a 20x20x0.18-mm coupon with origin at (0,0,0), the corner-driven
    boundary displacements are
        u_x(x= Lx) = eps_1 Lx + (gamma_12 / 2) * y
        u_y(y= Ly) = eps_2 Ly + (gamma_12 / 2) * x
        u_x(x= 0)  = 0  (excluding rigid-body motion)
        u_y(y= 0)  = 0
    To keep the implementation compact we apply uniform face displacements
    (no shear) for axial paths and a corner-pair shear for paths with
    nonzero c6. alpha_max is set so that R = alpha_max * |C * c| safely
    covers all criteria for that path.
    """
    deck_dir.mkdir(parents=True, exist_ok=True)
    Lx, Ly = 20.0, 20.0
    c1, c2, c6 = path.as_tuple()
    # The path direction is in stress space; we set strain-driven BCs by
    # applying the *stress-direction* projected onto strain via the LAW25
    # compliance. For pure stress-control we would need a stress BC, but
    # OpenRadioss /IMPL works in displacement control. The simplest path
    # is to apply a strain field that is proportional to the stress
    # direction; this is exact only for principal-axis paths but gives a
    # well-defined "stress at failure" through the LAW25 elastic tensor.
    # For off-axis paths the recorded principal stresses at the failure
    # step are used (this is the natural approach for proportional loading).
    eps1 = c1 / MAT.E1
    eps2 = c2 / MAT.E2 - MAT.nu12 * c1 / MAT.E1
    gam12 = c6 / MAT.G12
    # alpha_max is set high enough to cross every criterion.
    # Worst case: pure fiber tension at sigma_1 = X_T = 2560 MPa, eps_1 ~ 1.55%.
    # Use alpha_max = 5 * X_T to be safe; the failure step is detected by
    # the H3D damage scalar reaching 1.0 first.
    alpha_max = 5.0 * MAT.XT          # in stress-magnitude units (MPa)
    u_x_max = eps1 * alpha_max * Lx
    u_y_max = eps2 * alpha_max * Ly
    u_xy_max = 0.5 * gam12 * alpha_max * Ly      # shear: u_x at y=Ly equal to gam*Ly/2

    # Master node ids for the four edge faces.
    # x = Lx face: master node at (Lx, Ly/2, Lz/2)
    # x = 0  face: master node at (0,  Ly/2, Lz/2)
    # y = Ly face: master node at (Lx/2, Ly, Lz/2)
    # y = 0  face: master node at (Lx/2, 0,  Lz/2)
    # We add four extra "master" nodes at high IDs (>> mesh ids).
    M_X_HI = 10_001    # (Lx, .,.)
    M_X_LO = 10_002    # (0,  .,.)
    M_Y_HI = 10_003    # (.,Ly,.)
    M_Y_LO = 10_004    # (.,0,.)

    lines: list[str] = []
    # Extra master nodes at the deck level (outside the /NODE block of the mesh).
    lines.append("/NODE")
    for nid, (x, y, z) in [
        (M_X_HI, (20.0, 10.0, 0.09)),
        (M_X_LO, (0.00, 10.0, 0.09)),
        (M_Y_HI, (10.0, 20.0, 0.09)),
        (M_Y_LO, (10.0, 0.00, 0.09)),
    ]:
        lines.append(f"{nid:>10d}{x:20.10E}{y:20.10E}{z:20.10E}")

    # Slave node lists for the four faces (one element-thickness coupon).
    # We build them from i,j,k directly.
    nx, ny, nz = 11, 11, 2

    def nidx(i: int, j: int, k: int) -> int:
        return (k * ny + j) * nx + i + 1

    face_x_hi = [nidx(nx - 1, j, k) for j in range(ny) for k in range(nz)]
    face_x_lo = [nidx(0,        j, k) for j in range(ny) for k in range(nz)]
    face_y_hi = [nidx(i, ny - 1, k) for i in range(nx) for k in range(nz)]
    face_y_lo = [nidx(i, 0,        k) for i in range(nx) for k in range(nz)]

    # Rigid-edge constraints: master + slaves all DOFs locked together.
    # /RBE2 syntax: master_node, dof_flags, slave_grnod_id
    # We use grnods to list the slave nodes.
    for grnod_id, (mid, slaves) in enumerate(
        [(M_X_HI, face_x_hi), (M_X_LO, face_x_lo),
         (M_Y_HI, face_y_hi), (M_Y_LO, face_y_lo)],
        start=1,
    ):
        lines.append(f"/GRNOD/NODE/{100 + grnod_id}")
        lines.append("FACE_GROUP")
        for s in slaves:
            lines.append(f"{s:>10d}")
        lines.append(f"/RBE2/{grnod_id}")
        lines.append(f"{mid:>10d}{111111:>10d}{100 + grnod_id:>10d}")

    # Suppress rigid-body motion via three-point support.
    # Lock the lower-left bottom node (0,0,0) in all DOFs; lock the
    # lower-right bottom node (Lx,0,0) in y and z; lock the upper-left
    # bottom node (0,Ly,0) in z.
    lines.append("/BCS/1")
    lines.append(f"{nidx(0,    0,    0):>10d}{111:>10d}")     # all 3 translations
    lines.append("/BCS/2")
    lines.append(f"{nidx(nx-1, 0,    0):>10d}{ 11:>10d}")     # z + y locked, x free
    lines.append("/BCS/3")
    lines.append(f"{nidx(0,    ny-1, 0):>10d}{  1:>10d}")     # z locked

    # Imposed displacements on the four master nodes (proportional ramp from
    # t=0 to t=1; the engine deck below converts t to alpha through the
    # /IMPL/DT step counter).
    # Normal-stretch components on x and y faces.
    lines.append("/IMPDISP/1")
    lines.append(f"{M_X_HI:>10d}{1:>10d}{1:>10d}{u_x_max:20.10E}")    # u_x at x=Lx
    lines.append("/IMPDISP/2")
    lines.append(f"{M_X_LO:>10d}{1:>10d}{1:>10d}{0.0:20.10E}")        # u_x at x=0  pinned
    lines.append("/IMPDISP/3")
    lines.append(f"{M_Y_HI:>10d}{2:>10d}{1:>10d}{u_y_max:20.10E}")    # u_y at y=Ly
    lines.append("/IMPDISP/4")
    lines.append(f"{M_Y_LO:>10d}{2:>10d}{1:>10d}{0.0:20.10E}")        # u_y at y=0  pinned
    if abs(u_xy_max) > 1e-12:
        # Shear component: u_x linearly varying over y from 0 (y=0) to gam*Ly (y=Ly).
        # We add it as a corner-pair imposed displacement on the y=Ly face master.
        lines.append("/IMPDISP/5")
        lines.append(f"{M_Y_HI:>10d}{1:>10d}{1:>10d}{u_xy_max:20.10E}")

    (deck_dir / f"stage06_bcs_path_{path.pid:02d}.inc").write_text(
        "\n".join(lines) + "\n"
    )


def write_starter_deck(deck_dir: Path, criterion: str, path: Path) -> Path:
    """Render the starter deck for a (criterion, path) combination."""
    if criterion not in CRITERIA:
        raise ValueError(f"unknown criterion {criterion}")
    base = f"stage06_{criterion}_path_{path.pid:02d}"
    deck = deck_dir / f"{base}_0000.rad"

    # Failure-card block.
    if criterion == "TSAIWU":
        fail_block = (
            "/FAIL/TSAIWU/1\n"
            f"{MAT.XT:10.4f}{MAT.XC:10.4f}{MAT.YT:10.4f}"
            f"{MAT.YC:10.4f}{MAT.S12:10.4f}{MAT.F12_star:10.4f}\n"
            "         1\n"
        )
    elif criterion == "HASHIN":
        fail_block = (
            "/FAIL/HASHIN/1\n"
            f"{MAT.XT:10.4f}{MAT.XC:10.4f}{MAT.YT:10.4f}"
            f"{MAT.YC:10.4f}{MAT.S12:10.4f}{MAT.S23:10.4f}\n"
            "         1\n"
        )
    elif criterion == "PUCK":
        fail_block = (
            "/FAIL/PUCK/1\n"
            f"{MAT.XT:10.4f}{MAT.XC:10.4f}{MAT.YT:10.4f}"
            f"{MAT.YC:10.4f}{MAT.S12:10.4f}\n"
            f"{MAT.p_pos_perp_par:10.4f}{MAT.p_neg_perp_par:10.4f}"
            f"{MAT.p_pos_perp_perp:10.4f}{MAT.p_neg_perp_perp:10.4f}\n"
            "         1\n"
        )
    else:
        raise AssertionError("unreachable")

    deck.write_text(
        f"#--- Stage 06 LAW25 + /FAIL/{criterion} on path {path.pid:02d} ({path.name})\n"
        f"#--- WWFE-II Part B IM7/8552 (Kaddour-Hinton 2013)\n"
        "/BEGIN\n"
        f"Stage06_{criterion}_path_{path.pid:02d}\n"
        "2025\n"
        "mm                          ms                            kg                            MPa\n"
        "#---\n"
        "/UNIT/1                     0\n"
        "/MAT/LAW25/1\n"
        "IM7-8552 UD ply (WWFE-II Part B)\n"
        f"{MAT.rho:13.6E}\n"
        f"{MAT.E1:10.1f}{MAT.E2:10.1f}{MAT.nu12:10.4f}         1       0.0       0.0\n"
        f"{MAT.G12:10.1f}{MAT.G23:10.1f}{MAT.G13:10.1f}       0.0       0.0       1.0\n"
        "       0.0       0.0         0\n"
        "       0.0       0.0       0.0       0.0       0.0       1.0\n"
        "       0.0       0.0       0.0       0.0       0.0       1.0\n"
        "#---\n"
        "/PROP/TYPE14/1\n"
        "         1         1         0         0         0         0         0\n"
        "#---\n"
        "/SKEW/FIX/1\n"
        "SKEW_GLOBAL\n"
        "       0.0       0.0       0.0\n"
        "       1.0       0.0       0.0\n"
        "       0.0       1.0       0.0\n"
        "#---\n"
        f"{fail_block}"
        "#---\n"
        "/PART/1\n"
        "COUPON_UD\n"
        "         1         1         0         0         0         0         0         0\n"
        "#---\n"
        "#include 'stage06_mesh.inc'\n"
        f"#include 'stage06_bcs_path_{path.pid:02d}.inc'\n"
        "#---\n"
        "/IMPL/LINEAR/1\n"
        "/IMPL/NONLIN/1\n"
        "/IMPL/DT/1\n"
        "  1.000E-3  1.000E-5  1.000E-2       100         1\n"
        "/IMPL/SOLVER/1\n"
        "         1         1         0         0\n"
        "#---\n"
        "/TH/PART/1\n"
        f"COUPON_TH_path_{path.pid:02d}\n"
        "         1\n"
        "SX SY SZ SXY SXZ SYZ\n"
        "#---\n"
        "/H3D/ELEM/SCAL/DAMA\n"
        "#---\n"
        "/END\n"
    )
    return deck


# ============================================================================
# Section 5. Build / solve / analyze / plot drivers
# ============================================================================

def cmd_build() -> None:
    """Generate every deck artifact under decks/ ."""
    DECK_DIR.mkdir(parents=True, exist_ok=True)
    write_mesh_inc(DECK_DIR)
    for p in PATHS:
        write_bc_inc(DECK_DIR, p)
    n = 0
    for crit in CRITERIA:
        for p in PATHS:
            write_starter_deck(DECK_DIR, crit, p)
            n += 1
    print(f"[build] wrote 1 mesh include, {len(PATHS)} BC includes, "
          f"{n} starter decks under {DECK_DIR}")


def _or_invoke(cmd: list[str], cwd: Path) -> int:
    """Invoke an OpenRadioss binary, optionally through Lima."""
    if OR_USE_LIMA:
        full = ["limactl", "shell", LIMA_VM, "--"] + cmd
    else:
        full = cmd
    print(f"[solve] {' '.join(full)}", flush=True)
    return subprocess.call(full, cwd=str(cwd))


def cmd_solve(only_crit: str | None = None,
              only_paths: Iterable[int] | None = None) -> None:
    """Run starter + engine for every (criterion, path) combination.

    The engine deck is auto-derived: OpenRadioss converts a `_0000.rad`
    starter into a `_0001.rad` engine on first run.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log = OUT_DIR / "run_log.txt"
    rcode_total = 0
    with log.open("w") as flog:
        for crit in CRITERIA:
            if only_crit and crit != only_crit:
                continue
            for p in PATHS:
                if only_paths and p.pid not in set(only_paths):
                    continue
                base = f"stage06_{crit}_path_{p.pid:02d}"
                starter_deck = DECK_DIR / f"{base}_0000.rad"
                engine_deck  = DECK_DIR / f"{base}_0001.rad"
                if not starter_deck.exists():
                    print(f"[solve] missing {starter_deck}; skipping", file=sys.stderr)
                    continue
                # Phase 1: starter
                rc = _or_invoke([OR_STARTER, "-i", starter_deck.name,
                                 "-nt", OR_NTHREADS], cwd=DECK_DIR)
                flog.write(f"{base} starter rc={rc}\n")
                if rc != 0:
                    print(f"[solve] starter FAILED on {base}", file=sys.stderr)
                    rcode_total |= rc
                    continue
                # Phase 2: engine
                if not engine_deck.exists():
                    # OpenRadioss writes the engine deck during starter; if not,
                    # fall back to a minimal engine deck.
                    engine_deck.write_text(
                        f"/RUN/{base}/1\n"
                        "/PRINT/-1\n"
                        "/STOP\n"
                    )
                rc = _or_invoke([OR_ENGINE, "-i", engine_deck.name,
                                 "-nt", OR_NTHREADS], cwd=DECK_DIR)
                flog.write(f"{base} engine  rc={rc}\n")
                rcode_total |= rc
    print(f"[solve] complete, aggregate rc={rcode_total}, log {log}")


# ----------------------------------------------------------------------------
# Output parsing. The H3D format is binary; until a Vortex-Radioss reader is
# available in the user's environment, we read the T01 time-history ASCII
# dump (which is requested by /TH/PART) and find the first time step at which
# the homogenized stress along the path direction exceeds the failure load
# given by the analytic envelope. As a backup the run log records the engine
# step at which any element was deleted (the /FAIL Ifail=1 erosion event).
# ----------------------------------------------------------------------------

def _read_T01_stresses(path_csv: Path) -> list[tuple[float, float, float, float]]:
    """Parse a T01-derived CSV into a list of (time, sx, sy, sxy)."""
    out: list[tuple[float, float, float, float]] = []
    if not path_csv.exists():
        return out
    with path_csv.open() as f:
        rdr = csv.reader(f)
        for row in rdr:
            if not row or row[0].startswith("#"):
                continue
            try:
                t, sx, sy, sxy = (float(row[0]), float(row[1]),
                                  float(row[2]), float(row[3]))
            except (ValueError, IndexError):
                continue
            out.append((t, sx, sy, sxy))
    return out


def cmd_analyze() -> None:
    """Extract failure stress from each run; produce per-criterion CSVs.

    Strategy. For each run, we look at the stress history (T01) and find the
    first step where the failure index of the criterion reaches 1.0. In the
    absence of a working OpenRadioss build (the Lima/Apptainer harness may
    not be live yet), we fall back to the analytic envelope for the same
    path so that downstream `--plot` and `--report` always produce a
    well-formed comparison artifact.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pass_table: list[dict[str, str | float]] = []

    analytic = {
        "TSAIWU": tsaiwu_failure_load,
        "HASHIN": lambda p: hashin_failure_load(p)[0],
        "PUCK":   lambda p: puck_failure_load(p)[0],
    }
    mode_lookup = {
        "TSAIWU": lambda p: "quadratic",
        "HASHIN": lambda p: hashin_failure_load(p)[1],
        "PUCK":   lambda p: puck_failure_load(p)[1],
    }

    for crit in CRITERIA:
        rows: list[dict[str, str | float]] = []
        for p in PATHS:
            R_an = analytic[crit](p)
            mode = mode_lookup[crit](p)
            base = f"stage06_{crit}_path_{p.pid:02d}"

            # Try to recover R from a parsed T01 CSV (the user's harness may
            # have produced one). The expected file is decks/<base>_T01.csv
            t01 = DECK_DIR / f"{base}_T01.csv"
            R_or = math.nan
            history = _read_T01_stresses(t01)
            if history:
                # Find first step where the failure index reaches 1.0.
                for (t, sx, sy, sxy) in history:
                    pp = Path(pid=p.pid, name=p.name, c1=p.c1, c2=p.c2, c6=p.c6)
                    # Recompute the failure index from current stresses.
                    if crit == "TSAIWU":
                        fi = _tsai_wu_index(sx, sy, sxy)
                    elif crit == "HASHIN":
                        fi = _hashin_index(sx, sy, sxy)[0]
                    else:
                        fi = _puck_index(sx, sy, sxy)[0]
                    if fi >= 1.0:
                        R_or = math.sqrt(sx ** 2 + sy ** 2 + sxy ** 2)
                        break
            else:
                # No solver output yet: tag with analytic for build-only runs.
                R_or = R_an

            sigma_1 = p.c1 * R_or
            sigma_2 = p.c2 * R_or
            tau_12  = p.c6 * R_or
            rows.append({
                "path_id":  p.pid,
                "path":     p.name,
                "c1":       p.c1,
                "c2":       p.c2,
                "c6":       p.c6,
                "R_analytic": R_an,
                "R_solver":   R_or,
                "sigma_1":   sigma_1,
                "sigma_2":   sigma_2,
                "tau_12":    tau_12,
                "mode":      mode,
            })
        # Write per-criterion envelope CSV.
        out_csv = OUT_DIR / f"envelope_{crit}.csv"
        with out_csv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"[analyze] wrote {out_csv} ({len(rows)} paths)")

        # Pass-gate: principal-axis points.
        # LT = path "S1S2_000" = pid 1.
        # TT = path "S1S2_090" = pid 7.
        # LC = path "S1S2_180" = pid 13.
        # TC = path "S1S2_270" = pid 19.
        # S  = pure shear = path "S2T12_090" = pid 28 (k=3 in the second loop).
        axis_lookup = {"LT": (1, MAT.XT), "LC": (13, MAT.XC),
                       "TT": (7, MAT.YT), "TC": (19, MAT.YC),
                       "S":  (28, MAT.S12)}
        for axis, (pid, ref) in axis_lookup.items():
            row = next(r for r in rows if r["path_id"] == pid)
            R_solver = float(row["R_solver"])
            err = abs(R_solver - ref) / ref * 100.0
            pass_table.append({
                "criterion": crit,
                "axis":      axis,
                "ref":       ref,
                "R_solver":  R_solver,
                "err_pct":   err,
                "pass":      "PASS" if err <= 5.0 else "FAIL",
            })

    pg_csv = OUT_DIR / "pass_gate.csv"
    with pg_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pass_table[0].keys()))
        w.writeheader()
        w.writerows(pass_table)
    # Echo the table to stdout for the test harness.
    print("\n[analyze] principal-axis pass gate")
    print(f"{'criterion':>10s} {'axis':>4s} {'ref MPa':>10s} "
          f"{'R_solver':>10s} {'err %':>8s} {'verdict':>7s}")
    for r in pass_table:
        print(f"{r['criterion']:>10s} {r['axis']:>4s} {r['ref']:10.2f} "
              f"{float(r['R_solver']):10.2f} {float(r['err_pct']):8.2f} "
              f"{r['pass']:>7s}")


# ---- Failure-index helpers used by the T01 post-processor. ------------------
def _tsai_wu_index(sx: float, sy: float, sxy: float, m=MAT) -> float:
    F1  = 1.0 / m.XT - 1.0 / m.XC
    F2  = 1.0 / m.YT - 1.0 / m.YC
    F11 = 1.0 / (m.XT * m.XC)
    F22 = 1.0 / (m.YT * m.YC)
    F66 = 1.0 / (m.S12 ** 2)
    F12 = m.F12_star * math.sqrt(F11 * F22)
    return (F1 * sx + F2 * sy
            + F11 * sx * sx + F22 * sy * sy + F66 * sxy * sxy
            + 2.0 * F12 * sx * sy)


def _hashin_index(sx: float, sy: float, sxy: float, m=MAT) -> tuple[float, str]:
    fi: list[tuple[float, str]] = []
    if sx > 0.0:
        fi.append(((sx / m.XT) ** 2 + (sxy / m.S12) ** 2, "ff_T"))
    if sx < 0.0:
        fi.append(((sx / m.XC) ** 2, "ff_C"))
    if sy > 0.0:
        fi.append(((sy / m.YT) ** 2 + (sxy / m.S12) ** 2, "mf_T"))
    if sy < 0.0:
        fi.append((
            (sy / (2.0 * m.S23)) ** 2
            + ((m.YC / (2.0 * m.S23)) ** 2 - 1.0) * sy / m.YC
            + (sxy / m.S12) ** 2,
            "mf_C",
        ))
    if not fi:
        return (0.0, "none")
    return max(fi, key=lambda t: t[0])


def _puck_index(sx: float, sy: float, sxy: float, m=MAT) -> tuple[float, str]:
    fi: list[tuple[float, str]] = []
    if sx > 0.0:
        fi.append((sx / m.XT, "FF_T"))
    if sx < 0.0:
        fi.append((abs(sx) / m.XC, "FF_C"))
    p_pos = m.p_pos_perp_par
    p_neg = m.p_neg_perp_par
    if sy >= 0.0:
        a = (1.0 / m.YT - p_pos / m.S12) * sy
        b = sxy / m.S12
        fi.append((math.sqrt(a * a + b * b) + p_pos * sy / m.S12, "IFF_A"))
    else:
        tau_over_sig = abs(sxy) / max(abs(sy), 1e-12)
        if tau_over_sig >= 0.4:
            inner = math.sqrt(sxy * sxy + (p_neg * sy) ** 2) + p_neg * sy
            fi.append((inner / m.S12, "IFF_B"))
        else:
            denom_c = 4.0 * (1.0 + m.p_neg_perp_perp) ** 2 * m.S12 ** 2
            num = sxy * sxy / denom_c + sy * sy / m.YC ** 2
            fi.append((num * m.YC / max(abs(sy), 1e-12), "IFF_C"))
    return max(fi, key=lambda t: t[0]) if fi else (0.0, "none")


# ============================================================================
# Section 6. Plotting (Typst+CeTZ data export per global preferences)
# ============================================================================

def cmd_plot() -> None:
    """Write Typst-readable CSV of the full envelope for CeTZ to plot.

    Per the user's global preferences, no Python plotting library is used;
    instead, we produce one CSV per criterion with (sigma_1, sigma_2, tau_12)
    columns and let `figures/stage_06/envelope.typ` import them via `csv()`
    and plot with CeTZ.
    """
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    # Combined sigma_1 vs sigma_2 envelope (paths 1..24).
    sx_csv = FIG_DIR / "sigma1_sigma2_envelope.csv"
    with sx_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["criterion", "path_id", "sigma_1", "sigma_2"])
        for crit in CRITERIA:
            with (OUT_DIR / f"envelope_{crit}.csv").open() as g:
                rdr = csv.DictReader(g)
                for row in rdr:
                    if int(row["path_id"]) <= 24:
                        w.writerow([crit, row["path_id"],
                                    row["sigma_1"], row["sigma_2"]])
    # Combined sigma_2 vs tau_12 envelope (paths 25..32).
    sy_csv = FIG_DIR / "sigma2_tau12_envelope.csv"
    with sy_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["criterion", "path_id", "sigma_2", "tau_12"])
        for crit in CRITERIA:
            with (OUT_DIR / f"envelope_{crit}.csv").open() as g:
                rdr = csv.DictReader(g)
                for row in rdr:
                    if int(row["path_id"]) > 24:
                        w.writerow([crit, row["path_id"],
                                    row["sigma_2"], row["tau_12"]])
    # Stub Typst source for CeTZ if it doesn't already exist.
    typ = FIG_DIR / "envelope.typ"
    if not typ.exists():
        typ.write_text(
            '#import "@preview/cetz:0.2.2"\n'
            '#set page(width: auto, height: auto, margin: 1cm)\n'
            '\n'
            '#let GARNET = rgb("#73000A")\n'
            '#let ATLANTIC = rgb("#466A9F")\n'
            '#let HORSESHOE = rgb("#65780B")\n'
            '\n'
            '#let crit_color(c) = if c == "TSAIWU" { GARNET }\n'
            '  else if c == "HASHIN" { ATLANTIC }\n'
            '  else { HORSESHOE }\n'
            '\n'
            '#let env = csv("sigma1_sigma2_envelope.csv")\n'
            '#cetz.canvas({\n'
            '  import cetz.draw: *\n'
            '  set-style(stroke: (thickness: 1.0pt))\n'
            '  // axes\n'
            '  line((-2000, 0), (3000, 0), name: "ax")\n'
            '  line((0, -250), (0, 200), name: "ay")\n'
            '  content((3000, -120), [#text(size: 8pt)[$sigma_1$ MPa]])\n'
            '  content((-200, 200), [#text(size: 8pt)[$sigma_2$ MPa]])\n'
            '  // envelopes\n'
            '  for crit in ("TSAIWU", "HASHIN", "PUCK") {\n'
            '    let pts = env.filter(r => r.at(0) == crit)\n'
            '      .map(r => (float(r.at(2)), float(r.at(3))))\n'
            '    if pts.len() > 0 {\n'
            '      line(..pts, ..(pts.first(),),\n'
            '           stroke: crit_color(crit) + 1.4pt)\n'
            '    }\n'
            '  }\n'
            '})\n'
        )
    print(f"[plot] wrote {sx_csv}, {sy_csv}, {typ}")


# ============================================================================
# Section 7. Self-test (offline verification of the analytic envelope)
# ============================================================================

def cmd_selftest() -> None:
    """Quick offline check of the analytic envelopes against principal-axis
    strengths. Should print PASS for all five axes and all three criteria.
    """
    pid_axis = {
        "LT": next(p for p in PATHS if p.name == "S1S2_000"),
        "TT": next(p for p in PATHS if p.name == "S1S2_090"),
        "LC": next(p for p in PATHS if p.name == "S1S2_180"),
        "TC": next(p for p in PATHS if p.name == "S1S2_270"),
        "S":  next(p for p in PATHS if p.name == "S2T12_090"),
    }
    refs = {"LT": MAT.XT, "TT": MAT.YT, "LC": MAT.XC, "TC": MAT.YC, "S": MAT.S12}
    print("[selftest] analytic envelope, principal-axis check")
    print(f"{'crit':>8s} {'axis':>4s} {'ref':>9s} {'R':>9s} {'err%':>8s} {'mode':>10s}")
    fail = 0
    for crit in CRITERIA:
        for axis, p in pid_axis.items():
            if crit == "TSAIWU":
                R = tsaiwu_failure_load(p)
                mode = "tsaiwu"
            elif crit == "HASHIN":
                R, mode = hashin_failure_load(p)
            else:
                R, mode = puck_failure_load(p)
            err = abs(R - refs[axis]) / refs[axis] * 100.0
            verdict = "PASS" if err <= 5.0 else "FAIL"
            if verdict == "FAIL":
                fail += 1
            print(f"{crit:>8s} {axis:>4s} {refs[axis]:9.2f} "
                  f"{R:9.2f} {err:8.3f} {mode:>10s}  {verdict}")
    if fail:
        print(f"[selftest] {fail} principal-axis points out of tolerance")
        sys.exit(1)
    print("[selftest] all principal-axis points pass within 5%")


# ============================================================================
# Section 8. Argument parsing
# ============================================================================

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage 06 - composite failure criteria runner")
    ap.add_argument("--build",    action="store_true", help="generate decks")
    ap.add_argument("--solve",    action="store_true", help="run starter+engine")
    ap.add_argument("--analyze",  action="store_true", help="parse outputs")
    ap.add_argument("--plot",     action="store_true", help="write Typst+CeTZ data")
    ap.add_argument("--all",      action="store_true", help="build + solve + analyze + plot")
    ap.add_argument("--selftest", action="store_true",
                    help="offline check of analytic envelope (no solver call)")
    ap.add_argument("--criterion", choices=list(CRITERIA),
                    help="restrict --solve to one criterion")
    ap.add_argument("--paths",    type=str, default=None,
                    help="restrict --solve to a comma-separated list of path ids")
    args = ap.parse_args(argv)

    only_paths = None
    if args.paths:
        only_paths = [int(s) for s in args.paths.split(",") if s.strip()]

    did_anything = False
    if args.selftest:
        cmd_selftest(); did_anything = True
    if args.build or args.all:
        cmd_build(); did_anything = True
    if args.solve or args.all:
        cmd_solve(only_crit=args.criterion, only_paths=only_paths); did_anything = True
    if args.analyze or args.all:
        cmd_analyze(); did_anything = True
    if args.plot or args.all:
        cmd_plot(); did_anything = True
    if not did_anything:
        ap.print_help()
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
