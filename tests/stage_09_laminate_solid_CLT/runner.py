#!/usr/bin/env python3
"""
Stage 09 runner — solid laminate, classical lamination theory verification.

Pipeline:
  1. Compute analytic CLT A-matrix (Jones 1999 Ch. 4) and write to CSV.
  2. Build GMSH structured-hex mesh, one named volume per ply layer.
  3. Convert to Abaqus .inp (meshio), then to OpenRadioss .rad via inp2rad.
  4. Template per-ply /PROP/TYPE14 + /MAT/LAW25 cards with rotated frames.
  5. Template starter + engine decks for each (layup, refinement, case).
  6. Invoke OpenRadioss starter + engine inside Lima/Apptainer.
  7. Convert .anim -> .vtkhdf (Kitware openradioss-to-vtkhdf).
  8. Read VTKHDF, extract per-element stress; integrate through thickness.
  9. Compute FEM A-matrix and compare to CLT analytic, per-component 2 percent.
 10. Write results CSV and stage summary.

Subcommands (CLI):
  python runner.py clt              # Step 1 only — analytic CLT to CSV.
  python runner.py mesh             # Step 2 — GMSH mesh build.
  python runner.py decks            # Step 4-5 — template OpenRadioss decks.
  python runner.py run              # Step 6 — invoke starter+engine.
  python runner.py post             # Step 7-9 — recover A^FEM and compare.
  python runner.py all              # End-to-end.
  python runner.py check            # Sanity checks only (no FEM).

Solid HEXA8 elements only. /MAT/LAW25 + /PROP/TYPE14. /IMPL/LINEAR.
SI units everywhere (m, kg, s, Pa, N).
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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

STAGE_DIR = Path(__file__).resolve().parent
REPO_ROOT = STAGE_DIR.parent.parent  # FEA_AP-PLY/
RUNS_DIR = STAGE_DIR / "runs"
RESULTS_DIR = STAGE_DIR / "results"
FIGURES_DIR = STAGE_DIR / "figures"
MESH_DIR = STAGE_DIR / "mesh"
DECKS_DIR = STAGE_DIR / "decks"

for _d in (RUNS_DIR, RESULTS_DIR, FIGURES_DIR, MESH_DIR, DECKS_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Material card — IM7/8552, SI units
# Sources:
#   Soden, Hinton, Kaddour 1998 (closely related IM7/8551-7 system).
#   Camanho, Maimi, Davila 2007 (IM7/8552 baseline).
#   CMH-17 Volume 3 (2012).
# Same card as stage 07 / stage 13. Strengths kept consistent so /MAT/LAW25
# is well-formed; in this stage they are not exercised (eps0 = 1e-4).
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class IM78552:
    rho: float = 1580.0          # kg/m^3
    E1: float = 171.42e9         # Pa
    E2: float = 9.08e9           # Pa
    E3: float = 9.08e9           # Pa
    nu12: float = 0.32
    nu13: float = 0.32
    nu23: float = 0.50
    G12: float = 5.29e9          # Pa
    G13: float = 5.29e9          # Pa
    XT: float = 2323.5e6         # Pa
    XC: float = 1200.1e6
    YT: float = 62.3e6
    YC: float = 199.8e6
    SL: float = 92.3e6

    @property
    def G23(self) -> float:
        return self.E2 / (2.0 * (1.0 + self.nu23))


# ---------------------------------------------------------------------------
# Layup definitions
# ---------------------------------------------------------------------------

PLY_THICKNESS_M = 0.18e-3  # 0.18 mm in metres (Kok ap_ply_model_creation default)

@dataclass(frozen=True)
class Layup:
    name: str
    thetas_deg: Tuple[float, ...]   # bottom (z = 0) -> top (z = t)

    @property
    def n_plies(self) -> int:
        return len(self.thetas_deg)

    @property
    def total_thickness_m(self) -> float:
        return self.n_plies * PLY_THICKNESS_M


LAYUP_A = Layup(name="A_crossply", thetas_deg=(0.0, 90.0, 90.0, 0.0))
LAYUP_B = Layup(name="B_quasiiso",
                thetas_deg=(0.0, 45.0, -45.0, 90.0,
                            90.0, -45.0, 45.0, 0.0))


# ---------------------------------------------------------------------------
# Geometry / mesh / loading parameters
# ---------------------------------------------------------------------------

L_X_M = 100.0e-3   # plate side (x), 100 mm
L_Y_M = 100.0e-3   # plate side (y), 100 mm
EPS0 = 1.0e-4       # applied uniform strain (small, linear-elastic)
GAMMA0 = 1.0e-4     # applied shear strain
N_X = 20            # in-plane elements along x (baseline)
N_Y = 20            # in-plane elements along y
INTERIOR_CROP_M = 5.0e-3  # crop 5 mm from each in-plane edge for A-matrix recovery

REFINEMENTS = {
    "ref1": 1,   # one element per ply (baseline)
    "ref2": 2,   # two elements per ply (first refinement)
    "ref4": 4,   # optional fallback
}

CASES = ("case1_xx", "case2_yy", "case3_xy")


# ---------------------------------------------------------------------------
# Section 7 — analytic CLT A-matrix (Jones 1999 Ch. 4)
# ---------------------------------------------------------------------------

def reduced_stiffness_Q(mat: IM78552) -> np.ndarray:
    """Jones (1999) eq. (2.66): plane-stress reduced stiffness in ply axes.

    Returns 3x3 Q in (1, 2, 6) ordering, units Pa.
    """
    nu21 = mat.nu12 * mat.E2 / mat.E1
    denom = 1.0 - mat.nu12 * nu21
    Q = np.zeros((3, 3))
    Q[0, 0] = mat.E1 / denom
    Q[1, 1] = mat.E2 / denom
    Q[0, 1] = mat.nu12 * mat.E2 / denom
    Q[1, 0] = Q[0, 1]
    Q[2, 2] = mat.G12
    return Q


def rotated_Q_bar(Q: np.ndarray, theta_deg: float) -> np.ndarray:
    """Jones (1999) eq. (2.84-2.85): rotated reduced stiffness Q_bar(theta).

    theta is measured counter-clockwise from the laminate x-axis to the
    ply 1-axis. Returns 3x3 Q_bar in (xx, yy, xy) ordering.
    """
    th = math.radians(theta_deg)
    c = math.cos(th)
    s = math.sin(th)
    c2 = c * c
    s2 = s * s
    c4 = c2 * c2
    s4 = s2 * s2
    cs = c * s

    Q11, Q22, Q12, Q66 = Q[0, 0], Q[1, 1], Q[0, 1], Q[2, 2]

    Qb = np.zeros((3, 3))
    Qb[0, 0] = Q11 * c4 + 2.0 * (Q12 + 2.0 * Q66) * s2 * c2 + Q22 * s4
    Qb[1, 1] = Q11 * s4 + 2.0 * (Q12 + 2.0 * Q66) * s2 * c2 + Q22 * c4
    Qb[0, 1] = (Q11 + Q22 - 4.0 * Q66) * s2 * c2 + Q12 * (s4 + c4)
    Qb[1, 0] = Qb[0, 1]
    Qb[2, 2] = (Q11 + Q22 - 2.0 * Q12 - 2.0 * Q66) * s2 * c2 + Q66 * (s4 + c4)
    Qb[0, 2] = (Q11 - Q12 - 2.0 * Q66) * cs * c2 + (Q12 - Q22 + 2.0 * Q66) * s2 * cs
    Qb[1, 2] = (Q11 - Q12 - 2.0 * Q66) * s2 * cs + (Q12 - Q22 + 2.0 * Q66) * cs * c2
    Qb[2, 0] = Qb[0, 2]
    Qb[2, 1] = Qb[1, 2]
    return Qb


def laminate_A_matrix(mat: IM78552, layup: Layup) -> np.ndarray:
    """Jones (1999) eq. (4.24a): A_ij = sum_k Q_bar_ij(theta_k) * t_k.

    Returns 3x3 A (units N/m) for the given layup with all plies of equal
    thickness PLY_THICKNESS_M.
    """
    Q = reduced_stiffness_Q(mat)
    A = np.zeros((3, 3))
    for theta in layup.thetas_deg:
        A += rotated_Q_bar(Q, theta) * PLY_THICKNESS_M
    return A


def apparent_moduli_from_A(A: np.ndarray, t_total: float) -> Dict[str, float]:
    """Daniel & Ishai 2006 Ch. 5: apparent in-plane moduli from A.

    Valid for symmetric laminates (B = 0) with negligible A_16, A_26.
    """
    A11, A22, A12, A66 = A[0, 0], A[1, 1], A[0, 1], A[2, 2]
    Ex = (A11 * A22 - A12 * A12) / (t_total * A22)
    Ey = (A11 * A22 - A12 * A12) / (t_total * A11)
    nu_xy = A12 / A22
    G_xy = A66 / t_total
    return dict(Ex=Ex, Ey=Ey, nu_xy=nu_xy, G_xy=G_xy)


def write_A_csv(A: np.ndarray, path: Path, header: str = "") -> None:
    """Write a 3x3 A-matrix to CSV in row-major (xx, yy, xy) ordering."""
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        if header:
            w.writerow([f"# {header}"])
        w.writerow(["", "xx", "yy", "xy"])
        for i, row_label in enumerate(("xx", "yy", "xy")):
            w.writerow([row_label] + [f"{A[i, j]:.10e}" for j in range(3)])


def cmd_clt(args: argparse.Namespace) -> int:
    """Step 1 — write analytic CLT A-matrix CSVs for both layups."""
    mat = IM78552()
    Q = reduced_stiffness_Q(mat)

    Q_path = RESULTS_DIR / "Q_ply_axes.csv"
    write_A_csv(Q, Q_path, header="Ply-axes reduced stiffness Q [Pa]")
    print(f"[clt] wrote {Q_path}")

    for layup in (LAYUP_A, LAYUP_B):
        A = laminate_A_matrix(mat, layup)
        path = RESULTS_DIR / f"cltA_{layup.name}_reference.csv"
        write_A_csv(A, path, header=f"CLT A-matrix [N/m] for layup {layup.name}, "
                                    f"plies {layup.thetas_deg} deg, "
                                    f"t_ply = {PLY_THICKNESS_M} m")
        print(f"[clt] {layup.name}  A_11 = {A[0,0]:.4e} N/m  A_22 = {A[1,1]:.4e}  "
              f"A_12 = {A[0,1]:.4e}  A_66 = {A[2,2]:.4e}  "
              f"A_16 = {A[0,2]:.4e}  A_26 = {A[1,2]:.4e}")

        mods = apparent_moduli_from_A(A, layup.total_thickness_m)
        mod_path = RESULTS_DIR / f"apparent_moduli_{layup.name}_CLT.csv"
        with mod_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["quantity", "value", "unit"])
            w.writerow(["Ex", f"{mods['Ex']:.6e}", "Pa"])
            w.writerow(["Ey", f"{mods['Ey']:.6e}", "Pa"])
            w.writerow(["nu_xy", f"{mods['nu_xy']:.6e}", "-"])
            w.writerow(["G_xy", f"{mods['G_xy']:.6e}", "Pa"])
        print(f"[clt] wrote {mod_path}  Ex = {mods['Ex']/1e9:.3f} GPa")

    # Sanity: cross-ply and quasi-iso symmetry to 1e-10 relative.
    A_A = laminate_A_matrix(mat, LAYUP_A)
    A_B = laminate_A_matrix(mat, LAYUP_B)
    rel = lambda x, y: abs(x - y) / max(abs(x), abs(y), 1e-30)
    assert rel(A_A[0, 0], A_A[1, 1]) < 1e-10, "cross-ply A11 != A22"
    assert abs(A_A[0, 2]) < 1e-8 * abs(A_A[0, 0]), "cross-ply A16 != 0"
    assert rel(A_B[0, 0], A_B[1, 1]) < 1e-10, "quasi-iso A11 != A22"
    assert abs(A_B[0, 2]) < 1e-8 * abs(A_B[0, 0]), "quasi-iso A16 != 0"
    A66_pred = (A_B[0, 0] - A_B[0, 1]) / 2.0
    assert rel(A_B[2, 2], A66_pred) < 1e-10, \
        f"quasi-iso A66 != (A11-A12)/2: {A_B[2,2]} vs {A66_pred}"
    print("[clt] symmetry sanity checks passed")
    return 0


# ---------------------------------------------------------------------------
# Step 2 — GMSH structured hex mesh, one named volume per ply
# Implementation note: this stage uses the GMSH Python API. If gmsh is not
# installed, the runner stops cleanly with a diagnostic. The mesh is also
# written as Abaqus .inp via meshio for downstream inp2rad conversion.
# ---------------------------------------------------------------------------

def cmd_mesh(args: argparse.Namespace) -> int:
    layup = LAYUP_A if args.layup == "A" else LAYUP_B
    n_per_ply = REFINEMENTS[args.refinement]
    out_msh = MESH_DIR / f"{layup.name}_{args.refinement}.msh"
    out_inp = MESH_DIR / f"{layup.name}_{args.refinement}.inp"

    try:
        import gmsh  # type: ignore
    except ImportError:
        print("[mesh] gmsh Python API not available; install with: pip install gmsh",
              file=sys.stderr)
        return 1

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"stage09_{layup.name}_{args.refinement}")

    # Build N_PLIES * n_per_ply layers stacked in z.
    # Use a structured 2D rectangle on the x-y plane, then extrude in z
    # in one extrude call per layer so each layer's volume gets a
    # distinct physical group.
    rect_tag = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, L_X_M, L_Y_M)
    gmsh.model.occ.synchronize()

    # Transfinite the base rectangle.
    surface_tag = rect_tag
    base_curves = [c[1] for c in gmsh.model.getBoundary([(2, surface_tag)],
                                                        oriented=False, recursive=False)]
    for c in base_curves:
        gmsh.model.mesh.setTransfiniteCurve(c, N_X + 1)
    gmsh.model.mesh.setTransfiniteSurface(surface_tag)
    gmsh.model.mesh.setRecombine(2, surface_tag)

    # Per-ply extrude. For each ply k and each sub-element s in [0, n_per_ply):
    #   layer thickness = PLY_THICKNESS_M / n_per_ply.
    cur_surface = surface_tag
    ply_volumes = []  # list of (ply_index, theta_deg, [vol_tags])
    for ply_idx, theta in enumerate(layup.thetas_deg):
        sub_vols = []
        for sub in range(n_per_ply):
            dz = PLY_THICKNESS_M / n_per_ply
            ext = gmsh.model.occ.extrude(
                [(2, cur_surface)], 0, 0, dz,
                numElements=[1], heights=[1.0], recombine=True
            )
            gmsh.model.occ.synchronize()
            # ext is a list of dim-tags; the first is the new top surface,
            # the second is the swept volume, the rest are the side surfaces.
            new_top = next(t for (d, t) in ext if d == 2 and t != cur_surface)
            new_vol = next(t for (d, t) in ext if d == 3)
            sub_vols.append(new_vol)
            cur_surface = new_top
        ply_volumes.append((ply_idx, theta, sub_vols))

    gmsh.model.occ.synchronize()

    # Physical groups so the .inp has one set per ply layer.
    for ply_idx, theta, vols in ply_volumes:
        ptag = 100 + ply_idx
        gmsh.model.addPhysicalGroup(3, vols, ptag)
        gmsh.model.setPhysicalName(3, ptag,
                                   f"ply_{ply_idx:02d}_theta_{int(round(theta)):+04d}")

    # Tag boundary surfaces (loaded edges + traction-free top/bottom).
    eps = 1e-9
    bbox_x0, bbox_x1 = 0.0, L_X_M
    bbox_y0, bbox_y1 = 0.0, L_Y_M
    bbox_z0, bbox_z1 = 0.0, layup.total_thickness_m

    def faces_in(box):
        x0, x1, y0, y1, z0, z1 = box
        return [t for (d, t) in gmsh.model.getEntitiesInBoundingBox(
            x0 - eps, y0 - eps, z0 - eps, x1 + eps, y1 + eps, z1 + eps, 2)]

    g_x0 = gmsh.model.addPhysicalGroup(2, faces_in(
        (bbox_x0, bbox_x0, bbox_y0, bbox_y1, bbox_z0, bbox_z1)), 200)
    gmsh.model.setPhysicalName(2, g_x0, "edge_x0")
    g_x1 = gmsh.model.addPhysicalGroup(2, faces_in(
        (bbox_x1, bbox_x1, bbox_y0, bbox_y1, bbox_z0, bbox_z1)), 201)
    gmsh.model.setPhysicalName(2, g_x1, "edge_x1")
    g_y0 = gmsh.model.addPhysicalGroup(2, faces_in(
        (bbox_x0, bbox_x1, bbox_y0, bbox_y0, bbox_z0, bbox_z1)), 202)
    gmsh.model.setPhysicalName(2, g_y0, "edge_y0")
    g_y1 = gmsh.model.addPhysicalGroup(2, faces_in(
        (bbox_x0, bbox_x1, bbox_y1, bbox_y1, bbox_z0, bbox_z1)), 203)
    gmsh.model.setPhysicalName(2, g_y1, "edge_y1")
    g_z0 = gmsh.model.addPhysicalGroup(2, faces_in(
        (bbox_x0, bbox_x1, bbox_y0, bbox_y1, bbox_z0, bbox_z0)), 204)
    gmsh.model.setPhysicalName(2, g_z0, "face_z0")
    g_z1 = gmsh.model.addPhysicalGroup(2, faces_in(
        (bbox_x0, bbox_x1, bbox_y0, bbox_y1, bbox_z1, bbox_z1)), 205)
    gmsh.model.setPhysicalName(2, g_z1, "face_z1")

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(1)  # linear HEXA8
    gmsh.write(str(out_msh))
    print(f"[mesh] wrote {out_msh}")

    gmsh.finalize()

    # Convert .msh -> .inp via meshio so downstream inp2rad can read it.
    try:
        import meshio  # type: ignore
        m = meshio.read(out_msh)
        meshio.write(out_inp, m)
        print(f"[mesh] wrote {out_inp}")
    except ImportError:
        print("[mesh] WARNING: meshio not installed; skipping .inp export",
              file=sys.stderr)

    # Persist a per-layer manifest so the deck templater knows which element
    # set goes with which theta.
    manifest = dict(
        layup=layup.name,
        thetas_deg=list(layup.thetas_deg),
        n_per_ply=n_per_ply,
        ply_thickness_m=PLY_THICKNESS_M,
        total_thickness_m=layup.total_thickness_m,
        L_x_m=L_X_M, L_y_m=L_Y_M,
        n_x=N_X, n_y=N_Y,
    )
    (MESH_DIR / f"{layup.name}_{args.refinement}_manifest.json").write_text(
        json.dumps(manifest, indent=2)
    )
    return 0


# ---------------------------------------------------------------------------
# Step 4-5 — OpenRadioss deck templating
# ---------------------------------------------------------------------------

STARTER_HEADER = """\
#RADIOSS STARTER
/BEGIN
{run_name}
      2024         0
                 g                 Pa
                 m                 s
/UNIT/1
SI units
kg                  m                  s
/IMPL/LINEAR
"""

ENGINE_DECK = """\
#RADIOSS ENGINE
/RUN/{run_name}/1
1.0
/IMPL/LINEAR
/PRINT/-1
/TFILE
0.05
/ANIM/DT
0.0  1.0
/ANIM/BRICK/STRESS
/ANIM/BRICK/STRAIN
/ANIM/BRICK/TENS/STRESS/ALL
/ANIM/BRICK/TENS/STRAIN/ALL
/H3D/DT
0.0  1.0
/STOP
"""


def render_law25_block(mat_id: int, name: str, mat: IM78552) -> str:
    """/MAT/LAW25 CRASURV card. SI units (kg, m, s, Pa).

    The LAW25 line ordering follows Altair help
    `mat_law25_composite_starter_r.htm`. Strengths are written so the deck
    is well-formed; this stage stays linear-elastic so the failure surface
    is not exercised.
    """
    return (
        f"/MAT/LAW25/{mat_id}\n"
        f"{name}\n"
        f"   {mat.rho:.6e}\n"
        f"   {mat.E1:.6e}   {mat.E2:.6e}   {mat.nu12:.6e}\n"
        f"   {mat.G12:.6e}   {mat.G23:.6e}   {mat.G13:.6e}\n"
        f"   {mat.E3:.6e}   {mat.nu23:.6e}\n"
        f"   {mat.XT:.6e}   {mat.XC:.6e}\n"
        f"   {mat.YT:.6e}   {mat.YC:.6e}   {mat.SL:.6e}\n"
    )


def render_prop_type14_block(prop_id: int, mat_id: int, theta_deg: float,
                             ply_idx: int) -> str:
    """/PROP/TYPE14 with orthotropic frame from theta.

    First in-plane vector V = (cos theta, sin theta, 0). OpenRadioss builds
    the second in-plane axis automatically (Iframe = 0). Iorth = 1 so the
    frame co-rotates with the element (per the audit, row 7).
    """
    th = math.radians(theta_deg)
    Vx, Vy, Vz = math.cos(th), math.sin(th), 0.0
    return (
        f"/PROP/TYPE14/{prop_id}\n"
        f"ply_{ply_idx:02d}_theta_{int(round(theta_deg)):+04d}\n"
        f"# Ihbe Ismstr Itetra4 Iframe Iorth Ip Icpre\n"
        f"   2     0      0       0      1     0   0\n"
        f"# Vx, Vy, Vz (orthotropic frame, theta = {theta_deg:.2f} deg)\n"
        f"   {Vx:.10e}   {Vy:.10e}   {Vz:.10e}\n"
        f"# pointer to /MAT/LAW25\n"
        f"   {mat_id}\n"
    )


def render_bc_block(case: str, layup: Layup) -> str:
    """Boundary conditions for case 1, 2, or 3.

    Case 1 (xx): u_x = 0 at edge_x0, u_x = eps0 * Lx at edge_x1, free in y.
    Case 2 (yy): u_y = 0 at edge_y0, u_y = eps0 * Ly at edge_y1, free in x.
    Case 3 (xy): u_x = 0 and u_y = 0 at edge_y0; u_x = gamma0 * Ly and
                 u_y = 0 at edge_y1; u_y = 0 at edge_x0 and edge_x1.

    All cases also pin one corner in u_z and one corner-pair in the
    transverse in-plane DOF to remove rigid-body modes and the in-plane
    rotation. The deck templater materializes the BC blocks against
    GMSH physical groups 200..205 (edge_x0..face_z1).
    """
    if case == "case1_xx":
        u_xL = EPS0 * L_X_M
        body = (
            "# /BCS for case1_xx — uniaxial extension in x\n"
            "/BCS/100\nfix_x_at_x0\n"
            "  1 1 1 1 0 0 0 0\n"
            "  edge_x0\n"
            f"/IMPDISP/100\nimp_x_at_x1   amplitude {u_xL:.10e}\n"
            "  1 0 0 0 0 0\n"
            "  edge_x1\n"
            "/BCS/101\nfix_z_face_z0\n"
            "  0 0 1 0 0 0 0 0\n"
            "  face_z0\n"
            "# remove in-plane rotation: pin u_y at one corner of edge_x0\n"
            "/BCS/102\npin_corner_y\n"
            "  0 1 0 0 0 0 0 0\n"
            "  corner_x0_y0\n"
        )
    elif case == "case2_yy":
        u_yL = EPS0 * L_Y_M
        body = (
            "# /BCS for case2_yy — uniaxial extension in y\n"
            "/BCS/100\nfix_y_at_y0\n"
            "  0 1 0 0 0 0 0 0\n"
            "  edge_y0\n"
            f"/IMPDISP/100\nimp_y_at_y1   amplitude {u_yL:.10e}\n"
            "  0 1 0 0 0 0\n"
            "  edge_y1\n"
            "/BCS/101\nfix_z_face_z0\n"
            "  0 0 1 0 0 0 0 0\n"
            "  face_z0\n"
            "/BCS/102\npin_corner_x\n"
            "  1 0 0 0 0 0 0 0\n"
            "  corner_x0_y0\n"
        )
    elif case == "case3_xy":
        u_top = GAMMA0 * L_Y_M
        body = (
            "# /BCS for case3_xy — pure in-plane shear\n"
            "/BCS/100\nfix_xy_at_y0\n"
            "  1 1 0 0 0 0 0 0\n"
            "  edge_y0\n"
            f"/IMPDISP/100\nimp_x_at_y1   amplitude {u_top:.10e}\n"
            "  1 0 0 0 0 0\n"
            "  edge_y1\n"
            "/BCS/101\nfix_y_at_y1\n"
            "  0 1 0 0 0 0 0 0\n"
            "  edge_y1\n"
            "/BCS/102\nfix_y_edges_x\n"
            "  0 1 0 0 0 0 0 0\n"
            "  edge_x0\n"
            "/BCS/103\nfix_y_edges_x\n"
            "  0 1 0 0 0 0 0 0\n"
            "  edge_x1\n"
            "/BCS/104\nfix_z_face_z0\n"
            "  0 0 1 0 0 0 0 0\n"
            "  face_z0\n"
        )
    else:
        raise ValueError(f"unknown case {case}")
    body += (
        "/FUNCT/1\nramp_0_to_1\n"
        "  0.0   0.0\n"
        "  1.0   1.0\n"
    )
    return body


def cmd_decks(args: argparse.Namespace) -> int:
    """Step 4-5 — write deck templates for every (layup, refinement, case)."""
    layups = (LAYUP_A, LAYUP_B)
    refs = (args.refinement,) if args.refinement else ("ref1", "ref2")
    mat = IM78552()

    n_decks = 0
    for layup in layups:
        for ref in refs:
            for case in CASES:
                run_name = f"stage_09_{layup.name}_{ref}_{case}"
                run_dir = RUNS_DIR / run_name
                run_dir.mkdir(parents=True, exist_ok=True)

                # Per-ply /MAT and /PROP blocks (one of each per ply layer).
                mat_blocks = []
                prop_blocks = []
                for ply_idx, theta in enumerate(layup.thetas_deg):
                    mat_id = ply_idx + 1
                    prop_id = ply_idx + 1
                    mat_blocks.append(render_law25_block(
                        mat_id, f"IM7_8552_ply_{ply_idx:02d}", mat))
                    prop_blocks.append(render_prop_type14_block(
                        prop_id, mat_id, theta, ply_idx))

                bc_block = render_bc_block(case, layup)

                starter = (
                    STARTER_HEADER.format(run_name=run_name)
                    + "# --- MESH BLOCK INSERTED FROM inp2rad OUTPUT ---\n"
                    + f"#INCLUDE {layup.name}_{ref}_mesh.rad\n"
                    + "# --- per-ply materials ---\n"
                    + "".join(mat_blocks)
                    + "# --- per-ply properties (with rotated orthotropic frame) ---\n"
                    + "".join(prop_blocks)
                    + "# --- boundary conditions ---\n"
                    + bc_block
                    + "/END\n"
                )

                (run_dir / f"{run_name}_0000.rad").write_text(starter)
                (run_dir / f"{run_name}_0001.rad").write_text(
                    ENGINE_DECK.format(run_name=run_name)
                )
                n_decks += 1

    print(f"[decks] templated {n_decks} starter+engine pairs in {RUNS_DIR}")
    return 0


# ---------------------------------------------------------------------------
# Step 6 — invoke OpenRadioss starter + engine inside Lima/Apptainer
# Same shape as stages 1-8: one wrapper script in REPO_ROOT/scripts/.
# ---------------------------------------------------------------------------

def run_openradioss(deck_dir: Path, run_name: str) -> int:
    """Invoke OpenRadioss starter then engine on a deck pair.

    Looks for an env var OR_RUNNER which should point to a script that
    forwards arguments into Lima. If unset, falls back to assuming
    starter/engine are on PATH (Linux native).
    """
    or_runner = os.environ.get("OR_RUNNER")
    starter_in = deck_dir / f"{run_name}_0000.rad"
    engine_in = deck_dir / f"{run_name}_0001.rad"

    if or_runner:
        cmd_starter = [or_runner, "starter", "-i", str(starter_in), "-nt", "1"]
        cmd_engine = [or_runner, "engine", "-i", str(engine_in), "-nt", "1"]
    else:
        cmd_starter = ["starter_linux64_gf", "-i", str(starter_in), "-nt", "1"]
        cmd_engine = ["engine_linux64_gf", "-i", str(engine_in), "-nt", "1"]

    print(f"[run] {' '.join(cmd_starter)}")
    rc = subprocess.call(cmd_starter, cwd=str(deck_dir))
    if rc != 0:
        print(f"[run] starter failed rc={rc}", file=sys.stderr)
        return rc
    print(f"[run] {' '.join(cmd_engine)}")
    rc = subprocess.call(cmd_engine, cwd=str(deck_dir))
    if rc != 0:
        print(f"[run] engine failed rc={rc}", file=sys.stderr)
    return rc


def cmd_run(args: argparse.Namespace) -> int:
    layups = (LAYUP_A, LAYUP_B)
    refs = (args.refinement,) if args.refinement else ("ref1", "ref2")
    rc_max = 0
    for layup in layups:
        for ref in refs:
            for case in CASES:
                run_name = f"stage_09_{layup.name}_{ref}_{case}"
                run_dir = RUNS_DIR / run_name
                rc = run_openradioss(run_dir, run_name)
                rc_max = max(rc_max, abs(rc))
    return rc_max


# ---------------------------------------------------------------------------
# Step 7-9 — post-process: extract A^FEM, compare to A^CLT
# ---------------------------------------------------------------------------

def _read_fem_stress_vtkhdf(vtkhdf_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (centroids, stress_global_3x3_per_elem, volumes) at final state.

    vtkhdf_path points to the converted Kitware VTKHDF file. The schema
    (cell-data 'Stress' as 6-component or 9-component, cell type
    VTK_HEXAHEDRON) is read with PyVista.
    """
    import pyvista as pv  # type: ignore

    grid = pv.read(str(vtkhdf_path))
    if grid.n_cells == 0:
        raise RuntimeError(f"{vtkhdf_path}: no cells found")
    # Last time step: pyvista.read returns the last frame for multiblock VTKHDF.

    centroids = grid.cell_centers().points
    # Stress field: try 'Stress' then fall back to 'STRESS' (OpenRadioss naming).
    sname = next((n for n in grid.cell_data.keys()
                  if n.lower() in ("stress", "stresses", "sigma")), None)
    if sname is None:
        raise RuntimeError(f"{vtkhdf_path}: no stress field found in cell data: "
                           f"{list(grid.cell_data.keys())}")
    s_raw = np.asarray(grid.cell_data[sname])
    # OpenRadioss writes (sxx, syy, szz, sxy, syz, sxz) per cell.
    if s_raw.shape[1] == 6:
        sigma = np.zeros((s_raw.shape[0], 3, 3))
        sigma[:, 0, 0] = s_raw[:, 0]
        sigma[:, 1, 1] = s_raw[:, 1]
        sigma[:, 2, 2] = s_raw[:, 2]
        sigma[:, 0, 1] = sigma[:, 1, 0] = s_raw[:, 3]
        sigma[:, 1, 2] = sigma[:, 2, 1] = s_raw[:, 4]
        sigma[:, 0, 2] = sigma[:, 2, 0] = s_raw[:, 5]
    elif s_raw.shape[1] == 9:
        sigma = s_raw.reshape(-1, 3, 3)
    else:
        raise RuntimeError(f"{vtkhdf_path}: unexpected stress component count "
                           f"{s_raw.shape}")

    volumes = np.asarray(grid.compute_cell_sizes(volume=True)
                         .cell_data["Volume"])
    return centroids, sigma, volumes


def recover_A_from_fem(layup: Layup, refinement: str
                       ) -> Tuple[np.ndarray, Dict[str, float]]:
    """Build A^FEM from three load-case .vtkhdf files for (layup, refinement).

    Returns (A_FEM, diagnostics).
    """
    L_x, L_y, t = L_X_M, L_Y_M, layup.total_thickness_m
    interior_mask_xy = lambda c: (
        (c[:, 0] >= INTERIOR_CROP_M) & (c[:, 0] <= L_x - INTERIOR_CROP_M)
        & (c[:, 1] >= INTERIOR_CROP_M) & (c[:, 1] <= L_y - INTERIOR_CROP_M)
    )

    case_strain = {
        "case1_xx": np.array([EPS0, 0.0, 0.0]),
        "case2_yy": np.array([0.0, EPS0, 0.0]),
        "case3_xy": np.array([0.0, 0.0, GAMMA0]),
    }
    N_cols = []  # each column is N for one case
    for case in CASES:
        run_name = f"stage_09_{layup.name}_{refinement}_{case}"
        run_dir = RUNS_DIR / run_name
        # Expect a .vtkhdf converted from OpenRadioss .anim by Kitware tool.
        vtk_files = list(run_dir.glob("*.vtkhdf"))
        if not vtk_files:
            raise FileNotFoundError(
                f"no .vtkhdf in {run_dir}; run openradioss-to-vtkhdf first")
        centroids, sigma, volumes = _read_fem_stress_vtkhdf(vtk_files[0])

        # Through-thickness integration: N_ij(x,y) ~ sum_e sigma_ij^e * (V_e / A_xy_e).
        # For a structured grid each element has equal in-plane footprint
        # h_xy^2; the "column" of elements at (x,y) sums to N_ij = sum_e sigma_ij^e
        # * h_z^e. Approximate h_z^e from V_e / h_xy^2; for a uniform mesh
        # h_xy = L_x / N_x.
        h_xy = L_x / N_X
        cell_dz = volumes / (h_xy * h_xy)

        # Bin by integer (i, j) in-plane index from centroid.
        ix = np.clip(np.floor(centroids[:, 0] / h_xy).astype(int), 0, N_X - 1)
        iy = np.clip(np.floor(centroids[:, 1] / h_xy).astype(int), 0, N_Y - 1)

        # Crop interior.
        keep = interior_mask_xy(centroids)
        ix_k, iy_k = ix[keep], iy[keep]
        sigma_k = sigma[keep]
        dz_k = cell_dz[keep]

        # Sum sigma_ij * dz over each in-plane column to get N_ij(x, y).
        n_cols_x = N_X
        n_cols_y = N_Y
        Nxx = np.zeros((n_cols_x, n_cols_y))
        Nyy = np.zeros((n_cols_x, n_cols_y))
        Nxy = np.zeros((n_cols_x, n_cols_y))
        for k in range(sigma_k.shape[0]):
            i, j = ix_k[k], iy_k[k]
            Nxx[i, j] += sigma_k[k, 0, 0] * dz_k[k]
            Nyy[i, j] += sigma_k[k, 1, 1] * dz_k[k]
            Nxy[i, j] += sigma_k[k, 0, 1] * dz_k[k]

        # Mask out columns that fall outside the cropped interior in x or y.
        col_x = (np.arange(n_cols_x) + 0.5) * h_xy
        col_y = (np.arange(n_cols_y) + 0.5) * h_xy
        mx = (col_x >= INTERIOR_CROP_M) & (col_x <= L_x - INTERIOR_CROP_M)
        my = (col_y >= INTERIOR_CROP_M) & (col_y <= L_y - INTERIOR_CROP_M)
        sub_xx = Nxx[np.ix_(mx, my)]
        sub_yy = Nyy[np.ix_(mx, my)]
        sub_xy = Nxy[np.ix_(mx, my)]

        N_mean = np.array([sub_xx.mean(), sub_yy.mean(), sub_xy.mean()])
        N_std = np.array([sub_xx.std(), sub_yy.std(), sub_xy.std()])
        rel = N_std / (np.abs(N_mean) + 1e-30)
        print(f"[post] {layup.name} {refinement} {case}: N_mean = {N_mean}, "
              f"rel_std = {rel}")
        N_cols.append(N_mean)

    # Stack columns and divide by applied strain to get A^FEM.
    N_mat = np.column_stack(N_cols)                      # 3 x 3
    eps_mat = np.column_stack([case_strain[c] for c in CASES])  # 3 x 3 (diagonal)
    A_fem = N_mat @ np.linalg.inv(eps_mat)
    diagnostics = dict(applied_strains=eps_mat.tolist(),
                       N_columns=N_mat.tolist())
    return A_fem, diagnostics


def compare_A(A_fem: np.ndarray, A_clt: np.ndarray, tol: float = 0.02
              ) -> Tuple[bool, Dict[Tuple[int, int], float]]:
    """Per-component 2 percent comparison.

    For (i, j) in {(0,0), (1,1), (0,1), (1,0), (2,2)} use relative error.
    For (0, 2), (1, 2) (CLT-zero) use error normalized to max diag.
    """
    rel = {}
    ok = True
    diag_scale = max(abs(A_clt[0, 0]), abs(A_clt[1, 1]), abs(A_clt[2, 2]), 1e-30)
    for i in range(3):
        for j in range(3):
            if abs(A_clt[i, j]) > 1e-3 * diag_scale:
                e = abs(A_fem[i, j] - A_clt[i, j]) / abs(A_clt[i, j])
            else:
                e = abs(A_fem[i, j]) / diag_scale
            rel[(i, j)] = e
            if e > tol:
                ok = False
    return ok, rel


def cmd_post(args: argparse.Namespace) -> int:
    mat = IM78552()
    summary_rows = []
    summary_rows.append(["layup", "refinement", "max_rel_err",
                         "A11_FEM", "A22_FEM", "A12_FEM", "A66_FEM",
                         "A11_CLT", "A22_CLT", "A12_CLT", "A66_CLT", "PASS"])
    for layup in (LAYUP_A, LAYUP_B):
        A_clt = laminate_A_matrix(mat, layup)
        write_A_csv(A_clt, RESULTS_DIR / f"cltA_{layup.name}_reference.csv",
                    header=f"CLT A-matrix [N/m] for layup {layup.name}")
        for ref in ("ref1", "ref2"):
            try:
                A_fem, diag = recover_A_from_fem(layup, ref)
            except FileNotFoundError as exc:
                print(f"[post] skipping {layup.name} {ref}: {exc}",
                      file=sys.stderr)
                continue
            write_A_csv(A_fem, RESULTS_DIR / f"cltA_{layup.name}_{ref}_FEM.csv",
                        header=f"FEM A-matrix [N/m] for layup {layup.name}, ref {ref}")
            ok, rel = compare_A(A_fem, A_clt, tol=0.02)
            cmp_path = RESULTS_DIR / f"cltA_{layup.name}_{ref}_compare.csv"
            with cmp_path.open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["i", "j", "A_FEM", "A_CLT", "abs_err", "rel_err",
                            "tol", "pass"])
                for i in range(3):
                    for j in range(3):
                        w.writerow([i, j,
                                    f"{A_fem[i, j]:.6e}",
                                    f"{A_clt[i, j]:.6e}",
                                    f"{abs(A_fem[i, j] - A_clt[i, j]):.6e}",
                                    f"{rel[(i, j)]:.6e}",
                                    f"{0.02:.2f}",
                                    "PASS" if rel[(i, j)] <= 0.02 else "FAIL"])
            print(f"[post] {layup.name} {ref}: max rel err = {max(rel.values()):.4e}, "
                  f"{'PASS' if ok else 'FAIL'}")
            summary_rows.append([
                layup.name, ref, f"{max(rel.values()):.4e}",
                f"{A_fem[0,0]:.4e}", f"{A_fem[1,1]:.4e}",
                f"{A_fem[0,1]:.4e}", f"{A_fem[2,2]:.4e}",
                f"{A_clt[0,0]:.4e}", f"{A_clt[1,1]:.4e}",
                f"{A_clt[0,1]:.4e}", f"{A_clt[2,2]:.4e}",
                "PASS" if ok else "FAIL",
            ])
    sum_path = RESULTS_DIR / "stage_09_summary.csv"
    with sum_path.open("w", newline="") as f:
        w = csv.writer(f)
        for row in summary_rows:
            w.writerow(row)
    print(f"[post] wrote {sum_path}")
    return 0


# ---------------------------------------------------------------------------
# Sanity-check command (Section 10 last paragraph)
# ---------------------------------------------------------------------------

def cmd_check(args: argparse.Namespace) -> int:
    mat = IM78552()
    A_A = laminate_A_matrix(mat, LAYUP_A)
    A_B = laminate_A_matrix(mat, LAYUP_B)

    # Cross-ply symmetry.
    rel = lambda x, y: abs(x - y) / max(abs(x), abs(y), 1e-30)
    print(f"[check] LAYUP_A A11 = {A_A[0,0]:.6e}, A22 = {A_A[1,1]:.6e}, "
          f"rel diff = {rel(A_A[0,0], A_A[1,1]):.2e}")
    assert rel(A_A[0, 0], A_A[1, 1]) < 1e-10
    assert abs(A_A[0, 2]) < 1e-8 * abs(A_A[0, 0])
    assert abs(A_A[1, 2]) < 1e-8 * abs(A_A[0, 0])

    # Quasi-iso symmetry.
    print(f"[check] LAYUP_B A11 = {A_B[0,0]:.6e}, A22 = {A_B[1,1]:.6e}, "
          f"rel diff = {rel(A_B[0,0], A_B[1,1]):.2e}")
    assert rel(A_B[0, 0], A_B[1, 1]) < 1e-10
    assert abs(A_B[0, 2]) < 1e-8 * abs(A_B[0, 0])
    assert abs(A_B[1, 2]) < 1e-8 * abs(A_B[0, 0])
    A66_pred = (A_B[0, 0] - A_B[0, 1]) / 2.0
    print(f"[check] LAYUP_B A66 = {A_B[2,2]:.6e}, predicted (A11-A12)/2 = "
          f"{A66_pred:.6e}, rel diff = {rel(A_B[2,2], A66_pred):.2e}")
    assert rel(A_B[2, 2], A66_pred) < 1e-10

    # Apparent moduli for layup B should be in-plane isotropic.
    mods = apparent_moduli_from_A(A_B, LAYUP_B.total_thickness_m)
    print(f"[check] LAYUP_B Ex = {mods['Ex']/1e9:.4f} GPa, "
          f"Ey = {mods['Ey']/1e9:.4f} GPa, "
          f"nu_xy = {mods['nu_xy']:.4f}, "
          f"G_xy = {mods['G_xy']/1e9:.4f} GPa")
    assert rel(mods['Ex'], mods['Ey']) < 1e-10, "quasi-iso Ex != Ey"

    print("[check] all sanity checks passed")
    return 0


# ---------------------------------------------------------------------------
# CLI plumbing
# ---------------------------------------------------------------------------

def cmd_all(args: argparse.Namespace) -> int:
    rc = cmd_check(args)
    if rc:
        return rc
    rc = cmd_clt(args)
    if rc:
        return rc
    for layup_name in ("A", "B"):
        for ref in ("ref1", "ref2"):
            ns = argparse.Namespace(layup=layup_name, refinement=ref)
            rc = cmd_mesh(ns)
            if rc:
                return rc
    rc = cmd_decks(argparse.Namespace(refinement=None))
    if rc:
        return rc
    rc = cmd_run(argparse.Namespace(refinement=None))
    if rc:
        print("[all] OpenRadioss step failed; skipping post-process",
              file=sys.stderr)
        return rc
    return cmd_post(args)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("clt", help="write analytic CLT A-matrix CSVs").set_defaults(
        func=cmd_clt)

    p_mesh = sub.add_parser("mesh", help="build GMSH structured-hex mesh")
    p_mesh.add_argument("--layup", choices=("A", "B"), required=True)
    p_mesh.add_argument("--refinement", choices=tuple(REFINEMENTS),
                        default="ref1")
    p_mesh.set_defaults(func=cmd_mesh)

    p_decks = sub.add_parser("decks", help="template OpenRadioss decks")
    p_decks.add_argument("--refinement", choices=tuple(REFINEMENTS), default=None,
                         help="(optional) limit to one refinement")
    p_decks.set_defaults(func=cmd_decks)

    p_run = sub.add_parser("run", help="invoke OpenRadioss starter+engine")
    p_run.add_argument("--refinement", choices=tuple(REFINEMENTS), default=None)
    p_run.set_defaults(func=cmd_run)

    sub.add_parser("post", help="recover A^FEM and compare to A^CLT").set_defaults(
        func=cmd_post)

    sub.add_parser("check", help="sanity checks (no FEM)").set_defaults(
        func=cmd_check)

    sub.add_parser("all", help="end-to-end").set_defaults(func=cmd_all)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
