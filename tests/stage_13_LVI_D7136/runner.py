"""
Stage 13 - Drop-weight low-velocity impact (ASTM D7136) on a 24-ply
quasi-isotropic IM7/8552 laminate. Author: J.C. Vaught.

Pipeline:
  1.  Build a Gmsh-driven HEXA8 mesh of the laminate (one element per
      ply through-thickness, refined under the impact site) and a
      hemispherical impactor surface for /INTER/TYPE7 contact.
  2.  Convert the mesh to Abaqus .inp via meshio, then to OpenRadioss
      .rad via the inp2rad Python converter shipped with
      OpenRadioss/Tools.
  3.  Render the starter (_0000.rad) and engine (_0001.rad) decks from
      Jinja2 templates, including 24 /MAT/LAW25 + /FAIL/HASHIN ply
      cards, 23 /INTER/TYPE2 cohesive interfaces, /INTER/TYPE7
      impactor contact, /RBODY rigid impactor, /BCS support fixture,
      /INIVEL initial drop velocity, /GRAV body force.
  4.  Invoke OpenRadioss starter + engine through Lima + Apptainer on
      macOS or natively on Linux.
  5.  Parse T01 time history with vortex-radioss; extract peak
      contact force F_max, peak time, contact duration.
  6.  Convert .anim to .vtkhdf with the Kitware
      openradioss-to-vtkhdf tool; read the cohesive-interface damage
      scalar; project to the back face; integrate failed-tie area to
      compute the projected delamination area.
  7.  Compute Olsson 2001 large-mass closed-form peak force as a
      secondary cross-check.
  8.  Compare to Lopes-Camanho 2009 reference (peak force ~ 8.0 kN
      within 10%, projected delam area ~ 1500 mm^2 within 20%).
  9.  Verify that the stage-14 state-file handoff is on disk and
      emit state_handoff_manifest.json so stage 14's runner can
      consume the damaged coupon.

This runner does not require OpenRadioss to be installed on the host
in order to *render* decks; it only needs a working OpenRadioss
binary if --execute is passed.

References (all in /Volumes/MacShare/Code/FEA_AP-PLY/references/):
  - Lopes2009LVIPart2 (a.k.a. Lopes2009Ballistic): canonical FEM
    benchmark; comparison figures cited in spec.md section 7.1.
  - Olsson2001LargeMass:    secondary analytical peak-force check.
  - DaviesOlsson2004:       review confirming the analytical regime.
  - SodenHintonKaddour1998: material card source.
  - ASTM_D7136:             test standard.
  - Turon2007MeshSize:      cohesive mesh-size scaling rule.
"""

from __future__ import annotations

import argparse
import dataclasses as dc
import hashlib
import json
import math
import os
import platform
import shutil
import struct
import subprocess
import sys
from pathlib import Path
from typing import Iterable

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

STAGE_DIR = Path(__file__).resolve().parent
REPO_ROOT = STAGE_DIR.parent.parent
OUT_DIR = STAGE_DIR / "out"
DECK_DIR = STAGE_DIR / "deck"
MESH_DIR = STAGE_DIR / "mesh"
TEMPLATES_DIR = STAGE_DIR / "templates"

# --------------------------------------------------------------------------- #
# Geometry / problem parameters (SI units: m, kg, s, N, Pa)
# --------------------------------------------------------------------------- #

@dc.dataclass(frozen=True)
class CouponGeometry:
    length_m: float = 0.150          # 150 mm along x
    width_m: float = 0.100           # 100 mm along y
    thickness_m: float = 0.004       # 4 mm total
    n_plies: int = 24                # 24-ply quasi-iso [0/+45/-45/90]_3s
    fixture_window_x_m: float = 0.125  # 125 mm support window length
    fixture_window_y_m: float = 0.075  # 75 mm support window width
    impact_zone_x_m: float = 0.030    # 30 mm refinement zone
    impact_zone_y_m: float = 0.030

    @property
    def ply_thickness_m(self) -> float:
        return self.thickness_m / self.n_plies  # 1.667e-4 m


@dc.dataclass(frozen=True)
class Impactor:
    diameter_m: float = 0.016        # 16 mm hemispherical tip
    mass_kg: float = 5.0
    impact_energy_J: float = 30.0
    initial_gap_m: float = 5.0e-5    # 0.05 mm
    density_kg_per_m3: float = 7850.0  # steel (rigid; only used for inertia)

    @property
    def initial_velocity_mps(self) -> float:
        # v0 = sqrt(2 * E_imp / m); negative z direction at injection time.
        return math.sqrt(2.0 * self.impact_energy_J / self.mass_kg)


@dc.dataclass(frozen=True)
class IM78552Card:
    """Soden-Hinton-Kaddour 1998 card, cross-checked Lopes 2009 Tab. 1."""
    rho: float = 1580.0          # kg / m^3
    E1: float = 161.0e9          # Pa  (fibre direction)
    E2: float = 11.38e9
    E3: float = 11.38e9
    G12: float = 5.17e9
    G13: float = 5.17e9
    G23: float = 3.98e9
    nu12: float = 0.32
    nu13: float = 0.32
    nu23: float = 0.43
    Xt: float = 2806.0e6         # Pa, fibre tension
    Xc: float = 1400.0e6         # Pa, fibre compression
    Yt: float = 60.0e6           # Pa, matrix tension (transverse)
    Yc: float = 185.0e6          # Pa, matrix compression
    S12: float = 90.0e6          # Pa, in-plane shear
    S23: float = 90.0e6
    sigma_crush: float = 850.0e6  # Pa, crush
    Gf_1t: float = 81.5e3        # J / m^2  (= 81.5 N/mm)
    Gf_1c: float = 106.3e3
    Gf_2t: float = 0.28e3
    Gf_2c: float = 0.79e3


@dc.dataclass(frozen=True)
class CohesiveCard:
    """/INTER/TYPE2 inter-laminar cohesive parameters."""
    GIc: float = 250.0           # J / m^2 (= 0.25 N/mm)
    GIIc: float = 800.0          # J / m^2 (= 0.80 N/mm)
    eta_BK: float = 1.45         # Benzeggagh-Kenane mode-mixity
    Kn: float = 1.0e15           # Pa / m (= 1e6 N/mm/mm^2 = 1e6 N/mm^3)
    Ks: float = 1.0e15
    sigma_n: float = 60.0e6      # Pa
    tau_s: float = 90.0e6


@dc.dataclass(frozen=True)
class ContactCard:
    """/INTER/TYPE7 impactor-laminate contact parameters."""
    friction: float = 0.30
    stiffness_factor: float = 0.10
    gap_min_m: float = 5.0e-5
    Inacti: int = 6


@dc.dataclass(frozen=True)
class StackingSequence:
    """[0/+45/-45/90]_3s, 24 plies, ply-angles in degrees."""

    @staticmethod
    def angles_deg() -> list[int]:
        repeat = [0, 45, -45, 90]
        s_top = repeat * 3              # plies 1..12, top to mid
        s_bot = list(reversed(s_top))   # symmetric mirror, plies 13..24
        return s_top + s_bot


# --------------------------------------------------------------------------- #
# Pass criteria from Lopes-Camanho 2009 + brief
# --------------------------------------------------------------------------- #

@dc.dataclass(frozen=True)
class PassCriteria:
    F_peak_ref_N: float = 8.0e3          # Lopes 2009 Fig. 4
    F_peak_tol: float = 0.10
    t_F_peak_ref_s: float = 1.45e-3
    t_F_peak_tol: float = 0.15
    contact_duration_ref_s: float = 3.0e-3
    contact_duration_tol: float = 0.15
    delam_area_ref_m2: float = 1.5e-3    # 1500 mm^2 = 1.5e-3 m^2
    delam_area_tol: float = 0.20


# --------------------------------------------------------------------------- #
# Olsson 2001 large-mass closed-form peak force (secondary check)
# --------------------------------------------------------------------------- #

def olsson_2001_large_mass_peak_force(
    g: CouponGeometry,
    imp: Impactor,
    mat: IM78552Card,
) -> float:
    """
    Olsson 2001 large-mass impact peak force ELASTIC upper bound
    (Composites Part A, vol 32 pp 1207-1215, eq. 5-7; review by
    Davies-Olsson 2004).

    Energy balance at peak indentation (zero KE):
        E_imp = (1/2) * F^2 / k_b + (3/5) * F * (F/k_alpha)^(2/3)
                ^---- plate bending     ^---- Hertzian indentation
    Solved implicitly for F.

    k_b: static central-point stiffness of a clamped rectangular plate
    (Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells, Table
    35; for the 75 x 125 mm window, b/a = 5/3, the deflection
    coefficient is alpha_T ~ 0.0088, giving k_b = D / (alpha_T * a^2)
    where a is the short side and
       D = E_x_QI * t^3 / (12 (1 - nu^2))
    with E_x_QI ~ 60 GPa for IM7/8552 quasi-iso laminate (CLT).

    k_alpha: Hertzian indentation stiffness for a hemispherical steel
    indenter on a transverse-isotropic CFRP (Olsson 2001 eq. 6,
    k_alpha = (4/3) * E_eff * sqrt(R) with E_eff ~ E_2 of the laminate).

    IMPORTANT: this is the ELASTIC peak. The closed form has no damage.
    For a 4 mm IM7/8552 QI laminate at 30 J the elastic upper bound is
    in the 17-22 kN range, while the FEM-with-damage and the
    experiment (Lopes-Camanho 2009 Fig. 4) settle near 8 kN because
    delamination limits the load before the elastic peak is reached.
    The function is therefore reported as a SECONDARY UPPER-BOUND
    sanity check only, per Davies-Olsson 2004's own review of the
    formula's range of applicability. Olsson 2001 explicitly notes
    that the elastic peak is a "delamination-threshold load" upper
    bound; the FEM is expected to come in below it.

    Returns F_peak_elastic_upper_bound in newtons.
    """
    # Plate bending stiffness for a clamped rectangle (Timoshenko Tab. 35)
    Ex_QI = 6.0e10                         # IM7/8552 quasi-iso, ~60 GPa
    nu_QI = 0.30                           # quasi-iso effective Poisson
    D = Ex_QI * g.thickness_m**3 / (12.0 * (1.0 - nu_QI**2))
    a_short = g.fixture_window_y_m         # 0.075 m, short side
    alpha_T = 0.0088                       # Timoshenko Tab. 35 b/a~1.67
    k_b = D / (alpha_T * a_short**2)       # N / m

    # Hertzian indentation stiffness
    R = imp.diameter_m / 2.0               # 0.008 m
    E_eff = mat.E2                         # transverse modulus
    k_alpha = (4.0 / 3.0) * E_eff * math.sqrt(R)  # N / m^(3/2)

    # Implicit energy balance, bisection
    E_imp = imp.impact_energy_J

    def lhs(F: float) -> float:
        return 0.5 * F * F / k_b + 0.6 * F * (F / k_alpha) ** (2.0 / 3.0)

    lo, hi = 1.0, 1.0e7
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if lhs(mid) < E_imp:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --------------------------------------------------------------------------- #
# Deck rendering
# --------------------------------------------------------------------------- #

def _ensure_dirs() -> None:
    for d in (OUT_DIR, DECK_DIR, MESH_DIR, TEMPLATES_DIR):
        d.mkdir(parents=True, exist_ok=True)


def _write_gmsh_geo(g: CouponGeometry, imp: Impactor) -> Path:
    """Emit a Gmsh .geo script that builds the 24-ply HEXA8 stack and
    the hemispherical impactor surface mesh.

    The script uses transfinite extrusion to guarantee one element per
    ply through-thickness and a structured 1-mm in-plane mesh in the
    impact window with a graded transition outside.
    """
    geo_path = MESH_DIR / "stage_13_lvi.geo"

    # In-plane partitioning: impact window centred at origin.
    # We model the full coupon (no symmetry) per spec §4.1.
    Lx = g.length_m * 1e3   # mm for Gmsh readability
    Ly = g.width_m * 1e3
    Iz = g.impact_zone_x_m * 1e3  # 30 mm
    Iy = g.impact_zone_y_m * 1e3
    nz = g.n_plies                 # 24 layers
    h_ply_mm = g.ply_thickness_m * 1e3  # 0.1667 mm per layer

    impactor_R_mm = imp.diameter_m * 1e3 / 2.0  # 8 mm
    impactor_z0_mm = (g.thickness_m + imp.initial_gap_m) * 1e3 + impactor_R_mm

    geo = f"""// Stage 13 LVI mesh - autogenerated by runner.py
SetFactory("OpenCASCADE");
// All units mm; runner converts to m on .rad export.

cl_fine    = 1.0;
cl_medium  = 2.0;
cl_coarse  = 4.0;

// Laminate base rectangle (z=0) with concentric impact-zone subrectangle.
Point(1) = {{-{Lx/2}, -{Ly/2}, 0, cl_coarse}};
Point(2) = {{ {Lx/2}, -{Ly/2}, 0, cl_coarse}};
Point(3) = {{ {Lx/2},  {Ly/2}, 0, cl_coarse}};
Point(4) = {{-{Lx/2},  {Ly/2}, 0, cl_coarse}};

Line(1) = {{1, 2}}; Line(2) = {{2, 3}};
Line(3) = {{3, 4}}; Line(4) = {{4, 1}};
Line Loop(1) = {{1, 2, 3, 4}};

Point(5) = {{-{Iz/2}, -{Iy/2}, 0, cl_fine}};
Point(6) = {{ {Iz/2}, -{Iy/2}, 0, cl_fine}};
Point(7) = {{ {Iz/2},  {Iy/2}, 0, cl_fine}};
Point(8) = {{-{Iz/2},  {Iy/2}, 0, cl_fine}};
Line(5) = {{5, 6}}; Line(6) = {{6, 7}};
Line(7) = {{7, 8}}; Line(8) = {{8, 5}};
Line Loop(2) = {{5, 6, 7, 8}};

Plane Surface(1) = {{1, 2}};   // outer ring with hole = inner zone
Plane Surface(2) = {{2}};      // refined inner zone

// Extrude through-thickness with one HEXA8 layer per ply.
out_outer[] = Extrude {{0, 0, {g.thickness_m * 1e3}}} {{
    Surface{{1}}; Layers{{ {nz} }}; Recombine;
}};
out_inner[] = Extrude {{0, 0, {g.thickness_m * 1e3}}} {{
    Surface{{2}}; Layers{{ {nz} }}; Recombine;
}};

// Hemispherical impactor (rigid; surface only for contact).
impactor_centre_z = {impactor_z0_mm};
SetFactory("OpenCASCADE");
Sphere(100) = {{0, 0, impactor_centre_z, {impactor_R_mm}, -Pi/2, 0}};
// The hemisphere's flat face becomes the upper boundary of the rigid
// body; only the lower hemisphere contacts the laminate.

// Physical groups for OpenRadioss /SUBSET / /GRBRIC tagging.
Physical Volume("LAMINATE_INNER") = {{out_inner[1]}};
Physical Volume("LAMINATE_OUTER") = {{out_outer[1]}};
Physical Volume("IMPACTOR")       = {{100}};
Physical Surface("LAM_TOP")       = {{out_inner[0], out_outer[0]}};
Physical Surface("LAM_BOTTOM")    = {{1, 2}};

Mesh.RecombineAll = 1;
Mesh.Algorithm    = 8;     // Frontal-Delaunay for quad
Mesh.Algorithm3D  = 1;     // Delaunay
Mesh.ElementOrder = 1;     // linear, HEXA8
"""
    geo_path.write_text(geo)
    return geo_path


def _run_gmsh(geo: Path) -> Path:
    msh = geo.with_suffix(".msh")
    cmd = ["gmsh", "-3", "-format", "msh22", "-o", str(msh), str(geo)]
    print(f"[stage 13] running: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"gmsh failed:\n{res.stdout}\n{res.stderr}")
    return msh


def _msh_to_inp(msh: Path) -> Path:
    """Convert .msh to Abaqus .inp via meshio so inp2rad can chew it."""
    import meshio  # imported lazily so the renderer works without the dep

    mesh = meshio.read(str(msh))
    inp = msh.with_suffix(".inp")
    meshio.write(str(inp), mesh, file_format="abaqus")
    return inp


def _inp_to_rad(inp: Path) -> Path:
    """Run OpenRadioss inp2rad. Falls back to a stub message if the
    converter binary is not on PATH; the renderer continues."""
    rad = inp.with_suffix(".rad")
    inp2rad = shutil.which("inp2rad.py") or shutil.which("inp2rad")
    if inp2rad is None:
        print("[stage 13] WARN: inp2rad not found; mesh .rad block "
              "will be templated with a placeholder. The user must "
              "run inp2rad inside the Lima VM before --execute.")
        rad.write_text("# placeholder for inp2rad output\n")
        return rad
    cmd = [inp2rad, str(inp), "-o", str(rad)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"inp2rad failed:\n{res.stdout}\n{res.stderr}")
    return rad


def _render_starter(
    g: CouponGeometry,
    imp: Impactor,
    mat: IM78552Card,
    coh: CohesiveCard,
    contact: ContactCard,
    mesh_rad_block: str,
) -> Path:
    """Render the OpenRadioss starter deck (_0000.rad).

    A concise, human-readable, fully-templated text deck.  Real
    production workflow would use Jinja2; we keep it as f-strings to
    avoid an extra dep.
    """
    angles = StackingSequence.angles_deg()  # 24-element list
    unique_angles = sorted(set(angles))     # [-45, 0, 45, 90]

    # /SKEW cards: one per unique ply angle, rotated about z-axis.
    skew_cards = []
    for k, ang in enumerate(unique_angles, start=1):
        c = math.cos(math.radians(ang))
        s = math.sin(math.radians(ang))
        skew_cards.append(
            f"/SKEW/FIX/{k}\nply_{ang:+d}_skew\n"
            f"# OX OY OZ X1 Y1 Z1 X2 Y2 Z2\n"
            f"  0.0 0.0 0.0  {c} {s} 0.0  {-s} {c} 0.0\n"
        )
    skew_block = "\n".join(skew_cards)
    angle_to_skew = {ang: k for k, ang in enumerate(unique_angles, start=1)}

    # /MAT/LAW25 + /FAIL/HASHIN: one per unique ply angle (material is
    # identical; only orientation changes via /SKEW, but OpenRadioss
    # convention is to bind material -> property -> skew so we replicate
    # the material 4 times and let the property carry /SKEW).
    mat_cards = []
    for k, ang in enumerate(unique_angles, start=1):
        mat_cards.append(f"""/MAT/LAW25/{k}
IM7_8552_{ang:+d}deg
# rho_I             E11           E22           NU12          NU23
  {mat.rho:.3e}     {mat.E1:.3e}  {mat.E2:.3e}  {mat.nu12}    {mat.nu23}
# G12               G23           EPS_F1        EPS_F2        EPS_M
  {mat.G12:.3e}     {mat.G23:.3e} 0.0           0.0           0.0
# Iform Sigy   B    N    FMAX
  1     0.0    0.0  0.0  0.0
# WPLAREF WPMAX Cbeta_T Cbeta_C Csigma
  0.0     0.0   0.0     0.0     0.0
/FAIL/HASHIN/{k}
# Ifail_sh  Ifail_so   ratio_max
  0          2          1.0
# sigma1t        sigma1c        sigma2t        sigma2c        sigma_crush
  {mat.Xt:.3e}   {mat.Xc:.3e}   {mat.Yt:.3e}   {mat.Yc:.3e}   {mat.sigma_crush:.3e}
# S12            S23
  {mat.S12:.3e}  {mat.S23:.3e}
# Gf_1t          Gf_1c          Gf_2t          Gf_2c
  {mat.Gf_1t:.3e} {mat.Gf_1c:.3e} {mat.Gf_2t:.3e} {mat.Gf_2c:.3e}
""")
    mat_block = "\n".join(mat_cards)

    # /PROP/TYPE14 + /PART per ply (24 plies; mat & skew chosen by angle).
    prop_part_cards = []
    for ply_idx, ang in enumerate(angles, start=1):
        mat_id = angle_to_skew[ang]
        skew_id = angle_to_skew[ang]
        prop_id = 100 + ply_idx
        part_id = ply_idx
        grbric_id = 1000 + ply_idx
        prop_part_cards.append(f"""/PROP/TYPE14/{prop_id}
ply_{ply_idx:02d}_solid_prop
# Ihbe Ismstr Iframe Iorth Iplas Itetra Itetra4 Iint
  24   2      0      1     0     0       0       0
# qa qb h dn LFAC SKEW_ID
  1.1 0.05 0.0 0.0 0.0 {skew_id}
/PART/{part_id}
ply_{ply_idx:02d}_part
# prop_ID  mat_ID  subset_ID
  {prop_id} {mat_id} 0
# group of bricks for this ply (set by mesh extrusion layer ply_idx)
/GRBRIC/PART/{grbric_id}
ply_{ply_idx:02d}_bricks
  {part_id}
""")
    prop_part_block = "\n".join(prop_part_cards)

    # 23 inter-ply /INTER/TYPE2 pairs.
    inter_type2_cards = []
    for k in range(1, g.n_plies):  # interfaces 1..23
        ml = 2000 + k     # main surface (ply k bottom)
        sl = 3000 + k     # secondary surface (ply k+1 top)
        iid = 7000 + k
        inter_type2_cards.append(f"""/INTER/TYPE2/{iid}
ply_{k:02d}_to_ply_{k+1:02d}_cohesive_tie
# Spotflag Idel Igap Inacti Iform Imult
  25       1    0    0      2     0
# Stmin Stfac Fric  Gap   Tstart Tstop
  0.0   1.0   0.0   0.0   0.0    0.0
# main_grnod (ply k bottom face nodes)  / secondary_grnod (ply k+1 top face nodes)
  {ml}                                   {sl}
# bilinear cohesive failure block (energy-based, Spotflag=25):
# alpha_n  alpha_t  N  Imass  sigma_n        sigma_t        GIc            GIIc        eta_BK
  1.0      1.0      1  0      {coh.sigma_n:.3e} {coh.tau_s:.3e} {coh.GIc:.3e} {coh.GIIc:.3e} {coh.eta_BK}
""")
    inter_type2_block = "\n".join(inter_type2_cards)

    # Impactor /RBODY + /MAT/LAW1 (rigid steel just for inertia consistency)
    impactor_block = f"""/MAT/LAW1/100
steel_rigid_impactor_dummy
# rho E nu
  {imp.density_kg_per_m3} 2.10e11 0.30
/PROP/TYPE14/200
impactor_solid_prop
  24 2 0 0 0 0 0 0
  1.1 0.05 0.0 0.0 0.0 0
/PART/100
impactor_part
  200 100 0
/RBODY/9001
hemispherical_impactor_rigid
# main_node ICOG Surf_ID Skew Ispher
  99001     1    0       0    0
# Mass        Jxx        Jyy         Jzz         Jxy Jxz Jyz
  {imp.mass_kg} 4.6e-4    4.6e-4      4.6e-4      0.0 0.0 0.0
# secondary_grnod (impactor brick nodes)
  4001
/INIVEL/TRA/1
impactor_drop
  0.0  0.0  -{imp.initial_velocity_mps:.6f}
# applied to grnod
  4001
/GRAV/1
gravity
# Ngrav Ndir Tstart Skew_ID
  -9.81 3    0.0    0
# applied to grnod (everyone)
  9999
"""

    # Impactor-laminate contact /INTER/TYPE7
    contact_block = f"""/INTER/TYPE7/8001
impactor_to_laminate_contact
# Igap Multimp Inacti Fpenmax Iform Idel sens_id
  0    0       {contact.Inacti}      0.0    2     0    0
# Stmin Stmax dtmin dtstart Tstart Tstop
  0.0   0.0   0.0   0.0     0.0    0.0
# main_surf_id (impactor lower hemisphere) / secondary_grnod (laminate upper face)
  4002                                       4003
# Stfac fric gapmin Iform2 sens_id
  {contact.stiffness_factor} {contact.friction} {contact.gap_min_m} 0 0
"""

    # Support BC /BCS — uz=0 on lower face of laminate outside the
    # 75 x 125 mm window (Lopes-Camanho 2009 §3.1 simplification).
    bcs_block = f"""/BCS/1
support_frame_uz_clamped
# trans_x trans_y trans_z rot_x rot_y rot_z
  0       0       1       0     0     0
# skew_id grnod_id
  0       5001
# (grnod 5001 = laminate bottom-face nodes outside fixture window)
"""

    # Output requests for stage 14 chaining (see spec §9.3).
    state_block = f"""/STATE/BRICK/FULL
# Tstart Tstop Tdt
  4.5e-3 5.0e-3 5.0e-4
/STATE/INTER/FULL
# Tstart Tstop Tdt
  4.5e-3 5.0e-3 5.0e-4
"""

    th_block = """/TH/RBODY/1
impactor_kinematics_TH
# var_list
  X Y Z VX VY VZ AX AY AZ FX FY FZ
  9001
/TH/INTER/2
contact_force_TH
# var_list
  FX FY FZ E
  8001
/TH/PART/3
energy_TH
# var_list
  IE KE HE
  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 100
"""

    deck = f"""# OpenRadioss starter deck for Stage 13 LVI (ASTM D7136)
# Generated by runner.py from FEA_AP-PLY repo
# Author: J.C. Vaught
/BEGIN
stage_13_LVI_D7136
       0                     1
                              1.0       1.0        1.0        1.0
/UNIT/1
m  kg  s

# ---------------- mesh imported from inp2rad ----------------
{mesh_rad_block}

# ---------------- ply orientations as /SKEW ----------------
{skew_block}

# ---------------- materials and Hashin failure ----------------
{mat_block}

# ---------------- properties and parts (24 plies) ----------------
{prop_part_block}

# ---------------- impactor rigid body ----------------
{impactor_block}

# ---------------- impactor-laminate contact ----------------
{contact_block}

# ---------------- inter-laminar cohesive interfaces (23) ----------------
{inter_type2_block}

# ---------------- support fixture BCs ----------------
{bcs_block}

# ---------------- state file output for stage 14 ----------------
{state_block}

# ---------------- time history requests ----------------
{th_block}

/END
"""
    p = DECK_DIR / "stage_13_LVI_D7136_0000.rad"
    p.write_text(deck)
    return p


def _render_engine() -> Path:
    deck = """# OpenRadioss engine deck for Stage 13 LVI (ASTM D7136)
/RUN/stage_13_LVI_D7136/1
        5.0e-3        0
/PRINT/-100
/STOP
/TFILE/1
        1.0e-6
/ANIM/DT
        0.0       2.5e-5
/ANIM/BRICK/STRESS/ALL
/ANIM/BRICK/DAMA
/ANIM/BRICK/EPSP
/ANIM/INTER
/H3D/DT
        0.0       2.5e-5
/H3D/COMPRESS
/DT/BRICK/CST
        0.667     5.0e-9
/INISTA/STATE_OUT
        4.5e-3
/END
"""
    p = DECK_DIR / "stage_13_LVI_D7136_0001.rad"
    p.write_text(deck)
    return p


# --------------------------------------------------------------------------- #
# OpenRadioss invocation (Lima + Apptainer on macOS, native on Linux)
# --------------------------------------------------------------------------- #

def _is_macos() -> bool:
    return platform.system() == "Darwin"


def _invoke_starter(starter_deck: Path) -> int:
    if _is_macos():
        cmd = [
            "limactl", "shell", "apptainer", "--",
            "/OpenRadioss/exec/starter_linuxa64",
            "-i", str(starter_deck), "-nt", "1",
        ]
    else:
        cmd = [
            "starter_linux64_gf", "-i", str(starter_deck), "-nt", "1",
        ]
    print(f"[stage 13] starter: {' '.join(cmd)}")
    return subprocess.call(cmd, cwd=DECK_DIR)


def _invoke_engine(engine_deck: Path, n_proc: int = 8) -> int:
    if _is_macos():
        cmd = [
            "limactl", "shell", "apptainer", "--",
            "mpirun", "-np", str(n_proc),
            "/OpenRadioss/exec/engine_linuxa64_ompi",
            "-i", str(engine_deck),
        ]
    else:
        cmd = [
            "mpirun", "-np", str(n_proc),
            "engine_linux64_gf_ompi", "-i", str(engine_deck),
        ]
    print(f"[stage 13] engine: {' '.join(cmd)}")
    return subprocess.call(cmd, cwd=DECK_DIR)


# --------------------------------------------------------------------------- #
# T01 time-history parsing (vortex-radioss preferred; fallback to ASCII)
# --------------------------------------------------------------------------- #

def _parse_t01(t01_path: Path) -> dict[str, list[tuple[float, float]]]:
    """Parse the T01 binary time-history file. Tries vortex-radioss
    first (the maintained Python tool, audit row B); falls back to
    OpenRadioss's bundled `T01ascii` text dump if vortex-radioss is
    not installed."""
    try:
        from vortex_radioss import T01Reader  # type: ignore
        rdr = T01Reader(str(t01_path))
        out: dict[str, list[tuple[float, float]]] = {}
        for series in rdr.list_series():
            ts = rdr.times()
            ys = rdr.values(series)
            out[series] = list(zip(ts, ys))
        return out
    except Exception:
        pass

    # Fallback: invoke T01ascii then parse the ASCII columns.
    ascii_path = t01_path.with_suffix(".T01.txt")
    t01ascii = shutil.which("T01ascii")
    if t01ascii:
        subprocess.run([t01ascii, str(t01_path), str(ascii_path)],
                       check=False, capture_output=True)
    if not ascii_path.exists():
        raise RuntimeError(
            f"Cannot parse {t01_path}: neither vortex-radioss nor "
            "T01ascii available."
        )

    out = {}
    current = None
    times: list[float] = []
    values: list[float] = []
    for line in ascii_path.read_text().splitlines():
        if line.startswith("#") and "series:" in line:
            if current is not None:
                out[current] = list(zip(times, values))
            current = line.split("series:")[-1].strip()
            times, values = [], []
            continue
        toks = line.split()
        if len(toks) >= 2:
            try:
                t = float(toks[0])
                v = float(toks[1])
                times.append(t)
                values.append(v)
            except ValueError:
                continue
    if current is not None:
        out[current] = list(zip(times, values))
    return out


# --------------------------------------------------------------------------- #
# Force-time analysis
# --------------------------------------------------------------------------- #

def _resultant_force(series: dict[str, list[tuple[float, float]]],
                     prefix: str) -> list[tuple[float, float]]:
    """Combine FX, FY, FZ components of a /TH/INTER block into the
    Euclidean resultant time-history."""
    fx = dict(series.get(f"{prefix}_FX", []))
    fy = dict(series.get(f"{prefix}_FY", []))
    fz = dict(series.get(f"{prefix}_FZ", []))
    ts = sorted(set(fx) | set(fy) | set(fz))
    out = []
    for t in ts:
        ax = fx.get(t, 0.0)
        ay = fy.get(t, 0.0)
        az = fz.get(t, 0.0)
        out.append((t, math.sqrt(ax * ax + ay * ay + az * az)))
    return out


def _analyze_force_time(force_t: list[tuple[float, float]]) -> dict:
    if not force_t:
        return dict(F_peak_N=float("nan"), t_peak_s=float("nan"),
                    contact_duration_s=float("nan"))
    ts, fs = zip(*force_t)
    F_peak = max(fs)
    t_peak = ts[fs.index(F_peak)]
    # Contact duration: continuous run of F > 0.05 * F_peak around t_peak.
    threshold = 0.05 * F_peak
    in_contact = [(t, f) for t, f in force_t if f > threshold]
    if in_contact:
        contact_duration = in_contact[-1][0] - in_contact[0][0]
    else:
        contact_duration = float("nan")
    return dict(F_peak_N=F_peak, t_peak_s=t_peak,
                contact_duration_s=contact_duration)


# --------------------------------------------------------------------------- #
# Delamination area extraction
# --------------------------------------------------------------------------- #

def _anim_to_vtkhdf(anim_dir: Path) -> Path:
    """Convert the most-recent .anim frame to .vtkhdf using the
    Kitware openradioss-to-vtkhdf converter (audit row B)."""
    converter = (shutil.which("openradioss-to-vtkhdf")
                 or shutil.which("openradioss_to_vtkhdf.py"))
    if converter is None:
        raise RuntimeError(
            "openradioss-to-vtkhdf is not on PATH; install per "
            "audit row B (gitlab.kitware.com/keu-public/openradioss-to-vtkhdf)."
        )
    last = sorted(anim_dir.glob("*A0*"))[-1]
    out = anim_dir / (last.name + ".vtkhdf")
    cmd = [converter, str(last), "-o", str(out)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"openradioss-to-vtkhdf failed:\n"
                           f"{res.stdout}\n{res.stderr}")
    return out


def _projected_delamination_area_m2(vtkhdf: Path) -> tuple[float, list[dict]]:
    """Read the cohesive-interface damage scalar from the .vtkhdf and
    integrate the projected (xy-plane) area of failed ties on each
    interface, returning (total_back_face_projected_area_m2,
    per_interface_breakdown).

    "Failed" = damage scalar >= 0.99 on the cohesive surface element.
    """
    try:
        import pyvista as pv  # type: ignore
    except ImportError as e:
        raise RuntimeError("pyvista is required for delamination "
                           "extraction; pip install pyvista") from e

    grid = pv.read(str(vtkhdf))
    blocks = grid if hasattr(grid, "n_blocks") else [grid]

    per_interface = []
    total_area_m2 = 0.0
    for block in blocks:
        if not hasattr(block, "cell_data"):
            continue
        if "INTER_DAM" not in block.cell_data and "DAMA" not in block.cell_data:
            continue
        scalar = (block.cell_data.get("INTER_DAM")
                  or block.cell_data["DAMA"])
        # Only cohesive-tie cells will carry INTER_DAM in OpenRadioss .vtkhdf.
        failed_mask = scalar >= 0.99
        if not failed_mask.any():
            continue
        # Project to xy-plane: per-cell area times indicator.
        sizes = block.compute_cell_sizes(area=True, volume=False)
        cell_areas = sizes.cell_data["Area"]  # m^2
        block_area = float((cell_areas * failed_mask).sum())
        per_interface.append(dict(
            block_name=getattr(block, "name", "?"),
            failed_cells=int(failed_mask.sum()),
            total_cells=int(scalar.size),
            area_m2=block_area,
        ))
        total_area_m2 += block_area
    return total_area_m2, per_interface


# --------------------------------------------------------------------------- #
# Pass / fail report
# --------------------------------------------------------------------------- #

def _within_tol(measured: float, ref: float, tol: float) -> tuple[bool, float]:
    if math.isnan(measured) or math.isnan(ref) or ref == 0:
        return False, float("nan")
    rel = abs(measured - ref) / abs(ref)
    return rel <= tol, rel


def _emit_pass_report(
    F_peak_N: float, t_peak_s: float, contact_duration_s: float,
    delam_area_m2: float, F_olsson_N: float,
    crit: PassCriteria,
) -> dict:
    pass_F, rel_F = _within_tol(F_peak_N, crit.F_peak_ref_N, crit.F_peak_tol)
    pass_t, rel_t = _within_tol(t_peak_s, crit.t_F_peak_ref_s,
                                crit.t_F_peak_tol)
    pass_dur, rel_dur = _within_tol(contact_duration_s,
                                    crit.contact_duration_ref_s,
                                    crit.contact_duration_tol)
    pass_dam, rel_dam = _within_tol(delam_area_m2, crit.delam_area_ref_m2,
                                    crit.delam_area_tol)
    # Olsson 2001 elastic upper bound: FEM peak must be BELOW this
    # (one-sided inequality, not a percent tolerance). If F_FEM > F_Olsson
    # the FEM is unphysical because the elastic energy bound is violated.
    pass_olsson_bound = (not math.isnan(F_peak_N)
                         and not math.isnan(F_olsson_N)
                         and F_peak_N < F_olsson_N)
    return dict(
        peak_force_N=dict(measured=F_peak_N,
                          ref_lopes_camanho_2009=crit.F_peak_ref_N,
                          olsson_2001_elastic_upper_bound=F_olsson_N,
                          fem_below_olsson_bound=pass_olsson_bound,
                          tolerance=crit.F_peak_tol,
                          rel_error_vs_lopes=rel_F, pass_=pass_F),
        peak_time_s=dict(measured=t_peak_s,
                         ref=crit.t_F_peak_ref_s,
                         tolerance=crit.t_F_peak_tol,
                         rel_error=rel_t, pass_=pass_t),
        contact_duration_s=dict(measured=contact_duration_s,
                                ref=crit.contact_duration_ref_s,
                                tolerance=crit.contact_duration_tol,
                                rel_error=rel_dur, pass_=pass_dur),
        delam_area_m2=dict(measured=delam_area_m2,
                           ref_lopes_camanho_2009=crit.delam_area_ref_m2,
                           tolerance=crit.delam_area_tol,
                           rel_error=rel_dam, pass_=pass_dam),
        overall_pass=all([pass_F, pass_dam, pass_olsson_bound]),
        primary_criteria="peak force within 10% of Lopes-Camanho 2009;"
                         " projected delam area within 20%;"
                         " FEM peak below Olsson 2001 elastic upper bound",
        secondary_criteria="peak time and contact duration"
                           " (informational, 15% tolerance)",
    )


# --------------------------------------------------------------------------- #
# Stage 14 handoff
# --------------------------------------------------------------------------- #

def _sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _emit_stage14_manifest(starter_deck: Path) -> Path:
    """Write state_handoff_manifest.json so Stage 14's runner can
    consume the damaged-coupon state file (see spec §9.3)."""
    sta = OUT_DIR / "stage_13_LVI_D7136_S0001.sta"
    bin_ = OUT_DIR / "stage_13_LVI_D7136_S0001"
    artefacts = []
    for p in (sta, bin_, starter_deck):
        if p.exists():
            artefacts.append(dict(
                role={
                    "stage_13_LVI_D7136_S0001.sta": "ascii_state_fragment",
                    "stage_13_LVI_D7136_S0001":     "binary_restart",
                    "stage_13_LVI_D7136_0000.rad":  "mesh_and_ids",
                }.get(p.name, "unknown"),
                path=str(p), sha256=_sha256_of(p),
                size_bytes=p.stat().st_size,
            ))
    manifest = dict(
        producer="tests/stage_13_LVI_D7136/runner.py",
        spec="tests/stage_13_LVI_D7136/spec.md",
        state_capture_time_s=4.5e-3,
        consumer="tests/stage_14_CAI_D7137/runner.py (stage 14)",
        notes="Stage 14 must /INISTA-load these into a fresh starter "
              "deck after replacing impactor and contact cards with "
              "the D7137 anti-buckling fixture and compressive load "
              "as per master_plan.md row 14.",
        artefacts=artefacts,
    )
    out = OUT_DIR / "state_handoff_manifest.json"
    out.write_text(json.dumps(manifest, indent=2))
    return out


# --------------------------------------------------------------------------- #
# Top-level entry point
# --------------------------------------------------------------------------- #

def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Stage 13 (D7136 LVI) runner.")
    parser.add_argument("--render-only", action="store_true",
                        help="Build mesh and render decks; do not run.")
    parser.add_argument("--execute", action="store_true",
                        help="Invoke OpenRadioss starter+engine.")
    parser.add_argument("--analyze", action="store_true",
                        help="Parse T01 + .anim and emit pass report.")
    parser.add_argument("--all", action="store_true",
                        help="Render, execute, and analyze.")
    parser.add_argument("--n-proc", type=int, default=8)
    args = parser.parse_args(list(argv) if argv else None)

    if not (args.render_only or args.execute or args.analyze or args.all):
        args.all = True

    _ensure_dirs()

    g = CouponGeometry()
    imp = Impactor()
    mat = IM78552Card()
    coh = CohesiveCard()
    contact = ContactCard()
    crit = PassCriteria()

    if args.render_only or args.all:
        geo = _write_gmsh_geo(g, imp)
        try:
            msh = _run_gmsh(geo)
            inp = _msh_to_inp(msh)
            rad_mesh = _inp_to_rad(inp)
            mesh_block = rad_mesh.read_text()
        except Exception as e:
            print(f"[stage 13] mesh pipeline incomplete ({e}); "
                  "rendering decks with placeholder mesh block.")
            mesh_block = ("# WARNING: mesh block placeholder — "
                          "rerun with gmsh + meshio + inp2rad available.")
        starter = _render_starter(g, imp, mat, coh, contact, mesh_block)
        engine = _render_engine()
        F_olsson = olsson_2001_large_mass_peak_force(g, imp, mat)
        print(f"[stage 13] starter -> {starter}")
        print(f"[stage 13] engine  -> {engine}")
        print(f"[stage 13] Olsson 2001 large-mass peak: "
              f"F_olsson = {F_olsson:.3e} N "
              f"({F_olsson * 1e-3:.2f} kN)")
        print(f"[stage 13] drop velocity v0 = "
              f"{imp.initial_velocity_mps:.4f} m/s "
              f"({imp.initial_velocity_mps:.4f} m/s)")

    if args.execute or args.all:
        starter_deck = DECK_DIR / "stage_13_LVI_D7136_0000.rad"
        engine_deck = DECK_DIR / "stage_13_LVI_D7136_0001.rad"
        rc = _invoke_starter(starter_deck)
        if rc != 0:
            print(f"[stage 13] starter exited {rc}", file=sys.stderr)
            return rc
        rc = _invoke_engine(engine_deck, n_proc=args.n_proc)
        if rc != 0:
            print(f"[stage 13] engine exited {rc}", file=sys.stderr)
            return rc

    if args.analyze or args.all:
        t01_candidates = sorted(DECK_DIR.glob("*T01*"))
        if not t01_candidates:
            print("[stage 13] no T01 found; skipping analysis.")
            return 1
        t01 = t01_candidates[-1]
        series = _parse_t01(t01)
        force_t = _resultant_force(series, "INTER_8001")

        # CSV: contact force vs. time.
        force_csv = OUT_DIR / "force_time.csv"
        with force_csv.open("w") as f:
            f.write("t_s,F_contact_N\n")
            for t, ff in force_t:
                f.write(f"{t:.6e},{ff:.6e}\n")

        analysis = _analyze_force_time(force_t)
        F_olsson = olsson_2001_large_mass_peak_force(g, imp, mat)

        # Delamination area from .anim -> .vtkhdf
        try:
            vtkhdf = _anim_to_vtkhdf(DECK_DIR)
            delam_area_m2, per_iface = _projected_delamination_area_m2(vtkhdf)
        except Exception as e:
            print(f"[stage 13] delamination extraction failed: {e}")
            delam_area_m2 = float("nan")
            per_iface = []

        # CSV: per-interface area
        delam_csv = OUT_DIR / "delam_area_per_interface.csv"
        with delam_csv.open("w") as f:
            f.write("interface_id,failed_cells,total_cells,area_m2\n")
            for k, info in enumerate(per_iface):
                f.write(f"{k},{info['failed_cells']},"
                        f"{info['total_cells']},{info['area_m2']:.6e}\n")

        report = _emit_pass_report(
            F_peak_N=analysis["F_peak_N"],
            t_peak_s=analysis["t_peak_s"],
            contact_duration_s=analysis["contact_duration_s"],
            delam_area_m2=delam_area_m2,
            F_olsson_N=F_olsson,
            crit=crit,
        )
        (OUT_DIR / "pass_report.json").write_text(
            json.dumps(report, indent=2))

        # Manifest for stage 14 handoff
        starter_deck = DECK_DIR / "stage_13_LVI_D7136_0000.rad"
        manifest_path = _emit_stage14_manifest(starter_deck)

        print(f"[stage 13] pass report -> {OUT_DIR / 'pass_report.json'}")
        print(f"[stage 13] handoff manifest -> {manifest_path}")
        print(json.dumps(report, indent=2))
        return 0 if report["overall_pass"] else 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
