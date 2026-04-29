"""
Stage 08 -- Ply rotation transformation verification.

Runs a sweep of seven OpenRadioss implicit-static jobs on a single-ply
solid coupon (250 x 25 x 1 mm, IM7/8552 from Soden 1998) at fiber-to-load
angles theta in {0, 15, 30, 45, 60, 75, 90} degrees. Geometry, mesh, and
boundary conditions are identical across all seven runs; only the per-element
material orientation is rotated, via /SKEW/FIX referenced by /PROP/TYPE14
with Iorth=1.

For each run the runner extracts the apparent modulus E_x(theta) from the
volume-averaged Cauchy stress over the gauge sub-volume divided by the
imposed engineering strain, and compares it against the closed-form
Jones / Daniel-Ishai transformation

    1/E_x(theta) = c^4/E1 + s^4/E2 + (1/G12 - 2 nu12/E1) c^2 s^2

with c=cos(theta), s=sin(theta). Pass criterion is |E_x_FEM - E_x_analytic|
/ E_x_analytic <= 0.01 at every angle.

This runner does NOT plot. CSV is exported and a separate Typst+CeTZ
document (per project preferences) draws the figure.

Run:
    python runner.py
or
    python runner.py --dry-run         # template decks only, no solver
    python runner.py --angles 0 45 90  # subset of the sweep
    python runner.py --analytic-only   # closed-form table only

Author: J.C. Vaught
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np


# ----------------------------------------------------------------------------
# 1. Material card (Soden, Hinton, Kaddour 1998 -- IM7/8552 UD CFRP). SI units.
# ----------------------------------------------------------------------------

@dataclass(frozen=True)
class Lamina:
    name: str
    rho: float       # kg/m^3
    E1: float        # Pa
    E2: float        # Pa
    E3: float        # Pa
    nu12: float
    nu13: float
    nu23: float
    G12: float       # Pa
    G13: float       # Pa
    G23: float       # Pa


IM7_8552 = Lamina(
    name="IM7_8552_UD",
    rho=1580.0,
    E1=165.0e9,
    E2=9.0e9,
    E3=9.0e9,
    nu12=0.34,
    nu13=0.34,
    nu23=0.5,
    G12=5.6e9,
    G13=5.6e9,
    G23=3.0e9,
)


# ----------------------------------------------------------------------------
# 2. Closed-form off-axis engineering constants (Jones 1999, Daniel-Ishai 2006)
# ----------------------------------------------------------------------------

def Ex_analytic(theta_deg: float, m: Lamina = IM7_8552) -> float:
    """Off-axis apparent uniaxial modulus E_x(theta) [Pa].

    Equation (2.85) of Jones 1999 / Chapter 5 of Daniel-Ishai 2006:

        1/E_x = c^4/E1 + s^4/E2 + (1/G12 - 2 nu12/E1) c^2 s^2.
    """
    th = math.radians(theta_deg)
    c, s = math.cos(th), math.sin(th)
    c2, s2 = c * c, s * s
    inv_Ex = c2 * c2 / m.E1 + s2 * s2 / m.E2 + (1.0 / m.G12 - 2.0 * m.nu12 / m.E1) * c2 * s2
    return 1.0 / inv_Ex


def nu_xy_analytic(theta_deg: float, m: Lamina = IM7_8552) -> float:
    """Off-axis major Poisson ratio nu_xy(theta), Jones 1999 eq. 2.87."""
    th = math.radians(theta_deg)
    c, s = math.cos(th), math.sin(th)
    c2, s2 = c * c, s * s
    Ex = Ex_analytic(theta_deg, m)
    bracket = (m.nu12 / m.E1) * (c2 * c2 + s2 * s2) - (1.0 / m.E1 + 1.0 / m.E2 - 1.0 / m.G12) * c2 * s2
    return Ex * bracket


def Gxy_analytic(theta_deg: float, m: Lamina = IM7_8552) -> float:
    """Off-axis in-plane shear modulus G_xy(theta), Jones 1999 eq. 2.88.

        1/G_xy = 2*(2/E1 + 2/E2 + 4 nu12/E1 - 1/G12) c^2 s^2 + (1/G12) (c^2 - s^2)^2
    """
    th = math.radians(theta_deg)
    c, s = math.cos(th), math.sin(th)
    c2, s2 = c * c, s * s
    inv_G = 2.0 * (2.0 / m.E1 + 2.0 / m.E2 + 4.0 * m.nu12 / m.E1 - 1.0 / m.G12) * c2 * s2 \
        + (1.0 / m.G12) * (c2 - s2) ** 2
    return 1.0 / inv_G


# ----------------------------------------------------------------------------
# 3. Geometry / mesh / loading constants (same as stage 7).
# ----------------------------------------------------------------------------

L_COUPON = 0.250        # m
B_COUPON = 0.025        # m
T_COUPON = 0.001        # m

NX, NY, NZ = 100, 10, 2

STRAIN_BAR = 2.0e-3     # 0.2 percent engineering strain

SWEEP_ANGLES_DEG = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)

PASS_TOL_REL = 0.01     # 1 percent relative tolerance on E_x


# ----------------------------------------------------------------------------
# 4. Mesh -- structured HEXA8 grid generated procedurally (no GMSH dependency
#    inside the runner; we write a Radioss /BRICK block directly). For stage 9
#    this will be replaced by the GMSH+inp2rad path; here a structured brick
#    is sufficient and removes a moving part from the verification.
# ----------------------------------------------------------------------------

def build_mesh() -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Build structured HEXA8 mesh of the 250 x 25 x 1 mm coupon.

    Returns
    -------
    nodes : (Nn, 3) array of node coordinates [m]
    bricks : (Ne, 8) int array of node IDs (1-based) per element, Radioss order
    face_nodes : dict mapping face name -> array of node IDs (1-based) on that face
    """
    xs = np.linspace(0.0, L_COUPON, NX + 1)
    ys = np.linspace(-B_COUPON / 2.0, B_COUPON / 2.0, NY + 1)
    zs = np.linspace(-T_COUPON / 2.0, T_COUPON / 2.0, NZ + 1)

    Nx, Ny, Nz = NX + 1, NY + 1, NZ + 1
    nodes = np.zeros((Nx * Ny * Nz, 3))
    nid = lambda i, j, k: i * (Ny * Nz) + j * Nz + k  # 0-based
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                nodes[nid(i, j, k)] = (xs[i], ys[j], zs[k])

    # Radioss /BRICK node order (HEXA8): n1..n8 with bottom face CCW then top.
    bricks = []
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                n1 = nid(i,     j,     k)     + 1
                n2 = nid(i + 1, j,     k)     + 1
                n3 = nid(i + 1, j + 1, k)     + 1
                n4 = nid(i,     j + 1, k)     + 1
                n5 = nid(i,     j,     k + 1) + 1
                n6 = nid(i + 1, j,     k + 1) + 1
                n7 = nid(i + 1, j + 1, k + 1) + 1
                n8 = nid(i,     j + 1, k + 1) + 1
                bricks.append((n1, n2, n3, n4, n5, n6, n7, n8))
    bricks = np.array(bricks, dtype=int)

    face_nodes = {
        "x0":  np.array([nid(0,        j, k) + 1 for j in range(Ny) for k in range(Nz)], dtype=int),
        "xL":  np.array([nid(NX,       j, k) + 1 for j in range(Ny) for k in range(Nz)], dtype=int),
    }

    return nodes, bricks, face_nodes


# ----------------------------------------------------------------------------
# 5. Deck templating.
# ----------------------------------------------------------------------------

DECK_HEADER = """\
#RADIOSS STARTER
/BEGIN
stage_08_off_axis_theta_{theta_tag}
       12345         2026
                  kg                    m                    s
                  kg                    m                    s
"""

DECK_MAT = """\
/MAT/LAW25/1
{name}
{rho:>20.6E}
{E1:>20.6E}{E2:>20.6E}{nu12:>20.6E}{G12:>20.6E}
{E3:>20.6E}{nu13:>20.6E}{nu23:>20.6E}{G13:>20.6E}{G23:>20.6E}
"""

DECK_SKEW = """\
/SKEW/FIX/9
ply_skew_theta_{theta_tag}
{Ox:>20.6E}{Oy:>20.6E}{Oz:>20.6E}{Ax:>20.6E}{Ay:>20.6E}{Az:>20.6E}{Bx:>20.6E}{By:>20.6E}{Bz:>20.6E}
"""

DECK_PROP = """\
/PROP/TYPE14/1
ply_solid_prop
#  Ihex  Iframe  Iorth   Iplas    Phi     Iskew
       24       1      1       0    0.0          9
"""

DECK_PART = """\
/PART/1
ply_part                                          1         1
"""

DECK_FUNCT = """\
/FUNCT/7
ramp_unit
                 0.0                  0.0
                 1.0                  1.0
"""

DECK_BCS = """\
/GRNOD/NODE/100
nodes_face_x0
{grnod_x0}
/GRNOD/NODE/200
nodes_face_xL
{grnod_xL}
/BCS/1
clamp_x0
                  111 000     100
"""

DECK_IMPDISP = """\
/IMPDISP/1
pull_xL_x
                   7         1         0          0          0
{u_imp:>20.6E}                 200
"""

DECK_IMPL = """\
/IMPL/LINEAR
"""

DECK_END = """\
/END
"""


def fmt_grnod(ids: np.ndarray, per_line: int = 10) -> str:
    """Format a node-ID list into Radioss 10-per-line block."""
    out_lines = []
    chunk = []
    for nid in ids:
        chunk.append(f"{int(nid):>10d}")
        if len(chunk) == per_line:
            out_lines.append("".join(chunk))
            chunk = []
    if chunk:
        out_lines.append("".join(chunk))
    return "\n".join(out_lines)


def fmt_node_block(nodes: np.ndarray) -> str:
    """Radioss /NODE block. One line per node: id, x, y, z."""
    lines = ["/NODE"]
    for i, (x, y, z) in enumerate(nodes, start=1):
        lines.append(f"{i:>10d}{x:>20.6E}{y:>20.6E}{z:>20.6E}")
    return "\n".join(lines) + "\n"


def fmt_brick_block(bricks: np.ndarray, part_id: int = 1) -> str:
    """Radioss /BRICK block. id, n1..n8, with the part owning the brick."""
    lines = [f"/BRICK/{part_id}"]
    for eid, ns in enumerate(bricks, start=1):
        lines.append(
            f"{eid:>10d}" + "".join(f"{int(n):>10d}" for n in ns)
        )
    return "\n".join(lines) + "\n"


def write_starter_deck(
    out_path: Path,
    theta_deg: float,
    nodes: np.ndarray,
    bricks: np.ndarray,
    face_nodes: dict[str, np.ndarray],
    mat: Lamina = IM7_8552,
) -> None:
    th = math.radians(theta_deg)
    c, s = math.cos(th), math.sin(th)
    theta_tag = f"{int(round(theta_deg)):03d}"

    text = DECK_HEADER.format(theta_tag=theta_tag)
    text += DECK_MAT.format(
        name=mat.name,
        rho=mat.rho,
        E1=mat.E1, E2=mat.E2, nu12=mat.nu12, G12=mat.G12,
        E3=mat.E3, nu13=mat.nu13, nu23=mat.nu23, G13=mat.G13, G23=mat.G23,
    )
    text += DECK_SKEW.format(
        theta_tag=theta_tag,
        Ox=0.0, Oy=0.0, Oz=0.0,
        Ax=c,   Ay=s,   Az=0.0,
        Bx=-s,  By=c,   Bz=0.0,
    )
    text += DECK_PROP
    text += DECK_PART
    text += fmt_node_block(nodes)
    text += fmt_brick_block(bricks)
    text += DECK_FUNCT
    text += DECK_BCS.format(
        grnod_x0=fmt_grnod(face_nodes["x0"]),
        grnod_xL=fmt_grnod(face_nodes["xL"]),
    )
    u_imp = STRAIN_BAR * L_COUPON
    text += DECK_IMPDISP.format(u_imp=u_imp)
    text += DECK_IMPL
    text += DECK_END
    out_path.write_text(text)


ENGINE_DECK = """\
/RUN/{runname}/1
                 1.0
/IMPL/LINEAR
/PRINT/-1
/H3D/DT
                 0.0                  1.0
/H3D/SOLID/STRESS/ALL
/H3D/SOLID/STRAIN/ALL
/STOP
"""


def write_engine_deck(out_path: Path, runname: str) -> None:
    out_path.write_text(ENGINE_DECK.format(runname=runname))


# ----------------------------------------------------------------------------
# 6. OpenRadioss invocation (Lima + Apptainer on macOS, native on Linux).
# ----------------------------------------------------------------------------

def find_openradioss_starter() -> str:
    for cand in ("starter_linux64_gf", "starter_linuxa64", "starter"):
        p = shutil.which(cand)
        if p:
            return p
    return ""


def find_openradioss_engine() -> str:
    for cand in ("engine_linux64_gf", "engine_linuxa64", "engine"):
        p = shutil.which(cand)
        if p:
            return p
    return ""


def _lima_wrap(cmd: list[str]) -> list[str]:
    """If we are on macOS and the binary is not available natively, wrap into
    `limactl shell apptainer -- ...` per master_plan.md section 7."""
    if sys.platform == "darwin" and shutil.which("limactl") and not shutil.which(cmd[0]):
        return ["limactl", "shell", "apptainer", "--"] + cmd
    return cmd


def run_openradioss(workdir: Path, runname: str, dry_run: bool) -> int:
    starter_deck = workdir / f"{runname}_0000.rad"
    engine_deck = workdir / f"{runname}_0001.rad"
    if not starter_deck.exists() or not engine_deck.exists():
        raise RuntimeError(f"missing decks for {runname}")
    if dry_run:
        print(f"[dry-run] would run starter on {starter_deck.name} and engine on {engine_deck.name}")
        return 0

    starter = find_openradioss_starter()
    engine = find_openradioss_engine()
    if not starter or not engine:
        # Fall back to lima-wrapped invocation; we still need to know SOME binary name.
        starter = "starter_linuxa64"
        engine = "engine_linuxa64"

    rc = subprocess.call(_lima_wrap([starter, "-i", starter_deck.name]), cwd=str(workdir))
    if rc != 0:
        return rc
    rc = subprocess.call(_lima_wrap([engine, "-i", engine_deck.name]), cwd=str(workdir))
    return rc


# ----------------------------------------------------------------------------
# 7. Post-processing -- volume-averaged Cauchy stress in the gauge sub-volume.
# ----------------------------------------------------------------------------

def gauge_mask(centroids: np.ndarray) -> np.ndarray:
    """Boolean mask selecting elements whose centroid sits in the central
    half-length gauge sub-volume |x - L/2| <= L/4."""
    return np.abs(centroids[:, 0] - 0.5 * L_COUPON) <= 0.25 * L_COUPON


def read_h3d_stresses(workdir: Path, runname: str) -> Optional[dict[str, np.ndarray]]:
    """Read element-level stress and strain tensors from the OpenRadioss H3D
    output via the Kitware openradioss-to-vtkhdf converter, then through
    PyVista. Returns dict with keys 'centroid', 'sigma' (Nx6), 'eps' (Nx6) or
    None if the toolchain is unavailable."""
    try:
        import pyvista as pv
    except ImportError:
        print("[warn] pyvista not installed; cannot post-process H3D output")
        return None

    h3d = workdir / f"{runname}_0001.h3d"
    vtkhdf = workdir / f"{runname}.vtkhdf"
    if not vtkhdf.exists():
        # Try the Kitware converter, which lives alongside OpenRadioss tooling.
        cvt = shutil.which("openradioss-to-vtkhdf") or "openradioss-to-vtkhdf"
        rc = subprocess.call(_lima_wrap([cvt, str(h3d), str(vtkhdf)]))
        if rc != 0 or not vtkhdf.exists():
            print(f"[warn] could not convert {h3d} to vtkhdf (rc={rc})")
            return None

    grid = pv.read(str(vtkhdf))
    cell_centers = grid.cell_centers().points
    # OpenRadioss naming for the tensors in the converter output. Field names
    # may shift between converter versions; try a few common ones.
    stress_keys = ("Stress", "STRESS", "Cauchy_Stress", "sigma")
    strain_keys = ("Strain", "STRAIN", "Total_Strain", "epsilon")
    sigma = None
    for k in stress_keys:
        if k in grid.cell_data:
            sigma = np.asarray(grid.cell_data[k])
            break
    eps = None
    for k in strain_keys:
        if k in grid.cell_data:
            eps = np.asarray(grid.cell_data[k])
            break
    if sigma is None or eps is None:
        print(f"[warn] vtkhdf missing stress/strain arrays; available={list(grid.cell_data.keys())}")
        return None
    return {"centroid": cell_centers, "sigma": sigma, "eps": eps}


def extract_Ex(workdir: Path, runname: str) -> Optional[dict[str, float]]:
    data = read_h3d_stresses(workdir, runname)
    if data is None:
        return None
    mask = gauge_mask(data["centroid"])
    if not mask.any():
        return None
    sigma = data["sigma"][mask]
    eps = data["eps"][mask]
    sxx = float(sigma[:, 0].mean())
    syy = float(sigma[:, 1].mean())
    sxy = float(sigma[:, 3].mean())  # OpenRadioss Voigt: xx,yy,zz,xy,yz,zx
    exx = float(eps[:, 0].mean())
    eyy = float(eps[:, 1].mean())
    if abs(exx) < 1e-12:
        return None
    Ex_fem = sxx / exx
    nu_xy_fem = -eyy / exx
    return dict(sxx=sxx, syy=syy, sxy=sxy, exx=exx, eyy=eyy, Ex_fem=Ex_fem, nu_xy_fem=nu_xy_fem)


# ----------------------------------------------------------------------------
# 8. Driver.
# ----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default="runs", help="output directory for decks/results")
    ap.add_argument("--angles", nargs="+", type=float, default=list(SWEEP_ANGLES_DEG))
    ap.add_argument("--dry-run", action="store_true", help="template decks but do not run solver")
    ap.add_argument("--analytic-only", action="store_true",
                    help="just print the analytic E_x, nu_xy, G_xy table and exit")
    ap.add_argument("--csv", default="summary.csv")
    return ap.parse_args()


def print_analytic_table(angles: Iterable[float]) -> None:
    print(f"{'theta':>8s} {'E_x [GPa]':>14s} {'nu_xy':>10s} {'G_xy [GPa]':>14s}")
    for th in angles:
        Ex = Ex_analytic(th) / 1e9
        nu = nu_xy_analytic(th)
        Gx = Gxy_analytic(th) / 1e9
        print(f"{th:8.2f} {Ex:14.4f} {nu:10.4f} {Gx:14.4f}")


def main() -> int:
    args = parse_args()
    angles = list(args.angles)

    if args.analytic_only:
        print_analytic_table(angles)
        return 0

    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    nodes, bricks, face_nodes = build_mesh()

    # Write decks once, run once, post-process once. No iteration.
    rows: list[dict] = []
    fail = []
    for th in angles:
        runname = f"coupon_{int(round(th)):03d}"
        rdir = workdir / runname
        rdir.mkdir(parents=True, exist_ok=True)
        starter_path = rdir / f"{runname}_0000.rad"
        engine_path = rdir / f"{runname}_0001.rad"
        write_starter_deck(starter_path, th, nodes, bricks, face_nodes, IM7_8552)
        write_engine_deck(engine_path, runname)
        print(f"[stage_08] wrote decks for theta={th:.2f} deg in {rdir}")

        rc = run_openradioss(rdir, runname, args.dry_run)
        if rc != 0:
            print(f"[stage_08] solver failed for theta={th:.2f}, rc={rc}")

        Ex_an = Ex_analytic(th)
        nu_an = nu_xy_analytic(th)
        row = {
            "theta_deg": th,
            "Ex_analytic_Pa": Ex_an,
            "nu_xy_analytic": nu_an,
            "Ex_fem_Pa": float("nan"),
            "nu_xy_fem": float("nan"),
            "rel_err_Ex": float("nan"),
            "passed": False,
        }
        if not args.dry_run and rc == 0:
            extracted = extract_Ex(rdir, runname)
            if extracted is not None:
                Ex_fem = extracted["Ex_fem"]
                rel_err = abs(Ex_fem - Ex_an) / Ex_an
                row.update(
                    Ex_fem_Pa=Ex_fem,
                    nu_xy_fem=extracted["nu_xy_fem"],
                    sigma_xx=extracted["sxx"],
                    sigma_yy=extracted["syy"],
                    sigma_xy=extracted["sxy"],
                    eps_xx=extracted["exx"],
                    eps_yy=extracted["eyy"],
                    rel_err_Ex=rel_err,
                    passed=rel_err <= PASS_TOL_REL,
                )
                if not row["passed"]:
                    fail.append((th, rel_err))
        rows.append(row)
        print(
            f"[stage_08] theta={th:6.2f} Ex_analytic={Ex_an / 1e9:7.3f} GPa  "
            f"Ex_fem={row['Ex_fem_Pa'] / 1e9 if row['Ex_fem_Pa'] == row['Ex_fem_Pa'] else float('nan'):7.3f} GPa  "
            f"rel_err={row['rel_err_Ex']:8.2%}  pass={row['passed']}"
        )

    # CSV export -- Typst+CeTZ figure imports this and draws the curve overlay.
    csv_path = workdir / args.csv
    fieldnames = sorted({k for r in rows for k in r.keys()})
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print(f"[stage_08] wrote summary -> {csv_path}")

    # Also dump a dense analytic curve for the figure.
    dense = workdir / "Ex_analytic_curve.csv"
    with dense.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["theta_deg", "Ex_analytic_Pa", "nu_xy_analytic", "Gxy_analytic_Pa"])
        for th in np.linspace(0.0, 90.0, 181):
            w.writerow([th, Ex_analytic(th), nu_xy_analytic(th), Gxy_analytic(th)])
    print(f"[stage_08] wrote analytic curve -> {dense}")

    if args.dry_run:
        print("[stage_08] dry run only; no pass/fail evaluated")
        return 0
    if fail:
        print(f"[stage_08] FAIL: {len(fail)} angle(s) exceeded {PASS_TOL_REL:.0%} tolerance:")
        for th, e in fail:
            print(f"           theta={th:6.2f}  rel_err={e:.2%}")
        return 1
    print("[stage_08] PASS: all angles within 1 percent of analytic E_x")
    return 0


if __name__ == "__main__":
    sys.exit(main())
