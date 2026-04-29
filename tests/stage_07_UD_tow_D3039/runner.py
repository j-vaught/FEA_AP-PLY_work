"""
Stage 07 runner. ASTM D3039 unidirectional tow tension on IM7/8552 in
OpenRadioss, implicit static, /MAT/LAW25 on /PROP/TYPE14 solid HEXA8.

Author: J.C. Vaught
Date: 2026-04-29

Usage
-----
    python runner.py                    # both 7A (0 deg) and 7B (45 deg)
    python runner.py --only 7A          # on-axis only
    python runner.py --only 7B          # off-axis only
    python runner.py --dry-run          # template decks but don't launch
    python runner.py --refine           # rerun with h/2 width refinement
    python runner.py --no-lima          # call starter/engine on PATH
                                        # directly, skipping Lima shell

What this script does
---------------------
1. Templates the OpenRadioss starter / engine .rad decks for both runs
   (psi = 0 deg and psi = 45 deg) using the IM7/8552 Soden 1998 card.
2. Generates the structured HEXA8 mesh inline (no GMSH dependency for
   this stage).
3. Invokes the OpenRadioss starter and engine binaries inside the
   user's documented Lima Apptainer environment (master plan section
   7) unless --no-lima is passed.
4. Parses the resulting T01 time-history file, extracts apparent
   modulus from a least-squares fit of sigma_xx vs eps_xx in the
   gauge window (spec section 9).
5. Compares to closed-form: E1 for 7A; the off-axis transformation
   formula at theta = 45 deg for 7B.
6. Writes a pass/fail CSV summary and a free-form log.

Hard requirements: solid elements only, SI base. The deck uses the
OpenRadioss "Mg, mm, s" unit set (mass in megagrams, length in
millimeters, time in seconds), which gives stress in MPa, density in
Mg/mm^3 = g/mm^3 * 1e-3 *no wait*: Mg/mm^3 = 1e9 kg/m^3, so density of
1.58 g/cm^3 = 1580 kg/m^3 = 1.58e-9 Mg/mm^3. We use that conversion
explicitly in the MAT_PARAMS dict below.
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
from typing import Iterable

# ----------------------------------------------------------------------
# 1. Material card. IM7/8552 from Soden, Hinton, Kaddour 1998
#    (and Kaddour-Hinton 2013 for nu23). Units in the deck are Mg, mm, s
#    so stress is MPa, density is Mg/mm^3.
# ----------------------------------------------------------------------

MAT_PARAMS: dict[str, float] = {
    # density: 1580 kg/m^3 = 1.58e-9 Mg/mm^3
    "rho": 1.58e-9,
    # stiffness, MPa
    "E1": 171_000.0,
    "E2": 9_080.0,
    "E3": 9_080.0,
    "nu12": 0.32,
    "nu13": 0.32,
    "nu23": 0.50,
    "G12": 5_290.0,
    "G13": 5_290.0,
    # G23 derived: E2 / (2 * (1 + nu23)) = 9080 / 3 = 3026.67
    "G23": 9_080.0 / (2.0 * (1.0 + 0.50)),
    # strengths, MPa  (used by LAW25 even in linear regime; not exercised
    # at the 0.2 percent strain we apply here, but kept so the deck is
    # reusable downstream).
    "Xt": 2_326.0,
    "Xc": 1_700.0,    # Soden 1998 IM7/8552 longitudinal compressive
    "Yt": 64.7,
    "Yc": 200.0,      # Soden 1998 IM7/8552 transverse compressive
    "S":  92.3,
}


def Ex_offaxis(theta_deg: float, p: dict[str, float]) -> float:
    """Closed-form off-axis Young's modulus for plane-stress orthotropic
    UD ply, Jones 1999 / Daniel-Ishai 2006:

        1/E_x = c^4/E1 + s^4/E2 + (1/G12 - 2 nu12/E1) c^2 s^2

    Returns E_x in the same units as the inputs (MPa here).
    """
    th = math.radians(theta_deg)
    c2, s2 = math.cos(th) ** 2, math.sin(th) ** 2
    c4, s4 = c2 * c2, s2 * s2
    E1, E2, G12, n12 = p["E1"], p["E2"], p["G12"], p["nu12"]
    inv = c4 / E1 + s4 / E2 + (1.0 / G12 - 2.0 * n12 / E1) * c2 * s2
    return 1.0 / inv


# ----------------------------------------------------------------------
# 2. Per-run configuration
# ----------------------------------------------------------------------

@dataclass
class RunCfg:
    name: str           # "7A" or "7B"
    theta_deg: float    # 0 or 45
    L: float            # coupon length, mm
    w: float            # coupon width, mm
    t: float            # coupon thickness, mm
    nx: int             # elements along length
    ny: int             # elements across width
    nz: int             # elements through thickness
    Lg: float           # gauge window length, mm (centered at L/2)
    eps_max: float      # peak applied longitudinal strain
    n_steps: int        # implicit load steps


def default_runs() -> list[RunCfg]:
    return [
        RunCfg(name="7A", theta_deg=0.0,
               L=250.0, w=15.0, t=1.0,
               nx=80, ny=8,  nz=1,
               Lg=100.0, eps_max=2.0e-3, n_steps=5),
        RunCfg(name="7B", theta_deg=45.0,
               L=250.0, w=25.0, t=1.0,
               nx=80, ny=12, nz=1,
               Lg=100.0, eps_max=2.0e-3, n_steps=5),
    ]


# ----------------------------------------------------------------------
# 3. Mesh generation. Structured HEXA8, lexicographic node numbering.
#    Node id = 1 + i + (nx+1)*j + (nx+1)*(ny+1)*k, 1-based.
#    Brick connectivity follows the OpenRadioss /BRICK convention:
#    n1..n4 form the bottom face (z = z_k) ccw seen from +z,
#    n5..n8 form the top face (z = z_{k+1}) ccw seen from +z.
# ----------------------------------------------------------------------

def gen_mesh(cfg: RunCfg) -> tuple[list[tuple[int, float, float, float]],
                                   list[tuple[int, list[int]]],
                                   list[int],   # gripA face nodes (x=0)
                                   list[int],   # gripB face nodes (x=L)
                                   list[int],   # gauge window node line at x=x1
                                   list[int]]:  # gauge window node line at x=x2
    nx, ny, nz = cfg.nx, cfg.ny, cfg.nz
    L, w, t = cfg.L, cfg.w, cfg.t
    Lg = cfg.Lg

    # node coords
    nodes: list[tuple[int, float, float, float]] = []
    nid_of = {}
    nid = 0
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nid += 1
                x = L * i / nx
                y = w * j / ny
                z = t * k / nz
                nodes.append((nid, x, y, z))
                nid_of[(i, j, k)] = nid

    # bricks
    bricks: list[tuple[int, list[int]]] = []
    eid = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                eid += 1
                n1 = nid_of[(i,     j,     k)]
                n2 = nid_of[(i + 1, j,     k)]
                n3 = nid_of[(i + 1, j + 1, k)]
                n4 = nid_of[(i,     j + 1, k)]
                n5 = nid_of[(i,     j,     k + 1)]
                n6 = nid_of[(i + 1, j,     k + 1)]
                n7 = nid_of[(i + 1, j + 1, k + 1)]
                n8 = nid_of[(i,     j + 1, k + 1)]
                bricks.append((eid, [n1, n2, n3, n4, n5, n6, n7, n8]))

    # node sets
    gripA = [nid_of[(0,  j, k)] for k in range(nz + 1) for j in range(ny + 1)]
    gripB = [nid_of[(nx, j, k)] for k in range(nz + 1) for j in range(ny + 1)]

    # gauge window: pick i indices closest to x1 = (L-Lg)/2 and x2 = (L+Lg)/2
    x1, x2 = (L - Lg) / 2.0, (L + Lg) / 2.0
    i1 = round(x1 / L * nx)
    i2 = round(x2 / L * nx)
    gauge_x1 = [nid_of[(i1, j, k)] for k in range(nz + 1) for j in range(ny + 1)]
    gauge_x2 = [nid_of[(i2, j, k)] for k in range(nz + 1) for j in range(ny + 1)]

    return nodes, bricks, gripA, gripB, gauge_x1, gauge_x2


# ----------------------------------------------------------------------
# 4. Deck templating. We deliberately keep this as direct string
#    formatting rather than Jinja2 to avoid an extra dependency for
#    this stage. The OpenRadioss starter free-format reader tolerates
#    mild whitespace variation; column-strict cards (LAW25, PROP) are
#    written with explicit field widths.
# ----------------------------------------------------------------------

def fmt_node_block(nodes) -> str:
    return "\n".join(f"{nid:10d}{x:20.7f}{y:20.7f}{z:20.7f}"
                     for (nid, x, y, z) in nodes)


def fmt_brick_block(part_id: int, bricks) -> str:
    lines = []
    for eid, nlist in bricks:
        # /BRICK card: eid, part, n1..n8
        lines.append(f"{eid:10d}{part_id:10d}" +
                     "".join(f"{n:10d}" for n in nlist))
    return "\n".join(lines)


def fmt_grnod_block(node_ids: Iterable[int]) -> str:
    out = []
    row = []
    for nid in node_ids:
        row.append(f"{nid:10d}")
        if len(row) == 10:
            out.append("".join(row))
            row = []
    if row:
        out.append("".join(row))
    return "\n".join(out)


def write_starter(path: Path, cfg: RunCfg, mat: dict[str, float]) -> None:
    nodes, bricks, gripA, gripB, gauge_x1, gauge_x2 = gen_mesh(cfg)
    u_max = cfg.eps_max * cfg.L

    nodes_block = fmt_node_block(nodes)
    bricks_block = fmt_brick_block(part_id=1, bricks=bricks)
    gripA_block = fmt_grnod_block(gripA)
    gripB_block = fmt_grnod_block(gripB)
    gauge_x1_block = fmt_grnod_block(gauge_x1)
    gauge_x2_block = fmt_grnod_block(gauge_x2)

    deck = f"""#RADIOSS STARTER
/BEGIN
stage7_{cfg.name}
      2024         0
                  Mg                  mm                   s
                  Mg                  mm                   s
/UNIT/1
mass length time
Mg mm s
#---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|
/MAT/LAW25/1
IM7_8552_Soden1998
{mat['rho']:.6e} {mat['E1']:.4f} {mat['E2']:.4f} {mat['E3']:.4f} {mat['nu12']:.4f} {mat['nu13']:.4f} {mat['nu23']:.4f}
{mat['G12']:.4f} {mat['G23']:.4f} {mat['G13']:.4f}      1      0
{mat['Xt']:.4f} {mat['Xc']:.4f} {mat['Yt']:.4f} {mat['Yc']:.4f} {mat['S']:.4f}
0.0 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0
#
/PROP/TYPE14/1
solid_UD_psi_{cfg.theta_deg:.1f}
   1   1   0   0   0   0
0.0 0.0 0.0 1.0 0.0 0.0   0   1   0   0   0   0
{cfg.theta_deg:.4f}
#
/SKEW/FIX/1
global
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
#
/NODE
{nodes_block}
#
/BRICK/1
{bricks_block}
#
/PART/1
ud_coupon
1 1 0
#
/GRNOD/NODE/1
gripA_face
{gripA_block}
/GRNOD/NODE/2
gripB_face
{gripB_block}
/GRNOD/NODE/3
gauge_x1
{gauge_x1_block}
/GRNOD/NODE/4
gauge_x2
{gauge_x2_block}
#
/BCS/1
gripA_clamped
111 111   1
1
#
/BCS/2
gripB_lateral_only
011 111   1
2
#
/IMPDISP/1
gripB_pull_x
2 1 1
   1   1   1   0.0   1.0   {u_max:.6f}
#
/FUNCT/1
ramp
0.0 0.0
1.0 1.0
#
/IMPL/PRINT/N
1
/IMPL/SOLVER/1
0  0  3  0  0
/IMPL/NONLIN/SMDISP
{cfg.n_steps}  0.0  1.0  1.0e-3  1.0e-6  20  0
/IMPL/DT/STOP
1.0e-12   1.0
/IMPL/DTINI
{1.0 / cfg.n_steps:.6f}
#
/TH/NODE/1
gripB_reaction
2 0 0
DX FX
/TH/NODE/2
gauge_x1_th
3 0 0
DX
/TH/NODE/3
gauge_x2_th
4 0 0
DX
#
/ANIM/BRICK/TENS/STRESS
/ANIM/BRICK/TENS/STRAIN
/ANIM/DT
0.2
#
/END
"""
    path.write_text(deck)


def write_engine(path: Path, cfg: RunCfg) -> None:
    deck = f"""#RADIOSS ENGINE
/RUN/stage7_{cfg.name}/1
1.0
/PRINT/-100
/HIS/DT
0.05
/STOP
0.0  1.0e-3  0.0  0.0
/END
"""
    path.write_text(deck)


# ----------------------------------------------------------------------
# 5. Solver invocation. Defaults to Lima Apptainer wrapper per master
#    plan section 7. The user can override the executable paths and
#    the wrapper command via env vars.
# ----------------------------------------------------------------------

def _resolve_solver_command(use_lima: bool) -> tuple[list[str], list[str]]:
    """Return (starter_cmd_prefix, engine_cmd_prefix). Each is a list of
    argv tokens that will have ['-i', deck_path] appended.
    """
    starter_bin = os.environ.get("OR_STARTER", "starter_linuxa64")
    engine_bin = os.environ.get("OR_ENGINE", "engine_linuxa64")
    if use_lima:
        # documented Lima Apptainer path: limactl shell apptainer -- <bin>
        lima_inst = os.environ.get("OR_LIMA_INSTANCE", "apptainer")
        prefix = ["limactl", "shell", lima_inst, "--"]
        return prefix + [starter_bin], prefix + [engine_bin]
    return [starter_bin], [engine_bin]


def run_or(cfg: RunCfg, work_dir: Path, use_lima: bool, dry_run: bool) -> None:
    starter_deck = work_dir / f"stage7_{cfg.name}_0000.rad"
    engine_deck  = work_dir / f"stage7_{cfg.name}_0001.rad"
    if dry_run:
        print(f"[dry-run] decks templated to {work_dir}; skipping solver")
        return
    starter_cmd, engine_cmd = _resolve_solver_command(use_lima)
    starter_argv = starter_cmd + ["-i", str(starter_deck)]
    engine_argv  = engine_cmd  + ["-i", str(engine_deck)]
    print(f"[stage7][{cfg.name}] starter: {' '.join(starter_argv)}")
    subprocess.run(starter_argv, check=True, cwd=work_dir)
    print(f"[stage7][{cfg.name}] engine : {' '.join(engine_argv)}")
    subprocess.run(engine_argv, check=True, cwd=work_dir)


# ----------------------------------------------------------------------
# 6. T01 time-history parsing. We support two backends: the
#    vortex-radioss Python package (preferred), and a fallback ASCII
#    T01 export that the OpenRadioss starter can produce when /TH is
#    requested with text output. The fallback is implemented inline.
# ----------------------------------------------------------------------

def parse_th(work_dir: Path, cfg: RunCfg) -> dict:
    """Parse the time-history file produced by OpenRadioss.

    Returns a dict with keys:
        t           list[float]  pseudo-time samples
        F_Bx        list[float]  summed reaction at grip-B in x (N)
        ux_x1       list[float]  through-width-mean ux at x = x1 (mm)
        ux_x2       list[float]  through-width-mean ux at x = x2 (mm)
    """
    # Backend 1: vortex-radioss
    try:
        from vortex_radioss import T01  # type: ignore
        t01_path = next(work_dir.glob(f"stage7_{cfg.name}*T01"))
        th = T01(str(t01_path))
        t = list(th.time)
        F_Bx = list(th.node_group("gripB_reaction").FX_sum)
        ux_x1 = list(th.node_group("gauge_x1_th").DX_mean)
        ux_x2 = list(th.node_group("gauge_x2_th").DX_mean)
        return dict(t=t, F_Bx=F_Bx, ux_x1=ux_x1, ux_x2=ux_x2)
    except Exception:
        pass

    # Backend 2: ASCII T01 fallback. The OpenRadioss "T01 ASCII" mode
    # writes one column per requested channel. We expect channel
    # ordering: time, gripB_reaction (DX, FX summed across the group),
    # gauge_x1 (DX averaged), gauge_x2 (DX averaged).
    ascii_path = work_dir / f"stage7_{cfg.name}.txt"
    if not ascii_path.exists():
        raise RuntimeError(
            f"could not locate T01 output for run {cfg.name}; "
            f"looked for vortex-radioss T01 and ASCII fallback "
            f"{ascii_path}. Re-run with /HIS/ASCII or install vortex-radioss."
        )
    t, F_Bx, ux_x1, ux_x2 = [], [], [], []
    with ascii_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 5:
                continue
            # cols: time, ux_gripB_mean, FX_gripB_sum, ux_x1_mean, ux_x2_mean
            t.append(float(cols[0]))
            F_Bx.append(float(cols[2]))
            ux_x1.append(float(cols[3]))
            ux_x2.append(float(cols[4]))
    return dict(t=t, F_Bx=F_Bx, ux_x1=ux_x1, ux_x2=ux_x2)


# ----------------------------------------------------------------------
# 7. Modulus extraction. Linear regression of sigma_xx on eps_xx in the
#    gauge window. Returns slope (= apparent E in MPa), R^2.
# ----------------------------------------------------------------------

def _linregress(x: list[float], y: list[float]) -> tuple[float, float, float]:
    n = len(x)
    if n < 2:
        return float("nan"), float("nan"), float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x)
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    syy = sum((yi - my) ** 2 for yi in y)
    if sxx == 0:
        return float("nan"), float("nan"), float("nan")
    slope = sxy / sxx
    intercept = my - slope * mx
    r2 = (sxy * sxy) / (sxx * syy) if syy > 0 else 1.0
    return slope, intercept, r2


def extract_modulus(th: dict, cfg: RunCfg) -> dict:
    A = cfg.w * cfg.t                 # mm^2
    F = th["F_Bx"]                    # N (Mg.mm/s^2 in the deck units; N in SI)
    ux1 = th["ux_x1"]                 # mm
    ux2 = th["ux_x2"]                 # mm
    sigma = [Fi / A for Fi in F]                  # MPa
    eps   = [(u2 - u1) / cfg.Lg for u1, u2 in zip(ux1, ux2)]
    slope, intercept, r2 = _linregress(eps, sigma)
    return dict(sigma=sigma, eps=eps, E_FEM=slope, intercept=intercept, R2=r2)


# ----------------------------------------------------------------------
# 8. Driver
# ----------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", choices=["7A", "7B"], default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--refine", action="store_true",
                        help="rerun with h/2 width refinement")
    parser.add_argument("--no-lima", action="store_true",
                        help="invoke starter/engine on PATH directly "
                             "(no Lima Apptainer wrapper)")
    parser.add_argument("--workdir", default=str(here / "runs"))
    parser.add_argument("--results", default=str(here / "results"))
    args = parser.parse_args(argv)

    work_root = Path(args.workdir).resolve()
    res_root  = Path(args.results).resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    res_root.mkdir(parents=True, exist_ok=True)

    runs = default_runs()
    if args.only:
        runs = [r for r in runs if r.name == args.only]
    if args.refine:
        for r in runs:
            r.ny *= 2
            r.nx = max(r.nx, 80)

    summary_rows: list[dict] = []
    log_lines: list[str] = []

    for cfg in runs:
        rdir = work_root / cfg.name
        rdir.mkdir(parents=True, exist_ok=True)
        write_starter(rdir / f"stage7_{cfg.name}_0000.rad", cfg, MAT_PARAMS)
        write_engine (rdir / f"stage7_{cfg.name}_0001.rad", cfg)
        run_or(cfg, rdir, use_lima=not args.no_lima, dry_run=args.dry_run)

        if args.dry_run:
            log_lines.append(f"[{cfg.name}] dry-run, no extraction.")
            continue

        th = parse_th(rdir, cfg)
        ext = extract_modulus(th, cfg)
        E_FEM = ext["E_FEM"]
        R2 = ext["R2"]

        E_ref = (MAT_PARAMS["E1"] if cfg.theta_deg == 0.0
                 else Ex_offaxis(cfg.theta_deg, MAT_PARAMS))
        rel_err = abs(E_FEM - E_ref) / E_ref
        passed = (rel_err <= 0.02) and (R2 >= 0.9999)

        msg = (
            f"[{cfg.name}] theta={cfg.theta_deg:>5.1f} deg  "
            f"E_FEM={E_FEM:10.2f} MPa  E_ref={E_ref:10.2f} MPa  "
            f"rel_err={rel_err*100:6.3f} %  R^2={R2:.6f}  "
            f"{'PASS' if passed else 'FAIL'}"
        )
        print(msg)
        log_lines.append(msg)
        summary_rows.append(dict(
            run=cfg.name,
            theta_deg=cfg.theta_deg,
            E_FEM_MPa=E_FEM,
            E_analytic_MPa=E_ref,
            rel_err=rel_err,
            R2=R2,
            passed=passed,
        ))

    # write CSV summary
    csv_path = res_root / "stage7_summary.csv"
    if summary_rows:
        with csv_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
            writer.writeheader()
            for row in summary_rows:
                writer.writerow(row)
        print(f"[stage7] wrote {csv_path}")

    # write log
    log_path = res_root / "stage7_log.txt"
    log_path.write_text("\n".join(log_lines) + "\n")

    # exit code: 0 only if every run passed
    if not summary_rows:
        return 0 if args.dry_run else 1
    return 0 if all(r["passed"] for r in summary_rows) else 1


if __name__ == "__main__":
    sys.exit(main())
