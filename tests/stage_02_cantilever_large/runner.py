"""Stage 02 runner: cantilever large-deflection geometric nonlinearity.

Author: J.C. Vaught
Project: FEA_AP-PLY
Stage:   02 of 16

Skeleton driver. Templates a per-alpha OpenRadioss starter+engine deck,
invokes ``starter`` and ``engine`` inside Lima/Apptainer, parses the T01
time-history for tip-centroid displacement, and compares against the
Bisshopp-Drucker (1945) elastica reference solution.

Pass criterion (per spec.md section 8): tip-displacement components
within 2 percent of the elliptic-integral reference at alpha in {1, 3, 5}.

This is a SKELETON. Blocks marked TODO are filled in once Stage 01 has
established the OpenRadioss toolchain end-to-end on the macOS host.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from scipy import optimize, special  # type: ignore[import-not-found]

# Constants -- spec.md sections 2 and 5 (SI: kg, m, s, Pa, N).
L_BEAM, B_BEAM, H_BEAM = 1.000, 0.0250, 0.0030
E_MOD, NU, RHO = 2.0e11, 0.30, 7850.0
I_ZZ = B_BEAM * H_BEAM ** 3 / 12.0
EI = E_MOD * I_ZZ                       # approx 11.25 N*m^2
PASS_ALPHAS = (1.0, 3.0, 5.0)             # spec.md sec 8 pass/fail gate
PLOT_ALPHAS = (0.05, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0)  # 0.05 gates the linear limit
PASS_TOL = 0.02


@dataclass(frozen=True)
class ElasticaRef:
    alpha: float
    theta_L: float
    delta_y_over_L: float
    delta_x_over_L: float


def _alpha_of_theta(theta_L: float) -> float:
    """Map tip rotation theta_L (rad) to alpha = P L^2 / EI via elliptic integrals."""
    k = math.sqrt((1.0 + math.sin(theta_L)) / 2.0)
    sin_phi0 = max(-1.0, min(1.0, 1.0 / (k * math.sqrt(2.0))))
    phi0 = math.asin(sin_phi0)
    big_K = float(special.ellipk(k * k))
    F_phi0 = float(special.ellipkinc(phi0, k * k))
    # Belendez et al. 2002: sqrt(alpha) = K(k) - F(phi_0, k); see spec.md sec 7.2.
    return (big_K - F_phi0) ** 2


def elastica_reference(alpha: float) -> ElasticaRef:
    """Closed-form Bisshopp-Drucker elastica solution for tip displacement.

    Returns delta_y / L (vertical sag, positive in load direction) and
    delta_x / L (horizontal foreshortening, positive toward the clamp).
    Uses x_tip/L = sqrt(2 sin theta_L)/sqrt(alpha) and
    delta_y/L = 1 - 2(E(k) - E(phi_0, k))/sqrt(alpha). See spec.md sec 7.3.
    """
    if alpha <= 0.0:
        return ElasticaRef(alpha, 0.0, 0.0, 0.0)
    theta_L = optimize.brentq(lambda th: _alpha_of_theta(th) - alpha,
                              1.0e-6, math.pi / 2.0 - 1.0e-6, xtol=1.0e-10)
    k = math.sqrt((1.0 + math.sin(theta_L)) / 2.0)
    phi0 = math.asin(max(-1.0, min(1.0, 1.0 / (k * math.sqrt(2.0)))))
    big_K = float(special.ellipk(k * k));   big_E = float(special.ellipe(k * k))
    F_phi0 = float(special.ellipkinc(phi0, k * k))
    E_phi0 = float(special.ellipeinc(phi0, k * k))
    sqrt_alpha = big_K - F_phi0
    delta_x_over_L = 1.0 - math.sqrt(2.0 * math.sin(theta_L)) / sqrt_alpha
    delta_y_over_L = 1.0 - 2.0 * (big_E - E_phi0) / sqrt_alpha
    return ElasticaRef(alpha, theta_L, delta_y_over_L, delta_x_over_L)


def template_deck(alpha: float, work_dir: Path) -> tuple[Path, Path]:
    """Render starter + engine .rad decks for one alpha into ``work_dir``.

    TODO: load Jinja2 templates ./templates/{starter,engine}.rad.j2 and
    render with the substitution dict below. Currently writes placeholders.
    """
    P = alpha * EI / L_BEAM ** 2
    subs = dict(alpha=alpha, P_total=P, L_beam=L_BEAM, b_beam=B_BEAM, h_beam=H_BEAM,
                E_mod=E_MOD, nu=NU, rho=RHO, Nx=80, Ny=6, Nz=3,
                ismstr=11, icpre=1, impl_solver=2,        # see spec.md section 6.1
                dt_init=0.05, dt_min=1.0e-4, dt_max=0.5)
    work_dir.mkdir(parents=True, exist_ok=True)
    starter = work_dir / f"cantilever_alpha_{alpha:0.2f}_0000.rad"
    engine = work_dir / f"cantilever_alpha_{alpha:0.2f}_0001.rad"
    starter.write_text(f"# placeholder starter, alpha={alpha}\n# subs={subs}\n")
    engine.write_text(f"# placeholder engine, alpha={alpha}\n")
    return starter, engine


def run_openradioss(starter: Path, engine: Path, lima_vm: str = "apptainer") -> None:
    """Invoke OpenRadioss starter + engine inside the Lima ARM64 VM."""
    for label, deck in (("starter", starter), ("engine", engine)):
        cmd = ["limactl", "shell", lima_vm, "--", f"{label}_linuxa64",
               "-i", deck.name, "-nt", "1"]
        proc = subprocess.run(cmd, cwd=deck.parent, capture_output=True,
                              text=True, check=False)
        log = deck.parent / f"{label}.log"
        log.write_text(proc.stdout + "\n--- STDERR ---\n" + proc.stderr)
        if proc.returncode != 0:
            raise RuntimeError(f"OpenRadioss {label} failed: see {log}")


def parse_tip_displacement(work_dir: Path, alpha: float) -> tuple[float, float]:
    """Return (delta_x, delta_y) in metres from the T01 time-history.

    TODO: use vortex-radioss (vortex_radioss.read_th) or
    OpenRadioss/Tools/th_to_csv.py to extract the tip-centroid node
    displacement at the final pseudo-time.
    """
    th_file = work_dir / f"cantilever_alpha_{alpha:0.2f}T01"
    if not th_file.exists():
        raise FileNotFoundError(f"T01 not found: {th_file}")
    raise NotImplementedError("T01 parsing not yet wired up; see spec.md section 9 step 5")


@dataclass
class SweepRow:
    alpha: float
    dy_fem: float; dx_fem: float
    dy_ref: float; dx_ref: float
    err_y: float; err_x: float
    passed: bool


def run_sweep(work_root: Path, alphas: tuple[float, ...] = PLOT_ALPHAS) -> list[SweepRow]:
    rows: list[SweepRow] = []
    for alpha in alphas:
        case = work_root / f"alpha_{alpha:0.2f}"
        starter, engine = template_deck(alpha, case)
        run_openradioss(starter, engine)
        dx_m, dy_m = parse_tip_displacement(case, alpha)
        ref = elastica_reference(alpha)
        dy, dx = abs(dy_m) / L_BEAM, abs(dx_m) / L_BEAM
        err_y = abs(dy - ref.delta_y_over_L) / max(ref.delta_y_over_L, 1.0e-12)
        err_x = abs(dx - ref.delta_x_over_L) / max(ref.delta_x_over_L, 1.0e-12)
        gated = alpha in PASS_ALPHAS
        passed = (not gated) or (err_y <= PASS_TOL and err_x <= PASS_TOL)
        rows.append(SweepRow(alpha, dy, dx, ref.delta_y_over_L, ref.delta_x_over_L,
                             err_y, err_x, passed))
    return rows


def write_csv(rows: list[SweepRow], path: Path) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["alpha", "dy_fem_over_L", "dx_fem_over_L",
                    "dy_ref_over_L", "dx_ref_over_L", "err_y", "err_x", "passed"])
        for r in rows:
            w.writerow([r.alpha, r.dy_fem, r.dx_fem, r.dy_ref, r.dx_ref,
                        r.err_y, r.err_x, int(r.passed)])


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage 02 cantilever large-deflection runner")
    parser.add_argument("--work-dir", type=Path, default=Path("./work"))
    parser.add_argument("--ref-only", action="store_true",
                        help="Skip FEM, just print elastica reference table")
    args = parser.parse_args()
    if args.ref_only:
        print(f"{'alpha':>6} {'theta_L':>10} {'dy/L':>10} {'dx/L':>10}")
        for a in PLOT_ALPHAS:
            r = elastica_reference(a)
            print(f"{r.alpha:>6.2f} {r.theta_L:>10.4f} "
                  f"{r.delta_y_over_L:>10.4f} {r.delta_x_over_L:>10.4f}")
        return 0
    if shutil.which("limactl") is None:
        raise SystemExit("limactl not on PATH; OpenRadioss requires Lima per master plan section 7")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    rows = run_sweep(args.work_dir)
    write_csv(rows, args.work_dir / "tip_displacement_vs_alpha.csv")
    print(f"{'alpha':>6} {'dy/L FEM':>10} {'dy/L ref':>10} {'err_y':>8} "
          f"{'dx/L FEM':>10} {'dx/L ref':>10} {'err_x':>8} verdict")
    overall = True
    for r in rows:
        gated = r.alpha in PASS_ALPHAS
        verdict = ("PASS" if r.passed else "FAIL") if gated else "info"
        if gated and not r.passed:
            overall = False
        print(f"{r.alpha:>6.2f} {r.dy_fem:>10.4f} {r.dy_ref:>10.4f} {r.err_y:>8.4f} "
              f"{r.dx_fem:>10.4f} {r.dx_ref:>10.4f} {r.err_x:>8.4f} {verdict}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
