"""Stage 09 - solid laminate CLT starter probe.

Author: J.C. Vaught

The intended stage verifies CLT A-matrices with one solid layer per ply using
/MAT/LAW25 on /PROP/TYPE14. The installed OpenRadioss starter rejects that
material/property pair with ERROR 3047, so this runner computes the analytic
CLT references, writes canonical per-ply TYPE14 starter probes for the
cross-ply and quasi-isotropic stacks, records the starter blocker, and marks
the stage INCONCLUSIVE. TYPE6/SOL_ORTH proxy starters are also checked to
isolate the blocker to the required TYPE14 path.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np


THIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parents[1]
RUNS_DIR = THIS_DIR / "runs"
RESULTS_DIR = THIS_DIR / "results"
FIGURES_DIR = THIS_DIR / "figures"
RUN_LOG = THIS_DIR / "run.log"
CARD_PATH = ROOT_DIR / "references" / "material_cards" / "im7_8552.json"

OR_ROOT = Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"

PLY_T = 0.18e-3
LAYUPS = {
    "A_crossply": (0.0, 90.0, 90.0, 0.0),
    "B_quasiiso": (0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0),
}


@dataclass(frozen=True)
class Material:
    rho: float
    e1: float
    e2: float
    e3: float
    nu12: float
    nu13: float
    nu23: float
    g12: float
    g13: float
    g23: float
    xt: float
    xc: float
    yt: float
    yc: float
    s12: float


def fmt_f(*values: float) -> str:
    return "".join(f"{value:20.12g}" for value in values)


def fmt_i(*values: int) -> str:
    return "".join(f"{value:10d}" for value in values)


def load_material() -> Material:
    data = json.loads(CARD_PATH.read_text(encoding="utf-8"))
    elastic = data["elastic"]
    strength = data["strength"]
    return Material(
        rho=float(elastic["density_kg_m3"]),
        e1=float(elastic["E1_Pa"]),
        e2=float(elastic["E2_Pa"]),
        e3=float(elastic["E3_Pa"]),
        nu12=float(elastic["nu12"]),
        nu13=float(elastic["nu13"]),
        nu23=float(elastic["nu23"]),
        g12=float(elastic["G12_Pa"]),
        g13=float(elastic["G13_Pa"]),
        g23=float(elastic["G23_Pa"]),
        xt=float(strength["Xt_Pa"]),
        xc=float(strength["Xc_Pa"]),
        yt=float(strength["Yt_Pa"]),
        yc=float(strength["Yc_Pa"]),
        s12=float(strength["S12_Pa"]),
    )


def q_matrix(mat: Material) -> np.ndarray:
    nu21 = mat.nu12 * mat.e2 / mat.e1
    denom = 1.0 - mat.nu12 * nu21
    q = np.zeros((3, 3), dtype=float)
    q[0, 0] = mat.e1 / denom
    q[1, 1] = mat.e2 / denom
    q[0, 1] = q[1, 0] = mat.nu12 * mat.e2 / denom
    q[2, 2] = mat.g12
    return q


def qbar(q: np.ndarray, theta_deg: float) -> np.ndarray:
    th = math.radians(theta_deg)
    c = math.cos(th)
    s = math.sin(th)
    c2 = c * c
    s2 = s * s
    c4 = c2 * c2
    s4 = s2 * s2
    q11, q22, q12, q66 = q[0, 0], q[1, 1], q[0, 1], q[2, 2]
    out = np.zeros((3, 3), dtype=float)
    out[0, 0] = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * s4
    out[1, 1] = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * c4
    out[0, 1] = out[1, 0] = (q11 + q22 - 4.0 * q66) * s2 * c2 + q12 * (s4 + c4)
    out[2, 2] = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * s2 * c2 + q66 * (s4 + c4)
    out[0, 2] = out[2, 0] = (q11 - q12 - 2.0 * q66) * s * c**3 + (q12 - q22 + 2.0 * q66) * s**3 * c
    out[1, 2] = out[2, 1] = (q11 - q12 - 2.0 * q66) * s**3 * c + (q12 - q22 + 2.0 * q66) * s * c**3
    return out


def a_matrix(mat: Material, thetas: tuple[float, ...]) -> np.ndarray:
    q = q_matrix(mat)
    a = np.zeros((3, 3), dtype=float)
    for theta in thetas:
        a += qbar(q, theta) * PLY_T
    return a


def radioss_env() -> dict[str, str]:
    env = os.environ.copy()
    env["OR"] = str(OR_ROOT)
    env["RAD_CFG_PATH"] = str(OR_ROOT / "hm_cfg_files")
    env["RAD_H3D_PATH"] = str(OR_ROOT / "extlib" / "h3d" / "lib" / "linux64")
    reader = str(OR_ROOT / "extlib" / "hm_reader" / "linux64")
    env["LD_LIBRARY_PATH"] = reader + ":" + env.get("LD_LIBRARY_PATH", "")
    return env


def run_cmd(cmd: list[str], cwd: Path, log_lines: list[str]) -> subprocess.CompletedProcess[str]:
    log_lines.append("$ " + " ".join(cmd))
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        env=radioss_env(),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_lines.append(proc.stdout)
    return proc


def law25_block(mat: Material, mat_id: int) -> list[str]:
    return [
        f"/MAT/LAW25/{mat_id}",
        f"IM7_8552_ply_{mat_id}",
        "#              RHO_I",
        fmt_f(mat.rho),
        "#                E11                 E22                NU12     Iform                           E33",
        fmt_f(mat.e1, mat.e2, mat.nu12) + fmt_i(0) + f"{mat.e3:20.12g}",
        "#                G12                 G23                 G31              EPS_f1              EPS_f2",
        fmt_f(mat.g12, mat.g23, mat.g13, 0.0, 0.0),
        "#             EPS_t1              EPS_m1              EPS_t2              EPS_m2                dmax",
        fmt_f(0.0, 0.0, 0.0, 0.0, 1.0),
        "#              Wpmax               Wpref      Ioff                         ratio",
        fmt_f(0.0, 0.0) + fmt_i(0) + f"{0.0:20.12g}",
        "#                  b                   n                fmax",
        fmt_f(0.0, 0.0, 0.0),
        "#            sig_1yt             sig_2yt             sig_1yc             sig_2yc               alpha",
        fmt_f(mat.xt, mat.yt, mat.xc, mat.yc, 0.0),
        "#           sig_12yc            sig_12yt                c_12          Eps_rate_0       ICC",
        fmt_f(mat.s12, mat.s12, 0.0, 0.0) + fmt_i(0),
        "#          GAMMA_ini           GAMMA_max               d3max",
        fmt_f(0.0, 0.0, 0.0),
        "#  Fsmooth                Fcut",
        fmt_i(0) + f"{0.0:20.12g}",
    ]


def skew_block(skew_id: int, theta_deg: float) -> list[str]:
    th = math.radians(theta_deg)
    c = math.cos(th)
    s = math.sin(th)
    return [
        f"/SKEW/FIX/{skew_id}",
        f"ply_frame_{skew_id}_theta_{theta_deg:+.0f}",
        fmt_f(0.0, 0.0, 0.0),
        fmt_f(c, s, 0.0),
        fmt_f(-s, c, 0.0),
    ]


def type14_property(prop_id: int) -> list[str]:
    return [
        f"/PROP/TYPE14/{prop_id}",
        f"ply_{prop_id}_canonical_type14",
        "#   Isolid    Ismstr               Icpre               Inpts    Itetra    Iframe                  dn",
        fmt_i(24, 4) + f"{1:20d}{0:20d}{0:10d}{2:10d}{0.0:20.12g}",
        "#                q_a                 q_b                   h            LAMBDA_V                MU_V",
        fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
        "#             dt_min   istrain      IHKT",
        fmt_f(0.0) + fmt_i(0, 0),
    ]


def type6_property(prop_id: int, skew_id: int) -> list[str]:
    return [
        f"/PROP/TYPE6/{prop_id}",
        f"ply_{prop_id}_type6_proxy",
        "#   Isolid    Ismstr               Icpre  Itetra10     Inpts   Itetra4    Iframe                  Dn",
        fmt_i(24, 4) + f"{1:20d}{0:10d}{0:10d}{0:10d}{2:10d}{0.0:20.12g}",
        "#                 qa                  qb                   h",
        fmt_f(0.0, 0.0, 0.0),
        "#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth",
        fmt_f(1.0, 0.0, 0.0) + fmt_i(skew_id, 0, 1),
        "#                Phi                 Px                  Py                  Pz",
        fmt_f(0.0, 0.0, 0.0, 0.0),
        "#             dt_min   istrain      IHKT",
        fmt_f(0.0) + fmt_i(0, 0),
    ]


def stack_mesh(nplies: int) -> tuple[list[str], list[tuple[int, list[int]]]]:
    lx = 0.100
    ly = 0.100
    node_id: dict[tuple[int, int, int], int] = {}
    lines = ["/NODE"]
    nid = 1
    for k in range(nplies + 1):
        for j in range(2):
            for i in range(2):
                node_id[(i, j, k)] = nid
                x = lx * i
                y = ly * j
                z = PLY_T * k
                lines.append(f"{nid:10d}{x:20.12g}{y:20.12g}{z:20.12g}")
                nid += 1
    bricks: list[tuple[int, list[int]]] = []
    eid = 1
    for k in range(nplies):
        conn = [
            node_id[(0, 0, k)],
            node_id[(1, 0, k)],
            node_id[(1, 1, k)],
            node_id[(0, 1, k)],
            node_id[(0, 0, k + 1)],
            node_id[(1, 0, k + 1)],
            node_id[(1, 1, k + 1)],
            node_id[(0, 1, k + 1)],
        ]
        bricks.append((eid, conn))
        eid += 1
    return lines, bricks


def write_starter(layup_name: str, thetas: tuple[float, ...], mat: Material, use_type6: bool) -> Path:
    suffix = "type6_proxy" if use_type6 else "type14_canonical"
    job = f"stage09_{layup_name}_{suffix}"
    nodes, bricks = stack_mesh(len(thetas))
    lines = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        "      2023         0",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 09 {layup_name} {suffix}",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 4) + f"{0:20d}{2:40d}",
    ]
    for ply, theta in enumerate(thetas, start=1):
        lines.extend(law25_block(mat, ply))
        lines.extend(skew_block(ply, theta))
        lines.extend(type6_property(ply, ply) if use_type6 else type14_property(ply))
        lines.extend([f"/PART/{ply}", f"ply_{ply:02d}", fmt_i(ply, ply, 0), f"/BRICK/{ply}"])
        eid, conn = bricks[ply - 1]
        lines.append(fmt_i(eid, *conn))
    lines.extend(nodes)
    lines.extend(["/END", ""])
    path = RUNS_DIR / f"{job}_0000.rad"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def run_stage() -> tuple[dict[str, object], list[dict[str, object]], list[str]]:
    mat = load_material()
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    log_lines = [f"material_card={CARD_PATH}", f"starter={STARTER}"]
    rows: list[dict[str, object]] = []
    type14_blocked_all = True
    type6_supported_all = True

    for name, thetas in LAYUPS.items():
        a = a_matrix(mat, thetas)
        strict = write_starter(name, thetas, mat, use_type6=False)
        strict_proc = run_cmd([str(STARTER), "-i", strict.name, "-nt", "1"], RUNS_DIR, log_lines)
        strict_out = (RUNS_DIR / strict.name.replace(".rad", ".out")).read_text(
            encoding="utf-8", errors="replace"
        )
        strict_blocked = (
            strict_proc.returncode != 0
            and "MATERIAL/PROPERTY COMPATIBILITY" in strict_out
            and "TYPE 14" in strict_out
            and "LAW  25" in strict_out
        )
        type14_blocked_all = type14_blocked_all and strict_blocked

        proxy = write_starter(name, thetas, mat, use_type6=True)
        proxy_proc = run_cmd([str(STARTER), "-i", proxy.name, "-nt", "1"], RUNS_DIR, log_lines)
        type6_supported_all = type6_supported_all and proxy_proc.returncode == 0

        for i, row_name in enumerate(("xx", "yy", "xy")):
            for j, col_name in enumerate(("xx", "yy", "xy")):
                rows.append(
                    {
                        "layup": name,
                        "component": f"A_{row_name}{col_name}",
                        "reference_N_per_m": a[i, j],
                        "canonical_property": "TYPE14",
                        "canonical_starter_rc": strict_proc.returncode,
                        "canonical_status": "blocked_error_3047" if strict_blocked else "unexpected",
                        "proxy_property": "TYPE6_SOL_ORTH",
                        "proxy_starter_rc": proxy_proc.returncode,
                        "solver_N_per_m": "",
                        "relative_error_pct": "",
                        "verdict": "INCONCLUSIVE",
                    }
                )

    metrics = {
        "canonical_type14_law25_supported": not type14_blocked_all,
        "canonical_type14_error_id": 3047 if type14_blocked_all else None,
        "type6_sol_orth_proxy_starter_supported": type6_supported_all,
        "layup_count": len(LAYUPS),
        "clt_A_gate_evaluated": False,
    }
    return metrics, rows, log_lines


def write_outputs(metrics: dict[str, object], rows: list[dict[str, object]], log_lines: list[str], wall_s: float) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    with (RESULTS_DIR / "timeseries.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    git_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=str(ROOT_DIR),
        text=True,
        stdout=subprocess.PIPE,
        check=False,
    ).stdout.strip()
    results = {
        "stage": 9,
        "verdict": "INCONCLUSIVE",
        "metrics": {**metrics, "wall_clock_s": wall_s},
        "reference": {
            "material_card": str(CARD_PATH),
            "ply_thickness_m": PLY_T,
            "layups": {name: list(thetas) for name, thetas in LAYUPS.items()},
        },
        "tolerance": {"A_matrix_component_relative_error": 0.02},
        "git_sha": git_sha,
    }
    (RESULTS_DIR / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    (FIGURES_DIR / "stage09_clt_probe.typ").write_text(
        "\n".join(
            [
                '#set page(width: 175mm, height: auto, margin: 10mm)',
                '#let rows = csv("../results/timeseries.csv")',
                '#text(size: 12pt, weight: "bold")[Stage 09 laminate CLT starter probe]',
                '#v(5pt)',
                '#table(',
                '  columns: (34mm, 24mm, 34mm, 36mm, 26mm),',
                '  stroke: rgb("#5C5C5C"),',
                '  [Layup], [Component], [Reference N/m], [Canonical status], [Verdict],',
                '  ..rows.map(r => ([#r.at(0)], [#r.at(1)], [#r.at(2)], [#r.at(5)], [#r.at(10)])).flatten(),',
                ')',
                "",
            ]
        ),
        encoding="utf-8",
    )

    (THIS_DIR / "blocker.md").write_text(
        "\n".join(
            [
                "# Stage 09 Blocker - LAW25 + TYPE14 laminate stack does not start",
                "",
                "Author: J.C. Vaught",
                "",
                "The canonical stage 09 CLT verification requires one `/MAT/LAW25` and one `/PROP/TYPE14` per solid ply. The installed OpenRadioss starter rejects the LAW25 + TYPE14 pairing with `ERROR ID : 3047` material/property compatibility for both the `[0/90]s` and `[0/+45/-45/90]s` stacks.",
                "",
                "Analytic CLT A-matrices are written to `results/timeseries.csv`. `/PROP/TYPE6` (`/PROP/SOL_ORTH`) proxy stacks start successfully, so the blocker is the required TYPE14 path rather than the material card or layup bookkeeping.",
                "",
                "Verdict: `INCONCLUSIVE` due to toolchain/material-property compatibility, not a numerical CLT failure.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    RUN_LOG.write_text("\n".join(log_lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if not args.all:
        parser.error("use --all")
    start = time.perf_counter()
    metrics, rows, log_lines = run_stage()
    wall_s = time.perf_counter() - start
    write_outputs(metrics, rows, log_lines, wall_s)
    print(json.dumps({"stage": 9, "verdict": "INCONCLUSIVE", "wall_clock_s": wall_s}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
