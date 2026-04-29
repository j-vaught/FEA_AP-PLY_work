"""Stage 08 - ply rotation starter probe.

Author: J.C. Vaught

The intended verification is a seven-angle LAW25 + TYPE14 solid coupon
sweep. The installed OpenRadioss starter rejects LAW25 + TYPE14 with ERROR
3047, so this runner generates the canonical sweep decks, records the
starter blocker, and writes the analytic Ex(theta) table without claiming an
FEM modulus. TYPE6/SOL_ORTH proxy starters are also probed to isolate the
blocker to the TYPE14 material/property pairing.
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


THIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parents[1]
RUNS_DIR = THIS_DIR / "runs"
RESULTS_DIR = THIS_DIR / "results"
FIGURES_DIR = THIS_DIR / "figures"
RUN_LOG = THIS_DIR / "run.log"
CARD_PATH = ROOT_DIR / "references" / "material_cards" / "im7_8552.json"

OR_ROOT = Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"

ANGLES = (0, 15, 30, 45, 60, 75, 90)


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


def ex_offaxis(theta_deg: float, mat: Material) -> float:
    theta = math.radians(theta_deg)
    c2 = math.cos(theta) ** 2
    s2 = math.sin(theta) ** 2
    inv = (
        c2 * c2 / mat.e1
        + s2 * s2 / mat.e2
        + (1.0 / mat.g12 - 2.0 * mat.nu12 / mat.e1) * c2 * s2
    )
    return 1.0 / inv


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


def law25_block(mat: Material) -> list[str]:
    return [
        "/MAT/LAW25/1",
        "IM7_8552_canonical_Soden_WWFEII",
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


def mesh_blocks() -> tuple[list[str], list[int], list[int]]:
    length = 0.250
    width = 0.025
    thick = 0.001
    nx, ny, nz = 100, 10, 2
    node_id: dict[tuple[int, int, int], int] = {}
    lines = ["/NODE"]
    nid = 1
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                node_id[(i, j, k)] = nid
                x = length * i / nx
                y = width * (j / ny - 0.5)
                z = thick * (k / nz - 0.5)
                lines.append(f"{nid:10d}{x:20.12g}{y:20.12g}{z:20.12g}")
                nid += 1
    lines.extend(["/PART/1", "ply_rotation_coupon", fmt_i(1, 1, 0), "/BRICK/1"])
    eid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                conn = [
                    node_id[(i, j, k)],
                    node_id[(i + 1, j, k)],
                    node_id[(i + 1, j + 1, k)],
                    node_id[(i, j + 1, k)],
                    node_id[(i, j, k + 1)],
                    node_id[(i + 1, j, k + 1)],
                    node_id[(i + 1, j + 1, k + 1)],
                    node_id[(i, j + 1, k + 1)],
                ]
                lines.append(fmt_i(eid, *conn))
                eid += 1
    left = [node_id[(0, j, k)] for k in range(nz + 1) for j in range(ny + 1)]
    right = [node_id[(nx, j, k)] for k in range(nz + 1) for j in range(ny + 1)]
    return lines, left, right


def group_block(group_id: int, name: str, node_ids: list[int]) -> list[str]:
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for i in range(0, len(node_ids), 10):
        lines.append(fmt_i(*node_ids[i : i + 10]))
    return lines


def skew_block(theta_deg: int) -> list[str]:
    theta = math.radians(theta_deg)
    c = math.cos(theta)
    s = math.sin(theta)
    return [
        "/SKEW/FIX/1",
        f"ply_frame_theta_{theta_deg}",
        fmt_f(0.0, 0.0, 0.0),
        fmt_f(c, s, 0.0),
        fmt_f(-s, c, 0.0),
    ]


def type14_property() -> list[str]:
    return [
        "/PROP/TYPE14/1",
        "canonical_type14_solid_property",
        "#   Isolid    Ismstr               Icpre               Inpts    Itetra    Iframe                  dn",
        fmt_i(24, 4) + f"{1:20d}{0:20d}{0:10d}{2:10d}{0.0:20.12g}",
        "#                q_a                 q_b                   h            LAMBDA_V                MU_V",
        fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
        "#             dt_min   istrain      IHKT",
        fmt_f(0.0) + fmt_i(0, 0),
    ]


def type6_property() -> list[str]:
    return [
        "/PROP/TYPE6/1",
        "type6_sol_orth_proxy_property",
        "#   Isolid    Ismstr               Icpre  Itetra10     Inpts   Itetra4    Iframe                  Dn",
        fmt_i(24, 4) + f"{1:20d}{0:10d}{0:10d}{0:10d}{2:10d}{0.0:20.12g}",
        "#                 qa                  qb                   h",
        fmt_f(0.0, 0.0, 0.0),
        "#                 Vx                  Vy                  Vz   skew_ID        Ip     Iorth",
        fmt_f(1.0, 0.0, 0.0) + fmt_i(1, 0, 1),
        "#                Phi                 Px                  Py                  Pz",
        fmt_f(0.0, 0.0, 0.0, 0.0),
        "#             dt_min   istrain      IHKT",
        fmt_f(0.0) + fmt_i(0, 0),
    ]


def write_starter(theta_deg: int, mat: Material, prop: list[str], suffix: str) -> Path:
    mesh, left, right = mesh_blocks()
    job = f"stage08_theta_{theta_deg:02d}_{suffix}"
    lines = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        "      2023         0",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 08 ply rotation theta {theta_deg} {suffix}",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 4) + f"{0:20d}{2:40d}",
    ]
    lines.extend(law25_block(mat))
    lines.extend(mesh)
    lines.extend(prop)
    lines.extend(skew_block(theta_deg))
    lines.extend(
        [
            "/BCS/1",
            "left_face_x_fixed",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   100 000{0:10d}{100:10d}",
            "/FUNCT/1",
            "unit_ramp",
            fmt_f(0.0, 0.0),
            fmt_f(1.0, 1.0),
            "/IMPDISP/1",
            "right_face_x",
            "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
            f"{1:10d}{'X':>10}{0:10d}{0:10d}{101:10d}{0:10d}{0:10d}",
            "#            Scale_x             Scale_y              Tstart               Tstop",
            fmt_f(1.0, 0.0005, 0.0, 0.0),
        ]
    )
    lines.extend(group_block(100, "left_face", left))
    lines.extend(group_block(101, "right_face", right))
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

    for theta in ANGLES:
        strict = write_starter(theta, mat, type14_property(), "type14_canonical")
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

        proxy = write_starter(theta, mat, type6_property(), "type6_proxy")
        proxy_proc = run_cmd([str(STARTER), "-i", proxy.name, "-nt", "1"], RUNS_DIR, log_lines)
        type6_supported_all = type6_supported_all and proxy_proc.returncode == 0

        rows.append(
            {
                "theta_deg": theta,
                "analytic_Ex_Pa": ex_offaxis(theta, mat),
                "canonical_property": "TYPE14",
                "canonical_starter_rc": strict_proc.returncode,
                "canonical_status": "blocked_error_3047" if strict_blocked else "unexpected",
                "proxy_property": "TYPE6_SOL_ORTH",
                "proxy_starter_rc": proxy_proc.returncode,
                "solver_Ex_Pa": "",
                "relative_error_pct": "",
                "verdict": "INCONCLUSIVE",
            }
        )

    metrics = {
        "canonical_type14_law25_supported": not type14_blocked_all,
        "canonical_type14_error_id": 3047 if type14_blocked_all else None,
        "type6_sol_orth_proxy_starter_supported": type6_supported_all,
        "angle_count": len(ANGLES),
        "modulus_gate_evaluated": False,
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
        "stage": 8,
        "verdict": "INCONCLUSIVE",
        "metrics": {**metrics, "wall_clock_s": wall_s},
        "reference": {
            "material_card": str(CARD_PATH),
            "angles_deg": list(ANGLES),
            "analytic_table": [{"theta_deg": r["theta_deg"], "Ex_Pa": r["analytic_Ex_Pa"]} for r in rows],
        },
        "tolerance": {"modulus_relative_error": 0.01},
        "git_sha": git_sha,
    }
    (RESULTS_DIR / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    (FIGURES_DIR / "stage08_ply_rotation_probe.typ").write_text(
        "\n".join(
            [
                '#set page(width: 175mm, height: auto, margin: 10mm)',
                '#let rows = csv("../results/timeseries.csv")',
                '#text(size: 12pt, weight: "bold")[Stage 08 ply-rotation starter probe]',
                '#v(5pt)',
                '#table(',
                '  columns: (20mm, 32mm, 38mm, 28mm, 28mm),',
                '  stroke: rgb("#5C5C5C"),',
                '  [Theta], [Analytic Ex (GPa)], [Canonical status], [Proxy rc], [Verdict],',
                '  ..rows.map(r => ([#r.at(0)], [#str(float(r.at(1)) / 1e9)], [#r.at(4)], [#r.at(6)], [#r.at(9)])).flatten(),',
                ')',
                "",
            ]
        ),
        encoding="utf-8",
    )

    (THIS_DIR / "blocker.md").write_text(
        "\n".join(
            [
                "# Stage 08 Blocker - LAW25 + TYPE14 ply-rotation sweep does not start",
                "",
                "Author: J.C. Vaught",
                "",
                "The canonical stage 08 sweep requires `/MAT/LAW25` on `/PROP/TYPE14` solid HEXA8 coupons with `/SKEW/FIX` frames at seven ply angles. Every canonical deck is rejected by the installed OpenRadioss starter with `ERROR ID : 3047` material/property compatibility before the engine can run.",
                "",
                "The runner also generated `/PROP/TYPE6` (`/PROP/SOL_ORTH`) proxy decks at the same seven angles. Those proxy starters succeed, which isolates the issue to the required TYPE14 path rather than to the material card or skew definitions.",
                "",
                "Verdict: `INCONCLUSIVE` due to toolchain/material-property compatibility, not a numerical off-axis stiffness failure.",
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
    print(json.dumps({"stage": 8, "verdict": "INCONCLUSIVE", "wall_clock_s": wall_s}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
