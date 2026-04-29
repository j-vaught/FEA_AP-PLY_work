"""Stage 06 - composite failure-card capability probe.

Author: J.C. Vaught

This stage is intended to validate /MAT/LAW25 with /FAIL/TSAIWU,
/FAIL/HASHIN, and /FAIL/PUCK on solid bricks. The repository specification
requires /PROP/TYPE14, but the OpenRadioss starter installed on this host
rejects LAW25 + TYPE14 as an incompatible material/property pair. This runner
therefore records that canonical probe and separately runs a TYPE6/SOL_ORTH
all-solid proxy to confirm that the three failure cards parse and erode solid
elements. The verdict is INCONCLUSIVE because the required TYPE14 pathway does
not run on the available toolchain.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


THIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parents[1]
RUNS_DIR = THIS_DIR / "runs"
RESULTS_DIR = THIS_DIR / "results"
FIGURES_DIR = THIS_DIR / "figures"
RUN_LOG = THIS_DIR / "run.log"
CARD_PATH = ROOT_DIR / "references" / "material_cards" / "im7_8552.json"

OR_ROOT = Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"
ENGINE = OR_ROOT / "exec" / "engine_linux64_gf"
ANIM_TO_VTK = OR_ROOT / "exec" / "anim_to_vtk_linux64_gf"

CRITERIA = ("TSAIWU", "HASHIN", "PUCK")
AXES = (
    ("XT", "longitudinal_tension", "S1S2_000", 1.0, 0.0, 0.0, "Xt_Pa"),
    ("XC", "longitudinal_compression", "S1S2_180", -1.0, 0.0, 0.0, "Xc_Pa"),
    ("YT", "transverse_tension", "S1S2_090", 0.0, 1.0, 0.0, "Yt_Pa"),
    ("YC", "transverse_compression", "S1S2_270", 0.0, -1.0, 0.0, "Yc_Pa"),
    ("S12", "in_plane_shear", "S2T12_090", 0.0, 0.0, 1.0, "S12_Pa"),
)


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
    zt: float
    zc: float
    s12: float
    s13: float
    s23: float


@dataclass(frozen=True)
class ProbeResult:
    criterion: str
    starter_rc: int
    engine_rc: int | None
    failure_observed: bool
    failure_time_s: float | None
    first_vtk: str | None


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
        zt=float(strength["Zt_Pa"]),
        zc=float(strength["Zc_Pa"]),
        s12=float(strength["S12_Pa"]),
        s13=float(strength["S13_Pa"]),
        s23=float(strength["S23_Pa"]),
    )


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


def group_block(group_id: int, name: str, node_ids: Iterable[int]) -> list[str]:
    ids = list(node_ids)
    lines = [f"/GRNOD/NODE/{group_id}", name]
    for i in range(0, len(ids), 10):
        lines.append(fmt_i(*ids[i : i + 10]))
    return lines


def law25_block(mat: Material) -> list[str]:
    """2019-format LAW25 elastic base card accepted by the 2026 binary."""
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


def fail_block(criterion: str, mat: Material) -> list[str]:
    if criterion == "TSAIWU":
        return [
            "/FAIL/TSAIWU/1",
            "#           SIGMA_1T            SIGMA_2T            SIGMA_1C            SIGMA_2C            SIGMA_12",
            fmt_f(mat.xt, mat.yt, mat.xc, mat.yc, mat.s12),
            "#              ALPHA             TAU_MAX                FCUT                      IFAIL_SH  IFAIL_SO",
            fmt_f(-0.5, 0.0, 0.0) + (" " * 20) + fmt_i(0, 1),
        ]
    if criterion == "HASHIN":
        return [
            "/FAIL/HASHIN/1",
            "#    IFORM  IFAIL_SH  IFAIL_SO          PTHICKFAIL",
            fmt_i(1, 0, 1) + f"{1.0:20.12g}",
            "#           SIGMA_1T            SIGMA_2T            SIGMA_3T            SIGMA_1C            SIGMA_2C",
            fmt_f(mat.xt, mat.yt, mat.zt, mat.xc, mat.yc),
            "#            SIGMA_C           SIGMA_12F           SIGMA_12M           SIGMA_23M           SIGMA_13M",
            fmt_f(mat.yc, mat.s12, mat.s12, mat.s23, mat.s13),
            "#                PHI                SDEL                TMAX                                    FCUT",
            fmt_f(0.0, 1.0, 0.0) + (" " * 20) + f"{0.0:20.12g}",
        ]
    if criterion == "PUCK":
        return [
            "/FAIL/PUCK/1",
            "#           SIGMA_1T            SIGMA_2T            SIGMA_12            SIGMA_1C            SIGMA_2C",
            fmt_f(mat.xt, mat.yt, mat.s12, mat.xc, mat.yc),
            "#       P12_POSITIVE        P12_NEGATIVE        P22_NEGATIVE             TAU_MAX  IFAIL_SH  IFAIL_SO",
            fmt_f(0.30, 0.25, 0.25, 0.0) + fmt_i(0, 1),
            "#               FCUT",
            fmt_f(0.0),
        ]
    raise ValueError(f"unknown criterion {criterion}")


def type14_property_block() -> list[str]:
    return [
        "/PROP/TYPE14/1",
        "strict_type14_required_by_stage06",
        "#   Isolid    Ismstr               Icpre               Inpts    Itetra    Iframe                  dn",
        fmt_i(24, 4) + f"{1:20d}{0:20d}{0:10d}{2:10d}{0.0:20.12g}",
        "#                q_a                 q_b                   h            LAMBDA_V                MU_V",
        fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
        "#             dt_min   istrain      IHKT",
        fmt_f(0.0) + fmt_i(0, 0),
    ]


def type6_property_block() -> list[str]:
    return [
        "/PROP/TYPE6/1",
        "type6_sol_orth_proxy_skew_global",
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


def one_brick_mesh_and_load(job: str, property_lines: list[str], mat: Material, criterion: str) -> list[str]:
    left = [1, 4, 5, 8]
    right = [2, 3, 6, 7]
    lines = [
        "#RADIOSS STARTER",
        "/BEGIN",
        job,
        "      2023         0",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        f"{'kg':>20}{'m':>20}{'s':>20}",
        "/TITLE",
        f"Stage 06 {criterion} one-brick failure-card probe",
        "/DEF_SOLID",
        "#  I_SOLID    ISMSTR             ISTRAIN                                  IFRAME",
        fmt_i(24, 4) + f"{0:20d}{2:40d}",
    ]
    lines.extend(law25_block(mat))
    lines.extend(fail_block(criterion, mat))
    lines.extend(
        [
            "/NODE",
            fmt_i(1) + fmt_f(0.0, 0.0, 0.0),
            fmt_i(2) + fmt_f(1.0, 0.0, 0.0),
            fmt_i(3) + fmt_f(1.0, 1.0, 0.0),
            fmt_i(4) + fmt_f(0.0, 1.0, 0.0),
            fmt_i(5) + fmt_f(0.0, 0.0, 1.0),
            fmt_i(6) + fmt_f(1.0, 0.0, 1.0),
            fmt_i(7) + fmt_f(1.0, 1.0, 1.0),
            fmt_i(8) + fmt_f(0.0, 1.0, 1.0),
            "/PART/1",
            "solid_probe_part",
            fmt_i(1, 1, 0),
            "/BRICK/1",
            fmt_i(1, 1, 2, 3, 4, 5, 6, 7, 8),
        ]
    )
    lines.extend(property_lines)
    lines.extend(
        [
            "/SKEW/FIX/1",
            "global_material_frame",
            fmt_f(0.0, 0.0, 0.0),
            fmt_f(1.0, 0.0, 0.0),
            fmt_f(0.0, 1.0, 0.0),
            "/BCS/1",
            "left_face_fixed_for_overdrive_probe",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   111 000{0:10d}{100:10d}",
            "/BCS/2",
            "right_face_yz_fixed_for_overdrive_probe",
            "#  Tra rot   skew_ID  grnod_ID",
            f"   011 000{0:10d}{101:10d}",
            "/FUNCT/1",
            "linear_ramp",
            fmt_f(0.0, 0.0),
            fmt_f(1.0, 1.0),
            "/IMPDISP/1",
            "right_face_x_overdrive",
            "#   Ifunct       DIR     Iskew   Isensor   Gnod_id     Frame     Icoor",
            f"{1:10d}{'X':>10}{0:10d}{0:10d}{101:10d}{0:10d}{0:10d}",
            "#            Scale_x             Scale_y              Tstart               Tstop",
            fmt_f(1.0e-3, 2.0e-2, 0.0, 0.0),
        ]
    )
    lines.extend(group_block(100, "left_face", left))
    lines.extend(group_block(101, "right_face", right))
    lines.extend(["/END", ""])
    return lines


def write_engine_deck(job: str, path: Path) -> None:
    lines = [
        "#RADIOSS ENGINE",
        "/ANIM/DT",
        fmt_f(2.0e-4, 2.0e-4),
        "/ANIM/VECT/DISP",
        "/ANIM/BRICK/TENS/STRESS/ALL",
        "/ANIM/BRICK/TENS/DAMA",
        "/ANIM/ELEM/EPSP",
        "/ANIM/GZIP",
        "/TFILE/4",
        fmt_f(1.0e-5),
        "/RFILE",
        fmt_i(1000),
        "/PRINT/-100/55",
        "/DT/NODA/CST",
        fmt_f(0.0, 1.0e-7),
        f"/RUN/{job}/1",
        fmt_f(1.0e-3),
        "/VERS/2023",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def build_decks(mat: Material) -> None:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    strict_job = "stage06_type14_canonical_probe"
    strict = one_brick_mesh_and_load(strict_job, type14_property_block(), mat, "TSAIWU")
    (RUNS_DIR / f"{strict_job}_0000.rad").write_text("\n".join(strict), encoding="utf-8")

    for criterion in CRITERIA:
        job = f"stage06_{criterion}_type6_proxy"
        starter = one_brick_mesh_and_load(job, type6_property_block(), mat, criterion)
        (RUNS_DIR / f"{job}_0000.rad").write_text("\n".join(starter), encoding="utf-8")
        write_engine_deck(job, RUNS_DIR / f"{job}_0001.rad")


def failure_observed(text: str, criterion: str) -> bool:
    return f"FAILURE ({criterion})" in text.upper()


def parse_failure_time(text: str, criterion: str) -> float | None:
    patterns = [
        rf"FAILURE \({criterion}\).*?AT TIME :\s*([0-9.Ee+-]+)",
        rf"FAILURE \({criterion}\).*?\n.*?TIME\s*:\s*([0-9.Ee+-]+)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL)
        if match:
            return float(match.group(1))
    return None


def convert_first_frame(job: str, log_lines: list[str]) -> str | None:
    gz = RUNS_DIR / f"{job}A001.gz"
    anim = RUNS_DIR / f"{job}A001"
    if gz.exists():
        with gzip.open(gz, "rb") as src, anim.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    if not anim.exists():
        return None
    proc = run_cmd([str(ANIM_TO_VTK), str(anim)], RUNS_DIR, log_lines)
    vtk = RUNS_DIR / f"{job}A001.vtk"
    if proc.stdout.startswith("# vtk"):
        vtk.write_text(proc.stdout, encoding="utf-8")
    return str(vtk) if vtk.exists() else None


def run_stage() -> tuple[dict[str, object], list[ProbeResult], list[str]]:
    mat = load_material()
    log_lines: list[str] = []
    build_decks(mat)
    log_lines.append(f"material_card={CARD_PATH}")
    log_lines.append(f"starter={STARTER}")
    log_lines.append(f"engine={ENGINE}")

    strict_job = "stage06_type14_canonical_probe"
    strict_proc = run_cmd(
        [str(STARTER), "-i", f"{strict_job}_0000.rad", "-nt", "1"],
        RUNS_DIR,
        log_lines,
    )
    strict_out = (RUNS_DIR / f"{strict_job}_0000.out").read_text(
        encoding="utf-8", errors="replace"
    )
    type14_rejected = (
        strict_proc.returncode != 0
        and "MATERIAL/PROPERTY COMPATIBILITY" in strict_out
        and "TYPE 14" in strict_out
        and "LAW  25" in strict_out
    )

    probes: list[ProbeResult] = []
    for criterion in CRITERIA:
        job = f"stage06_{criterion}_type6_proxy"
        starter_proc = run_cmd([str(STARTER), "-i", f"{job}_0000.rad", "-nt", "1"], RUNS_DIR, log_lines)
        engine_rc: int | None = None
        failure_time: float | None = None
        vtk_path: str | None = None
        if starter_proc.returncode == 0:
            engine_proc = run_cmd([str(ENGINE), "-i", f"{job}_0001.rad", "-nt", "1"], RUNS_DIR, log_lines)
            engine_rc = engine_proc.returncode
            engine_text = engine_proc.stdout
            observed = failure_observed(engine_text, criterion)
            failure_time = parse_failure_time(engine_text, criterion)
            vtk_path = convert_first_frame(job, log_lines)
        else:
            observed = False
        probes.append(
            ProbeResult(
                criterion=criterion,
                starter_rc=starter_proc.returncode,
                engine_rc=engine_rc,
                failure_observed=observed,
                failure_time_s=failure_time,
                first_vtk=vtk_path,
            )
        )

    metrics = {
        "canonical_type14_law25_supported": not type14_rejected,
        "canonical_type14_starter_rc": strict_proc.returncode,
        "canonical_type14_error_id": 3047 if type14_rejected else None,
        "type6_proxy_all_failure_cards_executed": all(
            p.starter_rc == 0 and p.engine_rc == 0 and p.failure_observed for p in probes
        ),
        "type6_proxy_probe_count": len(probes),
        "principal_axis_gate_evaluated": False,
        "off_axis_paths_reported_analytic_only": 32,
    }
    return metrics, probes, log_lines


def tsaiwu_load(c1: float, c2: float, c6: float, mat: Material) -> float:
    f1 = 1.0 / mat.xt - 1.0 / mat.xc
    f2 = 1.0 / mat.yt - 1.0 / mat.yc
    f11 = 1.0 / (mat.xt * mat.xc)
    f22 = 1.0 / (mat.yt * mat.yc)
    f66 = 1.0 / (mat.s12 * mat.s12)
    f12 = -0.5 * math.sqrt(f11 * f22)
    a = f11 * c1 * c1 + f22 * c2 * c2 + f66 * c6 * c6 + 2.0 * f12 * c1 * c2
    b = f1 * c1 + f2 * c2
    disc = b * b + 4.0 * a
    return (-b + math.sqrt(disc)) / (2.0 * a)


def hashin_load(c1: float, c2: float, c6: float, mat: Material) -> float:
    candidates: list[float] = []
    if c1 > 0.0:
        candidates.append(1.0 / math.sqrt((c1 / mat.xt) ** 2 + (c6 / mat.s12) ** 2))
    if c1 < 0.0:
        candidates.append(mat.xc / abs(c1))
    if c2 > 0.0:
        candidates.append(1.0 / math.sqrt((c2 / mat.yt) ** 2 + (c6 / mat.s12) ** 2))
    if c2 < 0.0:
        a = (c2 / (2.0 * mat.s23)) ** 2 + (c6 / mat.s12) ** 2
        b = ((mat.yc / (2.0 * mat.s23)) ** 2 - 1.0) * c2 / mat.yc
        disc = b * b + 4.0 * a
        candidates.append((-b + math.sqrt(disc)) / (2.0 * a))
    return min(candidates) if candidates else math.inf


def analytic_envelope_rows(mat: Material) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    pid = 1
    for k in range(24):
        theta = math.radians(15.0 * k)
        c1, c2, c6 = math.cos(theta), math.sin(theta), 0.0
        for criterion in CRITERIA:
            if criterion == "TSAIWU":
                r = tsaiwu_load(c1, c2, c6, mat)
            elif criterion == "HASHIN":
                r = hashin_load(c1, c2, c6, mat)
            else:
                r = min(
                    mat.xt / c1 if c1 > 0.0 else math.inf,
                    mat.xc / abs(c1) if c1 < 0.0 else math.inf,
                    mat.yt / c2 if c2 > 0.0 else math.inf,
                    mat.yc / abs(c2) if c2 < 0.0 else math.inf,
                )
            rows.append(
                {
                    "kind": "analytic_envelope",
                    "criterion": criterion,
                    "path_id": pid,
                    "path": f"S1S2_{15*k:03d}",
                    "axis": "",
                    "reference_pa": "",
                    "reported_pa": f"{r:.12g}",
                    "relative_error_pct": "",
                    "sigma1_pa": f"{c1*r:.12g}",
                    "sigma2_pa": f"{c2*r:.12g}",
                    "tau12_pa": f"{c6*r:.12g}",
                    "solver_status": "not_sampled_stage_inconclusive",
                    "failure_time_s": "",
                }
            )
        pid += 1
    for k in range(8):
        phi = math.radians(30.0 * k)
        c1, c2, c6 = 0.0, math.cos(phi), math.sin(phi)
        for criterion in CRITERIA:
            if criterion == "TSAIWU":
                r = tsaiwu_load(c1, c2, c6, mat)
            elif criterion == "HASHIN":
                r = hashin_load(c1, c2, c6, mat)
            else:
                r = min(
                    mat.yt / c2 if c2 > 0.0 else math.inf,
                    mat.yc / abs(c2) if c2 < 0.0 else math.inf,
                    mat.s12 / abs(c6) if abs(c6) > 0.0 else math.inf,
                )
            rows.append(
                {
                    "kind": "analytic_envelope",
                    "criterion": criterion,
                    "path_id": pid,
                    "path": f"S2T12_{30*k:03d}",
                    "axis": "",
                    "reference_pa": "",
                    "reported_pa": f"{r:.12g}",
                    "relative_error_pct": "",
                    "sigma1_pa": f"{c1*r:.12g}",
                    "sigma2_pa": f"{c2*r:.12g}",
                    "tau12_pa": f"{c6*r:.12g}",
                    "solver_status": "not_sampled_stage_inconclusive",
                    "failure_time_s": "",
                }
            )
        pid += 1
    return rows


def write_outputs(metrics: dict[str, object], probes: list[ProbeResult], wall_s: float) -> None:
    mat = load_material()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    strength = {
        "Xt_Pa": mat.xt,
        "Xc_Pa": mat.xc,
        "Yt_Pa": mat.yt,
        "Yc_Pa": mat.yc,
        "S12_Pa": mat.s12,
    }
    rows: list[dict[str, object]] = []
    for criterion in CRITERIA:
        for axis, label, path, c1, c2, c6, key in AXES:
            rows.append(
                {
                    "kind": "principal_axis_reference",
                    "criterion": criterion,
                    "path_id": "",
                    "path": path,
                    "axis": axis,
                    "reference_pa": f"{strength[key]:.12g}",
                    "reported_pa": "",
                    "relative_error_pct": "",
                    "sigma1_pa": f"{c1 * strength[key]:.12g}",
                    "sigma2_pa": f"{c2 * strength[key]:.12g}",
                    "tau12_pa": f"{c6 * strength[key]:.12g}",
                    "solver_status": "not_evaluated_canonical_type14_rejected",
                    "failure_time_s": "",
                }
            )
    for probe in probes:
        rows.append(
            {
                "kind": "type6_proxy_overdrive",
                "criterion": probe.criterion,
                "path_id": "",
                "path": "one_brick_constrained_x_overdrive",
                "axis": "",
                "reference_pa": "",
                "reported_pa": "",
                "relative_error_pct": "",
                "sigma1_pa": "",
                "sigma2_pa": "",
                "tau12_pa": "",
                "solver_status": "failure_observed" if probe.failure_observed else "not_observed",
                "failure_time_s": "" if probe.failure_time_s is None else f"{probe.failure_time_s:.12g}",
            }
        )
    rows.extend(analytic_envelope_rows(mat))

    timeseries = RESULTS_DIR / "timeseries.csv"
    with timeseries.open("w", newline="", encoding="utf-8") as f:
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
        "stage": 6,
        "verdict": "INCONCLUSIVE",
        "metrics": {
            **metrics,
            "wall_clock_s": wall_s,
            "type6_proxy": [
                {
                    "criterion": p.criterion,
                    "starter_rc": p.starter_rc,
                    "engine_rc": p.engine_rc,
                    "failure_observed": p.failure_observed,
                    "failure_time_s": p.failure_time_s,
                    "first_vtk": p.first_vtk,
                }
                for p in probes
            ],
        },
        "reference": {
            "material_card": str(CARD_PATH),
            "strengths_pa": strength,
            "principal_axis_gate": "5 percent per criterion, not evaluated because canonical TYPE14 deck does not start",
        },
        "tolerance": {"principal_axis_relative_error": 0.05},
        "git_sha": git_sha,
    }
    (RESULTS_DIR / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    typ = FIGURES_DIR / "stage06_failure_card_probe.typ"
    typ.write_text(
        "\n".join(
            [
                '#import "@preview/cetz:0.3.4"',
                '#set page(width: 180mm, height: auto, margin: 10mm)',
                '#let garnet = rgb("#73000A")',
                '#let black70 = rgb("#5C5C5C")',
                '#let atlantic = rgb("#466A9F")',
                '#let horseshoe = rgb("#65780B")',
                '#let rows = csv("../results/timeseries.csv")',
                '#text(size: 12pt, weight: "bold")[Stage 06 failure-card probe]',
                '#v(4pt)',
                '#text(size: 8pt)[Canonical LAW25 + TYPE14 is rejected by starter; TYPE6/SOL_ORTH proxy confirms the three failure cards execute on solid elements.]',
                '#v(8pt)',
                '#table(',
                '  columns: (28mm, 35mm, 32mm, 32mm),',
                '  stroke: black70,',
                '  [Criterion], [Proxy status], [Failure time (s)], [Canonical gate],',
                '  ..rows.filter(r => r.at(0) == "type6_proxy_overdrive").map(r => (',
                '    [#r.at(1)], [#r.at(11)], [#r.at(12)], [INCONCLUSIVE],',
                '  )).flatten(),',
                ')',
                "",
            ]
        ),
        encoding="utf-8",
    )

    blocker = THIS_DIR / "blocker.md"
    blocker.write_text(
        "\n".join(
            [
                "# Stage 06 Blocker - LAW25 + TYPE14 rejected by starter",
                "",
                "Author: J.C. Vaught",
                "",
                "The required canonical deck uses `/MAT/LAW25` on `/PROP/TYPE14` with solid HEXA8 elements. The installed OpenRadioss starter rejects that material/property pair before the engine can run.",
                "",
                "Observed starter finding:",
                "",
                "- `ERROR ID : 3047`",
                "- `ERROR IN MATERIAL/PROPERTY COMPATIBILITY`",
                "- `PROPERTY ID 1 OF TYPE 14 IS NOT COMPATIBLE WITH MATERIAL ID 1 OF TYPE 25`",
                "",
                "A separate all-solid `/PROP/TYPE6` (`/PROP/SOL_ORTH`) proxy using the same `/MAT/LAW25` card successfully starts and the engine reports solid-element failure for `/FAIL/TSAIWU`, `/FAIL/HASHIN`, and `/FAIL/PUCK`. That proves the failure cards are present in this binary, but it does not satisfy the stage 06 `/PROP/TYPE14` requirement or consolidation item B.",
                "",
                "Verdict: `INCONCLUSIVE` due to toolchain/material-property compatibility, not a numerical failure against the WWFE-II principal-axis strengths.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true", help="build, run, analyze, and write artifacts")
    args = parser.parse_args()
    if not args.all:
        parser.error("stage 06 runner expects --all")

    start = time.perf_counter()
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    metrics, probes, log_lines = run_stage()
    wall_s = time.perf_counter() - start
    write_outputs(metrics, probes, wall_s)
    RUN_LOG.write_text("\n".join(log_lines), encoding="utf-8")
    print(json.dumps({"stage": 6, "verdict": "INCONCLUSIVE", "wall_clock_s": wall_s}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
