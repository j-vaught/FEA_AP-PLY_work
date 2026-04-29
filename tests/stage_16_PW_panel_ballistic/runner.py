#!/usr/bin/env python3
# coding: utf-8
"""
Stage 16 runner: Pseudo-Woven Panel Ballistic (UofSC P1-TWT, Vakili Rad 2020).

This is the FEA_AP-PLY project's terminal-goal driver. It executes the Kok
ap_ply_model_creation -> GMSH -> Abaqus .inp -> inp2rad -> OpenRadioss .rad
pipeline, then runs a V50 bracketing sweep (seven shots), parses the per-shot
residual velocity from the OpenRadioss T01 time-history, fits the Recht-Ipson
1963 perforation model, and reports V50 with a PASS / FAIL / INCONCLUSIVE
verdict against the Vakili Rad 2020 published value with a 7 percent
tolerance.

This script is intentionally self-contained except for two external dependencies
that must be installed in the virtualenv per master_plan.md section 7:

    pip install numpy scipy gmsh meshio jinja2

and the Lima + Apptainer + OpenRadioss + Tools/inp2rad + ap_ply_model_creation
toolchain per master_plan.md section 7. The runner detects whether each piece
is present and falls back to a dry-run mode if not, so the script is safe to
import or invoke for spec / inspection purposes on a workstation without the
full toolchain.

Author: J.C. Vaught
Date: 2026-04-29
Stage: 16 (terminal goal)
Master plan reference: plan/master_plan.md section 3 row 16
Spec: tests/stage_16_PW_panel_ballistic/spec.md
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

# ----------------------------------------------------------------------------
# Module-level constants pulled from spec.md sections 2 through 7.
# Single source of truth for material data is material_card_im7_8552.json
# (loaded at runtime); this module only carries the geometry / velocity-sweep
# constants that drive the bracketing scheme.
# ----------------------------------------------------------------------------

V50_TOLERANCE_PERCENT = 7.0
V50_TOLERANCE_PERCENT_BC_DEGRADED = 9.0
HOURGLASS_RATIO_GATE = 0.05
SIMULATION_DURATION_S = 0.5e-3
TIME_HISTORY_DT_S = 1.0e-6
ANIMATION_DT_S = 5.0e-6
DEFAULT_TIMESTEP_S = 4.4e-8

# Velocity-sweep multipliers about V50_published (per spec section 9 step 4).
V_SWEEP_MULTIPLIERS: Tuple[float, ...] = (0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20)

# BC-insensitivity tolerances (spec section 6).
BC_PEAK_FORCE_TOL_PERCENT = 5.0
BC_TIME_OF_PEAK_TOL_PERCENT = 10.0
BC_DURATION_TOL_PERCENT = 10.0

# CFL parameters for the wall-clock estimator (spec section 3.2).
CFRP_FIBER_DIRECTION_WAVE_SPEED_M_S = 10250.0  # sqrt(E1 / rho), IM7/8552, fiber axis.
CFL_SAFETY_FACTOR = 0.9

# Configurations the runner knows how to drive.
KNOWN_PANELS: Tuple[str, ...] = ("P1_TWT", "P2_WTW", "P3_Control")

# Default file locations relative to the runner's directory.
RUNNER_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = RUNNER_DIR / "kok_config.json"
DEFAULT_VAKILI_RAD_TABLE = RUNNER_DIR / "vakili_rad_2020_v50_table.json"
DEFAULT_MATERIAL_CARD = RUNNER_DIR / "material_card_im7_8552.json"
DEFAULT_BUILD_CACHE = RUNNER_DIR / "build"
DEFAULT_OUTPUT_DIR = RUNNER_DIR / "runs" / "stage_16"
DEFAULT_FIGURES_DIR = RUNNER_DIR / "figures"

# OpenRadioss / Lima invocation patterns (spec section 2.4).
KOK_PREPROCESSOR_CMD = (
    "python ap_ply_model_creation/ap_ply_model.py"
    " --config {config} --output-inp {inp_path}"
)
INP2RAD_CMD = (
    "python OpenRadioss/Tools/inp2rad/inp2rad.py {inp_path} --output-rad {rad_path}"
)
STARTER_CMD = "starter_linuxa64 -i {starter_deck} -nt {cores}"
ENGINE_CMD = "engine_linuxa64 -i {engine_deck} -nt {cores}"

# Logging.
log = logging.getLogger("stage_16")


# ----------------------------------------------------------------------------
# Dataclasses for results.
# ----------------------------------------------------------------------------


@dataclass
class ShotResult:
    """Single-shot output, one row of v50_sweep.csv."""

    v_test_m_per_s: float
    v_residual_m_per_s: float
    f_peak_n: float
    delta_back_face_mm: float
    n_penetrated_plies: int
    hourglass_energy_j: float
    internal_energy_j: float
    hourglass_ratio: float
    contact_duration_us: float
    bc_option: str
    status: str  # "OK" / "FAIL_HOURGLASS" / "FAIL_TIMEOUT" / "DRY_RUN"


@dataclass
class V50Summary:
    """Final per-panel summary, written to v50_summary.json."""

    panel_configuration: str
    v50_simulated_m_per_s: Optional[float] = None
    v50_uncertainty_m_per_s: Optional[float] = None
    v50_published_m_per_s: Optional[float] = None
    absolute_percent_error: Optional[float] = None
    tolerance_percent: float = V50_TOLERANCE_PERCENT
    bc_option: str = "A"
    bc_insensitivity_passed: Optional[bool] = None
    hourglass_max_ratio: Optional[float] = None
    hourglass_check_passed: Optional[bool] = None
    ranking_check_passed: Optional[bool] = None
    verdict: str = "INCONCLUSIVE"
    shots: List[ShotResult] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


# ----------------------------------------------------------------------------
# Configuration loaders.
# ----------------------------------------------------------------------------


def load_json(path: Path) -> Dict[str, Any]:
    """Load a JSON file. Returns empty dict if the file does not exist (the
    runner emits an actionable warning rather than crashing during dry-run /
    spec-inspection use)."""
    if not path.exists():
        log.warning("config file not found: %s", path)
        return {}
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: Path, data: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=4, sort_keys=False, default=_jsonable)


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "__dataclass_fields__"):
        return asdict(obj)
    if isinstance(obj, (set, tuple)):
        return list(obj)
    raise TypeError(f"Type {type(obj).__name__} is not JSON serializable")


def load_v50_published(table_path: Path, panel: str) -> Tuple[Optional[float], Optional[float], List[str]]:
    """Pull V50_published for the given panel from the Vakili Rad 2020 JSON
    table. Returns (v50, uncertainty, notes). If the published number is the
    string '[UNVERIFIED]', returns None for v50 and adds a note explaining the
    INCONCLUSIVE downstream behaviour."""
    table = load_json(table_path)
    notes: List[str] = []
    panels = table.get("panels", {})
    rec = panels.get(panel, {})
    v50 = rec.get("v50_m_per_s")
    unc = rec.get("v50_uncertainty_m_per_s")
    if isinstance(v50, str) and "UNVERIFIED" in v50.upper():
        notes.append(
            f"V50_published for {panel} is [UNVERIFIED] in {table_path.name};"
            " the runner will produce a simulated V50 but cannot evaluate the"
            " pass criterion. Re-extract the paywalled Vakili Rad 2020 *Composites"
            " Part B* paper Table 5 and update the JSON."
        )
        v50 = None
    if isinstance(unc, str) and "UNVERIFIED" in unc.upper():
        unc = None
    return (
        float(v50) if isinstance(v50, (int, float)) else None,
        float(unc) if isinstance(unc, (int, float)) else None,
        notes,
    )


# ----------------------------------------------------------------------------
# Geometry build.
# ----------------------------------------------------------------------------


def build_geometry(
    config: Dict[str, Any],
    build_cache: Path,
    panel: str,
    lima_vm: Optional[str],
    dry_run: bool,
) -> Tuple[Path, List[str]]:
    """Run the Kok preprocessor then inp2rad. Caches the .rad on disk so that
    subsequent V50 sweep shots reuse the same geometry. Returns the path to
    the converted OpenRadioss .rad fragment, plus a list of notes."""
    notes: List[str] = []
    build_cache.mkdir(parents=True, exist_ok=True)
    inp_path = build_cache / f"{panel.lower()}_section.inp"
    rad_path = build_cache / f"{panel.lower()}_section.rad"

    if rad_path.exists():
        notes.append(f"geometry build cache hit: {rad_path.name}")
        log.info("geometry cache hit, skipping rebuild: %s", rad_path)
        _sanity_check_inp_vs_rad(inp_path, rad_path, notes)
        return rad_path, notes

    log.info("building Kok AP-PLY geometry for %s", panel)
    config_for_kok = dict(config)
    config_for_kok.setdefault("ap_ply_layup", {})["configuration_name"] = panel

    kok_config_path = build_cache / f"kok_config_{panel.lower()}.json"
    write_json(kok_config_path, config_for_kok)

    kok_cmd = KOK_PREPROCESSOR_CMD.format(
        config=str(kok_config_path), inp_path=str(inp_path)
    )
    inp2rad_cmd = INP2RAD_CMD.format(
        inp_path=str(inp_path), rad_path=str(rad_path)
    )

    _run_in_lima(kok_cmd, lima_vm=lima_vm, dry_run=dry_run, notes=notes)
    _run_in_lima(inp2rad_cmd, lima_vm=lima_vm, dry_run=dry_run, notes=notes)

    if not dry_run:
        _sanity_check_inp_vs_rad(inp_path, rad_path, notes)
    else:
        notes.append("DRY-RUN: skipped Kok preprocessor + inp2rad invocation.")
    return rad_path, notes


def _sanity_check_inp_vs_rad(inp_path: Path, rad_path: Path, notes: List[str]) -> None:
    """Spec section 2.4 + risk 3: inp2rad is beta status; verify element and node
    count consistency. We count *NODE and *ELEMENT keywords in the .inp and the
    /NODE and /BRICK or /TETRA blocks in the .rad. Bookkeeping mismatch is a
    hard fail (runner aborts) because a silent /SKEW orientation drop would
    corrupt every per-element fiber direction in the laminate."""
    if not inp_path.exists() or not rad_path.exists():
        notes.append(
            f"sanity check skipped: {inp_path.name} or {rad_path.name} missing"
        )
        return
    n_nodes_inp, n_elem_inp = _count_inp_entities(inp_path)
    n_nodes_rad, n_elem_rad = _count_rad_entities(rad_path)
    notes.append(
        f"inp2rad sanity: nodes inp={n_nodes_inp} rad={n_nodes_rad}; "
        f"elem inp={n_elem_inp} rad={n_elem_rad}"
    )
    if n_nodes_inp and n_nodes_rad and n_nodes_inp != n_nodes_rad:
        raise RuntimeError(
            f"inp2rad node-count mismatch: inp {n_nodes_inp} vs rad {n_nodes_rad};"
            " converter is beta-status (audit row C); fix or fall back to manual"
            " conversion before continuing."
        )
    if n_elem_inp and n_elem_rad and n_elem_inp != n_elem_rad:
        raise RuntimeError(
            f"inp2rad element-count mismatch: inp {n_elem_inp} vs rad {n_elem_rad};"
            " converter dropped elements; abort before V50 sweep."
        )


def _count_inp_entities(path: Path) -> Tuple[int, int]:
    n_nodes = n_elem = 0
    section = None
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("*NODE"):
                section = "node"
                continue
            if stripped.startswith("*ELEMENT"):
                section = "elem"
                continue
            if stripped.startswith("*"):
                section = None
                continue
            if section == "node":
                n_nodes += 1
            elif section == "elem":
                n_elem += 1
    return n_nodes, n_elem


def _count_rad_entities(path: Path) -> Tuple[int, int]:
    n_nodes = n_elem = 0
    section = None
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("/NODE"):
                section = "node"
                continue
            if stripped.startswith(("/BRICK", "/HEXA", "/TETRA")):
                section = "elem"
                continue
            if stripped.startswith("/"):
                section = None
                continue
            if section == "node":
                n_nodes += 1
            elif section == "elem":
                n_elem += 1
    return n_nodes, n_elem


# ----------------------------------------------------------------------------
# Deck templating.
# ----------------------------------------------------------------------------


STARTER_TEMPLATE = """\
#-RADIOSS STARTER : Stage 16 P{panel} ballistic V50 sweep
#-Vakili Rad 2020 *Composites Part B* 203, 108478
/BEGIN
Stage 16 PW panel ballistic | {panel} | v_impact = {v_impact_m_per_s:.4f} m/s | BC = {bc_option}
2025         0
                  Mg              mm               s
                  Mg              mm               s
/UNIT/1
Mg mm s
#---
/INCLUDE
{geometry_rad_fragment}
#---
# LAW25 IM7G/8552-1 orthotropic, Soden 1998 + WWFE-II strengths
/MAT/LAW25/1
IM7G_8552_1_LAW25_orthotropic
{rho_kg_m3:.6e} {E1_Pa:.6e} {E2_Pa:.6e} {nu12:.6e} {G12_Pa:.6e}
{Xt_Pa:.6e} {Xc_Pa:.6e} {Yt_Pa:.6e} {Yc_Pa:.6e} {Sl_Pa:.6e}
#---
# /FAIL/HASHIN with Ifail = 2 (element erosion)
/FAIL/HASHIN/1
1 2 1.0 0.0 0.0 1.0e-4
#---
# LAW1 neat 8552 epoxy isotropic, resin pockets
/MAT/LAW1/2
neat_8552_isotropic
{rho_resin_kg_m3:.6e} {E_resin_Pa:.6e} {nu_resin:.6e}
#---
# Rigid steel projectile (LAW2 minimal) and /RBODY
/MAT/LAW2/3
steel_projectile_rigid
{rho_proj_kg_m3:.6e} {E_proj_Pa:.6e} {nu_proj:.6e}
#---
/PROP/TYPE14/1
TYPE14_HEXA8_HEPH_for_LAW25
0 24 0 0 0 0
#---
/PROP/TYPE14/2
TYPE14_HEXA8_HEPH_for_LAW1
0 24 0 0 0 0
#---
/INTER/TYPE2/1
inter_laminar_cohesive_ply_to_ply
{tau_n_max_Pa:.6e} {tau_s_max_Pa:.6e} {GIc_J_m2:.6e} {GIIc_J_m2:.6e} {eta_BK:.6e} 4
#---
/INTER/TYPE7/1
projectile_to_laminate_general_contact
0.30 0 0 0
#---
/RBODY/1
projectile_rbody
{m_proj_kg:.6e}
#---
/IMPVEL/1
projectile_initial_velocity
1 0.0 0.0 -{v_impact_m_per_s:.6e} 0.0 0.0 0.0
#---
/BCS/1
clamped_cut_edges {bc_option}
{bc_card_block}
#---
/TH/PART/1
projectile reaction_force back_face_node
/TH/NODE/1
projectile_centroid back_face_centre
#---
/H3D/ELEM/DAMA
/H3D/ELEM/EROS
/ANIM/ELEM/DAMA
/ANIM/ELEM/EROS
/ANIM/DT
{anim_dt_s:.6e}
#---
/END
"""

ENGINE_TEMPLATE = """\
#-RADIOSS ENGINE : Stage 16 P{panel} V_imp={v_impact_m_per_s:.4f} m/s, BC={bc_option}
/RUN/STAGE16/1
{panel}_v_imp_{v_impact_int}/
/STOP/ETIME
{simulation_duration_s:.6e}
/DT/BRICK/CST/1
0.9 {target_dt_s:.6e}
/DT/INTER/DEL/1
0.9 {target_dt_s:.6e}
/TFILE
{th_dt_s:.6e}
/PRINT/-1
/END
"""


def template_decks(
    output_dir: Path,
    panel: str,
    v_impact_m_per_s: float,
    bc_option: str,
    geometry_rad_fragment: Path,
    material_card: Dict[str, Any],
    cohesive_card: Dict[str, Any],
    projectile_card: Dict[str, Any],
) -> Tuple[Path, Path]:
    """Render the starter and engine decks. Returns (starter_path, engine_path)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    v_int = int(round(v_impact_m_per_s))
    base = f"{panel.lower()}_v{v_int:04d}_bc{bc_option}"
    starter_path = output_dir / f"{base}_0000.rad"
    engine_path = output_dir / f"{base}_0001.rad"

    bc_card_block = _bc_card_block(bc_option)

    starter_text = STARTER_TEMPLATE.format(
        panel=panel,
        v_impact_m_per_s=v_impact_m_per_s,
        bc_option=bc_option,
        geometry_rad_fragment=str(geometry_rad_fragment),
        rho_kg_m3=material_card.get("rho_kg_m3", 1570.0),
        E1_Pa=material_card.get("E1_Pa", 165e9),
        E2_Pa=material_card.get("E2_Pa", 8.4e9),
        nu12=material_card.get("nu12", 0.34),
        G12_Pa=material_card.get("G12_Pa", 5.6e9),
        Xt_Pa=material_card.get("Xt_Pa", 2560e6),
        Xc_Pa=material_card.get("Xc_Pa", 1590e6),
        Yt_Pa=material_card.get("Yt_Pa", 73e6),
        Yc_Pa=material_card.get("Yc_Pa", 185e6),
        Sl_Pa=material_card.get("Sl_Pa", 90e6),
        rho_resin_kg_m3=material_card.get("rho_resin_kg_m3", 1300.0),
        E_resin_Pa=material_card.get("E_resin_Pa", 4.67e9),
        nu_resin=material_card.get("nu_resin", 0.36),
        rho_proj_kg_m3=projectile_card.get("rho_kg_m3", 7850.0),
        E_proj_Pa=projectile_card.get("E_Pa", 200e9),
        nu_proj=projectile_card.get("nu", 0.30),
        tau_n_max_Pa=cohesive_card.get("tau_n_max_Pa", 60e6),
        tau_s_max_Pa=cohesive_card.get("tau_s_max_Pa", 90e6),
        GIc_J_m2=cohesive_card.get("GIc_J_m2", 277.0),
        GIIc_J_m2=cohesive_card.get("GIIc_J_m2", 788.0),
        eta_BK=cohesive_card.get("eta_BK", 1.45),
        m_proj_kg=projectile_card.get("m_kg", 13.4e-3),
        bc_card_block=bc_card_block,
        anim_dt_s=ANIMATION_DT_S,
    )
    engine_text = ENGINE_TEMPLATE.format(
        panel=panel,
        v_impact_m_per_s=v_impact_m_per_s,
        bc_option=bc_option,
        v_impact_int=v_int,
        simulation_duration_s=SIMULATION_DURATION_S,
        target_dt_s=DEFAULT_TIMESTEP_S,
        th_dt_s=TIME_HISTORY_DT_S,
    )

    starter_path.write_text(starter_text)
    engine_path.write_text(engine_text)
    return starter_path, engine_path


def _bc_card_block(bc_option: str) -> str:
    """Spec section 4.1: Option A clamped, Option B absorbing."""
    if bc_option.upper() == "A":
        return (
            "1 1 1 1 1 1 cut_edge_north\n"
            "1 1 1 1 1 1 cut_edge_south\n"
            "1 1 1 1 1 1 cut_edge_east\n"
            "1 1 1 1 1 1 cut_edge_west"
        )
    return (
        "/DAMP/INTER/1 absorbing_north 1.07e4\n"
        "/DAMP/INTER/2 absorbing_south 1.07e4\n"
        "/DAMP/INTER/3 absorbing_east  1.07e4\n"
        "/DAMP/INTER/4 absorbing_west  1.07e4"
    )


# ----------------------------------------------------------------------------
# Solver invocation (Lima / Apptainer wrapping).
# ----------------------------------------------------------------------------


def _run_in_lima(
    cmd: str,
    lima_vm: Optional[str],
    dry_run: bool,
    notes: List[str],
    timeout_s: Optional[float] = None,
) -> Tuple[int, str, str]:
    """Invoke a command inside the Lima Apptainer VM if `lima_vm` is set, else
    execute on the host."""
    if lima_vm:
        wrapped = f"limactl shell {shlex.quote(lima_vm)} -- bash -lc {shlex.quote(cmd)}"
    else:
        wrapped = cmd

    log.info("invoke: %s", wrapped)
    if dry_run:
        notes.append(f"DRY-RUN: would invoke `{wrapped}`")
        return 0, "", ""

    try:
        completed = subprocess.run(
            wrapped,
            shell=True,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )
    except subprocess.TimeoutExpired:
        notes.append(f"TIMEOUT: `{wrapped}` after {timeout_s} s")
        return 124, "", "TIMEOUT"
    if completed.returncode != 0:
        notes.append(
            f"non-zero exit ({completed.returncode}) on `{wrapped}`:\n"
            f"{completed.stderr[-2000:]}"
        )
    return completed.returncode, completed.stdout, completed.stderr


def run_starter_engine(
    starter_deck: Path,
    engine_deck: Path,
    cores: int,
    lima_vm: Optional[str],
    dry_run: bool,
    notes: List[str],
) -> Tuple[int, Path]:
    """Run starter then engine. Returns (return_code, working_dir)."""
    work_dir = starter_deck.parent
    starter_cmd = STARTER_CMD.format(
        starter_deck=starter_deck.name, cores=cores
    )
    engine_cmd = ENGINE_CMD.format(
        engine_deck=engine_deck.name, cores=cores
    )

    rc, _, _ = _run_in_lima(
        f"cd {shlex.quote(str(work_dir))} && {starter_cmd}",
        lima_vm=lima_vm, dry_run=dry_run, notes=notes,
    )
    if rc != 0 and not dry_run:
        return rc, work_dir
    rc, _, _ = _run_in_lima(
        f"cd {shlex.quote(str(work_dir))} && {engine_cmd}",
        lima_vm=lima_vm, dry_run=dry_run, notes=notes,
    )
    return rc, work_dir


# ----------------------------------------------------------------------------
# T01 time-history parsing and per-shot diagnostic extraction.
# ----------------------------------------------------------------------------


def parse_time_history(work_dir: Path, panel: str, v_int: int, bc_option: str) -> Dict[str, Any]:
    """Pull the projectile centroid x-velocity vs time, the contact force, the
    back-face displacement, and the energy bookkeeping out of the OpenRadioss
    T01 binary. The actual decode requires `vortex-radioss`; we wrap it here
    so the runner gracefully degrades if the package is absent."""
    th_path = work_dir / f"{panel.lower()}_v{v_int:04d}_bc{bc_option}_T01"
    result: Dict[str, Any] = {
        "v_residual_m_per_s": 0.0,
        "f_peak_n": 0.0,
        "delta_back_face_mm": 0.0,
        "n_penetrated_plies": 0,
        "hourglass_energy_j": 0.0,
        "internal_energy_j": 0.0,
        "contact_duration_us": 0.0,
        "status": "DRY_RUN",
    }
    if not th_path.exists():
        log.warning("T01 not found: %s; returning dry-run defaults", th_path)
        return result

    try:
        from vortex_radioss import VortexT01  # type: ignore
    except ImportError:
        log.warning(
            "vortex-radioss not installed; cannot decode %s. install with"
            " `pip install vortex-radioss` or clone github.com/Vortex-CAE/Vortex-Radioss",
            th_path,
        )
        result["status"] = "PARSE_FAILED_NO_VORTEX"
        return result

    try:
        th = VortexT01(str(th_path))  # type: ignore
        t = th.time
        v_proj = th.node_velocity("projectile_centroid", axis="z")
        f_contact = th.contact_force("projectile_to_laminate")
        d_back = th.node_displacement("back_face_centre", axis="z")
        e_int = th.global_quantity("internal_energy")
        e_hg = th.global_quantity("hourglass_energy")
        n_eros = th.global_quantity("n_eroded_elements")
    except Exception as exc:  # pragma: no cover, vortex-radioss internals
        log.error("T01 decode failed: %s", exc)
        result["status"] = "PARSE_FAILED"
        return result

    v_residual = float(v_proj[-1]) if len(v_proj) else 0.0
    if v_residual < 0.0:
        v_residual = 0.0  # rebound; clamp to zero per Recht-Ipson convention
    f_peak = float(max(abs(x) for x in f_contact)) if len(f_contact) else 0.0
    delta_max = float(max(abs(x) for x in d_back)) if len(d_back) else 0.0

    # Contact duration: full-width at 10% of peak.
    contact_duration_us = _contact_duration_us(t, f_contact, threshold=0.10 * f_peak)

    e_int_min = float(min(x for x in e_int if x > 0)) if any(x > 0 for x in e_int) else 1.0
    hg_max = float(max(e_hg)) if len(e_hg) else 0.0

    n_penetrated = _ply_count_from_erosion(n_eros)

    result.update(
        v_residual_m_per_s=v_residual,
        f_peak_n=f_peak,
        delta_back_face_mm=delta_max * 1000.0,  # m -> mm
        n_penetrated_plies=n_penetrated,
        hourglass_energy_j=hg_max,
        internal_energy_j=e_int_min,
        contact_duration_us=contact_duration_us,
        status="OK",
    )
    return result


def _contact_duration_us(t: Sequence[float], f: Sequence[float], threshold: float) -> float:
    if not len(t) or not len(f) or threshold <= 0:
        return 0.0
    above = [(ti, fi) for ti, fi in zip(t, f) if abs(fi) >= threshold]
    if not above:
        return 0.0
    return (above[-1][0] - above[0][0]) * 1.0e6  # s -> us


def _ply_count_from_erosion(n_eros_history: Sequence[float]) -> int:
    """Heuristic: penetrated plies are estimated as the saturation level of
    eroded elements divided by the per-ply impact-zone element count. The
    runner's geometry dump should provide the per-ply count; pending that
    plumbing this returns a zero-floor estimate."""
    if not len(n_eros_history):
        return 0
    n_max = int(max(n_eros_history))
    elements_per_ply_estimate = 20000  # spec section 3.1 impact-zone divided by 24 plies
    return min(24, max(0, n_max // elements_per_ply_estimate))


# ----------------------------------------------------------------------------
# Recht-Ipson 1963 V50 fit.
# ----------------------------------------------------------------------------


def recht_ipson_fit(
    v_test: Sequence[float], v_residual: Sequence[float]
) -> Tuple[float, float, float]:
    """Fit V_r^2 = a (V_i^2 - V_50^2) for V_i >= V_50 (over-match shots only).
    Returns (V_50_simulated, V_50_uncertainty, a)."""
    try:
        import numpy as np
        from scipy.optimize import curve_fit
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "scipy and numpy required for Recht-Ipson fit;"
            " pip install scipy numpy"
        ) from exc

    v_test = np.asarray(v_test, dtype=float)
    v_res = np.asarray(v_residual, dtype=float)
    over = v_res > 0.0
    if over.sum() < 2:
        # Not enough over-match shots to fit two parameters; fall back to the
        # bracket midpoint between the highest-velocity sub-V50 shot and the
        # lowest-velocity over-match shot.
        sub = ~over
        if over.any() and sub.any():
            v50_est = 0.5 * (v_test[sub].max() + v_test[over].min())
            return float(v50_est), float(v50_est * 0.05), 1.0
        return float(v_test.mean()), float(v_test.std() or v_test.mean() * 0.05), 1.0

    def _model(vi, a, v50):
        return np.sqrt(np.clip(a * (vi ** 2 - v50 ** 2), 0.0, None))

    p0 = [1.0, float(v_test[over].min())]
    bounds = ([0.1, 0.0], [10.0, float(v_test.max())])
    popt, pcov = curve_fit(_model, v_test[over], v_res[over], p0=p0, bounds=bounds)
    a_fit, v50_fit = popt
    sigma_v50 = float(np.sqrt(pcov[1, 1])) if pcov is not None else float("nan")
    return float(v50_fit), sigma_v50, float(a_fit)


# ----------------------------------------------------------------------------
# Hourglass-ratio gate.
# ----------------------------------------------------------------------------


def hourglass_check(shots: Sequence[ShotResult]) -> Tuple[float, bool]:
    """Spec section 3.3: E_hg / E_int <= 0.05."""
    ratios = [
        s.hourglass_energy_j / s.internal_energy_j
        for s in shots
        if s.internal_energy_j > 0
    ]
    if not ratios:
        return 0.0, False
    max_ratio = max(ratios)
    return max_ratio, bool(max_ratio <= HOURGLASS_RATIO_GATE)


# ----------------------------------------------------------------------------
# BC-insensitivity gate (spec section 6).
# ----------------------------------------------------------------------------


def bc_insensitivity_check(
    shot_a: ShotResult, shot_b: ShotResult, notes: List[str]
) -> bool:
    if shot_a.f_peak_n <= 0 or shot_b.f_peak_n <= 0:
        notes.append("BC insensitivity: zero peak force on one of the shots")
        return False
    delta_f = abs(shot_a.f_peak_n - shot_b.f_peak_n) / max(shot_a.f_peak_n, shot_b.f_peak_n)
    delta_t = (
        abs(shot_a.contact_duration_us - shot_b.contact_duration_us)
        / max(shot_a.contact_duration_us, shot_b.contact_duration_us, 1e-9)
    )
    notes.append(
        f"BC insensitivity: dF_peak = {delta_f * 100.0:.2f}%,"
        f" dT_duration = {delta_t * 100.0:.2f}%"
    )
    return (
        delta_f * 100.0 <= BC_PEAK_FORCE_TOL_PERCENT
        and delta_t * 100.0 <= BC_DURATION_TOL_PERCENT
    )


# ----------------------------------------------------------------------------
# V50 sweep driver.
# ----------------------------------------------------------------------------


def run_v50_sweep(
    panel: str,
    v50_published: Optional[float],
    v50_publication_uncertainty: Optional[float],
    config: Dict[str, Any],
    material_card: Dict[str, Any],
    cohesive_card: Dict[str, Any],
    projectile_card: Dict[str, Any],
    geometry_rad: Path,
    bc_option: str,
    output_dir: Path,
    cores: int,
    lima_vm: Optional[str],
    dry_run: bool,
) -> List[ShotResult]:
    """Loop over V_SWEEP_MULTIPLIERS, template decks, run starter + engine,
    parse T01, and assemble per-shot ShotResult records."""
    if v50_published is None:
        log.warning(
            "V50_published is [UNVERIFIED]; sweep centred on placeholder 110 m/s."
        )
        v50_centre = 110.0
    else:
        v50_centre = float(v50_published)

    velocities = [v50_centre * m for m in V_SWEEP_MULTIPLIERS]

    shots: List[ShotResult] = []
    for v in velocities:
        v_int = int(round(v))
        starter_deck, engine_deck = template_decks(
            output_dir=output_dir / f"{panel}_v{v_int:04d}_bc{bc_option}",
            panel=panel,
            v_impact_m_per_s=v,
            bc_option=bc_option,
            geometry_rad_fragment=geometry_rad,
            material_card=material_card,
            cohesive_card=cohesive_card,
            projectile_card=projectile_card,
        )
        notes: List[str] = []
        rc, work_dir = run_starter_engine(
            starter_deck=starter_deck,
            engine_deck=engine_deck,
            cores=cores,
            lima_vm=lima_vm,
            dry_run=dry_run,
            notes=notes,
        )
        th = parse_time_history(work_dir, panel, v_int, bc_option)
        shot_status = th["status"]
        if rc != 0 and not dry_run:
            shot_status = "FAIL_SOLVER"
        hg_ratio = (
            th["hourglass_energy_j"] / th["internal_energy_j"]
            if th["internal_energy_j"] > 0 else 0.0
        )
        shots.append(
            ShotResult(
                v_test_m_per_s=v,
                v_residual_m_per_s=float(th["v_residual_m_per_s"]),
                f_peak_n=float(th["f_peak_n"]),
                delta_back_face_mm=float(th["delta_back_face_mm"]),
                n_penetrated_plies=int(th["n_penetrated_plies"]),
                hourglass_energy_j=float(th["hourglass_energy_j"]),
                internal_energy_j=float(th["internal_energy_j"]),
                hourglass_ratio=hg_ratio,
                contact_duration_us=float(th["contact_duration_us"]),
                bc_option=bc_option,
                status=shot_status,
            )
        )
    return shots


# ----------------------------------------------------------------------------
# CSV emission for downstream Typst / CeTZ figures.
# ----------------------------------------------------------------------------


def write_v50_sweep_csv(path: Path, shots: Sequence[ShotResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "v_test_m_per_s", "v_residual_m_per_s", "F_peak_N",
                "delta_back_face_mm", "n_penetrated_plies",
                "hourglass_energy_J", "internal_energy_J", "hourglass_ratio",
                "contact_duration_us", "BC_option", "status",
            ]
        )
        for s in shots:
            writer.writerow(
                [
                    f"{s.v_test_m_per_s:.4f}",
                    f"{s.v_residual_m_per_s:.4f}",
                    f"{s.f_peak_n:.4e}",
                    f"{s.delta_back_face_mm:.4f}",
                    s.n_penetrated_plies,
                    f"{s.hourglass_energy_j:.4e}",
                    f"{s.internal_energy_j:.4e}",
                    f"{s.hourglass_ratio:.4e}",
                    f"{s.contact_duration_us:.2f}",
                    s.bc_option,
                    s.status,
                ]
            )


def write_recht_ipson_csv(
    path: Path, shots: Sequence[ShotResult], v50: float, a_fit: float
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["v_i_m_per_s", "v_r_simulated_m_per_s", "v_r_recht_ipson_m_per_s"])
        for s in shots:
            v_r_model = (
                math.sqrt(max(a_fit * (s.v_test_m_per_s ** 2 - v50 ** 2), 0.0))
            )
            writer.writerow(
                [f"{s.v_test_m_per_s:.4f}",
                 f"{s.v_residual_m_per_s:.4f}",
                 f"{v_r_model:.4f}"]
            )


# ----------------------------------------------------------------------------
# Entry point.
# ----------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Stage 16 runner: V50 sweep on a UofSC pseudo-woven AP-PLY panel."
            " Project terminal goal."
        )
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--vakili-rad-table", type=Path, default=DEFAULT_VAKILI_RAD_TABLE,
    )
    parser.add_argument(
        "--material-card", type=Path, default=DEFAULT_MATERIAL_CARD,
    )
    parser.add_argument("--panel", choices=KNOWN_PANELS, default="P1_TWT")
    parser.add_argument("--all-configs", action="store_true",
                        help="run P1_TWT, P2_WTW, P3_Control sequentially.")
    parser.add_argument("--cores", type=int, default=16)
    parser.add_argument("--mpi-hosts", type=str, default="")
    parser.add_argument("--mass-scaling", type=float, default=1.0,
                        help="elemental mass-scaling factor; >1 trades inertia"
                             " overshoot for shorter wall clock")
    parser.add_argument("--cohesive-mode", choices=("type2", "law83"),
                        default="type2")
    parser.add_argument("--bc-option", choices=("A", "B", "both"), default="A")
    parser.add_argument("--lima-vm", type=str, default="apptainer",
                        help="Lima VM name; pass empty string to run on host.")
    parser.add_argument("--build-cache", type=Path, default=DEFAULT_BUILD_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figures", type=Path, default=DEFAULT_FIGURES_DIR)
    parser.add_argument("--dry-run", action="store_true",
                        help="template decks but do not invoke OpenRadioss.")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)-7s %(name)s | %(message)s",
        datefmt="%H:%M:%S",
    )

    config = load_json(args.config)
    material_card = load_json(args.material_card).get("LAW25", {})
    cohesive_card = load_json(args.material_card).get("INTER_TYPE2", {})
    projectile_card = load_json(args.material_card).get("PROJECTILE", {})
    if args.mass_scaling > 1.0:
        log.warning(
            "mass scaling = %.2f enabled; expect <5%% inertial overshoot in"
            " contact force per Stage 14 audit.",
            args.mass_scaling,
        )

    panels_to_run: List[str] = (
        list(KNOWN_PANELS) if args.all_configs else [args.panel]
    )
    bc_options: List[str] = (
        ["A", "B"] if args.bc_option == "both" else [args.bc_option]
    )

    overall_results: Dict[str, V50Summary] = {}

    for panel in panels_to_run:
        log.info("=" * 60)
        log.info("Stage 16 | panel %s | BCs=%s", panel, bc_options)
        log.info("=" * 60)

        v50_pub, v50_unc, pub_notes = load_v50_published(args.vakili_rad_table, panel)
        summary = V50Summary(
            panel_configuration=panel,
            v50_published_m_per_s=v50_pub,
            bc_option="+".join(bc_options),
            notes=list(pub_notes),
        )

        try:
            geometry_rad, geom_notes = build_geometry(
                config=config,
                build_cache=args.build_cache,
                panel=panel,
                lima_vm=(args.lima_vm or None),
                dry_run=args.dry_run,
            )
            summary.notes.extend(geom_notes)
        except Exception as exc:  # surface upstream pipeline failures here
            summary.notes.append(f"geometry build failed: {exc}")
            summary.verdict = "FAIL_GEOMETRY"
            overall_results[panel] = summary
            continue

        # Run the BC-insensitivity check on the V50_centre shot only,
        # iff both options are requested OR if the user picked Option A
        # (we still want a Section 6 record for the run log).
        bc_check_passed: Optional[bool] = None
        if "A" in bc_options and "B" in bc_options:
            log.info("BC-insensitivity: running V50_centre shot under both BCs")
            bc_a_shots = run_v50_sweep(
                panel, v50_pub, v50_unc, config,
                material_card, cohesive_card, projectile_card,
                geometry_rad, "A", args.output, args.cores,
                args.lima_vm or None, args.dry_run,
            )
            bc_b_shots = run_v50_sweep(
                panel, v50_pub, v50_unc, config,
                material_card, cohesive_card, projectile_card,
                geometry_rad, "B", args.output, args.cores,
                args.lima_vm or None, args.dry_run,
            )
            # Use only the v_centre = V_SWEEP_MULTIPLIERS[3] = 1.00x shot for
            # the comparison.
            shot_a_centre = bc_a_shots[3] if len(bc_a_shots) > 3 else bc_a_shots[0]
            shot_b_centre = bc_b_shots[3] if len(bc_b_shots) > 3 else bc_b_shots[0]
            bc_check_passed = bc_insensitivity_check(
                shot_a_centre, shot_b_centre, summary.notes
            )
            shots = bc_a_shots + bc_b_shots
            summary.tolerance_percent = (
                V50_TOLERANCE_PERCENT_BC_DEGRADED if not bc_check_passed
                else V50_TOLERANCE_PERCENT
            )
            if not bc_check_passed:
                summary.notes.append(
                    "BC insensitivity FAIL; tolerance relaxed to"
                    f" {V50_TOLERANCE_PERCENT_BC_DEGRADED}%."
                )
        else:
            shots = []
            for bc in bc_options:
                shots.extend(run_v50_sweep(
                    panel, v50_pub, v50_unc, config,
                    material_card, cohesive_card, projectile_card,
                    geometry_rad, bc, args.output, args.cores,
                    args.lima_vm or None, args.dry_run,
                ))
            bc_check_passed = None
            summary.notes.append(
                "BC-insensitivity check skipped (single BC option requested)."
            )

        summary.shots = shots
        summary.bc_insensitivity_passed = bc_check_passed

        # Hourglass gate.
        hg_max, hg_ok = hourglass_check(shots)
        summary.hourglass_max_ratio = hg_max
        summary.hourglass_check_passed = hg_ok
        if not hg_ok:
            summary.notes.append(
                f"hourglass ratio max {hg_max:.4f} exceeds gate"
                f" {HOURGLASS_RATIO_GATE:.2f}; V50 reported but pass criterion"
                " is not met."
            )

        # Recht-Ipson fit + pass criterion.
        if any(s.status == "OK" for s in shots):
            v50_sim, v50_sigma, a_fit = recht_ipson_fit(
                [s.v_test_m_per_s for s in shots if s.status == "OK"],
                [s.v_residual_m_per_s for s in shots if s.status == "OK"],
            )
            summary.v50_simulated_m_per_s = v50_sim
            summary.v50_uncertainty_m_per_s = v50_sigma

            if v50_pub is not None:
                err = abs(v50_sim - v50_pub) / v50_pub * 100.0
                summary.absolute_percent_error = err
                pass_v50 = err <= summary.tolerance_percent
                if pass_v50 and (hg_ok if hg_ok is not None else True):
                    summary.verdict = "PASS"
                else:
                    summary.verdict = "FAIL"
            else:
                summary.verdict = "INCONCLUSIVE"
                summary.notes.append(
                    "no V50_published available; cannot evaluate pass criterion."
                )

            # CSV emission for Typst / CeTZ.
            args.figures.mkdir(parents=True, exist_ok=True)
            write_v50_sweep_csv(
                args.figures / f"v50_sweep_{panel}.csv", shots
            )
            write_recht_ipson_csv(
                args.figures / f"recht_ipson_fit_{panel}.csv",
                shots, v50_sim, a_fit,
            )
        else:
            summary.verdict = "FAIL_NO_SHOTS"

        # JSON summary.
        write_json(args.output / f"v50_summary_{panel}.json", asdict(summary))

        log.info("PANEL %s VERDICT: %s", panel, summary.verdict)
        if summary.v50_simulated_m_per_s is not None:
            log.info(
                "  V50_sim = %.2f +/- %.2f m/s | V50_pub = %s | err = %s%%",
                summary.v50_simulated_m_per_s,
                summary.v50_uncertainty_m_per_s or float("nan"),
                f"{v50_pub:.2f}" if v50_pub is not None else "[UNVERIFIED]",
                f"{summary.absolute_percent_error:.2f}"
                if summary.absolute_percent_error is not None else "n/a",
            )
        for note in summary.notes:
            log.info("  note: %s", note)

        overall_results[panel] = summary

    # Optional ranking gate when --all-configs.
    if args.all_configs and all(
        overall_results[p].v50_simulated_m_per_s is not None
        for p in KNOWN_PANELS
    ):
        ranking_sim = sorted(
            KNOWN_PANELS,
            key=lambda p: overall_results[p].v50_simulated_m_per_s or 0.0,
            reverse=True,
        )
        # Per Vakili Rad 2020 the published ranking is P1_TWT > P2_WTW > P3_Control.
        ranking_pub = ["P1_TWT", "P2_WTW", "P3_Control"]
        ranking_passed = ranking_sim == ranking_pub
        for p in KNOWN_PANELS:
            overall_results[p].ranking_check_passed = ranking_passed
        log.info(
            "ranking sim=%s pub=%s : %s",
            ranking_sim, ranking_pub, "PASS" if ranking_passed else "FAIL",
        )

    # Final aggregate JSON.
    write_json(
        args.output / "v50_summary_all.json",
        {p: asdict(s) for p, s in overall_results.items()},
    )

    overall_pass = all(s.verdict == "PASS" for s in overall_results.values())
    return 0 if overall_pass else 1


if __name__ == "__main__":
    sys.exit(main())
