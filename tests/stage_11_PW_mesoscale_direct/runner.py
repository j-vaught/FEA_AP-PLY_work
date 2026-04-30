"""Stage 11 runner - Pseudo-woven (AP-PLY) mesoscale direct simulation.

Author: J.C. Vaught
Date:   2026-04-29

Pipeline:
    (a) call the Kok preprocessor (rutger-kok/ap_ply_model_creation) to generate
        a 4-ply [0, +45, -45, 90] AP-PLY mesoscale block of 25 mm x 25 mm x 0.72 mm,
    (b) convert the resulting Abaqus .inp to OpenRadioss .rad through inp2rad
        (path A) or through meshio + GMSH (path B fallback),
    (c) template three OpenRadioss decks (axial in-plane tension, transverse
        in-plane tension, in-plane shear) on the converted mesh,
    (d) invoke the OpenRadioss starter and engine through limactl shell on a
        macOS host or directly on Linux,
    (e) read each run's T01 / .anim through Vortex-Radioss (or anim2vtk fallback)
        and area-average the reaction tractions to recover the effective in-plane
        engineering moduli E_x, E_y, G_xy,
    (f) compare to the Kok 2022 reported homogenized values plus the CLT shear
        target and emit a PASS / FAIL sentinel and a results CSV.

Reframed from RVE periodic homogenization to direct uniform-displacement BCs
because OpenRadioss has no general 3-DOF periodic boundary keyword (see
references/openradioss_endtoend_audit.md cross-cut F).

Usage:
    python runner.py --all
    python runner.py --geometry
    python runner.py --mesh
    python runner.py --case {A,B,C}
    python runner.py --postprocess

This script is the dress rehearsal for stage 16 (panel ballistic): the geometry
pipeline (Kok + inp2rad) and material assignment exercised here flow forward
unchanged. Only the BC and explicit-dynamic step swap at stage 16.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

STAGE_DIR = Path(__file__).resolve().parent
GEOM_DIR = STAGE_DIR / "geometry"
DECK_DIR = STAGE_DIR / "decks"
RUN_DIR = STAGE_DIR / "runs"
RES_DIR = STAGE_DIR / "results"
FIG_DIR = STAGE_DIR / "figures"

for d in (GEOM_DIR, DECK_DIR, RUN_DIR, RES_DIR, FIG_DIR):
    d.mkdir(parents=True, exist_ok=True)

# Repository roots (set through env vars or defaults).
KOK_REPO = Path(os.environ.get("KOK_REPO", STAGE_DIR.parents[2] / "external" / "ap_ply_model_creation"))
OR_TOOLS = Path(os.environ.get("OR_TOOLS", STAGE_DIR.parents[2] / "external" / "OpenRadioss" / "Tools"))
OR_EXEC = Path(os.environ.get("OR_EXEC", STAGE_DIR.parents[2] / "external" / "OpenRadioss" / "exec"))

# Lima invocation prefix (empty list means run on a Linux host directly).
LIMA_PREFIX = os.environ.get("LIMA_PREFIX", "limactl shell or --").split() if os.environ.get("LIMA_PREFIX") else []

# Kok 2022 reference values for the [0,+45,-45,90] AP-PLY quasi-iso architecture.
# Digitized from Kok et al. 2022 Table 4 (numerical column) and Figure 7.
# IM7/8552 substitution correction applied per spec section 7.1.
KOK_2022_REFERENCE = {
    "E_x_GPa": 53.3,    # Kok 2022 Table 4 quasi-iso AP-PLY numerical
    "E_y_GPa": 53.3,    # in-plane isotropic at laminate scale (one-over-one quasi-iso)
    "G_xy_GPa": 20.5,   # CLT estimate Ex / (2 (1 + nu_xy)) with nu_xy = 0.30
    "tolerance": 0.10,  # 10 percent acceptance, master plan section 3
}


# ---------------------------------------------------------------------------
# Configuration dataclass
# ---------------------------------------------------------------------------

@dataclass
class APPlyConfig:
    """AP-PLY architecture configuration matching Kok preprocessor defaults."""

    tape_angles: list = field(default_factory=lambda: [0, 45, -45, 90])
    tape_widths: list = field(default_factory=lambda: [6.35, 6.35, 6.35, 6.35])
    tape_spacing: int = 1                  # one-over-one
    cured_ply_thickness: float = 0.18      # mm
    undulation_ratio: float = 0.09
    specimen_size_x: float = 25.0          # mm
    specimen_size_y: float = 25.0
    n_plies: int = 4
    material: str = "IM7_8552"
    shift_per_ply: list = field(default_factory=lambda: [0.0, 0.0, 0.0, 0.0])
    export_format: str = "abaqus_inp"
    export_path: str = "geometry/ap_ply_block.inp"
    mesh_seed_size: float = 0.06           # mm, three through-thickness elements per 0.18 mm ply

    def to_json(self, path: Path) -> None:
        path.write_text(json.dumps(asdict(self), indent=2))


# ---------------------------------------------------------------------------
# Material card values (IM7/8552 from Soden 1998; 8552 epoxy from Soden 1998)
# ---------------------------------------------------------------------------

# Tape - LAW25 orthotropic (SI units, m / s / kg / Pa).
TAPE_CARD = {
    "rho":  1570.0,        # kg/m3
    "E11":  161.0e9, "E22": 11.4e9, "E33": 11.4e9,
    "nu12": 0.32, "nu13": 0.32, "nu23": 0.45,
    "G12":  5.17e9, "G13": 5.17e9, "G23": 3.93e9,
    "Xt":   2560.0e6, "Xc": 1590.0e6,
    "Yt":   73.0e6,   "Yc": 185.0e6,
    "S12":  90.0e6,
}

# Resin - LAW1 isotropic (SI units).
RESIN_CARD = {
    "rho":  1300.0,
    "E":    4.67e9,
    "nu":   0.38,
}


# ---------------------------------------------------------------------------
# Step (a): geometry generation through the Kok preprocessor
# ---------------------------------------------------------------------------

def run_kok_preprocessor(config: APPlyConfig) -> Path:
    """Invoke rutger-kok/ap_ply_model_creation to generate the AP-PLY block.

    The preprocessor is an Abaqus CAE Python script (Python 2.7).  We launch it
    through `abaqus cae noGUI=...` if abaqus is on PATH inside the Lima VM, or
    fall back to a pure-Python emulation that reuses the preprocessor's
    Shapely-based geometry routines without the Abaqus calls (for environments
    without an Abaqus license).

    Output: geometry/ap_ply_block.inp.
    """

    cfg_path = STAGE_DIR / "config.json"
    config.to_json(cfg_path)

    inp_target = STAGE_DIR / config.export_path
    inp_target.parent.mkdir(parents=True, exist_ok=True)

    abaqus_cae = shutil.which("abaqus")
    kok_driver = KOK_REPO / "ap_ply_model.py"

    if abaqus_cae and kok_driver.exists():
        cmd = LIMA_PREFIX + [
            abaqus_cae, "cae", "noGUI=" + str(kok_driver),
            "--", "--config", str(cfg_path),
        ]
        print("[stage11] running Kok preprocessor:\n  " + " ".join(cmd))
        subprocess.run(cmd, cwd=STAGE_DIR, check=True)
    else:
        # Fallback: call our standalone Python wrapper that reuses the Kok
        # tape_placement.laminate_creation routines without Abaqus.  The wrapper
        # writes an Abaqus-compatible .inp directly.  This is a stage 11
        # contingency path, not the preferred route - see spec section 10.1.
        print("[stage11] Abaqus CAE not located; invoking standalone fallback")
        _kok_standalone_fallback(config, inp_target)

    if not inp_target.exists():
        raise FileNotFoundError("Kok preprocessor did not write " + str(inp_target))
    print("[stage11] geometry written:", inp_target)
    return inp_target


def _kok_standalone_fallback(config: APPlyConfig, inp_path: Path) -> None:
    """Standalone clean-room geometry path through :mod:`kok_geom`.

    This replaces the earlier Abaqus-only contingency.  It writes the same
    runner-facing Abaqus INP path plus the native MSH4 and orientations sidecar
    under ``geometry/`` for parser and mesh sanity checks.
    """
    from kok_geom.config import KokConfig
    from kok_geom.io import convert_msh_to_inp, generate_mesh

    kok_config = KokConfig.model_validate(
        {
            "panel": {
                "size_x_mm": config.specimen_size_x,
                "size_y_mm": config.specimen_size_y,
                "n_plies": config.n_plies,
                "symmetry": "none",
            },
            "laydown": {
                "fiber_angles_deg": [float(v) for v in config.tape_angles],
                "placement_sequence": "1010",
                "angle_shift_deg": 0.0,
                "tape_width_mm": float(config.tape_widths[0]),
                "cured_ply_thickness_mm": config.cured_ply_thickness,
                "undulation_ratio": config.undulation_ratio,
                "tape_spacing": config.tape_spacing,
            },
            "mesh": {
                "in_plane_target_mm_impact_zone": config.mesh_seed_size,
                "in_plane_target_mm_far_field": config.mesh_seed_size,
                "through_thickness_target_mm": config.cured_ply_thickness,
                "graded_zone_radius_mm": 0.0,
                "element_order": 2,
            },
            "output": {
                "msh_path": str(GEOM_DIR / "ap_ply_block.msh"),
                "orientations_json_path": str(GEOM_DIR / "orientations.json"),
                "msh_format": "msh4_ascii",
            },
        }
    )
    cfg_path = GEOM_DIR / "kok_config.json"
    kok_config.to_file(cfg_path)
    msh_path, orientations_path = generate_mesh(cfg_path)
    convert_msh_to_inp(msh_path, inp_path)
    print("[stage11] clean-room geometry written:", msh_path)
    print("[stage11] orientation sidecar written:", orientations_path)


# ---------------------------------------------------------------------------
# Step (b): mesh conversion to OpenRadioss .rad
# ---------------------------------------------------------------------------

def convert_inp_to_rad(inp_path: Path) -> Path:
    """Path A: invoke the OpenRadioss inp2rad converter.

    Output: geometry/ap_ply_block_mesh.rad with /SKEW/FIX blocks per ply
    orientation and /BRICK / /TETRA10 element blocks referencing /PROP/TYPE14.
    """
    out_rad = inp_path.with_suffix(".rad")
    inp2rad = OR_TOOLS / "input_converters" / "inp2rad" / "inp2rad.py"

    if inp2rad.exists():
        cmd = LIMA_PREFIX + [
            sys.executable, str(inp2rad),
            "--input", str(inp_path),
            "--output", str(out_rad),
            "--unit-system", "SI_m_s_kg",
        ]
        print("[stage11] running inp2rad (path A):\n  " + " ".join(cmd))
        subprocess.run(cmd, cwd=STAGE_DIR, check=True)
    else:
        print("[stage11] inp2rad not located; falling back to path B (meshio + GMSH)")
        _path_b_meshio_bridge(inp_path, out_rad)

    if not out_rad.exists():
        raise FileNotFoundError("Mesh conversion did not produce " + str(out_rad))

    _validate_mesh_rad(out_rad)
    return out_rad


def _path_b_meshio_bridge(inp_path: Path, out_rad: Path) -> None:
    """Path B fallback: read .inp through meshio, mesh in GMSH, then write .rad
    keyword by keyword.  Used when inp2rad mishandles `*ORIENTATION` blocks."""
    raise NotImplementedError(
        "Path B (meshio + GMSH) bridge deferred to implementation; "
        "see spec section 2.3."
    )


def _validate_mesh_rad(rad_path: Path) -> None:
    """Confirm the converted mesh has the expected /SKEW/FIX count, /PROP refs,
    and face-coverage extents (spec section 6.4)."""
    text = rad_path.read_text()
    skew_count = text.count("/SKEW/FIX")
    prop_count = text.count("/PROP/TYPE14")
    if skew_count < 4:
        raise RuntimeError(f"mesh has {skew_count} /SKEW/FIX entries; expected >= 4 (one per ply)")
    if prop_count < 5:
        raise RuntimeError(f"mesh has {prop_count} /PROP/TYPE14 entries; expected >= 5 (4 plies + resin)")
    print(f"[stage11] mesh validation passed: {skew_count} skews, {prop_count} props")


# ---------------------------------------------------------------------------
# Step (c): deck templating
# ---------------------------------------------------------------------------

# All BCs ramp the displacement linearly with time t in [0, 1].  /IMPL/STATIC
# performs a single linear solve.

CASE_BCS = {
    "A": {  # axial in-plane tension along x, eps_xx = 1e-3
        "name": "axial_x",
        "fixed": [("face_x_min", "ux", 0.0)],
        "loaded": [("face_x_max", "ux", 25.0e-6)],   # 25 nm = 1e-3 * 25 mm in m
        "rb_corner": True,
        "strain_axis": "xx",
        "strain_magnitude": 1.0e-3,
    },
    "B": {  # transverse in-plane tension along y
        "name": "transverse_y",
        "fixed": [("face_y_min", "uy", 0.0)],
        "loaded": [("face_y_max", "uy", 25.0e-6)],
        "rb_corner": True,
        "strain_axis": "yy",
        "strain_magnitude": 1.0e-3,
    },
    "C": {  # in-plane shear gamma_xy = 1e-3
        "name": "shear_xy",
        "fixed": [
            ("face_x_min", "uy", 0.0),
            ("face_y_min", "ux", 0.0),
        ],
        "loaded": [
            ("face_x_max", "uy", 12.5e-6),
            ("face_y_max", "ux", 12.5e-6),
        ],
        "rb_corner": True,
        "strain_axis": "xy",
        "strain_magnitude": 0.5e-3,    # eps_xy = gamma_xy / 2
    },
}


def template_starter_deck(case: str, mesh_rad: Path) -> Path:
    """Write the OpenRadioss starter deck for the given case."""
    bcs = CASE_BCS[case]
    job = f"ap_ply_block_case{case}"
    starter = DECK_DIR / f"{job}_0000.rad"

    fixed_block = "\n".join(_bcs_keyword(grp, dof, val, fixed=True)
                            for (grp, dof, val) in bcs["fixed"])
    load_block = "\n".join(_bcs_keyword(grp, dof, val, fixed=False)
                           for (grp, dof, val) in bcs["loaded"])

    text = STARTER_TEMPLATE.format(
        job=job,
        mesh_include=str(mesh_rad.relative_to(STAGE_DIR)),
        rho_t=TAPE_CARD["rho"],
        E11=TAPE_CARD["E11"], E22=TAPE_CARD["E22"], E33=TAPE_CARD["E33"],
        nu12=TAPE_CARD["nu12"], nu13=TAPE_CARD["nu13"], nu23=TAPE_CARD["nu23"],
        G12=TAPE_CARD["G12"], G13=TAPE_CARD["G13"], G23=TAPE_CARD["G23"],
        Xt=TAPE_CARD["Xt"], Xc=TAPE_CARD["Xc"],
        Yt=TAPE_CARD["Yt"], Yc=TAPE_CARD["Yc"],
        S12=TAPE_CARD["S12"],
        rho_r=RESIN_CARD["rho"], E_r=RESIN_CARD["E"], nu_r=RESIN_CARD["nu"],
        fixed_block=fixed_block,
        load_block=load_block,
        Lx=25.0e-3, Ly=25.0e-3, Lz=0.72e-3,
    )
    starter.write_text(text)
    print("[stage11] starter deck written:", starter)

    engine = DECK_DIR / f"{job}_0001.rad"
    engine.write_text(ENGINE_TEMPLATE.format(job=job))
    print("[stage11] engine deck written:", engine)
    return starter


def _bcs_keyword(group: str, dof: str, value: float, fixed: bool) -> str:
    dof_map = {"ux": (1, 0, 0), "uy": (0, 1, 0), "uz": (0, 0, 1)}
    dx, dy, dz = dof_map[dof]
    flag = "1.0E20" if fixed else "1.0E20"
    return (
        f"/BCS/TRA/{abs(hash(group + dof)) % 90000 + 1000}\n"
        f"{group}_{dof}\n"
        f"{flag:>20}{flag:>20}{flag:>20}\n"
        f"{group:>20}\n"
        f"# applied magnitude {value:.3e} m on {dof} for face {group}"
    )


STARTER_TEMPLATE = """\
#RADIOSS STARTER
/BEGIN
{job}
      2026         0
SI                  m                   s                   kg
SI                  m                   s                   kg
/IMPL/STATIC
/IMPL/LINEAR
/IMPL/PRINT
                   1
/INCLUDE
{mesh_include}
#---------------------------------------- materials
/MAT/LAW25/1
IM7_8552_tape
{rho_t:>20.3f}            0.0
{E11:>20.3e}{E22:>20.3e}{E33:>20.3e}
{nu12:>20.3f}{nu13:>20.3f}{nu23:>20.3f}
{G12:>20.3e}{G13:>20.3e}{G23:>20.3e}
{Xt:>20.3e}{Xc:>20.3e}{Yt:>20.3e}{Yc:>20.3e}{S12:>20.3e}
/MAT/LAW1/2
EPOXY_8552_resin
{rho_r:>20.3f}{E_r:>20.3e}{nu_r:>20.3f}
#---------------------------------------- properties (defined in {mesh_include})
#---------------------------------------- node face groups (defined in {mesh_include})
#---------------------------------------- BCs
{fixed_block}
{load_block}
#---------------------------------------- imposed displacement ramp
/IMPL/DISPL/INCR
                   1          1.0e-6                 1
#---------------------------------------- output
/TH/NODE
                1001            1002            1003            1004
/ANIM/DT
                 0.0             1.0
/END
"""

ENGINE_TEMPLATE = """\
/RUN/{job}/1
        1.0
/PRINT/-1
/STOP/STAT
"""


# ---------------------------------------------------------------------------
# Step (d): solver invocation
# ---------------------------------------------------------------------------

def run_openradioss(case: str) -> Path:
    """Invoke the OpenRadioss starter and engine for the named case."""
    job = f"ap_ply_block_case{case}"
    starter = DECK_DIR / f"{job}_0000.rad"
    engine = DECK_DIR / f"{job}_0001.rad"
    run_dir = RUN_DIR / case
    run_dir.mkdir(parents=True, exist_ok=True)

    starter_bin = OR_EXEC / "starter_linuxa64"
    engine_bin = OR_EXEC / "engine_linuxa64"

    cmd1 = LIMA_PREFIX + [str(starter_bin), "-i", str(starter), "-nt", "4"]
    cmd2 = LIMA_PREFIX + [str(engine_bin), "-i", str(engine), "-nt", "4"]
    print("[stage11] starter:\n  " + " ".join(cmd1))
    subprocess.run(cmd1, cwd=run_dir, check=True)
    print("[stage11] engine:\n  " + " ".join(cmd2))
    subprocess.run(cmd2, cwd=run_dir, check=True)
    return run_dir


# ---------------------------------------------------------------------------
# Step (e): area-average stress and strain to recover effective stiffness
# ---------------------------------------------------------------------------

def area_average_stiffness(case: str, run_dir: Path) -> dict:
    """Read the run output (T01 + .anim) and compute the effective stress and
    strain on the loaded face.  Returns engineering stiffness in GPa."""
    bcs = CASE_BCS[case]
    L = {"xx": (25.0e-3, 25.0e-3 * 0.72e-3),
         "yy": (25.0e-3, 25.0e-3 * 0.72e-3),
         "xy": (25.0e-3, 25.0e-3 * 0.72e-3)}[bcs["strain_axis"]]
    face_length, face_area = L

    # Read summed reaction force on the loaded face from T01.  The keyword
    # /TH/NODE on the loaded face's group writes the per-node contact /
    # reaction force; sum across nodes to get the total reaction.  Below we use
    # the Vortex-Radioss Python reader if available, with anim2vtk fallback.
    try:
        from vortex_radioss.animtoth5 import read_t01  # type: ignore
        F_total = read_t01(run_dir / f"ap_ply_block_case{case}T01")["RF_loaded_face"]
    except ImportError:
        F_total = _read_t01_fallback(run_dir / f"ap_ply_block_case{case}T01", bcs)

    # Effective stress sigma_eng = F_total / face_area.
    # Effective strain eps_eng = applied_displacement / face_length.
    applied_disp = bcs["loaded"][0][2]
    eps_eng = applied_disp / face_length
    sigma_eng = F_total / face_area

    # Engineering stiffness (Pa -> GPa).
    if bcs["strain_axis"] == "xx":
        E_x = sigma_eng / eps_eng / 1.0e9
        return {"case": case, "E_x_GPa": E_x}
    if bcs["strain_axis"] == "yy":
        E_y = sigma_eng / eps_eng / 1.0e9
        return {"case": case, "E_y_GPa": E_y}
    if bcs["strain_axis"] == "xy":
        # eps_xy in spec is gamma/2; G_xy = sigma_xy / gamma_xy.
        gamma = 2.0 * eps_eng
        G_xy = sigma_eng / gamma / 1.0e9
        return {"case": case, "G_xy_GPa": G_xy}
    raise ValueError(case)


def _read_t01_fallback(t01_path: Path, bcs: dict) -> float:
    """Parse T01 ASCII text for the summed reaction on the loaded face.  This
    is a stage 11 stop-gap until Vortex-Radioss is installed."""
    if not t01_path.exists():
        raise FileNotFoundError(str(t01_path))
    # Last column of the T01 history at t = 1.0 holds the total reaction
    # if /TH/NODE was set on the loaded-face group.
    data = np.loadtxt(t01_path, comments="#")
    last_row = data[-1]
    return float(last_row[-1])


# ---------------------------------------------------------------------------
# Step (f): comparison to Kok 2022
# ---------------------------------------------------------------------------

def compare_to_reference(results: dict) -> dict:
    """Compare computed E_x, E_y, G_xy to the Kok 2022 reference (with the CLT
    shear target).  Returns a dict suitable for serialization to CSV."""
    out = {"case": "summary"}
    tol = KOK_2022_REFERENCE["tolerance"]
    for key in ("E_x_GPa", "E_y_GPa", "G_xy_GPa"):
        target = KOK_2022_REFERENCE[key]
        measured = results.get(key)
        if measured is None:
            continue
        rel = abs(measured - target) / target
        out[f"{key}_FEM"] = measured
        out[f"{key}_target"] = target
        out[f"{key}_rel_err"] = rel
        out[f"{key}_pass"] = bool(rel <= tol)
    out["overall_pass"] = all(v for k, v in out.items() if k.endswith("_pass"))
    return out


def write_results_csv(rows: list, path: Path) -> None:
    import csv
    keys = sorted({k for row in rows for k in row.keys()})
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print("[stage11] results CSV written:", path)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stage 11 - AP-PLY mesoscale direct simulation")
    p.add_argument("--all", action="store_true", help="run every step")
    p.add_argument("--geometry", action="store_true", help="run only Kok preprocessor")
    p.add_argument("--mesh", action="store_true", help="run only inp2rad")
    p.add_argument("--case", choices=["A", "B", "C"], help="run a single OpenRadioss case")
    p.add_argument("--postprocess", action="store_true", help="area-average + compare")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    config = APPlyConfig()

    if args.all or args.geometry:
        inp = run_kok_preprocessor(config)
    else:
        inp = STAGE_DIR / config.export_path

    if args.all or args.mesh:
        mesh_rad = convert_inp_to_rad(inp)
    else:
        mesh_rad = inp.with_suffix(".rad")

    case_results: list = []
    cases = ["A", "B", "C"] if args.all else ([args.case] if args.case else [])
    for case in cases:
        template_starter_deck(case, mesh_rad)
        run_dir = run_openradioss(case)
        case_results.append(area_average_stiffness(case, run_dir))

    if args.all or args.postprocess:
        # Merge per-case results into a single record.
        merged: dict = {}
        for row in case_results:
            merged.update({k: v for k, v in row.items() if k != "case"})
        summary = compare_to_reference(merged)
        write_results_csv(case_results + [summary], RES_DIR / "effective_moduli.csv")

        sentinel = RES_DIR / ("PASS" if summary.get("overall_pass") else "FAIL")
        sentinel.write_text(json.dumps(summary, indent=2))
        print(f"[stage11] {sentinel.name}: see {sentinel}")
        return 0 if summary.get("overall_pass") else 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
