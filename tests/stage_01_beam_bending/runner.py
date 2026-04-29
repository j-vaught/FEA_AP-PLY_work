"""Stage 01 — Linear-elastic 3-pt and 4-pt bending of a solid prismatic beam.

Toolchain smoke test for the FEA_AP-PLY single-tool OpenRadioss pipeline.
Pipeline per spec.md §9: GMSH HEXA8 mesh -> Jinja2 .rad render (TODO deck.rad.j2)
-> inp2rad -> starter+engine in Lima -> .anim to .vtkhdf -> PyVista headless probe
-> 1% pass/fail vs Euler-Bernoulli. Author: J.C. Vaught.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import pathlib
import subprocess
import sys

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

STAGE_DIR = pathlib.Path(__file__).resolve().parent
TEMPLATE_PATH = STAGE_DIR / "deck.rad.j2"            # TODO: not yet written; see spec.md §6
RESULTS_CSV = STAGE_DIR / "results.csv"

# Lima VM name (matches master_plan.md §7 example "or" / "apptainer"); user override via CLI.
DEFAULT_LIMA_VM = "apptainer"
# Path to OpenRadioss exec dir *inside* the Lima VM.
DEFAULT_OR_EXEC = "/opt/openradioss/exec"
# Path to Kitware converter inside the Lima VM.
DEFAULT_VTKHDF_CONVERTER = "/opt/openradioss-to-vtkhdf/openradioss_to_vtkhdf.py"


@dataclasses.dataclass(frozen=True)
class BeamCase:
    """Geometry, material, and load definition for stage 01 (SI base units)."""
    L: float = 0.200          # m, length
    b: float = 0.020          # m, width
    h: float = 0.010          # m, height
    E: float = 6.89e10        # Pa, Young's modulus (6061-T6 per brief)
    nu: float = 0.33          # -, Poisson's ratio
    rho: float = 2700.0       # kg/m^3, density
    P: float = 1000.0         # N, total applied load (split as P/2,P/2 for 4-pt)

    @property
    def I(self) -> float:
        """Second moment of area about y-axis: I = b h^3 / 12 [m^4]."""
        return self.b * self.h**3 / 12.0

    def delta_3pt(self) -> float:
        """Euler-Bernoulli midspan deflection, 3-pt: P L^3 / (48 E I) [m]."""
        return self.P * self.L**3 / (48.0 * self.E * self.I)

    def delta_4pt(self) -> float:
        """Euler-Bernoulli midspan deflection, 4-pt with a = L/3:
        P a (3 L^2 - 4 a^2) / (24 E I) [m].
        """
        a = self.L / 3.0
        return self.P * a * (3 * self.L**2 - 4 * a**2) / (24.0 * self.E * self.I)


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def generate_mesh(case: BeamCase, mesh_level: str, out_inp: pathlib.Path) -> None:
    """Step 1. Build a structured HEXA8 mesh via the GMSH Python API.
    mesh_level in {M0, M1, M2} sets (n_x, n_y, n_z) per spec.md §3.
    TODO: gmsh transfinite hex over [0,L]x[-b/2,b/2]x[-h/2,h/2], tag the four required
    node sets (NSET_SUPPORT_LEFT/RIGHT, NSET_LOAD_*, NODE_MIDSPAN), write Abaqus .inp.
    """
    raise NotImplementedError("generate_mesh: GMSH Python API mesh generation not yet wired.")


def render_deck(case: BeamCase, load_case: str, mesh_inp: pathlib.Path,
                out_starter: pathlib.Path, out_engine: pathlib.Path) -> None:
    """Step 2. Render OpenRadioss .rad deck from Jinja2 template at TEMPLATE_PATH.
    TODO (deck.rad.j2 not yet written; see spec.md §6 for the keyword skeleton): fill
    {L, b, h, E, nu, rho, P, mesh_inp, load_case, force_per_node} into _0000.rad/_0001.rad.
    """
    raise NotImplementedError("render_deck: Jinja2 template not yet written; see spec.md §6.")


def _lima(vm: str, *cmd: str) -> None:
    """Run a command inside the Lima VM via `limactl shell <vm> -- <cmd...>`."""
    subprocess.run(["limactl", "shell", vm, "--", *cmd], check=True)


def inp_to_rad(mesh_inp: pathlib.Path, mesh_rad: pathlib.Path, lima_vm: str) -> None:
    """Step 3. Convert Abaqus .inp -> OpenRadioss .rad via inp2rad (OpenRadioss/Tools)."""
    _lima(lima_vm, "python3",
          "/opt/openradioss/Tools/input_converters/inp2rad/inp2rad.py",
          str(mesh_inp), "-o", str(mesh_rad))


def run_openradioss(starter_rad: pathlib.Path, engine_rad: pathlib.Path,
                    lima_vm: str, or_exec: str, n_threads: int) -> None:
    """Step 4. Invoke OpenRadioss starter then engine inside Lima.
    Implicit linear requires a MUMPS-linked build — see spec.md §10 risk #1.
    """
    _lima(lima_vm, f"{or_exec}/starter_linuxa64", "-i", str(starter_rad), "-nt", str(n_threads))
    _lima(lima_vm, f"{or_exec}/engine_linuxa64", "-i", str(engine_rad), "-nt", str(n_threads))


def anim_to_vtkhdf(anim_path: pathlib.Path, vtkhdf_path: pathlib.Path,
                   lima_vm: str, converter: str) -> None:
    """Step 5. Convert OpenRadioss .anim -> .vtkhdf via the Kitware converter."""
    _lima(lima_vm, "python3", converter, str(anim_path), "-o", str(vtkhdf_path))


def extract_midspan_uz(vtkhdf_path: pathlib.Path, case: BeamCase) -> float:
    """Step 6. Read .vtkhdf with PyVista headless, return u_z at (L/2, 0, -h/2) [m].
    Field name pinned by the Kitware converter version (typ. "Displacement" or "DISP").
    """
    import pyvista as pv  # lazy import so the rest of the module is import-light.
    pv.OFF_SCREEN = True
    grid = pv.read(str(vtkhdf_path))
    pid = grid.find_closest_point((case.L / 2.0, 0.0, -case.h / 2.0))
    disp = grid.point_data.get("Displacement", grid.point_data.get("DISP"))
    if disp is None:
        raise KeyError(f"Displacement field not found; available: {list(grid.point_data.keys())}")
    return float(disp[pid][2])  # z-component, metres


def evaluate(case: BeamCase, load_case: str, uz_fem: float) -> tuple[float, float, bool]:
    """Step 7. Compute relative error vs Euler-Bernoulli closed form; PASS if <= 1%."""
    delta_ref = -case.delta_3pt() if load_case == "3pt" else -case.delta_4pt()
    rel_err = abs(uz_fem - delta_ref) / abs(delta_ref)
    return delta_ref, rel_err, rel_err <= 0.01


def write_results_csv(rows: list[dict]) -> None:
    """Append rows to results.csv (created on first call)."""
    new_file = not RESULTS_CSV.exists()
    with RESULTS_CSV.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        if new_file:
            writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Stage 01 — beam bending toolchain runner.")
    p.add_argument("--load-case", choices=["3pt", "4pt", "both"], default="both")
    p.add_argument("--mesh", choices=["M0", "M1", "M2"], default="M1")
    p.add_argument("--solver", choices=["implicit", "explicit"], default="implicit",
                   help="implicit = /IMPL/LINEAR (needs MUMPS-linked build); "
                        "explicit = dynamic relaxation fallback (spec.md §10 risk #1).")
    p.add_argument("--lima-vm", default=DEFAULT_LIMA_VM)
    p.add_argument("--or-exec", default=DEFAULT_OR_EXEC)
    p.add_argument("--vtkhdf-converter", default=DEFAULT_VTKHDF_CONVERTER)
    p.add_argument("--n-threads", type=int, default=4)
    p.add_argument("--dry-run", action="store_true",
                   help="Print closed-form deflections and exit; no solver call.")
    args = p.parse_args(argv)

    case = BeamCase()
    print(f"[stage_01] L={case.L} m, b={case.b} m, h={case.h} m, "
          f"E={case.E:.3e} Pa, I={case.I:.6e} m^4, P={case.P} N")
    print(f"[stage_01] closed-form deflections (Euler-Bernoulli):")
    print(f"  3-pt midspan = {case.delta_3pt()*1e3:.4f} mm")
    print(f"  4-pt midspan = {case.delta_4pt()*1e3:.4f} mm")
    if args.dry_run:
        return 0

    cases = ["3pt", "4pt"] if args.load_case == "both" else [args.load_case]
    rows: list[dict] = []
    overall_pass = True

    for lc in cases:
        work = STAGE_DIR / f"run_{lc}_{args.mesh}_{args.solver}"
        work.mkdir(exist_ok=True)
        mesh_inp = work / "mesh.inp"
        mesh_rad = work / "mesh.rad"
        starter = work / f"bend_{lc}_0000.rad"
        engine = work / f"bend_{lc}_0001.rad"
        anim = work / f"bend_{lc}_0001.anim"
        vtkhdf = work / f"bend_{lc}.vtkhdf"

        generate_mesh(case, args.mesh, mesh_inp)
        render_deck(case, lc, mesh_inp, starter, engine)
        inp_to_rad(mesh_inp, mesh_rad, args.lima_vm)
        run_openradioss(starter, engine, args.lima_vm, args.or_exec, args.n_threads)
        anim_to_vtkhdf(anim, vtkhdf, args.lima_vm, args.vtkhdf_converter)

        uz = extract_midspan_uz(vtkhdf, case)
        delta_ref, rel_err, passed = evaluate(case, lc, uz)
        overall_pass = overall_pass and passed
        print(f"[{lc}/{args.mesh}] u_z_FEM = {uz*1e3:+.4f} mm, "
              f"closed-form = {delta_ref*1e3:+.4f} mm, "
              f"rel_err = {rel_err*100:.3f}%  ->  {'PASS' if passed else 'FAIL'}")
        rows.append({
            "load_case": lc, "mesh": args.mesh, "solver": args.solver,
            "uz_fem_mm": uz * 1e3, "delta_ref_mm": delta_ref * 1e3,
            "rel_err": rel_err, "pass": int(passed),
        })

    write_results_csv(rows)
    return 0 if overall_pass else 1


if __name__ == "__main__":
    sys.exit(main())
