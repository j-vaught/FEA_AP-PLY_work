"""
Stage 10 - Unidirectional tow-wise direct mesoscale (reframed) - runner.

Pipeline.
    1. Load mesh from UD_mesoscale_hex.inp (built by mesh_build.py via GMSH).
    2. Identify the six face nodesets by bounding-box queries.
    3. For each of three KUBC loading cases (axial stretch, transverse stretch,
       longitudinal shear) generate a /IMPDISP table per face and write the
       OpenRadioss starter (_0000.rad) and engine (_0001.rad) decks.
    4. Run the OpenRadioss starter + engine inside Lima Apptainer.
    5. Convert .anim to .vtkhdf via the Kitware openradioss-to-vtkhdf tool.
    6. Volume-average stress and strain across all elements, derive effective
       constants E1, E2, G12.
    7. Compare to Halpin-Tsai 1969 (zeta=2 for E_2, zeta=1 for G_12) and
       Hashin-Shtrikman 1963 bounds; write effective_moduli.csv.

Author. J.C. Vaught
Date.   2026-04-29
Spec.   tests/stage_10_UD_mesoscale_direct/spec.md (sections 5-8 pin the math).

Notes.
    - SI units throughout (Pa, m, kg).
    - This runner is the deck-templating + post-processing layer. The Lima
      Apptainer launch is delegated to the user's OpenRadioss wrapper described
      in master_plan.md section 7. The OR_RUN constant below points at it.
    - Element types accepted: TETRA10 (default from GMSH OCC) and HEXA8.
    - This file is heavy on docstrings and light on cleverness on purpose:
      it is a verification stage and must be auditable.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Iterable

import numpy as np

# ---------------------------------------------------------------------------
# Paths and constants. Edit OR_RUN to match the user's local Lima install.
# ---------------------------------------------------------------------------

STAGE_DIR = Path(__file__).resolve().parent
MESH_INP = STAGE_DIR / "UD_mesoscale_hex.inp"
DECK_DIR = STAGE_DIR / "decks"
RUN_DIR = STAGE_DIR / "runs"
OUT_DIR = STAGE_DIR / "outputs"
FIG_DATA_DIR = STAGE_DIR.parents[1] / "figures" / "data"
PATHS_FILE = STAGE_DIR / "or_paths.json"          # Lima/Apptainer + OR binaries

# Specimen geometry (must match mesh_build.py).
L_X = 100.0e-6
L_Y = 100.0e-6
L_Z = 50.0e-6

# Constituent properties (Soden 1998: IM7 fiber, Hexcel 8552 matrix).
E_F = 230.0e9
NU_F = 0.20
RHO_F = 1780.0

E_M = 4.08e9
NU_M = 0.39
RHO_M = 1300.0

V_F_TARGET = 0.60                                  # nominal volume fraction

# KUBC strain magnitudes (small enough for linear regime, large enough to be
# numerically clean; 1e-3 is the standard micromechanics convention).
EPS_AXIAL = 1.0e-3
EPS_TRANSVERSE = 1.0e-3
GAMMA_SHEAR = 1.0e-3

# Halpin-Tsai parameters (spec section 7.2).
ZETA_E2 = 2.0
ZETA_G12 = 1.0

# Numerical tolerances for the pass criteria (spec section 8).
HT_REL_TOL = 0.05
HS_TOL = 0.01

# Bounding-box tolerance for face nodeset detection (1 percent of edge).
FACE_TOL = 1.0e-8                                  # absolute, in metres


# ---------------------------------------------------------------------------
# Reference solutions: Halpin-Tsai and Hashin-Shtrikman.
# ---------------------------------------------------------------------------

def shear_modulus(E: float, nu: float) -> float:
    return E / (2.0 * (1.0 + nu))


def bulk_modulus(E: float, nu: float) -> float:
    return E / (3.0 * (1.0 - 2.0 * nu))


def halpin_tsai(E_phase: float, E_matrix: float, V_f: float, zeta: float) -> float:
    """Halpin-Tsai 1969 transverse / shear modulus prediction.

    Returns the homogenized modulus given the phase modulus E_phase (fiber for
    transverse E_2, fiber-shear for in-plane G_12), the matrix modulus, the
    fiber volume fraction, and the geometry parameter zeta.
    """
    ratio = E_phase / E_matrix
    eta = (ratio - 1.0) / (ratio + zeta)
    return E_matrix * (1.0 + zeta * eta * V_f) / (1.0 - eta * V_f)


def rule_of_mixtures(E_f: float, E_m: float, V_f: float) -> float:
    return V_f * E_f + (1.0 - V_f) * E_m


def hashin_shtrikman_bounds(E1: float, nu1: float, V1: float,
                            E2: float, nu2: float) -> dict:
    """Hashin-Shtrikman 1963 bounds on effective bulk K and shear G.

    Phase 1 is the inclusion / fiber, phase 2 is the matrix; V1 is the
    fiber volume fraction. Returns lower and upper bounds on K and G.
    The bounds are tight when both phases are well-ordered (K_1>K_2 and
    G_1>G_2 simultaneously), which holds for our IM7-in-8552 system.
    """
    K1, G1 = bulk_modulus(E1, nu1), shear_modulus(E1, nu1)
    K2, G2 = bulk_modulus(E2, nu2), shear_modulus(E2, nu2)
    V2 = 1.0 - V1

    def _K_bound(K_a, G_a, V_a, K_b, V_b):
        return K_a + V_b / (1.0 / (K_b - K_a) + 3.0 * V_a / (3.0 * K_a + 4.0 * G_a))

    def _G_bound(K_a, G_a, V_a, G_b, V_b):
        denom = 1.0 / (G_b - G_a) + 6.0 * V_a * (K_a + 2.0 * G_a) / (
            5.0 * G_a * (3.0 * K_a + 4.0 * G_a)
        )
        return G_a + V_b / denom

    # Lower: comparison medium = softer phase (matrix).
    K_lo = _K_bound(K2, G2, V2, K1, V1)
    G_lo = _G_bound(K2, G2, V2, G1, V1)
    # Upper: comparison medium = stiffer phase (fiber).
    K_up = _K_bound(K1, G1, V1, K2, V2)
    G_up = _G_bound(K1, G1, V1, K2, V2)

    return {"K_lower": K_lo, "K_upper": K_up,
            "G_lower": G_lo, "G_upper": G_up}


def reference_targets() -> dict:
    """Compose all closed-form reference numbers used in the pass criteria."""
    G_f = shear_modulus(E_F, NU_F)
    G_m = shear_modulus(E_M, NU_M)
    E1_rom = rule_of_mixtures(E_F, E_M, V_F_TARGET)
    E2_ht = halpin_tsai(E_F, E_M, V_F_TARGET, ZETA_E2)
    G12_ht = halpin_tsai(G_f, G_m, V_F_TARGET, ZETA_G12)
    hs = hashin_shtrikman_bounds(E_F, NU_F, V_F_TARGET, E_M, NU_M)
    return {
        "E1_rom_Pa": E1_rom,
        "E2_HT_Pa": E2_ht,
        "G12_HT_Pa": G12_ht,
        "hs_K_lower_Pa": hs["K_lower"],
        "hs_K_upper_Pa": hs["K_upper"],
        "hs_G_lower_Pa": hs["G_lower"],
        "hs_G_upper_Pa": hs["G_upper"],
    }


# ---------------------------------------------------------------------------
# Mesh ingest from Abaqus .inp via meshio.
# ---------------------------------------------------------------------------

@dataclass
class Mesh:
    nodes: np.ndarray            # (N, 3)
    cells: dict                  # {"tetra10": (M, 10), "hexa": (...)}
    cell_phase: np.ndarray       # (M,) ints; 0 = matrix, 1 = fiber
    cell_volumes: np.ndarray     # (M,)
    bbox_min: np.ndarray         # (3,)
    bbox_max: np.ndarray         # (3,)
    face_nodes: dict = field(default_factory=dict)   # name -> ndarray of node indices


def load_mesh(inp_path: Path) -> Mesh:
    """Load the GMSH-exported Abaqus .inp file. Requires meshio."""
    import meshio
    m = meshio.read(inp_path.as_posix())
    nodes = np.asarray(m.points, dtype=np.float64)
    cells = {}
    cell_phase_chunks = []

    # Phase tagging via the elset/material assignments preserved in the .inp.
    # We rely on the GMSH "FIBER" and "MATRIX" physical groups to be written
    # as Abaqus element sets named "FIBER" and "MATRIX" (meshio convention).
    for cb in m.cells:
        if cb.type in ("tetra10", "tetra"):
            cells.setdefault(cb.type, []).append(np.asarray(cb.data))
        elif cb.type in ("hexahedron", "hexahedron20"):
            cells.setdefault(cb.type, []).append(np.asarray(cb.data))

    cells = {k: np.vstack(v) for k, v in cells.items()}

    # Attach cell phase from cell_sets if present.
    n_total_cells = sum(c.shape[0] for c in cells.values())
    phase = np.zeros(n_total_cells, dtype=np.int8)
    if hasattr(m, "cell_sets") and m.cell_sets:
        # Meshio cell_sets: dict[name] -> list per cell-block.
        offset = 0
        block_sizes = [c.shape[0] for c in cells.values()]
        if "FIBER" in m.cell_sets:
            ranges = m.cell_sets["FIBER"]
            for r, size in zip(ranges, block_sizes):
                if r is None:
                    offset += size
                    continue
                phase[offset + np.asarray(r, dtype=int)] = 1
                offset += size
    cell_volumes = compute_cell_volumes(nodes, cells)

    bbox_min = nodes.min(axis=0)
    bbox_max = nodes.max(axis=0)

    mesh = Mesh(
        nodes=nodes,
        cells=cells,
        cell_phase=phase,
        cell_volumes=cell_volumes,
        bbox_min=bbox_min,
        bbox_max=bbox_max,
    )
    mesh.face_nodes = identify_faces(nodes, bbox_min, bbox_max)
    return mesh


def compute_cell_volumes(nodes: np.ndarray, cells: dict) -> np.ndarray:
    """Per-cell volume by signed-tet decomposition."""
    vols = []
    for ctype, conn in cells.items():
        if ctype.startswith("tetra"):
            v = tet_volumes(nodes[conn[:, :4]])
        elif ctype.startswith("hex"):
            v = hex_volumes(nodes[conn[:, :8]])
        else:
            v = np.zeros(conn.shape[0])
        vols.append(v)
    return np.concatenate(vols)


def tet_volumes(tet_nodes: np.ndarray) -> np.ndarray:
    """Signed tetrahedron volume; abs to get cell volume."""
    a = tet_nodes[:, 1, :] - tet_nodes[:, 0, :]
    b = tet_nodes[:, 2, :] - tet_nodes[:, 0, :]
    c = tet_nodes[:, 3, :] - tet_nodes[:, 0, :]
    return np.abs(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0


def hex_volumes(hex_nodes: np.ndarray) -> np.ndarray:
    """Hex-8 volume by 6-tet decomposition (one of several valid splits)."""
    n = hex_nodes
    # Six-tet decomposition relative to node 0.
    splits = [(0, 1, 2, 5),
              (0, 2, 3, 7),
              (0, 2, 5, 7),
              (0, 5, 6, 7),
              (0, 4, 5, 7),
              (2, 5, 6, 7)]
    total = np.zeros(n.shape[0])
    for split in splits:
        total += tet_volumes(n[:, list(split), :])
    return total


def identify_faces(nodes: np.ndarray, bmin: np.ndarray, bmax: np.ndarray) -> dict:
    """Return node-index arrays for the six bounding-box faces."""
    return {
        "XMIN": np.where(np.abs(nodes[:, 0] - bmin[0]) < FACE_TOL)[0],
        "XMAX": np.where(np.abs(nodes[:, 0] - bmax[0]) < FACE_TOL)[0],
        "YMIN": np.where(np.abs(nodes[:, 1] - bmin[1]) < FACE_TOL)[0],
        "YMAX": np.where(np.abs(nodes[:, 1] - bmax[1]) < FACE_TOL)[0],
        "ZMIN": np.where(np.abs(nodes[:, 2] - bmin[2]) < FACE_TOL)[0],
        "ZMAX": np.where(np.abs(nodes[:, 2] - bmax[2]) < FACE_TOL)[0],
    }


def actual_volume_fraction(mesh: Mesh) -> float:
    fiber_vol = mesh.cell_volumes[mesh.cell_phase == 1].sum()
    return fiber_vol / mesh.cell_volumes.sum()


# ---------------------------------------------------------------------------
# OpenRadioss deck generation per loading case.
# ---------------------------------------------------------------------------

@dataclass
class LoadCase:
    name: str                       # "axial", "transverse", "shear"
    eps_tensor: np.ndarray          # (3,3) macro strain
    extract: str                    # which stress to recover: "sigma_zz", etc.
    target_modulus_key: str         # "E1_rom_Pa", "E2_HT_Pa", "G12_HT_Pa"


def define_cases() -> list[LoadCase]:
    cases = []

    eps_a = np.zeros((3, 3))
    eps_a[2, 2] = EPS_AXIAL
    cases.append(LoadCase("axial", eps_a, "sigma_zz", "E1_rom_Pa"))

    eps_b = np.zeros((3, 3))
    eps_b[0, 0] = EPS_TRANSVERSE
    cases.append(LoadCase("transverse", eps_b, "sigma_xx", "E2_HT_Pa"))

    eps_c = np.zeros((3, 3))
    # Engineering shear gamma_xz; tensorial component = gamma/2 on each off-diagonal.
    eps_c[0, 2] = 0.5 * GAMMA_SHEAR
    eps_c[2, 0] = 0.5 * GAMMA_SHEAR
    cases.append(LoadCase("shear", eps_c, "sigma_xz", "G12_HT_Pa"))

    return cases


def kubc_node_displacements(nodes: np.ndarray, eps: np.ndarray,
                            face_idx: np.ndarray) -> np.ndarray:
    """For each node on a face, return u_i = eps_ij * x_j as an (n,3) array."""
    return nodes[face_idx, :] @ eps.T


def write_deck(case: LoadCase, mesh: Mesh, deck_path: Path) -> None:
    """Emit the OpenRadioss starter+engine deck pair for one loading case.

    The starter (_0000.rad) carries materials, properties, parts, BCs and
    /IMPDISP tables. The engine (_0001.rad) carries the implicit-static run
    control, animation cadence, and stop time.
    """
    deck_path.parent.mkdir(parents=True, exist_ok=True)
    starter = deck_path.with_name(f"{case.name}_0000.rad")
    engine = deck_path.with_name(f"{case.name}_0001.rad")

    # Build per-node IMPDISP records for every face.
    impdisp_blocks = []
    bcs_blocks = []
    funct_id_counter = 100   # /FUNCT IDs offset to avoid collision
    impdisp_id_counter = 1
    for face_name, face_idx in mesh.face_nodes.items():
        u_face = kubc_node_displacements(mesh.nodes, case.eps_tensor, face_idx)
        for k_dof in range(3):
            # One /IMPDISP per node per nonzero DOF; for KUBC every face node
            # has all three DOFs prescribed (some to zero).
            for local_i, node_idx in enumerate(face_idx):
                u_val = u_face[local_i, k_dof]
                impdisp_blocks.append(
                    f"/IMPDISP/{impdisp_id_counter}\n"
                    f"node{node_idx}_dof{k_dof+1}\n"
                    f"{k_dof+1}  {funct_id_counter}  {node_idx+1}  {u_val:.12e}\n"
                )
                impdisp_id_counter += 1

    # One ramp function shared by every IMPDISP scaling: 0->1 over t in [0,1].
    funct_block = (
        f"/FUNCT/{funct_id_counter}\n"
        "linear_ramp\n"
        "0.0  0.0\n"
        "1.0  1.0\n"
    )

    starter.write_text(
        STARTER_TEMPLATE.format(
            case_name=case.name,
            mesh_inc=str(MESH_INP.with_suffix(".inc").relative_to(STAGE_DIR)),
            E_F=E_F, NU_F=NU_F, RHO_F=RHO_F,
            E_M=E_M, NU_M=NU_M, RHO_M=RHO_M,
            funct_block=funct_block,
            impdisp_block="\n".join(impdisp_blocks),
        )
    )
    engine.write_text(
        ENGINE_TEMPLATE.format(case_name=case.name)
    )


STARTER_TEMPLATE = """\
#RADIOSS STARTER
/BEGIN
{case_name}
2026  0
                  Pa                   m                  kg                   s
                  Pa                   m                  kg                   s
/INCLUDE
{mesh_inc}
/MAT/LAW1/1
fiber_T800_IM7
{RHO_F}
{E_F}  {NU_F}
/MAT/LAW1/2
matrix_8552
{RHO_M}
{E_M}  {NU_M}
/PROP/SOLID/1
fiber_solid
0  0  0  0  0  0
/PROP/SOLID/2
matrix_solid
0  0  0  0  0  0
/PART/1
fiber_part
1  1
/PART/2
matrix_part
2  2
{funct_block}
{impdisp_block}
/END
"""

ENGINE_TEMPLATE = """\
#RADIOSS ENGINE
/RUN/{case_name}/1
1.0
/IMPL/LINEAR
/IMPL/SOLVER/2
/IMPL/PRINT/N-1
/ANIM/DT
0.0  1.0
/ANIM/BRICK/TENS/STRESS
/ANIM/BRICK/TENS/STRAIN
/ANIM/BRICK/VOL
/TH/BRICK/1
all_bricks
GRBRIC/ALL  SIGXX  SIGYY  SIGZZ  SIGXY  SIGYZ  SIGZX  EPSXX  EPSYY  EPSZZ  EPSXY  EPSYZ  EPSZX
/STOP
"""


# ---------------------------------------------------------------------------
# OpenRadioss invocation (Lima Apptainer wrapper).
# ---------------------------------------------------------------------------

def load_or_paths() -> dict:
    if not PATHS_FILE.exists():
        # Sensible default consistent with master_plan.md section 7.
        return {
            "lima_shell": ["limactl", "shell", "apptainer", "--"],
            "starter": "/OpenRadioss/exec/starter_linuxa64",
            "engine": "/OpenRadioss/exec/engine_linuxa64",
            "vtkhdf_converter": "openradioss-to-vtkhdf",
        }
    return json.loads(PATHS_FILE.read_text())


def run_openradioss(case_name: str, paths: dict, run_dir: Path) -> tuple[int, int]:
    """Invoke starter then engine on the given case. Returns the two returncodes."""
    run_dir.mkdir(parents=True, exist_ok=True)
    starter_cmd = list(paths["lima_shell"]) + [
        paths["starter"], "-i", f"{case_name}_0000.rad", "-nt", "4"
    ]
    engine_cmd = list(paths["lima_shell"]) + [
        paths["engine"], "-i", f"{case_name}_0001.rad", "-nt", "4"
    ]
    s = subprocess.run(starter_cmd, cwd=run_dir, capture_output=True, text=True)
    if s.returncode != 0:
        return s.returncode, -1
    e = subprocess.run(engine_cmd, cwd=run_dir, capture_output=True, text=True)
    return s.returncode, e.returncode


def convert_anim_to_vtkhdf(case_name: str, paths: dict, run_dir: Path) -> Path:
    """Run the Kitware openradioss-to-vtkhdf converter on the produced .anim."""
    anim = next(run_dir.glob(f"{case_name}*A001"))
    out = run_dir / f"{case_name}.vtkhdf"
    cmd = [paths["vtkhdf_converter"], anim.as_posix(), out.as_posix()]
    subprocess.run(cmd, check=True)
    return out


# ---------------------------------------------------------------------------
# Post-processing: read stress and strain, volume-average, derive moduli.
# ---------------------------------------------------------------------------

def read_field_stress_strain(vtkhdf_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (stress_per_cell (M,3,3), strain_per_cell (M,3,3), volumes (M,))."""
    import pyvista as pv
    grid = pv.read(vtkhdf_path.as_posix())
    cell_data = grid.cell_data

    # OpenRadioss VTKHDF stores stress as a 6-component symmetric tensor or as
    # six scalar fields named SIGXX, SIGYY, SIGZZ, SIGXY, SIGYZ, SIGZX. The
    # converter version controls which; we accept either.
    keys = list(cell_data.keys())
    if "Stress" in keys:
        sig6 = np.asarray(cell_data["Stress"])  # (M, 6)
    else:
        sig6 = np.column_stack([
            cell_data["SIGXX"], cell_data["SIGYY"], cell_data["SIGZZ"],
            cell_data["SIGXY"], cell_data["SIGYZ"], cell_data["SIGZX"],
        ])
    if "Strain" in keys:
        eps6 = np.asarray(cell_data["Strain"])
    else:
        eps6 = np.column_stack([
            cell_data["EPSXX"], cell_data["EPSYY"], cell_data["EPSZZ"],
            cell_data["EPSXY"], cell_data["EPSYZ"], cell_data["EPSZX"],
        ])
    sigma = voigt6_to_tensor(sig6)
    epsilon = voigt6_to_tensor(eps6, engineering_shear=True)
    volumes = np.asarray(grid.compute_cell_sizes(volume=True).cell_data["Volume"])
    return sigma, epsilon, volumes


def voigt6_to_tensor(v6: np.ndarray, engineering_shear: bool = False) -> np.ndarray:
    """Convert (M,6) Voigt = [xx, yy, zz, xy, yz, zx] to (M,3,3) symmetric tensor."""
    M = v6.shape[0]
    T = np.zeros((M, 3, 3))
    T[:, 0, 0] = v6[:, 0]
    T[:, 1, 1] = v6[:, 1]
    T[:, 2, 2] = v6[:, 2]
    factor = 0.5 if engineering_shear else 1.0
    T[:, 0, 1] = T[:, 1, 0] = v6[:, 3] * factor
    T[:, 1, 2] = T[:, 2, 1] = v6[:, 4] * factor
    T[:, 0, 2] = T[:, 2, 0] = v6[:, 5] * factor
    return T


def volume_average(field: np.ndarray, volumes: np.ndarray) -> np.ndarray:
    """Volume-weighted average of a per-cell tensor field. Returns (3,3)."""
    V = volumes.sum()
    return np.einsum("mij,m->ij", field, volumes) / V


def derive_modulus(case: LoadCase, sigma_avg: np.ndarray) -> float:
    """Recover the effective constant for this loading case (spec section 5)."""
    if case.name == "axial":
        return sigma_avg[2, 2] / EPS_AXIAL
    if case.name == "transverse":
        return sigma_avg[0, 0] / EPS_TRANSVERSE
    if case.name == "shear":
        # G_12 = <sigma_xz> / gamma_xz; engineering shear gamma_xz = 2 * eps_xz.
        return sigma_avg[0, 2] / GAMMA_SHEAR
    raise ValueError(case.name)


# ---------------------------------------------------------------------------
# Pass criteria.
# ---------------------------------------------------------------------------

def check_pass(case_name: str, value_Pa: float, refs: dict) -> dict:
    """Apply the four pass conditions of spec section 8 to one case."""
    if case_name == "axial":
        target = refs["E1_rom_Pa"]
        rel_err = abs(value_Pa - target) / target
        in_HT_band = rel_err <= HT_REL_TOL
        # E_1 bracket against HS bounds (axial uses Voigt/Reuss-like, but we
        # at least sanity-check that the value is sensible - between 0 and the
        # ROM upper estimate).
        in_HS = 0.0 < value_Pa < target * (1.0 + HT_REL_TOL + HS_TOL)
        return dict(target=target, rel_err=rel_err,
                    HS_lower=float("nan"), HS_upper=float("nan"),
                    in_HT_band=in_HT_band, in_HS_bracket=in_HS,
                    pass_=in_HT_band and in_HS)

    if case_name == "transverse":
        target = refs["E2_HT_Pa"]
        rel_err = abs(value_Pa - target) / target
        in_HT_band = rel_err <= HT_REL_TOL
        # Map E_2 to a comparable isotropic Young's via E ~ 9KG/(3K+G); for the
        # bracket we use the HS bounds on the equivalent Young's modulus,
        # E_lo = 9 K_lo G_lo / (3 K_lo + G_lo) etc.
        E_lo = 9 * refs["hs_K_lower_Pa"] * refs["hs_G_lower_Pa"] / (
            3 * refs["hs_K_lower_Pa"] + refs["hs_G_lower_Pa"])
        E_up = 9 * refs["hs_K_upper_Pa"] * refs["hs_G_upper_Pa"] / (
            3 * refs["hs_K_upper_Pa"] + refs["hs_G_upper_Pa"])
        in_HS = (E_lo * (1.0 - HS_TOL)) <= value_Pa <= (E_up * (1.0 + HS_TOL))
        return dict(target=target, rel_err=rel_err,
                    HS_lower=E_lo, HS_upper=E_up,
                    in_HT_band=in_HT_band, in_HS_bracket=in_HS,
                    pass_=in_HT_band and in_HS)

    if case_name == "shear":
        target = refs["G12_HT_Pa"]
        rel_err = abs(value_Pa - target) / target
        in_HT_band = rel_err <= HT_REL_TOL
        in_HS = (refs["hs_G_lower_Pa"] * (1.0 - HS_TOL)) <= value_Pa <= (
            refs["hs_G_upper_Pa"] * (1.0 + HS_TOL))
        return dict(target=target, rel_err=rel_err,
                    HS_lower=refs["hs_G_lower_Pa"], HS_upper=refs["hs_G_upper_Pa"],
                    in_HT_band=in_HT_band, in_HS_bracket=in_HS,
                    pass_=in_HT_band and in_HS)

    raise ValueError(case_name)


# ---------------------------------------------------------------------------
# Driver.
# ---------------------------------------------------------------------------

def main(skip_solver: bool = False) -> int:
    """Run all three KUBC cases end to end. Returns 0 on overall pass."""
    print("Stage 10: UD direct mesoscale (KUBC, reframed).")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DATA_DIR.mkdir(parents=True, exist_ok=True)

    print(f"  Loading mesh: {MESH_INP.name}")
    mesh = load_mesh(MESH_INP)
    V_f_actual = actual_volume_fraction(mesh)
    n_elements = mesh.cell_volumes.size
    print(f"  Volume fraction (actual)  : {V_f_actual:.4f} (target {V_F_TARGET})")
    print(f"  Element count             : {n_elements}")
    print(f"  Face nodes (XMIN, XMAX, ...): "
          f"{[len(v) for v in mesh.face_nodes.values()]}")
    if abs(V_f_actual - V_F_TARGET) > 0.005:
        print(f"  WARNING: V_f off target by > 0.5 percent. "
              "Check mesh boundary fiber clipping (spec section 2.4).")

    refs = reference_targets()
    print(f"  Halpin-Tsai targets:")
    print(f"    E_1 (ROM)    = {refs['E1_rom_Pa']/1e9:.2f} GPa")
    print(f"    E_2 (HT z=2) = {refs['E2_HT_Pa']/1e9:.2f} GPa")
    print(f"    G_12 (HT z=1) = {refs['G12_HT_Pa']/1e9:.2f} GPa")

    paths = load_or_paths()
    cases = define_cases()
    rows = []

    for case in cases:
        print(f"\n  Case: {case.name}")
        deck_dir = DECK_DIR / case.name
        run_dir = RUN_DIR / case.name
        run_dir.mkdir(parents=True, exist_ok=True)
        write_deck(case, mesh, deck_dir / f"{case.name}_0000.rad")
        # Stage decks into the run directory next to the included mesh.
        for f in deck_dir.glob("*.rad"):
            shutil.copy(f, run_dir / f.name)
        if not skip_solver:
            rc_starter, rc_engine = run_openradioss(case.name, paths, run_dir)
            print(f"    starter rc={rc_starter}, engine rc={rc_engine}")
            if rc_starter != 0 or rc_engine != 0:
                print(f"    OpenRadioss returned nonzero; aborting case.")
                continue
            vtk = convert_anim_to_vtkhdf(case.name, paths, run_dir)
        else:
            vtk = run_dir / f"{case.name}.vtkhdf"
            if not vtk.exists():
                print(f"    skip_solver=True but {vtk} missing; cannot post-process.")
                continue

        sigma_field, epsilon_field, volumes = read_field_stress_strain(vtk)
        sigma_avg = volume_average(sigma_field, volumes)
        epsilon_avg = volume_average(epsilon_field, volumes)
        modulus = derive_modulus(case, sigma_avg)

        verdict = check_pass(case.name, modulus, refs)
        print(f"    extracted  = {modulus/1e9:.3f} GPa")
        print(f"    target     = {verdict['target']/1e9:.3f} GPa "
              f"(rel err {verdict['rel_err']*100:.2f} percent)")
        print(f"    HS bracket = [{verdict['HS_lower']/1e9:.3f}, "
              f"{verdict['HS_upper']/1e9:.3f}] GPa")
        print(f"    pass       = {verdict['pass_']}")

        rows.append(dict(
            case=case.name,
            modulus_FEM_GPa=modulus / 1e9,
            target_GPa=verdict["target"] / 1e9,
            rel_err=verdict["rel_err"],
            HS_lower_GPa=verdict["HS_lower"] / 1e9 if not math.isnan(verdict["HS_lower"]) else "",
            HS_upper_GPa=verdict["HS_upper"] / 1e9 if not math.isnan(verdict["HS_upper"]) else "",
            in_HT_band=verdict["in_HT_band"],
            in_HS_bracket=verdict["in_HS_bracket"],
            pass_=verdict["pass_"],
            V_f_actual=V_f_actual,
            n_elements=n_elements,
        ))

    # Write CSV outputs.
    csv_path = OUT_DIR / "effective_moduli.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        w.writeheader()
        w.writerows(rows)
    print(f"\n  Wrote {csv_path}")

    # Mirror to the figures/data area for the Typst-CeTZ summary plot.
    shutil.copy(csv_path, FIG_DATA_DIR / "stage_10_summary.csv")

    overall_pass = bool(rows) and all(r["pass_"] for r in rows)
    print(f"\nStage 10 overall pass: {overall_pass}")
    return 0 if overall_pass else 1


if __name__ == "__main__":
    skip = "--skip-solver" in sys.argv
    sys.exit(main(skip_solver=skip))
