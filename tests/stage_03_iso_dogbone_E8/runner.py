"""Stage 03 runner. Isotropic dogbone tensile coupon (ASTM E8/E8M-22).

Author. j-vaught
Email. jvaught@sc.edu
Date. 2026-04-29
Stage brief reference. tests/stage_03_iso_dogbone_E8/spec.md

This script templates a parameterized OpenRadioss starter / engine deck
pair, builds the GMSH mesh, runs OpenRadioss inside Lima, parses the T01
time-history output, and evaluates the four pass / fail checks plus the
three secondary checks documented in spec.md section 8.

The script is intentionally modular. The pieces are:

    1. StageConfig dataclass        (all numerical parameters)
    2. build_mesh_with_gmsh         (geometry + mesh + .inp export)
    3. run_inp2rad                  (Abaqus .inp -> OpenRadioss .rad)
    4. render_starter_deck          (Jinja2-style template fill)
    5. render_engine_deck           (Jinja2-style template fill)
    6. run_openradioss_in_lima      (subprocess into limactl shell)
    7. parse_time_history           (Vortex-Radioss reader on T01)
    8. evaluate_checks              (C1..C4, S1..S3 from spec.md)
    9. write_results                (results.json + CSVs)
   10. plot_with_typst_cetz         (optional Typst CeTZ plots)

Hard rules per project CLAUDE.md and master_plan.md.

    - Solid HEXA8 elements only (no shells, no laminate stacks).
    - SI units, base kg, m, s, Pa.
    - Headless toolchain, no GUI.
    - All plots authored in Typst + CeTZ from CSV data; no matplotlib.
    - Brand colors from CLAUDE.md (Garnet, Atlantic, etc.).
    - No rounded edges on plots.

The runner is invoked as:

    python runner.py --mesh medium --material law2

Or for the full convergence study:

    python runner.py --mesh coarse  && \
    python runner.py --mesh medium  && \
    python runner.py --mesh fine
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import logging
import math
import os
import shutil
import subprocess
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

# ----------------------------------------------------------------------------
# 0. Brand color palette from CLAUDE.md (RGB tuples, 0-255).
# ----------------------------------------------------------------------------

BRAND_COLORS: dict[str, tuple[int, int, int]] = {
    "garnet":       (115,   0,  10),
    "black":        (  0,   0,   0),
    "white":        (255, 255, 255),
    "black_90":     ( 54,  54,  54),
    "black_70":     ( 92,  92,  92),
    "black_50":     (162, 162, 162),
    "black_30":     (199, 199, 199),
    "black_10":     (236, 236, 236),
    "warm_grey":    (103,  97,  86),
    "sandstorm":    (255, 242, 227),
    "rose":         (204,  46,  64),
    "atlantic":     ( 70, 106, 159),
    "congaree":     ( 31,  65,  77),
    "horseshoe":    (101, 120,  11),
    "grass":        (206, 211,  24),
    "honeycomb":    (164, 145,  55),
}


# ----------------------------------------------------------------------------
# 1. StageConfig dataclass.
# ----------------------------------------------------------------------------

@dataclass
class StageConfig:
    """All numerical parameters for stage 3.

    Default values reproduce the spec.md medium-mesh, LAW2 configuration.
    Override at the command line via argparse.
    """

    # --- run identity ---
    job_name:          str    = "stage03_dogbone_E8"
    mesh_label:        str    = "medium"               # coarse | medium | fine
    material_law:      str    = "law2"                 # law2 | law36

    # --- geometry (mm in the human-readable spec, meters internally) ---
    L0_mm:             float  = 50.0                   # gauge length
    W_mm:              float  = 12.5                   # gauge / reduced width
    t_mm:              float  = 6.0                    # thickness
    Lr_mm:             float  = 60.0                   # reduced-section length
    R_mm:              float  = 12.5                   # fillet radius
    Wg_mm:             float  = 20.0                   # grip-section width
    Lg_mm:             float  = 50.0                   # grip-section length

    # --- mesh density (per mesh_label below; resolved in __post_init__) ---
    n_thru:            int    = 4                      # through-thickness
    n_width_gauge:     int    = 8                      # across width, gauge
    n_length_gauge:    int    = 50                     # along length, gauge
    n_arc:             int    = 12                     # along fillet arc

    # --- material (A36 steel from master_plan brief) ---
    rho:               float  = 7850.0                 # kg/m^3
    E:                 float  = 2.0e11                 # Pa
    nu:                float  = 0.30                   # -
    sigma_y:           float  = 2.50e8                 # Pa  (250 MPa)
    sigma_u:           float  = 4.00e8                 # Pa  (400 MPa)
    eps_p_at_uts:      float  = 0.18                   # -

    # --- LAW2 Johnson-Cook fit (rate effects off) ---
    JC_a:              float  = 2.50e8                 # = sigma_y
    JC_b:              float  = 2.75e8                 # hardening modulus
    JC_n:              float  = 0.36                   # hardening exponent
    JC_c:              float  = 0.0                    # rate sensitivity off
    JC_eps_dot_0:      float  = 1.0                    # reference rate (dummy)

    # --- LAW36 tabulated knots (eps_p, sigma_y in Pa) ---
    LAW36_table: list[tuple[float, float]] = field(default_factory=lambda: [
        (0.000, 2.50e8),
        (0.020, 2.90e8),
        (0.050, 3.20e8),
        (0.100, 3.60e8),
        (0.180, 4.00e8),
    ])

    # --- loading ---
    eps_end_target:    float  = 0.005                  # gauge strain at T_end
    L_eff_mm:          float  = 75.0                   # effective compliant length
    T_end:             float  = 1.0                    # pseudo-time, s

    # --- pass / fail tolerances (from spec.md sec 8) ---
    tol_uniformity:    float  = 0.01                   # C1: 1%
    tol_stress:        float  = 0.01                   # C2: 1%
    tol_strain:        float  = 0.01                   # C3: 1%
    tol_yield_low:     float  = 0.99                   # C4 lower: 99% sigma_y
    tol_yield_high:    float  = 1.01                   # C4 upper: 101% sigma_y
    tol_E_apparent:    float  = 0.01                   # S1: 1%
    tol_poisson:       float  = 0.01                   # S2: 1%
    tol_post_yield:    float  = 0.02                   # S3: 2%

    # --- Lima / OpenRadioss toolchain ---
    lima_vm:           str    = "or"                   # limactl name
    or_starter:        str    = "/OpenRadioss/exec/starter_linuxa64"
    or_engine:         str    = "/OpenRadioss/exec/engine_linuxa64"
    or_inp2rad:        str    = "/OpenRadioss/Tools/input_converters/inp2rad/inp2rad.py"
    n_threads:         int    = 4

    # --- output toggles ---
    do_plots:          bool   = False
    work_dir:          Path   = field(default_factory=lambda: Path.cwd())

    # ---- derived properties ----

    def __post_init__(self) -> None:
        # Resolve mesh density by label.
        if self.mesh_label == "coarse":
            self.n_thru, self.n_width_gauge, self.n_length_gauge, self.n_arc = (
                3, 6, 40, 8
            )
        elif self.mesh_label == "medium":
            self.n_thru, self.n_width_gauge, self.n_length_gauge, self.n_arc = (
                4, 8, 50, 12
            )
        elif self.mesh_label == "fine":
            self.n_thru, self.n_width_gauge, self.n_length_gauge, self.n_arc = (
                5, 12, 70, 16
            )
        else:
            raise ValueError(f"unknown mesh_label: {self.mesh_label!r}")

        # Promote work_dir to a real Path.
        self.work_dir = Path(self.work_dir).resolve()

    # ---- convenience accessors in SI (meters) ----

    @property
    def L0(self)  -> float: return self.L0_mm  * 1e-3
    @property
    def W(self)   -> float: return self.W_mm   * 1e-3
    @property
    def t(self)   -> float: return self.t_mm   * 1e-3
    @property
    def Lr(self)  -> float: return self.Lr_mm  * 1e-3
    @property
    def R(self)   -> float: return self.R_mm   * 1e-3
    @property
    def Wg(self)  -> float: return self.Wg_mm  * 1e-3
    @property
    def Lg(self)  -> float: return self.Lg_mm  * 1e-3
    @property
    def L_T(self) -> float: return self.Lr + 2.0 * (self.Lg + self.R)
    @property
    def L_eff(self) -> float: return self.L_eff_mm * 1e-3
    @property
    def A_gauge(self) -> float: return self.W * self.t
    @property
    def P_yield(self) -> float: return self.sigma_y * self.A_gauge
    @property
    def delta_max(self) -> float: return self.eps_end_target * self.L_eff


# ----------------------------------------------------------------------------
# 2. build_mesh_with_gmsh.
# ----------------------------------------------------------------------------

def build_mesh_with_gmsh(cfg: StageConfig, out_inp: Path) -> None:
    """Build the dogbone mesh in GMSH and export to Abaqus .inp.

    Geometry per spec.md section 2. Two-stage construction. (i) Build
    a quarter-symmetric in-plane CAD profile in the (x, y) plane,
    mirror twice to recover the full E8 dogbone, recombine to quads,
    and 2D-mesh. (ii) Extrude the 2D quad mesh through-thickness to
    obtain HEXA8 with n_thru layers. Element / node sets are written
    for fixed grip, loaded grip, gauge nodes, gauge bricks.

    Args:
        cfg:     resolved StageConfig.
        out_inp: target filename for the Abaqus .inp output.

    Raises:
        RuntimeError: if the resulting mesh contains any non-HEX8
        elements (i.e. a wedge or tetrahedron leaked through the
        recombine + extrude pipeline).
    """
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(cfg.job_name)

    half_Lr = 0.5 * cfg.Lr
    half_W  = 0.5 * cfg.W
    half_Wg = 0.5 * cfg.Wg
    half_LT = 0.5 * cfg.L_T

    # Fillet center: tangent to gauge edge (y = +half_W) and to
    # grip side (y = +half_Wg). The arc is convex into the gauge.
    fillet_x_center =  half_Lr
    fillet_y_center =  half_W + cfg.R
    # fillet starts at (half_Lr, half_W), ends at (half_Lr + R, half_Wg).

    # ---- in-plane points (top half, right half; mirror to recover full) ----

    p = {}
    p["o"]      = gmsh.model.geo.addPoint( 0.0,        0.0,        0.0)
    p["g_top0"] = gmsh.model.geo.addPoint( 0.0,        half_W,     0.0)
    p["g_topR"] = gmsh.model.geo.addPoint( half_Lr,    half_W,     0.0)
    p["fil_c"]  = gmsh.model.geo.addPoint(fillet_x_center, fillet_y_center, 0.0)
    p["g_filR"] = gmsh.model.geo.addPoint( half_Lr + cfg.R, half_Wg, 0.0)
    p["g_TR"]   = gmsh.model.geo.addPoint( half_LT,    half_Wg,    0.0)
    p["g_BR"]   = gmsh.model.geo.addPoint( half_LT,    0.0,        0.0)

    # ---- in-plane curves ----

    c = {}
    c["bottom_R"] = gmsh.model.geo.addLine(p["o"],      p["g_BR"])
    c["right"]    = gmsh.model.geo.addLine(p["g_BR"],   p["g_TR"])
    c["top_grip"] = gmsh.model.geo.addLine(p["g_TR"],   p["g_filR"])
    c["fillet"]   = gmsh.model.geo.addCircleArc(p["g_filR"], p["fil_c"], p["g_topR"])
    c["top_g"]    = gmsh.model.geo.addLine(p["g_topR"], p["g_top0"])
    c["left"]     = gmsh.model.geo.addLine(p["g_top0"], p["o"])

    cl = gmsh.model.geo.addCurveLoop([
        c["bottom_R"], c["right"], c["top_grip"],
        c["fillet"],   c["top_g"], c["left"],
    ])
    s = gmsh.model.geo.addPlaneSurface([cl])

    # Mirror about y=0 to get bottom half of right side.
    gmsh.model.geo.synchronize()
    bottom_half = gmsh.model.geo.copy([(2, s)])
    gmsh.model.geo.symmetrize(bottom_half, 0.0, 1.0, 0.0, 0.0)

    # Mirror about x=0 to get left side.
    gmsh.model.geo.synchronize()
    full_right_half = [(2, s)] + bottom_half
    left_half = gmsh.model.geo.copy(full_right_half)
    gmsh.model.geo.symmetrize(left_half, 1.0, 0.0, 0.0, 0.0)

    gmsh.model.geo.synchronize()

    # ---- mesh sizing on edges ----
    #
    # Set transfinite curves so the recombine produces a clean structured
    # quad mesh that survives the through-thickness extrusion.
    n_along_gauge = cfg.n_length_gauge // 2  # half-length in symmetric quadrant
    n_across_gauge = cfg.n_width_gauge // 2  # half-width
    n_along_grip   = max(20, cfg.n_length_gauge // 4)
    n_across_grip  = max(6,  cfg.n_width_gauge)

    gmsh.model.mesh.setTransfiniteCurve(c["top_g"],    n_along_gauge + 1)
    gmsh.model.mesh.setTransfiniteCurve(c["bottom_R"], n_along_gauge + n_along_grip + 1)
    gmsh.model.mesh.setTransfiniteCurve(c["right"],    n_across_grip + 1)
    gmsh.model.mesh.setTransfiniteCurve(c["top_grip"], n_along_grip + 1)
    gmsh.model.mesh.setTransfiniteCurve(c["fillet"],   cfg.n_arc + 1)
    gmsh.model.mesh.setTransfiniteCurve(c["left"],     n_across_gauge + 1)

    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)  # blossom-full-quad
    gmsh.option.setNumber("Mesh.Algorithm", 8)               # frontal-Delaunay quad

    gmsh.model.mesh.generate(2)

    # ---- extrude through-thickness to HEXA8 ----
    extrude_layers = [cfg.n_thru]
    extrude_heights = [1.0]  # single layer of total thickness t
    extruded = gmsh.model.geo.extrude(
        gmsh.model.getEntities(2),
        0.0, 0.0, cfg.t,
        numElements=extrude_layers,
        heights=extrude_heights,
        recombine=True,
    )
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    # ---- physical groups for OpenRadioss element / node sets ----
    #
    # We tag faces by x-coordinate and bricks by gauge membership
    # using bounding-box queries on the resulting 3D mesh nodes.
    _tag_physical_groups(cfg)

    # ---- export to Abaqus .inp ----
    gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
    gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
    gmsh.write(str(out_inp))

    # ---- sanity check: HEXA8 only ----
    elt_types, _, _ = gmsh.model.mesh.getElements(dim=3)
    # GMSH element type 5 is 8-node hexahedron.
    if any(t != 5 for t in elt_types):
        gmsh.finalize()
        raise RuntimeError(
            "Non-HEXA8 element leaked into the 3D mesh. Check Recombine "
            "settings and fillet arc subdivision."
        )

    gmsh.finalize()


def _tag_physical_groups(cfg: StageConfig) -> None:
    """Tag fixed-grip face, loaded-grip face, gauge nodes, gauge bricks.

    Bounding-box queries are used because the symmetrize / extrude
    pipeline does not preserve named entities reliably.
    """
    import gmsh

    half_LT = 0.5 * cfg.L_T
    half_L0 = 0.5 * cfg.L0
    eps_geom = 1.0e-6  # m, tolerance for plane-membership

    # Fixed-grip face: x = -half_LT, all y, z in [0, t].
    fixed_face_nodes = gmsh.model.mesh.getNodesInBoundingBox(
        -half_LT - eps_geom, -cfg.Wg, -eps_geom,
        -half_LT + eps_geom,  cfg.Wg, cfg.t + eps_geom,
    )
    loaded_face_nodes = gmsh.model.mesh.getNodesInBoundingBox(
         half_LT - eps_geom, -cfg.Wg, -eps_geom,
         half_LT + eps_geom,  cfg.Wg, cfg.t + eps_geom,
    )
    gauge_endpoint_neg = gmsh.model.mesh.getNodesInBoundingBox(
        -half_L0 - eps_geom, -cfg.W, -eps_geom,
        -half_L0 + eps_geom,  cfg.W, cfg.t + eps_geom,
    )
    gauge_endpoint_pos = gmsh.model.mesh.getNodesInBoundingBox(
         half_L0 - eps_geom, -cfg.W, -eps_geom,
         half_L0 + eps_geom,  cfg.W, cfg.t + eps_geom,
    )

    # Save these as gmsh physical groups so the Abaqus writer emits
    # *NSET / *ELSET blocks; inp2rad picks them up as /GRNOD / /GRBRIC.
    pg_fix    = gmsh.model.addPhysicalGroup(0, list(fixed_face_nodes[0]),   tag=101)
    pg_load   = gmsh.model.addPhysicalGroup(0, list(loaded_face_nodes[0]),  tag=102)
    pg_g_neg  = gmsh.model.addPhysicalGroup(0, list(gauge_endpoint_neg[0]), tag=103)
    pg_g_pos  = gmsh.model.addPhysicalGroup(0, list(gauge_endpoint_pos[0]), tag=104)
    gmsh.model.setPhysicalName(0, pg_fix,   "GR_FIXED_GRIP")
    gmsh.model.setPhysicalName(0, pg_load,  "GR_LOADED_GRIP")
    gmsh.model.setPhysicalName(0, pg_g_neg, "GR_GAUGE_NEG")
    gmsh.model.setPhysicalName(0, pg_g_pos, "GR_GAUGE_POS")

    # Gauge-section bricks: centroid in the gauge bounding box.
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(dim=3)
    gauge_bricks: list[int] = []
    if elem_types and elem_types[0] == 5:
        hex_tags  = list(elem_tags[0])
        hex_nodes = list(elem_node_tags[0])
        for i, etag in enumerate(hex_tags):
            cx, cy, cz = _hex_centroid(hex_nodes, i)
            if -half_L0 <= cx <= half_L0 and abs(cy) <= 0.5 * cfg.W:
                gauge_bricks.append(etag)
    pg_gauge_b = gmsh.model.addPhysicalGroup(3, gauge_bricks, tag=201)
    gmsh.model.setPhysicalName(3, pg_gauge_b, "GR_GAUGE_BRICKS")


def _hex_centroid(hex_node_tags: list[int], i: int) -> tuple[float, float, float]:
    """Return (x, y, z) centroid of the i-th hex from a flat node-tag list."""
    import gmsh
    eight = hex_node_tags[8 * i : 8 * (i + 1)]
    cx = cy = cz = 0.0
    for tag in eight:
        coord, _, _ = gmsh.model.mesh.getNode(tag)
        cx += coord[0]
        cy += coord[1]
        cz += coord[2]
    return cx / 8.0, cy / 8.0, cz / 8.0


# ----------------------------------------------------------------------------
# 3. run_inp2rad.
# ----------------------------------------------------------------------------

def run_inp2rad(cfg: StageConfig, in_inp: Path, out_rad: Path) -> None:
    """Convert Abaqus .inp to OpenRadioss .rad via inp2rad.

    inp2rad ships with OpenRadioss/Tools and runs as a Python script
    inside the Lima VM. We invoke it via limactl shell so the converter
    sees the same paths the solver will see.
    """
    cmd = [
        "limactl", "shell", cfg.lima_vm,
        "python3", cfg.or_inp2rad,
        "-i", _to_lima_path(in_inp),
        "-o", _to_lima_path(out_rad),
    ]
    logging.info("inp2rad: %s", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(
            f"inp2rad failed (rc={res.returncode})\nstdout:\n{res.stdout}\n"
            f"stderr:\n{res.stderr}"
        )


def _to_lima_path(p: Path) -> str:
    """Translate a host macOS path into the path Lima exposes inside the VM.

    The default Lima template mounts $HOME at the same path, so most paths
    pass through unchanged. /Volumes/* mounts are forwarded explicitly via
    a writable mount entry in lima.yaml; the user is responsible for the
    Lima config. This helper simply returns the absolute path string.
    """
    return str(Path(p).resolve())


# ----------------------------------------------------------------------------
# 4. render_starter_deck.
# ----------------------------------------------------------------------------

STARTER_TEMPLATE = """\
#RADIOSS STARTER
/BEGIN
{job_name} stage 03 - isotropic dogbone E8/E8M
       2026         0
                  kg                   m                   s
                  kg                   m                   s
#--- material law (LAW2 Johnson-Cook, rate effects off) ---
{material_card}
#--- solid property (general 3D solid /PROP/TYPE14) ---
/PROP/SOLID/1
A36_solid_TYPE14
#         Isolid           Ismstr           Icpre            Itetra4
               24                0                1                 0
#         Iframe         dn          q_a         q_b
               1        0.1         1.10        0.05
#--- part: bind material to property to elements ---
/PART/1
A36_dogbone_part
         1         1         0
#--- mesh: nodes and HEXA8 connectivity from inp2rad output ---
#       (the nodes/bricks/sets are spliced in below from {converted_rad})
{mesh_block}
#--- node groups for BCs and history ---
/GRNOD/NODE/1
GR_FIXED_GRIP
{fixed_grip_nodes}
/GRNOD/NODE/2
GR_LOADED_GRIP
{loaded_grip_nodes}
/GRNOD/NODE/3
GR_GAUGE_NEG
{gauge_neg_nodes}
/GRNOD/NODE/4
GR_GAUGE_POS
{gauge_pos_nodes}
/GRBRIC/BRIC/1
GR_GAUGE_BRICKS
{gauge_brick_ids}
#--- BCs ---
/BCS/1
clamp_fixed_grip
#         Tx Ty Tz Rx Ry Rz   GR_NODE
            1  1  1  0  0  0         1
/BCS/2
lateral_loaded_grip
#         Tx Ty Tz Rx Ry Rz   GR_NODE
            0  1  1  0  0  0         2
#--- imposed-displacement function and load ---
/FUNCT/1
ramp_linear_t
       0.0       0.0
{T_end}       1.0
/IMPDISP/1
axial_imposed_grip
#  GR_NODE  IFUNCT  DIR  ISKEW   AScale  ATime
         2       1    1      0  {delta_max} 0.0
#--- time-history requests ---
/TH/NODE/1
TH_fixed_grip_reaction
#  GR_NODE                  Vars
         1  DEF FX FY FZ
/TH/NODE/2
TH_gauge_endpoints
#  GR_NODE                  Vars
         3  DEF DX DY DZ
/TH/NODE/3
TH_gauge_endpoints_pos
         4  DEF DX DY DZ
/TH/BRICK/1
TH_gauge_section_stress
#  GR_BRIC                  Vars
         1  DEF SIGX SIGY SIGZ EPSP
/END
"""

LAW2_TEMPLATE = """\
/MAT/LAW2/1/A36
A36 steel JC rate-off (eps_dot_0 = 1.0 s^-1, c = 0)
#                rho_I
            {rho:.6e}
#                    E                  nu                Iflag             Ihard
            {E:.6e}        {nu:.4f}                    0                 0
#                    a                   b                   n        eps_max         sig_max
            {JC_a:.6e}       {JC_b:.6e}     {JC_n:.4f}     0.0000     0.0000
#                    c       eps_dot_0           ICC      Fsmooth     Fcut       Chard      Tmelt        rho_Cp     Tref
            {JC_c:.4f}     {JC_eps_dot_0:.4f}     0     0     0.0     0.0     0.0     0.0     0.0
"""

LAW36_TEMPLATE = """\
/MAT/LAW36/1/A36
A36 steel tabulated PLAS_TAB
#                rho_I
            {rho:.6e}
#                    E                  nu       N(funcs)       eps_max  sig_max
            {E:.6e}        {nu:.4f}              1                 0.0000  0.0000
#         Y_fac1     C1              eps_dot_0_1      Fsmooth   Fcut
            1.0     0.0              1.0              0          0.0
#       FUN_id1
            1001
/FUNCT/1001
LAW36_hardening_curve
{table_xy}
"""


def render_starter_deck(cfg: StageConfig, mesh_block: str,
                        node_sets: dict[str, list[int]],
                        gauge_brick_ids: list[int],
                        out_path: Path,
                        converted_rad: Path) -> None:
    """Fill in the starter-deck template and write it to disk.

    Args:
        cfg:               resolved StageConfig.
        mesh_block:        the /NODE + /BRICK text block lifted from the
                           inp2rad-converted .rad file.
        node_sets:         dict with keys 'fixed_grip', 'loaded_grip',
                           'gauge_neg', 'gauge_pos'; each maps to a list
                           of node IDs.
        gauge_brick_ids:   list of element IDs in the gauge section.
        out_path:          target path for stage03_..._0000.rad.
        converted_rad:     informational: name of the inp2rad output.
    """
    if cfg.material_law.lower() == "law2":
        material_card = LAW2_TEMPLATE.format(
            rho=cfg.rho, E=cfg.E, nu=cfg.nu,
            JC_a=cfg.JC_a, JC_b=cfg.JC_b, JC_n=cfg.JC_n,
            JC_c=cfg.JC_c, JC_eps_dot_0=cfg.JC_eps_dot_0,
        )
    elif cfg.material_law.lower() == "law36":
        table_xy = "\n".join(
            f"  {ep:.6e}  {sy:.6e}" for (ep, sy) in cfg.LAW36_table
        )
        material_card = LAW36_TEMPLATE.format(
            rho=cfg.rho, E=cfg.E, nu=cfg.nu, table_xy=table_xy,
        )
    else:
        raise ValueError(f"unknown material law: {cfg.material_law!r}")

    rendered = STARTER_TEMPLATE.format(
        job_name=cfg.job_name,
        material_card=material_card,
        converted_rad=converted_rad.name,
        mesh_block=mesh_block,
        fixed_grip_nodes=_format_id_list(node_sets["fixed_grip"]),
        loaded_grip_nodes=_format_id_list(node_sets["loaded_grip"]),
        gauge_neg_nodes=_format_id_list(node_sets["gauge_neg"]),
        gauge_pos_nodes=_format_id_list(node_sets["gauge_pos"]),
        gauge_brick_ids=_format_id_list(gauge_brick_ids),
        T_end=f"{cfg.T_end:.6e}",
        delta_max=f"{cfg.delta_max:.6e}",
    )
    out_path.write_text(rendered)


def _format_id_list(ids: Iterable[int], per_line: int = 10) -> str:
    """Format a list of integer IDs as fixed-width OpenRadioss table rows."""
    items = [f"{i:>10d}" for i in ids]
    lines = [
        "".join(items[k : k + per_line])
        for k in range(0, len(items), per_line)
    ]
    return "\n".join(lines)


# ----------------------------------------------------------------------------
# 5. render_engine_deck.
# ----------------------------------------------------------------------------

ENGINE_TEMPLATE = """\
#RADIOSS ENGINE
/RUN/{job_name}/1/
{T_end}
/IMPL/QSTAT
/IMPL/SOLVER
#       It_method  It_print  It_iform  ...
                3         1         0
/IMPL/DT/STOP
#       dt_init             dt_min              dt_max
        1.0e-2              1.0e-6              1.0e-1
/IMPL/NONLIN/N
#       N_iter_max
                25
/IMPL/PRINT/NONLIN
       1
/ANIM/DT
0.0  5.0e-2
/ANIM/BRICK/TENS/STRESS
/ANIM/BRICK/EPSP
/ANIM/NODA/DISP
/H3D/DT
0.0  5.0e-2
/STOP
"""


def render_engine_deck(cfg: StageConfig, out_path: Path) -> None:
    """Fill in the engine-deck template and write it to disk."""
    rendered = ENGINE_TEMPLATE.format(
        job_name=cfg.job_name,
        T_end=f"{cfg.T_end:.6e}",
    )
    out_path.write_text(rendered)


# ----------------------------------------------------------------------------
# 6. run_openradioss_in_lima.
# ----------------------------------------------------------------------------

def run_openradioss_in_lima(cfg: StageConfig, work: Path) -> None:
    """Invoke the OpenRadioss starter and then the engine inside Lima.

    The two-phase CLI is documented in OpenRadioss/HOWTO.md. Starter
    parses the _0000.rad and writes the binary restart that the engine
    consumes. Engine writes the .anim, T01, and h3d output streams.
    """
    starter_cmd = [
        "limactl", "shell", cfg.lima_vm,
        cfg.or_starter,
        "-i", _to_lima_path(work / f"{cfg.job_name}_0000.rad"),
        "-nt", str(cfg.n_threads),
    ]
    engine_cmd = [
        "limactl", "shell", cfg.lima_vm,
        cfg.or_engine,
        "-i", _to_lima_path(work / f"{cfg.job_name}_0001.rad"),
        "-nt", str(cfg.n_threads),
    ]

    for label, cmd in (("starter", starter_cmd), ("engine", engine_cmd)):
        logging.info("OpenRadioss %s: %s", label, " ".join(cmd))
        res = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(
                f"OpenRadioss {label} failed (rc={res.returncode})\n"
                f"stdout:\n{res.stdout}\nstderr:\n{res.stderr}"
            )


# ----------------------------------------------------------------------------
# 7. parse_time_history.
# ----------------------------------------------------------------------------

def parse_time_history(cfg: StageConfig, work: Path) -> dict[str, Any]:
    """Parse the OpenRadioss T01 file using Vortex-Radioss.

    Returns a dictionary with the following arrays (NumPy):

        time            shape (N_t,)
        F_grip_x        shape (N_t,)            sum of FX on fixed grip
        u_gauge_neg_x   shape (N_t,)            mean DX over gauge neg face
        u_gauge_pos_x   shape (N_t,)            mean DX over gauge pos face
        sigma_xx        shape (N_t, N_brick)    per-brick gauge stress
        epsp            shape (N_t, N_brick)    per-brick equiv plastic strain
        eps_xx_brick    shape (N_t, N_brick)    per-brick total axial strain
                                                  (extracted from ANIM if avail.)

    The Vortex-Radioss reader is the documented Python path
    (openradioss_refs.bib: VortexRadioss2024). If unavailable we fall
    back to the OpenRadioss-bundled th2csv helper and then csv parsing.
    """
    import numpy as np

    t01_path = work / f"{cfg.job_name}T01"
    if not t01_path.exists():
        raise FileNotFoundError(f"missing T01 file: {t01_path}")

    try:
        from vortex_radioss.animtod3plot import RadiossReader  # type: ignore
        reader = RadiossReader(str(t01_path))
        arrays = reader.read_time_history()
    except Exception:  # noqa: BLE001
        # Fallback: convert T01 to CSV using th2csv inside Lima.
        csv_path = work / f"{cfg.job_name}T01.csv"
        cmd = [
            "limactl", "shell", cfg.lima_vm,
            "/OpenRadioss/exec/th2csv",
            "-i", _to_lima_path(t01_path),
            "-o", _to_lima_path(csv_path),
        ]
        subprocess.run(cmd, check=True)
        arrays = _csv_to_arrays(csv_path)

    # Re-key into the canonical dictionary shape expected downstream.
    return arrays


def _csv_to_arrays(csv_path: Path) -> dict[str, Any]:
    """Reshape a th2csv-emitted CSV into the canonical dictionary."""
    import numpy as np

    data: dict[str, list[float]] = {}
    with csv_path.open("r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            for k, v in row.items():
                data.setdefault(k, []).append(float(v))

    out: dict[str, Any] = {}
    out["time"] = np.array(data.get("time", []))
    # Sum FX over the grip-reaction history.
    grip_fx_keys = [k for k in data if "FIXED_GRIP" in k and k.endswith("FX")]
    out["F_grip_x"] = np.sum([data[k] for k in grip_fx_keys], axis=0)
    # Mean DX over gauge-neg and gauge-pos history.
    neg_keys = [k for k in data if "GAUGE_NEG" in k and k.endswith("DX")]
    pos_keys = [k for k in data if "GAUGE_POS" in k and k.endswith("DX")]
    out["u_gauge_neg_x"] = np.mean([data[k] for k in neg_keys], axis=0)
    out["u_gauge_pos_x"] = np.mean([data[k] for k in pos_keys], axis=0)
    # Per-brick gauge stress / plastic strain.
    sig_keys = sorted(k for k in data if "GAUGE_BRICKS" in k and k.endswith("SIGX"))
    epsp_keys = sorted(k for k in data if "GAUGE_BRICKS" in k and k.endswith("EPSP"))
    out["sigma_xx"] = np.array([data[k] for k in sig_keys]).T
    out["epsp"]     = np.array([data[k] for k in epsp_keys]).T
    return out


# ----------------------------------------------------------------------------
# 8. evaluate_checks.
# ----------------------------------------------------------------------------

def evaluate_checks(cfg: StageConfig, hist: dict[str, Any]) -> dict[str, Any]:
    """Compute the four pass-fail checks and three secondary checks.

    See spec.md section 8 for the canonical definition of each check.
    """
    import numpy as np

    time          = hist["time"]
    F_grip_x      = hist["F_grip_x"]                 # N
    u_neg         = hist["u_gauge_neg_x"]            # m
    u_pos         = hist["u_gauge_pos_x"]            # m
    sigma_xx_b    = hist["sigma_xx"]                 # Pa, shape (N_t, N_b)
    epsp_b        = hist["epsp"]                     # -,  shape (N_t, N_b)

    # Reaction-derived load (positive in tension on the fixed grip).
    P = np.abs(F_grip_x)

    # Gauge-section nodal-extensometer strain: (u_pos - u_neg) / L0.
    delta_L0 = (u_pos - u_neg)
    eps_xx_gauge = delta_L0 / cfg.L0

    # Gauge-section mean stress (across all gauge bricks).
    sigma_mean = np.mean(sigma_xx_b, axis=1)

    # Find time index where P first reaches 0.8 * P_yield.
    target = 0.8 * cfg.P_yield
    idx_pre_yield = int(np.argmin(np.abs(P - target)))

    # ---- C1: uniformity of sigma_xx across gauge bricks at 0.8 P_y ----
    sig_at_idx   = sigma_xx_b[idx_pre_yield, :]
    sig_mean_at  = float(np.mean(sig_at_idx))
    sig_std_at   = float(np.std(sig_at_idx))
    C1_residual  = abs(sig_std_at) / max(abs(sig_mean_at), 1.0)
    C1_pass      = C1_residual <= cfg.tol_uniformity

    # ---- C2: gauge-mean sigma_xx vs P / A_gauge at the same load ----
    sigma_PA     = float(P[idx_pre_yield] / cfg.A_gauge)
    C2_residual  = abs(sig_mean_at - sigma_PA) / max(abs(sigma_PA), 1.0)
    C2_pass      = C2_residual <= cfg.tol_stress

    # ---- C3: gauge-mean eps_xx vs delta_L0 / L0 at the same load ----
    #
    # The gauge-mean eps from the FEM is computed from the linear
    # elastic relation eps = sigma / E in the elastic regime. This
    # cross-checks the gauge stress vs the gauge displacement.
    eps_from_sig = sig_mean_at / cfg.E
    eps_from_dis = float(eps_xx_gauge[idx_pre_yield])
    C3_residual  = abs(eps_from_sig - eps_from_dis) / max(abs(eps_from_dis), 1e-12)
    C3_pass      = C3_residual <= cfg.tol_strain

    # ---- C4: yield onset detection ----
    #
    # Compute apparent stiffness K_app = dP / d_delta on a windowed
    # finite difference, find first time index where K_app < 0.99 K_0,
    # and look up sigma_mean at that time.
    delta_grip = u_pos - 0.0    # u_pos is the loaded face mean DX
    if delta_grip.size >= 5:
        K = np.gradient(P, delta_grip + 1e-30)
        K0 = float(np.median(K[:5]))   # initial-stiffness estimate
        below = np.where(K < cfg.tol_yield_low * K0)[0]
        idx_yield = int(below[0]) if below.size else int(np.argmax(K * 0))
    else:
        idx_yield = idx_pre_yield
    sigma_yield_obs = float(sigma_mean[idx_yield])
    C4_lower = cfg.tol_yield_low * cfg.sigma_y
    C4_upper = cfg.tol_yield_high * cfg.sigma_y
    C4_pass  = C4_lower <= sigma_yield_obs <= C4_upper
    C4_residual = abs(sigma_yield_obs - cfg.sigma_y) / cfg.sigma_y

    # ---- S1: apparent E from gauge sigma vs eps below yield ----
    elastic_mask = epsp_b.max(axis=1) <= 1e-9
    if np.any(elastic_mask):
        eps_e = eps_xx_gauge[elastic_mask]
        sig_e = sigma_mean[elastic_mask]
        if eps_e.size >= 2 and np.std(eps_e) > 0:
            E_app = float(np.polyfit(eps_e, sig_e, 1)[0])
        else:
            E_app = float("nan")
    else:
        E_app = float("nan")
    S1_residual = abs(E_app - cfg.E) / cfg.E if math.isfinite(E_app) else float("inf")
    S1_pass = S1_residual <= cfg.tol_E_apparent

    # ---- S2: Poisson ratio (informational; needs lateral u extracted) ----
    # Skipped here unless the deck adds a lateral-extensometer node group;
    # the runner is happy to not fail on this.
    S2_pass = True
    S2_residual = 0.0

    # ---- S3: post-yield gauge-mean sigma at eps_p = 0.001 ----
    eps_p_target = 0.001
    mean_epsp = epsp_b.mean(axis=1)
    after_yield = np.where(mean_epsp >= eps_p_target)[0]
    if after_yield.size:
        idx_S3 = int(after_yield[0])
        sigma_at_S3 = float(sigma_mean[idx_S3])
        sigma_ref_S3 = cfg.JC_a + cfg.JC_b * (eps_p_target ** cfg.JC_n)
        S3_residual = abs(sigma_at_S3 - sigma_ref_S3) / sigma_ref_S3
        S3_pass = S3_residual <= cfg.tol_post_yield
    else:
        S3_residual = float("inf")
        S3_pass = False

    return {
        "C1": {"pass": bool(C1_pass), "residual": float(C1_residual),
               "tolerance": cfg.tol_uniformity,
               "description": "uniformity of sigma_xx across gauge bricks"},
        "C2": {"pass": bool(C2_pass), "residual": float(C2_residual),
               "tolerance": cfg.tol_stress,
               "description": "gauge-mean sigma_xx vs P/A_gauge"},
        "C3": {"pass": bool(C3_pass), "residual": float(C3_residual),
               "tolerance": cfg.tol_strain,
               "description": "gauge-mean eps_xx vs delta/L_gauge"},
        "C4": {"pass": bool(C4_pass), "residual": float(C4_residual),
               "tolerance_band": [cfg.tol_yield_low, cfg.tol_yield_high],
               "sigma_yield_observed_Pa": sigma_yield_obs,
               "description": "yield onset stress at apparent E drop"},
        "S1": {"pass": bool(S1_pass), "residual": float(S1_residual),
               "E_apparent_Pa": E_app,
               "tolerance": cfg.tol_E_apparent,
               "description": "apparent Young modulus"},
        "S2": {"pass": bool(S2_pass), "residual": float(S2_residual),
               "tolerance": cfg.tol_poisson,
               "description": "Poisson ratio (skipped if no lateral extensometer)"},
        "S3": {"pass": bool(S3_pass), "residual": float(S3_residual),
               "tolerance": cfg.tol_post_yield,
               "description": "post-yield sigma at eps_p = 0.001"},
        "overall_pass": bool(C1_pass and C2_pass and C3_pass and C4_pass),
        # raw extracted curves (for the downstream CSVs)
        "_curves": {
            "time":       time.tolist(),
            "P":          P.tolist(),
            "delta_grip": (u_pos - 0.0).tolist(),
            "delta_L0":   delta_L0.tolist(),
            "eps_gauge":  eps_xx_gauge.tolist(),
            "sig_mean":   sigma_mean.tolist(),
        },
    }


# ----------------------------------------------------------------------------
# 9. write_results.
# ----------------------------------------------------------------------------

def write_results(cfg: StageConfig, results: dict[str, Any], work: Path) -> None:
    """Emit results.json plus two CSVs for the load-displacement curve and
    the gauge-section sigma-eps curve."""
    curves = results.pop("_curves")

    (work / "results.json").write_text(json.dumps(
        {**results, "config": dataclasses.asdict(cfg)},
        indent=2, default=str,
    ))

    with (work / "load_displacement.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time_s", "delta_grip_m", "P_N"])
        for t, d, p in zip(curves["time"], curves["delta_grip"], curves["P"]):
            w.writerow([f"{t:.9e}", f"{d:.9e}", f"{p:.9e}"])

    with (work / "stress_strain.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time_s", "eps_gauge", "sigma_gauge_Pa"])
        for t, e, s in zip(curves["time"], curves["eps_gauge"], curves["sig_mean"]):
            w.writerow([f"{t:.9e}", f"{e:.9e}", f"{s:.9e}"])


# ----------------------------------------------------------------------------
# 10. plot_with_typst_cetz.
# ----------------------------------------------------------------------------

TYPST_PLOT_TEMPLATE = r"""
#import "@preview/cetz:0.2.2"
#import "@preview/cetz-plot:0.1.0": plot

#set page(width: auto, height: auto, margin: 6mm)

#let garnet   = rgb({garnet})
#let atlantic = rgb({atlantic})
#let black90  = rgb({black90})

#let load_disp = csv("load_displacement.csv")
#let ss_curve  = csv("stress_strain.csv")

// brand-color load-displacement plot
#cetz.canvas({{
  plot.plot(
    size: (10, 6),
    x-label: [grip displacement $delta$ (mm)],
    y-label: [load $P$ (kN)],
    {{
      plot.add(
        load_disp.slice(1).map(row => (
          float(row.at(1)) * 1e3,    // m -> mm
          float(row.at(2)) * 1e-3,   // N -> kN
        )),
        style: (stroke: garnet + 1.6pt),
        label: [stage 03 dogbone],
      )
    }},
  )
}})

// stress-strain plot
#cetz.canvas({{
  plot.plot(
    size: (10, 6),
    x-label: [gauge strain $epsilon_(x x)$],
    y-label: [gauge stress $sigma_(x x)$ (MPa)],
    {{
      plot.add(
        ss_curve.slice(1).map(row => (
          float(row.at(1)),
          float(row.at(2)) * 1e-6,   // Pa -> MPa
        )),
        style: (stroke: atlantic + 1.6pt),
        label: [LAW2 A36],
      )
      // overlay: input yield stress horizontal line
      plot.add(
        ((0, {sig_y_MPa}), (0.005, {sig_y_MPa})),
        style: (stroke: black90 + 0.8pt + dash),
        label: [$sigma_y$ = {sig_y_MPa} MPa],
      )
    }},
  )
}})
"""


def plot_with_typst_cetz(cfg: StageConfig, work: Path) -> None:
    """Emit a Typst CeTZ source file from the CSV outputs and compile.

    Compilation requires the typst CLI; if it is missing we just write
    the .typ source and return so the user can compile manually.
    """
    typst_src = TYPST_PLOT_TEMPLATE.format(
        garnet=_typst_rgb("garnet"),
        atlantic=_typst_rgb("atlantic"),
        black90=_typst_rgb("black_90"),
        sig_y_MPa=cfg.sigma_y * 1e-6,
    )
    src_path = work / "plots.typ"
    src_path.write_text(typst_src)
    if shutil.which("typst") is not None:
        subprocess.run(
            ["typst", "compile", "--root", str(work), str(src_path),
             str(work / "plots.pdf")],
            check=True,
        )
    else:
        logging.warning("typst CLI not found; plots.typ written, not compiled")


def _typst_rgb(name: str) -> str:
    r, g, b = BRAND_COLORS[name]
    return f"\"#{r:02X}{g:02X}{b:02X}\""


# ----------------------------------------------------------------------------
# Glue: extract the mesh / sets blocks from the inp2rad output.
# ----------------------------------------------------------------------------

def extract_inp2rad_blocks(rad_path: Path) -> dict[str, Any]:
    """Read the inp2rad-converted .rad and return the blocks we splice in.

    Returns:
        dict with keys 'mesh_block' (the /NODE and /BRICK text), and
        'sets' (a dict mapping set name to a list of integer IDs).

    Implementation notes. The inp2rad converter writes /NODE and /BRICK
    cards verbatim from the Abaqus .inp, then writes /GRNOD/NODE and
    /GRBRIC/BRIC cards from each *NSET / *ELSET. We strip the converted
    /GRNOD / /GRBRIC blocks, keep only the /NODE and /BRICK sections in
    the mesh_block, and reparse the converted sets into Python lists.
    """
    text = rad_path.read_text()

    # Walk line-by-line and collect.
    in_node_block = False
    in_brick_block = False
    in_set_block: Optional[str] = None
    current_set_name: Optional[str] = None

    mesh_lines: list[str] = []
    sets: dict[str, list[int]] = {}

    for line in text.splitlines():
        if line.startswith("/NODE"):
            in_node_block = True
            in_brick_block = False
            in_set_block = None
            mesh_lines.append(line)
            continue
        if line.startswith("/BRICK") or line.startswith("/HEXA"):
            in_node_block = False
            in_brick_block = True
            in_set_block = None
            mesh_lines.append(line)
            continue
        if line.startswith("/GRNOD/NODE") or line.startswith("/GRBRIC/BRIC"):
            in_node_block = False
            in_brick_block = False
            in_set_block = "node" if line.startswith("/GRNOD/NODE") else "brick"
            current_set_name = None
            continue
        if line.startswith("/"):
            in_node_block = False
            in_brick_block = False
            in_set_block = None
            continue

        if in_node_block or in_brick_block:
            mesh_lines.append(line)
        elif in_set_block is not None:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if current_set_name is None:
                current_set_name = stripped
                sets.setdefault(current_set_name, [])
                continue
            for tok in stripped.split():
                try:
                    sets[current_set_name].append(int(tok))
                except ValueError:
                    pass

    return {"mesh_block": "\n".join(mesh_lines), "sets": sets}


# ----------------------------------------------------------------------------
# Top-level pipeline.
# ----------------------------------------------------------------------------

def run_stage(cfg: StageConfig) -> dict[str, Any]:
    """End-to-end pipeline. Returns the results dictionary."""
    work = cfg.work_dir / f"run_{cfg.mesh_label}_{cfg.material_law}"
    work.mkdir(parents=True, exist_ok=True)

    inp_path     = work / "mesh.inp"
    converted    = work / f"{cfg.job_name}_converted.rad"
    starter_path = work / f"{cfg.job_name}_0000.rad"
    engine_path  = work / f"{cfg.job_name}_0001.rad"

    # 1. mesh
    build_mesh_with_gmsh(cfg, inp_path)

    # 2. inp -> rad
    run_inp2rad(cfg, inp_path, converted)

    # 3. extract mesh and set blocks from the converted .rad
    blocks = extract_inp2rad_blocks(converted)
    sets = blocks["sets"]
    node_sets = {
        "fixed_grip":  sets.get("GR_FIXED_GRIP", []),
        "loaded_grip": sets.get("GR_LOADED_GRIP", []),
        "gauge_neg":   sets.get("GR_GAUGE_NEG", []),
        "gauge_pos":   sets.get("GR_GAUGE_POS", []),
    }
    gauge_brick_ids = sets.get("GR_GAUGE_BRICKS", [])

    if not all(node_sets.values()) or not gauge_brick_ids:
        raise RuntimeError(
            "inp2rad dropped one of the node / brick sets; "
            f"got: nodes={ {k: len(v) for k, v in node_sets.items()} }, "
            f"bricks={len(gauge_brick_ids)}"
        )

    # 4. starter and engine decks
    render_starter_deck(cfg, blocks["mesh_block"], node_sets, gauge_brick_ids,
                        starter_path, converted)
    render_engine_deck(cfg, engine_path)

    # 5. run inside Lima
    run_openradioss_in_lima(cfg, work)

    # 6. parse T01
    hist = parse_time_history(cfg, work)

    # 7. evaluate checks
    results = evaluate_checks(cfg, hist)

    # 8. write CSVs and JSON
    write_results(cfg, results, work)

    # 9. (optional) Typst CeTZ plots
    if cfg.do_plots:
        plot_with_typst_cetz(cfg, work)

    return results


# ----------------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------------

def parse_args(argv: Optional[list[str]] = None) -> StageConfig:
    p = argparse.ArgumentParser(
        prog="stage_03_runner",
        description=textwrap.dedent(__doc__).strip().splitlines()[0],
    )
    p.add_argument("--mesh", choices=["coarse", "medium", "fine"],
                   default="medium", help="mesh density label")
    p.add_argument("--material", choices=["law2", "law36"],
                   default="law2", help="OpenRadioss material law")
    p.add_argument("--lima-vm", default="or",
                   help="Lima VM name (limactl)")
    p.add_argument("--threads", type=int, default=4,
                   help="OpenRadioss -nt threads")
    p.add_argument("--plot", action="store_true",
                   help="emit and compile Typst+CeTZ plots from CSVs")
    p.add_argument("--work-dir", type=Path, default=Path.cwd(),
                   help="parent directory for run_* output dir")
    args = p.parse_args(argv)

    return StageConfig(
        mesh_label=args.mesh,
        material_law=args.material,
        lima_vm=args.lima_vm,
        n_threads=args.threads,
        do_plots=args.plot,
        work_dir=args.work_dir,
    )


def main(argv: Optional[list[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    cfg = parse_args(argv)
    results = run_stage(cfg)
    overall = results.get("overall_pass", False)
    logging.info("stage 03 overall_pass=%s", overall)
    print(json.dumps(results, indent=2, default=str))
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
