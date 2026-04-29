#!/usr/bin/env python3
"""Generate and run OpenRadioss LAW/PROP/FAIL starter compatibility probes."""

from __future__ import annotations

import csv
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
THIS = Path(__file__).resolve().parent
OR_ROOT = Path(os.environ.get("OR", "/mnt/storage/j-vaught/openradioss/OpenRadioss")).resolve()
STARTER = OR_ROOT / "exec" / "starter_linux64_gf"


def fmt_f(*values: float) -> str:
    return "".join(f"{value:20.12g}" for value in values)


def fmt_i(*values: int) -> str:
    return "".join(f"{value:10d}" for value in values)


def env() -> dict[str, str]:
    out = os.environ.copy()
    out["OR"] = str(OR_ROOT)
    out["RAD_CFG_PATH"] = str(OR_ROOT / "hm_cfg_files")
    out["RAD_H3D_PATH"] = str(OR_ROOT / "extlib" / "h3d" / "lib" / "linux64")
    reader = str(OR_ROOT / "extlib" / "hm_reader" / "linux64")
    out["LD_LIBRARY_PATH"] = reader + ":" + out.get("LD_LIBRARY_PATH", "")
    return out


@dataclass(frozen=True)
class Law:
    law: str
    name: str
    note: str
    cards: tuple[str, ...]


def law_cards() -> dict[str, Law]:
    common_comp = {
        "rho": 1.55e-9,
        "e1": 135000.0,
        "e2": 10000.0,
        "e3": 10000.0,
        "nu12": 0.30,
        "nu23": 0.40,
        "nu31": 0.0222222222,
        "g12": 5000.0,
        "g23": 3800.0,
        "g31": 5000.0,
        "xt": 1500.0,
        "xc": 1200.0,
        "yt": 50.0,
        "yc": 200.0,
        "zt": 50.0,
        "zc": 200.0,
        "s12": 90.0,
        "s23": 70.0,
        "s13": 90.0,
    }
    c = common_comp
    law12 = (
        "/MAT/LAW12/1",
        "probe_3d_comp",
        fmt_f(c["rho"]),
        fmt_f(c["e1"], c["e2"], c["e3"]),
        fmt_f(c["nu12"], c["nu23"], c["nu31"]),
        fmt_f(c["g12"], c["g23"], c["g31"]),
        fmt_f(c["xt"], c["yt"], c["zt"], 0.0),
        fmt_f(1.0, 1.0, 1.0),
        fmt_f(c["xt"], c["yt"], c["xc"], c["yc"]),
        fmt_f(c["s12"], c["s12"], c["s23"], c["s23"]),
        fmt_f(c["zt"], c["zc"], c["s13"], c["s13"]),
        fmt_f(0.0, 0.0, 0.0, 0.0) + fmt_i(0),
    )
    law14 = (
        "/MAT/LAW14/1",
        "probe_compso",
        fmt_f(c["rho"]),
        fmt_f(c["e1"], c["e2"], c["e3"]),
        fmt_f(c["nu12"], c["nu23"], c["nu31"]),
        fmt_f(c["g12"], c["g23"], c["g31"]),
        fmt_f(c["xt"], c["yt"], c["zt"], 0.0),
        fmt_f(1.0, 1.0, 1.0),
        fmt_f(c["xt"], c["yt"], c["xc"], c["yc"]),
        fmt_f(c["s12"], c["s12"], c["s23"], c["s23"]),
        fmt_f(0.0, 0.0, 0.0, 0.0) + fmt_i(0),
    )
    law15 = (
        "/MAT/LAW15/1",
        "probe_chang",
        fmt_f(c["rho"]),
        fmt_f(c["e1"], c["e2"], c["nu12"]),
        fmt_f(c["g12"], c["g23"], c["g31"]),
        fmt_f(1.0, 1.0, 1.0e30),
        fmt_f(1.0e30, 1.0) + fmt_i(0),
        fmt_f(c["xt"], c["yt"], c["xc"], c["yc"], 1.0),
        fmt_f(c["s12"], c["s12"], 0.0, 0.0) + fmt_i(0),
        fmt_f(1.0, 1.0e30, c["xt"], c["yt"], c["s12"]),
        fmt_i(0) + fmt_f(0.0, c["xc"], c["yc"]),
    )
    law25 = (
        "/MAT/LAW25/1",
        "probe_compsh_crasurv",
        fmt_f(c["rho"]),
        fmt_f(c["e1"], c["e2"], c["nu12"]) + fmt_i(0) + f"{c['e3']:20.12g}",
        fmt_f(c["g12"], c["g23"], c["g31"], 0.0, 0.0),
        fmt_f(0.0, 0.0, 0.0, 0.0, 1.0),
        fmt_f(0.0, 0.0) + fmt_i(0) + f"{0.0:20.12g}",
        fmt_f(0.0, 0.0, 0.0),
        fmt_f(c["xt"], c["yt"], c["xc"], c["yc"], 0.0),
        fmt_f(c["s12"], c["s12"], 0.0, 0.0) + fmt_i(0),
        fmt_f(0.0, 0.0, 0.0),
        fmt_i(0) + f"{0.0:20.12g}",
    )
    law28 = (
        "/MAT/LAW28/1",
        "probe_honeycomb",
        fmt_f(7.8e-5),
        fmt_f(6.429, 6.429, 6.429),
        fmt_f(2.5, 2.5, 2.5),
        fmt_i(1, 1, 1, 0) + fmt_f(0.0, 0.0, 0.0),
        fmt_f(0.0, 0.0, 0.0),
        fmt_i(1, 1, 1, 0) + fmt_f(0.0, 0.0, 0.0),
        fmt_f(0.0, 0.0, 0.0),
    )
    law42 = (
        "/MAT/LAW42/1",
        "probe_ogden",
        fmt_f(940.0),
        fmt_f(0.495, 0.0, 0.0, 0.0) + fmt_i(0, 0),
        fmt_f(827.6, 55.0, 0.0, 0.0, 0.0),
        "",
        fmt_f(2.0, -2.0, 0.0, 0.0, 0.0),
        "",
    )
    law53 = (
        "/MAT/LAW53/1",
        "probe_tsai_tab",
        fmt_f(c["rho"]),
        fmt_f(c["e1"], c["e2"]),
        fmt_f(c["g12"], c["g23"]),
        fmt_i(1, 1, 1, 1, 1),
        fmt_f(1.0, 1.0, 1.0, 1.0, 1.0),
    )
    law58 = (
        "/MAT/LAW58/1",
        "probe_fabric_a",
        fmt_f(c["rho"]),
        fmt_f(c["e1"], 0.0, c["e2"], 0.0, 1.0),
        fmt_f(c["g12"], c["g12"], 0.0),
        fmt_f(0.0, 0.0),
        fmt_i(1, 1) + fmt_f(1.0, 1.0),
    )
    law62 = (
        "/MAT/LAW62/1",
        "probe_visc_hyp",
        fmt_f(5.5e-8),
        fmt_f(0.0) + fmt_i(2, 4),
        fmt_f(6.9119e-5, 1.20126e-7),
        fmt_f(30.4979, -3.6961),
        fmt_f(0.34391, 0.10785, 0.070763, 0.077511),
        fmt_f(461.99, 14206.0, 199183.3, 8625495.0),
    )
    law123 = (
        "/MAT/LAW123/1",
        "probe_laminated_fracture_daimler_pinho",
        fmt_f(c["rho"]),
        fmt_f(c["e1"], c["e2"], c["e3"]),
        fmt_f(c["g12"], c["g31"], c["g23"]),
        fmt_f(c["nu12"], c["nu31"], c["nu23"]),
        fmt_f(0.1, 0.1, 0.1, 0.1, 0.1),
        fmt_f(c["xc"], c["xt"], c["yc"], c["yt"], c["s12"]),
        fmt_f(0.0, 0.0) + fmt_i(0) + f"{0.0:20.12g}",
        fmt_f(0.0, 1.0, 0.0),
    )
    law125 = (
        "/MAT/LAW125/1",
        "probe_laminated_composite",
        fmt_f(c["rho"]),
        fmt_f(c["e1"], c["e2"], c["e3"]) + (" " * 20) + fmt_i(1),
        fmt_f(c["g12"], c["g31"], c["g23"]),
        fmt_f(c["nu12"], c["nu31"], c["nu23"]),
        fmt_i(0) + f"{0.0:20.12g}" + fmt_i(0) + fmt_f(c["xt"], 0.0),
        fmt_i(0) + f"{0.0:20.12g}" + fmt_i(0) + fmt_f(c["xc"], 0.0),
        fmt_i(0) + f"{0.0:20.12g}" + fmt_i(0) + fmt_f(c["yt"], 0.0),
        fmt_i(0) + f"{0.0:20.12g}" + fmt_i(0) + fmt_f(c["yc"], 0.0),
        fmt_i(0) + f"{0.0:20.12g}" + fmt_i(0) + fmt_f(c["zt"], 0.0),
        fmt_i(0) + f"{0.0:20.12g}" + fmt_i(0) + fmt_f(c["zc"], 0.0),
        fmt_f(0.0, 0.0, 0.0, c["s12"], 0.0),
        fmt_i(0, 0, 0, 0),
        fmt_f(0.0, 0.0, 0.0, c["s13"], 0.0),
        fmt_i(0, 0, 0, 0),
        fmt_f(0.0, 0.0, 0.0, c["s23"], 0.0),
        fmt_i(0, 0, 0, 0),
        fmt_f(0.0, 0.0, 1.0),
        fmt_i(0) + f"{0.0:20.12g}",
        fmt_f(0.0),
    )
    law127 = (
        "/MAT/LAW127/1",
        "probe_enhanced_composite",
        fmt_f(c["rho"]),
        fmt_f(c["e1"], c["e2"], c["e3"]),
        fmt_f(c["g12"], c["g31"], c["g23"]),
        fmt_f(c["nu12"], c["nu31"], c["nu23"]),
        fmt_f(c["xt"], 0.0) + (" " * 10) + fmt_i(0) + f"{0.0:20.12g}",
        fmt_f(c["yt"], 0.0) + (" " * 10) + fmt_i(0) + f"{0.0:20.12g}",
        fmt_f(c["s12"], 0.0) + (" " * 10) + fmt_i(0) + f"{0.0:20.12g}",
        fmt_f(c["xc"], 0.0) + (" " * 10) + fmt_i(0) + f"{0.0:20.12g}",
        fmt_f(c["yc"], 0.0) + (" " * 10) + fmt_i(0) + f"{0.0:20.12g}",
        fmt_f(0.0),
        fmt_f(0.0, 0.0) + fmt_i(0, 0),
        fmt_f(0.0, 0.0, 0.0, 0.0, 1.0),
        (" " * 10) + fmt_i(0) + fmt_f(0.0, 0.0, 0.0),
        fmt_f(0.0, 0.0, 0.0, 0.0),
    )
    law128 = (
        "/MAT/LAW128/1",
        "probe_hill_visc_plast",
        fmt_f(7.8e-9),
        fmt_f(210000.0, 0.30, 200.0, 0.0),
        fmt_i(0) + (" " * 10) + fmt_f(0.0, 0.0),
        fmt_f(0.0, 0.0, 0.0, 0.0),
        fmt_f(0.0, 0.0, 0.0, 0.0),
        fmt_f(0.0, 0.0),
        fmt_f(1.0, 1.0, 1.0),
        fmt_f(0.5, 0.5, 0.5),
        fmt_f(1.5, 1.5, 1.5),
    )
    law132 = (
        "/MAT/LAW132/1",
        "probe_laminated_fracture_daimler_camanho",
        fmt_f(c["rho"]),
        fmt_f(c["e1"], c["e2"], c["e3"]),
        fmt_f(c["g12"], c["g31"], c["g23"]),
        fmt_f(c["nu12"], c["nu31"], c["nu23"]),
        fmt_f(0.1, 0.1, 0.1, 0.1, 0.1),
        fmt_f(c["xc"], c["xt"], c["yc"], c["yt"], c["s12"]),
        fmt_f(0.1, 0.1, c["xc"], c["xt"]),
        fmt_f(0.0, 0.0, 0.0, 0.0) + fmt_i(0),
        fmt_f(0.0, 0.0, 0.0),
        fmt_f(0.0, 0.0, 0.0),
        fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
        fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
        fmt_f(1.0, 0.0),
    )
    law158 = (
        "/MAT/LAW158/1",
        "probe_fabric_nl",
        fmt_f(c["rho"]),
        fmt_f(1.0, 1.0, 0.0, 0.0, 0.0),
        fmt_f(0.0) + (" " * 10) + fmt_i(0),
        fmt_i(1) + (" " * 10) + f"{1.0:20.12g}",
        fmt_i(1) + (" " * 10) + f"{1.0:20.12g}",
        fmt_i(1) + (" " * 10) + f"{1.0:20.12g}",
        fmt_i(0, 0),
    )
    laws = [
        Law("LAW1", "ELAST", "isotropic linear elastic baseline", ("/MAT/LAW1/1", "probe_elast", fmt_f(7.8e-9), fmt_f(210000.0, 0.30))),
        Law("LAW2", "PLAS_JOHNS", "Johnson-Cook isotropic plastic baseline", ("/MAT/LAW2/1", "probe_plas_johns", fmt_f(7.8e-9), fmt_f(210000.0, 0.30) + fmt_i(0), fmt_f(200.0, 210.0, 1.0, 0.0, 0.0), fmt_f(0.0, 1.0, 0.0, 0.0, 0.0, 0.0), fmt_f(0.0, 0.0, 0.0, 0.0))),
        Law("LAW12", "3D_COMP", "3D solid composite orthotropic", law12),
        Law("LAW14", "COMPSO", "legacy composite solid orthotropic", law14),
        Law("LAW15", "CHANG", "Chang composite shell material", law15),
        Law("LAW25", "COMPSH", "CRASURV/COMPSH composite shell card", law25),
        Law("LAW28", "HONEYCOMB", "orthotropic honeycomb, not a fiber composite ply", law28),
        Law("LAW42", "OGDEN", "hyperelastic rubber control", law42),
        Law("LAW53", "TSAI_TAB", "Tsai tabulated anisotropic law", law53),
        Law("LAW58", "FABR_A", "fabric anisotropic material", law58),
        Law("LAW62", "VISC_HYP", "visco-hyperelastic rubber control", law62),
        Law("LAW123", "LAMINATED_FRACTURE_DAIMLER_PINHO", "laminated composite fracture, Pinho", law123),
        Law("LAW125", "LAMINATED_COMPOSITE", "laminated composite with built-in failure", law125),
        Law("LAW127", "ENHANCED_COMPOSITE", "enhanced composite", law127),
        Law("LAW128", "HILL_VISC_PLAST", "Hill anisotropic viscoplastic metal", law128),
        Law("LAW132", "LAMINATED_FRACTURE_DAIMLER_CAMANHO", "laminated composite fracture, Camanho", law132),
        Law("LAW158", "FABR_NL", "nonlinear fabric material", law158),
    ]
    return {law.law: law for law in laws}


FAILS = {
    "base": (),
    "HASHIN": (
        "/FAIL/HASHIN/1",
        fmt_i(1, 0, 1) + f"{1.0:20.12g}",
        fmt_f(1500.0, 50.0, 50.0, 1200.0, 200.0),
        fmt_f(200.0, 90.0, 90.0, 70.0, 90.0),
        fmt_f(0.0, 1.0, 0.0) + (" " * 20) + f"{0.0:20.12g}",
    ),
    "PUCK": (
        "/FAIL/PUCK/1",
        fmt_f(1500.0, 50.0, 90.0, 1200.0, 200.0),
        fmt_f(0.30, 0.25, 0.25, 0.0) + fmt_i(0, 1),
        fmt_f(0.0),
    ),
    "TSAIWU": (
        "/FAIL/TSAIWU/1",
        fmt_f(1500.0, 50.0, 1200.0, 200.0, 90.0),
        fmt_f(-0.5, 0.0, 0.0) + (" " * 20) + fmt_i(0, 1),
    ),
    "CHANG": (
        "/FAIL/CHANG/1",
        fmt_f(1500.0, 50.0, 90.0, 1200.0, 200.0),
        fmt_f(0.0, 0.0) + fmt_i(0, 1) + f"{0.0:20.12g}",
    ),
    "USER1": (
        "/FAIL/USER1/1",
        "0",
    ),
    "JOHNSON": (
        "/FAIL/JOHNSON/1",
        fmt_f(0.1, 0.0, 0.0, 0.0, 0.0),
        fmt_f(1.0) + fmt_i(0, 1) + fmt_f(0.0, 0.0) + (" " * 10) + fmt_i(0),
        fmt_i(0),
    ),
}


def type14_prop() -> tuple[str, ...]:
    return (
        "/PROP/TYPE14/1",
        "probe_type14_solid",
        f"{0:10d}{0:10d}{0:20d}{0:20d}{0:20d}{0.0:20.12g}",
        fmt_f(0.0, 0.0, 0.0),
        fmt_f(0.0),
    )


def type6_prop() -> tuple[str, ...]:
    return (
        "/PROP/TYPE6/1",
        "probe_type6_sol_orth",
        f"{0:10d}{0:10d}{0:20d}{0:20d}{0:20d}{0.0:20.12g}",
        fmt_f(0.0, 0.0, 0.0),
        fmt_f(1.0, 0.0, 0.0) + fmt_i(0, 1, 0),
        fmt_f(0.0, 0.0, 0.0, 0.0),
        fmt_f(0.0) + fmt_i(0, 0),
    )


def type1_prop() -> tuple[str, ...]:
    return (
        "/PROP/TYPE1/1",
        "probe_type1_shell",
        fmt_i(0, 0, 0, 0) + (" " * 30) + f"{0.0:20.12g}",
        fmt_f(0.0, 0.0, 0.0, 0.0, 0.0),
        fmt_i(1, 0) + f"{1.0:20.12g}{0.0:20.12g}{0:20d}{0:10d}",
    )


PROPS = {
    "TYPE14": ("BRICK", type14_prop),
    "TYPE6": ("BRICK", type6_prop),
    "TYPE1": ("SHELL", type1_prop),
}


def function_block(law: str) -> tuple[str, ...]:
    # Several tabulated/fabric laws reference function 1 in their minimal cards.
    return (
        "/FUNCT/1",
        f"{law}_constant_function",
        fmt_f(-1.0e9, 1.0),
        fmt_f(1.0e9, 1.0),
    )


def geometry_block(element_kind: str) -> tuple[str, ...]:
    if element_kind == "BRICK":
        return (
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
            "probe_part",
            fmt_i(1, 1, 0),
            "/BRICK/1",
            fmt_i(1, 1, 2, 3, 4, 5, 6, 7, 8),
        )
    return (
        "/NODE",
        fmt_i(1) + fmt_f(0.0, 0.0, 0.0),
        fmt_i(2) + fmt_f(1.0, 0.0, 0.0),
        fmt_i(3) + fmt_f(1.0, 1.0, 0.0),
        fmt_i(4) + fmt_f(0.0, 1.0, 0.0),
        "/PART/1",
        "probe_part",
        fmt_i(1, 1, 0),
        "/SHELL/1",
        fmt_i(1, 1, 2, 3, 4),
    )


def load_block(element_kind: str) -> tuple[str, ...]:
    if element_kind == "BRICK":
        fixed_nodes = fmt_i(1, 4, 5, 8)
        loaded_nodes = fmt_i(2, 3, 6, 7)
    else:
        fixed_nodes = fmt_i(1, 4)
        loaded_nodes = fmt_i(2, 3)
    return (
        "/BCS/1",
        "probe_clamp",
        f"   111 111{0:10d}{1:10d}",
        "/GRNOD/NODE/1",
        "probe_fixed_nodes",
        fixed_nodes,
        "/CLOAD/1",
        "probe_unit_load",
        f"{99:10d}         X{0:10d}{0:10d}{2:10d}{1.0:30.12g}{1.0:20.12g}",
        "/GRNOD/NODE/2",
        "probe_loaded_nodes",
        loaded_nodes,
        "/FUNCT/99",
        "probe_constant_load",
        fmt_f(0.0, 1.0),
        fmt_f(1.0, 1.0),
    )


def starter_deck(law: Law, prop: str, fail: str) -> str:
    element_kind, prop_func = PROPS[prop]
    lines: list[str] = [
        "#RADIOSS STARTER",
        "/BEGIN",
        f"deck_{law.law}_{prop}_{fail}",
        "      2026         0",
        f"{'kg':>20}{'mm':>20}{'ms':>20}",
        f"{'kg':>20}{'mm':>20}{'ms':>20}",
        "/TITLE",
        f"{law.law} {law.name} {prop} {fail} compatibility probe",
        "/DEF_SOLID",
        f"{0:10d}{0:10d}{0:20d}{0:40d}",
        "/DEF_SHELL",
        fmt_i(0, 0, 0, 0, 0) + f"{1:30d}{0:10d}",
    ]
    lines.extend(law.cards)
    if fail != "base":
        lines.extend(FAILS[fail])
    lines.extend(geometry_block(element_kind))
    lines.extend(load_block(element_kind))
    lines.extend(prop_func())
    if law.law in {"LAW28", "LAW53", "LAW158"}:
        lines.extend(function_block(law.law))
    lines.extend(("/END", ""))
    return "\n".join(lines)


def classify(text: str, rc: int) -> tuple[str, str, str]:
    error_ids = re.findall(r"ERROR ID\s*:\s*(\d+)", text)
    warning_ids = re.findall(r"WARNING ID\s*:\s*(\d+)", text)
    if rc == 0 and "ERROR TERMINATION" not in text and not error_ids:
        return "OK", "", ",".join(warning_ids)
    if error_ids:
        return f"BLOCKED:{error_ids[0]}", error_ids[0], ",".join(warning_ids)
    if "ERROR TERMINATION" in text:
        return "BLOCKED:ERROR_TERMINATION", "", ",".join(warning_ids)
    return f"BLOCKED:RC{rc}", "", ",".join(warning_ids)


def run_one(deck_path: Path) -> tuple[int, str]:
    proc = subprocess.run(
        [str(STARTER), "-i", str(deck_path.name), "-check"],
        cwd=str(THIS),
        env=env(),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    listing = deck_path.with_suffix(".out")
    starter_listing = deck_path.with_name(deck_path.stem.removesuffix("_0000") + "_0000.out")
    text = proc.stdout
    if starter_listing.exists():
        text = starter_listing.read_text(encoding="utf-8", errors="replace")
        if listing != starter_listing:
            listing.write_text(text, encoding="utf-8")
    else:
        listing.write_text(text, encoding="utf-8")
    return proc.returncode, text


def clean_generated() -> None:
    for pattern in (
        "deck_LAW*.rad",
        "deck_LAW*.out",
        "deck_LAW*_0000.out",
        "deck_LAW*_0000.sta",
        "deck_LAW*_0000.rad",
        "deck_LAW*.ctl",
        "deck_LAW*.log",
        "results.csv",
        "starter_version.txt",
        "cfg_candidates.txt",
    ):
        for path in THIS.glob(pattern):
            if path.name == Path(__file__).name:
                continue
            path.unlink(missing_ok=True)


def starter_version() -> str:
    proc = subprocess.run(
        [str(STARTER), "-version"],
        cwd=str(THIS),
        env=env(),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return proc.stdout


def main() -> int:
    clean_generated()
    laws = law_cards()
    (THIS / "starter_version.txt").write_text(starter_version(), encoding="utf-8")
    rows: list[dict[str, str | int]] = []
    for law in laws.values():
        for prop in PROPS:
            for fail in FAILS:
                name = f"deck_{law.law}_{prop}_{fail}"
                deck_path = THIS / f"{name}_0000.rad"
                deck_path.write_text(starter_deck(law, prop, fail), encoding="utf-8")
                rc, text = run_one(deck_path)
                verdict, error_id, warnings = classify(text, rc)
                rows.append(
                    {
                        "law": law.law,
                        "name": law.name,
                        "prop": prop,
                        "fail": fail,
                        "return_code": rc,
                        "verdict": verdict,
                        "error_id": error_id,
                        "warnings": warnings,
                        "deck": deck_path.name,
                        "out": deck_path.with_suffix(".out").name,
                        "note": law.note,
                    }
                )
                print(f"{law.law:6s} {prop:6s} {fail:8s} {verdict}")
    with (THIS / "results.csv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
