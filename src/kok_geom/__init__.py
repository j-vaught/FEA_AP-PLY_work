"""Clean-room AP-PLY geometry generation for GMSH-OCC."""

from kok_geom.config import KokConfig
from kok_geom.geometry.laminate import Laminate, LaminateSolid
from kok_geom.geometry.ply import Ply, PlySolid
from kok_geom.geometry.tow import Tow, TowSolid, UndulationProfile
from kok_geom.io import generate_mesh

__all__ = [
    "KokConfig",
    "Laminate",
    "LaminateSolid",
    "Ply",
    "PlySolid",
    "Tow",
    "TowSolid",
    "UndulationProfile",
    "generate_mesh",
]

__version__ = "0.1.0"
