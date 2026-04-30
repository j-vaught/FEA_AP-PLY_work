"""Mesh generation helpers."""

from __future__ import annotations

import gmsh


def configure_msh_format(msh_format: str) -> None:
    if msh_format == "msh4_ascii":
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Binary", 0)
    elif msh_format == "msh4_binary":
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Binary", 1)
    elif msh_format == "msh2_ascii":
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
    else:
        raise ValueError(f"Unsupported msh_format: {msh_format}")


def generate_volume_mesh(*, target_size_m: float, element_order: int = 2) -> None:
    """Generate a 3-D mesh using TETRA10 by default."""

    if target_size_m <= 0.0:
        raise ValueError("target_size_m must be positive")
    gmsh.model.occ.synchronize()
    entities = gmsh.model.getEntities(0)
    gmsh.model.mesh.setSize(entities, target_size_m)
    gmsh.option.setNumber("Mesh.ElementOrder", element_order)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    if element_order == 2:
        # Netgen optimization rewrites the linear mesh, so promote after it.
        gmsh.model.mesh.setOrder(2)
