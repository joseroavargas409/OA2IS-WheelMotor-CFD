"""Extruded (2.5D) prism/hex meshes generated with gmsh for overset components.

A component is a planar region (outer loop minus solid holes) extruded along Z
between z0 and z1 with a fixed number of layers. Boundary patches:
  overset  : outer loop side walls (+ optionally the two end caps)
  <name>   : side walls of every hole (moving gear / carrier surfaces)
  caps     : end caps when they are walls (slab tests)
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import gmsh
import numpy as np

from .foam import env as foam_env


def _loop(points, lc, tag_list):
    pts = [gmsh.model.geo.addPoint(float(x), float(y), 0.0, lc) for x, y in points]
    lines = [gmsh.model.geo.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
    tag_list.extend(lines)
    return gmsh.model.geo.addCurveLoop(lines)


def build(path: Path, outer: np.ndarray, holes: dict[str, list[np.ndarray]], z0: float, z1: float,
          layers: int, lc_fine: float, lc_coarse: float, fine_dist: float,
          cap_patch: str = 'overset', recombine: bool = False, outer_patch: str = 'overset',
          zone: str = 'fluid'):
    """Write a gmsh v2 ASCII mesh with physical groups ready for gmshToFoam."""
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.model.add(path.stem)
    outer_lines, hole_lines = [], {k: [] for k in holes}
    loops = [_loop(outer, lc_coarse, outer_lines)]
    for name, polys in holes.items():
        for poly in polys:
            loops.append(_loop(poly, lc_fine, hole_lines[name]))
    surf = gmsh.model.geo.addPlaneSurface(loops)
    gmsh.model.geo.synchronize()
    # size field: fine near holes, coarse far away
    fd = gmsh.model.mesh.field.add('Distance')
    gmsh.model.mesh.field.setNumbers(fd, 'CurvesList', [l for v in hole_lines.values() for l in v])
    gmsh.model.mesh.field.setNumber(fd, 'Sampling', 4)
    th = gmsh.model.mesh.field.add('Threshold')
    gmsh.model.mesh.field.setNumber(th, 'InField', fd)
    gmsh.model.mesh.field.setNumber(th, 'SizeMin', lc_fine)
    gmsh.model.mesh.field.setNumber(th, 'SizeMax', lc_coarse)
    gmsh.model.mesh.field.setNumber(th, 'DistMin', 0.0)
    gmsh.model.mesh.field.setNumber(th, 'DistMax', fine_dist)
    gmsh.model.mesh.field.setAsBackgroundMesh(th)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.Algorithm', 8 if recombine else 6)
    gmsh.option.setNumber('Mesh.RecombinationAlgorithm', 1)
    gmsh.model.geo.translate([(2, surf)], 0, 0, z0)
    if recombine:
        gmsh.model.geo.mesh.setRecombine(2, surf)
    ext = gmsh.model.geo.extrude([(2, surf)], 0, 0, z1 - z0, numElements=[layers], recombine=True)
    gmsh.model.geo.synchronize()
    # ext: [top surface, volume, side surfaces in loop-line order...]
    top, vol, sides = ext[0][1], ext[1][1], [e[1] for e in ext[2:]]
    n_outer = len(outer_lines)
    gmsh.model.addPhysicalGroup(3, [vol], name=zone)
    side_outer = sides[:n_outer]
    k = n_outer
    groups = {outer_patch: list(side_outer)}
    for name, lines in hole_lines.items():
        groups.setdefault(name, []).extend(sides[k:k + len(lines)])
        k += len(lines)
    groups.setdefault(cap_patch, []).extend([surf, top])
    for name, tags in groups.items():
        gmsh.model.addPhysicalGroup(2, tags, name=name)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.write(str(path))
    n = len(gmsh.model.mesh.getElements(3)[1][0]) if gmsh.model.mesh.getElements(3)[1] else 0
    gmsh.finalize()
    return n


def to_foam(msh: Path, case: Path, log: Path):
    """Convert with gmshToFoam into case/constant/polyMesh."""
    with open(log, 'w') as f:
        subprocess.run(['gmshToFoam', str(msh.resolve()), '-case', str(case)], check=True, stdout=f,
                       stderr=subprocess.STDOUT, env=foam_env())
