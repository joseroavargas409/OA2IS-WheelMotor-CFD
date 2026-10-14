#!/usr/bin/env python3
"""Export the PROCESSED CFD geometry (after the repairs of docs/GEOMETRIA.md)
as STL surfaces, rest pose t = 0:

  geometry/cfd/wall_<patch>.stl         body-fitted walls (fixed spindle/rings; axisymmetric rotating hub, shaft, sleeve)
  geometry/cfd/moving_<body>.stl        immersed moving bodies (suns, planets, carriers)
  geometry/cfd/static_ring<Stage>.stl   static immersed ring-gear tooth bands
  geometry/cfd/cavity_envelope.stl      outer envelope of the oil cavity
  geometry/cfd/geometry_summary.json    volumes, z ranges, provenance notes

usage: python3 tools/export_cfd_geometry.py [out_dir]
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import trimesh

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oa2is_cfd.model import REPAIRS, Z_CLOSURE, Z_COVER, R_HUB, default_case   # noqa: E402

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / 'geometry' / 'cfd')
OUT.mkdir(parents=True, exist_ok=True)


def prism_mesh(p):
    poly = p.polygon()
    polys = list(poly.geoms) if poly.geom_type == 'MultiPolygon' else [poly]
    meshes = []
    for g in polys:
        if g.area < 1e-10:
            continue
        m = trimesh.creation.extrude_polygon(g, p.z1 - p.z0)
        m.apply_translation([0, 0, p.z0])
        meshes.append(m)
    return meshes


def main():
    gc = default_case()
    groups = defaultdict(list)
    notes = defaultdict(set)
    for p in gc.background_solids(for_volume=True):
        groups['wall_' + p.patch] += prism_mesh(p)
        notes['wall_' + p.patch].add(p.note)
    for c in gc.components():
        for p in c['solids']:
            groups['moving_' + c['body']] += prism_mesh(p)
            notes['moving_' + c['body']].add(p.note)
    for z in gc.ring_zones():
        for p in z['solids']:
            groups['static_' + z['body']] += prism_mesh(p)
            notes['static_' + z['body']].add(p.note)
    env = trimesh.creation.cylinder(radius=R_HUB, height=Z_CLOSURE - Z_COVER, sections=180)
    env.apply_translation([0, 0, 0.5 * (Z_COVER + Z_CLOSURE)])
    groups['cavity_envelope'] = [env]
    summary = dict(units='m', frame='model frame: z = wheel axis (inboard +), y up; rest pose t = 0',
                   cavity=dict(z_cover=Z_COVER, z_closure=Z_CLOSURE, r_hub=R_HUB),
                   fluid_volume_geometric_L=gc.fluid_volume() * 1e3, repairs=REPAIRS, parts={})
    for name, ms in sorted(groups.items()):
        m = trimesh.util.concatenate(ms)
        m.export(OUT / f'{name}.stl')
        summary['parts'][name] = dict(triangles=int(len(m.faces)), bbox=np.round(m.bounds, 4).tolist(),
                                      notes=sorted(n for n in notes[name] if n))
        print(f'{name:28s} {len(m.faces):7d} tris')
    (OUT / 'geometry_summary.json').write_text(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
