#!/usr/bin/env python3
"""Immersed-boundary pipeline: mesh (optional) + case for a scenario.

usage: python3 tools/build_ib_case.py <scenario.yaml> <case_dir> [--mesh DIR] [--res ib_coarse]
                                     [--maxco 1.0] [--duration T]
If --mesh DIR has no mesh yet it is built there.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oa2is_cfd import foam                                                 # noqa: E402
from oa2is_cfd.ib import IB_PRESETS, build_ib_mesh, initial_oil_ib, write_ib_case   # noqa: E402
from oa2is_cfd.model import default_case                                   # noqa: E402
from oa2is_cfd.scenario import Scenario                                    # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('scenario')
ap.add_argument('case')
ap.add_argument('--mesh', default=None)
ap.add_argument('--res', default='ib_coarse', choices=sorted(IB_PRESETS))
ap.add_argument('--maxco', type=float, default=None, help='override max Courant (default 1.0)')
ap.add_argument('--duration', type=float, default=None, help='override scenario duration [s]')
a = ap.parse_args()
sc = Scenario.load(a.scenario)
if a.maxco:
    sc.numerics['maxCo'] = a.maxco
if a.duration:
    sc.duration = a.duration
gc = default_case()
mesh = Path(a.mesh or Path(a.case).parent / f'mesh_{a.res}')
if not (mesh / 'constant/polyMesh/faces').exists():
    print('building mesh', mesh, flush=True)
    print(json.dumps(build_ib_mesh(mesh, gc, IB_PRESETS[a.res], log=lambda s: print(s, flush=True))['total_cells']))
write_ib_case(Path(a.case), mesh, sc, gc)
foam.run(a.case, 'checkMesh', '-constant', check=False)
cm = (Path(a.case) / 'log.checkMesh').read_text()
print('checkMesh:', 'Mesh OK' if 'Mesh OK' in cm else 'FAILED (see log.checkMesh)')
print(json.dumps(initial_oil_ib(Path(a.case), sc, gc), indent=1))
