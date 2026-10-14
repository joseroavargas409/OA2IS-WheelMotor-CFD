#!/usr/bin/env python3
"""Create a runnable case from a mesh and a scenario.

usage: python3 tools/setup_case.py <mesh_dir> <scenario.yaml> <case_dir> [--np 2]
Steps: write dictionaries, zoneID (setFields), checkMesh, initial oil (flat
surface normal to g_eff(0) holding the scenario litres), decomposePar.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oa2is_cfd import foam                                          # noqa: E402
from oa2is_cfd.case import set_initial_oil, write_case              # noqa: E402
from oa2is_cfd.model import default_case                            # noqa: E402
from oa2is_cfd.scenario import Scenario                             # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('mesh')
ap.add_argument('scenario')
ap.add_argument('case')
ap.add_argument('--np', type=int, default=1)
a = ap.parse_args()
case = Path(a.case)
sc = Scenario.load(a.scenario)
stats = json.loads((Path(a.mesh) / 'mesh_stats.json').read_text())
backlash = tuple(stats.get('backlash', (0.004, 0.003)))
zones = write_case(case, Path(a.mesh), sc, backlash=backlash, n_procs=a.np)
(case / 'mesh_stats.json').write_text((Path(a.mesh) / 'mesh_stats.json').read_text())
cd = case / 'system/controlDict'
cd.write_text(cd.read_text().replace('writeFormat binary;', 'writeFormat ascii;'))
foam.run(case, 'setFields')
foam.run(case, 'checkMesh', '-constant', check=False)
cm = (case / 'log.checkMesh').read_text()
print('checkMesh:', 'Mesh OK' if 'Mesh OK' in cm else 'FAILED checks (see log.checkMesh)')
for line in cm.splitlines():
    if any(k in line for k in ('cells:', 'non-orthogonality Max', 'Max skewness', 'Max aspect', '***')):
        print('  ', line.strip())
g = sc.g_eff(0.0)
n = -g / (g @ g) ** 0.5
gc = default_case(backlash)
level = gc.level_for_volume(sc.oil_litres, normal=n)
alpha, C, V = set_initial_oil(case, sc, level, n)
cd.write_text(cd.read_text().replace('writeFormat ascii;', 'writeFormat binary;'))
meta = dict(level=level, up=list(map(float, n)), litres=sc.oil_litres, zones=zones)
(case / 'initial_oil.json').write_text(json.dumps(meta, indent=1))
print(f'initial free surface: n.x = {level:.4f} m, n = {n.round(4).tolist()}, {sc.oil_litres} L')
if a.np > 1:
    foam.run(case, 'decomposePar', '-force')
    print('decomposed for', a.np, 'processes')
