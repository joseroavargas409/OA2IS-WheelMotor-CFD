#!/usr/bin/env python3
"""Verify that OpenFOAM moves every cellZone exactly as the analytic planetary
kinematics (sun spin, planet spin + orbit, carrier rotation, static rings).

usage: python3 tests/test_kinematics_openfoam.py <case_dir> [t_end]
Runs moveDynamicMesh on a copy of the case up to t_end and compares the cell
centres of every zone with the analytic pose of its rigid body.
"""
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oa2is_cfd import foam                                   # noqa: E402
from oa2is_cfd.case import read_vector_field, zone_names     # noqa: E402
from oa2is_cfd.kinematics import stages                      # noqa: E402


def rotor_angle(table, t):
    """Integral of a piecewise-linear omega(t) table from 0 to t."""
    ts, ws = np.array([r[0] for r in table]), np.array([r[1] for r in table])
    grid = np.unique(np.r_[ts[ts < t], t, 0.0])
    w = np.interp(grid, ts, ws)
    return float(np.trapezoid(w, grid))


def cell_zone_ids(case):
    txt = (Path(case) / 'constant/polyMesh/cellZones').read_text()
    out = {}
    for name in zone_names(case):
        i = txt.index('\n' + name + '\n')
        j = txt.index('(', txt.index('cellLabels', i))
        k = txt.index(')', j)
        out[name] = np.array(txt[j + 1:k].split(), dtype=np.int64)
    return out


def main(case, t_end=0.02):
    case = Path(case)
    work = Path(tempfile.mkdtemp(prefix='kin_'))
    for d in ('constant', 'system'):
        shutil.copytree(case / d, work / d)
    shutil.copytree(case / '0', work / '0')
    sc = json.loads((case / 'scenario.json').read_text())
    cd = work / 'system/controlDict'
    t = cd.read_text()
    for k, v in (('endTime', t_end), ('deltaT', t_end / 2), ('writeInterval', t_end / 2)):
        t = __import__('re').sub(rf'\n{k} [^;]*;', f'\n{k} {v:g};', t)
    t = t.replace('adjustTimeStep yes;', 'adjustTimeStep no;').replace('writeFormat binary;', 'writeFormat ascii;')
    t = t.replace('writeControl adjustableRunTime;', 'writeControl timeStep;').replace(f'writeInterval {t_end/2:g};', 'writeInterval 1;')
    cd.write_text(t)
    foam.run(work, 'moveDynamicMesh', log='move')
    foam.run(work, 'postProcess', '-func', 'writeCellCentres', '-time', f'0,{t_end:g}', log='cc')
    C0 = read_vector_field(work / '0/C')
    C1 = read_vector_field(work / f'{t_end:g}/C')
    from oa2is_cfd.scenario import Scenario
    s = Scenario(**{k: v for k, v in sc.items() if k in Scenario.__dataclass_fields__})
    th = rotor_angle(s.rotor_omega_table(), t_end)
    st = {x.id: x for x in stages(tuple(sc['backlash']))}
    zones = cell_zone_ids(work)
    R = lambda a: np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    report = {}
    for name, ids in zones.items():
        p0 = C0[ids, :2]
        if name.startswith('sun'):
            S = st[name[3:]]
            exp = p0 @ R(S.input * th).T
        elif name.startswith('planet'):
            S = st[name[6:-1]]
            j = int(name[-1])
            c = np.array(S.planet_center(j, 0))
            spun = (p0 - c) @ R(S.planet_rel * th).T + c
            exp = spun @ R(S.output * th).T
        elif name.startswith('carrier'):
            S = st['High' if 'High' in name else 'Low']
            exp = p0 @ R(S.output * th).T
        else:
            exp = p0
        err = np.abs(C1[ids, :2] - exp).max()
        report[name] = float(err)
    ok = all(v < 1e-6 for v in report.values())
    print(json.dumps(dict(t_end=t_end, rotor_angle_rad=th, max_error_m=report, passed=ok), indent=1))
    shutil.rmtree(work)
    return ok


if __name__ == '__main__':
    sys.exit(0 if main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 0.02) else 1)
