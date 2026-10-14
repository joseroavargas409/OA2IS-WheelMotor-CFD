#!/usr/bin/env python3
"""Check, inside the CFD results, that the immersed bodies are where the
analytic kinematics says they should be.

For every written time, chi (solid fraction written by interIBFoamOA2IS) is
sampled on the planet orbit circle at the mid-plane of each stage. The three
planets appear as three arcs with chi > 0.5; the arc mid-points are the
measured orbit angles and are compared with the analytic carrier angle
(integral of the scenario's omega(t) table x stage ratio).

usage: python3 tests/test_ib_motion.py <case_dir> [--csv out.csv]
Pass criterion: |error| < 1.5 cell sizes along the orbit.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oa2is_cfd.ibpost import Fields, load_scenario      # noqa: E402
from oa2is_cfd.kinematics import stages                 # noqa: E402


def rotor_angle(table, t):
    ts, ws = np.array([r[0] for r in table]), np.array([r[1] for r in table])
    grid = np.unique(np.r_[ts[ts < t], t, 0.0])
    return float(np.trapezoid(np.interp(grid, ts, ws), grid))


def arcs(mask):
    """Mid-angles (index units) of circular runs of True."""
    n = len(mask)
    if mask.all() or not mask.any():
        return []
    start = np.argmin(mask)          # a False index: rotate so runs do not wrap
    m = np.roll(mask, -start)
    d = np.diff(np.r_[0, m.astype(int), 0])
    s, e = np.where(d == 1)[0], np.where(d == -1)[0]
    return [((a + b - 1) / 2 + start) % n for a, b in zip(s, e)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('--csv')
    a = ap.parse_args()
    case = Path(a.case)
    sc = load_scenario(case)
    meta = json.loads((case / 'scenario.json').read_text())
    lc = json.loads((case / 'mesh_stats.json').read_text())['resolution']['lc_moving']
    sts = stages(tuple(meta['backlash']))
    F = Fields(case)
    table = sc.rotor_omega_table()
    rows = []
    N = 1440
    for t in F.times:
        m = F.mesh(t)
        if 'chi' not in m.cell_data:
            continue
        theta = rotor_angle(table, t)
        for st in sts:
            ang = np.linspace(0, 2 * math.pi, N, endpoint=False)
            pts = np.c_[st.orbit_r * np.cos(ang), st.orbit_r * np.sin(ang), np.full(N, st.z + 0.002)]
            chi = pv.PolyData(pts).sample(m).point_data['chi']
            mids = arcs(chi > 0.5)
            meas = sorted((np.array(mids) * 2 * math.pi / N) % (2 * math.pi))
            exp = sorted([(j * 2 * math.pi / 3 + st.carrier_angle(theta)) % (2 * math.pi) for j in range(3)])
            if len(meas) != 3:
                rows.append(dict(t=t, stage=st.id, n_arcs=len(meas), err_deg=np.nan, err_m=np.nan))
                continue
            err = []
            for e in exp:
                d = [(mm - e + math.pi) % (2 * math.pi) - math.pi for mm in meas]
                err.append(min(d, key=abs))
            emax = max(abs(x) for x in err)
            rows.append(dict(t=t, stage=st.id, n_arcs=3, carrier_deg=math.degrees(st.carrier_angle(theta)),
                             err_deg=math.degrees(emax), err_m=emax * st.orbit_r))
    import pandas as pd
    df = pd.DataFrame(rows)
    if a.csv:
        df.to_csv(a.csv, index=False)
    print(df.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    ok = bool((df.n_arcs == 3).all() and (df.err_m < 1.5 * lc).all())
    print(f'max error {df.err_m.max() * 1e3:.2f} mm (cell size {lc * 1e3:.0f} mm) ->', 'PASSED' if ok else 'FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
