#!/usr/bin/env python3
"""Regenerated gears: no interpenetration and positive clearance over a full
mesh cycle of both planetary stages (sun-planet and planet-ring pairs), using
the analytic kinematics (same poses the immersed-boundary solver uses).

usage: python3 tests/test_gears.py        (exit code 0 = pass)
"""
import math
import sys
from pathlib import Path

import numpy as np
from shapely.geometry import Polygon

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oa2is_cfd.gears import external_profile, ring_space_profile, rotate   # noqa: E402
from oa2is_cfd.model import circle, default_case                         # noqa: E402


def main(n=60):
    gc = default_case()
    ok = True
    for st in gc.stages:
        sun0 = external_profile(st.sun_spec())
        pl0 = external_profile(st.planet_spec())
        ring = Polygon(circle(st.pitch + 3 * st.module, h=0.002), [ring_space_profile(st.ring_spec())]).buffer(0)
        # one tooth-mesh period of the sun relative to the carrier
        period = 2 * math.pi / st.ns / abs(st.input - st.output)
        worst = dict(ov_sp=0.0, ov_pr=0.0, gap_sp=1.0, gap_pr=1.0)
        for theta in np.linspace(0, period, n):
            sun = Polygon(rotate(sun0, st.sun_angle(theta)))
            for j in range(3):
                c = np.array(st.planet_center(j, theta))
                pl = Polygon(rotate(pl0, st.planet_angle(j, theta)) + c)
                worst['ov_sp'] = max(worst['ov_sp'], sun.intersection(pl).area)
                worst['ov_pr'] = max(worst['ov_pr'], ring.intersection(pl).area)
                worst['gap_sp'] = min(worst['gap_sp'], sun.distance(pl))
                worst['gap_pr'] = min(worst['gap_pr'], ring.distance(pl))
        print(f"stage {st.id}: backlash {st.backlash * 1e3:.1f} mm | overlap sun-planet {worst['ov_sp'] * 1e6:.3f} mm2, "
              f"planet-ring {worst['ov_pr'] * 1e6:.3f} mm2 | min gap sun-planet {worst['gap_sp'] * 1e3:.2f} mm, "
              f"planet-ring {worst['gap_pr'] * 1e3:.2f} mm")
        ok &= worst['ov_sp'] < 1e-9 and worst['ov_pr'] < 1e-9 and worst['gap_sp'] > 0 and worst['gap_pr'] > 0
    print('PASSED' if ok else 'FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
