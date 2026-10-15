#!/usr/bin/env python3
"""Checks of the rotating hub-cover plug sensors (no CFD run needed):
  * electrode stencil: 9 points, long axis along the requested axis
  * plug ring: all ring points in fluid at rest (outside every solid of the model)
  * plug angle kinematics: theta_hub(t) = ratio_out * integral(omega_rotor)
usage: python3 tests/test_plug_sensors.py
"""
import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'tools'))
from oa2is_cfd.ib import chi_at_rest, electrode_stencil, plug_ring   # noqa: E402
from oa2is_cfd.model import default_case                              # noqa: E402
from oa2is_cfd.scenario import Scenario                               # noqa: E402

off = np.array(electrode_stencil(2))
assert off.shape == (9, 3) and np.allclose(off.sum(0), 0)
assert np.isclose(np.ptp(off[:, 2]), 0.025) and np.isclose(np.ptp(off[:, 0]), 0.012)

sc = Scenario.load(ROOT / 'scenarios/07_real_acarreo.yaml')
g = sc.rotating_sensors[0]
ang, roff, pts = plug_ring(g)
assert len(pts) == g['n_ring'] * 9
solid = chi_at_rest(default_case(), np.array(pts))
assert not solid.any(), 'plug ring points inside a solid'

# kinematics against the closed form for the 07 table (ramp 0-0.05 s, then constant)
from plug_sensors import part_angle                                  # noqa: E402
with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    meta = dict(sc.__dict__, ratios=dict(inp=1.0, mid=0.2, out=0.03066666666666667))
    (d / 'scenario.json').write_text(json.dumps(meta, default=float))
    t = np.array([0.05, 0.35])
    th = part_angle(d, 'out', t)
    w = sc.rotor_sign() * sc.motor_rpm[1][1] * 2 * math.pi / 60 * 0.03066666666666667
    exact = np.array([0.5 * w * 0.05, 0.5 * w * 0.05 + w * 0.30])
    assert np.allclose(th, exact, rtol=1e-6), (th, exact)
print('plug sensors OK: stencil, ring in fluid, hub angle at 0.35 s = %.2f deg' % math.degrees(th[1]))
