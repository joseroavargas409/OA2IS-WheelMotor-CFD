"""Capacitive oil/air sensor model applied to CFD probe data.

Python port of the sensor model of komatsu-oil-sensors-v4-1
(dashboard/sensor-model.js, MIT License, (c) 2026 OA2IS simulation contributors;
copied unchanged in third_party/komatsu-oil-sensors-v4-1/sensor-model.js).
What is kept from the original, for each sensor:

* sensitive volume = 9 points (centre, +-6 mm across, +-6.25 / +-12.5 mm along the
  coaxial electrode, which is vertical); points inside a solid are excluded and
  a sensor whose 9 points are all inside solids is 'unavailable' (never 'air');
* oil fraction = mean of the CFD volume fraction alpha.oil over the valid points;
* coaxial capacitance (a = 4 mm, b = 6 mm, L = 25 mm, eps_oil 2.2, eps_air 1.0,
  8 pF stray) with a residual oil film (target 50 um x oil fraction, filling
  time constant 0.1 s, draining 0.7 s);
* first-order electronics (tau = 0.08 s) and Gaussian noise 0.01 pF from the
  same xorshift32 generator and seeds (930 + (i+1)*2654435761);
* signal = (C - C_dry) / (C_wet - C_dry), state dry < 0.2 <= mixed < 0.8 <= wet,
  wet flag signal > 0.6.

Difference from the original: the oil field comes from the OpenFOAM VOF
solution (probes of alpha.oil and of the immersed-body fraction chi) instead
of SPH particles, so the 'bulk/detached' split (particle connectivity) and the
level column are not computed.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

EPS0 = 8.8541878128e-12
A, B, LENGTH, EPS_OIL, EPS_GAS, STRAY = .004, .006, .025, 2.2, 1.0, 8e-12
TAU_S = 0.08          # FIXED.sensor_tau_s
NOISE_PF = 0.01       # FIXED.sensor_noise_pf
SEED = 930            # FIXED.sensor_seed


def coax_capacitance(oil_fraction, film=0.0):
    """pF; identical to coaxCapacitance() in sensor-model.js."""
    film = max(0.0, min(.00015, film))
    mix = EPS_GAS + (EPS_OIL - EPS_GAS) * min(1.0, max(0.0, oil_fraction))
    inv = (math.log((A + film) / A) / EPS_OIL + math.log((B - film) / (A + film)) / mix
           + math.log(B / (B - film)) / EPS_OIL)
    return (STRAY + 2 * math.pi * EPS0 * LENGTH / inv) * 1e12


C_DRY, C_WET = coax_capacitance(0.0), coax_capacitance(1.0)


class _XorShift:
    def __init__(self, i):
        self.a = (SEED + (i + 1) * 2654435761) & 0xFFFFFFFF

    def __call__(self):
        a = self.a
        a ^= (a << 13) & 0xFFFFFFFF
        a ^= a >> 17
        a ^= (a << 5) & 0xFFFFFFFF
        self.a = a & 0xFFFFFFFF
        return (self.a + .5) / 4294967296


def run(times, alpha, solid, ids, chi_threshold=0.5):
    """times [n]; alpha, solid [n, n_sensors, 9] (solid: bool, True = point in a
    solid or outside the fluid mesh). Returns a long DataFrame."""
    n, ns, _ = alpha.shape
    filt = np.zeros(ns)
    film = np.zeros(ns)
    rng = [_XorShift(i) for i in range(ns)]
    rows = []
    last = None
    for k in range(n):
        t = float(times[k])
        dt = 0.0 if last is None else max(0.0, t - last)
        g = 1.0 if last is None else 1 - math.exp(-dt / TAU_S)
        for i in range(ns):
            ok = ~solid[k, i]
            valid = int(ok.sum())
            if valid == 0:
                rows.append(dict(t=t, sensor=ids[i], oil_fraction=np.nan, air_fraction=np.nan,
                                 film_thickness_um=np.nan, capacitance_pf=np.nan, ideal_capacitance_pf=np.nan,
                                 signal=np.nan, wet=False, state='unavailable', solid_excluded_fraction=1.0))
                continue
            oil = float(np.clip(alpha[k, i][ok], 0, 1).mean())
            target_film = .00005 * oil
            if last is None:
                film[i] = target_film
            else:
                tau_f = .1 if target_film > film[i] else .7
                film[i] += (target_film - film[i]) * (1 - math.exp(-dt / tau_f))
            target = coax_capacitance(oil, film[i])
            filt[i] = target if last is None else filt[i] + g * (target - filt[i])
            u1, u2 = rng[i](), rng[i]()
            noise = math.sqrt(-2 * math.log(u1)) * math.cos(2 * math.pi * u2) * NOISE_PF
            meas = filt[i] + noise
            sig = min(1.0, max(0.0, (meas - C_DRY) / (C_WET - C_DRY)))
            rows.append(dict(t=t, sensor=ids[i], oil_fraction=oil, air_fraction=1 - oil,
                             film_thickness_um=film[i] * 1e6, capacitance_pf=meas, ideal_capacitance_pf=target,
                             signal=sig, wet=sig > .6, state='dry' if sig < .2 else ('mixed' if sig < .8 else 'wet'),
                             solid_excluded_fraction=1 - valid / 9))
        last = t
    return pd.DataFrame(rows)
