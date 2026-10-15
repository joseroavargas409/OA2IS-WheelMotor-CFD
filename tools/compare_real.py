#!/usr/bin/env python3
"""Comparison of the real-cycle scenarios (06-09) from their results/ folders.

Writes results/validation/real_cycle_summary.csv and
results/validation/real_cycle.png (oil pool tilt, oil lifted above the axis,
gear wetting, and the level the hub-cover plug ring would infer).
usage: python3 tools/compare_real.py
"""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402
import pandas as pd               # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / 'results'
CASES = ['06_real_reposo_nivel', '07_real_acarreo', '08_real_subida_cargado', '09_real_bajada_retardo']
T0 = 0.15                                     # statistics window start (after spin-up and first transient)


def main():
    rows = []
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True)
    for c in CASES:
        d = RES / c
        if not (d / 'summary.json').exists():
            continue
        s = json.loads((d / 'summary.json').read_text())
        fm = pd.read_csv(d / 'field_metrics.csv')
        ring = pd.read_csv(d / 'plug_ring.csv')
        plugs = pd.read_csv(d / 'plug_sensors_summary.csv').set_index('sensor')
        w = fm[fm.t >= T0]
        # tilt of the oil pool: angle of the oil centroid from straight down (+ = forward)
        fm['pool_angle_deg'] = [math.degrees(math.atan2(x, -y)) for x, y in zip(fm.oil_centroid_x, fm.oil_centroid_y)]
        wr = ring[ring.t >= T0]
        lh = pd.read_csv(d / 'log_history.csv')
        seg = (lh.wall_s.diff() < 0).cumsum()           # the solver clock restarts when a run is continued
        wall = lh.groupby(seg).wall_s.max().sum()
        rows.append(dict(
            case=c, t_end=s['t_end'], steps=s['steps'], wall_min=round(wall / 60, 1),
            dt_mean_ms=round(s['dt_mean'] * 1e3, 4), Co_max=s['Co_max'], Co_p95=round(s['Co_p95'], 3),
            oil_L=round(s['oil_L_initial'], 2), oil_drift_max_pct=s['oil_drift_max_pct'],
            oil_mean_speed=round(w.oil_mean_speed.mean(), 3), oil_above_axis_pct=round(100 * w.oil_above_axis_frac.mean(), 2),
            wetted_gear_frac=round(w.wetted_gear_frac.mean(), 3), wetted_high=round(w.wetted_frac_stageHigh.mean(), 3),
            wetted_low=round(w.wetted_frac_stageLow.mean(), 3), cover_pocket_L=round(w.oil_L_cover_pocket.mean(), 2),
            oil_in_bodies_L=round(w.oil_excess_in_bodies_L.mean(), 2),
            pool_angle_deg=round(fm[fm.t >= T0].pool_angle_deg.mean(), 2),
            ring_level_y=round(wr.level_y_from_ring.mean(), 4), ring_wet_arc_deg=round(wr.wet_arc_deg.mean(), 1),
            P6_wet_mixed=round(plugs.loc['P6', 'frac_wet'] + plugs.loc['P6', 'frac_mixed'], 3),
            P7_wet_mixed=round(plugs.loc['P7', 'frac_wet'] + plugs.loc['P7', 'frac_mixed'], 3)))
        axs[0, 0].plot(fm.t, fm.pool_angle_deg, label=c)
        axs[0, 1].plot(fm.t, 100 * fm.oil_above_axis_frac, label=c)
        axs[1, 0].plot(fm.t, fm.wetted_gear_frac, label=c)
        axs[1, 1].plot(ring.t, ring.level_y_from_ring, label=c)
    axs[0, 0].set_ylabel('inclinacion del pozo de aceite [deg]\n(+ hacia adelante)')
    axs[0, 1].set_ylabel('aceite sobre el eje [%]')
    axs[1, 0].set_ylabel('fraccion mojada de engranajes')
    axs[1, 1].set_ylabel('nivel inferido por el anillo de tapones [m]')
    axs[1, 1].axhline(0.270 * math.sin(math.radians(240)), color='k', ls=':', lw=.8)
    for ax in axs[1]:
        ax.set_xlabel('t [s]')
    axs[0, 0].legend(fontsize=7)
    fig.suptitle('Escenarios del ciclo real (malla ib_test): 06 reposo, 07 acarreo, 08 subida cargado, 09 bajada con retardo')
    fig.tight_layout()
    out = RES / 'validation'
    out.mkdir(exist_ok=True)
    fig.savefig(out / 'real_cycle.png', dpi=110)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'real_cycle_summary.csv', index=False)
    print(df.T.to_string())


if __name__ == '__main__':
    main()
