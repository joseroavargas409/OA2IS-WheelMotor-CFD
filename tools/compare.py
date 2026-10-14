#!/usr/bin/env python3
"""Cross-case validation figures and numbers (results/validation/).

Reads the per-case outputs written by tools/postprocess.py (results/<case>/)
and produces:
  stopped_vs_moving.png      01_reposo vs 02_operacion_normal (same mesh, same oil)
  mesh_sensitivity.png       ib_test vs ib_mid, nominal case, common time window
  dt_sensitivity.png         Co_max 1.0 vs 0.5, nominal case, common time window
  inertial_response.png      oil-pool angle vs g_eff angle (03 inclined, 04 transient)
  inventory_all.png          oil inventory drift of every case
  validation_summary.json    every number quoted in docs/VALIDACION.md

usage: python3 tools/compare.py --test results/test_02_operacion_normal ... (see defaults)
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oa2is_cfd.scenario import Scenario   # noqa: E402

R = ROOT / 'results'
METRICS = [('wetted_gear_frac', 'fraccion mojada zona engranajes [-]'),
           ('oil_above_axis_frac', 'aceite sobre el eje [-]'),
           ('oil_mean_speed', 'velocidad media del aceite [m/s]'),
           ('oil_centroid_y', 'centroide y del aceite [m]')]


def load(name):
    d = R / name
    if not (d / 'summary.json').exists():
        return None
    return dict(name=name, dir=d, summary=json.loads((d / 'summary.json').read_text()),
                met=pd.read_csv(d / 'field_metrics.csv'), inv=pd.read_csv(d / 'inventory.csv'),
                log=pd.read_csv(d / 'log_history.csv'), sens=pd.read_csv(d / 'sensors.csv'))


def pool_angle(met):
    """Angle [deg] of the oil centroid around the wheel axis, measured from -y
    towards +x (0 = oil pool at the bottom)."""
    return np.degrees(np.arctan2(met.oil_centroid_x, -met.oil_centroid_y))


def compare_pair(a, b, la, lb, fname, title, t0=0.0, t1=None):
    t1 = t1 or min(a['met'].t.max(), b['met'].t.max())
    fig, ax = plt.subplots(1, len(METRICS), figsize=(16, 3.6))
    out = {}
    for k, (m, lab) in enumerate(METRICS):
        for c, l in ((a, la), (b, lb)):
            d = c['met'][(c['met'].t >= t0) & (c['met'].t <= t1 + 1e-9)]
            ax[k].plot(d.t, d[m], marker='o', ms=3, label=l)
        ax[k].set(xlabel='t [s]', title=lab)
        da = a['met'].set_index('t')[m]
        db = b['met'].set_index('t')[m]
        common = sorted(set(np.round(da.index, 6)) & set(np.round(db.index, 6)))
        common = [t for t in common if t0 <= t <= t1 + 1e-9]
        if common:
            va = np.array([da[np.isclose(da.index, t)].iloc[0] for t in common])
            vb = np.array([db[np.isclose(db.index, t)].iloc[0] for t in common])
            win = np.array(common) >= max(t0, 0.05)       # after the numerical spin-up
            if win.any():
                out[m] = dict(mean_a=float(va[win].mean()), mean_b=float(vb[win].mean()),
                              mean_abs_diff=float(np.abs(va[win] - vb[win]).mean()),
                              rel_diff_of_means=float(abs(va[win].mean() - vb[win].mean()) /
                                                      max(abs(va[win].mean()), 1e-12)))
    ax[0].legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(R / 'validation' / fname, dpi=110)
    plt.close(fig)
    return dict(window=[t0, t1], metrics=out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prefix', default='test_')
    a = ap.parse_args()
    (R / 'validation').mkdir(parents=True, exist_ok=True)
    P = a.prefix
    cases = {n: load(n) for n in [P + '01_reposo', P + '02_operacion_normal', P + '03_operacion_inclinada',
                                   P + '04_transitorio', P + '05_bajo_inventario', 'dt_02_Co05',
                                   'mid_02_operacion_normal']}
    S = {}
    c01, c02 = cases[P + '01_reposo'], cases[P + '02_operacion_normal']
    if c01 and c02:
        S['stopped_vs_moving'] = compare_pair(c01, c02, 'detenido (01)', 'en movimiento (02)', 'stopped_vs_moving.png',
                                              'Mecanismos detenidos vs en movimiento (malla ib_test, 95 L)')
        for lab, c in (('stopped', c01), ('moving', c02)):
            sm = c['sens'][c['sens'].t > 0.05].groupby('sensor').agg(alpha=('alpha', 'mean'), wet=('oil_present', 'mean'),
                                                                     U=('U_mag', 'mean'))
            S.setdefault('sensors_stopped_vs_moving', {})[lab] = sm.round(4).to_dict()
        S['stopped_max_speed'] = float(c01['met'].max_speed_fluid.max())
    cm = cases['mid_02_operacion_normal']
    if c02 and cm:
        S['mesh'] = compare_pair(c02, cm, 'ib_test 146k', 'ib_mid 477k', 'mesh_sensitivity.png',
                                 'Sensibilidad de malla, caso nominal', t1=min(c02['met'].t.max(), cm['met'].t.max()))
    cd = cases['dt_02_Co05']
    if c02 and cd:
        S['dt'] = compare_pair(c02, cd, 'Co_max 1.0', 'Co_max 0.5', 'dt_sensitivity.png',
                               'Sensibilidad al paso de tiempo, caso nominal, malla ib_test')
    # inertial response: pool angle vs g_eff angle
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
    for k, (n, ref) in enumerate(((P + '03_operacion_inclinada', P + '02_operacion_normal'),
                                  (P + '04_transitorio', P + '02_operacion_normal'))):
        c = cases[n]
        if not c:
            continue
        sc = Scenario.load(ROOT / 'scenarios' / f'{n[len(P):]}.yaml')
        ang = pool_angle(c['met'])
        ax[k].plot(c['met'].t, ang, 'o-', ms=3, label='angulo del pozo de aceite (CFD)')
        if cases[ref]:
            ax[k].plot(cases[ref]['met'].t, pool_angle(cases[ref]['met']), 's--', ms=3, alpha=.6, label='nominal 02 (CFD)')
        tt = np.linspace(0, c['met'].t.max(), 200)
        g = np.array([sc.g_eff(t) for t in tt])
        gang = np.degrees(np.arctan2(g[:, 0], -g[:, 1]))
        ax[k].plot(tt, gang, 'k-', lw=1, label='direccion de g_eff (cuasi-estatico)')
        ax[k].set(xlabel='t [s]', ylabel='angulo desde -y hacia +x [deg]', title=n)
        ax[k].legend(fontsize=7)
        S.setdefault('inertial', {})[n] = dict(pool_angle_deg=dict(zip(map(float, c['met'].t.round(4)), map(float, ang.round(2)))),
                                              g_eff_angle_end_deg=float(gang[-1]),
                                              oil_centroid_z=dict(zip(map(float, c['met'].t.round(4)),
                                                                      map(float, c['met'].oil_centroid_z.round(4)))))
    fig.tight_layout()
    fig.savefig(R / 'validation' / 'inertial_response.png', dpi=110)
    plt.close(fig)
    # inventory + numerics table
    fig, ax = plt.subplots(figsize=(8, 3.6))
    tab = []
    for n, c in cases.items():
        if not c:
            continue
        ax.plot(c['inv'].t, c['inv'].drift_pct, label=n)
        s = c['summary']
        tab.append(dict(case=n, cells=s['cells'], t_end=s['t_end'], steps=s['steps'], wall_h=s['wall_s'] / 3600,
                        dt_mean=s['dt_mean'], Co_max=s['Co_max'], Co_p95=s['Co_p95'],
                        p_rgh_final_res_max=s['p_rgh_res_final_max'], oil_L0=s['oil_L_initial'],
                        oil_L_end=s['oil_L_final'], drift_max_pct=s['oil_drift_max_pct'], states=s['states']))
    ax.set(xlabel='t [s]', ylabel='deriva de inventario [%]', title='Conservacion del aceite, todos los casos')
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(R / 'validation' / 'inventory_all.png', dpi=110)
    plt.close(fig)
    S['runs'] = tab
    (R / 'validation' / 'validation_summary.json').write_text(json.dumps(S, indent=1, default=float))
    print(pd.DataFrame(tab).drop(columns=['states']).to_string(index=False, float_format=lambda v: f'{v:.4g}'))


if __name__ == '__main__':
    main()
