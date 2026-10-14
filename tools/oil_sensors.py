#!/usr/bin/env python3
"""Oil/air sensors (capacitive model of komatsu-oil-sensors-v4-1) on a CFD case.

The case must have been run with the 9-point electrode stencil of every sensor
as probe locations (sensor_layout.json in the case, written when the case is
prepared). Reads postProcessing/sensors/*/{alpha.oil,chi}, applies
oa2is_cfd.oil_sensor.run() and writes to results/<case>/:

  oil_sensors.csv / .parquet   every time step x sensor: oil/air fraction, film,
                               capacitance [pF], signal, wet flag, state
  oil_sensors_summary.csv      per sensor: time fractions dry/mixed/wet/unavailable,
                               mean oil fraction, mean and std capacitance
  plots/oil_sensors_signal.png, plots/oil_sensors_states.png

usage: python3 tools/oil_sensors.py <case_dir> [--out results/NAME] [--t0 0.05] [--history CASE]

--history CASE: a case run earlier on the SAME mesh (e.g. the one this case was
restarted from); its written fields (every write interval) are sampled at the
same 180 points and prepended, so the sensors cover the whole simulated time.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oa2is_cfd import oil_sensor   # noqa: E402


def read_probes(path: Path, n_probes: int):
    """Return times [n], values [n, n_probes] (NaN for probes not found)."""
    lines = path.read_text().splitlines()
    heads = [l for l in lines if re.match(r'#\s*Probe \d+ \(', l)]
    found = [i for i, l in enumerate(heads) if 'Not Found' not in l]
    rows = []
    for l in lines:
        if not l.strip() or l.startswith('#'):
            continue
        rows.append([float(x) for x in l.split()])
    a = np.array(rows)
    vals = np.full((len(a), n_probes), np.nan)
    if len(a):
        if a.shape[1] == 1 + len(heads):          # OpenFOAM writes every probe, flags missing ones in the header
            vals[:, :len(heads)] = a[:, 1:]
            missing = [i for i, l in enumerate(heads) if 'Not Found' in l]
            vals[:, missing] = np.nan
        else:                                     # only the found probes are written
            vals[:, found] = a[:, 1:1 + len(found)]
    return (a[:, 0] if len(a) else np.zeros(0)), vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('--out')
    ap.add_argument('--t0', type=float, default=0.0, help='ignore samples before t0 in the summary')
    ap.add_argument('--history', help='earlier case on the same mesh: sample its written fields before this case starts')
    a = ap.parse_args()
    case = Path(a.case).resolve()
    out = Path(a.out or ROOT / 'results' / case.name)
    (out / 'plots').mkdir(parents=True, exist_ok=True)
    L = json.loads((case / 'sensor_layout.json').read_text())
    S, OFF = L['sensors'], L['offsets']
    ns, no = len(S), len(OFF)
    T, AL, CH = [], [], []
    for d in sorted(case.glob('postProcessing/sensors/*'), key=lambda p: float(p.name)):
        t, al = read_probes(d / 'alpha.oil', ns * no)
        t2, ch = read_probes(d / 'chi', ns * no)
        n = min(len(t), len(t2))
        T.append(t[:n]); AL.append(al[:n]); CH.append(ch[:n])
    if a.history:
        from oa2is_cfd.ibpost import Fields
        F = Fields(a.history)
        pts = np.array([[s['position'][i] + o[i] for i in range(3)] for s in S for o in OFF])
        m0 = F.mesh(F.times[0])
        cells = np.array([m0.find_containing_cell(p) for p in pts])
        t_start = min(x[0] for x in T if len(x))
        for tt in F.times:
            if tt >= t_start - 1e-9:
                continue
            m = F.mesh(tt)
            av = np.full(len(pts), np.nan)
            cv = np.full(len(pts), np.nan)
            ok = cells >= 0
            av[ok] = m.cell_data['alpha.oil'][cells[ok]]
            cv[ok] = m.cell_data['chi'][cells[ok]] if 'chi' in m.cell_data else 0.0
            T.insert(0, np.array([tt])); AL.insert(0, av[None]); CH.insert(0, cv[None])
    t = np.concatenate(T)
    al = np.concatenate(AL)
    ch = np.concatenate(CH)
    t, idx = np.unique(t, return_index=True)
    al, ch = al[idx], ch[idx]
    alpha = al.reshape(len(t), ns, no)
    solid = np.isnan(alpha) | (np.nan_to_num(ch.reshape(len(t), ns, no), nan=1.0) > 0.5)
    ids = [s['id'] for s in S]
    df = oil_sensor.run(t, np.nan_to_num(alpha), solid, ids)
    names = {s['id']: s['name'] for s in S}
    pos = {s['id']: s['position'] for s in S}
    df['name'] = df.sensor.map(names)
    for k, c in enumerate('xyz'):
        df[c] = df.sensor.map(lambda s: pos[s][k])
    df.to_csv(out / 'oil_sensors.csv', index=False)
    df.to_parquet(out / 'oil_sensors.parquet', index=False)
    w = df[df.t >= a.t0]
    summ = w.groupby('sensor', sort=False).agg(
        name=('name', 'first'),
        frac_wet=('state', lambda s: (s == 'wet').mean()), frac_mixed=('state', lambda s: (s == 'mixed').mean()),
        frac_dry=('state', lambda s: (s == 'dry').mean()), frac_unavailable=('state', lambda s: (s == 'unavailable').mean()),
        oil_fraction_mean=('oil_fraction', 'mean'), capacitance_pf_mean=('capacitance_pf', 'mean'),
        capacitance_pf_std=('capacitance_pf', 'std'), film_um_mean=('film_thickness_um', 'mean')).reset_index()
    summ.to_csv(out / 'oil_sensors_summary.csv', index=False)
    # plots
    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for ax, grp in zip(axs, (ids[:8], ids[8:])):
        for s in grp:
            g = df[df.sensor == s]
            ax.plot(g.t, g.capacitance_pf, lw=.8, label=s)
        ax.axhline(oil_sensor.C_DRY, color='k', lw=.6, ls=':')
        ax.axhline(oil_sensor.C_WET, color='k', lw=.6, ls=':')
        ax.set_ylabel('capacitancia [pF]')
        ax.legend(ncol=6, fontsize=7)
    axs[0].set_title('Sensores capacitivos (modelo komatsu-oil-sensors-v4-1) sobre la solucion CFD: '
                     'S1-S8 montajes OA2IS, V1-V12 sondas V4.1')
    axs[1].set_xlabel('t [s]')
    fig.tight_layout()
    fig.savefig(out / 'plots/oil_sensors_signal.png', dpi=110)
    plt.close(fig)
    code = {'dry': 0, 'mixed': 1, 'wet': 2, 'unavailable': -1}
    M = df.pivot_table(index='sensor', columns='t', values='state', aggfunc='first').reindex(ids)
    Z = M.apply(lambda col: col.map(code)).values.astype(float)
    fig, ax = plt.subplots(figsize=(12, 5))
    cmap = matplotlib.colors.ListedColormap(['#9a9a9a', '#f4e3c1', '#8fb8de', '#08306b'])
    tt = M.columns.values.astype(float)
    edges = np.r_[tt, tt[-1] + (tt[-1] - tt[-2] if len(tt) > 1 else 1e-3)]
    ax.pcolormesh(edges, np.arange(len(ids) + 1) - .5, Z, cmap=cmap, vmin=-1.5, vmax=2.5, shading='flat')
    ax.invert_yaxis()
    ax.set_yticks(range(len(ids)))
    ax.set_yticklabels([f'{s}  {names[s][:38]}' for s in ids], fontsize=7)
    ax.set_xlabel('t [s]')
    ax.set_title('Estado de cada sensor: gris no disponible (dentro de un solido), crema seco, celeste mixto, azul mojado')
    fig.tight_layout()
    fig.savefig(out / 'plots/oil_sensors_states.png', dpi=110)
    plt.close(fig)
    print(summ.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
