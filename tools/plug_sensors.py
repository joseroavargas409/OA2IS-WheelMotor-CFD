#!/usr/bin/env python3
"""Oil/air sensors at the hub-cover magnetic plugs (rotating with the hub).

The case records a ring of probes, fixed in space, on the plug circle of the
hub cover (sensor_layout.json -> rotating[], function object ring_<group>),
each ring point with the 9-point electrode volume of the capacitive sensor of
komatsu-oil-sensors-v4-1 (electrode along the wheel axis, as a sensor screwed
into a plug port). The hub turns, so at time t plug k is at

    phi_k(t) = phi_k0 + theta_hub(t),   theta_hub(t) = int_0^t ratio_out * omega_rotor dt

(the same kinematics the solver applies to the hub walls). The ring values are
interpolated linearly in angle at phi_k(t) and fed to oa2is_cfd.oil_sensor.run().

Writes to results/<case>/:
  plug_sensors.csv / .parquet   every time step x plug: angle, position, oil fraction,
                                capacitance, signal, state (same columns as oil_sensors.csv)
  plug_sensors_summary.csv      per plug: time fractions wet/mixed/dry, mean oil fraction
  plug_ring.csv                 oil fraction along the whole plug circle (<= 400 times), wetted arc and
                                the pool level it implies (what one plug sees over a full wheel turn)
  plots/plug_sensors.png        capacitance and states of P1-P8
  plots/plug_ring.png           oil fraction around the plug circle vs time, with the plug tracks

usage: python3 tools/plug_sensors.py <case_dir> [--out results/NAME] [--t0 0.05]
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oa2is_cfd import oil_sensor                 # noqa: E402
from oa2is_cfd.ibpost import load_scenario        # noqa: E402
sys.path.insert(0, str(ROOT / 'tools'))
from oil_sensors import read_probes               # noqa: E402


def part_angle(case: Path, part: str, t: np.ndarray) -> np.ndarray:
    """Rotation angle [rad] of a part (ratio key of scenario.json: inp/mid/out) at times t."""
    meta = json.loads((case / 'scenario.json').read_text())
    k = meta['ratios'][part]
    W = np.array(load_scenario(case).rotor_omega_table(), float)
    tf = np.linspace(0, max(t.max(), 1e-6), 200001)
    w = k * np.interp(tf, W[:, 0], W[:, 1])
    th = np.r_[0, np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(tf))]
    return np.interp(t, tf, th)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('--out')
    ap.add_argument('--t0', type=float, default=0.05)
    a = ap.parse_args()
    case = Path(a.case).resolve()
    out = Path(a.out or ROOT / 'results' / case.name)
    (out / 'plots').mkdir(parents=True, exist_ok=True)
    L = json.loads((case / 'sensor_layout.json').read_text())
    import pandas as pd
    for g in L.get('rotating', []):
        ang = np.array(g['ring_angles'])
        n, no = len(ang), len(g['offsets'])
        T, AL, CH = [], [], []
        for d in sorted(case.glob(f"postProcessing/{g['function_object']}/*"), key=lambda p: float(p.name)):
            t, al = read_probes(d / 'alpha.oil', n * no)
            t2, ch = read_probes(d / 'chi', n * no)
            m = min(len(t), len(t2))
            T.append(t[:m]); AL.append(al[:m]); CH.append(ch[:m])
        t = np.concatenate(T)
        t, idx = np.unique(t, return_index=True)
        al = np.concatenate(AL)[idx].reshape(len(t), n, no)
        ch = np.concatenate(CH)[idx].reshape(len(t), n, no)
        th = part_angle(case, g['part'], t)
        plugs = g['plugs']
        ids = [p['id'] for p in plugs]
        A = np.zeros((len(t), len(plugs), no))
        S = np.zeros((len(t), len(plugs), no), bool)
        phi = np.zeros((len(t), len(plugs)))
        rows = np.arange(len(t))
        for k, p in enumerate(plugs):
            phi[:, k] = np.radians(p['angle_deg']) + th
            f = np.mod(phi[:, k], 2 * np.pi) / (2 * np.pi) * n
            i0 = np.floor(f).astype(int) % n
            i1 = (i0 + 1) % n
            w = (f - np.floor(f))[:, None]
            a0, a1 = al[rows, i0], al[rows, i1]
            c0, c1 = ch[rows, i0], ch[rows, i1]
            A[:, k] = np.nan_to_num((1 - w) * a0 + w * a1)
            S[:, k] = np.isnan(a0) | np.isnan(a1) | (np.nan_to_num((1 - w) * c0 + w * c1, nan=1.0) > 0.5)
        df = oil_sensor.run(t, A, S, ids)
        name = {p['id']: p['name'] for p in plugs}
        k_of = {p['id']: k for k, p in enumerate(plugs)}
        ti = {v: i for i, v in enumerate(t)}
        r, z = g['radius'], g['z']
        ph = np.array([phi[ti[tt], k_of[s]] for tt, s in zip(df.t, df.sensor)])
        df['name'] = df.sensor.map(name)
        df['angle_rest_deg'] = df.sensor.map({p['id']: p['angle_deg'] for p in plugs})
        df['hub_angle_deg'] = np.degrees(np.array([th[ti[tt]] for tt in df.t]))
        df['angle_deg'] = np.mod(np.degrees(ph), 360)
        df['x'], df['y'], df['z'] = r * np.cos(ph), r * np.sin(ph), z
        df.to_csv(out / 'plug_sensors.csv', index=False)
        df.to_parquet(out / 'plug_sensors.parquet', index=False)
        w = df[df.t >= a.t0]
        summ = w.groupby('sensor', sort=False).agg(
            name=('name', 'first'), angle_rest_deg=('angle_rest_deg', 'first'),
            y_start=('y', 'first'), y_end=('y', 'last'),
            frac_wet=('state', lambda s: (s == 'wet').mean()), frac_mixed=('state', lambda s: (s == 'mixed').mean()),
            frac_dry=('state', lambda s: (s == 'dry').mean()), frac_unavailable=('state', lambda s: (s == 'unavailable').mean()),
            oil_fraction_mean=('oil_fraction', 'mean'), capacitance_pf_mean=('capacitance_pf', 'mean'),
            capacitance_pf_std=('capacitance_pf', 'std')).reset_index()
        summ.to_csv(out / 'plug_sensors_summary.csv', index=False)
        # whole ring (every written step is enough: 1 row per probe-file time)
        ring = np.nanmean(np.where(ch > 0.5, np.nan, al), axis=2)
        step = max(1, len(t) // 400)
        rg = pd.DataFrame(ring[::step], columns=[f'{d:.0f}' for d in np.degrees(ang)])
        rg.insert(0, 't', t[::step])
        rg.insert(1, 'hub_angle_deg', np.degrees(th[::step]))
        # level inferred from the ring: a plug sweeps the whole circle in one wheel turn,
        # so the wetted arc 2*beta (centred at the bottom) gives the pool height on the
        # plug circle: y_level = -r cos(beta)  (beta = pi * wetted fraction of the circle)
        wf = np.mean(ring[::step] > 0.5, axis=1)
        rg.insert(2, 'wet_arc_deg', 360 * wf)
        rg.insert(3, 'level_y_from_ring', -r * np.cos(np.pi * wf))
        rg.to_csv(out / 'plug_ring.csv', index=False)
        lv = rg[rg.t >= a.t0] if (rg.t >= a.t0).any() else rg
        print('ring level estimate [m]: start %.4f  mean(t>=t0) %.4f  end %.4f  (plug level y = %.4f)' %
              (rg.level_y_from_ring.iloc[0], lv.level_y_from_ring.mean(), rg.level_y_from_ring.iloc[-1],
               r * np.sin(np.radians(240))))
        # plots
        fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
        for s in ids:
            gg = df[df.sensor == s]
            axs[0].plot(gg.t, gg.capacitance_pf, lw=.9, label=f"{s} ({gg.angle_rest_deg.iloc[0]} deg)")
            axs[1].plot(gg.t, gg.oil_fraction, lw=.9, label=s)
        axs[0].axhline(oil_sensor.C_DRY, color='k', lw=.6, ls=':')
        axs[0].axhline(oil_sensor.C_WET, color='k', lw=.6, ls=':')
        axs[0].set_ylabel('capacitancia [pF]')
        axs[0].legend(ncol=4, fontsize=7)
        axs[1].set_ylabel('fraccion de aceite')
        axs[1].set_xlabel('t [s]')
        axs[0].set_title(f'{case.name}: sensores en los tapones de la tapa de maza (giran con la maza; '
                         f'nivel: P6 240 deg / P7 300 deg en reposo)')
        fig.tight_layout()
        fig.savefig(out / 'plots/plug_sensors.png', dpi=110)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(10, 5.5))
        dg = np.degrees(ang)
        te = np.r_[t[::step], t[-1] + 1e-4]
        de = np.r_[dg, 360.0]
        pc = ax.pcolormesh(te, de, ring[::step].T, cmap='Blues', vmin=0, vmax=1, shading='flat')
        for k, s in enumerate(ids):
            ax.plot(t, np.mod(np.degrees(phi[:, k]), 360), '.', ms=1, color='crimson')
            ax.text(t[-1], np.mod(np.degrees(phi[-1, k]), 360), ' ' + s, fontsize=7, color='crimson', va='center')
        ax.set_ylabel('angulo en la tapa (marco fijo) [deg]  (270 = abajo)')
        ax.set_xlabel('t [s]')
        ax.set_title('Aceite a lo largo del circulo de tapones (r = %.3f m) y trayectoria de cada tapon' % r)
        fig.colorbar(pc, ax=ax, label='fraccion de aceite')
        fig.tight_layout()
        fig.savefig(out / 'plots/plug_ring.png', dpi=110)
        plt.close(fig)
        print(summ.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
