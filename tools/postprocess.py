#!/usr/bin/env python3
"""Post-process one interIBFoamOA2IS case into results/<name>/.

Outputs (all computed from the case's own solver output):
  log_history.csv            dt, Courant, interface Courant, p_rgh residuals, continuity, wall time
  inventory.csv              volIntegrate(alpha.oil) every step (L) + drift vs initial
  sensors.csv / .parquet     virtual OA2IS sensors S1..S8 every step (alpha, presence, U, p, p_rgh, chi)
  system_state.csv/.parquet  per-step state from the sensors (NORMAL / ALERTA / CRITICO) + CFD inventory
  field_metrics.csv          per written time: oil distribution metrics from the full field
  plots/*.png                inventory, numerics, sensor time series
  frames/frame_NNNN.png      cut-plane frames of every written time            (--frames)
  animation.mp4              H.264 video of the frames                         (--mp4)
  summary.json               headline numbers used in docs/VALIDACION.md

usage: python3 tools/postprocess.py <case_dir> [--out results/NAME] [--frames] [--mp4] [--every N]
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oa2is_cfd import ibpost as P          # noqa: E402
from oa2is_cfd.render import make_mp4, render_frame   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('--out')
    ap.add_argument('--frames', action='store_true')
    ap.add_argument('--mp4', action='store_true')
    ap.add_argument('--every', type=int, default=1)
    a = ap.parse_args()
    case = Path(a.case).resolve()
    sc = P.load_scenario(case)
    out = Path(a.out or Path(__file__).resolve().parents[1] / 'results' / case.name)
    (out / 'plots').mkdir(parents=True, exist_ok=True)

    log = P.parse_log(case)
    log.to_csv(out / 'log_history.csv', index=False)
    oil = P.oil_volume(case)
    oil['drift_L'] = oil.oil_L - oil.oil_L.iloc[0]
    oil['drift_pct'] = 100 * oil.drift_L / oil.oil_L.iloc[0]
    oil.to_csv(out / 'inventory.csv', index=False)
    sens = P.sensors(case, sc)
    sens.to_csv(out / 'sensors.csv', index=False)
    sens.to_parquet(out / 'sensors.parquet', index=False)
    state = P.system_state(sens, oil, sc.oil_litres)
    state.to_csv(out / 'system_state.csv', index=False)
    state.to_parquet(out / 'system_state.parquet', index=False)

    F = P.Fields(case)
    times = [t for t in F.times]
    met = pd.DataFrame([F.metrics(t) for t in times])
    met.to_csv(out / 'field_metrics.csv', index=False)

    # ---------------------------------------------------------------- plots
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    ax[0].plot(oil.t, oil.oil_L)
    ax[0].set(xlabel='t [s]', ylabel='aceite integrado [L]', title='Inventario (volIntegrate alpha)')
    ax[0].ticklabel_format(useOffset=False)
    ax[1].plot(oil.t, oil.drift_pct)
    ax[1].set(xlabel='t [s]', ylabel='deriva [%]', title='Conservacion de masa')
    fig.tight_layout()
    fig.savefig(out / 'plots/inventario.png', dpi=110)
    plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.6))
    ax[0].semilogy(log.t, log.dt)
    ax[0].set(xlabel='t [s]', ylabel='dt [s]', title='Paso de tiempo')
    ax[1].plot(log.t, log.Co, label='Co max')
    ax[1].plot(log.t, log.alphaCo, label='Co interfase max', alpha=.7)
    ax[1].legend(fontsize=8)
    ax[1].set(xlabel='t [s]', title='Courant')
    ax[2].semilogy(log.t, log.p_rgh_res0, label='p_rgh residuo inicial')
    ax[2].semilogy(log.t, log.p_rgh_res_final, label='p_rgh residuo final')
    ax[2].semilogy(log.t, log.cont_local, label='error continuidad local')
    ax[2].legend(fontsize=8)
    ax[2].set(xlabel='t [s]', title='Residuos')
    fig.tight_layout()
    fig.savefig(out / 'plots/numerica.png', dpi=110)
    plt.close(fig)
    fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    for sid, g in sens.groupby('sensor'):
        ax[0].plot(g.t, g.alpha, label=f'{sid}')
        ax[1].plot(g.t, g.U_mag, label=f'{sid}')
    ax[0].set(ylabel='alpha aceite', title=f'Sensores virtuales OA2IS - {sc.name}')
    ax[1].set(xlabel='t [s]', ylabel='|U| [m/s]')
    ax[0].legend(ncol=8, fontsize=7)
    fig.tight_layout()
    fig.savefig(out / 'plots/sensores.png', dpi=110)
    plt.close(fig)

    # --------------------------------------------------------------- frames
    if a.frames or a.mp4:
        fd = out / 'frames'
        fd.mkdir(exist_ok=True)
        for f in fd.glob('frame_*.png'):
            f.unlink()
        o = oil.set_index('t').oil_L
        for k, t in enumerate(times[::a.every]):
            render_frame(F, t, fd / f'frame_{k:04d}.png', sc, oil_L=float(np.interp(t, o.index, o.values)))
        if a.mp4:
            make_mp4(fd, out / 'animation.mp4', fps=8)

    s = dict(case=str(case), scenario=sc.name, mechanisms=sc.mechanisms, oil_litres_target=sc.oil_litres,
             cells=json.loads((case / 'mesh_stats.json').read_text()).get('total_cells'),
             t_end=float(log.t.iloc[-1]), steps=int(len(log)), wall_s=float(log.wall_s.iloc[-1]),
             dt_mean=float(log.dt.mean()), Co_max=float(log.Co.max()), Co_p95=float(log.Co.quantile(.95)),
             p_rgh_res_final_max=float(log.p_rgh_res_final.max()),
             oil_L_initial=float(oil.oil_L.iloc[0]), oil_L_final=float(oil.oil_L.iloc[-1]),
             oil_drift_max_pct=float(oil.drift_pct.abs().max()),
             states=state.state.value_counts().to_dict(),
             final_metrics=met.iloc[-1].to_dict(),
             sensor_alpha_mean_last_half=sens[sens.t > sens.t.max() / 2].groupby('sensor').alpha.mean().to_dict())
    (out / 'summary.json').write_text(json.dumps(s, indent=1, default=float))
    print(json.dumps({k: s[k] for k in ('scenario', 't_end', 'steps', 'oil_drift_max_pct', 'Co_max', 'states')}, default=float))


if __name__ == '__main__':
    main()
