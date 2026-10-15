#!/usr/bin/env python3
"""Cut-plane frames + MP4 with the oil/air sensors drawn in their state.

usage: python3 tools/render_with_sensors.py <oil_sensors.csv> <out_dir> <case> [<case> ...] [--plugs plug_sensors.csv]
--plugs: also draw the hub-cover plug sensors (positions rotate with the hub) in an
extra panel cut parallel to the hub cover.
Cases are concatenated in time (e.g. the original run and its continuation):
for each case only the written times before the next case starts are used.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oa2is_cfd.ibpost import Fields, load_scenario, oil_volume   # noqa: E402
from oa2is_cfd.render import make_mp4, render_frame             # noqa: E402


def main():
    args = sys.argv[1:]
    plugs = None
    if '--plugs' in args:
        i = args.index('--plugs')
        plugs = pd.read_csv(args[i + 1])
        del args[i:i + 2]
    csv, out, cases = Path(args[0]), Path(args[1]), [Path(c) for c in args[2:]]
    out.mkdir(parents=True, exist_ok=True)
    for f in out.glob('frame_*.png'):
        f.unlink()
    df = pd.read_csv(csv)
    pos = df.groupby('sensor')[['x', 'y', 'z']].first()
    ts = np.sort(df.t.unique())
    F = [Fields(c) for c in cases]
    starts = [min(t for t in f.times if t > 0) if len(f.times) > 1 else 0 for f in F]
    sc = load_scenario(cases[-1])
    oils = pd.concat([oil_volume(c) for c in cases]).drop_duplicates('t').sort_values('t')
    k = 0
    for i, f in enumerate(F):
        nxt = starts[i + 1] - 0.0099 if i + 1 < len(F) else 1e9
        for t in f.times:
            if t > nxt or (i > 0 and t < starts[i] - 1e-9):
                continue
            tk = ts[ts <= t + 1e-9].max() if (ts <= t + 1e-9).any() else ts.min()
            st = df[np.isclose(df.t, tk)].set_index('sensor').state
            sensors = [(s, tuple(pos.loc[s]), st.get(s, 'unavailable')) for s in pos.index]
            extra = {}
            if plugs is not None:
                pt = np.sort(plugs.t.unique())
                tp = pt[np.argmin(np.abs(pt - t))]
                g = plugs[plugs.t == tp]
                sensors += [(r.sensor, (r.x, r.y, r.z), r.state) for r in g.itertuples()]
                extra = dict(cover_z=float(g.z.iloc[0]), hub_angle_deg=float(g.hub_angle_deg.iloc[0]))
            render_frame(f, t, out / f'frame_{k:04d}.png', sc,
                         oil_L=float(np.interp(t, oils.t, oils.oil_L)), sensors=sensors, **extra)
            k += 1
    make_mp4(out, out.parent / 'animation_sensores.mp4', fps=8)
    print(k, 'frames')


if __name__ == '__main__':
    main()
