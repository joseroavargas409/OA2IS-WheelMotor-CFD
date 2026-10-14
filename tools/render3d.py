#!/usr/bin/env python3
"""3D view of a CFD solution: oil free surface (alpha = 0.5 iso-surface,
blue), moving immersed bodies at their solved pose (chi = 0.5, orange) and a
transparent housing (cavity envelope). Frames for every written time
(or every N-th) and an optional MP4.

usage: xvfb-run -a python3 tools/render3d.py <case_dir> <out_dir> [--every 1] [--mp4]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oa2is_cfd.ibpost import Fields, load_scenario   # noqa: E402
from oa2is_cfd.render import make_mp4               # noqa: E402

pv.OFF_SCREEN = True


def frame(F, t, out, sc, env, fixed):
    m = F.mesh(t)
    pm = m.cell_data_to_point_data()
    p = pv.Plotter(off_screen=True, window_size=(1400, 900))
    p.set_background('white')
    p.add_mesh(env, color='#9aa3ad', opacity=0.08)
    p.add_mesh(env.extract_feature_edges(), color='#555555', line_width=1)
    for f in fixed:
        p.add_mesh(f, color='#b8bec6', opacity=0.25)
    if 'chi' in pm.point_data:
        bodies = pm.contour([0.5], scalars='chi')
        if bodies.n_points:
            bodies = bodies.clip(normal='x', origin=(0, 0, 0), invert=True)   # cut-away x > 0
        if bodies.n_points:
            p.add_mesh(bodies, color='#d9822b', smooth_shading=True)
    oil = pm.contour([0.5], scalars='alpha.oil')
    if oil.n_points:
        p.add_mesh(oil, color='#1f5fa8', opacity=0.85, smooth_shading=True)
    rpm = np.interp(t, [r[0] for r in sc.motor_rpm], [r[1] for r in sc.motor_rpm]) if sc.mechanisms != 'stopped' else 0
    acc = np.interp(t, [r[0] for r in sc.accel], [r[1] for r in sc.accel])
    p.add_text(f'{sc.name}   t = {t:.3f} s   motor {rpm:.0f} rpm   a_veh {acc:+.1f} m/s2\n'
               'azul: superficie libre del aceite (alpha = 0,5)   naranjo: cuerpos inmersos en su pose resuelta (corte x > 0)',
               font_size=10, color='black')
    p.camera_position = [(2.5, 0.9, -0.45), (0, -0.08, -1.36), (0, 1, 0)]
    p.screenshot(out)
    p.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('out')
    ap.add_argument('--every', type=int, default=1)
    ap.add_argument('--mp4', action='store_true')
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    F = Fields(a.case)
    sc = load_scenario(a.case)
    g = ROOT / 'geometry' / 'cfd'
    env = pv.read(g / 'cavity_envelope.stl')
    fixed = [pv.read(g / f'{n}.stl').clip(normal='x', origin=(0, 0, 0), invert=True)
             for n in ('wall_spindle', 'wall_ringHigh', 'wall_ringLow')]
    for k, t in enumerate(F.times[::a.every]):
        frame(F, t, str(out / f'frame_{k:04d}.png'), sc, env, fixed)
    if a.mp4:
        make_mp4(out, out.parent / 'animation_3d.mp4', fps=8)


if __name__ == '__main__':
    main()
