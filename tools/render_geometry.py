#!/usr/bin/env python3
"""3D picture of the processed CFD geometry (geometry/cfd/*.stl).
usage: xvfb-run -a python3 tools/render_geometry.py docs/img/geometria_cfd.png"""
import glob
import os
import sys
from pathlib import Path

import pyvista as pv

pv.OFF_SCREEN=True
d = str(Path(__file__).resolve().parents[1] / 'geometry' / 'cfd')
p=pv.Plotter(off_screen=True, window_size=(1600,1000))
col={'wall_spindle':'#b0b4ba','wall_seal':'#b0b4ba','wall_ringHigh':'#8c96a3','wall_ringLow':'#8c96a3','static_ringHigh':'#7d8794','static_ringLow':'#7d8794',
     'wall_hubWall':'#d9c27a','wall_carrierLowRing':'#4c9f70','wall_intermediate':'#c0504d','wall_shaft':'#c0504d'}
for f in sorted(glob.glob(d+'/*.stl')):
    n=os.path.basename(f)[:-4]
    if n=='cavity_envelope': continue
    m=pv.read(f)
    if n.startswith('wall_hubWall'):
        m=m.clip(normal='x', origin=(0,0,0)); 
    c=col.get(n, '#e07b39' if 'planet' in n else ('#c0504d' if 'sun' in n else ('#4c9f70' if 'carrier' in n else '#999999')))
    if not ('planet' in n or 'sun' in n):
        m=m.clip(normal='x', origin=(0,0,0))
    p.add_mesh(m, color=c, smooth_shading=False, opacity=1.0)
env=pv.read(d+'/cavity_envelope.stl')
p.add_mesh(env.extract_feature_edges(), color='black', line_width=1)
p.camera_position=[(2.2,1.0,0.3),(0,0,-1.35),(0,1,0)]
p.add_text('Geometria CFD procesada (t = 0): engranajes regenerados, porta-satelites, corona, spindle y cubo\ncorte x>0 en carcasa/spindle/coronas', font_size=10)
p.screenshot(sys.argv[1])
