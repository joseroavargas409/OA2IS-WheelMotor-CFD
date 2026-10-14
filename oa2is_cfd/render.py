"""Frames and animations rendered from the CFD fields (no synthetic data).

Each frame shows three cuts through the solved field at one written time:
  (a) plane z = high-stage mid plane, (b) z = low-stage mid plane,
  (c) vertical plane x = 0 (axial section, y up).
Colours: oil volume fraction alpha.oil (white = air, dark blue = oil),
moving immersed bodies (chi > 0.5) dark grey, fixed solids (outside the fluid
mesh) light steel grey. A header gives time, motor speed, vehicle acceleration
and the integrated oil inventory of that frame.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402
from matplotlib.collections import PolyCollection    # noqa: E402
from matplotlib.colors import LinearSegmentedColormap   # noqa: E402

from .ibpost import STAGE_Z, Fields, load_scenario   # noqa: E402

OIL_CMAP = LinearSegmentedColormap.from_list('oil', ['#ffffff', '#9ecae1', '#2171b5', '#08306b'])
BODY = np.array([0.30, 0.30, 0.32, 1.0])
STEEL = '#c9ccd1'


def _polys(sl, axes):
    """Polygons of a pyvista slice (PolyData) projected on `axes` (i, j)."""
    faces = sl.faces
    pts = sl.points[:, axes]
    out, i = [], 0
    while i < len(faces):
        n = faces[i]
        out.append(pts[faces[i + 1:i + 1 + n]])
        i += n + 1
    return out


def _panel(ax, m, normal, origin, axes, title, xlim, ylim, quiver=False):
    sl = m.slice(normal=normal, origin=origin)
    if sl.n_cells == 0:
        ax.set_title(title + ' (vacio)')
        return
    a = np.clip(sl.cell_data['alpha.oil'], 0, 1)
    chi = sl.cell_data['chi'] if 'chi' in sl.cell_data else np.zeros_like(a)
    col = OIL_CMAP(a)
    col[chi > 0.5] = BODY
    pc = PolyCollection(_polys(sl, axes), facecolors=col, edgecolors='none', antialiased=False)
    ax.add_collection(pc)
    if quiver:
        c = sl.cell_centers().points[:, axes]
        U = sl.cell_data['U'][:, axes]
        keep = (chi < 0.5) & (a > 0.5)
        idx = np.where(keep)[0][::max(1, keep.sum() // 350)]
        if len(idx):
            ax.quiver(c[idx, 0], c[idx, 1], U[idx, 0], U[idx, 1], color='#d95f02', scale=40, width=0.003)
    ax.set_facecolor(STEEL)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=10)
    ax.tick_params(labelsize=7)


STATE_COLOR = {'wet': '#08306b', 'mixed': '#6baed6', 'dry': '#fdd49e', 'unavailable': '#7f7f7f'}


def _sensors(ax, sensors, axes, sel):
    """Draw sensors [(id, (x, y, z), state)] whose coordinate sel(pos) is True."""
    for sid, pos, state in sensors:
        if not sel(pos):
            continue
        u, v = pos[axes[0]], pos[axes[1]]
        ax.scatter([u], [v], s=46, marker='o', c=STATE_COLOR.get(state, 'w'), edgecolors='#d7301f', linewidths=1.4,
                   zorder=5)
        ax.annotate(sid, (u, v), xytext=(4, 3), textcoords='offset points', fontsize=6.5, color='#d7301f', zorder=6)


def render_frame(F: Fields, t: float, out: Path, sc=None, oil_L=None, quiver=True, dpi=110, sensors=None):
    m = F.mesh(t)
    sc = sc or load_scenario(F.case)
    fig = plt.figure(figsize=(14, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.45])
    lim = (-0.55, 0.55)
    ax0, ax1, ax2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
    _panel(ax0, m, 'z', (0, 0, STAGE_Z['High'] + 0.002), (0, 1),
           'Etapa alta (z = -1.13 m), vista desde interior', lim, lim, quiver)
    _panel(ax1, m, 'z', (0, 0, STAGE_Z['Low'] + 0.002), (0, 1),
           'Etapa baja (z = -1.49 m)', lim, lim, quiver)
    _panel(ax2, m, 'x', (0.0005, 0, 0), (2, 1),
           'Corte axial x = 0 (z horizontal, y vertical)', (-1.75, -0.98), lim, quiver)
    if sensors:
        _sensors(ax0, sensors, (0, 1), lambda p: p[2] > -1.31)
        _sensors(ax1, sensors, (0, 1), lambda p: p[2] <= -1.31)
        _sensors(ax2, sensors, (2, 1), lambda p: abs(p[0]) < 0.07)
    rpm = np.interp(t, [r[0] for r in sc.motor_rpm], [r[1] for r in sc.motor_rpm]) \
        if sc.mechanisms != 'stopped' else 0.0
    acc = np.interp(t, [r[0] for r in sc.accel], [r[1] for r in sc.accel])
    hdr = (f'{sc.name} | t = {t:.3f} s | motor {rpm:.0f} rpm | a_veh {acc:+.1f} m/s2 | '
           f'pendiente {sc.pitch_deg:.1f} deg, peralte {sc.roll_deg:.1f} deg')
    if oil_L is not None:
        hdr += f' | aceite integrado {oil_L:.2f} L'
    fig.suptitle(hdr, fontsize=11)
    sm = plt.cm.ScalarMappable(cmap=OIL_CMAP)
    sm.set_array([0, 1])
    cb = fig.colorbar(sm, ax=fig.axes, fraction=0.015, pad=0.01)
    cb.set_label('fraccion de aceite alpha', fontsize=8)
    fig.text(0.01, 0.01, 'OA2IS-WheelMotor-CFD | OpenFOAM v1912 interIBFoamOA2IS (VOF + frontera inmersa) | '
             'gris oscuro: engranajes/porta-satelites moviles; gris claro: solidos fijos; flechas: U en aceite'
             + ('\ncirculos: sensores capacitivos (relleno: azul mojado, celeste mixto, crema seco, gris no disp.); '
                'etapa alta = sensores con z > -1.31, etapa baja = resto' if sensors else ''),
             fontsize=6.5)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


def make_mp4(frames_dir: Path, out: Path, fps=8):
    cmd = ['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(Path(frames_dir) / 'frame_%04d.png'),
           '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '20', str(out)]
    subprocess.run(cmd, check=True)
