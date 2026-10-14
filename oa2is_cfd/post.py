"""Post-processing of overset results: de-duplicated oil inventory, virtual
OA2IS sensors and field extraction.

Overset meshes overlap, so a plain volume integral counts the overlap twice.
Each cell gets an *active* flag:
  active = (cellType != HOLE) and (cell centre not inside the footprint of a
           higher-priority zone at that time)
Priority (high -> low): planets, suns, carriers, ring zones, background.
Footprints are the zone cylinders/annuli moved with the analytic kinematics.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pyvista as pv

from .kinematics import stages as make_stages
from .scenario import Scenario

PRIORITY = ('planet', 'sun', 'carrier', 'ring', 'background')


def zone_rank(name):
    for i, p in enumerate(PRIORITY):
        if name.startswith(p):
            return i
    return len(PRIORITY)


class Case:
    def __init__(self, case):
        self.case = Path(case)
        meta = json.loads((self.case / 'scenario.json').read_text())
        self.meta = meta
        self.sc = Scenario(**{k: v for k, v in meta.items() if k in Scenario.__dataclass_fields__})
        self.zones = meta['zones']                      # index = zoneID
        self.stats = json.loads((self.case / 'mesh_stats.json').read_text())
        self.zmeta = self.stats['zone_meta']
        self.st = {s.id: s for s in make_stages(tuple(meta['backlash']))}
        (self.case / 'case.foam').touch()
        decomposed = (self.case / 'processor0').exists() and not any(
            p.name.replace('.', '', 1).isdigit() and p.name != '0' for p in self.case.iterdir())
        self.reader = pv.POpenFOAMReader(str(self.case / 'case.foam'))
        if decomposed:
            self.reader.case_type = 'decomposed'
        self.reader.enable_all_cell_arrays()
        self.times = [t for t in self.reader.time_values]

    # ---------------------------------------------------------- kinematics
    def rotor_angle(self, t):
        tab = self.sc.rotor_omega_table()
        ts, ws = np.array([r[0] for r in tab]), np.array([r[1] for r in tab])
        grid = np.unique(np.r_[ts[ts < t], t, 0.0])
        return float(np.trapezoid(np.interp(grid, ts, ws), grid))

    def rotor_omega(self, t):
        tab = self.sc.rotor_omega_table()
        return float(np.interp(t, [r[0] for r in tab], [r[1] for r in tab]))

    def zone_center(self, name, t):
        c = np.array(self.zmeta[name]['center'], float)
        if name.startswith('planet') and not self.zmeta[name]['static']:
            S = self.st['High' if 'High' in name else 'Low']
            a = S.output * self.rotor_angle(t)
            ca, sa = math.cos(a), math.sin(a)
            c = np.array([ca * c[0] - sa * c[1], sa * c[0] + ca * c[1]])
        return c

    def footprint(self, name, t, pts):
        zm = self.zmeta[name]
        c = self.zone_center(name, t)
        r = np.hypot(pts[:, 0] - c[0], pts[:, 1] - c[1])
        return (r < zm['r_out']) & (r > zm['r_in']) & (pts[:, 2] > zm['z0']) & (pts[:, 2] < zm['z1'])

    # ---------------------------------------------------------- reading
    def read(self, t):
        self.reader.set_active_time_value(t)
        mb = self.reader.read()
        m = mb['internalMesh']
        return m

    def active_mask(self, m, t):
        cc = m.cell_centers().points
        zid = np.asarray(m.cell_data['zoneID']).astype(int)
        ct = np.asarray(m.cell_data['cellTypes']).astype(int) if 'cellTypes' in m.cell_data else np.zeros(len(zid), int)
        active = ct != 2
        ranks = {z: zone_rank(z) for z in self.zones}
        for zi, zname in enumerate(self.zones):
            sel = zid == zi
            if not sel.any():
                continue
            higher = [z for z in self.zones if ranks[z] < ranks[zname] and z in self.zmeta]
            cov = np.zeros(sel.sum(), bool)
            p = cc[sel]
            for h in higher:
                cov |= self.footprint(h, t, p)
            a = active[sel]
            a &= ~cov
            active[sel] = a
        return active, cc, zid, ct

    def inventory(self, t, m=None):
        m = m if m is not None else self.read(t)
        active, cc, zid, ct = self.active_mask(m, t)
        V = m.compute_cell_sizes(length=False, area=False, volume=True).cell_data['Volume']
        a = np.asarray(m.cell_data['alpha.oil'])
        return dict(time=t, oil_L=float(np.sum(a * V * active) * 1e3), fluid_L=float(np.sum(V * active) * 1e3),
                    oil_all_cells_L=float(np.sum(a * V) * 1e3))

    def sensors(self, t, sensors, radius=0.010, m=None):
        """Virtual probes: point values from the highest-priority active cell
        containing the point and a volume-averaged oil fraction within `radius`."""
        m = m if m is not None else self.read(t)
        active, cc, zid, ct = self.active_mask(m, t)
        V = m.compute_cell_sizes(length=False, area=False, volume=True).cell_data['Volume']
        a = np.asarray(m.cell_data['alpha.oil'])
        U = np.asarray(m.cell_data['U'])
        p = np.asarray(m.cell_data['p']) if 'p' in m.cell_data else np.full(len(a), np.nan)
        prgh = np.asarray(m.cell_data['p_rgh'])
        rows = []
        g = self.sc.g_eff(t)
        for s in sensors:
            x = np.asarray(s['position'], float)
            d = np.linalg.norm(cc - x, axis=1)
            near = np.where((d < radius) & active)[0]
            if len(near):
                w = V[near]
                avol = float(np.sum(a[near] * w) / np.sum(w))
                i = near[np.argmin(d[near])]
                ok = True
            else:
                cand = np.where(active)[0]
                i = cand[np.argmin(d[cand])]
                avol = float(a[i])
                ok = d[i] < 3 * radius
            rows.append(dict(time=t, sensor=s['id'], x=x[0], y=x[1], z=x[2], valid=bool(ok),
                             alpha_point=float(a[i]), alpha_vol=avol, oil_present=bool(avol > 0.5),
                             Ux=float(U[i, 0]), Uy=float(U[i, 1]), Uz=float(U[i, 2]), Umag=float(np.linalg.norm(U[i])),
                             p=float(p[i]), p_rgh=float(prgh[i]), zone=self.zones[zid[i]],
                             rotor_omega=self.rotor_omega(t), gx=float(g[0]), gy=float(g[1]), gz=float(g[2]),
                             pitch_deg=self.sc.pitch_deg, roll_deg=self.sc.roll_deg))
        return rows
