"""Post-processing of immersed-boundary (interIBFoamOA2IS) results.

Everything here is computed from the solver's own output: the written time
directories (alpha.oil, U, p_rgh, chi), the function-object files in
postProcessing/ (oilVolume, extremes, sensors) and the solver log.

Frames: model frame of oa2is_cfd (z = wheel axis, inboard positive; y up; for a
right-hand wheel x points forward). Litres = m^3 * 1000.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .scenario import Scenario

STAGE_Z = {'High': -1.13, 'Low': -1.49}
Z_COVER_POCKET = (-1.7345, -1.62)        # outboard cover pocket (hub cover side)
GEAR_R = 0.47                            # radius of the moving-gear region
BOTTOM_SENSORS = ('S1', 'S2', 'S3', 'S7')


# ------------------------------------------------------------------ helpers
def load_scenario(case: Path) -> Scenario:
    meta = json.loads((Path(case) / 'scenario.json').read_text())
    return Scenario(**{k: v for k, v in meta.items() if k in Scenario.__dataclass_fields__})


def solver_log(case: Path) -> Path:
    case = Path(case)
    for n in ('log.solver', 'log.overInterDyMFoam', 'log.interIBFoamOA2IS'):
        if (case / n).exists():
            return case / n
    raise FileNotFoundError('no solver log in %s' % case)


def parse_log(case: Path) -> pd.DataFrame:
    """One row per time step: t, dt, Co, interface Co, p_rgh residuals,
    continuity error, wall-clock time."""
    rows, cur = [], None
    co = aco = None
    for line in solver_log(case).read_text(errors='replace').splitlines():
        if line.startswith('Courant Number'):
            co = float(line.split('max:')[1])
        elif line.startswith('Interface Courant'):
            aco = float(line.split('max:')[1])
        elif line.startswith('deltaT ='):
            dt = float(line.split('=')[1])
            cur = dict(dt=dt, Co=co, alphaCo=aco, p_rgh_res0=np.nan, p_rgh_res_final=np.nan, n_p=0)
        elif line.startswith('Time =') and cur is not None:
            cur['t'] = float(line.split()[2])
        elif 'Solving for p_rgh' in line and cur is not None:
            r0 = float(re.search(r'Initial residual = ([^,]+)', line).group(1))
            rf = float(re.search(r'Final residual = ([^,]+)', line).group(1))
            if cur['n_p'] == 0:
                cur['p_rgh_res0'] = r0
            cur['p_rgh_res_final'] = rf
            cur['n_p'] += 1
        elif line.startswith('time step continuity errors') and cur is not None:
            cur['cont_local'] = float(re.search(r'sum local = ([^,]+)', line).group(1))
        elif line.startswith('ExecutionTime') and cur is not None and 't' in cur:
            cur['wall_s'] = float(line.split()[2])
            rows.append(cur)
            cur = None
    return pd.DataFrame(rows)


def oil_volume(case: Path) -> pd.DataFrame:
    """Oil inventory from the volIntegrate(alpha.oil) function object."""
    f = sorted(Path(case).glob('postProcessing/oilVolume/*/volFieldValue*.dat'))
    frames = []
    for p in f:
        d = pd.read_csv(p, sep=r'\s+', comment='#', header=None, names=['t', 'oil_m3'])
        frames.append(d)
    d = pd.concat(frames).drop_duplicates('t', keep='last').sort_values('t').reset_index(drop=True)
    d['oil_L'] = d['oil_m3'] * 1e3
    return d


def _probe_file(path: Path, ncomp: int):
    txt = path.read_text().splitlines()
    rows = []
    for line in txt:
        if not line or line.startswith('#'):
            continue
        if ncomp == 1:
            vals = line.split()
            rows.append([float(v) for v in vals])
        else:
            t = float(line.split()[0])
            vecs = re.findall(r'\(([^()]*)\)', line)
            rows.append([t] + [float(x) for v in vecs for x in v.split()])
    return np.array(rows)


def sensors(case: Path, sc: Scenario | None = None, wet=0.5) -> pd.DataFrame:
    """Virtual OA2IS sensors from the probes function object. One row per
    (time step, sensor): alpha, oil presence, velocity, pressure, chi."""
    case = Path(case)
    sc = sc or load_scenario(case)
    dirs = sorted(case.glob('postProcessing/sensors/*'), key=lambda p: float(p.name))
    # a case prepared for the capacitive sensors records 9 points per sensor
    # (sensor_layout.json); the centre point (offset 0) is the sensor position
    lay = case / 'sensor_layout.json'
    stride = len(json.loads(lay.read_text())['offsets']) if lay.exists() else 1
    out = []
    for d in dirs:
        hdr = [l for l in (d / 'alpha.oil').read_text().splitlines() if l.startswith('# Probe') and '(' in l]
        ok = [i for i, l in enumerate(hdr) if 'Not Found' not in l]
        a = _probe_file(d / 'alpha.oil', 1)
        U = _probe_file(d / 'U', 3)
        p = _probe_file(d / 'p', 1)
        prgh = _probe_file(d / 'p_rgh', 1)
        chi = _probe_file(d / 'chi', 1) if (d / 'chi').exists() else None
        n = min(len(a), len(U), len(p), len(prgh))
        allcols = a.shape[1] == 1 + len(hdr)          # OpenFOAM writes every probe column
        for si, s in enumerate(sc.sensors):
            i = si * stride
            if i not in ok:
                continue
            k = i if allcols else ok.index(i)
            u = U[:n, 1 + 3 * k:4 + 3 * k]
            df = pd.DataFrame(dict(t=a[:n, 0], sensor=s['id'], name=s['name'],
                                   x=s['position'][0], y=s['position'][1], z=s['position'][2],
                                   alpha=a[:n, 1 + k], Ux=u[:, 0], Uy=u[:, 1], Uz=u[:, 2],
                                   p=p[:n, 1 + k], p_rgh=prgh[:n, 1 + k],
                                   chi=(chi[:n, 1 + k] if chi is not None and len(chi) >= n else np.nan)))
            out.append(df)
    df = pd.concat(out).drop_duplicates(['t', 'sensor'], keep='last').sort_values(['t', 'sensor'])
    df['U_mag'] = np.sqrt(df.Ux ** 2 + df.Uy ** 2 + df.Uz ** 2)
    df['oil_present'] = df.alpha > wet
    return df.reset_index(drop=True)


def system_state(sens: pd.DataFrame, oil: pd.DataFrame, litres0: float) -> pd.DataFrame:
    """Per-time OA2IS system state built only from the virtual sensors
    (what a real installation could observe) + the CFD inventory for audit.
    Rule (demonstrative): bottom sensors S1,S2,S3,S7 wet -> NORMAL (>=3),
    ALERTA (1-2), CRITICO (0)."""
    piv = sens.pivot_table(index='t', columns='sensor', values='alpha')
    bottom = [s for s in BOTTOM_SENSORS if s in piv.columns]
    wet = (piv[bottom] > 0.5).sum(axis=1)
    st = pd.DataFrame({'t': piv.index, 'bottom_wet': wet.values,
                       'mean_alpha_bottom': piv[bottom].mean(axis=1).values,
                       'mean_alpha_all': piv.mean(axis=1).values})
    st['state'] = np.where(st.bottom_wet >= 3, 'NORMAL', np.where(st.bottom_wet >= 1, 'ALERTA', 'CRITICO'))
    o = oil.set_index('t')['oil_L']
    st['oil_L_cfd'] = np.interp(st.t, o.index.values, o.values)
    st['oil_L_ref'] = litres0
    return st


# ------------------------------------------------------------------ fields
class Fields:
    """Lazy access to written time directories through VTK's OpenFOAM reader."""

    def __init__(self, case: Path):
        import pyvista as pv
        self.case = Path(case)
        (self.case / 'case.foam').touch()
        self.reader = pv.POpenFOAMReader(str(self.case / 'case.foam'))
        decomposed = (self.case / 'processor0').exists() and not any(
            p.name.replace('.', '', 1).replace('e-', '', 1).isdigit() and p.name != '0' for p in self.case.iterdir())
        if decomposed:
            self.reader.case_type = 'decomposed'
        self.reader.enable_all_cell_arrays()
        self.reader.disable_all_patch_arrays()
        self.reader.enable_patch_array('internalMesh')
        self.times = list(self.reader.time_values)
        self._vol = None

    def mesh(self, t):
        self.reader.set_active_time_value(t)
        m = self.reader.read()['internalMesh']
        if self._vol is None:
            self._vol = np.abs(m.compute_cell_sizes(length=False, area=False)['Volume'])
            self._cent = m.cell_centers().points
        m.cell_data['V'] = self._vol
        return m

    def metrics(self, t) -> dict:
        m = self.mesh(t)
        a = np.clip(m.cell_data['alpha.oil'], 0, 1)
        chi = m.cell_data['chi'] if 'chi' in m.cell_data else np.zeros_like(a)
        U = m.cell_data['U']
        V = self._vol
        C = self._cent
        oil = a * V
        r = np.hypot(C[:, 0], C[:, 1])
        fl = (1 - chi) * V
        gear = (r < GEAR_R) & (np.abs(C[:, 2] - STAGE_Z['High']) < 0.06) | \
               (r < GEAR_R) & (np.abs(C[:, 2] - STAGE_Z['Low']) < 0.06)
        res = dict(t=t, oil_L=oil.sum() * 1e3,
                   # oil exceeding the free (1 - chi) part of partially covered cells
                   oil_excess_in_bodies_L=float((np.maximum(0.0, a - (1 - chi)) * V).sum() * 1e3),
                   fluid_L=fl.sum() * 1e3)
        tot = oil.sum()
        for k, c in zip('xyz', C.T):
            res[f'oil_centroid_{k}'] = float((oil * c).sum() / tot)
        res['oil_above_axis_frac'] = float(oil[C[:, 1] > 0].sum() / tot)
        res['wetted_gear_frac'] = float((a * fl)[gear].sum() / max(fl[gear].sum(), 1e-30))
        for sid, z in STAGE_Z.items():
            band = np.abs(C[:, 2] - z) < 0.06
            res[f'oil_L_stage{sid}'] = float(oil[band].sum() * 1e3)
            res[f'wetted_frac_stage{sid}'] = float((a * fl)[band & (r < GEAR_R)].sum() /
                                                   max(fl[band & (r < GEAR_R)].sum(), 1e-30))
        pocket = (C[:, 2] < Z_COVER_POCKET[1])
        res['oil_L_cover_pocket'] = float(oil[pocket].sum() * 1e3)
        um = np.linalg.norm(U, axis=1)
        res['oil_mean_speed'] = float((oil * um).sum() / tot)
        res['oil_kinetic_J'] = float((0.5 * 900.0 * oil * um ** 2).sum())
        res['max_speed_fluid'] = float(um[chi < 0.5].max())
        return res


def time_label(t):
    return f'{t:.3f} s'
