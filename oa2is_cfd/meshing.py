"""Build the multi-zone (overset) mesh of the gear case from oa2is_cfd.model.

Zones (cellZones, one per mesh region):
  componentes móviles : sunHigh, planetHigh0..2, carrierHighIn, carrierHighOut,
                        sunLow,  planetLow0..2,  carrierLowIn,  carrierLowOut
  fondo fijo          : background (hub walls, spindle, ring gears, sleeve, shaft)

Each zone is meshed with oa2is_cfd.layered (2D gmsh mesh with every solid
outline embedded, extruded through Z levels that contain every solid end).
Zones are merged with mergeMeshes; components come first so that the first
patch of the merged mesh is an overset patch (OpenFOAM requirement).
"""
from __future__ import annotations

import json
import math
import shutil
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from . import foam
from .layered import Solid, extrude, merge, mesh2d, write_polymesh
from .model import GearCase, Z_CLOSURE, Z_COVER, circle


@dataclass
class Resolution:
    name: str = 'coarse'
    # background
    bg_lc: float = 0.016            # bulk 2D size
    bg_ring_lc: float = 0.006       # on ring-gear teeth
    bg_wall_lc: float = 0.012       # on other embedded outlines
    bg_grading: float = 0.03
    bg_dz_stage: float = 0.012      # axial size around the two stages
    bg_dz_far: float = 0.025
    stage_halfwidth: float = 0.15   # axial half-width of the 'stage' band
    # components
    gear_lc: dict = None            # per stage, on gear teeth
    comp_lc_max: float = 0.012
    comp_grading: float = 0.02
    gear_layers: int = 6            # across the 0.102 m face width
    comp_dz: float = 0.007          # axial size outside the gear face
    carrier_lc: float = 0.012
    carrier_dz: float = 0.008
    z_snap: float = 0.003           # merge axial break levels closer than this
    ring_lc: float = 0.006          # on ring teeth (static ring zones)

    def __post_init__(self):
        if self.gear_lc is None:
            self.gear_lc = {'High': 0.005, 'Low': 0.0045}


PRESETS = {
    # low-cost configuration (what runs on a 2-core workstation in hours)
    'coarse': Resolution(name='coarse', bg_lc=0.024, bg_ring_lc=0.008, bg_wall_lc=0.018, bg_dz_stage=0.020,
                         bg_dz_far=0.035, gear_lc={'High': 0.006, 'Low': 0.005}, comp_lc_max=0.014,
                         gear_layers=5, comp_dz=0.007, carrier_lc=0.016, carrier_dz=0.010, ring_lc=0.007),
    # reference configuration (~2-3x finer; for a workstation/cluster)
    'medium': Resolution(name='medium'),
    'fine': Resolution(name='fine', bg_lc=0.008, bg_ring_lc=0.003, bg_wall_lc=0.005, bg_dz_stage=0.005,
                       bg_dz_far=0.012, gear_lc={'High': 0.0025, 'Low': 0.002}, comp_lc_max=0.006,
                       gear_layers=16, comp_dz=0.004, carrier_lc=0.004, carrier_dz=0.003, ring_lc=0.003),
}


def cluster(breaks, tol, keep=()):
    """Merge break levels closer than `tol` (keeps the ends and any `keep` level).
    Solid ends that fall between levels are resolved to the nearest level by the
    cell-midpoint classification of oa2is_cfd.layered.extrude."""
    b = sorted(set(np.round(np.asarray(breaks, float), 7)))
    out = [b[0]]
    for z in b[1:]:
        if z - out[-1] < tol:
            if z == b[-1] or z in keep:
                out[-1] = z
            continue
        out.append(z)
    return out


def subdivide(breaks, dz_of, tol=0.0):
    """Insert levels so that no interval exceeds dz_of(z_mid)."""
    b = cluster(breaks, tol) if tol > 0 else sorted(set(np.round(np.asarray(breaks, float), 7)))
    out = [b[0]]
    for z0, z1 in zip(b[:-1], b[1:]):
        if z1 - z0 < 1e-6:
            continue
        n = max(1, int(math.ceil((z1 - z0) / dz_of(0.5 * (z0 + z1)) - 1e-9)))
        out += list(np.linspace(z0, z1, n + 1)[1:])
    return out


def _solids(prisms, zlo, zhi):
    out = []
    for p in prisms:
        if p.z1 <= zlo or p.z0 >= zhi:
            continue
        out.append(Solid(p.patch, p.outline, max(p.z0, zlo), min(p.z1, zhi), p.holes))
    return out


def _outlines(prisms):
    loops = []
    for p in prisms:
        loops.append(np.asarray(p.outline))
        loops += [np.asarray(h) for h in p.holes]
    return loops


def build_component(c: dict, res: Resolution, stage_id: str):
    zlo, zhi = c['z0'], c['z1']
    solids = _solids(c['solids'], zlo, zhi)
    is_gear = c['name'].startswith(('sun', 'planet', 'ring'))
    lc = (res.ring_lc if c['name'].startswith('ring') else res.gear_lc[stage_id]) if is_gear else res.carrier_lc
    emb = []
    cx, cy = c['center']
    inside = lambda loop: np.max(np.hypot(loop[:, 0] - cx, loop[:, 1] - cy)) < c['r_out'] - 1e-4
    for p in c['solids']:
        is_teeth = p.patch.startswith(('sun', 'planet', 'ring')) and (len(p.outline) > 200 or
                                                                     any(len(hh) > 200 for hh in p.holes))
        for loop in [p.outline] + list(p.holes):
            loop = np.asarray(loop)
            if inside(loop):
                emb.append((loop, lc if (is_teeth or not is_gear) else 1.5 * lc, res.comp_grading))
    outer = circle(c['r_out'], c['center'], h=res.comp_lc_max)
    holes = [(circle(c['r_in'], c['center'], h=lc), lc)] if c['r_in'] > 0 else []
    xy, el = mesh2d(outer, emb, res.comp_lc_max, holes=holes, lc_max=res.comp_lc_max)
    brk = [zlo, zhi] + [z for s in solids for z in (s.z0, s.z1)]
    if is_gear:
        zc = 0.5 * (zlo + zhi)
        hw = .059 if c['name'].startswith('ring') else .051
        face = (zc - hw, zc + hw)
        dz = lambda z: (face[1] - face[0]) / res.gear_layers if face[0] < z < face[1] else res.comp_dz
    else:
        dz = lambda z: res.carrier_dz
    zl = subdivide(brk, dz, tol=res.z_snap)
    patch = 'overset' + c['name'][0].upper() + c['name'][1:]
    m = extrude(xy, el, zl, solids, patch, (patch, patch))
    return m, patch


def build_background(gc: GearCase, res: Resolution, comps):
    prisms = gc.background_solids()
    emb = []
    for p in prisms:
        teeth = p.patch.startswith('ring') and len(p.holes) and len(p.holes[0]) > 200
        for loop in [p.outline] + list(p.holes):
            loop = np.asarray(loop)
            if np.max(np.hypot(loop[:, 0], loop[:, 1])) > 0.5361:   # outlines outside the domain
                continue
            emb.append((loop, res.bg_ring_lc if teeth else res.bg_wall_lc, res.bg_grading))
    outer = circle(0.536, h=res.bg_wall_lc)
    xy, el = mesh2d(outer, emb, res.bg_wall_lc, lc_max=res.bg_lc)
    zst = [s.z for s in gc.stages]
    dz = lambda z: res.bg_dz_stage if min(abs(z - a) for a in zst) < res.stage_halfwidth else res.bg_dz_far
    brk = [Z_COVER, Z_CLOSURE] + [z for p in prisms for z in (p.z0, p.z1) if Z_COVER < z < Z_CLOSURE]
    brk += [z for c in comps for z in (c['z0'], c['z1'])]
    zl = subdivide(brk, dz, tol=res.z_snap)
    solids = _solids(prisms, Z_COVER, Z_CLOSURE)
    m = extrude(xy, el, zl, solids, 'hubWall', ('hubCover', 'inboardClosure'))
    return m


def build_all(case: Path, gc: GearCase, res: Resolution, log=print, only=None):
    case = Path(case)
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    # discretise circular outlines at the local mesh size (embedded points are kept by gmsh)
    gcc = GearCase(gc.stages, h=res.carrier_lc)
    comps = gcc.components() + gcc.ring_zones()
    all_comps = comps
    if only is not None:
        comps = [c for c in comps if c['name'] in only]
    gc_bg = GearCase(gc.stages, h=res.bg_wall_lc)
    meshes, stats = [], {}
    for c in comps:
        t0 = time.time()
        stage = 'High' if 'High' in c['name'] else 'Low'
        m, patch = build_component(c, res, stage)
        meshes.append((c['name'], m))
        stats[c['name']] = dict(cells=m['ncells'], seconds=round(time.time() - t0, 1), body=c['body'])
        log(f"  {c['name']:16s} {m['ncells']:8d} cells ({time.time() - t0:.1f}s)")
    t0 = time.time()
    mb = build_background(gc_bg, res, all_comps)
    meshes.append(('background', mb))
    stats['background'] = dict(cells=mb['ncells'], seconds=round(time.time() - t0, 1), body='fixed')
    log(f"  {'background':16s} {mb['ncells']:8d} cells ({time.time() - t0:.1f}s)")
    types = {}
    for _, m in meshes:
        for p, *_ in m['patches']:
            types[p] = 'overset' if p.startswith('overset') else 'wall'
    merged, zranges = merge(meshes)
    write_polymesh(case, merged, None, types, zones_override=zranges)
    foam.write(case, 'system/controlDict', 'application none; startFrom startTime; startTime 0; stopAt endTime; '
                                           'endTime 0; deltaT 1; writeControl timeStep; writeInterval 1;')
    stats['total_cells'] = sum(v['cells'] for k, v in stats.items() if isinstance(v, dict))
    stats['resolution'] = asdict(res)
    stats['zones'] = [n for n, _ in meshes]
    stats['backlash'] = [st.backlash for st in gc.stages]
    # footprint of each overset zone (rest pose) for post-processing de-duplication
    stats['zone_meta'] = {c['name']: dict(body=c['body'], center=list(map(float, c['center'])), r_in=float(c['r_in']),
                                          r_out=float(c['r_out']), z0=float(c['z0']), z1=float(c['z1']),
                                          static=bool(c.get('static', False))) for c in comps}
    (case / 'mesh_stats.json').write_text(json.dumps(stats, indent=1))
    return stats
