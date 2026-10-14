"""CFD geometry of the wheel-motor gear case (final drive), built from the
wheel-motor-inspector rev. 03 model (third_party/wheel-motor-inspector).

World frame (same as the original): metres, rotation axis = +Z, the gear case
lies at Z < -0.89 (outboard end at the hub cover). All solids inside the oil
cavity are prisms along Z, so each one is an outline (polygon) plus a Z range.

Every dimension below is copied from the original sources; the comment gives
the call it comes from (v2 = dist/assembly-v2.js, v3 = dist/assembly-v3.js).
Geometric repairs needed for a valid CFD domain are listed in REPAIRS and in
docs/GEOMETRY.md.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

from .gears import external_profile, ring_space_profile, rotate
from .kinematics import Stage, stages as make_stages

TAU = 2 * math.pi

# Gear-case (oil cavity) axial extent
Z_COVER = -1.75 + 0.031 / 2          # v2 hub: ring(hub,.562,.14,.031,-1.025) at hub base -.725 -> inner face
Z_CLOSURE = -0.865 - 0.055 / 2       # v2 bearingSeats ring(.537,.369,.055) at -.865: outboard face (-0.8925)
R_HUB = 0.536                        # v2 hub: ring(hub,.565,.536,1.95)
R_HUB_FLANGE = 0.516                 # v2 hub: ring(hub,.600,.516,.042,-.995) -> world [-1.741,-1.699]
Z_HUB_FLANGE_END = -1.699


@dataclass
class Prism:
    """A solid: outline (+holes) extruded over [z0, z1]. `body` names the rigid
    body it belongs to (for kinematics); `patch` names the wall patch."""
    patch: str
    body: str
    outline: np.ndarray
    z0: float
    z1: float
    holes: list = field(default_factory=list)
    note: str = ''

    def polygon(self) -> Polygon:
        return Polygon(np.asarray(self.outline), [np.asarray(h) for h in self.holes]).buffer(0)


def circle(r, c=(0.0, 0.0), h=0.004, nmin=48):
    n = max(nmin, int(math.ceil(TAU * r / h)))
    t = np.linspace(0, TAU, n, endpoint=False)
    return np.c_[c[0] + r * np.cos(t), c[1] + r * np.sin(t)]


def annulus(name, body, ro, ri, z0, z1, h=0.004, note=''):
    return Prism(name, body, circle(ro, h=h), min(z0, z1), max(z0, z1), [circle(ri, h=h)] if ri > 0 else [], note)


def ring_at(name, body, ro, ri, depth, zc, h=0.004, note=''):
    """Original ringGeometry(ro, ri, depth) centred at world zc."""
    return annulus(name, body, ro, ri, zc - depth / 2, zc + depth / 2, h, note)


def hexagon(r, c=(0.0, 0.0), phase=0.0):
    """THREE.CylinderGeometry with 6 radial segments -> hexagonal prism, circumradius r."""
    t = phase + np.arange(6) * TAU / 6
    return np.c_[c[0] + r * np.sin(t), c[1] + r * np.cos(t)]


REPAIRS = [
    'Gear teeth regenerated as exact involutes with backlash (original interpenetrates: mirrored ring profile, zero backlash).',
    'Spindle ring support at low stage outboard (v2 ring(spindle,.526,.5,.018) at -1.557) trimmed to [-1.559,-1.548]: '
    'original overlaps the rotating low-carrier ring at [-1.566,-1.561].',
    'Thrust rings at -1.22 (v2 thrust: ring(.099,.062,.007,.516), ring(.110,.061,.016,.52)) omitted: the original keeps '
    'them static while they intersect the rotating high carrier plate and intermediate flange.',
    'Interior of the intermediate sleeve (r<.049) and of the shaft/spindle bores treated as closed (sealed) volumes.',
    'Radii closer than ~1-3 mm snapped for meshability: spindle tube outer .532->.536 (closes the 4 mm tube/hub gap), '
    'low ring outer .507->.510, side rings .482/.466->.483/.467 and .499->.500, plate/short tube .484/.468->.483/.467, '
    'spindle collar .281->.280, intermediate .061->.059.',
    'Oil cavity closed at the hub cover inner face (z=-1.7345) and at the outboard face of the hub-bearing seat '
    '(z=-0.8925); regions beyond belong to the thrust cap and hub-bearing cavity.',
]


@dataclass
class GearCase:
    stages: list[Stage]
    h: float = 0.004          # polygon resolution for circles [m]

    # ------------------------------------------------------------ background
    def background_solids(self, for_volume: bool = False) -> list[Prism]:
        """Fixed solids and axisymmetric rotating solids (walls in the background mesh)."""
        S = []
        h = self.h
        # hub flange step at the cover (rotating with the hub)
        S.append(annulus('hubWall', 'hub', R_HUB + 0.01, R_HUB_FLANGE, Z_COVER - 0.01, Z_HUB_FLANGE_END, h,
                         'v2 hub flange ring(.600,.516,.042,-.995)'))
        # spindle (fixed): torque tube, short tube, plate, main tube + bore, collars
        S.append(annulus('spindle', 'spindle', R_HUB + .01, .510, -1.552, -1.006, h, 'v2 ring(spindle,.532,.510,.546,-1.279); outer snapped to hub (4 mm gap closed)'))
        S.append(annulus('spindle', 'spindle', .483, .467, -1.077, -1.009, h, 'v2 ring(spindle,.484,.468,.068,-1.043) snapped 1 mm'))
        S.append(annulus('spindle', 'spindle', .483, .157, -1.015, -0.995, h, 'v2 ring(spindle,.484,.157,.020,-1.005) snapped 1 mm'))
        S.append(annulus('spindle', 'spindle', .280, 0.0, -0.990, Z_CLOSURE + 0.01, h,
                         'v2 ring(spindle,.280,.157,1.36) main tube; bore sealed (assumption)'))
        S.append(annulus('spindle', 'spindle', .302, 0.0, -0.993, -0.967, h, 'v2 ring(spindle,.302,.157,.026,-.67)'))
        S.append(annulus('spindle', 'spindle', .280, 0.0, -0.910, Z_CLOSURE + 0.01, h, 'v2 ring(spindle,.281,.157,.09,-.555) snapped 1 mm'))
        # ring-gear supports (fixed, on the spindle)
        S.append(annulus('spindle', 'spindle', .526, .483, -1.072, -1.054, h, 'v2 ring(spindle,.526,.483,.018,high+.067)'))
        S.append(annulus('spindle', 'spindle', .526, .500, -1.432, -1.414, h, 'v2 ring(spindle,.526,.5,.018,low+.067)'))
        S.append(annulus('spindle', 'spindle', .526, .500, -1.559, -1.548, h,
                         'v2 ring(spindle,.526,.5,.018,low-.067) TRIMMED (was [-1.566,-1.548])'))
        # seal and hub-bearing near the inboard closure (fixed)
        S.append(annulus('seal', 'spindle', .389, .280, -0.941, -0.919, h, 'v2 seals ring(.389,.28,.022,-.930)'))
        S.append(annulus('seal', 'spindle', .368, .280, -0.901, Z_CLOSURE + 0.01, h, 'v2 hubBearing outboard (solid envelope)'))
        # ring gears (fixed). The toothed band lives in a static overset zone
        # (ring_zones); the background keeps the smooth ring body beyond it.
        for st in self.stages:
            zc = st.z
            R = .490 if st.id == 'High' else .510   # low ring outer .507 snapped to the torque tube (.510)
            # In the background the ring body starts at the tooth-TIP radius: the
            # tooth spaces exist only in the ring zone. (A background body starting
            # at the root would leave background 'fluid' inside the teeth that
            # overset hole-cutting cannot remove, because that region touches a
            # real background wall.)
            r_in = self.ring_zone_r_out(st) if for_volume else st.pitch - st.module
            S.append(annulus('ring' + st.id, 'ring' + st.id, R, r_in, zc - .059, zc + .059, h,
                             'v2 gearGeometry(nr,pitch,.118,0,true): body from the tooth-tip radius'))
            for s in (-1, 1):
                ro, ri = (.483, .467) if st.id == 'High' else (.500, .483)   # .482/.466 and .499/.483 snapped 1 mm
                S.append(ring_at('ring' + st.id, 'ring' + st.id, ro, ri, .018, zc + s * .068, h,
                                 'v2 ring(ringP,pitch+.057,pitch+.041,.018,+-.068) snapped'))
        # low carrier outer ring: axisymmetric, joined to the hub (rotates with hub)
        S.append(annulus('carrierLowRing', 'carrierLow', .536 + .01, .381, -1.593, -1.561, h,
                         'v2 ring(carrierP,.536,.381,.032,-.087) low stage'))
        # intermediate sleeve (omega/5), solid incl. sealed bore
        for ro, z0, z1, src in [(.059, -1.558, -1.218, 'ring(.059,.049,.34,0)'), (.081, -1.228, -1.206, 'ring(.081,.049,.022,.171)'),
                                (.064, -1.5725, -1.5575, 'ring(.064,.049,.015,-.177)'), (.059, -1.742, -1.548, 'ring(.061,.049,.194,-.257) snapped .061->.059')]:
            S.append(annulus('intermediate', 'intermediate', ro, 0, max(z0, Z_COVER - .01), z1, h, 'v2 intermediate ' + src))
        # rotor shaft (omega): ends at -1.205
        S.append(annulus('shaft', 'shaft', .046, 0, -1.205, -0.985, h, 'v2 shaft cylinder(.046,.046,2.545)'))
        return S

    # ------------------------------------------------------------ ring zones
    @staticmethod
    def ring_zone_r_out(st):
        """Outer radius of the static toothed-band zone: 6 mm beyond the tooth root."""
        return st.pitch + 1.25 * st.module + .006

    def ring_zones(self) -> list[dict]:
        """Static overset zones holding the regenerated internal teeth."""
        Z = []
        for st in self.stages:
            zc = st.z
            space = ring_space_profile(st.ring_spec())
            Z.append(dict(name='ring' + st.id, body='ring' + st.id, center=(0.0, 0.0),
                          r_out=self.ring_zone_r_out(st), r_in=st.pitch - st.module - .010,
                          z0=zc - .065, z1=zc + .065, static=True,
                          solids=[Prism('ring' + st.id, 'ring' + st.id, circle(self.ring_zone_r_out(st) + .01, h=self.h),
                                        zc - .059, zc + .059, [space], 'v2 internal gear regenerated involute')]))
        return Z

    # ------------------------------------------------------------ components
    def components(self) -> list[dict]:
        """Moving component regions (overset meshes or immersed bodies).
        Each: name, body, centre (rest pose), outer radius, inner radius,
        z0, z1, solids (Prisms in world coordinates, rest pose)."""
        C = []
        h = self.h
        for st in self.stages:
            zc = st.z
            # sun
            sun = rotate(external_profile(st.sun_spec()), st.sun_angle(0))
            core = .046 if st.id == 'High' else .059   # shaft / intermediate inside the sun bore
            C.append(dict(name='sun' + st.id, body='sun' + st.id, center=(0.0, 0.0),
                          r_out=st.sun_r + st.module + .015, r_in=0.0, z0=zc - .065, z1=zc + .065,
                          solids=[Prism('sun' + st.id, 'sun' + st.id, sun, zc - .051, zc + .051, note='v2 gearGeometry(ns,sunR,.102,bore) regenerated'),
                                  annulus('sun' + st.id, 'sun' + st.id, core, 0, zc - .065, zc + .065, h,
                                          'shaft (high) / intermediate (low) inside the sun bore, same speed')]))
            # planets
            for j in range(3):
                c = st.planet_center(j, 0)
                pl = rotate(external_profile(st.planet_spec()), st.planet_angle(j, 0)) + np.array(c)
                hub = circle(.044, c, h=h / 2)
                C.append(dict(name=f'planet{st.id}{j}', body=f'planet{st.id}{j}', center=c,
                              r_out=st.planet_r + st.module + .015, r_in=0.0, z0=zc - .065, z1=zc + .065,
                              solids=[Prism(f'planet{st.id}', f'planet{st.id}{j}', pl, zc - .051, zc + .051, note='v2 gearGeometry(np,planetR,.102,.036) regenerated'),
                                      Prism(f'planet{st.id}', f'planet{st.id}{j}', hub, zc - .065, zc + .065,
                                            note='v2 thrust washers ring(.044,.036,.014,+-.06) + pin: solid core r<=.044')]))
            # carrier plates
            win = [circle(.080, (.198 * math.cos(math.pi / 3 + i * TAU / 3), .198 * math.sin(math.pi / 3 + i * TAU / 3)), h=h)
                   for i in range(3)]
            pins = lambda z0, z1: [Prism('carrier' + st.id, 'carrier' + st.id, circle(.023, st.planet_center(j, 0), h=h / 2), z0, z1,
                                         note='v2 pin cylinder(.023,.023,.221)') for j in range(3)]
            heads = [Prism('carrier' + st.id, 'carrier' + st.id, hexagon(.034, st.planet_center(j, 0)), zc - .125, zc - .109,
                           note='v2 pin head cylinder(.034,.034,.016,6)') for j in range(3)]
            plate = lambda z0, z1: Prism('carrier' + st.id, 'carrier' + st.id, circle(.384, h=h), z0, z1, [circle(.082, h=h)] + win,
                                         'v2 ring(carrierP,.384,.082,.025,+-.087,holes)')
            if st.id == 'High':
                zin = (zc + .0675, zc + .106)
                C.append(dict(name='carrierHighIn', body='carrierHigh', center=(0, 0), r_out=.40, r_in=.064,
                              z0=zin[0], z1=zin[1], solids=[plate(zc + .0745, zc + .0995)] + pins(*zin)))
                zout = (zc - .137, zc - .0675)
                C.append(dict(name='carrierHighOut', body='carrierHigh', center=(0, 0), r_out=.40, r_in=0.0,
                              z0=zout[0], z1=zout[1],
                              solids=[plate(zc - .0995, zc - .0745)] + pins(zc - .1105, zout[1]) + heads + [
                                  annulus('carrierHigh', 'carrierHigh', .086, 0, zc - .1045, zc - .0695, h, 'v2 ring(carrierP,.086,.049,.035,-.087) + struts'),
                                  annulus('carrierHigh', 'carrierHigh', .059, 0, zout[0], zc - .088, h, 'intermediate sleeve (same speed); '
                                          'flange r.081 and bore r.049 lie inside the hub ring / sleeve')]))
            else:
                zin = (zc + .0675, zc + .137)
                C.append(dict(name='carrierLowIn', body='carrierLow', center=(0, 0), r_out=.40, r_in=.070,
                              z0=zin[0], z1=zin[1], solids=[plate(zc + .0745, zc + .0995)] + pins(zin[0], zc + .1105)))
                zout = (zc - .137, zc - .0675)
                C.append(dict(name='carrierLowOut', body='carrierLow', center=(0, 0), r_out=.40, r_in=.072,
                              z0=zout[0], z1=zout[1],
                              solids=[plate(zc - .0995, zc - .0745)] + pins(zc - .1105, zout[1]) + heads + [
                                  annulus('carrierLow', 'carrierLow', .40 + .01, .381, zc - .103, zc - .071, h, 'v2 carrier ring .381-.536 (inner part)')]))
        return C

    # ------------------------------------------------------------ volumes
    def outer_radius(self, z):
        return R_HUB_FLANGE if z < Z_HUB_FLANGE_END else R_HUB

    def z_breaks(self):
        zs = {Z_COVER, Z_CLOSURE, Z_HUB_FLANGE_END}
        for p in self.all_solids():
            zs.update([p.z0, p.z1])
        return sorted(z for z in zs if Z_COVER <= z <= Z_CLOSURE)

    def all_solids(self):
        out = list(self.background_solids(for_volume=True))
        for c in self.components() + self.ring_zones():
            out += c['solids']
        return out

    def fluid_slices(self, level=None, normal=(0.0, 1.0, 0.0)):
        """Fluid cross-section area per z interval. With `level`, only the part
        below the plane n.x = level (n = unit 'up' = -g/|g|). Returns (z0, z1, area)."""
        zb = self.z_breaks()
        polys = [(s, s.polygon()) for s in self.all_solids()]
        n = np.asarray(normal, float)
        n = n / np.linalg.norm(n)
        out = []
        for z0, z1 in zip(zb[:-1], zb[1:]):
            zm = 0.5 * (z0 + z1)
            disc = Point(0, 0).buffer(self.outer_radius(zm), 512)
            act = [p for s, p in polys if s.z0 < zm < s.z1]
            fluid = disc.difference(unary_union(act)) if act else disc
            if level is not None:
                nxy = n[:2]
                h = level - n[2] * zm
                if np.linalg.norm(nxy) < 1e-9:      # gravity along the axis
                    fluid = fluid if h > 0 else Polygon()
                else:
                    m = nxy / np.linalg.norm(nxy)
                    hh = h / np.linalg.norm(nxy)
                    t = np.array([-m[1], m[0]])
                    p0 = m * hh
                    big = 10.0
                    half = Polygon([p0 + t * big, p0 - t * big, p0 - t * big - m * big, p0 + t * big - m * big])
                    fluid = fluid.intersection(half)
            out.append((z0, z1, fluid.area))
        return out

    def fluid_volume(self, **kw):
        return sum((z1 - z0) * a for z0, z1, a in self.fluid_slices(**kw))

    def level_for_volume(self, litres, normal=(0.0, 1.0, 0.0), tol_l=0.01):
        """Plane offset `level` (n.x = level) holding `litres` below it."""
        target = litres / 1000.0
        lo, hi = -1.0, 1.0
        lo, hi = -2.5, 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            v = self.fluid_volume(level=mid, normal=normal)
            if abs(v - target) * 1000 < tol_l:
                return mid
            lo, hi = (mid, hi) if v < target else (lo, mid)
        return 0.5 * (lo + hi)


def default_case(backlash=(0.004, 0.003)) -> GearCase:
    return GearCase(make_stages(backlash))
