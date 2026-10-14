"""Involute spur-gear profiles for the CFD geometry.

Why this module exists
----------------------
wheel-motor-inspector (rev. 03) builds gears with ``gearGeometry`` in
``dist/model.js``. Its internal ring is the *mirror* of an external profile
(``rad = 2*pitch - r``) and its teeth carry zero backlash. Sliced at the
mid-plane, the rest pose interpenetrates (high stage: planet/sun 62 mm^2,
planet/ring 32 mm^2 at the optimal phase; see docs/GEOMETRY.md). Interpenetrating
solids cannot bound a CFD fluid domain, so these profiles are regenerated here
with the SAME tooth counts, pitch radii, centre positions, tooth phases and
face widths as the original, using:

* exact involute flanks from the base circle (pressure angle 20 deg),
* addendum 1.0 m, dedendum 1.25 m (same circles as the original),
* internal ring teeth as the complement of an involute *space*,
* an explicit, configurable circumferential backlash (numerical clearance).

Tooth i of every gear is centred at polar angle i*2*pi/N in the gear frame,
the same convention as the original generator, so its motion phases apply
unchanged.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

PRESSURE_ANGLE = math.radians(20.0)


def inv(alpha: float) -> float:
    return math.tan(alpha) - alpha


@dataclass(frozen=True)
class GearSpec:
    teeth: int
    pitch_radius: float          # m
    internal: bool = False
    backlash: float = 0.0        # circumferential backlash removed from THIS gear's tooth, at pitch circle [m]
    addendum: float = 1.0        # in modules
    dedendum: float = 1.25       # in modules

    @property
    def module(self) -> float:
        return 2.0 * self.pitch_radius / self.teeth

    @property
    def base_radius(self) -> float:
        return self.pitch_radius * math.cos(PRESSURE_ANGLE)


def _half_angle_at(r: float, rb: float, half_pitch_angle: float) -> float:
    """Half angular thickness (from the tooth centre line) of an involute
    tooth at radius r, given its half angular thickness at the pitch circle.
    Below the base circle the flank is continued radially."""
    rp = rb / math.cos(PRESSURE_ANGLE)
    if r <= rb:
        return half_pitch_angle + inv(PRESSURE_ANGLE)
    a = math.acos(rb / r)
    return half_pitch_angle + inv(PRESSURE_ANGLE) - inv(a)


def _tooth_polygon(n_teeth: int, rp: float, rb: float, r_root: float, r_tip: float,
                   half_pitch_angle: float, n_flank: int = 10, n_land: int = 3) -> np.ndarray:
    """Closed outline (counter-clockwise) of an external-gear-shaped star whose
    teeth have the given half angular thickness at the pitch circle."""
    pts = []
    pitch = 2 * math.pi / n_teeth
    r0 = max(r_root, min(rb, r_tip))
    radii = np.linspace(r0, r_tip, n_flank)
    for i in range(n_teeth):
        c = i * pitch
        # root land from previous gap centre to start of the left flank
        h_root = _half_angle_at(r_root, rb, half_pitch_angle)
        for t in np.linspace(c - pitch / 2, c - h_root, n_land, endpoint=False):
            pts.append((r_root * math.cos(t), r_root * math.sin(t)))
        if r0 > r_root:  # radial segment below the base circle
            t = c - h_root
            pts.append((r_root * math.cos(t), r_root * math.sin(t)))
        for r in radii:  # left (leading) flank going up
            t = c - _half_angle_at(r, rb, half_pitch_angle)
            pts.append((r * math.cos(t), r * math.sin(t)))
        h_tip = _half_angle_at(r_tip, rb, half_pitch_angle)
        for t in np.linspace(c - h_tip, c + h_tip, n_land + 2)[1:-1]:
            pts.append((r_tip * math.cos(t), r_tip * math.sin(t)))
        for r in radii[::-1]:  # right flank going down
            t = c + _half_angle_at(r, rb, half_pitch_angle)
            pts.append((r * math.cos(t), r * math.sin(t)))
        # root land (its first point closes the radial segment below the base circle)
        for t in np.linspace(c + h_root, c + pitch / 2, n_land, endpoint=False):
            pts.append((r_root * math.cos(t), r_root * math.sin(t)))
    return np.asarray(pts)


def external_profile(g: GearSpec, **kw) -> np.ndarray:
    """Outline of an external gear (solid inside)."""
    m = g.module
    s = math.pi * m / 2 - g.backlash          # tooth thickness at pitch circle
    half = s / (2 * g.pitch_radius)
    return _tooth_polygon(g.teeth, g.pitch_radius, g.base_radius,
                          g.pitch_radius - g.dedendum * m, g.pitch_radius + g.addendum * m, half, **kw)


def internal_space_profile(g: GearSpec, **kw) -> np.ndarray:
    """Outline of the fluid-side *space* of an internal ring gear: a star whose
    'teeth' are the ring's tooth spaces. Ring solid = annulus(R_out) - this.
    The ring's own tooth (between spaces) is centred at angle (i+1/2)*pitch in
    this outline, so call with the result rotated by -pitch/2 to centre ring
    tooth i at i*pitch (done in ``ring_space_profile``)."""
    m = g.module
    e = math.pi * m / 2 + g.backlash          # space width = pitch - ring tooth thickness
    half = e / (2 * g.pitch_radius)
    # space reaches the ring root (rp + dedendum) and opens at the ring tip (rp - addendum)
    return _tooth_polygon(g.teeth, g.pitch_radius, g.base_radius,
                          g.pitch_radius - g.addendum * m, g.pitch_radius + g.dedendum * m, half, **kw)


def ring_space_profile(g: GearSpec, **kw) -> np.ndarray:
    """Internal-gear space outline with ring tooth i centred at i*2*pi/N.

    The original generator mirrors (r -> 2*rp - r) an external outline whose
    tooth i is centred at i*2*pi/N; the mirrored tooth reaches the smallest
    radius (rp - m), i.e. it becomes ring material. Ring teeth are therefore
    centred at i*2*pi/N and spaces at (i+1/2)*2*pi/N. We keep that convention."""
    return rotate(internal_space_profile(g, **kw), math.pi / g.teeth)


def rotate(pts: np.ndarray, angle: float, center=(0.0, 0.0)) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    p = np.asarray(pts) - center
    return np.c_[c * p[:, 0] - s * p[:, 1], s * p[:, 0] + c * p[:, 1]] + center
