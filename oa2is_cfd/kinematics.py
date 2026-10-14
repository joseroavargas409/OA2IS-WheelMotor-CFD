"""Planetary kinematics of the wheel-motor-inspector rev. 03 model.

All speeds are expressed as ratios of the rotor (input) speed, exactly as the
original ``motions`` array: rotor/shaft = 1, high-stage carrier = intermediate
sleeve = low-stage sun = 1/5, low-stage carrier = hub = 1/(5*(1+127/23)).
Rings and spindle are fixed. Planet spin is given relative to its carrier, as
in the original (``planetRelative = -(ns/np)*(input-output)``).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from .gears import GearSpec

# Values copied from third_party/wheel-motor-inspector/dist/assembly-v2.js
LAYOUT = dict(motor=.85, high=-1.13, low=-1.49, hub=-.725, hubLength=1.95,
              spindle=-.31, spindleLength=1.36, bearings=[-.865, .115])
STAGES_SRC = [
    dict(id='High', ns=18, np=27, nr=72, input=1.0, z=LAYOUT['high'], pitch=.425),
    dict(id='Low', ns=23, np=52, nr=127, input=.2, z=LAYOUT['low'], pitch=.442),
]
MODEL_RATIO = 5 * (1 + 127 / 23)
GEAR_FACE_WIDTH = .102      # sun/planet extrusion depth in the original
RING_FACE_WIDTH = .118      # ring extrusion depth in the original


@dataclass
class Stage:
    id: str
    ns: int
    np: int
    nr: int
    input: float
    z: float
    pitch: float
    backlash: float = 0.0
    module: float = field(init=False)
    sun_r: float = field(init=False)
    planet_r: float = field(init=False)
    orbit_r: float = field(init=False)
    output: float = field(init=False)
    planet_rel: float = field(init=False)
    sun_phase: float = field(init=False)

    def __post_init__(self):
        self.module = 2 * self.pitch / self.nr
        self.sun_r = self.module * self.ns / 2
        self.planet_r = self.module * self.np / 2
        self.orbit_r = self.sun_r + self.planet_r
        self.output = self.input / (1 + self.nr / self.ns)
        self.planet_rel = -(self.ns / self.np) * (self.input - self.output)
        self.sun_phase = math.pi / self.ns if self.np % 2 else 0.0

    # --- gear specs (backlash split equally between the two members of a pair)
    def sun_spec(self):
        return GearSpec(self.ns, self.sun_r, backlash=self.backlash / 2)

    def planet_spec(self):
        return GearSpec(self.np, self.planet_r, backlash=self.backlash / 2)

    def ring_spec(self):
        return GearSpec(self.nr, self.pitch, internal=True, backlash=self.backlash / 2)

    def planet_slot_angle(self, j: int) -> float:
        return j * 2 * math.pi / 3

    def planet_phase(self, j: int) -> float:
        a = self.planet_slot_angle(j)
        return (1 - self.nr / self.np) * a - math.pi / self.np

    # --- poses as a function of rotor angle theta [rad]
    def sun_angle(self, theta):
        return self.input * theta + self.sun_phase

    def carrier_angle(self, theta):
        return self.output * theta

    def planet_center(self, j, theta):
        a = self.planet_slot_angle(j) + self.carrier_angle(theta)
        return self.orbit_r * math.cos(a), self.orbit_r * math.sin(a)

    def planet_angle(self, j, theta):
        """Absolute orientation of planet j (carrier rotation + relative spin + phase)."""
        return self.carrier_angle(theta) + self.planet_rel * theta + self.planet_phase(j)


def stages(backlash=(0.0, 0.0)):
    return [Stage(**s, backlash=b) for s, b in zip(STAGES_SRC, backlash)]
