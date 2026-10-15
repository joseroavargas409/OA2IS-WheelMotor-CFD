"""Operating scenarios (YAML) -> physical inputs for the CFD case.

Frames
------
Body frame of the truck: x forward, y up, z to the truck's left (right-handed).
Model frame = frame of the wheel-motor model (axis +Z points inboard).
* right-side wheel: model frame == body frame (inboard = +left).
* left-side wheel : model frame = body frame rotated 180 deg about y
                    (x_model = -x_body, z_model = -z_body).
The CFD is solved in the model frame attached to the truck (non-rotating; the
mechanisms rotate inside it). Road grade (pitch) and cross slope (roll) rotate
gravity; vehicle accelerations enter as an apparent body force -a. Both are
combined ONCE into the effective gravity g_eff(t) = R^T g - a(t), applied with
OpenFOAM's tabulatedAccelerationSource (g = g0 - a_tab), so inertial effects
are not counted twice. Vehicle rotation (yaw rate) is neglected.

Wheel rotation for forward travel (v > 0): omega_z,model = -v/R (right wheel),
+v/R (left wheel); the motor turns 32.6087 times faster in the same sense.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

G = 9.81


def walther_nu(T_c, ref):
    """ASTM D341 (Walther) kinematic viscosity [m2/s] at T_c from two reference
    points ref = [[T1_c, nu1_cSt], [T2_c, nu2_cSt]]."""
    (t1, n1), (t2, n2) = ref
    f = lambda n: math.log10(math.log10(n + 0.7))
    T1, T2, T = t1 + 273.15, t2 + 273.15, T_c + 273.15
    b = (f(n1) - f(n2)) / (math.log10(T2) - math.log10(T1))
    a = f(n1) + b * math.log10(T1)
    nu_cst = 10 ** (10 ** (a - b * math.log10(T))) - 0.7
    return nu_cst * 1e-6


@dataclass
class Scenario:
    name: str
    description: str
    duration: float
    write_interval: float
    motor_rpm: list            # [[t, rpm_magnitude], ...] (forward travel positive)
    side: str = 'right'
    oil_litres: float = 95.0
    oil_rho: float = 900.0
    oil_temperature_c: float = 60.0
    oil_nu: float | None = 2.64238e-4
    oil_viscosity_ref: list | None = None     # Walther reference points (cSt), overrides oil_nu
    air_rho: float = 1.06
    air_nu: float = 1.9e-5
    sigma: float = 0.030
    pitch_deg: float = 0.0     # road grade, nose up positive
    roll_deg: float = 0.0      # cross slope, left side up positive
    accel: list = field(default_factory=lambda: [[0.0, 0.0, 0.0]])   # [[t, a_forward, a_left], ...] m/s2 (body)
    numerics: dict = field(default_factory=dict)
    sensors: list = field(default_factory=list)
    mechanisms: str = 'moving'                 # 'moving' or 'stopped' (comparison case)
    mesh: str = 'coarse'
    initial_level: str = 'flat'                # flat free surface normal to g_eff(0)
    oil_level_y_flat: float | None = None      # if set: inventory = fluid volume below this height (level ground, rest)
    rotating_sensors: list = field(default_factory=list)   # sensors fixed to a rotating part (e.g. hub-cover plugs)
    capacitive_stencil: bool = False           # record the 9-point electrode volume of every sensor (oil_sensor.py)

    @staticmethod
    def load(path):
        d = yaml.safe_load(Path(path).read_text())
        return Scenario(**d)

    # ---------------------------------------------------------- physics
    def nu_oil(self):
        if self.oil_viscosity_ref:
            return walther_nu(self.oil_temperature_c, self.oil_viscosity_ref)
        return self.oil_nu

    def body_to_model(self, v):
        v = np.asarray(v, float)
        return v if self.side == 'right' else np.array([-v[0], v[1], -v[2]])

    def gravity_body(self):
        th, ph = math.radians(self.pitch_deg), math.radians(self.roll_deg)
        # world 'down' expressed in the pitched + rolled body frame
        return G * np.array([-math.sin(th), -math.cos(th) * math.cos(ph), -math.cos(th) * math.sin(ph)])

    def g0_model(self):
        return self.body_to_model(self.gravity_body())

    def accel_model(self):
        """[(t, a_model[3])] vehicle acceleration in the model frame."""
        return [(row[0], self.body_to_model([row[1], 0.0, row[2]])) for row in self.accel]

    def g_eff(self, t):
        ts = [r[0] for r in self.accel]
        a = np.array([np.interp(t, ts, [r[k] for r in self.accel]) for k in (1, 2)])
        return self.g0_model() - self.body_to_model([a[0], 0.0, a[1]])

    def rotor_sign(self):
        """Sign of the rotor angular velocity about +Z(model) for forward travel."""
        return -1.0 if self.side == 'right' else 1.0

    def rotor_omega_table(self):
        """[(t, omega_rotor [rad/s])], signed about +Z(model)."""
        if self.mechanisms == 'stopped':
            return [(0.0, 0.0), (max(1.0, self.duration), 0.0)]
        s = self.rotor_sign()
        return [(t, s * rpm * 2 * math.pi / 60.0) for t, rpm in self.motor_rpm]
