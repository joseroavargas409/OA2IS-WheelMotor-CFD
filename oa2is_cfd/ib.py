"""Immersed-boundary (volume penalisation) variant of the gear-case model.

Mesh: one static layered mesh of the whole oil cavity. Fixed parts (ring gears
with their regenerated teeth, spindle, hub, seals, sleeve, shaft) are
body-fitted walls. Moving parts (suns, planets, carriers) are immersed bodies
handled by the solver interIBFoamOA2IS from constant/immersedBodiesDict, which
is generated here from the SAME prisms (oa2is_cfd.model.GearCase.components)
and the SAME kinematics (oa2is_cfd.kinematics) as the overset variant.
"""
from __future__ import annotations

import json
import math
import shutil
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import Polygon

from . import foam
from .layered import extrude, mesh2d, write_polymesh
from .meshing import _solids, subdivide
from .model import GearCase, Z_CLOSURE, Z_COVER, circle
from .scenario import Scenario


@dataclass
class IBResolution:
    name: str = 'ib_coarse'
    lc_moving: float = 0.009        # 2D size where bodies move (r < r_moving)
    r_moving: float = 0.47
    lc_ring: float = 0.006          # on the (fixed) ring teeth
    lc_wall: float = 0.012
    lc_outer: float = 0.016
    dz_band: float = 0.013          # axial size in the stage bands
    dz_far: float = 0.03
    band_halfwidth: float = 0.15
    z_snap: float = 0.004


IB_PRESETS = {
    'ib_test': IBResolution(name='ib_test', lc_moving=0.018, lc_ring=0.010, lc_wall=0.018, lc_outer=0.024,
                            dz_band=0.025, dz_far=0.04, z_snap=0.006),
    'ib_mid': IBResolution(name='ib_mid', lc_moving=0.012, lc_ring=0.008, lc_wall=0.014, lc_outer=0.020,
                           dz_band=0.018, dz_far=0.035, z_snap=0.005),
    'ib_coarse': IBResolution(),
    'ib_medium': IBResolution(name='ib_medium', lc_moving=0.006, lc_ring=0.004, lc_wall=0.009, lc_outer=0.012,
                              dz_band=0.007, dz_far=0.02),
    'ib_fine': IBResolution(name='ib_fine', lc_moving=0.004, lc_ring=0.003, lc_wall=0.006, lc_outer=0.010,
                            dz_band=0.004, dz_far=0.015),
}


def fixed_solids(gc: GearCase):
    """Body-fitted fixed solids of the IB mesh. The ring-gear bodies start 6 mm
    beyond the tooth root; the toothed band is a STATIC immersed body (its
    outline crosses the other stage's outlines in plan view, which a single
    layered 2D mesh cannot conform to without degenerate cells)."""
    return list(gc.background_solids(for_volume=True))


def build_ib_mesh(case: Path, gc: GearCase, res: IBResolution, log=print):
    case = Path(case)
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    t0 = time.time()
    g = GearCase(gc.stages, h=res.lc_wall)
    prisms = fixed_solids(g)
    emb = []
    for p in prisms:
        teeth = any(len(hh) > 200 for hh in p.holes)
        for loop in [p.outline] + list(p.holes):
            loop = np.asarray(loop)
            if np.max(np.hypot(loop[:, 0], loop[:, 1])) > 0.5361:
                continue
            emb.append((loop, res.lc_ring if teeth else res.lc_wall, 0.012))
    import oa2is_cfd.layered as L
    L.LOOP_MIN_SEG, old = 0.45, L.LOOP_MIN_SEG     # embedded profiles resampled at ~0.45 lc
    try:
        xy, el = mesh2d(circle(0.536, h=res.lc_wall), emb, res.lc_wall, lc_max=res.lc_outer,
                        size_points=[(0.0, 0.0, res.r_moving, res.lc_moving)], node=True)
    finally:
        L.LOOP_MIN_SEG = old
    log(f'  2D: {len(el)} elements ({time.time() - t0:.1f}s)')
    zst = [s.z for s in gc.stages]
    dz = lambda z: res.dz_band if min(abs(z - a) for a in zst) < res.band_halfwidth else res.dz_far
    moving = [p for c in g.components() for p in c['solids']]
    brk = [Z_COVER, Z_CLOSURE] + [z for p in prisms + moving for z in (p.z0, p.z1) if Z_COVER < z < Z_CLOSURE]
    zl = subdivide(brk, dz, tol=res.z_snap)
    m = extrude(xy, el, zl, _solids(prisms, Z_COVER, Z_CLOSURE), 'hubWall', ('hubCover', 'inboardClosure'),
                remove_pockets=True)
    log(f'  3D: {m["ncells"]} cells, {len(zl) - 1} layers, {m["pockets_removed"]} cells in closed pockets '
        f'removed ({time.time() - t0:.1f}s)')
    write_polymesh(case, m, 'fluid', {})
    foam.write(case, 'system/controlDict', 'application none; startFrom startTime; startTime 0; stopAt endTime; '
                                           'endTime 0; deltaT 1; writeControl timeStep; writeInterval 1;')
    stats = dict(total_cells=m['ncells'], layers=len(zl) - 1, elements2d=len(el), resolution=asdict(res),
                 pocket_cells_removed=m['pockets_removed'],
                 backlash=[st.backlash for st in gc.stages], method='immersed_boundary')
    (case / 'mesh_stats.json').write_text(json.dumps(stats, indent=1))
    return stats


# ------------------------------------------------------------- bodies
def _polar_table(outline, center, n=4096):
    """r(theta) of a star-shaped outline about `center` (rest pose), by
    periodic interpolation of the vertices in polar coordinates."""
    d = np.asarray(outline, float) - np.asarray(center, float)
    th = np.mod(np.arctan2(d[:, 1], d[:, 0]), 2 * math.pi)
    r = np.hypot(d[:, 0], d[:, 1])
    o = np.argsort(th, kind='stable')
    t = np.linspace(0, 2 * math.pi, n, endpoint=False)
    return np.interp(t, th[o], r[o], period=2 * math.pi)


def _is_circle(loop, rel=2e-3):
    loop = np.asarray(loop, float)
    c = loop.mean(axis=0)
    r = np.hypot(*(loop - c).T)
    return (len(loop) >= 12 and (r.max() - r.min()) < rel * r.mean()), c, r.mean()


def bodies_dict(gc: GearCase, sc: Scenario, backlash=None):
    """immersedBodiesDict entries for every moving rigid body."""
    st = {s.id: s for s in gc.stages}
    W = sc.rotor_omega_table()
    tab = lambda k: '(' + ' '.join(f'({t:.6g} {k * w:.8g})' for t, w in W) + ')'
    bodies = {}
    for c in gc.components():
        name = c['body']
        S = st['High' if 'High' in name else 'Low']
        if name.startswith('sun'):
            orbit, spin, ctr = 0.0, S.input, (0.0, 0.0)
        elif name.startswith('planet'):
            orbit, spin, ctr = S.output, S.planet_rel, tuple(c['center'])
        else:
            orbit, spin, ctr = 0.0, S.output, (0.0, 0.0)
        b = bodies.setdefault(name, dict(center=ctr, orbit=orbit, spin=spin, prims=[]))
        for p in c['solids']:
            z = (p.z0, p.z1)
            ok, cc, rr = _is_circle(p.outline)
            if ok:
                rin = 0.0
                hole_rest = []
                for h in p.holes:
                    okh, ch, rh = _is_circle(h)
                    if okh and np.hypot(*(ch - cc)) < 1e-6:
                        rin = rh
                    else:
                        hole_rest.append((ch, rh))
                b['prims'].append(dict(type='annulus', z=z, c=tuple(cc - np.array(ctr)), rIn=rin, rOut=rr))
                for ch, rh in hole_rest:
                    b['prims'].append(dict(type='hole', z=z, c=tuple(ch - np.array(ctr)), rOut=rh))
            else:
                # hexagon pin heads -> area-equivalent disc; gears -> polar table
                outline = np.asarray(p.outline)
                if len(outline) <= 8:
                    cc = outline.mean(axis=0)
                    area = Polygon(outline).area
                    b['prims'].append(dict(type='annulus', z=z, c=tuple(cc - np.array(ctr)), rIn=0.0,
                                           rOut=math.sqrt(area / math.pi)))
                else:
                    b['prims'].append(dict(type='polar', z=z, r=_polar_table(outline, ctr)))
        b['W'] = (orbit, spin)
    from .gears import ring_space_profile
    for S_ in gc.stages:
        sp = ring_space_profile(S_.ring_spec())
        rr = _polar_table(sp, (0.0, 0.0))
        bodies['ring' + S_.id] = dict(center=(0.0, 0.0), W=(0.0, 0.0), prims=[
            dict(type='polarOut', z=(S_.z - .059, S_.z + .059), r=rr, rOut=gc.ring_zone_r_out(S_) + 0.01)])
    out = []
    for name, b in bodies.items():
        prims = []
        for i, p in enumerate(b['prims']):
            body = f'type {p["type"]}; z ({p["z"][0]:.7g} {p["z"][1]:.7g});'
            if 'c' in p:
                body += f' c ({p["c"][0]:.8g} {p["c"][1]:.8g});'
            if p['type'] in ('annulus', 'hole'):
                body += f' rIn {p.get("rIn", 0):.7g}; rOut {p["rOut"]:.7g};'
            if p['type'] == 'polarOut':
                body += f' rOut {p["rOut"]:.7g};'
            if p['type'] in ('polar', 'polarOut'):
                body += f' r List<scalar> {len(p["r"])} (' + ' '.join(f'{x:.6g}' for x in p['r']) + ');'
            prims.append(f'            p{i} {{ {body} }}')
        orbit, spin = b['W']
        out.append(f'''    {name}
    {{
        center ({b["center"][0]:.8g} {b["center"][1]:.8g});
        orbit {tab(orbit)};
        spin {tab(spin)};
        primitives
        {{
{chr(10).join(prims)}
        }}
    }}''')
    return 'penalty 1000;\ndisplaceOil true;\nsubSamples 2;\nbodies\n{\n' + '\n'.join(out) + '\n}'


def chi_at_rest(gc: GearCase, C):
    """Binary solid indicator of cell centres C (N,3) for the rest pose."""
    pts = shapely.points(C[:, :2])
    chi = np.zeros(len(C), bool)
    for c in gc.components() + gc.ring_zones():
        for p in c['solids']:
            poly = p.polygon()
            shapely.prepare(poly)
            sel = (C[:, 2] > p.z0) & (C[:, 2] < p.z1)
            if sel.any():
                chi[sel] |= shapely.contains(poly, pts[sel])
    return chi


# ------------------------------------------------------------- case
def write_ib_case(case: Path, mesh_dir: Path, sc: Scenario, gc: GearCase):
    from .case import ROT_WALLS, _table, patch_names
    case, mesh_dir = Path(case), Path(mesh_dir)
    if case.exists():
        shutil.rmtree(case)
    (case / 'constant').mkdir(parents=True)
    shutil.copytree(mesh_dir / 'constant/polyMesh', case / 'constant/polyMesh')
    (case / 'constant/polyMesh/cellZones').unlink(missing_ok=True)
    shutil.copy(mesh_dir / 'mesh_stats.json', case / 'mesh_stats.json')
    patches = patch_names(case)
    st = {s.id: s for s in gc.stages}
    ratios = dict(inp=1.0, mid=st['High'].output, out=st['Low'].output)
    W = sc.rotor_omega_table()
    num = dict(maxCo=1.0, maxAlphaCo=0.5, maxDeltaT=2e-3, nOuter=1, nCorr=3, writeFormat='binary')
    num.update(sc.numerics)

    def bf(kind):
        out = []
        for name, typ in patches:
            if kind == 'U':
                if name in ROT_WALLS:
                    k = {'in': 1.0, 'mid': ratios['mid'], 'out': ratios['out']}[ROT_WALLS[name]]
                    out.append(f'{name} {{ type rotatingWallVelocity; origin (0 0 0); axis (0 0 1); '
                               f'omega {_table([(t, k * w) for t, w in W])}; value uniform (0 0 0); }}')
                else:
                    out.append(f'{name} {{ type noSlip; }}')
            elif kind == 'p_rgh':
                out.append(f'{name} {{ type fixedFluxPressure; value uniform 0; }}')
            else:
                out.append(f'{name} {{ type zeroGradient; }}')
        return 'boundaryField\n{\n    ' + '\n    '.join(out) + '\n}'

    foam.write(case, '0/U', 'dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0 0 0);\n' + bf('U'), 'volVectorField')
    foam.write(case, '0/p_rgh', 'dimensions [1 -1 -2 0 0 0 0];\ninternalField uniform 0;\n' + bf('p_rgh'), 'volScalarField')
    foam.write(case, '0/alpha.oil', 'dimensions [0 0 0 0 0 0 0];\ninternalField uniform 0;\n' + bf('a'), 'volScalarField')
    foam.write(case, 'constant/transportProperties', f'''phases (oil air);
oil {{ transportModel Newtonian; nu {sc.nu_oil():.6g}; rho {sc.oil_rho:g}; }}
air {{ transportModel Newtonian; nu {sc.air_nu:.6g}; rho {sc.air_rho:g}; }}
sigma {sc.sigma:g};''')
    foam.write(case, 'constant/turbulenceProperties', 'simulationType laminar;')
    foam.write(case, 'constant/dynamicMeshDict', 'dynamicFvMesh staticFvMesh;')
    g0 = sc.g0_model()
    foam.write(case, 'constant/g', f'dimensions [0 1 -2 0 0 0 0];\nvalue {foam.vec(g0)};', 'uniformDimensionedVectorField')
    rows = sc.accel_model()
    if len(rows) == 1:
        rows = rows + [(max(10.0, sc.duration * 2), rows[0][1])]
    (case / 'constant/acceleration-terms.dat').write_text(
        '(\n' + '\n'.join(f'    ({t:.6g} ({foam.vec(a)} (0 0 0) (0 0 0)))' for t, a in rows) + '\n)\n')
    foam.write(case, 'constant/fvOptions', 'vehicleAcceleration\n{\n    type tabulatedAccelerationSource;\n    active yes;\n'
                                           '    timeDataFileName "constant/acceleration-terms.dat";\n}')
    foam.write(case, 'constant/immersedBodiesDict', bodies_dict(gc, sc, None))
    probes = ' '.join(foam.vec(s['position']) for s in sc.sensors)
    foam.write(case, 'system/controlDict', '\n'.join([
        'application interIBFoamOA2IS;', 'libs ("libfvOptions.so" "libsampling.so" "libfieldFunctionObjects.so");',
        'startFrom latestTime;', 'startTime 0;', 'stopAt endTime;', f'endTime {sc.duration:g};', 'deltaT 1e-4;',
        'writeControl adjustableRunTime;', f'writeInterval {sc.write_interval:g};', 'purgeWrite 0;',
        f'writeFormat {num["writeFormat"]};', 'writePrecision 8;', 'writeCompression off;', 'timeFormat general;',
        'timePrecision 8;', 'runTimeModifiable yes;', 'adjustTimeStep yes;', f'maxCo {num["maxCo"]};',
        f'maxAlphaCo {num["maxAlphaCo"]};', f'maxDeltaT {num["maxDeltaT"]};',
        'functions', '{',
        '    oilVolume { type volFieldValue; fields (alpha.oil); operation volIntegrate; regionType all;',
        '                writeFields false; writeControl timeStep; writeInterval 1; log false; }',
        '    extremes { type fieldMinMax; fields (U alpha.oil chi); writeControl timeStep; writeInterval 1; log false; }',
        f'    sensors {{ type probes; fields (alpha.oil U p p_rgh chi); writeControl timeStep; writeInterval 1;',
        f'              probeLocations ({probes}); }}',
        '}']))
    foam.write(case, 'system/fvSchemes', '''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(rhoPhi,U) Gauss linearUpwind grad(U);
    div(phi,alpha) Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes { default limited 0.5; }''')
    foam.write(case, 'system/fvSolution', f'''solvers
{{
    "alpha.oil.*" {{ nAlphaCorr 2; nAlphaSubCycles 1; cAlpha 1; MULESCorr yes; nLimiterIter 5;
                    solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0; }}
    "pcorr.*" {{ solver PCG; preconditioner DIC; tolerance 1e-7; relTol 0; }}
    p_rgh {{ solver GAMG; smoother DIC; tolerance 1e-7; relTol 0.01; }}
    p_rghFinal {{ $p_rgh; relTol 0; }}
    "U.*" {{ solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0; }}
}}
PIMPLE
{{
    momentumPredictor no; nOuterCorrectors {num["nOuter"]}; nCorrectors {num["nCorr"]};
    nNonOrthogonalCorrectors 1; pRefPoint (0.2 0.45 -1.3); pRefValue 0;
}}
relaxationFactors {{ equations {{ ".*" 1; }} }}''')
    meta = dict(sc.__dict__, g0_model=list(map(float, g0)), ratios=ratios,
                backlash=[s.backlash for s in gc.stages], method='immersed_boundary')
    (case / 'scenario.json').write_text(json.dumps(meta, indent=1, default=float))
    return patches


def initial_oil_ib(case: Path, sc: Scenario, gc: GearCase):
    """Flat free surface normal to g_eff(0) holding the scenario litres in the
    fluid space (outside fixed AND moving solids); alpha = 0 inside bodies."""
    from .case import read_scalar_field, read_vector_field
    case = Path(case)
    g = sc.g_eff(0.0)
    n = -g / np.linalg.norm(g)
    level = gc.level_for_volume(sc.oil_litres, normal=n)
    cd = case / 'system/controlDict'
    txt = cd.read_text()
    cd.write_text(txt.replace('writeFormat binary;', 'writeFormat ascii;'))
    foam.run(case, 'postProcess', '-func', 'writeCellCentres', '-time', '0', log='cellCentres')
    foam.run(case, 'postProcess', '-func', 'writeCellVolumes', '-time', '0', log='cellVolumes')
    C = read_vector_field(case / '0/C')
    V = read_scalar_field(case / '0/V')
    chi = chi_at_rest(gc, C)
    # level from the actual mesh (exact inventory in the discrete fluid volume)
    h = np.cbrt(V)
    proj = C @ n
    vol = lambda L: float(np.sum(np.clip(0.5 - (proj - L) / h, 0, 1) * (~chi) * V))
    lo, hi = proj.min() - 0.05, proj.max() + 0.05
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if vol(mid) < sc.oil_litres / 1000 else (lo, mid)
    level_geom, level = level, 0.5 * (lo + hi)
    s = proj - level
    alpha = np.clip(0.5 - s / h, 0, 1) * (~chi)
    t = (case / '0/alpha.oil').read_text()
    t = t.replace('internalField uniform 0;', 'internalField nonuniform List<scalar>\n%d\n(\n%s\n)\n;' %
                  (len(alpha), '\n'.join('%.6g' % a for a in alpha)))
    (case / '0/alpha.oil').write_text(t)
    cd.write_text(txt)
    meta = dict(level=float(level), level_geometric=float(level_geom), up=list(map(float, n)), litres=sc.oil_litres,
                oil_in_field_L=float(np.sum(alpha * V) * 1e3), fluid_L=float(np.sum(V * (~chi)) * 1e3))
    (case / 'initial_oil.json').write_text(json.dumps(meta, indent=1))
    return meta
