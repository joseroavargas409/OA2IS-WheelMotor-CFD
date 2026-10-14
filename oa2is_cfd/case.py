"""Write an overInterDyMFoam case (OpenFOAM v1912) for a scenario.

Mechanisms (multiSolidBodyMotionSolver, one entry per moving cellZone):
  sun<S>            rotatingMotion about Z, omega = ratio_sun * omega_rotor(t)
  planet<S><j>      multiMotion: spin about its rest-pose centre (relative
                    rate) then orbit about Z with the carrier
  carrier<S>In/Out  rotatingMotion about Z at the carrier rate
  ring<S>, background : stationary
Axisymmetric rotating walls in the background use rotatingWallVelocity with
the same omega tables (hub, hub cover, low-carrier ring at the output rate;
intermediate sleeve at 1/5; rotor shaft at 1).
"""
from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

import numpy as np

from . import foam
from .kinematics import stages as make_stages
from .scenario import Scenario

MOVING_PREFIX = ('sun', 'planet', 'carrier')
ROT_WALLS = {'hubWall': 'out', 'hubCover': 'out', 'carrierLowRing': 'out', 'intermediate': 'mid', 'shaft': 'in'}


def _table(rows):
    return 'table (' + ' '.join(f'({t:.6g} {w:.8g})' for t, w in rows) + ')'


def zone_names(case):
    txt = (Path(case) / 'constant/polyMesh/cellZones').read_text()
    return re.findall(r'\n(\w+)\s*\n\{\s*\n\s*type\s+cellZone', txt)


def patch_names(case):
    b = (Path(case) / 'constant/polyMesh/boundary').read_text()
    return re.findall(r'\n\s+(\w+)\s*\n\s*\{\s*\n\s*type\s+(\w+);', b)


def motion_entries(sc: Scenario, backlash):
    st_hi, st_lo = make_stages(backlash)
    W = sc.rotor_omega_table()
    scale = lambda k: [(t, k * w) for t, w in W]
    rot = lambda origin, k: ('solidBodyMotionFunction rotatingMotion; rotatingMotionCoeffs { origin %s; axis (0 0 1); omega %s; }'
                             % (foam.vec(origin), _table(scale(k))))
    E = {}
    for st in (st_hi, st_lo):
        E['sun' + st.id] = rot((0, 0, 0), st.input)
        for j in range(3):
            c = st.planet_center(j, 0)
            # multiMotion composes TR = TR_0 * TR_1 (septernion); verified by
            # tests/test_kinematics_openfoam.py against the analytic pose.
            E[f'planet{st.id}{j}'] = ('solidBodyMotionFunction multiMotion; multiMotionCoeffs { orbit { %s } spin { %s } }'
                                     % (rot((0, 0, 0), st.output), rot((c[0], c[1], 0), st.planet_rel)))
        for side in ('In', 'Out'):
            E[f'carrier{st.id}{side}'] = rot((0, 0, 0), st.output)
    ratios = dict(inp=1.0, mid=st_hi.output, out=st_lo.output)
    return E, ratios, W


def write_case(case: Path, mesh_dir: Path, sc: Scenario, backlash=(0.004, 0.003), n_procs=1):
    case, mesh_dir = Path(case), Path(mesh_dir)
    if case.exists():
        shutil.rmtree(case)
    (case / 'constant').mkdir(parents=True)
    shutil.copytree(mesh_dir / 'constant/polyMesh', case / 'constant/polyMesh')
    zones = zone_names(case)
    zones = ['background'] + [z for z in zones if z != 'background']
    patches = patch_names(case)
    E, ratios, W = motion_entries(sc, backlash)
    num = dict(maxCo=2.0, maxAlphaCo=0.6, maxDeltaT=1e-3, nOuter=2, nCorr=3, searchBox=(250, 250, 120),
               writeFormat='binary', purgeWrite=0)
    num.update(sc.numerics)

    # ---------------- boundary conditions
    def bf(kind):
        out = ['#includeEtc "caseDicts/setConstraintTypes"']
        for name, typ in patches:
            if typ == 'overset':
                v = {'U': 'uniform (0 0 0)'}.get(kind, 'uniform 0')
                out.append(f'{name} {{ type overset; value {v}; }}')
                continue
            if kind == 'U':
                if name.startswith(MOVING_PREFIX) and name != 'carrierLowRing':
                    out.append(f'{name} {{ type movingWallVelocity; value uniform (0 0 0); }}')
                elif name in ROT_WALLS:
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
    foam.write(case, '0/alpha.oil', 'dimensions [0 0 0 0 0 0 0];\ninternalField uniform 0;\n' + bf('alpha'), 'volScalarField')
    foam.write(case, '0/zoneID', 'dimensions [0 0 0 0 0 0 0];\ninternalField uniform 0;\n' + bf('alpha'), 'volScalarField')
    foam.write(case, '0/pointDisplacement', 'dimensions [0 1 0 0 0 0 0];\ninternalField uniform (0 0 0);\n'
               'boundaryField\n{\n    #includeEtc "caseDicts/setConstraintTypes"\n    ".*" { type calculated; value uniform (0 0 0); }\n}',
               'pointVectorField')

    # ---------------- physical properties
    foam.write(case, 'constant/transportProperties', f'''phases (oil air);
oil {{ transportModel Newtonian; nu {sc.nu_oil():.6g}; rho {sc.oil_rho:g}; }}
air {{ transportModel Newtonian; nu {sc.air_nu:.6g}; rho {sc.air_rho:g}; }}
sigma {sc.sigma:g};''')
    foam.write(case, 'constant/turbulenceProperties', 'simulationType laminar;')
    g0 = sc.g0_model()
    foam.write(case, 'constant/g', f'dimensions [0 1 -2 0 0 0 0];\nvalue {foam.vec(g0)};', 'uniformDimensionedVectorField')
    # vehicle acceleration (model frame) for tabulatedAccelerationSource: (t ((a) (Omega) (dOmega/dt)))
    rows = sc.accel_model()
    if len(rows) == 1:
        rows = rows + [(max(10.0, sc.duration * 2), rows[0][1])]
    acc = '(\n' + '\n'.join(f'    ({t:.6g} ({foam.vec(a)} (0 0 0) (0 0 0)))' for t, a in rows) + '\n)\n'
    (case / 'constant/acceleration-terms.dat').write_text(acc)
    foam.write(case, 'constant/fvOptions', '''vehicleAcceleration
{
    type            tabulatedAccelerationSource;
    active          yes;
    timeDataFileName "constant/acceleration-terms.dat";
}''')
    body = '\n    '.join(f'{z}\n    {{\n        {E[z]}\n    }}' for z in zones if z in E)
    foam.write(case, 'constant/dynamicMeshDict', f'''dynamicFvMesh dynamicOversetFvMesh;
solver multiSolidBodyMotionSolver;
multiSolidBodyMotionSolverCoeffs
{{
    {body}
}}''')

    # ---------------- system
    sb = num['searchBox']
    foam.write(case, 'system/controlDict', '\n'.join([
        'application overInterDyMFoamOA2IS;', 'libs ("liboverset.so" "libfvOptions.so");', 'startFrom latestTime;', 'startTime 0;',
        'stopAt endTime;', f'endTime {sc.duration:g};', 'deltaT 1e-5;', 'writeControl adjustableRunTime;',
        f'writeInterval {sc.write_interval:g};', f'purgeWrite {num["purgeWrite"]};', f'writeFormat {num["writeFormat"]};',
        'writePrecision 8;', 'writeCompression off;', 'timeFormat general;', 'timePrecision 8;', 'runTimeModifiable yes;',
        'adjustTimeStep yes;', f'maxCo {num["maxCo"]};', f'maxAlphaCo {num["maxAlphaCo"]};', f'maxDeltaT {num["maxDeltaT"]};',
        'functions', '{',
        '    oilVolumeAll { type volFieldValue; libs ("libfieldFunctionObjects.so"); fields (alpha.oil); operation volIntegrate;',
        '                   regionType all; writeFields false; writeControl timeStep; writeInterval 1; log false; }',
        '    extremes { type fieldMinMax; libs ("libfieldFunctionObjects.so"); fields (U alpha.oil); writeControl timeStep;',
        '               writeInterval 1; log false; }',
        '}']))
    foam.write(case, 'system/fvSchemes', f'''ddtSchemes {{ default Euler; }}
gradSchemes {{ default Gauss linear; }}
divSchemes
{{
    default none;
    div(rhoPhi,U) Gauss linearUpwind grad(U);
    div(phi,alpha) Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}}
laplacianSchemes {{ default Gauss linear limited 0.5; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default limited 0.5; }}
oversetInterpolation {{ method inverseDistance; searchBoxDivisions ({sb[0]} {sb[1]} {sb[2]}); }}
fluxRequired {{ default no; pcorr; p_rgh; alpha.oil; }}
oversetInterpolationSuppressed {{ grad(p_rgh); surfaceIntegrate(phiHbyA); }}''')
    foam.write(case, 'system/fvSolution', f'''solvers
{{
    "alpha.oil.*" {{ nAlphaCorr 2; nAlphaSubCycles 1; cAlpha 1; MULESCorr yes; nLimiterIter 5;
                    solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0; }}
    p_rgh {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-7; relTol 0.01; maxIter 300; }}
    p_rghFinal {{ $p_rgh; relTol 0; }}
    pcorr {{ solver PCG; preconditioner DIC; tolerance 1e-7; relTol 0; maxIter 300; }}
    pcorrFinal {{ $pcorr; relTol 0; }}
    "U.*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-7; relTol 0; }}
}}
PIMPLE
{{
    momentumPredictor no; correctPhi no; oversetAdjustPhi yes;
    nOuterCorrectors {num["nOuter"]}; nCorrectors {num["nCorr"]}; nNonOrthogonalCorrectors 0; ddtCorr true;
    pRefPoint (0.2 0.45 -1.3); pRefValue 0;
}}
relaxationFactors {{ equations {{ ".*" 1; }} }}''')
    sets = '\n'.join(f'    zoneToCell {{ zone {z}; fieldValues ( volScalarFieldValue zoneID {i} ); }}' for i, z in enumerate(zones))
    foam.write(case, 'system/setFieldsDict', 'defaultFieldValues ( volScalarFieldValue zoneID 0 );\nregions\n(\n' + sets + '\n);')
    if n_procs > 1:
        foam.write(case, 'system/decomposeParDict', f'numberOfSubdomains {n_procs};\nmethod scotch;')
    (case / 'scenario.json').write_text(json.dumps(dict(sc.__dict__, g0_model=list(g0), zones=zones,
                                                        ratios=ratios, backlash=list(backlash)), indent=1, default=float))
    return zones


def set_initial_oil(case: Path, sc: Scenario, level: float, normal):
    """alpha = 1 below the plane n.x = level (n = -g_eff/|g_eff|), smeared over one
    cell size to reduce the volume error of a binary initialisation."""
    case = Path(case)
    foam.run(case, 'postProcess', '-func', 'writeCellCentres', '-time', '0', log='cellCentres')
    foam.run(case, 'postProcess', '-func', 'writeCellVolumes', '-time', '0', log='cellVolumes')
    C = read_vector_field(case / '0/C')
    V = read_scalar_field(case / '0/V')
    n = np.asarray(normal, float) / np.linalg.norm(normal)
    d = np.cbrt(V)
    s = C @ n - level
    alpha = np.clip(0.5 - s / d, 0.0, 1.0)
    txt = (case / '0/alpha.oil').read_text()
    body = 'nonuniform List<scalar>\n%d\n(\n%s\n)\n' % (len(alpha), '\n'.join('%.6g' % a for a in alpha))
    txt = txt.replace('internalField uniform 0;', 'internalField ' + body + ';')
    (case / '0/alpha.oil').write_text(txt)
    return alpha, C, V


def _read_list(path, ncomp):
    s = Path(path).read_text()
    i = s.index('internalField')
    j = s.index('(', s.index('\n', i))
    n = int(s[i:j].split()[-1])
    k = s.index('\n)', j)
    data = s[j + 1:k].replace('(', ' ').replace(')', ' ').split()
    return np.array(data, float).reshape(n, ncomp) if ncomp > 1 else np.array(data, float)[:n]


def read_vector_field(path):
    return _read_list(path, 3)


def read_scalar_field(path):
    return _read_list(path, 1)
