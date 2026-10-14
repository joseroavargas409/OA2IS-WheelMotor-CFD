"""Layered (2.5D) conformal mesher writing OpenFOAM polyMesh directly.

Every solid inside the wheel-motor gear case is a prism along the rotation
axis Z: a planar outline (gear profile, carrier plate with windows, pin,
shaft, ring, spindle plate...) extruded over a Z interval. A mesh region is
therefore described exactly by

* an outer loop (boundary of the region),
* solids ``Solid(name, outline, z0, z1, holes)``,
* Z levels (``zlevels``) that must contain every solid z0/z1.

Algorithm
1. gmsh meshes the outer loop in 2D with every solid outline *embedded*, so
   element edges follow all outlines (quads with recombination, else tris).
2. Each 2D element is classified by its centroid against each solid outline
   (shapely). The element column is extruded through the z levels; a cell
   (element, layer) is fluid unless its 2D centroid is inside a solid active
   over that layer.
3. Faces: internal between two fluid cells; wall patch (named after the
   solid) between fluid and solid cells; outer lateral faces and the z caps go
   to user-chosen patches. Output: points/faces/owner/neighbour/boundary and
   one cellZone, in OpenFOAM ASCII format.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import gmsh
import numpy as np
import shapely
from shapely.geometry import Polygon


@dataclass
class Solid:
    name: str
    outline: np.ndarray
    z0: float
    z1: float
    holes: list = field(default_factory=list)

    def polygon(self):
        return Polygon(np.asarray(self.outline), [np.asarray(h) for h in self.holes]).buffer(0)


# ---------------------------------------------------------------- 2D meshing
SIMPLIFY_FRAC = 0.03  # embedded curves simplified with tolerance SIMPLIFY_FRAC*lc (node=True path)
LOOP_MIN_SEG = 0.05   # embedded loops: drop vertices closer than this fraction of the local size


def clean_loop(pts, tol):
    """Drop consecutive points closer than `tol` (incl. the closing pair)."""
    pts = np.asarray(pts, dtype=float)
    keep = [0]
    for i in range(1, len(pts)):
        if np.hypot(*(pts[i] - pts[keep[-1]])) > tol:
            keep.append(i)
    if len(keep) > 2 and np.hypot(*(pts[keep[-1]] - pts[keep[0]])) <= tol:
        keep.pop()
    return pts[keep]


def _mesh2d(outer: np.ndarray, embedded: list, lc_outer: float, holes: list | None = None,
            recombine: bool = True, size_points: list[tuple] | None = None, lc_max: float | None = None,
            algorithm: int | None = None, node: bool = False):
    """2D mesh of the region bounded by `outer` (minus `holes`) whose element
    edges follow every `embedded` loop.

    embedded : list of (loop[N,2], lc_on_loop, grading_distance)
    holes    : list of (loop[N,2], lc_on_loop) removed from the region
    Size field: lc grows linearly from lc_on_loop on each loop to lc_max over
    grading_distance (Threshold fields combined with Min).
    Returns (nodes[N,2], list of CCW element node arrays (3 or 4 nodes))."""
    lc_max = lc_max or lc_outer
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.model.add('m2d')
    geo = gmsh.model.geo

    def loop_curves(pts, lc):
        pts = clean_loop(pts, LOOP_MIN_SEG * lc)
        p = [geo.addPoint(float(x), float(y), 0, lc) for x, y in pts]
        return [geo.addLine(p[i], p[(i + 1) % len(p)]) for i in range(len(p))]

    oc = loop_curves(outer, lc_outer)
    loops = [geo.addCurveLoop(oc)]
    groups = [(oc, lc_outer, 0.0)]
    for loop, lc in holes or []:
        hc = loop_curves(loop, lc)
        loops.append(geo.addCurveLoop(hc))
        groups.append((hc, lc, 0.0))
    surf = geo.addPlaneSurface(loops)
    emb = []
    # identical loops (e.g. the same circle used by several solids) must be embedded once
    uniq = {}
    for loop, lc, dist in embedded:
        loop = np.asarray(loop, float)
        key = (len(loop), round(float(loop[:, 0].mean()), 6), round(float(loop[:, 1].mean()), 6),
               round(float(np.hypot(*(loop[0] - loop.mean(axis=0)))), 6))
        if key in uniq:
            l0, lc0, d0 = uniq[key]
            uniq[key] = (l0, min(lc0, lc), max(d0, dist))
        else:
            uniq[key] = (loop, lc, dist)
    if not node:
        for loop, lc, dist in uniq.values():
            c = loop_curves(loop, lc)
            emb += c
            groups.append((c, lc, dist))
    else:
        # loops at different z may cross in plan view: split them at their
        # intersections so that gmsh gets a conforming network of curves
        import shapely
        from shapely.geometry import LineString
        from shapely.ops import unary_union
        loops = [(clean_loop(l, LOOP_MIN_SEG * lc), lc, dist) for l, lc, dist in uniq.values()]
        lines = [LineString(np.vstack([l, l[:1]])) for l, _, _ in loops]
        noded = unary_union(lines)
        segs = list(getattr(noded, 'geoms', [noded]))
        ptag = {}

        def pt(x, y, lc):
            k = (round(x, 9), round(y, 9))
            if k not in ptag:
                ptag[k] = geo.addPoint(x, y, 0, lc)
            return ptag[k]
        tree = shapely.STRtree(lines)
        bylc = {}
        for sg in segs:
            mid = sg.interpolate(0.5, normalized=True)
            j = int(tree.nearest(mid))
            lc, dist = loops[j][1], loops[j][2]
            # Douglas-Peucker: nearly collinear vertices create flat triangles
            sg = sg.simplify(SIMPLIFY_FRAC * lc, preserve_topology=False)
            cs = np.asarray(sg.coords)
            if len(cs) > 2:   # drop interior vertices too close to the segment ends (intersection nodes)
                tol = LOOP_MIN_SEG * lc
                keep = [0] + [i for i in range(1, len(cs) - 1)
                              if np.hypot(*(cs[i] - cs[0])) > tol and np.hypot(*(cs[i] - cs[-1])) > tol] + [len(cs) - 1]
                cs = cs[keep]
            for a, b in zip(cs[:-1], cs[1:]):
                if np.hypot(*(b - a)) < 1e-9:
                    continue
                c = geo.addLine(pt(a[0], a[1], lc), pt(b[0], b[1], lc))
                emb.append(c)
                bylc.setdefault((lc, dist), []).append(c)
        for (lc, dist), cs in bylc.items():
            groups.append((cs, lc, dist))
    geo.synchronize()
    if emb:
        gmsh.model.mesh.embed(1, emb, 2, surf)
    fields = []
    for curves, lc, dist in groups:
        fd = gmsh.model.mesh.field.add('Distance')
        gmsh.model.mesh.field.setNumbers(fd, 'CurvesList', curves)
        gmsh.model.mesh.field.setNumber(fd, 'Sampling', 6)
        th = gmsh.model.mesh.field.add('Threshold')
        gmsh.model.mesh.field.setNumber(th, 'InField', fd)
        gmsh.model.mesh.field.setNumber(th, 'SizeMin', lc)
        gmsh.model.mesh.field.setNumber(th, 'SizeMax', lc_max)
        gmsh.model.mesh.field.setNumber(th, 'DistMin', 0.0)
        gmsh.model.mesh.field.setNumber(th, 'DistMax', max(dist, 1e-6) if dist > 0 else 1e-6)
        fields.append(th)
    for (x, y, r, lc) in size_points or []:
        b = gmsh.model.mesh.field.add('Ball')
        for k, v in dict(XCenter=x, YCenter=y, ZCenter=0, Radius=r, VIn=lc, VOut=lc_max).items():
            gmsh.model.mesh.field.setNumber(b, k, v)
        fields.append(b)
    mn = gmsh.model.mesh.field.add('Min')
    gmsh.model.mesh.field.setNumbers(mn, 'FieldsList', fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(mn)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
    gmsh.option.setNumber('Mesh.Algorithm', algorithm or (8 if recombine else 6))
    gmsh.option.setNumber('Mesh.RecombinationAlgorithm', 1)
    if recombine:
        gmsh.model.mesh.setRecombine(2, surf)
    gmsh.model.mesh.generate(2)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xy = coords.reshape(-1, 3)[:, :2]
    idx = {int(t): i for i, t in enumerate(tags)}
    elems = []
    etypes, _, enodes = gmsh.model.mesh.getElements(2, surf)
    for et, en in zip(etypes, enodes):
        nn = {2: 3, 3: 4}[int(et)]
        for row in np.asarray(en, dtype=np.int64).reshape(-1, nn):
            elems.append(np.array([idx[int(t)] for t in row]))
    cons = set()
    for dim, tag in gmsh.model.getEntities(1):
        _, _, en = gmsh.model.mesh.getElements(1, tag)
        for arr in en:
            a = np.asarray(arr, dtype=np.int64).reshape(-1, 2)
            for u, v in a:
                iu, iv = idx[int(u)], idx[int(v)]
                cons.add((min(iu, iv), max(iu, iv)))
    gmsh.finalize()
    elems = improve_flat(xy, elems, cons)
    used = np.unique(np.concatenate(elems))
    remap = -np.ones(len(xy), dtype=np.int64)
    remap[used] = np.arange(len(used))
    xy = xy[used]
    elems = [remap[e] for e in elems]
    # gmsh may return a few tangled (inverted) elements near very short
    # embedded segments: detect them against the majority orientation.
    sgn = np.sign([_area(xy[e]) for e in elems])
    major = 1.0 if np.sum(sgn > 0) >= np.sum(sgn < 0) else -1.0
    n_inv = int(np.sum(sgn != major))
    if n_inv:
        raise TangledMesh(f'{n_inv} inverted 2D elements')
    if major < 0:
        elems = [e[::-1] for e in elems]
    return xy, fix_elements(xy, elems)


class TangledMesh(RuntimeError):
    pass


def mesh2d(*args, **kw):
    """Robust 2D meshing: tries quad-dominant frontal (8), frontal-Delaunay
    with recombination (6), then pure triangles (6, 5) until no element is
    inverted."""
    attempts = [dict(recombine=True, algorithm=8), dict(recombine=True, algorithm=6),
                dict(recombine=False, algorithm=6), dict(recombine=False, algorithm=5)]
    if not kw.pop('recombine', True):
        attempts = attempts[2:]
    err = None
    for a in attempts:
        try:
            return _mesh2d(*args, **kw, **a)
        except TangledMesh as e:
            err = e
    raise err


def _angles(p):
    n = len(p)
    out = []
    for i in range(n):
        u = p[i - 1] - p[i]
        v = p[(i + 1) % n] - p[i]
        c = u @ v / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-300)
        out.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return np.array(out)


def improve_flat(xy, elems, cons, qmax=160.0, sweeps=4):
    """Remove nearly flat elements (an interior angle > qmax degrees):
    quads are split from the obtuse vertex into two triangles; triangles get
    their longest edge flipped with the neighbour (unconstrained edges only)
    when that improves the worst angle. Prevents the near-singular cells that
    break fvc::reconstruct."""
    out = []
    for e in elems:
        e = np.asarray(e)
        if len(e) == 4:
            ang = _angles(xy[e])
            if ang.max() > qmax:
                k = int(np.argmax(ang))
                e = np.roll(e, -k)
                out += [e[[0, 1, 2]], e[[0, 2, 3]]]
                continue
        out.append(e)
    elems = out
    for _ in range(sweeps):
        edge = {}
        for ei, e in enumerate(elems):
            m = len(e)
            for j in range(m):
                a, b = int(e[j]), int(e[(j + 1) % m])
                edge.setdefault((min(a, b), max(a, b)), []).append(ei)
        changed = 0
        done = set()
        for ei, e in enumerate(elems):
            if len(e) != 3 or ei in done:
                continue
            ang = _angles(xy[e])
            if ang.max() <= qmax:
                continue
            k = int(np.argmax(ang))            # obtuse vertex; opposite edge is the longest
            a, b = int(e[(k + 1) % 3]), int(e[(k + 2) % 3])
            key = (min(a, b), max(a, b))
            if key in cons or len(edge.get(key, [])) != 2:
                continue
            ej = [x for x in edge[key] if x != ei][0]
            f = elems[ej]
            if len(f) != 3 or ej in done:
                continue
            c = int(e[k])
            d = [int(x) for x in f if int(x) not in (a, b)][0]
            t1, t2 = np.array([c, a, d]), np.array([c, d, b])
            if _area(xy[t1]) <= 0 or _area(xy[t2]) <= 0:
                t1, t2 = t1[::-1], t2[::-1]
            if _area(xy[t1]) <= 0 or _area(xy[t2]) <= 0:
                continue
            new_worst = max(_angles(xy[t1]).max(), _angles(xy[t2]).max())
            if new_worst < max(ang.max(), _angles(xy[f]).max()):
                elems[ei], elems[ej] = t1, t2
                done.update((ei, ej))
                changed += 1
        if not changed:
            break
    return elems


def _area(p):
    return 0.5 * np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1])


def fix_elements(xy, elems, min_area=1e-12):
    """Orient elements counter-clockwise; split non-convex/bow-tie quads into
    two triangles (choosing the diagonal that keeps both positive)."""
    out, nsplit = [], 0
    for e in elems:
        e = np.asarray(e)
        if _area(xy[e]) < 0:
            e = e[::-1]
        if len(e) == 4:
            p = xy[e]
            d = np.roll(p, -1, axis=0) - p
            cross = d[:, 0] * np.roll(d, -1, axis=0)[:, 1] - d[:, 1] * np.roll(d, -1, axis=0)[:, 0]
            if np.all(cross > 0):
                out.append(e)
                continue
            nsplit += 1
            for tri in ((e[[0, 1, 2]], e[[0, 2, 3]]), (e[[0, 1, 3]], e[[1, 2, 3]])):
                if all(_area(xy[t]) > min_area for t in tri):
                    out += list(tri)
                    break
            else:
                raise ValueError('degenerate quad cannot be split: %s' % xy[e].tolist())
        else:
            if _area(xy[e]) <= min_area:
                raise ValueError('degenerate triangle: %s' % xy[e].tolist())
            out.append(e)
    return out


# --------------------------------------------------------------- extrusion
def extrude(xy, elems, zlevels, solids: list[Solid], outer_patch: str, cap_patches=('caps', 'caps'),
            outer_tol=1e-7, outer_polygon: Polygon | None = None, remove_pockets: bool = False,
            pocket_patch: str = 'closedPocket'):
    """Build polyMesh arrays. Returns dict(points, faces, owner, neighbour, patches, ncells, cellmask)."""
    zl = np.asarray(zlevels, dtype=float)
    nz = len(zl) - 1
    n2 = len(xy)
    ne = len(elems)
    cent = np.array([xy[e].mean(axis=0) for e in elems])
    zmid = 0.5 * (zl[:-1] + zl[1:])
    # solid id per (elem, layer); -1 fluid. Later solids override earlier ones.
    sid = -np.ones((ne, nz), dtype=np.int64)
    pts = shapely.points(cent)
    for k, s in enumerate(solids):
        poly = s.polygon()
        shapely.prepare(poly)
        inside = shapely.contains(poly, pts) | shapely.touches(poly, pts)
        lay = (zmid > s.z0) & (zmid < s.z1)
        if inside.any() and lay.any():
            sid[np.ix_(inside, lay)] = k
    fluid = sid < 0
    # keep only the largest connected fluid region: closed pockets (e.g. a
    # gap sealed between two fixed parts) make the pressure equation singular
    pocket_id = len(solids)
    if remove_pockets:
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        idx = -np.ones((ne, nz), dtype=np.int64)
        idx[fluid] = np.arange(fluid.sum())
        rows, cols = [], []
        emap = {}
        for ei, e in enumerate(elems):
            m_ = len(e)
            for j in range(m_):
                a, b = int(e[j]), int(e[(j + 1) % m_])
                emap.setdefault((min(a, b), max(a, b)), []).append(ei)
        for lst in emap.values():
            if len(lst) == 2:
                e0, e1 = lst
                both = fluid[e0] & fluid[e1]
                rows.append(idx[e0][both]); cols.append(idx[e1][both])
        vert = fluid[:, :-1] & fluid[:, 1:]
        rows.append(idx[:, :-1][vert]); cols.append(idx[:, 1:][vert])
        r = np.concatenate(rows); c = np.concatenate(cols)
        nf = int(fluid.sum())
        g = coo_matrix((np.ones(len(r)), (r, c)), shape=(nf, nf))
        ncomp, lab = connected_components(g, directed=False)
        if ncomp > 1:
            main = np.argmax(np.bincount(lab))
            drop = np.zeros((ne, nz), bool)
            drop[fluid] = lab != main
            sid[drop] = pocket_id
            fluid = sid < 0
            pockets_removed = int(drop.sum())
        else:
            pockets_removed = 0
    else:
        pockets_removed = 0
    cid = -np.ones((ne, nz), dtype=np.int64)
    cid[fluid] = np.arange(fluid.sum())
    ncells = int(fluid.sum())
    P = np.vstack([np.c_[xy, np.full(n2, z)] for z in zl])      # node (i, level) -> i + level*n2

    # 2D edges: map edge -> list of (elem, local index)
    edge_map = {}
    for ei, e in enumerate(elems):
        m = len(e)
        for j in range(m):
            a, b = int(e[j]), int(e[(j + 1) % m])
            key = (a, b) if a < b else (b, a)
            edge_map.setdefault(key, []).append((ei, a, b))

    faces, owner, neigh, kind = [], [], [], []   # kind: -1 internal, else patch key

    def add(face, o, n, k):
        faces.append(face); owner.append(o); neigh.append(n); kind.append(k)

    patch_names = [outer_patch, cap_patches[0], cap_patches[1]] + [s.name for s in solids] + [pocket_patch]
    pid = {}
    for nm in patch_names:
        pid.setdefault(nm, len(pid))
    solid_pid = [pid[s.name] for s in solids] + [pid[pocket_patch]]

    # lateral faces (one per 2D edge per layer)
    for key, lst in edge_map.items():
        for l in range(nz):
            if len(lst) == 2:
                (e0, a0, b0), (e1, _, _) = lst
                c0, c1 = cid[e0, l], cid[e1, l]
                # quad oriented so its normal points out of e0 (CCW elem: edge a->b, outward = right side)
                q = [a0 + l * n2, b0 + l * n2, b0 + (l + 1) * n2, a0 + (l + 1) * n2]
                if c0 >= 0 and c1 >= 0:
                    if c0 < c1:
                        add(q, c0, c1, -1)
                    else:
                        add(q[::-1], c1, c0, -1)
                elif c0 >= 0:
                    add(q, c0, -1, solid_pid[sid[e1, l]])
                elif c1 >= 0:
                    add(q[::-1], c1, -1, solid_pid[sid[e0, l]])
            else:
                (e0, a0, b0), = lst
                c0 = cid[e0, l]
                if c0 >= 0:
                    q = [a0 + l * n2, b0 + l * n2, b0 + (l + 1) * n2, a0 + (l + 1) * n2]
                    add(q, c0, -1, pid[outer_patch])
    # horizontal faces at each level
    for ei, e in enumerate(elems):
        for lev in range(nz + 1):
            below = cid[ei, lev - 1] if lev > 0 else -1
            above = cid[ei, lev] if lev < nz else -1
            ring = [int(v) + lev * n2 for v in e]           # CCW seen from +z => normal +z
            if below >= 0 and above >= 0:
                add(ring, below, above, -1)
            elif below >= 0:
                k = pid[cap_patches[1]] if lev == nz else solid_pid[sid[ei, lev]]
                add(ring, below, -1, k)
            elif above >= 0:
                k = pid[cap_patches[0]] if lev == 0 else solid_pid[sid[ei, lev - 1]]
                add(ring[::-1], above, -1, k)
    owner = np.asarray(owner); neigh = np.asarray(neigh); kind = np.asarray(kind)
    internal = np.where(kind < 0)[0]
    internal = internal[np.lexsort((neigh[internal], owner[internal]))]
    order = [internal]
    patches = []
    start = len(internal)
    for nm, p in pid.items():
        f = np.where(kind == p)[0]
        f = f[np.argsort(owner[f], kind='stable')]
        if len(f):
            patches.append((nm, start, len(f)))
            start += len(f)
            order.append(f)
    order = np.concatenate(order)
    faces = [faces[i] for i in order]
    owner = owner[order]
    neigh = neigh[order][:len(internal)]
    # compact points (remove unused)
    used = np.unique(np.concatenate([np.asarray(f) for f in faces]))
    remap = -np.ones(len(P), dtype=np.int64)
    remap[used] = np.arange(len(used))
    P = P[used]
    faces = [remap[np.asarray(f)] for f in faces]
    return dict(points=P, faces=faces, owner=owner, neighbour=neigh, patches=patches,
                ncells=ncells, fluid=fluid, sid=sid, ne=ne, nz=nz, pockets_removed=pockets_removed)


# ------------------------------------------------------------------ merging
def merge(meshes: list[tuple[str, dict]]):
    """Concatenate independent zone meshes into one polyMesh WITHOUT merging
    coincident points/faces (OpenFOAM's mergeMeshes stitches coincident
    boundary faces, which ties a rotating zone to the fixed background).
    Returns (mesh dict, zones {name: (start, count)})."""
    pts, faces, own, nei = [], [], [], []
    internal_faces, patch_faces = [], {}
    p_off = c_off = 0
    zones = {}
    for name, m in meshes:
        nint = len(m['neighbour'])
        F = m['faces']
        for i in range(nint):
            internal_faces.append((m['owner'][i] + c_off, m['neighbour'][i] + c_off, F[i] + p_off))
        for pname, start, n in m['patches']:
            lst = patch_faces.setdefault(pname, [])
            for i in range(start, start + n):
                lst.append((m['owner'][i] + c_off, F[i] + p_off))
        pts.append(m['points'])
        zones[name] = (c_off, m['ncells'])
        p_off += len(m['points'])
        c_off += m['ncells']
    # internal faces are already upper-triangular within each block and blocks are
    # ordered by cell offset, so concatenation keeps the ordering valid
    out_faces = [f for _, _, f in internal_faces]
    out_own = [o for o, _, _ in internal_faces]
    out_nei = [n for _, n, _ in internal_faces]
    patches = []
    start = len(out_faces)
    for pname, lst in patch_faces.items():
        lst.sort(key=lambda t: t[0])
        patches.append((pname, start, len(lst)))
        out_faces += [f for _, f in lst]
        out_own += [o for o, _ in lst]
        start += len(lst)
    return dict(points=np.vstack(pts), faces=out_faces, owner=np.asarray(out_own), neighbour=np.asarray(out_nei),
                patches=patches, ncells=c_off), zones


# ------------------------------------------------------------------ writer
def _hdr(cls, obj, note=''):
    return ('FoamFile\n{\n    version 2.0;\n    format ascii;\n    class %s;\n    location "constant/polyMesh";\n'
            '    object %s;%s\n}\n' % (cls, obj, ('\n    note "%s";' % note) if note else ''))


def write_polymesh(case: Path, m: dict, zone: str, patch_types: dict[str, str], zones_override=None):
    d = Path(case) / 'constant' / 'polyMesh'
    d.mkdir(parents=True, exist_ok=True)
    P = m['points']
    with open(d / 'points', 'w') as f:
        f.write(_hdr('vectorField', 'points'))
        f.write('%d\n(\n' % len(P))
        f.write('\n'.join('(%.9g %.9g %.9g)' % tuple(p) for p in P))
        f.write('\n)\n')
    with open(d / 'faces', 'w') as f:
        f.write(_hdr('faceList', 'faces'))
        f.write('%d\n(\n' % len(m['faces']))
        f.write('\n'.join('%d(%s)' % (len(fc), ' '.join(map(str, fc))) for fc in m['faces']))
        f.write('\n)\n')
    note = 'nPoints:%d nCells:%d nFaces:%d nInternalFaces:%d' % (len(P), m['ncells'], len(m['faces']), len(m['neighbour']))
    for name, arr in (('owner', m['owner']), ('neighbour', m['neighbour'])):
        with open(d / name, 'w') as f:
            f.write(_hdr('labelList', name, note))
            f.write('%d\n(\n' % len(arr))
            f.write('\n'.join(map(str, arr)))
            f.write('\n)\n')
    with open(d / 'boundary', 'w') as f:
        f.write(_hdr('polyBoundaryMesh', 'boundary'))
        f.write('%d\n(\n' % len(m['patches']))
        for name, start, n in m['patches']:
            t = patch_types.get(name, 'wall')
            f.write('    %s\n    {\n        type %s;\n%s        nFaces %d;\n        startFace %d;\n    }\n' %
                    (name, t, '        inGroups List<word> 1(%s);\n' % t if t in ('wall', 'overset') else '', n, start))
        f.write(')\n')
    zones = zones_override or {zone: (0, m['ncells'])}
    with open(d / 'cellZones', 'w') as f:
        f.write(_hdr('regIOobject', 'cellZones'))
        f.write('%d\n(\n' % len(zones))
        for zname, (c0, n) in zones.items():
            f.write('%s\n{\n    type cellZone;\ncellLabels List<label> %d\n(\n' % (zname, n))
            f.write('\n'.join(map(str, range(c0, c0 + n))))
            f.write('\n)\n;\n}\n')
        f.write(')\n')
