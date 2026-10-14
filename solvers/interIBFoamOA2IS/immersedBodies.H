/*---------------------------------------------------------------------------*\
  OA2IS-WheelMotor-CFD : immersed moving rigid bodies for interIBFoamOA2IS
  License: GPL-3.0-or-later (derived work of OpenFOAM v1912, (C) OpenCFD Ltd.)

  Volume-penalisation immersed boundary (Brinkman type) for rigid bodies that
  are prisms along Z: each body is a union of primitives minus 'hole'
  primitives, given in the body frame at the rest pose:
    polar   : star-shaped outline r(theta) (gear profiles), table of N radii
    polarOut: solid outside r(theta) up to rOut (internal ring gears, static)
    annulus : rIn < |x - c| < rOut   (shafts, pins, plates, hubs)
    hole    : |x - c| < rOut removed from the body (carrier windows)
  Pose: the body centre c0 orbits about the Z axis (angle phi(t), carrier),
  the body spins about its own centre (angle psi(t), relative to the orbit
  frame). Both angular velocities are piecewise-linear tables
  ((t omega) ...) in rad/s, integrated exactly.

  Every time step:
    chi : solid volume fraction of each cell (2x2x2 sub-samples)
    Us  : rigid-body velocity at the cell centre
  Oil displacement: oil in the part of a cell newly covered by a body is moved
  to the neighbouring fluid cells (conservative, one pass), so the bodies push
  oil aside instead of swallowing it.
\*---------------------------------------------------------------------------*/

#ifndef immersedBodiesOA2IS_H
#define immersedBodiesOA2IS_H

#include "fvCFD.H"
#include "Tuple2.H"
#include "IOdictionary.H"

namespace Foam
{

class ibPrimitive
{
public:
    enum kinds { POLAR, POLAROUT, ANNULUS, HOLE };
    kinds kind;
    scalar z0, z1;
    vector2D c;
    scalar rIn, rOut;
    scalarField rTab;   // POLAR: r(theta_k), theta_k = 2 pi k/N
    scalar rMax;

    ibPrimitive(const dictionary& d)
    {
        const word t(d.get<word>("type"));
        const vector2D z(d.get<vector2D>("z"));
        z0 = z.x(); z1 = z.y();
        c = d.getOrDefault<vector2D>("c", vector2D(0, 0));
        rIn = d.getOrDefault<scalar>("rIn", 0);
        rOut = d.getOrDefault<scalar>("rOut", 0);
        if (t == "polar")
        {
            kind = POLAR;
            rTab = scalarField(d.lookup("r"));
            rMax = max(rTab);
        }
        else if (t == "polarOut")
        {
            // solid OUTSIDE the star-shaped outline r(theta), up to rOut (internal gears)
            kind = POLAROUT;
            rTab = scalarField(d.lookup("r"));
            rMax = rOut;
            rIn = min(rTab);      // nothing solid inside the tooth-tip circle
        }
        else if (t == "annulus") { kind = ANNULUS; rMax = rOut; }
        else if (t == "hole") { kind = HOLE; rMax = rOut; }
        else FatalErrorInFunction << "unknown primitive " << t << exit(FatalError);
    }

    //- point in body frame (x, y relative to the body centre), world z
    bool inside(scalar x, scalar y, scalar z) const
    {
        if (z <= z0 || z >= z1) return false;
        const scalar dx = x - c.x(), dy = y - c.y();
        const scalar r2 = dx*dx + dy*dy;
        if (r2 >= rMax*rMax) return false;
        if (kind == ANNULUS || kind == HOLE) return r2 > rIn*rIn;
        if (kind == POLAROUT && r2 <= rIn*rIn) return false;
        const label N = rTab.size();
        scalar th = ::atan2(dy, dx);
        if (th < 0) th += constant::mathematical::twoPi;
        const scalar s = th/constant::mathematical::twoPi*N;
        label i = label(s);
        const scalar f = s - i;
        i = i % N;
        const scalar rb = (1 - f)*rTab[i] + f*rTab[(i + 1) % N];
        return (kind == POLAR) ? (r2 < rb*rb) : (r2 > rb*rb);
    }
};


//- piecewise-linear omega(t) table with exact integral (angle)
class ibTable
{
public:
    scalarField t_, w_;
    ibTable() {}
    explicit ibTable(const List<Tuple2<scalar, scalar>>& tab)
    :
        t_(tab.size()), w_(tab.size())
    {
        forAll(tab, i) { t_[i] = tab[i].first(); w_[i] = tab[i].second(); }
    }
    scalar value(scalar t) const
    {
        if (t <= t_[0]) return w_[0];
        for (label i = 1; i < t_.size(); ++i)
        {
            if (t <= t_[i])
            {
                const scalar f = (t - t_[i-1])/max(t_[i] - t_[i-1], VSMALL);
                return w_[i-1] + f*(w_[i] - w_[i-1]);
            }
        }
        return w_.last();
    }
    //- integral from t_[0] (=0) to t
    scalar integral(scalar t) const
    {
        scalar s = 0;
        for (label i = 1; i < t_.size(); ++i)
        {
            const scalar a = t_[i-1], b = min(t_[i], t);
            if (b <= a) break;
            s += 0.5*(value(a) + value(b))*(b - a);
        }
        if (t > t_.last()) s += w_.last()*(t - t_.last());
        return s;
    }
};


class ibBody
{
public:
    word name;
    vector2D c0;
    ibTable orbitT, spinT;
    PtrList<ibPrimitive> prims;
    scalar zMin, zMax, reach;
    scalar rFree;    // radius around the body centre that is certainly fluid (rings)

    ibBody(const word& n, const dictionary& d)
    :
        name(n),
        c0(d.get<vector2D>("center")),
        orbitT(List<Tuple2<scalar, scalar>>(d.lookup("orbit"))),
        spinT(List<Tuple2<scalar, scalar>>(d.lookup("spin")))
    {
        const dictionary& pd = d.subDict("primitives");
        zMin = GREAT; zMax = -GREAT; reach = 0; rFree = GREAT;
        for (const entry& e : pd)
        {
            prims.append(new ibPrimitive(e.dict()));
            const ibPrimitive& p = prims.last();
            zMin = min(zMin, p.z0); zMax = max(zMax, p.z1);
            reach = max(reach, mag(p.c) + p.rMax);
            const bool centred = mag(p.c) < SMALL
                && (p.kind == ibPrimitive::POLAROUT || p.kind == ibPrimitive::ANNULUS);
            rFree = centred ? min(rFree, p.rIn) : 0;
        }
    }
};


class immersedBodies
{
    const fvMesh& mesh_;
    IOdictionary dict_;
    PtrList<ibBody> bodies_;
    scalar penalty_;
    bool displaceOil_;
    label nSub_;
    bool debug_;
    vectorField halfExt_;   // half extents of each cell bounding box

public:

    immersedBodies(const fvMesh& mesh)
    :
        mesh_(mesh),
        dict_
        (
            IOobject("immersedBodiesDict", mesh.time().constant(), mesh,
                     IOobject::MUST_READ, IOobject::NO_WRITE)
        ),
        penalty_(dict_.getOrDefault<scalar>("penalty", 1000)),
        displaceOil_(dict_.getOrDefault<bool>("displaceOil", true)),
        nSub_(dict_.getOrDefault<label>("subSamples", 2)),
        debug_(dict_.getOrDefault<bool>("debug", false)),
        halfExt_(mesh.nCells(), Zero)
    {
        const dictionary& bd = dict_.subDict("bodies");
        for (const entry& e : bd)
        {
            bodies_.append(new ibBody(e.keyword(), e.dict()));
        }
        const pointField& pts = mesh.points();
        const labelListList& cp = mesh.cellPoints();
        forAll(cp, ci)
        {
            boundBox bb(pts, cp[ci], false);
            halfExt_[ci] = 0.5*bb.span();
        }
        Info<< "immersedBodies: " << bodies_.size() << " bodies, penalty "
            << penalty_ << ", oil displacement " << displaceOil_ << endl;
    }

    scalar penalty() const { return penalty_; }

    //- Pose of a body at one instant (trigonometry evaluated once per step)
    struct pose
    {
        scalar cp, sp, cs, ss;
        vector2D c;
    };

    static pose makePose(const ibBody& b, scalar phi, scalar psi)
    {
        pose q;
        q.cp = ::cos(phi); q.sp = ::sin(phi);
        q.cs = ::cos(psi); q.ss = ::sin(psi);
        q.c = vector2D(q.cp*b.c0.x() - q.sp*b.c0.y(), q.sp*b.c0.x() + q.cp*b.c0.y());
        return q;
    }

    //- true if world point p is inside body b in pose q
    static bool insideBody(const ibBody& b, const point& p, const pose& q)
    {
        if (p.z() <= b.zMin || p.z() >= b.zMax) return false;
        // orbit frame coordinates
        const scalar qx =  q.cp*p.x() + q.sp*p.y();
        const scalar qy = -q.sp*p.x() + q.cp*p.y();
        const scalar dx = qx - b.c0.x(), dy = qy - b.c0.y();
        if (dx*dx + dy*dy > sqr(b.reach)) return false;
        const scalar lx =  q.cs*dx + q.ss*dy;
        const scalar ly = -q.ss*dx + q.cs*dy;
        bool solid = false;
        forAll(b.prims, k)
        {
            const ibPrimitive& pr = b.prims[k];
            if (pr.kind == ibPrimitive::HOLE)
            {
                if (pr.inside(lx, ly, p.z())) return false;
            }
            else if (!solid && pr.inside(lx, ly, p.z()))
            {
                solid = true;
            }
        }
        return solid;
    }

    void update(scalar t, volScalarField& chi, volVectorField& Us, volScalarField& alpha1)
    {
        const scalarField chiOld(chi.primitiveField());
        const vectorField& C = mesh_.C();
        scalarField& chiI = chi.primitiveFieldRef();
        vectorField& UsI = Us.primitiveFieldRef();
        chiI = 0;
        UsI = Zero;

        List<pose> Q(bodies_.size());
        List<scalar> wo(bodies_.size()), ws(bodies_.size());
        if (debug_) Info<< "ib: angles" << endl;
        forAll(bodies_, bi)
        {
            Q[bi] = makePose(bodies_[bi], bodies_[bi].orbitT.integral(t), bodies_[bi].spinT.integral(t));
            wo[bi] = bodies_[bi].orbitT.value(t);
            ws[bi] = bodies_[bi].spinT.value(t);
        }
        const label n = nSub_;
        const scalar w = 1.0/(n*n*n);
        forAll(C, ci)
        {
            const point& cc = C[ci];
            const vector& h = halfExt_[ci];
            const scalar hxy = ::sqrt(sqr(h.x()) + sqr(h.y()));
            forAll(bodies_, bi)
            {
                const ibBody& b = bodies_[bi];
                if (cc.z() + h.z() <= b.zMin || cc.z() - h.z() >= b.zMax) continue;
                // cheap rejection: cell farther from the body's CURRENT centre than its reach
                const pose& q = Q[bi];
                const scalar d2 = sqr(cc.x() - q.c.x()) + sqr(cc.y() - q.c.y());
                if (d2 > sqr(b.reach + hxy)) continue;
                if (b.rFree > hxy && d2 < sqr(b.rFree - hxy)) continue;
                scalar frac = 0;
                for (label i = 0; i < n; ++i)
                for (label j = 0; j < n; ++j)
                for (label k = 0; k < n; ++k)
                {
                    const point ps
                    (
                        cc.x() + h.x()*((2*i + 1.0)/n - 1),
                        cc.y() + h.y()*((2*j + 1.0)/n - 1),
                        cc.z() + h.z()*((2*k + 1.0)/n - 1)
                    );
                    if (insideBody(b, ps, q)) frac += w;
                }
                if (frac > chiI[ci])
                {
                    chiI[ci] = frac;
                    // v = w_orbit z x p + w_spin z x (p - c(t))
                    UsI[ci] = vector
                    (
                        -wo[bi]*cc.y() - ws[bi]*(cc.y() - q.c.y()),
                         wo[bi]*cc.x() + ws[bi]*(cc.x() - q.c.x()),
                         0
                    );
                }
            }
        }
        if (debug_) Info<< "ib: chi done" << endl;
        chi.correctBoundaryConditions();
        Us.correctBoundaryConditions();
        if (debug_) Info<< "ib: bc done" << endl;

        if (!displaceOil_) return;

        // conservative push-out of oil from cells newly covered by bodies
        scalarField& a = alpha1.primitiveFieldRef();
        const scalarField& V = mesh_.V();
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        scalarField excess(mesh_.nCells(), 0);
        forAll(a, ci)
        {
            const scalar dchi = chiI[ci] - chiOld[ci];
            if (dchi > 0 && a[ci] > 0)
            {
                const scalar fluidBefore = max(1 - chiOld[ci], SMALL);
                const scalar moved = min(a[ci], a[ci]*dchi/fluidBefore)*V[ci];
                excess[ci] = moved;
                a[ci] -= moved/V[ci];
            }
        }
        // capacity of receivers: empty fluid volume
        scalarField cap(mesh_.nCells(), 0);
        forAll(a, ci) { cap[ci] = max(0.0, (1 - chiI[ci]) - a[ci])*V[ci]; }
        scalarField sumCap(mesh_.nCells(), 0);
        forAll(own, f)
        {
            if (excess[own[f]] > 0) sumCap[own[f]] += cap[nei[f]];
            if (excess[nei[f]] > 0) sumCap[nei[f]] += cap[own[f]];
        }
        // proposed transfers per face (donor -> receiver), then scaled so that
        // no receiver gets more than its capacity; unsent oil stays put
        scalarField mON(own.size(), 0), mNO(own.size(), 0), demand(mesh_.nCells(), 0);
        forAll(own, f)
        {
            const label o = own[f], nb = nei[f];
            if (excess[o] > 0 && sumCap[o] > VSMALL)
            {
                mON[f] = min(excess[o], sumCap[o])*cap[nb]/sumCap[o];
                demand[nb] += mON[f];
            }
            if (excess[nb] > 0 && sumCap[nb] > VSMALL)
            {
                mNO[f] = min(excess[nb], sumCap[nb])*cap[o]/sumCap[nb];
                demand[o] += mNO[f];
            }
        }
        scalarField sc(mesh_.nCells(), 1.0), gain(mesh_.nCells(), 0), sent(mesh_.nCells(), 0);
        forAll(sc, ci) { if (demand[ci] > cap[ci]) sc[ci] = cap[ci]/max(demand[ci], VSMALL); }
        forAll(own, f)
        {
            const label o = own[f], nb = nei[f];
            const scalar a1 = mON[f]*sc[nb], a2 = mNO[f]*sc[o];
            gain[nb] += a1; sent[o] += a1;
            gain[o] += a2;  sent[nb] += a2;
        }
        forAll(a, ci)
        {
            a[ci] += (gain[ci] + max(0.0, excess[ci] - sent[ci]))/V[ci];
        }
        alpha1.correctBoundaryConditions();
    }
};

} // End namespace Foam

#endif
