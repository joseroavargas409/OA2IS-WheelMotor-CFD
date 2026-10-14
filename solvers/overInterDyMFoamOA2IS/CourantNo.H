/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

if (mesh.nInternalFaces())
{
    surfaceScalarField phiMask(localMin<scalar>(mesh).interpolate(cellMask));

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phiMask*phi))().internalField()
    );

    CoNum = 0.5*gMax
        (
            // OA2IS: interpolated (acceptor) cells are not solved; their face
            // fluxes are taken from donors and must not limit the time step.
            // Solved neighbours still see those faces in their own sums.
            sumPhi/mesh.V().field()*cellMask.primitiveField()*interpolatedCells.primitiveField()
        )*runTime.deltaTValue();

    if (getenv("OA2IS_DEBUG_CO"))
    {
        const scalarField coC
        (
            0.5*sumPhi/mesh.V().field()*cellMask.primitiveField()
           *interpolatedCells.primitiveField()*runTime.deltaTValue()
        );
        scalarField coL(coC);
        forAll(coL, i) { if (alpha1[i] <= 1e-3) coL[i] = 0; }
        const label im = findMax(coL);
        const labelIOList& zid = mesh.lookupObject<labelIOList>("zoneID");
        Pout<< "COMAX-LIQ " << coC[im] << " cell " << im << " at " << mesh.C()[im]
            << " V " << mesh.V()[im] << " zone " << zid[im] << " U " << U[im]
            << " alpha " << alpha1[im] << endl;
    }

    // OA2IS: separate Courant limits for liquid-bearing cells and for pure gas.
    // CoNum (checked against maxCo by setDeltaT.H) is the max over cells with
    // oil (alpha > 1e-3) and over gas cells scaled by maxCo/maxAirCo, i.e. the
    // gas only has to satisfy Co < maxAirCo (controlDict, default = maxCo).
    {
        const scalar maxCoD(runTime.controlDict().get<scalar>("maxCo"));
        const scalar maxAirCo(runTime.controlDict().getOrDefault<scalar>("maxAirCo", maxCoD));
        const scalarField coC
        (
            0.5*sumPhi/mesh.V().field()*cellMask.primitiveField()
           *interpolatedCells.primitiveField()*runTime.deltaTValue()
        );
        const scalarField& a1 = alpha1.primitiveField();
        scalar coLiq = 0, coGas = 0;
        forAll(coC, i)
        {
            if (a1[i] > 1e-3) coLiq = max(coLiq, coC[i]);
            else coGas = max(coGas, coC[i]);
        }
        reduce(coLiq, maxOp<scalar>());
        reduce(coGas, maxOp<scalar>());
        CoNum = max(coLiq, coGas*maxCoD/maxAirCo);
        Info<< "Courant Number liquid max: " << coLiq << " gas max: " << coGas << endl;
    }

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;

// ************************************************************************* //
