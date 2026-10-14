// OA2IS: immersed moving bodies (see immersedBodies.H)
Info<< "Creating immersed-body fields chi, Us\n" << endl;
volScalarField chi
(
    IOobject("chi", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(dimless, Zero),
    zeroGradientFvPatchScalarField::typeName
);
volVectorField Us
(
    IOobject("Us", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector(dimVelocity, Zero),
    zeroGradientFvPatchVectorField::typeName
);
immersedBodies ib(mesh);
// initial pose (no oil displacement at t0: the initial field is set outside the bodies)
{
    volScalarField alphaTmp("alphaTmp", alpha1);
    ib.update(runTime.value(), chi, Us, alphaTmp);
}
