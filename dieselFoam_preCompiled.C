    while (runTime.run())
    {
# 1 "/home/openfoam/OpenFOAM/OpenFOAM-1.6/src/finiteVolume/lnInclude/readPISOControls.H" 1
    dictionary piso = mesh.solutionDict().subDict("PISO");

    int nCorr(readInt(piso.lookup("nCorrectors")));

    int nNonOrthCorr =
        piso.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    bool momentumPredictor =
        piso.lookupOrDefault<Switch>("momentumPredictor", true);

    bool transonic =
        piso.lookupOrDefault<Switch>("transonic", false);

    int nOuterCorr =
        piso.lookupOrDefault<int>("nOuterCorrectors", 1);
# 68 "dieselFoam.C" 2
# 1 "/home/openfoam/OpenFOAM/OpenFOAM-1.6/src/finiteVolume/lnInclude/compressibleCourantNo.H" 1
# 33 "/home/openfoam/OpenFOAM/OpenFOAM-1.6/src/finiteVolume/lnInclude/compressibleCourantNo.H"
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

if (mesh.nInternalFaces())
{
    surfaceScalarField SfUfbyDelta =
        mesh.surfaceInterpolation::deltaCoeffs()*mag(phi)/fvc::interpolate(rho);

    CoNum = max(SfUfbyDelta/mesh.magSf())
        .value()*runTime.deltaT().value();

    meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
        .value()*runTime.deltaT().value();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;
# 69 "dieselFoam.C" 2
# 1 "/home/openfoam/OpenFOAM/OpenFOAM-1.6/src/finiteVolume/lnInclude/setDeltaT.H" 1
# 35 "/home/openfoam/OpenFOAM/OpenFOAM-1.6/src/finiteVolume/lnInclude/setDeltaT.H"
if (adjustTimeStep)
{
    scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaT().value(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " << runTime.deltaT().value() << endl;
}
# 70 "dieselFoam.C" 2

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info << "Evolving Spray" << endl;

        dieselSpray.evolve();

        Info << "Solving chemistry" << endl;

        chemistry.solve
        (
            runTime.value() - runTime.deltaT().value(),
            runTime.deltaT().value()
        );


        {
            volScalarField tk =
                Cmix*sqrt(turbulence->muEff()/rho/turbulence->epsilon());
            volScalarField tc = chemistry.tc();


            kappa = (runTime.deltaT() + tc)/(runTime.deltaT()+tc+tk);
        }

# 1 "../dieselEngineFoam/rhoEqn.H" 1
# 33 "../dieselEngineFoam/rhoEqn.H"
volScalarField Sevap
(
    IOobject
    (
        "Sevap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.0)
);

for (label i=0; i<Y.size(); i++)
{
    if (dieselSpray.isLiquidFuel()[i])
    {
        Sevap += dieselSpray.evaporationSource(i);
    }
}

{
    solve
    (
        fvm::ddt(rho)
      + fvc::div(phi)
      ==
        Sevap
    );
}
# 97 "dieselFoam.C" 2
# 1 "../dieselEngineFoam/UEqn.H" 1
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho*g
      + dieselSpray.momentumSource()
    );

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }
# 98 "dieselFoam.C" 2

        for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
        {
# 1 "../dieselEngineFoam/YEqn.H" 1
tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{

    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];

    for (label i=0; i<Y.size(); i++)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            solve
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
              ==
                dieselSpray.evaporationSource(i)
              + kappa*chemistry.RR(i),
                mesh.solver("Yi")
            );

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
# 102 "dieselFoam.C" 2
# 1 "../dieselEngineFoam/hEqn.H" 1
{
    solve
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
       DpDt
     + dieselSpray.heatTransferSource()
    );

    thermo.correct();
}
# 103 "dieselFoam.C" 2


            for (int corr=1; corr<=nCorr; corr++)
            {
# 1 "pEqn.H" 1
rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

if (transonic)
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rUA, p)
         ==
            Sevap
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi =
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, p)
         ==
            Sevap
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
}

# 1 "../dieselEngineFoam/rhoEqn.H" 1
# 33 "../dieselEngineFoam/rhoEqn.H"
volScalarField Sevap
(
    IOobject
    (
        "Sevap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.0)
);

for (label i=0; i<Y.size(); i++)
{
    if (dieselSpray.isLiquidFuel()[i])
    {
        Sevap += dieselSpray.evaporationSource(i);
    }
}

{
    solve
    (
        fvm::ddt(rho)
      + fvc::div(phi)
      ==
        Sevap
    );
}
# 67 "pEqn.H" 2
# 1 "/home/openfoam/OpenFOAM/OpenFOAM-1.6/src/finiteVolume/lnInclude/compressibleContinuityErrs.H" 1
# 33 "/home/openfoam/OpenFOAM/OpenFOAM-1.6/src/finiteVolume/lnInclude/compressibleContinuityErrs.H"
{
    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr =
        (fvc::domainIntegrate(mag(rho - thermo.rho()))/totalMass).value();

    scalar globalContErr =
        (fvc::domainIntegrate(rho - thermo.rho())/totalMass).value();

    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr
        << endl;
}
# 68 "pEqn.H" 2

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
# 108 "dieselFoam.C" 2
            }
        }

        turbulence->correct();

# 1 "../dieselEngineFoam/spraySummary.H" 1
    label Nparcels = dieselSpray.size();
    reduce(Nparcels, sumOp<label>());

    Info<< "\nNumber of parcels in system.... | "
        << Nparcels << endl
        << "Injected liquid mass........... | "
        << 1e6*dieselSpray.injectedMass(runTime.value()) << " mg" << endl
        << "Liquid Mass in system.......... | "
        << 1e6*dieselSpray.liquidMass() << " mg" << endl
        << "SMD, Dmax...................... | "
        << dieselSpray.smd()*1e6 << " mu, "
        << dieselSpray.maxD()*1e6 << " mu"
        << endl;

    scalar evapMass =
        dieselSpray.injectedMass(runTime.value())
    - dieselSpray.liquidMass();

    scalar gasMass = fvc::domainIntegrate(rho).value();

    if (dieselSpray.twoD())
    {
        gasMass *= 2.0*mathematicalConstant::pi/dieselSpray.angleOfWedge();
    }

    scalar addedMass = gasMass - gasMass0;

    Info<< "Added gas mass................. | " << 1e6*addedMass << " mg"
        << nl << "Evaporation Continuity Error... | "
        << 1e6*(addedMass - evapMass) << " mg" << endl;
# 114 "dieselFoam.C" 2

        rho = thermo.rho();

        if (runTime.write())
        {
            chemistry.dQ()().write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
