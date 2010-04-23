
#if PRINTVECTOR & 1000100 > 0
TS_TOGGLE(true);
#else
TS_TOGGLE(false);
#endif

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

List<bool> dieselSprayArray = dieselSpray.isLiquidFuel();

#ifdef PARALLELIZE
int chunk = Y.size() / omp_get_num_threads();
#pragma omp parallel for default(shared) private(i) schedule(static, chunk)
#endif
for (label i=0; i<Y.size(); i++)
{
    if (dieselSprayArray[i])
    {
        Sevap += dieselSpray.evaporationSource(i);
    }
}

{
	TS_START("rho > solve");
	#ifdef PARALLELIZE
	TS_START("rho > 1");
    fvm::ddt(rho);
	TS_END("rho > 1");

	TS_START("rho > 2");
    fvc::div(phi);
	TS_END("rho > 2");

	TS_START("rho > 3");
    fvm::ddt(rho) + fvc::div(phi);
	TS_END("rho > 3");
	#else
    solve
    (
        fvm::ddt(rho)
      + fvc::div(phi)
      ==
        Sevap
    );
	#endif
	TS_END("rho > solve");
}

TS_TOGGLE(false);

