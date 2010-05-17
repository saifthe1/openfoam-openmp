
if((PRINTVECTOR & PRINT_rhoEqn) > 0)
	TS_TOGGLE(true);
else
	TS_TOGGLE(false);

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

//#ifdef PARALLELIZE
//int chunk = Y.size() / omp_get_num_threads();
//#pragma omp parallel for default(shared) private(i) schedule(static, chunk)
//#endif
TS_START("rhoEqn.H > for");
for (label i=0; i<Y.size(); i++)
{
    if (dieselSprayArray[i])
    {
        Sevap += dieselSpray.evaporationSource(i);
    }
}
TS_END("rhoEqn.H > for");

{
	TS_START("rhoEqn.H > solve");
/*
	#ifdef PARALLELIZE
	TS_START("rhoEqn.H > 1");
    fvm::ddt(rho);
	TS_END("rhoEqn.H > 1");

	TS_START("rhoEqn.H > 2");
    fvc::div(phi);
	TS_END("rhoEqn.H > 2");

	TS_START("rhoEqn.H > 3");
    fvm::ddt(rho) + fvc::div(phi);
	TS_END("rhoEqn.H > 3");
	#else
*/
	TS_START("rhoEqn.H > solve()");
    solve
    (
        fvm::ddt(rho)
      + fvc::div(phi)
      ==
        Sevap
    );
	TS_END("rhoEqn.H > solve()");
//	#endif
	TS_END("rhoEqn.H > solve");
}

TS_TOGGLE(false);

