
if((PRINTVECTOR & PRINT_hEqn) > 1) 
	TS_TOGGLE(true);
else
	TS_TOGGLE(false);

{
	TS_START("hEqn.H > solve");
	#ifdef PARALLELIZE
	fvScalarMatrix *left1 = NULL,
				   *left2 = NULL,
				   *left3 = NULL;
	volScalarField *right = NULL;

	#pragma omp parallel sections
	{
		#pragma omp section
		{
			left1 = new fvScalarMatrix(fvm::ddt(rho, h)());
		}
		#pragma omp section
		{
			left2 = new fvScalarMatrix(mvConvection->fvmDiv(phi, h));
			left3 = new fvScalarMatrix(fvm::laplacian(turbulence->alphaEff(), h));
			right = new volScalarField(DpDt + dieselSpray.heatTransferSource());
		}
	}

	solve( *left1 + *left2 -(*left3) == *right );
	#else
    solve
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
       DpDt
     + dieselSpray.heatTransferSource()
    );
	#endif
	TS_END("hEqn.H > solve");

	TS_START("hEqn.H > thermo.correct()");
    thermo.correct();
	TS_END("hEqn.H > thermo.correct()");
}

TS_TOGGLE(false);
