
if((PRINTVECTOR & PRINT_hEqn) > 0) 
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
			TS_START("hEqn.H > solve > left1");
			left1 = new fvScalarMatrix(fvm::ddt(rho, h)());
			TS_END("hEqn.H > solve > left1");
		}
		#pragma omp section
		{
			TS_START("hEqn.H > solve > left2 + left3 + right");
			left2 = new fvScalarMatrix(mvConvection->fvmDiv(phi, h));
			left3 = new fvScalarMatrix(fvm::laplacian(turbulence->alphaEff(), h));
			right = new volScalarField(DpDt + dieselSpray.heatTransferSource());
			TS_END("hEqn.H > solve > left2 + left3 + right");
		}
	}

	TS_START("hEqn.H > solve()");
	solve( *left1 + *left2 - (*left3) == *right );
	TS_END("hEqn.H > solve()");
	#else
	TS_START("hEqn.H > solve()");
    solve
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
       DpDt
     + dieselSpray.heatTransferSource()
    );
	TS_END("hEqn.H > solve()");
	#endif
	TS_END("hEqn.H > solve");

    thermo.correct();
}

TS_TOGGLE(false);
