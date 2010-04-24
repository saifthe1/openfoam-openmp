
if((PRINTVECTOR & PRINT_UEqn) > 1)
	TS_TOGGLE(true);
else
	TS_TOGGLE(false);


#ifdef PARALLELIZE
fvVectorMatrix *t11 = NULL, *t12 = NULL;
volVectorField *t2 = NULL;

#pragma omp parallel sections default(shared)
{
    #pragma omp section
    {
		t11 = new fvVectorMatrix(fvm::ddt(rho, U) + fvm::div(phi, U));
		t2 = new volVectorField((rho*g + dieselSpray.momentumSource())());
    }
    #pragma omp section
	{
		t12 = new fvVectorMatrix(turbulence->divDevRhoReff(U));
	}
}
fvVectorMatrix UEqn(*t11 + *t12 == *t2);
#else


fvVectorMatrix UEqn
(
	fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
 ==
	rho*g
  + dieselSpray.momentumSource()
);
#endif
if (momentumPredictor)
{
	solve(UEqn == -fvc::grad(p));
}

TS_TOGGLE(false);
