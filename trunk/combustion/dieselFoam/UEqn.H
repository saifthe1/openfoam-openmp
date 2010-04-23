
#if PRINTVECTOR & 1000010 > 0
TS_TOGGLE(true);
#else
TS_TOGGLE(false);
#endif

#ifdef PARALLELIZE
fvVectorMatrix *t11 = NULL, *t12 = NULL;
volVectorField *t2 = NULL;

#pragma omp sections
{
    #pragma omp section
    {
        t11 = new fvVectorMatrix(fvm::ddt(rho, U) + fvm::div(phi, U));     // 0.78646 s
        t2 = new volVectorField((rho*g + dieselSpray.momentumSource())()); // 0.10644 s
    }
    #pragma omp section
        t12 = new fvVectorMatrix(turbulence->divDevRhoReff(U));            // 1.2071  s
}
fvVectorMatrix UEqn(*t11 + *t12 == *t2);                                   // 0.13733 s

// Totalt 2.2374 s
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
