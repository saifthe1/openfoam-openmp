
if((PRINTVECTOR & PRINT_UEqn) > 0)
	TS_TOGGLE(true);
else
	TS_TOGGLE(false);

TS_START("UEqn.H > UEqn");
#ifdef PARALLELIZE
fvVectorMatrix *left1 = NULL,
			   *left2 = NULL;
volVectorField *right = NULL;

#pragma omp parallel sections default(shared)
{
    #pragma omp section
    {
		TS_START("UEqn.H > UEqn > left1 + right");
		left1 = new fvVectorMatrix(fvm::ddt(rho, U) + fvm::div(phi, U));
		right = new volVectorField((rho*g + dieselSpray.momentumSource())());
		TS_END("UEqn.H > UEqn > left1 + right");
    }
    #pragma omp section
	{
		TS_START("UEqn.H > UEqn > left2");
		left2 = new fvVectorMatrix(turbulence->divDevRhoReff(U));
		TS_END("UEqn.H > UEqn > left2");
	}
}
TS_START("UEqn.H > UEqn()");
fvVectorMatrix UEqn(*left1 + *left2 == *right);
TS_END("UEqn.H > UEqn()");

delete left1;
delete left2;
delete right;
#else

TS_START("UEqn.H > UEqn()");
fvVectorMatrix UEqn
(
	fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
 ==
	rho*g
  + dieselSpray.momentumSource()
);
TS_END("UEqn.H > UEqn()");
#endif
TS_END("UEqn.H > UEqn");
if (momentumPredictor)
{
	TS_START("UEqn.H > solve()");
	solve(UEqn == -fvc::grad(p));
	TS_END("UEqn.H > solve()");
}
TS_TOGGLE(false);

