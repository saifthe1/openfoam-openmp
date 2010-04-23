
#if PRINTVECTOR & 1001000 > 0
TS_TOGGLE(true);
#else
TS_TOGGLE(false);
#endif

TS_START("pEqn.H > thermo.rho()");
rho = thermo.rho();
TS_END("pEqn.H > thermo.rho()");

TS_START("pEqn.H > 1.0/UEqn.A()");
volScalarField rUA = 1.0/UEqn.A();
TS_END("pEqn.H > 1.0/UEqn.A()");

TS_START("pEqn.H > rUA*UEqn.H()");
U = rUA*UEqn.H();
TS_END("pEqn.H > rUA*UEqn.H()");

TS_START("pEqn.H > if/else");


if(transonic)
{
	TS_START("pEqn.H > if > phi =");
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
    );
	TS_END("pEqn.H > if > phi =");

	TS_START("pEqn.H > if > for");
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
		TS_START("pEqn.H > if > for > pEqn");
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rUA, p)
         ==
            Sevap
        );
		TS_END("pEqn.H > if > for > pEqn");

		TS_START("pEqn.H > if > for > pEqn.solve()");
        pEqn.solve();
		TS_END("pEqn.H > if > for > pEqn.solve()");

		TS_START("pEqn.H > if > for > if");
        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
        }
		TS_END("pEqn.H > if > for > if");
    }
	TS_END("pEqn.H > if > for");
}
else
{
	TS_START("pEqn.H > else > phi =");
    phi =
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        );
	TS_END("pEqn.H > else > phi =");

	TS_START("pEqn.H > else > for");
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
		TS_START("pEqn.H > else > for > pEqn");
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, p)
         ==
            Sevap
        );
		TS_END("pEqn.H > else > for > pEqn");

		TS_START("pEqn.H > else > for > pEqn.solve()");
        pEqn.solve();
		TS_END("pEqn.H > else > for > pEqn.solve()");

		TS_START("pEqn.H > else > for > if");
        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
		TS_END("pEqn.H > else > for > if");
    }
	TS_END("pEqn.H > else > for");
}
TS_END("pEqn.H > if/else");

TS_START("pEqn.H > rhoEqn.H");
#include "rhoEqn.H"
// In case rho turned it off.
#if PRINTVECTOR & 1001000 > 0
TS_TOGGLE(true);
#else
TS_TOGGLE(false);
#endif

TS_END("pEqn.H > rhoEqn.H");

TS_START("pEqn.H > compressibleContinuityErrs.H");
#include "compressibleContinuityErrs.H"
TS_END("pEqn.H > compressibleContinuityErrs.H");

TS_START("pEqn.H > rUA*fvc::grad(p)");
U -= rUA*fvc::grad(p);
TS_END("pEqn.H > rUA*fvc::grad(p)");

TS_START("pEqn.H > U.correctBoundaryConditions()");
U.correctBoundaryConditions();
TS_END("pEqn.H > U.correctBoundaryConditions()");

TS_START("pEqn.H > fvc::DDt()");
DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
TS_END("pEqn.H > fvc::DDt()");
TS_TOGGLE(false);

