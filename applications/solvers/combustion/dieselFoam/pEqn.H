
if((PRINTVECTOR & PRINT_pEqn) > 0)
	TS_TOGGLE(true);
else
	TS_TOGGLE(false);

rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

TS_START("pEqn.H > if/else");
if(transonic)
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

        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
        }
    }
	TS_END("pEqn.H > if > for");
}
else
{
    phi =
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        );

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

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
	TS_END("pEqn.H > else > for");
}
TS_END("pEqn.H > if/else");

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
TS_TOGGLE(false);

