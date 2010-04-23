
#if PRINTVECTOR & 1010000 > 0
TS_TOGGLE(true);
#else
TS_TOGGLE(false);
#endif

{
	TS_START("hEqn.H > solve");
	#ifdef PARALLELIZE
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
