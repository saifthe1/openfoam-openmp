
if((PRINTVECTOR & PRINT_YEqn) > 0)
	TS_TOGGLE(true);
else
	TS_TOGGLE(false);

tmp<fv::convectionScheme<scalar> > mvConvection
(
	fv::convectionScheme<scalar>::New
	(
		mesh,
		fields,
		phi,
		mesh.divScheme("div(phi,Yi_h)")
	)
);
{

    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];

	/** Kan förmodligen parallelliseras
	 * const unsigned int chunk = Y.size() / omp_get_num_threads();
	 * #pragma omp for default(shared) private(i) schedule(static, chunk)
	 **/
	TS_START("YEqn.H > for");
    for (label i=0; i< Y.size(); i++)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];
			
			TS_START("YEqn.H > for > if > solve");
			#ifdef PARALLELIZE
			fvScalarMatrix *left1 = NULL,
						   *left2 = NULL,
						   *left3 = NULL;
			volScalarField *right = NULL;

			#pragma omp parallel sections default(shared)
			{
				#pragma omp section
				{
					TS_START("YEqn.H > for > if > left1 & right");
					left1 = new fvScalarMatrix(fvm::ddt(rho, Yi)());
					right = new volScalarField((dieselSpray.evaporationSource(i) + kappa*chemistry.RR(i))());
					TS_END("YEqn.H > for > if > left1 & right");
				}
				#pragma omp section
				{
					TS_START("YEqn.H > for > if > left2");
					left2 = new fvScalarMatrix(mvConvection->fvmDiv(phi, Yi));
					TS_END("YEqn.H > for > if > left2");
				}
				#pragma omp section
				{
					TS_START("YEqn.H > for > if > left3");
					left3 = new fvScalarMatrix(fvm::laplacian(turbulence->muEff(), Yi));
					TS_END("YEqn.H > for > if > left3");
				}
			}

			TS_START("YEqn.H > for > if > solve()");
			solve( *left1 + *left2 - *left3 == *right, mesh.solver("Yi") );
			TS_END("YEqn.H > for > if > solve()");
			
			delete left1;
			delete left2;
			delete left3;
			delete right;
			#else
			TS_START("YEqn.H > for > if > solve()");
            solve
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
              ==
                dieselSpray.evaporationSource(i)
              + kappa*chemistry.RR(i),
                mesh.solver("Yi")
            );
			TS_END("YEqn.H > for > if > solve()");
			#endif
			TS_END("YEqn.H > for > if > solve");
			
            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }
	TS_END("YEqn.H > for");
	Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
TS_TOGGLE(false);
