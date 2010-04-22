TS_TOGGLE(true);

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
	std::cerr << "\n\n Entering the loop..." << std::endl;
    for (label i=0; i< Y.size(); i++)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];
						
			fvScalarMatrix  *left1 = NULL,
							*left2 = NULL,
							*left3 = NULL;
			volScalarField *right = NULL;

			#pragma omp parallel sections
			{
				#pragma omp section
				{
					std::cerr << "\t\tEntering section 1..." << std::endl;
					TS_START("left1 & right");
					left1 = new fvScalarMatrix(fvm::ddt(rho, Yi)());
					right = new volScalarField((dieselSpray.evaporationSource(i) + kappa*chemistry.RR(i))());
					TS_END("left1 & right");
				}
				#pragma omp section
				{
					std::cerr << "\t\tEntering section 2..." << std::endl;
					TS_START("left2");
					left2 = new fvScalarMatrix(mvConvection->fvmDiv(phi, Yi));
					TS_END("left2");
				}
				#pragma omp section
				{
					std::cerr << "\t\tEntering section 3..." << std::endl;
					TS_START("left3");
					left3 = new fvScalarMatrix(fvm::laplacian(turbulence->muEff(), Yi));
					TS_END("left3");
				}
			}

			TS_START("solve()");
			solve( *left1 + *left2 - *left3 == *right, mesh.solver("Yi") );
			TS_END("solve()");
			
			/*
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
			*/

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }
	Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
TS_TOGGLE();
