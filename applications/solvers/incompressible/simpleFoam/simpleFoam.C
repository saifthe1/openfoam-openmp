/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow

\*---------------------------------------------------------------------------*/

#include <omp.h>

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	// * * MOVED OUT OF LOOP * * //
	dictionary simple = mesh.solutionDict().subDict("SIMPLE");

	int nNonOrthCorr =
		simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
	// * * * * * * * * * * * * * //

	Info<< "\nStarting time loop\n" << endl;
	const int startTime = runTime.startTime().value(),
			  endTime = runTime.endTime().value(),
			  deltaT = runTime.deltaT().value();
	#pragma omp parallel for default(shared) ordered schedule(static, 1)
	for(int bigLoopI = startTime; bigLoopI < endTime; bigLoopI += deltaT)
    {
		#pragma omp ordered
		{
			runTime++;
		}

        Info<< "Time = " << runTime.timeName() << nl << endl;

        /* * Replaced 
		 * #include "readSIMPLEControls.H"
         * #include "initConvergenceCheck.H"
		 */
		// readSIMPLEControls
		// initConvergenceCheck
		scalar eqnResidual = 1, maxResidual = 0;
		scalar convergenceCriterion = 0;

		simple.readIfPresent("convergence", convergenceCriterion);
        /* * End of replacement * */

        p.storePrevIter();

        // Pressure-velocity SIMPLE corrector
        {
			/* * Replaced
           	 * #include "UEqn.H"
             * #include "pEqn.H"
			 */
			// UEqn:
			tmp<fvVectorMatrix> *UEqn = NULL;
			
			UEqn = new tmp<fvVectorMatrix> 
			(
				fvm::div(phi, U)
			  + turbulence->divDevReff(U)
			);

			(*UEqn)().relax();

			eqnResidual = solve
			(
				(*UEqn)() == -fvc::grad(p)
			).initialResidual();

			maxResidual = max(eqnResidual, maxResidual);

			// pEqn:
			p.boundaryField().updateCoeffs();

			volScalarField AU = (*UEqn)().A();
			U = (*UEqn)().H()/AU;
			(*UEqn).clear();

			phi = fvc::interpolate(U) & mesh.Sf();
			adjustPhi(phi, U, p);


			for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
			{
				fvScalarMatrix *pEqn = NULL;
				pEqn = new fvScalarMatrix 
				(
					fvm::laplacian(1.0/AU, p) == fvc::div(phi)
				);
				
				pEqn->setReference(pRefCell, pRefValue);

				if (nonOrth == 0)
				{
					{
						eqnResidual = pEqn->solve().initialResidual();
						maxResidual = max(eqnResidual, maxResidual);
					}
				}
				else
				{
					pEqn->solve();
				}

				if (nonOrth == nNonOrthCorr)
				{
					phi -= pEqn->flux();
				}
			}

			p.relax();

			U -= fvc::grad(p)/AU;
			U.correctBoundaryConditions();
        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        /* * Replaced
		 * #include "convergenceCheck.H"
		 */
		if (maxResidual < convergenceCriterion)
		{
			Info<< "reached convergence criterion: " << convergenceCriterion << endl;
			// If this happens, we're screwed. =)
			runTime.writeAndEnd();
			Info<< "latestTime = " << runTime.timeName() << endl;
		}
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
