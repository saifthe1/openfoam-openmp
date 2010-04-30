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
#include "semaphore/semaphore.hpp"

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
	omp_set_num_threads(4);
	
	// * * MOVED OUT OF LOOP * * //
	dictionary simple = mesh.solutionDict().subDict("SIMPLE");

	int nNonOrthCorr =
		simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
	// * * * * * * * * * * * * * //

	Info<< "\nStarting time loop@" << runTime.startTime().value() << "\n" << endl;
	
	const unsigned int  startTime	= runTime.startTime().value(),
						endTime		= runTime.endTime().value(),
						deltaT		= runTime.deltaT().value();
	unsigned int		bigLoopI;
	int 				c			= -1;
				
	#pragma omp parallel for default(shared) private(bigLoopI) ordered schedule(static, 1)
	for(bigLoopI = startTime; bigLoopI < endTime; bigLoopI += deltaT)
    {
		// I am the one is set True for the first one into the loop.
		bool iato = false;
		#pragma omp critical 
		{
			if(c == -1) {
				iato = true;
				c++;
			}
		}
		if(!iato && c <= 2)
		{
			while(!iato && c <= 2) ;

		}

		tmp<fvVectorMatrix> *UEqn = NULL;
		fvScalarMatrix *pEqn = NULL;

		unsigned int sect = 0;
		// Sect#0
		runTime++;
	
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
			
			// Sect#3
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
			
			if(c <= 2)
			{
				ORDER_END(0);
			}

			for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
			{
				pEqn = new fvScalarMatrix 
				(
					fvm::laplacian(1.0/AU, p) == fvc::div(phi)
				);

				pEqn->setReference(pRefCell, pRefValue);

				if (nonOrth == 0)
				{
					eqnResidual = pEqn->solve().initialResidual();
					maxResidual = max(eqnResidual, maxResidual);
				}
				else
				{
					pEqn->solve();
				}

				if (nonOrth == nNonOrthCorr)
				{
					ORDER_START((bigLoopI-startTime)/deltaT, nonOrth);
					phi -= pEqn->flux();
					ORDER_START((bigLoopI-startTime)/deltaT, nonOrth);
				}
				delete pEqn;
			}

			ORDER_START(0, (bigLoopI-startTime)/deltaT);
			if(c <=2) c++;
			
			p.relax();
			U -= fvc::grad(p)/AU;
			U.correctBoundaryConditions();
			delete UEqn;
        }

		// Sect#5
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
