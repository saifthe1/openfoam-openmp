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
#include "/home/zut/openfoam-openmp/ordered/order.hpp"

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
	omp_set_num_threads(10);
	
	// * * MOVED OUT OF LOOP * * //
	dictionary simple = mesh.solutionDict().subDict("SIMPLE");

	int nNonOrthCorr =
		simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
	// * * * * * * * * * * * * * //

	Info<< "\nStarting time loop\n" << endl;
	const unsigned int  startTime	= runTime.startTime().value(),
						endTime		= runTime.endTime().value(),
						deltaT		= runTime.deltaT().value();
	unsigned int		bigLoopI;
						
	#pragma omp parallel for default(shared) private(bigLoopI) ordered schedule(static, 1)
	for(bigLoopI = startTime; bigLoopI < endTime; bigLoopI += deltaT)
    {
		unsigned int sect = 0;

		ORDER_START(sect, bigLoopI);
		runTime++;
		ORDER_END(sect++);
	
		ORDER_START(sect, bigLoopI);
        Info<< "Time = " << runTime.timeName() << nl << endl;
		ORDER_END(sect++);

        /* * Replaced 
		 * #include "readSIMPLEControls.H"
         * #include "initConvergenceCheck.H"
		 */
		// readSIMPLEControls
		// initConvergenceCheck
		scalar eqnResidual = 1, maxResidual = 0;
		scalar convergenceCriterion = 0;

		ORDER_START(sect, bigLoopI);
		simple.readIfPresent("convergence", convergenceCriterion);
		ORDER_END(sect++);
	
        /* * End of replacement * */

		ORDER_START(sect, bigLoopI);
        p.storePrevIter();
		ORDER_END(sect++);
	

        // Pressure-velocity SIMPLE corrector
        {
			/* * Replaced
           	 * #include "UEqn.H"
             * #include "pEqn.H"
			 */
			// UEqn:
			tmp<fvVectorMatrix> *UEqn = NULL;
			
			ORDER_START(sect, bigLoopI);
			UEqn = new tmp<fvVectorMatrix> 
			(
				fvm::div(phi, U)
			  + turbulence->divDevReff(U)
			);
			ORDER_END(sect++);
		

			ORDER_START(sect, bigLoopI);
			(*UEqn)().relax();
			ORDER_END(sect++);
		

			ORDER_START(sect, bigLoopI);
			eqnResidual = solve
			(
				(*UEqn)() == -fvc::grad(p)
			).initialResidual();
			ORDER_END(sect++);
		

			ORDER_START(sect, bigLoopI);
			maxResidual = max(eqnResidual, maxResidual);
			ORDER_END(sect++);
		

			// pEqn:
			ORDER_START(sect, bigLoopI);
			p.boundaryField().updateCoeffs();
			ORDER_END(sect++);
		

			volScalarField AU = (*UEqn)().A();
			U = (*UEqn)().H()/AU;
			(*UEqn).clear();

			ORDER_START(sect, bigLoopI);
			phi = fvc::interpolate(U) & mesh.Sf();
			ORDER_END(sect++);
		
			ORDER_START(sect, bigLoopI);
			adjustPhi(phi, U, p);
			ORDER_END(sect++);
		

			for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
			{
				fvScalarMatrix *pEqn = NULL;
				ORDER_START(sect, bigLoopI);
				pEqn = new fvScalarMatrix 
				(
					fvm::laplacian(1.0/AU, p) == fvc::div(phi)
				);
				ORDER_END(sect++);
			

				ORDER_START(sect, bigLoopI);
				pEqn->setReference(pRefCell, pRefValue);
				ORDER_END(sect++);
			

				if (nonOrth == 0)
				{
					ORDER_START(sect, bigLoopI);
					eqnResidual = pEqn->solve().initialResidual();
					maxResidual = max(eqnResidual, maxResidual);
					ORDER_END(sect++);
				
				}
				else
				{
					ORDER_START(sect, bigLoopI);
					pEqn->solve();
					ORDER_END(sect++);
				
				}

				if (nonOrth == nNonOrthCorr)
				{
					ORDER_START(sect, bigLoopI);
					phi -= pEqn->flux();
					ORDER_END(sect++);
				
				}
			}

			ORDER_START(sect, bigLoopI);
			p.relax();
			ORDER_END(sect++);
		

			ORDER_START(sect, bigLoopI);
			U -= fvc::grad(p)/AU;
			ORDER_END(sect++);
		
			ORDER_START(sect, bigLoopI);
			U.correctBoundaryConditions();
			ORDER_END(sect++);
		
        }

		ORDER_START(sect, bigLoopI);
        turbulence->correct();
		ORDER_END(sect++);
	

		ORDER_START(sect, bigLoopI);
        runTime.write();
		ORDER_END(sect++);
	

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
			ORDER_START(sect, bigLoopI);
			runTime.writeAndEnd();
			ORDER_END(sect++);
			
			Info<< "latestTime = " << runTime.timeName() << endl;
		}
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
