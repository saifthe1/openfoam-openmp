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
    dieselFoam

Description
    Solver for diesel spray and combustion.

\*---------------------------------------------------------------------------*/

#include <iomanip>

#include <omp.h>
	#define PARALLELIZE 1

#include "/home/zut/OpenFOAM/timestamp.hpp"

	#define PRINT_dieselFoam	0x01	// 0000 0001
	#define PRINT_hEqn			0x02	// 0000 0010
	#define PRINT_pEqn			0x04	// 0000 0100
	#define PRINT_rhoEqn		0x08	// 0000 1000
	#define PRINT_UEqn			0x10	// 0001 0000
	#define PRINT_YEqn			0x20	// 0010 0000
	#define PRINT_ALL			0xFF	// 1111 1111

	const unsigned int PRINTVECTOR = PRINT_ALL & ~PRINT_pEqn & ~PRINT_rhoEqn;


#include "fvCFD.H"
#include "hCombustionThermo.H"
#include "turbulenceModel.H"
#include "spray.H"
#include "psiChemistryModel.H"
#include "chemistrySolver.H"

#include "multivariateScheme.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
	TS_TOGGLE(false);
//	omp_set_num_threads(4);
	
	if((PRINTVECTOR & PRINT_dieselFoam) > 0)
		TS_TOGGLE(true);
	else 
		TS_TOGGLE(false);
	
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readGravitationalAcceleration.H"
    #include "readCombustionProperties.H"
    #include "createSpray.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	double	startTime	= runTime.startTime().value(),
			endTime		= runTime.endTime().value(),
			deltaT		= runTime.deltaT().value();

    Info << "\nStarting time loop\n" << endl;

	while(runTime.run())
	{
        #include "readPISOControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info << "Evolving Spray" << endl;

        dieselSpray.evolve();

        Info << "Solving chemistry" << endl;

        chemistry.solve
        (
            runTime.value() - runTime.deltaT().value(),
            runTime.deltaT().value()
        );

        // turbulent time scale
        {
            volScalarField tk =
                Cmix*sqrt(turbulence->muEff()/rho/turbulence->epsilon());
            volScalarField tc = chemistry.tc();

            // Chalmers PaSR model
            kappa = (runTime.deltaT() + tc)/(runTime.deltaT()+tc+tk);
        }

        #include "rhoEqn.H"
		if((PRINTVECTOR & PRINT_dieselFoam) > 0) TS_TOGGLE(true); else TS_TOGGLE(false);
		TS_START("UEqn.H");
        #include "UEqn.H"
		if((PRINTVECTOR & PRINT_dieselFoam) > 0) TS_TOGGLE(true); else TS_TOGGLE(false);
		TS_END("UEqn.H");

        for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
        {
			TS_START("YEqn.H");
            #include "YEqn.H"
			if((PRINTVECTOR & PRINT_dieselFoam) > 0) TS_TOGGLE(true); else TS_TOGGLE(false);
			TS_END("YEqn.H");

			TS_START("hEqn.H");
            #include "hEqn.H"
			if((PRINTVECTOR & PRINT_dieselFoam) > 0) TS_TOGGLE(true); else TS_TOGGLE(false);
			TS_END("hEqn.H");

            // --- PISO loop
            for (int corr=1; corr<=nCorr; corr++)
            {
                #include "pEqn.H"
            }
        }

		std::cout << "turbulence->correct()" << std::endl;
        turbulence->correct();

        #include "spraySummary.H"

		std::cout << "thermo.rho()" << std::endl;
        rho = thermo.rho();

		std::cout << "runTime.write()" << std::endl;
        if (runTime.write())
        {
			std::cout << "chemistry.dQ()().write();" << std::endl;
            chemistry.dQ()().write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
	
	std::cout	<< "*** RESULTS ***" << std::endl
				<< std::setiosflags(std::ios::left)
				<< std::setw(13) << "startTime" << " = " << startTime << std::endl
				<< std::setw(13) << "endTime" << " = " << endTime << std::endl
				<< std::setw(13) << "deltaT" << " = " << deltaT << std::endl
				<< std::setw(13) << "threads"  << " = " << getenv("OMP_NUM_THREADS") << std::endl
				<< std::setw(13) << "meshsize" << " = " << mesh.cells().size() << std::endl 
				<< std::endl;

	TS_TOGGLE(true);
	TS_PRINT();
	
    return 0;
}


// ************************************************************************* //
