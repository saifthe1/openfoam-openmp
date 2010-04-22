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
#include <omp.h>
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

    Info << "\nStarting time loop\n" << endl;
	
  //  int yttreforlooptot = 0;
  //  int inreforlooptot = 0;
#pragma omp parallel
{
    while (runTime.run())
    {
//skapa variablerna som sätts i critical-sektionen
   dictionary *piso;
    int *nCorr;
    int *nNonOrthCorr;
    bool *momentumPredictor;
    bool *transonic;
    int *nOuterCorr;
    
//pEqn.H anropar UEqn.A()
    fvVectorMatrix *UEqn;
// I rhoEqn.H
	volScalarField *Sevap;

#pragma omp single
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
        #include "UEqn.H"
}//pragma omp critical slut

//int yttreforloop = 0;
//int inreforloop = 0;

#pragma omp for
        for (label ocorr=1; ocorr <= *nOuterCorr; ocorr++)
        {
//Info<< nl << "ANTAL TRÅDAR I OPENMP " << omp_get_num_threads() << nl <<endl;
//yttreforloop++;
//yttreforlooptot++;
            #include "YEqn.H"
            #include "hEqn.H"

            // --- PISO loop
//#pragma omp for
            for (int corr=1; corr<=*nCorr; corr++)
            {
//inreforlooptot++;
//inreforloop++;
                #include "pEqn.H"
            }
        }
		
        turbulence->correct();

        #include "spraySummary.H"

        rho = thermo.rho();

        if (runTime.write())
        {
            chemistry.dQ()().write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
//	    << nl
//	    << "yttre for-loopen körs " << yttreforlooptot << " gånger"
//	    << nl
//	    << "inre for-loopen körs  " << inreforlooptot << " gånger"
            << nl << endl;
    }
}//pragma omp parallel slut

    Info<< "End\n" << endl;
std:: cout << omp_get_max_threads() << std:: endl;
    return 0;
}


// ************************************************************************* //
