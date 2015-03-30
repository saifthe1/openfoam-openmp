# DieselFoam #

```
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

    while (runTime.run())
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

        for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
        {
            #include "YEqn.H"
            #include "hEqn.H"

            // --- PISO loop
            for (int corr=1; corr<=nCorr; corr++)
            {
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
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
```

```
/**###################################################################
/*# förslag för parallellisering av dieselFoam.C                     #
/*####################################################################*/

#pragma omp parallel
{
    while( runTime.run() )
    {
        #pragma omp critical 
        {
            //Evolving spray och annat 
        } 

        //yttre for-loopen kör 1 gång och inre 2 ggr för varje iteration
        //av while-loopen 
        #pragma omp for 
        for () //yttre for-loopen
        {
            #include "YEqn.H"
            #include "hEqn.H" 

            #pragma omp for //inte bestämt om den här ska användas 
            for () //inre for-loopen
            {
                #include "pEqn.H" 
            } 
        } 
    } 
}

################################################
## så här ser den nya koden ut                ##
## parallelliseringen ligger just nu på YEqn.H##
## istället för på dieselFoam.C               ##
## vi får tyvärr att exekveringstiden ökar    ##
## när man ökar antalet trådar                ##
## vi har kört denna kod på en 4-kärnig dator ##
################################################
#include <omp.h>
//      #define PARALLELIZE 1

#include "/home/zut/OpenFOAM/timestamp.hpp"

        #define PRINT_dieselFoam        0x01    // 0000 0001
        #define PRINT_hEqn                      0x02    // 0000 0010
        #define PRINT_pEqn                      0x04    // 0000 0100
        #define PRINT_rhoEqn            0x08    // 0000 1000
        #define PRINT_UEqn                      0x10    // 0001 0000
        #define PRINT_YEqn                      0x20    // 0010 0000
        #define PRINT_ALL                       0xFF    // 1111 1111

        const unsigned int PRINTVECTOR = PRINT_YEqn;


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
        omp_set_num_threads(4);
        
        if((PRINTVECTOR & PRINT_dieselFoam) > 1)
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

        Info << "\tSolving rho:" << endl;
        #include "rhoEqn.H"
        Info << "\tSolving U:" << endl;
        #include "UEqn.H"

        Info << "\tEntering loop:" << endl;
        for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
        {
                Info << "\t\tSolving Y:" << endl;
            #include "YEqn.H"
                Info << "\t\tSolving h:" << endl;
            #include "hEqn.H"

            // --- PISO loop
                Info << "\t\tEntering inner loop:" << endl;
            for (int corr=1; corr<=nCorr; corr++)
            {
                        Info << "\t\t\tSolving p:" << endl;
                #include "pEqn.H"
            }
        }

        Info << "\tOutside of loops, doing other stuff:" << endl;
        turbulence->correct();

        #include "spraySummary.H"

        rho = thermo.rho();

        if (runTime.write())
        {
            chemistry.dQ()().write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
        
        TS_TOGGLE(true);
        TS_PRINT();
        
    return 0;
}

```