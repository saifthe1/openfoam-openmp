## rhoEqn.H ##
[code](http://code.google.com/p/openfoam-openmp/source/browse/trunk/rhoEqn.H)

Körs i roten av `while`-loopen, men anropas även av [pEqn.H](Equations#pEqn.H.md)

```
volScalarField Sevap
(
    IOobject
    (
        "Sevap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.0)
);

for (label i=0; i<Y.size(); i++)
{
    if (dieselSpray.isLiquidFuel()[i])
    {
        Sevap += dieselSpray.evaporationSource(i);
    }
}

{
    solve
    (
        fvm::ddt(rho)
      + fvc::div(phi)
      ==
        Sevap
    );
}
```

Det är bara solve som är muterande.

Parallellisering av `for`-loopen m.h.a. är pågående
```
#pragma omp parallel for      \  
  default(shared) private(i)  \  
  schedule(static,chunk)      \  
```


---


## UEqn.H ##
[code](http://code.google.com/p/openfoam-openmp/source/browse/trunk/UEqn.H)

Körs i roten av `while`-loopen

```
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
 ==
    rho*g
  + dieselSpray.momentumSource()
);

if (momentumPredictor)
{
    solve(UEqn == -fvc::grad(p));
}
```

Inga tilldelningar alls. Vi måste undersöka vidare vad de olika operatorerna och `solve()` gör, för att kunna förstå vart resultatet från ekvationen sparas.

Eftersom det andra stycket är beroende av det första, måste det hela utföras sekventiellt. Dock kan man med fördel dela upp beräkningarna som blir argument till `UEqn`'s konstruktor. Detta t.ex. med:
```
fvVectorMatrix *t11, *t12;
volVectorField *t2;

#pragma omp sections
{
    #pragma omp section
    {
        t11 = new fvVectorMatrix(fvm::ddt(rho, U) + fvm::div(phi, U));     // 0.78646 s
        t2 = new volVectorField((rho*g + dieselSpray.momentumSource())()); // 0.10644 s
    }
    #pragma omp section
        t12 = new fvVectorMatrix(turbulence->divDevRhoReff(U));            // 1.2071  s
}
fvVectorMatrix UEqn(*t11 + *t12 == *t2);                                   // 0.13733 s

// Totalt 2.2374 s
```

Stamps tagna under sekventiell körning

```
hashResult efter rensning av katalogerna (tidigare resultat)
med Pers parallellisering  d41d8cd98f00b204e9800998ecf8427e
med OpenFOAMs kod          d41d8cd98f00b204e9800998ecf8427e
```

---


## YEqn.H ##
[code](http://code.google.com/p/openfoam-openmp/source/browse/trunk/YEqn.H)

Körs i den yttre av `for`-looparna

```
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

    for (label i=0; i<Y.size(); i++)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

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
```

Förutom det trolleri som `solve()` åstakommer, ändrar den bara på lokala variabler samt `Y[inertIndex]`. Om man antar att varje slinga har ett unikt `inertIndex`, är det bara `solve()` kvar.


---


## hEqn.H ##
[code](http://code.google.com/p/openfoam-openmp/source/browse/trunk/hEqn.H)

Körs i den yttre av `for`-looparna

```
{
    solve
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
       DpDt
     + dieselSpray.heatTransferSource()
    );

    thermo.correct();
}
```

Det är endast `thermo.correct()` och `solve()` som anropas, kontrollera dessa.

Parallellisering av `solve()`'s parametrar:

```
fvScalarMatrix *left1, *left2, *left3;
volScalarField *right;

#pragma omp sections
{
    #pragma omp section
        left1 = new fvScalarMatrix(fvm::ddt(rho, h)());
    #pragma omp section
    {
        left2 = new fvScalarMatrix(mvConvection->fvmDiv(phi, h));
        left3 = new fvScalarMatrix(fvm::laplacian(turbulence->alphaEff(), h));
        right = new volScalarField(DpDt + dieselSpray.heatTransferSource());
    }
}

solve( *left1 + *left2 -(*left3) == *right );
```


---


## pEqn.H ##
[code](http://code.google.com/p/openfoam-openmp/source/browse/trunk/pEqn.H)

Körs i den innre `for`-loopen

```
rho = thermo.rho();

volScalarField A = UEqn.A();
U = UEqn.H()/A;

if (transonic)
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *((fvc::interpolate(U) & mesh.Sf()) - fvc::meshPhi(rho, U))
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho/A, p)
         ==
            Sevap
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi = fvc::interpolate(rho)
         *((fvc::interpolate(U) & mesh.Sf()) - fvc::meshPhi(rho, U));

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho/A, p)
         ==
            Sevap
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U -= fvc::grad(p)/A;
U.correctBoundaryConditions();

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
```

Ändrar på många globala variabler så som `rho`, `phi` och `U`. Inkluderar 2 ytterligare filer, som kan ha beroenden i sig.

Vi har även (för mig) nykomligen `DpDt` insmygandes på slutet. Vad är det för rackare?