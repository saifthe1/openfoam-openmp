# Plot #

Varje applikations-solver består oftast av flera ekvationer som körs sekventiellt. Dessa måste bibehålla sin ordning, och alla delar av ekvationerna måste som regel bibehålla sin sekventiella ordning. Dock finns visst utrymme till parallellism, om vissa delar av koden skulle vara oberoende relativt andra tidspunkter.

# Kod #

Den mycket triviala lösningen blir att byta ut huvudloopen, som ser ut som något i stil med:
```
while(runTime.run())
{
    runTime++;
    ...
}
```
Eller
```
while(runTime.loop())
{
    ...
}
```

Mot
```
const int startTime = runTime.startTime().value(),
          endTime = runTime.endTime().value(),
          deltaT = runTime.deltaT().value();
#pragma omp parallel for default(shared) ordered schedule(static, 1)
for(int bigLoopI = startTime; bigLoopI < endTime; bigLoopI += deltaT)
{
    ...
}
```

# Analys av implementering i simpleFoam #

Då simpleFoam har följande struktur:
```
while (runTime.loop())
{
	...

	tmp<fvVectorMatrix> UEqn
	(
		fvm::div(phi, U)	// Använder både U och phi
	  + turbulence->divDevReff(U)
	);
	
	...

	// Nollställer U och phi
	volScalarField AU = (*UEqn)().A();
	U = UEqn().H()/AU;
	phi = fvc::interpolate(U) & mesh.Sf();

	...

	for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
	{
		...
		if (nonOrth == nNonOrthCorr)
		{
			// Uppdatering inför NÄSTA iteration
			phi -= pEqn.flux(); 
		}
	}
	...
	// Uppdatering inför NÄSTA iteration
	U -= fvc::grad(p)/AU; 
	U.correctBoundaryConditions(); 
	
	...

}
```

Och då det inte är mycket kod "censurerat" i början och slutet av koden, ser man att man inte kan påbörja (alls mycket) av en iteration, förän den föregående är helt utförd.