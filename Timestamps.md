# Timestamps #
Tagna med vår nya ts klass!

|Filnamn|Lokalitet|Tidsåtgång|Körningar|
|:------|:--------|:-----------|:---------|
|dieselFoam.C|> while|27.012|1 |
|dieselFoam.C|> while > for|16.892|2 |
|dieselFoam.C|> while > for > for > pEqn|10.361|4 |
|pEqn.H|> if/else|08.333|4 |
|pEqn.H|> else > for|06.770|4 |
|pEqn.H|> else > for > pEqn.solve()|05.919|4 |
|dieselFoam.C|> while > UEqn.H|04.635|2 |
|dieselFoam.C|> while > for > YEqn.H|03.772|2 |
|YEqn.H|> for|03.743|2 |
|YEqn.H|> for > if > solve|03.703|8 |
|dieselFoam.C|> while > turbulence->correct()|03.644|2 |
|dieselFoam.C|> while > for > hEqn.H|02.756|2 |
|UEqn.H|> UEqn|02.359|2 |
|UEqn.H|> if|02.276|2 |
|pEqn.H|> else > phi =|01.563|4 |
|hEqn.H|> solve()|01.438|2 |
|hEqn.H|> thermo.correct()|01.318|2 |
|dieselFoam.C|> while > dieselSplay.evolve()|01.273|2 |
|pEqn.H|> rUA`*`UEqn.H()|00.592|4 |
|pEqn.H|> else > for > pEqn|00.591|4 |
|rhoEqn.h|> `*`|00.499|6 |
|pEqn.H|> rUA`*`fvc::grad(p)|00.459|4 |
|pEqn.H|> fvc::DDt()|00.416|4 |
|rhoEqn.h|> solve|00.400|6 |
|pEqn.H|> rhoEqn.H|00.332|4 |
|pEqn.H|> else > for > if|00.255|4 |
|dieselFoam.C|> while > compressibleCourantNo.H|00.174|2 |
|dieselFoam.C|> while > rhoEqn.H|00.167|2 |
|pEqn.H|> compressibleContinuityErrs.H|00.157|4 |
|dieselFoam.C|> while > turbulent time scale|00.142|2 |
|compressibleContinuityErrs.H|> globalContErr|00.082|4 |
|rhoEqn.h|> for|00.075|6 |
|rhoEqn.h|> for > `*`|00.075|30|
|compressibleContinuityErrs.H|> sumLocalContErr|00.060|4 |
|pEqn.H|> 1.0/UEqn.A()|00.059|4 |
|rhoEqn.h|> Sevap|00.024|6 |
|dieselFoam.C|> while > spraySummary.H|00.018|2 |
|dieselFoam.C|> while > chimistry.solve()|00.018|2 |
|compressibleContinuityErrs.H|> totalMass|00.014|4 |
|pEqn.H|> thermo.rho()|00.013|4 |
|dieselFoam.C|> while > thermo.rho()|00.009|2 |
|dieselFoam.C|> while > readPISOControls.H|00.000|2 |
|YEqn.H|> mvConvection|00.000|2 |
|dieselFoam.C|> while > setDeltaT.H|00.000|2 |
|pEqn.H|> U.correctBoundaryConditions()|00.000|4 |
|compressibleContinuityErrs.H|> cumulativeContErr|00.000|4 |