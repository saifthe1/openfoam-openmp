# Kod #
## Orginalkod ##

```
solve
(
    fvm::ddt(rho, Yi)                        // left1   \
  + mvConvection->fvmDiv(phi, Yi)            // left2    |--- left
  - fvm::laplacian(turbulence->muEff(), Yi)  // left3   /
  ==
    dieselSpray.evaporationSource(i)         // right1  \___  right
  + kappa*chemistry.RR(i),                   // right2  /
    mesh.solver("Yi")
);
```

## Parallelliserad kod ##
```
fvScalarMatrix  *left1 = NULL,
                *left2 = NULL,
                *left3 = NULL;
volScalarField *right = NULL;

#pragma omp parallel sections default(shared)
{
    #pragma omp section
    {
        TS_START("YEqn.H > for > if > left1 & right");
        left1 = new fvScalarMatrix(fvm::ddt(rho, Yi)());
        right = new volScalarField((dieselSpray.evaporationSource(i) + kappa*chemistry.RR(i))());
        TS_END("YEqn.H > for > if > left1 & right");
    }
    #pragma omp section
    {
        TS_START("YEqn.H > for > if > left2");
        left2 = new fvScalarMatrix(mvConvection->fvmDiv(phi, Yi));
        TS_END("YEqn.H > for > if > left2");
    }
    #pragma omp section
    {
        TS_START("YEqn.H > for > if > left3");
        left3 = new fvScalarMatrix(fvm::laplacian(turbulence->muEff(), Yi));
        TS_END("YEqn.H > for > if > left3");
    }
}

TS_START("YEqn.H > for > if > solve()");
solve( *left1 + *left2 - *left3 == *right, mesh.solver("Yi") );
TS_END("YEqn.H > for > if > solve()");
```

## Sammansättning ##

För att lättare kunna utföra tester har vi strukturerat den parallella och seriella biten såhär:

```
TS_START("YEqn.H > for > if > solve");
#ifdef PARALLELIZE
  // Parallelliserad kod
#else
  // Orginalkod
#endif
TS_END("YEqn.H > for > if > solve");
```

Vi har verifierat resultatet genom att med md5 hasha all utdata.

# Tidtagningar #

Koden har körts på en 4-kärnig dator. Namngivningen av stycken som är parallelliserade är inspirerad av CSS-selectors, och ser ut som t.ex:
```
YEqn.H > for > if > solve
```
Vilket innebär i `YEqn.H`, i en `for`-loop, i en `if`-sats, finns något vi kallar `solve`, och det är det vi har tagit tid på. `solve` är hela stycket kring solve()-funktionen. Studera kodstyckena ovan för att se hur det hela ter sig.

## Seriell körning med orginalkoden ##

```
dieselFoam/aachenBomb$ dieselFoam | tail -n 3

What                                               Time               Count   
YEqn.H > for > if > solve......................... 1.4973             12      
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 3

What                                               Time               Count   
YEqn.H > for > if > solve......................... 1.501              12      
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 3

What                                               Time               Count   
YEqn.H > for > if > solve......................... 1.5008             12   
```

## Seriell körning av den parallelliserade koden ##
Pragman utkodade och PARALLELIZE satt till 1.

```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.1055             12      
YEqn.H > for > if > left2......................... 0.19928            12      
YEqn.H > for > if > left3......................... 0.19718            12      
YEqn.H > for > if > solve......................... 1.6077             12      
YEqn.H > for > if > solve()....................... 1.1055             12      
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.10746            12      
YEqn.H > for > if > left2......................... 0.20146            12      
YEqn.H > for > if > left3......................... 0.20119            12      
YEqn.H > for > if > solve......................... 1.6214             12      
YEqn.H > for > if > solve()....................... 1.111              12      
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.10823            12      
YEqn.H > for > if > left2......................... 0.20387            12      
YEqn.H > for > if > left3......................... 0.20146            12      
YEqn.H > for > if > solve......................... 1.6246             12      
YEqn.H > for > if > solve()....................... 1.1109             12    
```

## Våran parallella kod ##
### omp\_set\_num\_threads(1) ###
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.10932            12      
YEqn.H > for > if > left2......................... 0.2044             12      
YEqn.H > for > if > left3......................... 0.20757            12      
YEqn.H > for > if > solve......................... 1.7351             12      
YEqn.H > for > if > solve()....................... 1.2134             12      
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.1062             12      
YEqn.H > for > if > left2......................... 0.19843            12      
YEqn.H > for > if > left3......................... 0.19709            12      
YEqn.H > for > if > solve......................... 1.6013             12      
YEqn.H > for > if > solve()....................... 1.0992             12      
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.10664            12      
YEqn.H > for > if > left2......................... 0.19885            12      
YEqn.H > for > if > left3......................... 0.19843            12      
YEqn.H > for > if > solve......................... 1.6013             12      
YEqn.H > for > if > solve()....................... 1.097              12      
```
### omp\_set\_num\_threads(2) ###
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.1776             12      
YEqn.H > for > if > left2......................... 0.3825             12      
YEqn.H > for > if > left3......................... 0.28764            12      
YEqn.H > for > if > solve......................... 1.672              12      
YEqn.H > for > if > solve()....................... 1.2061             12      
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.17727            12      
YEqn.H > for > if > left2......................... 0.38812            12      
YEqn.H > for > if > left3......................... 0.29323            12      
YEqn.H > for > if > solve......................... 1.6795             12      
YEqn.H > for > if > solve()....................... 1.2084             12      
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.18576            12      
YEqn.H > for > if > left2......................... 0.35751            12      
YEqn.H > for > if > left3......................... 0.28088            12      
YEqn.H > for > if > solve......................... 1.6396             12      
YEqn.H > for > if > solve()....................... 1.173              12      
```
### omp\_set\_num\_threads(3) ###
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.39329            12      
YEqn.H > for > if > left2......................... 0.43809            12      
YEqn.H > for > if > left3......................... 0.47465            12      
YEqn.H > for > if > solve......................... 1.7064             12      
YEqn.H > for > if > solve()....................... 1.2084             12      
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.30621            12      
YEqn.H > for > if > left2......................... 0.45615            12      
YEqn.H > for > if > left3......................... 0.42762            12      
YEqn.H > for > if > solve......................... 1.7287             12      
YEqn.H > for > if > solve()....................... 1.2068             12      
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.38802            12      
YEqn.H > for > if > left2......................... 0.45907            12      
YEqn.H > for > if > left3......................... 0.45635            12      
YEqn.H > for > if > solve......................... 1.7532             12      
YEqn.H > for > if > solve()....................... 1.2612             12      
```
### omp\_set\_num\_threads(4) ###
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count
YEqn.H > for > if > left1 & right................. 0.37724            12
YEqn.H > for > if > left2......................... 0.49329            12
YEqn.H > for > if > left3......................... 0.43352            12
YEqn.H > for > if > solve......................... 1.8343             12
YEqn.H > for > if > solve()....................... 1.2566             12
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count
YEqn.H > for > if > left1 & right................. 0.34974            12
YEqn.H > for > if > left2......................... 0.49357            12
YEqn.H > for > if > left3......................... 0.38694            12
YEqn.H > for > if > solve......................... 1.8047             12
YEqn.H > for > if > solve()....................... 1.2445             12
```
```
dieselFoam/aachenBomb$ dieselFoam | tail -n 7

What                                               Time               Count
YEqn.H > for > if > left1 & right................. 0.38617            12
YEqn.H > for > if > left2......................... 0.48187            12
YEqn.H > for > if > left3......................... 0.37564            12
YEqn.H > for > if > solve......................... 1.8635             12
YEqn.H > for > if > solve()....................... 1.2824             12
```

# Slutsatser #

## Förväntan ##
Tidtagningarna på orginalkoden ger en snittlig tid på `(1.4973 + 1.501 + 1.5008) / 3 = 1.4997` sekunder. Tiderna för den parallella koden, om den körs seriellt utan OpenMP ger
```
YEqn.H > for > if > left1 & right:  (0.1055  + 0.10746 + 0.10823) / 3 = 0.1071  (a)
YEqn.H > for > if > left2:          (0.19928 + 0.20146 + 0.20387) / 3 = 0.2015  (b)
YEqn.H > for > if > left3:          (0.19718 + 0.20119 + 0.20146) / 3 = 0.1999  (c)
YEqn.H > for > if > solve:          (1.6077  + 1.6214  + 1.6246)  / 3 = 1.6179  (d)
YEqn.H > for > if > solve():        (1.1055  + 1.111   + 1.1109)  / 3 = 1.1091  (e)
```
Då a, b och c kan köras parallelt helt oberoende av varandra, borde man med parallellisering ha möjlighet att åstakomma en körtid på `MAX(a,b,c) + e = 0.2015 + 1.1091 = 1.3106` sekunder. Dock har man ju givetvis någon form av overhead vid parallellisering, så resultatet kan förväntas bli en något längre tid i den verkliga implementeringen.

Under ideella förhållanden skulle den här parallelliseringen kunna ge `100*(1.4997 - 1.3106) / 1.4997 = 12.6092`% uppsnabbing, av just det här stycket. Då solve-stykcet i YEqn bidrar till ca 11% av den totala körtiden vid en kortare körning, innebär den här parallelliseringen en uppsnabbning på totalt 1.3% av den totala körtiden.

Dock bör det märkas väl att YEqn är den enda delen som ökar mest vid en längre körning relativt andra delar, samt att filen innehåller en loop som med största sannolikhet går att parallellisera.

## Resultat ##
Vi får sämre resultat ju fler trådar vi aktiverar. Detta är mycket märkligt, specielt i det fallet med 4 trådar, eftersom det ju bara är maximalt 3 trådar som kan utföra det här stycket.