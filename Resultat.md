# 3. Resultat #

Det praktiska arbetet har nästan uteslutande berört OpenFOAM och det är endast detta som tas upp i rapporten.

## 3.1 Parallellisering av OpenFOAM ##
OpenFOAM är mycket modulärt och parallellisering i de olika delarna medför väldigt skilda resultat. En bedömning av hur förändringarna i koden påverkar projektet i helhet behövs för flertalet aspekter för att avgöra vilken del som är lämpligast att parallellisera.

OpenFOAMs lösare är anpassade för att författare av dessa inte ska behöva några vidare programmeringskunskaper för att skriva dem. Således medför förändringar i dessa möjligen besvär för eventuella fysiker, matematiker eller kemister som uppdaterar applikationslösare. Dock kan det diskuteras huruvida ett fysikalisk beskrivning av ett flöde behöver "patchas".

Att parallellisera en lösare ger endast effekt för just den lösaren, och möjligen några ytterligare som brukar dennes ekvationer. Effekten av parallellisering i en lösare blir hur som helst begränsad till mycket få applikationer. Om en viktig del av OpenFOAMs kärna däremot skulle parallelliseras skulle en uppsnabbning av samtliga lösare uppnås.

Kontinuitetsproblem kan uppstå vid framtagandet av en parallelliserad version av OpenFOAMs kärna. En ny version skulle innebära förändringar i kärnan, och en förändring i kärnan skulle innebära att den parallella versionen skrivs över. OpenCFD uppdaterar dock inte fysikaliska omständigheter, och således behövs inga förändringar i lösare vid uppdatering av OpenFOAM. Därför skulle inte dessa problem uppstå vid parallellisering av en lösare och man slipper helt kontinuetitsproblem vid framtagandet av en parallelliserad version av en lösare.

För projektets ändamål valdes att parallellisera lösaren dieselFoam. Detta dels på grund av att parallellisera en lösare är mindre komplext än att parallellisera en algoritm i kärnan (kärnmekanism), som till exempel gauss-elimination, matrismultiplikation, etc.. Men framför allt valdes lösare för att det tar betydligt kortare tid - ca 30 sekunder istället för ett par timmar - att kompilera en lösare, än att kompilera om kärnan.

### 3.1.1 Analys ###

För att bedömma vilken form av parallellisering som skulle ge mest effekt är profilering och analys av koden en självklarhet. Parallellisering av en lösare är problematiskt då den är ämnad att lösa ett fysikaliskt problem, och därmed förenklat kan beskrivas som en serie uträkningar. Dessa uträkningar måste ofta utföras i korrekt ordning.
Det är nu av värde att titta på koden för dieselFoam:

```
	while(runTime.run()) //tidsloopen
	{

        dieselSpray.evolve();

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
    }

```
Ovanstående kodstycke ska ses som pseudokod. Den verkliga koden står att finna i appendix.
Vid en första anblick lägger man märke till den nästlade for-loopen. Iterationer av detta slag har stor potential vid parallellisering. Men det är inte tillräckligt att identifiera potential då den potentialen är beroende av vad tidsanalysen kommer fram till. Således tar vi analysverktyg till hjälp.

#### 3.1.1.1 GProf ####

Analys m.h.a. GProf gav oss en tabell med funktioner sorterade efter tidsåtgång i procent. Genomgång av dessa funktioner gav vid handen att det var kärnmekanismer som listades. Dessa var inte av intresse då strategin var att fokusera på en lösare istället för på kärnmekanismer. Dessutom tog dessa funktioner upp stor del av tiden därför att de anropas ofta, inte för att varje körning av funktionen tar mycket tid. Detta är inte vad vi söker efter när vi vill parallellisera.
GProf är något av industristandard när det kommer till profilering. Dock utför det inte alls det arbete som behövs för att utföra parallellisering i en lösare då den mäter tidsåtgången vid funktionsanrop. Mer om detta i avsnitt 4 Diskussion.


#### 3.1.1.2 Rad-för-rad timestamps ####

Istället användes en mer brute-force metod, med tidtagning varje instruktion i main()-funktionen.

### 3.1.2 Genomförande ###

Efter att ha fastställt hur stor tidsåtgång varje del av koden kräver, ses till vilket beroende som finns mellan olika delar av koden. Endast ett fåtal stycken utan beroenden finnes, och inga uppenmbara parallelliseringsmöjligheter hittas. Man kan argumentera för att det skulle vara lönsamt att låta de olika slingorna i YEqn köras parallelt, men det är endast 4 i varje simulation, vilket gjorde att detta förbisågs.

#### 3.1.2.1 Uppdelning av solve()'s parametrar ####

Solvefunktionen anropas med argument i form av ekvationer. Dessa ekvationer byggs upp av mindre delar som beräknas var för sig. I originalimplementationen beräknades ekvationsdelarna inom solvefunktionens anrop, detta visade sig vara en möjlighet till parallellisering.

Efter att ha analyserat beräkningarna som gjordes för att skapa ekvationsdelarna fastslogs det att de inte modifierade några gemensamma variabler utan bara den variabel de var till för att räkna ut. Genom att flytta ut beräkningen av ekvationsdelarna ur solveanropet gick det att beräkna ett antal av delarna parallellt.

**Kod**

```
solve
(
    fvm::ddt(rho, h)                            // left1 \
  + mvConvection->fvmDiv(phi, h)                // left2  | ------ Vänsterledet
  - fvm::laplacian(turbulence->alphaEff(), h)   // left3 /
 ==
   DpDt                                         // right1 \_______ Högerledet
 + dieselSpray.heatTransferSource()             // right2 /
);
```

#### 3.1.2.2 Parallellisering längsmed tidsaxeln ####

Varje simulation i OpenFOAM utför ett visst antal tidsiterationer, vilket ger ett visst potential till parallellism. För att testa detta tillvägagångssätt användes simpleFoam, vilket är en mycket minder komplex applikationslösare än t.ex. dieselFoam. Dock visade sig resultaten långt ifrån tillräckkliga eftersom endast små delar i början och slutet av varje tidsiteration i den algoritm som simpleFoam (SIMPLE) använder kunde parallelliseras. Tidigt i varje iteration behöver man läsa från en variabel som skrivs till i slutet av föregående iteration. Se:

```
while (runTime.loop())
{
        ...
        tmp<fvVectorMatrix> UEqn
        (
                fvm::div(phi, U) + turbulence->divDevReff(U)   // Använder både U
        );
        ...
        U = UEqn().H()/AU;                                     // Nollställer U
        ...
        U -= fvc::grad(p)/AU;
        U.correctBoundaryConditions();                         // Uppdatering inför NÄSTA iteration
        ...
}
```

Betrakta hur man genomgående och på ett mycket retfullt vis räknar ut saker i rätt ordning. Man beräknar `UEqn`, som använnds för att initiera `U` under varje iteration. Uträkning av `UEqn` kräver dessutom att `U` från föregående iteration är fullständigt utförd. Eftersom `U` inte beräknas förän i slutet av varje iteration och `UEqn` beräknas i ett tidigt skede av varje iteration. Finns inte mer än ett fåtal iterationer att tjäna på parallellisering längs tidsaxeln i fallet av simpleFoam.

Det finns dock ett möjligt potential; om man finner en lösare där varje instruktion i algoritmen endast påverkar någon eller några nästkommande instruktioner. Då projektet var under stor tidspress hann inte detta undersökas närmare.

## 3.2 Resultatredovisning ##

Resultaten härrör från parallellisering i dieselFoam genom uppdelning av solve()'s parametrar.

### 3.2.1 Kodbriefing ###

`[Kod från wiki]`

Då vi varit tvugna att i viss mån ändra på algoritmen för att uppnå parallellism, kan inte en körning med en ensam tråd i den parallelliserade version anses vara densamma som en körning med orginalkoden. Därav brukas förkompilationsdirektiv för att det ska vara lätt att växla mellan orginalkoden och den parallelliserade versionen samtidigt som man får korrekta tidsuppmätningar.

### 3.2.2 Miljöbeskrivning ###

Tester har gjorts med det färdiga Aachenbomb-exemplet som följer med OpenFOAM-distributionen.

#### 3.2.2.1 Tutorialbeskrivning ####

För att ekvationerna som de olika anropen till solve() löser ska ta en icke försumbar tid att ställa upp, måste problemet man försöker simulera vara tillräckligt stort. Samtidigt för inte körningen ta allt för lång tid för att testen skulle hinna bli behandlade i projektet.

För att öka problemstorleken har modellens storlek ökats, genom konfigurationsfilen constant/polyMesh/blockMeshDict. Dessutom har noggrannheten för varje ekvation ökats, genom konfigurationsfilen system/fvSolution. Därtill har man alltså minskat antalet tidsiterationer genom system/controlDict.

#### 3.2.2.1 Hårdvaruspecifikation ####

Beskrivning av hårdvara vi kör testerna på:

Intel core2 quad, GB, ubuntu 10.04.

### 3.2.3 Redovisning ###

Redovisning av testkörningarna, tabeller. Inget mer behövs?

### 3.2.4 Resultatsammanfattning ###

http://code.google.com/p/openfoam-openmp/wiki/YEqn