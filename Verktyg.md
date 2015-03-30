# 2 Verktyg #

Parallellprogrammering innebär djup analys och djup förståelse, ofta även analys av stora mängder kod. Därför är goda verktyg ett måste för att utföra sådant arbete. Dessutom är parallellisering en komplicerad process som kräver att man behärskar många tekniker inom ämnesområdet. Detta inte minst då trådhantering i olika operativsystem ter sig helt olika.

## 2.1 OpenMP ##

OpenMP är ett parallelliseringsramverk för språken C, C++ och Fortran, som underlättar många av de mer komplicerade uppgifterna under ett paralllelliseringsarbete.
"OpenMP is a portable, scalable model that gives shared-memory parallel programmers a simple and flexible interface for developing parallel applications for platforms ranging from the desktop to the supercomputer." (1)

En viktig insikt är att OpenMP inte är en produkt, utan endast en specifikation. Det är upp till varje kompilator för sig att implementera OpenMP. OpenMP-specifikationen utvecklas och godkänns av OpenMP Architecture Review Board, som består av representanter från många IT-jättar så som Microsoft, IBM, Cray, HP, AMD, Intel, Texas Instruments och Oracle. (1)

OpenMP opererar på enskilda kodblock, och inte på ett helt program. Ett parallellt block har ett tillhörande lag av trådar, som utför den inneslutna koden. Om inga direktiv ges vid skapandet av blocket, kommer samtliga trådar i laget att utföra koden. Det är lätt att styra hur trådar skall hantera kod inom ett parallellt block genom att skapa ett kontrollblock inom det parallelliserade blocket. Kontrollblock kan tillexempel vara till för att synkronisera, dela upp arbetet eller specifiera en uppgifts-orienterad struktur.

(1) http://openmp.org/wp/about-openmp/ 2010-05-06

### 2.1.1 Pragma ###

För C och C++ nyttjar sig OpenMP av förkompileringsdirektivet pragma för att skapa och styra sina parallelliserade block. Pragma är en förkortning av pragmatisk och är till för att låta C/C++-programmerare få tillgång till maskinspecifika egenskaper.

Detta är en stor fördel då kompilatorn ignorerar dessa pragmatiska direktiv om den inte vet vad de innebär. Således kan man kompilera kod som innehåller OpenMP-instruktioner på en kompilator som inte stödjer OpenMP och fortfarande få ett fungerande program som resultat.

### 2.1.2 Worksharing ###

En av de viktigaste egenskaperna med OpenMP är att dess möjlighet att tillhandahålla en otroligt väl skalande parallellism. Detta genom att de tidigare nämnda uppgifts-uppdelande kontrollblocken, och då specifikt möjligheten att låta varje iteration av en for-loop skötas parallellt. Se de två nedanstående kodstyckena som exemplifierar detta.

**Kodsektion A**
```
#pragma omp sections
{
    #pragma omp section
    {
        // Work
    }
    #pragma omp section
    {
        // Work
    }
}
```

**Kodsektion B**
```
#pragma omp for
for(int i=0; i<N; ++i)
{
    // Work
}
```

Om man antar att tidsåtgången för att skapa trådar är försumbar kan man anta att A skulle köras dubbelt så fort på en multicore-processor som en enkärnig processor. B skulle dock köras C gånger snabbare, där C är antalet processorkärnor på den arkitektur man exekverar programmet.

Då antalet trådar ej är kopplat till programmet, utan till arkitekturen den körs på, kan man använda B på vilken dator som helst, och förvänta sig en uppsnabbning med faktor C.

## 2.2 Profilering ##

Profilering motsvara dynamisk kodanalys. En insikt i vad detta innebär fås lättast genom att betrakta motsatsen - statisk kodanalys. Satisk kodanalys innebär att man helt enkelt analyserar koden, och betraktar hur programfödet ser ut. Dynamisk kodanalys innebär alltså att man _inte_ tittar på koden, utan bara betraktar hur programmet beter sig när det faktiskt är under exekvering.

Även de programmerare som tror sig aldrig ha profilerat ett program, kan nog säkerligen minnas ett tillfälle då denne tog tid på en programkörning. Detta innebär således att programmeraren har profilerat programmet då denne analyserade programexekveringen, om än på en väldigt detaljfattig nivå.

### 2.2.1 GProf ###

För att underlätta profileringsarbetet har som tidigare nämnts gprof använts. Verktyger är enkelt, snabbt och har väldigt många inställningar för att finjustera nogrannhet. Programmet mäter tidsåtgång och antal anrop för varje funktion, samtidigt som den spara information om anropsgrafen.

## 2.3 OpenFOAM ##

OpenFOAM är ett ramverk för att konstruera och definiera simuleringsmiljöer, och tillhandahåller verktyg för att skapa representationer av den fysiska miljön programmet skall köras i och för att visualisera utdata.  Dessutom ges goda möjligheter att representera en fysikalisk eller kemisk process genom att skapa s.k. lösare.

En fysikalisk eller kemisk process kan givetvis ske under många olika förhållanden och i många olika miljöer. Därför kan man med OpenFOAM lätt konfigurera flera olika saker inför en simulation. Förutom att definiera den fysiska miljön och vad som påverkar den, kan man även finjustera med hur stor precision beräkningar under simulationen skall utföras.

### 2.3.1 Programstruktur ###

![http://openfoam.com/docs/user/img/user0x.png](http://openfoam.com/docs/user/img/user0x.png)

(`http://openfoam.com/docs/user/img/user0x.png`)

Figuren ovan visar den verkliga arbetsgången för en simulering.
  * User Applications och Standard Applications kallas i texten lösare.
  * Utilities och Meshing Tools kallas verktyg.
  * OpenFOAM C++ Library kallas kärnmekańismer.
  * Solving kallas simuleringsfall.
  * Pre-processing och Post-processing kommer bara nämnas där det är absolut nödvändigt för förståelse.

#### 2.3.1.1 Lösare ####
Lösare är program som löser specifika problem inom flödesdynamiken. Dessa använder sig av verktyg som har hand om datamanipulering och algebraiska beräkningar.
Det finns en uppsjö av lösare som kommer med i programvaran. Denna rapport kommer behandla tre stycken lösare simpleFoam, dieselFoam och pisoFoam. Källkoden för dessa tre ges i appendix.

Det riktigt intressanta med lösare i OpenFOAM är att vem som helst med kunskaper i flödesdynamik ska kunna sätta ihop en egen. Grundkunskaper i programmering ska vara det enda kravet för att klara av detta. **OpenCFD säljer kurser**
simpleFoam tas som exempel för att illustrera förfarandet. Nedan följer ett utdrag från openfoamwiki.net
http://openfoamwiki.net/index.php/The_SIMPLE_algorithm_in_OpenFOAM -Wuilbert Lopez 2010-05-10 00.50
som förklarar vad simple algoritmen gör. Jämför den med koden som följer direkt efter.

> The SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) allows to couple the Navier-Stokes equations with an iterative procedure, which can be summed up as follows:

  1. Set the boundary conditions.
  1. Solve the discretized momentum equation to compute the intermediate velocity field.
  1. Compute the mass fluxes at the cells faces.
  1. Solve the pressure equation and apply under-relaxation.
  1. Correct the mass fluxes at the cell faces.
  1. Correct the velocities on the basis of the new pressure field.
  1. Update the boundary conditions.
  1. Repeat till convergence.

```
while (runTime.loop())
{
    ...
    p.storePrevIter(); //behövs för att kunna göra punkt 4
    // Pressure-velocity SIMPLE corrector
    {
2&3     #include "UEqn.H"
4&5     #include "pEqn.H"
    }
    ...
    turbulence->correct();
    ...
8   #include "convergenceCheck.H"
    ...
}
```
Punkterna 1-8 i algoritmen matchar punkterna i koden. Man ser att koden följer algoritmen väldigt väl. Man måste naturligtvis sätta sig in i OpenFOAM och vad den kräver för att uträtta en viss uppgift men som synes av koden är det lite som att bygga med Legoklossar. Alla delar finns där redan, dvs alla verktyg. En kort genomgång av koden är nu på sin plats.
Översta raden är en while-loop, i fortsättningen kallar vi den för tidsloopen. Man anger i en separat fil (se 2.3.1.4) hur länge man vill att simuleringen ska köras och tidsloopen sköter detta åt oss.
UEqn.H och pEqn.H sköter ekvationslösning åt oss. I dessa finns en funktion benämnd solve(). Den beskrivs närmare i 2.3.1.3.


#### 2.3.1.2 Verktyg ####

Alla problem som innebär en modell i 2-D eller högre behöver minst ett datagenereringsverktyg, blockMesh. Det tar en fil med beskrivning av modellen och bygger den. Det gör blockMesh m.h.a. kärnmekanismer.
Efter att en simulering har körts finns det naturligtvis verktyg för att analysera simuleringen. Det kanske viktigaste av dessa verktyg är paraFoam. Det visar modellen grafiskt och man kan via kontrollknappar stega sig igenom simuleringen. Om en del av simuleringen har gått snett är det relativt enkelt att se det i paraFoam. Det finns även alternativ för att kolla på en specifik variabel såsom tryck över hela modellen.

#### 2.3.1.3 Kärnmekanismer ####

Om ändringar i lösare är enkelt så är ändringar i OpenFOAM kärnan desto svårare. Både lösare och verktyg använder sig av kärnmekanismer. Ett exempel på en kärnmekanism är solve. Den får som argument ekvationer (i form av matriser) som beskriver någon egenskap i modellen (t.ex. hastighet och tryck) och beräknar lösningen.

#### 2.3.1.4 Simuleringsfall ####

Lösare och blockMeshverktyget behöver konfigurationsfiler för att uträtta något. Ett simuleringfall är då helt enkelt en mapp som innehåller dessa konfigurationsfiler. Mappen måste se ut på följande vis.

> system ska ha tre filer, controlDict, fvSolver och fvSchemes. I controlDict styr du den simulerade tiden. fvSolver ger dig möjlighet att förändra ekvationslösarnas tolerans. fvSchemes anger vilka ekvationslösare som ska användas.
Om exempelprogrammet kräver komplicerade kommandon kommer man finna ett script kallat Allrun. Man kan antingen köra det eller öppna det och lära sig hur exempelprogrammet ska köras.

### 2.3.2 Hur det funkar ###

För att ge en bättre förståelse skapar vi ett eget simuleringsfall. Vi kallar den Elektrogarden. Det är en simulering av Elektrogården på Chalmers Campus Johanneberg. Modellen ska förutom rudimentär form på gården också simulera vind, regn och väggarnas ytbeläggning.2010-05-09 problem när jag la till trappor, tunnel och glasvägg; se utdelat dokument Elektrogården -Wuilbert Lopez 2010-05-09 20.28

#### 2.3.2.1 Preprocessing ####

Strukturen för Elektrogarden är:
```
Elektrogarden
|
---- system
|
---- constant
---- polyMesh
---- blockMeshDict
```
visa vad som finns i blockMeshDict? eller appendix? -Wuilbert Lopez 5/9/10 12:46 PM

> lättare att göra en bild -Wuilbert Lopez 5/9/10 12:45 PM

#### 2.3.2.2 Running ####

kör man icoFoam.

#### 2.3.2.3 Postprocessing ####

paraFoam