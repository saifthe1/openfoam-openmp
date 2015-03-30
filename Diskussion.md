# 4 Diskussion #

## 4.1 Profileringsverktyg ##

Som omnämnt i kapitel 3.1.1.1 förefaller gprof sig ej vara lämpat för parallelliseringsarbete. En tidsåtgång per funktionsanrop kan till vara väldigt bra när man skall optimera ett program. Dock är det inte alls lika relevant när man vill parallellisera ett program.

Ett profileringsverktyg som operar på rad-för-rad basis istället för funktion-för-funktion basis skulle passat vårat projekt mycket bättre. Ett sådant verktyg skulle gärna profilera tiden kring varje block, såväl som varje rad. På så vis skulle det göra precis det jobb gprof gör, fast med större detaljrikedom.

## 4.2 Parallelliseringsarbete ##

### 4.2.1 Ramverk ###

Valet av OpenMP har fungerat bra under projektets gång. Det har varit lätt att använda och levt upp till de förväntningar som ställts på det. Det primära alternativet, POSIX Threads hade eventuellt kunnat ge bättre resultat, men med tanke på den begränsade tidstillgången hade det antagligen inte hjälpt.

Vid en eventuell parallellisering av kärnan hade den lägre overheaden kunnat ge bättre prestanda, men på den nivå som projektet rört sig hade resultatet knappast förbättrats i någon större utsträckning.

### 4.2.2 Testning ###

#### 4.2.2.1 Benchmarks ####

#### 4.2.2.2 Resultatverifiering ####

### 4.2.3 Collaborative work in parallelized programming ###

Funkar det bättre/sämre att programmera flera på samma projekt? Hur pass mycket mer "ömtålig" blir koden?

### 4.2.4 Approach efter hårdvara? ###

Cahcemissar är mindre kritiska i vissa system, osv.

## 4.3 OpenFOAM ##

OpenFOAM har en ovanlig kodstruktur, särskilt på lösarnivå eftersom den är tänkt för fysiker och inte programmerare. Detta gör att koden är svår att analysera eftersom den bryter mot många grundläggande principer för hur programkod bör skrivas.

Exempelvis finns det ett antal globala variabler som används utan dokumentation. Denna typ av kod gör det mycket svårt att analysera beroenden och utröna vilka kodstycken som modifierar samma platser i minnet eftersom den variabel som modifieras kan finnas i en importerad fil, ett stort antal steg längre in i programmet.

Utvecklarna har också valt att överlagra ett antal operatorer så som == och `*`, i vissa fall är det motiverat till exempel när de använt det för att implementera matrismultiplikation. I andra fall är det dock otydligt och svårläst och verkar ha gjorts enbart för att få koden att lika en matematisk ekvation.

Denna typ av "matematisering" gör koden svårhanterlig och är en stor anledning till att det eventuellt hade varit bättre att hitta ett annat projekt alternativt att parallelliser i kärnan istället.