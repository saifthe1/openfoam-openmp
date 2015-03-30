# 1. Introduktion #

Projektet började med en hel del instudering och efterforskning för att bättre kunna behandla ämnesområdet. Efterforskningar gjordes bland annat kring parallelliseringsramverk och profileringsverktyg. Instuderingen syftade till att ge en bättre förståelse för parallellprogrammering och en närmare bekantskap med programmeringsspråket.

Under slutfasen av projektet har vi i samråd med Volvo Technology tagit oss an parallellisering av OpenFOAM. Det är en programvara som kan hantera godtyckliga simuleringsproblem gällande flödesdynamik. Programmet utvecklas av företaget OpenCFD vars affärsidé är att sälja service och utbildning i OpenFOAM, samtidigt som man låter källkoden vara öppen.

## 1.1 Bakgrund ##

Sedan datorn uppfanns har dess utveckling som bekant följt Moores Lag: ”Antalet transistorer som kan placeras på en yta fördubblas vartannat år.” Detta har inneburit en motsvarande prestandaökning för processorer och andra datorkomponenter.

I många år har denna ökning hållit i sig och datorernas hastighet har ökat exponentiellt. I början på 2000-talet har denna trend stött på ett hinder, istället för att utrymmet på kiselchipet hindrar prestandaökning är det värmeutvecklingen från komponenterna som stoppar utvecklingen. Detta är ett problem som varit svårt att avhjälpa med konventionella tekniker och man har därför fått ta fram nya lösningar.

Istället för att öka antalet transistorer i varje processor har många tillverkare valt att sätta flera processorer i samma dator. Dessa delar sedan på minne och andra resurser i datorn men kan arbeta med varsin uppgift parallellt.

I servrar och stordatorer för forskning har detta realiserats genom flera helt självständiga processorer med egna cacheminnen och kommunikationsvägar till resten av datorn. I persondatorer är det vanligaste istället att man istället har flera processorkärnor, den del av processorn som utför själva beräkningarna, som delar på cacheminne och kommunikationsmöjligheter.

### 1.1.1 Parallellismens problematik ###

En vanlig enkärnig dator kan bara utföra en beräkning i taget, men genom att den arbetar så fort kan den rotera olika uppgifter så att den verkar arbeta med dem parallellt. Detta gör att du kan leva i illusionen av att Spotify spelar musik samtidigt som Firefox laddar en hemsida, trots att all kod körs sekventiellt.

En flerkärnig dator kan däremot utföra flera beräkningar samtidigt. Parallellt exekverad kod introducerar en osäkerhetsfaktor i och med att instruktionernas ordningsföljd inte är deterministisk på samma sätt som när koden körs på den enkärniga datorn. Varje förändring från seriellt exekverad kod till parallellt exekverad kod måste därför kontrolleras noga. Det är naturligtvis bättre att från början skriva ett program så att det kan köras parallellt.

## 1.2 Syfte ##

Projektet har syftat till att undersöka möjligheterna för mjukvaruutvecklare, som inte tidigare specialiserat sig i ämnet, att parallellisera kod skriven av någon annan. Detta för att kunna sammanställa tankar, problem och eventuella tips till nästa generations
utvecklare.

## 1.3 Problemformulering ##

För att ge relevanta och konkreta exempel på arbetsgången i ett parallelliseringsarbete har OpenFOAM parallelliserats. Under arbetet ämnades besvara frågor så som:

  * Vilka är de huvudsakliga problem eller principer involverade i parallellisering?
  * Finns det möjlighet att formulera en generell guide för parallellisering?
  * Vilka verktyg skulle behöva utvecklas för att göra det lättare att parallellisera?
  * Vilken del av arbetet är mest tidskrävande?
  * Hur långt klarar man sig utan förståelse för själva algoritmen som koden implementerar?

## 1.4 Fokus & Avgränsningar ##

Rapporten fokuserar på mjukvara och kommer inte diskutera hårdvara i andra syften än för att upprätthålla en korrekt och relevant diskussion kring prestanda - som ju faktiskt är själva syftet med att parallellisera. För att göra detta måste vi belysa vissa termer så som cache- och RAM-minnen, processorkärnor och registerhantering vid trådväxling.

Rapporten avser att ge tankar rörande parallellisering av kod utan att programmeraren behöver ha en djupare förståelse kring problemet. Därmed menas inte att ge direkta riktlinjer, och långt ifrån en komplett lösning. Problem och möjligheter kommer presenteras, men ingen totallösning kan ges.

Rapporten behandlar ej parallelliseringsarbete generellt, utan presenterar istället det arbete som under projektets gång utförts i OpenFOAM för att ge exempel på olika problem och lösningar specifika för just dess kodbas. De resultat som uppnås och de motgångar som uppstår kommer sedan användas för diskussion som förhoppningsvis kan appliceras på mer generella problem.

## 1.5 Metod ##

Valet av programmeringsspråk föll på C och C++. Dessa två språk är nära besläktade och kompileras till maskinkod. Jämfört med ett språk som Java som körs på en virtuell maskin ger detta bättre prestanda vilket gör att C och C++ ofta används till beräkningsintensiva uppgifter, just den typ av uppgifter där prestandaökningar genom parallellisering är mest intressant.

Detta val avklarat måste ett ramverk för parallellisering väljas. OpenMP valdes eftersom det ger ett tydligt och lättförstått stöd till programmeraren, det krävs endast några få rader för att parallellisera ett större stycke. Eftersom parallelliseringen sker med hjälp av så kallade pragman, kodstycken som enbart hanteras av kompilatorn om OpenMP-biblioteket är aktiverat, är det enkelt att testa koden i seriellt läge. Detta gör det enklare att avgöra om ett problem kommer ur själva logiken eller ur det parallella flödet.

Det främsta alternativet till OpenMP var Posix Threads som är ett sätt att på väldigt låg nivå hantera trådar i UNIX-system. Detta kan ge bättre parallelliseringsmöjligheter men kräver att man detaljstyr koden och trådningen mycket hårdare. (1)

För att profilera och analysera körning program bestämdes att använda gprof, ett profileringsverktyg utvecklat av organisationen GNU. Detta inte att förväxla med UNIX-varianten av samma program med samma namn. gprof presenterar tydligt och detaljerat vilka funktionsanrop som upptar mest tid under en körning.

(1) http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.13.1493&rep=rep1&type=pdf

## 1.6 Resultatsammanfattning ##

Det viktigaste resultatet från arbetet med OpenFOAM har varit insikten att vissa av de verktyg som används vid parallelliseringsarbete är långt ifrån fullkomliga. Det främsta exemplet är gprof, som i stor utsträckning kan ses som industristandard vid profilering. Dock var GProf inte till alls stor nytta i det här fallet eftersom det ger information endast om hur funktioner behandlas. Då vi inte alls behöver parallellisera en hel funktion utvecklades ett eget verktyg för att tillhandahålla profilering kring egendinierade sektioner.

Andra viktiga slutsatser har varit de svårigheter som finns vid uppskattning av en möjlig tidsvinst av parallellisering, och att i seriell kod hitta stycken som är parallelliseringsvänliga.

Ett av de största problemen som uppstått under projektets gång har varit sporadiska skillnader vid identiska testkörningar. Dessa uppstår inte på grånd av något fel eller egentligen något alls underligt, utan kan härledas till att vi har operativsystem och många andra processer körandes i bakgrunden. Man vill gärna låta sitt testprogram köras under en tidsrymd av ett fåtal sekunder, men gör man det så kommer variationen att vara så stor att en uppsnabbning blir omöjlig att mäta.

## 1.7 Disposition ##

Rapporten kommer under kapitel `2. Verktyg` att behandla de olika programvaror som använts för att utföra parallelliseringsarbetet. Kapitlet behandlar även OpenFOAM till sådan grad att läsaren skall kunna ta till sig resultat och diskussion på ett givande sätt.

Under kapitel `3. Genomförande och resultat` presenteras först projektets arbetsgång och sedan de resultat projektet slutligen resulterade i. Samtliga fakta presenteras i detta avsnitt.

Kapitel `4. Diskussion` bearbetar de resultat som presenterats i det föregående kapitlet. Här ledes läsaren in i en argumentation till varför nya verktyg bör tas fram och hur dessa skall fungera.