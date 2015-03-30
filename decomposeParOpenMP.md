# OpenMPI #

decomposePar används ihop med composePar för att dela upp indata (se blockMesh och decomposeParDict) och sedan sätta ihop det när det väl har blivit processat.

```
decomposePar
|
|
mpirun -np 4 simpleFoam -parallel
|
|
composePar
```

decomposePar delar upp modellen i n stycken delar beroende på hur man har valt sina regioner. Den lägger de n delarna i n stycken olika kataloger. Vi antar att OpenMPI sedan skickar katalogerna till de n olika noderna. Varje nod kör varsin del och OpenMPI sköter om de uppdateringar som måste göras.

# OpenMP #

Vi skulle behöva göra en egen decomposeParOMP (modifierad decomposePar) som inte delar upp det i kataloger. Vi behöver fortfarande dela upp modellen men vi behöver istället dela upp det i någon slags variabel eller struktur. Det skulle kunna se ut så här (pseudokod):
```
//case är modellen
// och ndelar är en vektor som innehåller de olika delarna vi har delat upp det i
ndelar = decomposeParOMP( case );

#pragma omp parallel for default(shared) //synkronisering av de olika variablerna måste undersökas
    for( int i = 0; i<n; i++)
    {
        //sätt igång solvern på en del av modellen
        simpleFoam( ndelar[i] );
    }
```