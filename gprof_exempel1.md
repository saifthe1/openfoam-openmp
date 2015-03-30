# Orginal #

```
#include <stdio.h>

/* Computes the length of Collatz sequences */

unsigned int step (unsigned int x)
{
    if (x % 2 == 0)
    {
        return (x / 2);
    }
    else
    {
        return (3 * x + 1);
    }
}

unsigned int nseq (unsigned int x0)
{
    unsigned int i = 1, x;
  
    if (x0 == 1 || x0 == 0)
      return i;

    x = step (x0);

    while (x != 1 && x != 0)
    {
        x = step (x);
        i++;
    }

    return i;
}

int main (void)
{
    unsigned int i, m = 0, im = 0;

    for (i = 1; i < 500000; i++)
    {
        unsigned int k = nseq (i);

        if (k > m)
        {
            m = k;
            im = i;
            printf ("sequence length = %u for %u\n", m, im);
        }
    }
    return 0;
}
```
# Profilering via gprof #
```
Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  us/call  us/call  name    
 55.92     68.92    68.92 39999980     1.72     3.07  nseq(unsigned int)
 43.71    122.78    53.87 1248562184     0.04     0.04  step(unsigned int)
  0.37    123.24     0.46                             main
```

# Försök till prestandaökning #
## Optimering ##

```
unsigned int step (unsigned int x)
{
    if (x % 2 == 0)     
    {           
        return (x / 2);
    }
    else
    {                              
        return (3 * x + 1);
    }                  
}                          
               
unsigned int nseq (unsigned int x0)
{                                 
    unsigned int i = 1, x;
                         
    if (x0 == 1 || x0 == 0)
      return i;

    x = step (x0);
 
    while (x != 1 && x != 0)
    {          
        x = step (x);
        i++;                          
    }
                                    
    return i;
}
```
Bytes mot
```
unsigned int step (unsigned int x)
{
    return (x & 1) == 0
        ? x / 2
        : 3 * x + 1;
}

unsigned int nseq (unsigned int x0)
{
    unsigned int i = 1;
    if (x0 == 1 || x0 == 0)
      return 1;

    for(unsigned int x = step(x0);
        x != 1 && x != 0;
        x = step (x), ++i
        )
    ;

    return i;
}
```

### Resultat ###
```
Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  us/call  us/call  name    
 54.93     64.94    64.94 39999980     1.62     2.94  nseq(unsigned int)
 44.64    117.70    52.77 1248562184     0.04     0.04  step(unsigned int)
  0.43    118.21     0.51                             main
```

## Parallellisering ##

```
int main (void)
{
    unsigned int i, m = 0, im = 0;

    for (i = 1; i < 500000; i++)
    {
        unsigned int k = nseq (i);

        if (k > m)
        {
            m = k;
            im = i;
            printf ("sequence length = %u for %u\n", m, im);
        }
    }
    return 0;
}
```
Bytes mot
```
int main (void)
{
    unsigned int i, m = 0, im = 0;

    #pragma omp parallel for private(i) default(shared)
    for (i = 1; i < 500000; i++)
    {
        unsigned int k = nseq (i);

        #pragma omp critical
        {
            if (k > m)
            {
                m = k;
                im = i;
                printf ("sequence length = %u for %u\n", m, im);
            }
        }
    }
    retu
```