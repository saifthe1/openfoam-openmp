# Introduction #

Försök att få YEqn.H:s tre olika segment att köra längre än 0.5 sekunder på Pers dator med 4 kärnor

# Details #

På wills dator finns ../dieselFoam/will.txt där jag sparat flera körningar av tre olika modellupplösningar

Jag kör med tre olika modellupplösningar
```
(41 41 100), (51 51 110) och (82 82 200)   
```

Den sista modellupplösningen kan inte exekveras till simuleringens slut utan klarar bara av de fyra första exekveringarna av YEqn. Därför presenterar jag bara fyra första exekveringarna för alla tre upplösningar.

Sammanfattningsvis ökar exekveringstiden i YEqn.H när man ökar
modellupplösningen, detta förvånar mig, jämför för den parallelliserade YEqn.H:
## (41 41 100) ##
```
What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.49478            4       
YEqn.H > for > if > left2......................... 0.48792            4       
YEqn.H > for > if > left3......................... 0.4694             4       
YEqn.H > for > if > solve......................... 1.09               3       
YEqn.H > for > if > solve()....................... 0.78614            4        
```
## (51 51 100) ##
```
YEqn.H > for > if > left1 & right................. 0.61952            4       
YEqn.H > for > if > left2......................... 0.85454            4       
YEqn.H > for > if > left3......................... 0.84286            4       
YEqn.H > for > if > solve......................... 2.1329             3       
YEqn.H > for > if > solve()....................... 1.6302             4      
```
## (82 82 200) ##
Följande utskrift fås innan datorn dödar processen
```
YEqn.H > for > if > left1 & right................. 83.251             4       
YEqn.H > for > if > left2......................... 58.074             4       
YEqn.H > for > if > left3......................... 77.121             4       
YEqn.H > for > if > solve......................... 116.88             3       
YEqn.H > for > if > solve()....................... 71.33              4        
```

Samma som ovan men för sekventiell kod:
## (41 41 100) ##
```
What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.18481            4       
YEqn.H > for > if > left2......................... 0.23309            4       
YEqn.H > for > if > left3......................... 0.16927            4       
YEqn.H > for > if > solve......................... 1.1186             3       
YEqn.H > for > if > solve()....................... 0.7932             4   
What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.60525            12      
YEqn.H > for > if > left2......................... 0.70461            12      
YEqn.H > for > if > left3......................... 0.50496            12      
YEqn.H > for > if > solve......................... 4.1823             12      
YEqn.H > for > if > solve()....................... 2.3659             12      
```
## (51 51 110) ##
```
What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.30688            4       
YEqn.H > for > if > left2......................... 0.33176            4       
YEqn.H > for > if > left3......................... 0.29096            4       
YEqn.H > for > if > solve......................... 1.9888             3       
YEqn.H > for > if > solve()....................... 1.488              4  
What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 0.93697            12      
YEqn.H > for > if > left2......................... 1.0057             12      
YEqn.H > for > if > left3......................... 0.86806            12      
YEqn.H > for > if > solve......................... 7.2211             12      
YEqn.H > for > if > solve()....................... 4.4087             12      
```
## (82 82 200) ##
Kraschar så jag har inte för Count 12
```
What                                               Time               Count   
YEqn.H > for > if > left1 & right................. 8.572              4       
YEqn.H > for > if > left2......................... 2.2497             4       
YEqn.H > for > if > left3......................... 63.709             4       
YEqn.H > for > if > solve......................... 121.96             3       
YEqn.H > for > if > solve()....................... 89.125             4     
```

Jag får kolla på meshes igen...