# Operatorer #

Många av klasserna har operatorer överlagrade, och den enda klassen vars överlagrade operatorer jag finner värd att nämna är tmp.

```
//	Dereference operator.
T & 	operator() ()
//	Const dereference operator.
const T & 	operator() () const
//	Const cast to the underlying type reference.
	operator const T & () const
//	Return object pointer.
T * 	operator-> ()
//	Return const object pointer.
const T * 	operator-> () const
//	Assignment operator. 
void 	operator= (const tmp< T > &)
```

Se upp för (), då detta kan se ut som ett funktionsanrop, men är en operator som returnerar T-objektet.

I övrigt finns tonvis med operatorer överlagrade för hela namespacet. Deklarationer hittas i [fvMatrixs header](http://foam.sourceforge.net/doc/Doxygen/html/fvMatrix_8H_source.html) och definitioner [fvMatrixs källkod](http://foam.sourceforge.net/doc/Doxygen/html/fvMatrix_8C_source.html).