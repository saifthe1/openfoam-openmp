#!/bin/sh

function ninjaProf()
{
	DIRECTORY="KandidatTestRun"
	if [ ! -d $DIRECTORY ]
	then
		mkdir $DIRECTORY
	fi

	i=1
	
	while [ -d $DIRECTORY/$i ]
	do
		let i=$i+1
	done
	
	DIRECTORY=$DIRECTORY/$i
	
	mkdir $DIRECTORY	

	$1 >& $DIRECTORY/log.$1
	mv gmon.out $DIRECTORY/gmon.$1
	gprof $FOAM_APPBIN/$1 $DIRECTORY/gmon.$1 | head -n 18 > $DIRECTORY/result.$1
	gprof -q $FOAM_APPBIN/$1 $DIRECTORY/gmon.$1 > $DIRECTORY/callgraph.$1
	cat $DIRECTORY/callgraph.$1 | \
		perl -e "while(<>) { while(/<[^<>]+>/) { s/<[^<>]+>//g } s/Foam:://g; print $_ }" \
		> $DIRECTORY/callgraph_clean.$1
}

