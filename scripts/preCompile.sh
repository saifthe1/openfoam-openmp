#!/bin/sh

function preCompile() 
{
		SOURCE=$1
		echo "LIB_SRC=$WM_PROJECT_DIR/src" > tmp.source
		. tmp.source

		cat Make/options | perl -e "my \$a=''; while(<>) {\$a.=\$_} \$a=~s/\\\\\n/ /g; \$a=~s/ +/ /g; \$a=~s/\(|\)|//g; \$a=~s/ = /=\"/g; \$a=~s/\n\s*\n/\"\n/g; \$a=~s/\s+$//; print \$a. \"\\\"\n\"" > tmp.source
		. tmp.source
		
		gcc	-E $EXE_INC $SOURCE

}
