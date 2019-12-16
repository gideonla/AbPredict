#!/bin/sh

#renumber a pdb file from 1->nres, ignoring chainID, TER statements and anything else taht doesn't have an ATOM card.

for pdb in "$@"
do

awk  'BEGIN {num=0;chain="L"} 	
	{
	if( substr($0,0,3)=="TER" ) chain="H";
	if( $1=="ATOM" ){
 		if( substr($0,14,3)=="N  " ) ++num;
		printf substr( $0, 1, 21 );
		printf "%s",chain
		print substr( $0,23,10000);
	}
	else print $0
	}' $pdb > $pdb.TMP

mv $pdb.TMP $pdb
done
