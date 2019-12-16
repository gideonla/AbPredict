#!/bin/sh

#renumber a pdb file from 1->nres, ignoring chainID, TER statements and anything else taht doesn't have an ATOM card.

pdb=$1
res_number_file=$2
name=( $(awk '{print $1}' $res_number_file) )
echo "${name[@]}"
 

awk  -v var="${name[*]}" 'BEGIN {num=1} 	
	{split(var,list," ")}
	{
	if( $1=="ATOM" ){
 		if( substr($0,14,3)=="N  " ) ++num;
		printf substr( $0, 1, 22 );
                printf "%4s",list[num]
		print substr( $0,28,10000);
	}
	else print $0
	}' $pdb > $pdb.TMP

#mv $pdb.TMP $pdb

