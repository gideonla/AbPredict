#!/bin/sh
 export PATH=/apps/RH7U2/gnu/hmmer/3.1b2/bin:$PATH
for i in `seq 1 3`; do
	/home/labs/fleishman/gideonla/.local/bin/ANARCI -i $i.fa --scheme chothia |cut -c1-10 |grep "^L" >$i.tmp  2>$i.error
	/home/labs/fleishman/gideonla/.local/bin/ANARCI -i $i.fa --scheme chothia |cut -c1-10 |grep "^H" >>$i.tmp
	cat $i.tmp|awk '{print $2$3}' >tmp
	mv tmp $i.tmp
done
	
