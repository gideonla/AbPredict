cat chain_break_PDBs| while read i ; do echo $i; pdbID=`echo $i|awk '{print $1}'`; resID=`echo $i|awk '{print $3}'`; sed -i "/ N.* [A-Z] $resID/i \TER" $pdbID.*; done
