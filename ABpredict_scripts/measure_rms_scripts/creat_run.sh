#!/bin/sh

#produces one pair of files, commandline and jobname. Commandline is used to execute jobname through lsf on the fleishman queue.
#USAGE:
rm -rf CDR_CO_rms
pdb=$1
rm command & rm job/* &
mkdir err out pdb job
for i in  ../pdb/*.pdb*; do name=`basename "$i" .pdb*`;./create_job_pymol.sh $name /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/measure_rms_scripts/measure_CDR_rms.py /home/labs/fleishman/gideonla/new_SpliceOutAntibody/2BRR.ppk_ideal.pdb.gz $pdb $i   ; done
