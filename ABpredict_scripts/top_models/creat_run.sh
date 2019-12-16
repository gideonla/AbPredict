#!/bin/sh

#produces one pair of files, commandline and jobname. Commandline is used to execute jobname through lsf on the fleishman queue.
#USAGE:
cdr_co_rms=$1
cluster_radius=$2

for i in `cat $cdr_co_rms |tr ":" " "|awk '{print $2}'`; do
readlink -f ../pdb/"$i"*;done > pdb_file_list 
sed -i "1d" pdb_file_list
rm command & rm job/* &
mkdir job err out
file_size=$(wc -l pdb_file_list|awk '{print $1}')
for ((i=$file_size; i>=0; i--)); do echo $i;./create_job_pymol.sh ./cluster_loop_conformations.py $cluster_radius $i; done
