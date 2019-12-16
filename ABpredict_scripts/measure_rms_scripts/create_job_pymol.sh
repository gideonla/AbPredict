
dir=`pwd`
date=`date +%s%N | awk '{print substr($0,9,length($0)-11)}'`
name=`basename $1 .pdb`
jobname=${dir}/job/job."$name"
commandname=${dir}/command
out=${dir}/out/out."$name"
err=${dir}/err/err."$name"

echo "#!/bin/bash" > $jobname
echo ". /usr/share/lsf/conf/profile.lsf" >> $jobname
echo ". /apps/RH7U2/Modules/default/init/bash"  >> $jobname
echo "module load pymol/1.8.4.0" >> $jobname

echo "cd $dir">> $jobname
echo -n "pymol -c ">> $jobname

for var in "${@:2}"
do
    echo -n " $var " >> $jobname 
done

echo "" >> $jobname

echo "bsub  -C 1024 -u /dev/null -R rusage[mem=1024] -L /bin/bash -G fleishman-wx-grp-lsf -q new-short -o $out -e $err /apps/RH6U4/blcr/0.8.5/bin/cr_run $jobname" >> $commandname
chmod +x $jobname

