
dir=`pwd`
date=`date +%s%N | awk '{print substr($0,9,length($0)-11)}'`
jobname=${dir}/job/job."${1:0:100}"
commandname=${dir}/command
out=${dir}/out/out."${1:0:100}"
err=${dir}/err/err."${1:0:100}"

echo "#!/bin/bash" > $jobname
echo ". /usr/share/lsf/conf/profile.lsf" >> $jobname
echo ". /apps/RH7U2/Modules/default/init/bash" >> $jobname
echo "module load gcc" >> $jobname



echo "cd $dir">> $jobname
printf "/home/labs/fleishman/WebServer/ABPREDICT/bin/rosetta_scripts.static.linuxgccrelease " >> $jobname

for var in "${@:2}"
do
    echo -n " $var " >> $jobname 
done

echo "" >> $jobname

echo "bsub  -C 2048 -u /dev/null -R rusage[mem=2048] -L /bin/bash -G fleishman-wx-grp-lsf -q fleishman -o $out -e $err $jobname" >> $commandname
chmod +x $jobname

