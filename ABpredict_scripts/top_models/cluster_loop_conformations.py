#!/usr/bin/python
import os
from os.path import basename
import re
import sys

#This script will align the pdbs to a template vh. tempalte vh is 2BRR.ppk.pdb

#file list
lines=[]
num_of_arg=len(sys.argv)
#print num_of_arg
if os.path.isfile("pdb_file_list") :
	pdb_file_list = open("pdb_file_list", "r")
	lines = pdb_file_list.read().split('\n')


rms_val = float(sys.argv[-2])
file_number = int(sys.argv[-1])
print rms_val
cluster_num=1
i=file_number-1

print "now loading "+lines[i]
print "i is :"+str(i)
cmd.load(lines[i])
print i
pdbs=cmd.get_object_list("all")
current_structure=pdbs[0]
last_current_res = {'resnums': []}
cmd.iterate(current_structure, 'resnums.append(resi)',space=last_current_res)
last_current_res=last_current_res.get('resnums')[-1]
print last_current_res
#	if not os.path.isdir("cluster_"+str(cluster_num)):
#	   	os.makedirs("cluster_"+str(cluster_num))
#	cmd.save("cluster_"+str(cluster_num)+"/"+current_structure+".pdb",current_structure)
j=i+1
f = open(str(file_number)+"_"+current_structure, 'w')
f.write(str(file_number)+" "+current_structure+"\n")
while j<len(lines):
	#print "now loading "+lines[j]
	if lines[j]=="":
		j=j+1
		continue
	cmd.load(lines[j])
	
	pdbs=cmd.get_object_list("all")
	target_structure=pdbs[1]
	print "target structure: "+target_structure
	last_target_res = {'resnums': []}
	cmd.iterate(target_structure, 'resnums.append(resi)',space=last_target_res)
	last_target_res=last_target_res.get('resnums')[-1]
	if last_current_res!=last_target_res:
		j=j+1
		cmd.delete(target_structure)
		continue
	pair_fit_val = cmd.pair_fit(current_structure+" and (name C+o)",target_structure+" and (name C+o)")
	if pair_fit_val<rms_val:
		print "rms val is %f and pair_fit val is %f"%(rms_val,pair_fit_val)
#		cmd.save("cluster_"+str(cluster_num)+"/"+target_structure+".pdb",target_structure)
		f = open(str(file_number)+"_"+current_structure, 'a')
		f.write(str(j+1)+" "+target_structure+"\n")
		lines[j] = ""
	cmd.delete(target_structure)
	j=j+1
cmd.delete(current_structure)
cluster_num=cluster_num+1
i=i+1
			
			
		
		



