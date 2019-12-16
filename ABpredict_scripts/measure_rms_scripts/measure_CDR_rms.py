#!/usr/bin/python
import os
from os.path import basename
import re
import sys
import zipfile


#This script will align the pdbs to a template vh. tempalte vh is 2BRR.ppk.pdb

#file list
lines=[]
num_of_arg=len(sys.argv)
#print num_of_arg
if os.path.isfile("pdb_file_list") :
	pdb_file_list = open("pdb_file_list", "r")
	lines = pdb_file_list.read().split('\n')
else:
	print(("template: "+str(sys.argv[num_of_arg-3])	)) #this will be 2BRR
	print(("target: "+str(sys.argv[num_of_arg-2]))) # this will be the target to be modeled. e.g 1AHW
	print(("model: "+str(sys.argv[num_of_arg-1]))) # this will be the target to be modeled. e.g 1AHW

lines.append(sys.argv[num_of_arg-3])
lines.append(sys.argv[num_of_arg-2])
lines.append(sys.argv[num_of_arg-1])

#load refrence structure
cmd.load(lines[0],"template")
cmd.load(lines[1],"target")
cmd.load(lines[2])
pdbs=cmd.get_object_list("all")
template=pdbs[0]
target=pdbs[1]
model=pdbs[2]


cmd.remove("het")

#cmd.select("template",template)
#cmd.select("target",target)
cmd.select("models",model)
	
	
cmd.select("template and resi 1-98",template+" and resi 1-98")
cmd.cealign(target+" and resi 1-98","template and resi 1-98") #align design to template

	###############################
	# select VL CDRs
	###############################

cmd.select("target_L1_first"," (br. (template and resi 22) around 1) and resn cys and target")
cmd.select("target_L1_last"," (br. (template and resi 35) around 1) and resn trp and target")
cmd.select("target_L2_first"," (br. (template and resi 48) around 0.8) and target")
cmd.select("target_L2_last"," (br. (template and resi 53 and name ca) around 1.3) and target")
cmd.select("target_L3_first"," (br. (template and resi 88) around 0.6) and target and resn cys")
cmd.select("target_L3_last"," (br. (template and resi 98) around 0.6) and target")


target_L1_first = {'resnums': []}
cmd.iterate("target_L1_first", 'resnums.append(resi)',space=target_L1_first)
#print (target_L1_first.get('resnums')[0])
target_L1_first=target_L1_first.get('resnums')[0]

target_L1_last = {'resnums': []}
cmd.iterate("target_L1_last", 'resnums.append(resi)',space=target_L1_last)
#print (target_L1_last.get('resnums')[0])
target_L1_last=target_L1_last.get('resnums')[0]


target_L2_first = {'resnums': []}
cmd.iterate("target_L2_first", 'resnums.append(resi)',space=target_L2_first)
target_L2_first=target_L2_first.get('resnums')[0]
print (target_L2_first)

target_L2_last = {'resnums': []}
cmd.iterate("target_L2_last", 'resnums.append(resi)',space=target_L2_last)
target_L2_last=target_L2_last.get('resnums')[0]
print (target_L2_last)


target_L3_first = {'resnums': []}
cmd.iterate("target_L3_first", 'resnums.append(resi)',space=target_L3_first)
target_L3_first=target_L3_first.get('resnums')[0]
print (target_L3_first)

target_L3_last = {'resnums': []}
cmd.iterate("target_L3_last", 'resnums.append(resi)',space=target_L3_last)
target_L3_last=target_L3_last.get('resnums')[0]
print (target_L3_last)

		
cmd.select("target_L1_CO","resi %s-%s and target and name C+O"%(target_L1_first,target_L1_last))
cmd.select("model_L1_CO","resi %s-%s and models and name C+O"%(target_L1_first,target_L1_last))
L1_pair_fit_val = cmd.pair_fit("target_L1_CO","model_L1_CO")

cmd.select("target_L2_CO","resi %s-%s and target and name C+O"%(target_L2_first,target_L2_last))
cmd.select("model_L2_CO","resi %s-%s and models and name C+O"%(target_L2_first,target_L2_last))
L2_pair_fit_val = cmd.pair_fit("target_L2_CO","model_L2_CO")

cmd.select("target_L3_CO","resi %s-%s and target and name C+O"%(target_L3_first,target_L3_last))
cmd.select("model_L3_CO","resi %s-%s and models and name C+O"%(target_L3_first,target_L3_last))
L3_pair_fit_val = cmd.pair_fit("target_L3_CO","model_L3_CO")

	###############################
	# select VH CDRs
	###############################
cmd.cealign(template+" and resi 120-","target and resi 120-190")
cmd.select("target_H1_first"," (br. (template and resi 126 and name CA) around 1) and resn cys and target")
cmd.select("target_H1_last"," (br. (template and resi 140 and name CA) around 1.5) and resn trp and target")
cmd.select("target_H2_first"," (br. (template and resi 151 and name CA) around 1.5) and target")
cmd.select("target_H2_last"," (br. (template and resi 164 and name ca) around 1.0) and target")
cmd.select("target_H3_first"," (br. (template and resi 200 and name CA) around 1.1) and target and resn cys")
cmd.select("target_H3_last"," (br. (template and resi 217 and name CA) around 1.5) and target and (resn trp+phe+leu)")


target_H1_first = {'resnums': []}
cmd.iterate("target_H1_first", 'resnums.append(resi)',space=target_H1_first)
#print (target_H1_first.get('resnums')[0])
target_H1_first=target_H1_first.get('resnums')[0]

target_H1_last = {'resnums': []}
cmd.iterate("target_H1_last", 'resnums.append(resi)',space=target_H1_last)
#print (target_H1_last.get('resnums')[0])
target_H1_last=target_H1_last.get('resnums')[0]


target_H2_first = {'resnums': []}
cmd.iterate("target_H2_first", 'resnums.append(resi)',space=target_H2_first)
target_H2_first=target_H2_first.get('resnums')[0]
print (target_H2_first)

target_H2_last = {'resnums': []}
cmd.iterate("target_H2_last", 'resnums.append(resi)',space=target_H2_last)
target_H2_last=target_H2_last.get('resnums')[0]
print (target_H2_last)


target_H3_first = {'resnums': []}
cmd.iterate("target_H3_first", 'resnums.append(resi)',space=target_H3_first)
target_H3_first=target_H3_first.get('resnums')[0]
print (target_H3_first)

target_H3_last = {'resnums': []}
cmd.iterate("target_H3_last", 'resnums.append(resi)',space=target_H3_last)
target_H3_last=target_H3_last.get('resnums')[0]
print (target_H3_last)

		
cmd.select("target_H1_CO","resi %s-%s and target and name C+O"%(target_H1_first,target_H1_last))
cmd.select("model_H1_CO","resi %s-%s and models and name C+O"%(target_H1_first,target_H1_last))
H1_pair_fit_val = cmd.pair_fit("target_H1_CO","model_H1_CO")

cmd.select("target_H2_CO","resi %s-%s and target and name C+O"%(target_H2_first,target_H2_last))
cmd.select("model_H2_CO","resi %s-%s and models and name C+O"%(target_H2_first,target_H2_last))
H2_pair_fit_val = cmd.pair_fit("target_H2_CO","model_H2_CO")

cmd.select("target_H3_CO","resi %s-%s and target and name C+O"%(target_H3_first,target_H3_last))
cmd.select("model_H3_CO","resi %s-%s and models and name C+O"%(target_H3_first,target_H3_last))
print(("target_H3_CO","resi %s-%s and target and name C+O"%(target_H3_first,target_H3_last)))
print(("model_H3_CO","resi %s-%s and models and name C+O"%(target_H3_first,target_H3_last)))
H3_pair_fit_val = cmd.pair_fit("target_H3_CO","model_H3_CO")


	###############################
	# measure whole VL/VH CO RMS
	###############################

#target_last_res={'resnums': []}
#cmd.iterate("target", 'resnums.append(resi)',space=target_last_res)
#target_last_res=target_last_res.get('resnums')[-1:]
whole_pair_fit_val = cmd.pair_fit("target and name C+O","models and name C+O")

target_dih_ang=cmd.get_dihedral("target_H1_first and name ca","target_H1_first and name o","target_L1_first and name ca","target_L1_first and name o")
model_dih_ang=cmd.get_dihedral("models and resi %s and name ca"%(target_H1_first),"models and resi %s and name o"%(target_H1_first),"models and resi %s and name ca"%(target_L1_first),"models and resi %s and name o"%(target_L1_first))

#get the energy from the pdb model
import gzip
try:
	f=open(sys.argv[num_of_arg-1], "r")
	for line in f:
		if re.search("pose", line):
			liner=line.rstrip()
			model_energy = liner.split(" ")[-1:]
			print(model_energy)
except:
	print ("not regular file")
try:
	f=gzip.open(sys.argv[num_of_arg-1],'rb')
	for line in f:
		if re.search("pose", line):
			liner=line.rstrip()
			model_energy = liner.split(" ")[-1:]
			print(model_energy)
except:
	print ("not zip file")


if not os.path.isdir("CDR_CO_rms"):
	try:   	
		os.makedirs("CDR_CO_rms")
	except OSError as err:		
		print ("oops, folder already exists")

f = open("CDR_CO_rms/"+model+"_CDR_CO_rms", 'w')
f.write("Model_name Energy whole_vl_vh  L1 L2 L3 H1 H2 H3 \n")
f.write("%s %s %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n"%(model,model_energy[0],whole_pair_fit_val,L1_pair_fit_val,L2_pair_fit_val,L3_pair_fit_val,H1_pair_fit_val,H2_pair_fit_val,H3_pair_fit_val,model_dih_ang))

f = open("../jsmol_cdr_label", 'w')
L1=((int(target_L1_first)+int(target_L1_last))/2)
f.write("L1:"+str(L1)+"\n")
L2=((int(target_L2_first)+int(target_L2_last))/2)
f.write("L2:"+str(L2)+"\n")
L3=((int(target_L3_first)+int(target_L3_last))/2)
f.write("L3:"+str(L3)+"\n")
H1=((int(target_H1_first)+int(target_H1_last))/2)
f.write("H1:"+str(H1)+"\n")
H2=((int(target_H2_first)+int(target_H2_last))/2)
f.write("H2:"+str(H2)+"\n")
H3=((int(target_H3_first)+int(target_H3_last))/2)
f.write("H3:"+str(H3)+"\n")


