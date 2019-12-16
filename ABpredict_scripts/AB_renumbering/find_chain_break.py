#!/usr/bin/python
import os
from os.path import basename
import csv


num_of_arg=len(sys.argv)
lines=[]
if os.path.isfile("pdb_file_list") :
	
	pdb_file_list = open("pdb_file_list", "r")
	lines = pdb_file_list.read().split('\n')
	#print lines
else:
	print "target"+str(sys.argv[num_of_arg-1])
	lines.append(sys.argv[num_of_arg-1])


f = open('chain_break_PDBs', 'wa')
for line in lines:
	try:
		if line !="":
			cmd.load(line)
			
			cmd.remove("hetatm")

			target=basename((line))

			target=target.split(".")[0]
			
			cmd.select("target",target)
			print target
			print line, target
	###############################
	# creat selection for target vl
	###############################
			residue={'residues':[]}
			cmd.iterate(target,'residues.append(resi)',space=residue)
			reduced_set=set(residue.get('residues'))
			reduced_set=set(residue.get('residues'))
			reduced_set=list(reduced_set)
			reduced_set = sorted(map(int,reduced_set))
			for i in range(0,len(reduced_set)-1):
				print i
				tmp=cmd.distance( "dist", "i. "+str(reduced_set[i])+" and name CA and "+target,"i. "+str(reduced_set[i+1])+" and name CA and "+target)
				print tmp
				if tmp > 4.5:
					string=target+" "+str(i+1)+" "+str((int(i)+2)) +" "+str(tmp)
					print string
					f.write(string+"\n")
#					os.remove(line) # - deletes file. don't want this option now (23Fab15)
					break
					
			cmd.delete(target)		
					
	except:
		continue

