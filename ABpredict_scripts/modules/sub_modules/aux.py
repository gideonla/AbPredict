import re
import sys
import pdb
import fileinput
def change_pdb_num(pdb_file,anarci_number_file):
	"""This function takes a pdb file name and an anarci 
	output file (a simple text file, each line has a residue number)
	and renumbers the pdb accordingly. I assume that the order of the chains 
	in the anarci file is light before heavy.
	"""
	anarci_lineList = [line.rstrip('\n') for line in open(anarci_number_file)]
	anarci_first_res_pos=[i for i,x in enumerate(anarci_lineList) if x == '1']
	vl_anarci_residues=anarci_lineList[0:anarci_first_res_pos[1]-1]
	vh_anarci_residues=anarci_lineList[anarci_first_res_pos[1]:]
	#count how many heavy and light chain residues there are in the pdb file
	pdb_file_str= open(pdb_file).read()
	num_of_VH_residues = len(re.findall("ATOM.*N .*[A-Z] H", pdb_file_str, re.IGNORECASE))
	num_of_VL_residues = len(re.findall("ATOM.*N .*[A-Z] L", pdb_file_str, re.IGNORECASE))
	# For some reason the anarci script starts it's numbering at '1' which is not always there in the PDB 
	# So in the following section I pop the un-needed residues from the anarci lists
	vl_delta=len(vl_anarci_residues)-num_of_VL_residues
	vh_delta=len(vh_anarci_residues)-num_of_VH_residues
	del vh_anarci_residues[:vh_delta] 
	del vl_anarci_residues[:vl_delta]
	combined_anarci_lists=vl_anarci_residues+vh_anarci_residues
	#change the actual pdb file
	num=-1
	#for line in (open(pdb_file)):
	for line in fileinput.input(pdb_file,inplace=True, backup='.bak'):
		line=line.rstrip("\n\r")
		if re.search("ATOM.*N .*[A-Z] [L/H]", line):			
			num+=1
		if line[0:4]=="ATOM":		
				print (line[0:22],end = '')
				if re.search('[a-zA-Z]', combined_anarci_lists[num]):
					print (combined_anarci_lists[num].rjust(5),end = '')
					#print (first_half)
				else:
					print (combined_anarci_lists[num].rjust(4)+" ",end = '')
				print (line[27:])
		else:
			print (line)
			
		


	#pdb.set_trace()


if __name__ == "__main__":
	change_pdb_num(sys.argv[1],sys.argv[2])

