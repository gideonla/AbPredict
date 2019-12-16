from ChainExtract import ChainExtract
from pprint import pprint
from sys import argv
import subprocess
import shutil
import os
import glob
import time
import re
import gzip
import pdb
from modules.modules import *
from modules.sub_modules.ABpredict import *
from modules.sub_modules.pdb_seq import *
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import warnings
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO
    from Bio import SeqIO
import logging

#set up logger
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG,filename='log.ABpredict',filemode='w')
logger = logging.getLogger(__name__)

#These are the sequences for the begining and trailing residues for the light and heavy variable domains. I use these for padding if the input sequence is too short
VH_starting_AAs="VQL"
VH_ending_AAs="VTVS"
VL_starting_AAs="IEMT"
VL_ending_AAs="GTKL"


def check_final_model_output(pdb_model,input_sequence):
	model_seq=pdbSeq2Fasta(pdb_model).split(" ")[3].replace("\n","")
	import re
	re.sub(r'\W+', '', model_seq)
	re.sub(r'\W+', '', input_sequence)
	print (model_seq)
	print (input_sequence)
	return model_seq==input_sequence



def main():
	parser = argparse.ArgumentParser(description='Runs ABpredict automatically.\nInput is an antibody fasta file',prog='ABpredict.py', usage='%(prog)s <fasta file>')
	parser.add_argument('fasta_file')
	parser.add_argument('--num_of_jobs', nargs='?', type=int,default=500,help='how many PDB models to generate (default: 500)')
	parser.add_argument('--MC', nargs='?', type=int,default=50,help='how many MC cycles (default: 50)')
	parser.add_argument('--clustering_radius', nargs='?', type=str,default=1,help='clustering radius of final models (default: 1)')
	args = parser.parse_args()
	fasta_file=args.fasta_file
	num_of_jobs=args.num_of_jobs
	MC_cyc=args.MC
	clustering_radius=args.clustering_radius
	
	logger.info("Starting run")
	cwd = os.getcwd()


	# setting up AbChains object
	chains = ChainExtract() #constructor
	

	if is_fasta(fasta_file):
		chains.parse_fasta(fasta_file) # pass fasta sequence to chains object
	else:
		logger.error("Input file is not fasta file. Exiting")
		sys.exit("Input file is not fasta file. Exiting")

	chains.parse_fasta(fasta_file)

	# Getting sequences using hmms
	logger.info("Extracting sequence between disulfides")
	heavy_seq_CC,light_seq_CC = extract_seqeunces(chains,"/home/labs/fleishman/gideonla/ABpredict/test_cut_length/ChainExtract/hmms/all_cys.hmm")

	heavy_seq_CW,light_seq_CF = extract_seqeunces(chains,"/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/hmmscan/kappa_C_F_heavy_smart_domains_aln_W.hmm")

	heavy_seq,light_seq = extract_seqeunces(chains,"/home/labs/fleishman/gideonla/ABpredict/test_cut_length/ChainExtract/database/igg-chains_cdhit.hmm")

	#I used Jakes Hmm to find the VL/VH domains. There is an error with the VH domain length. It's 1 aa too long
	# So I need to remove it here
	VH_seq_delta=heavy_seq.find(heavy_seq_CC)-20
	if (VH_seq_delta>0):
		heavy_seq=heavy_seq[1:]
		


	#If the input sequence is too short in the tail regions then we need to pad it 
	########################HEAVY CHAIN########################
	heavy_seq=check_variable_sequecne(heavy_seq,heavy_seq_CC, heavy_seq_CW, 21,9,VH_starting_AAs,VH_ending_AAs)
	########################LIGHT CHAIN########################
	light_seq=check_variable_sequecne(light_seq,light_seq_CC, light_seq_CF, 22,6,VL_starting_AAs,VL_ending_AAs)

	

	#calculate lengths of AB segments
	H3_len=len(heavy_seq)-len(heavy_seq_CC)-20
	L3_len=len(light_seq)-len(light_seq_CC)-21
	H1_H2_len =  (len(heavy_seq_CC)*4)-4 #The minus 4 and times 4 is to calculate the dihedral db length
	L1_L2_len =  (len(light_seq_CC)*4)-4 #The minus 4 and times 4 is to calculate the dihedral db length
	L3_len=(L3_len+1)*4
	H3_len=(H3_len+1)*4
	segment_lengths={"L1_L2":L1_L2_len,"L3":L3_len,"H1_H2":H1_H2_len,"H3":H3_len}
	logger.info(segment_lengths)


# Here I make the top 3 antibody structures that the user will download 
	if not os.path.isdir("top_models"):
		os.makedirs("top_models")
	os.chdir("top_models")

	logger.info("Calculating top models")
	f= open("top_models", 'rb')
	num_lines = sum(1 for line in open('top_models'))
	#The num of lines in the "top_models" should be exactly 3
	if (num_lines<3):
		logger.error("Failed to find top 3 models. Exiting")
		sys.exit("Failed to find top 3 models. Exiting")
	i=0
	for line in f:
		i=i+1
		for file in glob.glob("../pdb/"+str(line,"utf-8").rstrip()+"*"):
			shutil.copy(file, ".")
			file=file.split("/")[-1]
			print (file)
			os.rename(file, str(i)+".pdb")
#Check sequence of top 3 models and compare it to input sequnce. If there is a mismatch then fail. 
# We saw in the past that if the input seqeunce is too short and AbPredict does not pick up on this then the seuqeunce of the model
# does not match to the input
	for filename in glob.iglob('*.pdb'):
		f = open(filename,'r',encoding='utf-8', errors='ignore')
		pdb_files = f.readlines()
		#print (pdb_files)
		f.close()
		#if not check_final_model_output(pdb_files,str(light_seq)+str(heavy_seq)):
		#	logger.error("Final Model sequence does not match input sequence. Exiting")
		#	sys.exit("Final Model sequence does not match input sequence. Exiting")
			

		

	os.chdir(cwd)
	logger.info("Generating graphs")
	if not os.path.isdir("plots"):
		os.makedirs("plots")

	for filename in glob.glob(os.path.join("/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/plot", '*')):
		shutil.copy(filename, cwd+"/plots/")

	os.chdir("plots/")
	try:
		comm = '/apps/RH7U2/gnu/python/2.7.12/bin/python2 plot_rms_plots.py ../measure_rms/cdr_co_rms ../top_models/top_models'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	except:
		logger.error("Could not generate graphs. Exiting")
		sys.exit("Could not generate graphs. Exiting")

#In the next block of code I am changing the final pdbs. Seperating VL and VH and renumbering according to IMGT numbering
	os.chdir(cwd+"/top_models/")
	with open('pdb_file_list', 'a') as the_file:
    		the_file.write('1.pdb\n')
    		the_file.write('2.pdb\n')
    		the_file.write('3.pdb\n')
	try:
		comm = '/apps/RH7U2/gnu/python/2.7.12/bin/pymol -c /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_renumbering/find_chain_break.py'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	except:
		logger.error("Could not run \"find chain break\"")
		sys.exit("Could not run \"find chain break\". Exiting")
	try:
		comm = 'sh /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_renumbering/add_vl_vh_break.sh'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	except:
		logger.error("Could not add chain break to pdbs")
		sys.exit("Could not add chain break to pdbs. Exiting")

	try:
		comm = 'sh /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_renumbering/renumber_W_chain.sh *.pdb'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	except:
		logger.error("could not add L\H chain identifiers")
		sys.exit("could not add L\H chain identifiers. Exiting")
	try:
		comm = '/apps/RH7U2/gnu/python/2.7.12/bin/python2 /home/labs/fleishman/gideonla/scripts/pdbTools_0.2.1/pdb_seq.py 1.pdb >1.fa'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
		comm = '/apps/RH7U2/gnu/python/2.7.12/bin/python2 /home/labs/fleishman/gideonla/scripts/pdbTools_0.2.1/pdb_seq.py 2.pdb >2.fa'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
		comm = '/apps/RH7U2/gnu/python/2.7.12/bin/python2 /home/labs/fleishman/gideonla/scripts/pdbTools_0.2.1/pdb_seq.py 3.pdb >3.fa'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)

	except:
		logger.error("Could not get sequence from PDB")
		sys.exit("Could not get sequence from PDB. Exiting")

	try:
		comm = 'sh /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_renumbering/creat_numer_file.sh &>anarci.out'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	except:
		logger.error("Could not creat Chothia numbering file")
		sys.exit("Could not creat Chothia numbering file. Exiting")

	try:
		comm = 'sh /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_renumbering/change_pdb_num.sh 1.pdb 1.tmp'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
		comm = 'sh /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_renumbering/change_pdb_num.sh 2.pdb 2.tmp'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
		comm = 'sh /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_renumbering/change_pdb_num.sh 3.pdb 3.tmp'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)

	except:
		logger.error("Could not change PDB numbering")
		sys.exit("Could not change PDB numbering. Exiting")


	logger.info("Processing completed")

if __name__ == '__main__':
    main()











	#print(count)
		
