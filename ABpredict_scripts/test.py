from ChainExtract import ChainExtract
from pprint import pprint
from sys import argv
import subprocess
import shutil
import os
import glob

import re
import gzip
import pdb
from modules.modules import *
from modules.sub_modules.aux import *
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


	

	os.chdir(cwd+"/top_models/")


	try:
		for i in (range(1,4)):
			comm = '/apps/RH7U2/gnu/python/2.7.12/bin/python2 /home/labs/fleishman/gideonla/scripts/pdbTools_0.2.1/pdb_seq.py '+str(i)+'.pdb >'+str(i)+'.fa'
			com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	except:
		logger.error("Could not get sequence from PDB")
		sys.exit("Could not get sequence from PDB. Exiting")

	try:
		comm = 'sh /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_renumbering/creat_numer_file.sh &>anarci.out'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	except:
		logger.error("Could not create Chothia numbering file")
		sys.exit("Could not create Chothia numbering file. Exiting")

	try:
		for i in (range(1,4)):
			#import pdb
			#pdb.set_trace()
			change_pdb_num(str(i)+'.pdb',str(i)+'.tmp')
	except:
		logger.error("Could not change PDB numbering")
		sys.exit("Could not change PDB numbering. Exiting")


	logger.info("Processing completed")

if __name__ == '__main__':
    main()











	#print(count)
		
