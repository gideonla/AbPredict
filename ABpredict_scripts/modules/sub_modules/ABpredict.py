import sys
#sys.path.append('/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/modules/sub_modules/lsfe')
sys.path.append('/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/modules/sub_modules/')
from rosetta_job import *
from lsfe.job_executer import *
from pprint import pprint
from sys import argv
import subprocess
import shutil
import os
import glob
import sys
import re
import gzip
import argparse
import logging
import warnings
import pdb
import random


from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO
    from Bio import SeqIO


###constant definitions######
DB_files={"L1_L2":"/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_db_files/L1_L2.db","L3":"/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_db_files/L3.db","H1_H2":"/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_db_files/H1_H2.db","H3":"/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/AB_db_files/H3.db"}
FLAG_PATH="/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/required_eden_Rosetta_files/flags"
RUN_PATH=os.getcwd()
ROSETTA_OBJ=RosettaPaths("/home/labs/fleishman/WebServer/ABPREDICT/bin","/home/labs/fleishman/WebServer/ABPREDICT/database",1)


logger = logging.getLogger(__name__)


def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file



def calculate_segment_lengths(chains):
	heavy_str=str(chains.heavy_seq)
	heavy_cc_str=str(chains.heavy_cc['seq'])
	light_str=str(chains.light_seq)
	light_cc_str=str(chains.light_cc['seq'])
	H1_H2_len =  (len(heavy_cc_str)*4)-4
	H3_len=len(heavy_str)-(heavy_str.rfind(heavy_cc_str)+len(heavy_cc_str))
	L3_len=len(light_str)-(light_str.rfind(light_cc_str)+len(light_cc_str))
	L1_L2_len =  (len(light_cc_str)*4)-4
	L3_len = L3_len*4-4
	H3_len = H3_len*4-4
	logger.info("Segment lengths: {},{},{},{}, L1_L2,L3, H1_H2, and H3 respectively ".format(L1_L2_len,L3_len,H1_H2_len,H3_len))
	return {"L1_L2":L1_L2_len,"L3":L3_len,"H1_H2":H1_H2_len,"H3":H3_len}

def get_db_entries(segment_lengths):
	segment_pdb_entries={}
	for segment, path in DB_files.items():
		pdb_names=[]
		try:
			with open(path, 'r') as filename:
				logger.debug("opened file: {0}".format(filename))
				for line in filename:
				    	if (len(line.split()) == segment_lengths[segment]) :
                            			pdb_names.append(line.split()[-1])
			segment_pdb_entries.update({segment:pdb_names})
		except IOError:
			logger.error("file not found: {}".format(filename))
		logger.debug("The selected pdbs for the given segment :{}".format(segment))
		logger.debug(segment_pdb_entries[segment])
		if len(segment_pdb_entries[segment])==0:
			logger.error("No PDB's of the appropriate length were found for {}, exiting".format(segment))
			sys.exit("No PDB's of the appropriate length were found for {}, exiting".format(segment))
	return segment_pdb_entries

def create_jobs(segment_pdb_entries,flag_object,num_of_jobs=500):
	for i in range(0,num_of_jobs):
		L1_L2=(random.choice(segment_pdb_entries["L1_L2"]))
		L3=(random.choice(segment_pdb_entries["L3"]))
		H1_H2=(random.choice(segment_pdb_entries["H1_H2"]))
		H3=(random.choice(segment_pdb_entries["H3"]))
		flag_object.add_script_var("entry_H1_H2", H1_H2)
		flag_object.add_script_var("entry_L1_L2", L1_L2)
		flag_object.add_script_var("entry_H3", H3)
		flag_object.add_script_var("entry_L3", L3) 
		job_name=L1_L2+"_"+L3+"_"+H1_H2+"_"+H3+"_"
		flag_object.add_option("out:prefix",job_name)
		logger.debug("job names:")
		logger.debug(job_name)
		flag_object
		job=RosettaJob(job_name,RUN_PATH,flag_object,ROSETTA_OBJ)
		job.create_job()
	
def make_flags():
	parse_flags()


#the main module accepts the fasta file
def main():
	parser = argparse.ArgumentParser(description='Runs ABpredict automatically.\nInput is an antibody fasta file',prog='ABpredict.py', usage='%(prog)s <fasta file>')
	parser.add_argument('fasta_file')
	args = parser.parse_args()
	filename=args.fasta_file
#####Check file exists###########
	try:
		with open(filename, 'r') as input_file:
		    logger.info("opened file: {0}".format(filename))
		    logger.debug("file opened in read mode")
	except IOError:
		logger.error("file not found: {}".format(filename))
#####Check file is fasta format###########
	is_fasta(args.fasta_file)
	try:
		if (not is_fasta(args.fasta_file)):
			sys.exit('Input is not a fasta file')
	except SystemExit as e:
	    # this log will include traceback
		logger.exception('Input is not a fasta file, exiting')
	#pdb.set_trace()
	chains = extract_chains(filename)
	segment_lengths=calculate_segment_lengths(chains)
	segment_pdb_entries=get_db_entries(segment_lengths)
	#make flags file object
	flag_object=parse_flags(FLAG_PATH)
	flag_object.add_script_var("sequence",str(chains.light_seq)+str(chains.heavy_seq) )
	flag_object.write_flags_file()
	create_jobs(segment_pdb_entries,flag_object)
	job=JobsExecuter("command")
	job.run_jobs()#this will submit the jobs and wait for the jobs to finish



if __name__ == '__main__':
    main()



