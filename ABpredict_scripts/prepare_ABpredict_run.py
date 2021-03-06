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
	logger.debug(segment_lengths)

	#Prepare
	#pdb.set_trace()

	segment_pdb_entries=get_db_entries(segment_lengths)
	#make flags file object
	flag_object=parse_flags(FLAG_PATH)
	flag_object.add_script_var("sequence",str(light_seq)+str(heavy_seq) )
	flag_object.add_script_var("MC_cyc",str(MC_cyc) )
	flag_object.write_flags_file()
	create_jobs(segment_pdb_entries,flag_object,num_of_jobs)
	job=JobsExecuter("command",wait_interval=0.1)
	job.run_jobs()#this will submit the jobs and wait for the jobs to finish



	# get lowest energy scoring pdb
	pdb_files = os.listdir(cwd+"/pdb/")
	if (len(pdb_files)==0):
		logger.error("Tried to get lowest energy pdb but no models were created. Exiting")
		sys.exit("Tried to get lowest energy pdb but no models were created")

	top_model_energy=sys.float_info.max
	top_pdb = ""
	for pdb_file in pdb_files:
		try:
			f= open(cwd+"/pdb/"+pdb_file, 'rb')
			#print (f)
			for line in f:
				str_line=str(line)
				if re.search("pose", str_line):
					liner=str_line.rstrip("\n\r")
					model_energy = liner.split(" ")[-1:]
					if (float (model_energy[0][:-4])<top_model_energy):
						top_model_energy=float (model_energy[0][:-4])
						top_pdb = cwd+"/pdb/"+pdb_file
			print(top_pdb)
		except:
			print ("not regular file")
		try:
			f=gzip.open(cwd+"/pdb/"+pdb_file, 'rb')
			#print (f)
			for line in f:
				str_line=str(line)
				if re.search("pose", str_line):
					liner=str_line.rstrip()
					model_energy = liner.split(" ")[-1:]
					print(model_energy)
					if (float (model_energy[0][:-4])<top_model_energy):
						top_model_energy=float (model_energy[0][:-4])
					top_pdb = cwd+"/pdb/"+pdb_file
		except:
			print ("not zip file")
	#this is a hard coded number but from expirience if a strcture has a higher R.e.u than -400 it might mean there is a problem with this model
	if (top_model_energy>-400):
		logger.warning("Model energies are high. This might be a problem\n")

	if not os.path.isdir("measure_rms"):
		os.makedirs("measure_rms")

	for filename in glob.glob(os.path.join("/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/measure_rms_scripts", '*')):
		shutil.copy(filename, cwd+"/measure_rms/")
	os.chdir("measure_rms")
	logger.info("Calculating model RMSD")
	creat_run_cmd = './creat_run.sh ' +top_pdb
	com_proc = subprocess.run([creat_run_cmd], shell=True, check=True,stdout=subprocess.PIPE)

	comm = 'sh command'
	com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	job_ids=re.findall('\d+',str(com_proc.stdout))
	num_of_jobs = len(job_ids)
	print ("In measure RMS:",num_of_jobs,"jobs were submitted")
	job_status = ["DONE", "EXIT"]
	count=0
	while count<num_of_jobs:
		for i in job_ids:	
			job_com = subprocess.run(["bjobs "+str(i)], shell=True, check=True,stdout=subprocess.PIPE)
			print (job_com.stdout)
			if any(x in str(job_com.stdout) for x in job_status):
				count = count+1
				job_ids.remove(i)
		print("numer of jobs in measure rms:"+str(count))
	import time
	time.sleep(60) 
	comm = 'sh post_run_command' 
	com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	# Check that measure RMS finished corectly.
	num_lines = sum(1 for line in open('cdr_co_rms'))
	#The num of lines in the "cdr_co_rms" file should be num_of_jobs+1, if not these is a problem and we should fail
	if (num_lines<num_of_jobs):
		logger.error("Failed to calculate RMSD. Exiting")
		sys.exit("Failed to calculate RMSD. Exiting")
		
	os.chdir("../")

# Here I make the top 3 antibody structures that the user will download 
	if not os.path.isdir("top_models"):
		os.makedirs("top_models")

	for filename in glob.glob(os.path.join("/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/top_models", '*')):
		shutil.copy(filename, cwd+"/top_models/")


	models_sorted_path=cwd+"/measure_rms/cdr_co_rms"

	os.chdir(cwd+"/top_models/")

	try:	
		comm = './creat_run.sh '+models_sorted_path+' '+str(clustering_radius)
		print ("executing :	",comm)
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)

	except:
		logger.error("Failed to find top models. Exiting")
		sys.exit("Failed to find top models. Exiting")
	comm = 'sh command'
	com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	job_ids=re.findall('\d+',str(com_proc.stdout))
	num_of_jobs = len(job_ids)
	print (num_of_jobs,"jobs were submitted")
	job_status = ["DONE", "EXIT"]
	count=0
	while count<num_of_jobs:
		for i in job_ids:	
			job_com = subprocess.run(["bjobs "+str(i)], shell=True, check=True,stdout=subprocess.PIPE)
			print (job_com.stdout)
			if any(x in str(job_com.stdout) for x in job_status):
				count = count+1
				job_ids.remove(i)
		print(count)
	time.sleep(60) 
	comm = 'sh post_job_command'
	com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	
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
		f.close()
		if not check_final_model_output(pdb_files,str(light_seq)+str(heavy_seq)):
			logger.error("Final Model sequence does not match input sequence. Exiting")
			sys.exit("Final Model sequence does not match input sequence. Exiting")
			

		

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
			change_pdb_num(str(i)+'.pdb',str(i)+'.tmp')
	except:
		logger.error("Could not change PDB numbering")
		sys.exit("Could not change PDB numbering. Exiting")
	try:
		comm = '/apps/RH7U2/gnu/python/2.7.12/bin/pymol -c /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/label_CDR.py /home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/required_eden_Rosetta_files/2BRR.ppk_ideal.pdb.gz 1.pdb 2.pdb'
		com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
	except:
		logger.error("could not add L\H chain identifiers")
		sys.exit("could not add L\H chain identifiers. Exiting")


	logger.info("Processing completed")

if __name__ == '__main__':
	try:
		main()
	except Exception as e:
		logger.error(e.message)
		sys.exit(e.message)
		










	#print(count)
		
