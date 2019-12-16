from ChainExtract import ChainExtract
from pprint import pprint
from sys import argv
import subprocess
import shutil
import os
import glob
import sys
import re
import gzip
import pdb
from modules import *
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

#These are the pssm scoring begining and trailing residues for the light and heavy variable domains. I use these for padding if the input sequence is too short
VH_starting_AAs="EVQL"
VH_ending_AAs="VTVS"
VL_starting_AAs="IEMT"
VL_ending_AAs="GTKL"


logger.info("Starting run")
cwd = os.getcwd()

for filename in glob.glob(os.path.join("/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/required_eden_Rosetta_files/", '*')):
	shutil.copy(filename, os.getcwd())


# setting up AbChains object
chains = ChainExtract() #constructor
fasta_file=argv[1]

if is_fasta(fasta_file):
	chains.parse_fasta(fasta_file) # pass fasta sequence to chains object
else:
	logger.error("Input file is not fasta file. Exiting")
	sys.exit("Input file is not fasta file. Exiting")

chains.parse_fasta(fasta_file)

# Getting sequencees using hmms
logger.info("Extracting sequence between disulfides")
heavy_seq_CC,light_seq_CC = extract_seqeunces(chains,"/home/labs/fleishman/gideonla/ABpredict/test_cut_length/ChainExtract/hmms/all_cys.hmm")

heavy_seq_CW,light_seq_CF = extract_seqeunces(chains,"/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/hmmscan/kappa_C_F_heavy_smart_domains_aln_W.hmm")

heavy_seq,light_seq = extract_seqeunces(chains,"/home/labs/fleishman/gideonla/ABpredict/test_cut_length/ChainExtract/database/igg-chains_cdhit.hmm")




#If the input sequence is too short in the tail regions then we need to pad it 
########################HEAVY CHAIN########################
first_heavy_cys_position=(heavy_seq.find(heavy_seq_CC))
W_heavy_position=(heavy_seq.find(heavy_seq_CW)+len(heavy_seq_CW)-1)
H3_tail_length=len(heavy_seq)-W_heavy_position-1
if (first_heavy_cys_position<22):
	if (first_heavy_cys_position>16):
		temp_seq=VH_starting_AAs[:21-first_heavy_cys_position]
		temp_seq=temp_seq+heavy_seq
		heavy_seq=temp_seq
	else:
		sys.exit('Error: The length of the VH sequnece is too short. Please check input')
if (H3_tail_length<9):
	if (H3_tail_length>5):
		heavy_seq+=VH_ending_AAs[H3_tail_length-9:]
	else:
		sys.exit('Error: The length of the VH sequnece is too short. Please check input')
		
H3_len=len(heavy_seq)-len(heavy_seq_CC)-20



########################LIGHT CHAIN########################
check_variable_sequecne(light_seq,light_seq_CC, light_seq_CF, 22,6,VL_starting_AAs,VL_ending_AAs)
first_light_cys_position=(light_seq.find(light_seq_CC))
F_light_position=(light_seq.find(light_seq_CF)+len(light_seq_CF)-1)
L3_tail_length=len(light_seq)-F_light_position-1
if (first_light_cys_position<22):
	if (first_light_cys_position>16):
		temp_seq=VL_starting_AAs[:21-first_light_cys_position]
		temp_seq=temp_seq+light_seq
		light_seq=temp_seq
	else:
		sys.exit('Error: The length of the VL sequnece is too short. Please check input')
if (L3_tail_length<6):
	if (L3_tail_length>2):
		light_seq+=VL_ending_AAs[L3_tail_length-6:]
	else:
		sys.exit('Error: The length of the VL sequnece is too short. Please check input')
L3_len=len(light_seq)-len(light_seq_CC)-20

	
#pdb.set_trace()


#print (chains.input_seq)
print ("\n")
H1_H2_len =  (len(heavy_seq_CC)*4)-4 #The minus 4 and times 4 is to calculate the dihedral db length
L1_L2_len =  (len(light_seq_CC)*4)-4 #The minus 4 and times 4 is to calculate the dihedral db length
L3_len=(L3_len+1)*4
H3_len=(H3_len+1)*4


f = open('seq_flag','w')
f.write("-parser:script_vars sequence="+str(light_seq)+str(heavy_seq))
f.close()
#pdb.set_trace()
creat_run_cmd = './creat_run.sh ' +str(L1_L2_len)+" "+str(L3_len)+" "+str(H1_H2_len)+" "+str(H3_len)
print (creat_run_cmd)
hmmscan_proc = subprocess.run([creat_run_cmd], shell=True, check=True)
comm = 'sh command'
com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)
job_ids=re.findall('\d+',str(com_proc.stdout))
num_of_jobs = len(job_ids)
print (num_of_jobs,"jobs were submitted")
count=0
job_status = ["DONE", "EXIT"]
while count<num_of_jobs:
	for i in job_ids:	
		job_com = subprocess.run(["bjobs "+str(i)], shell=True, check=True,stdout=subprocess.PIPE)
		print (job_com.stdout)
		if any(x in str(job_com.stdout) for x in job_status):
			print ("yes")
			count = count+1
			job_ids.remove(i)
	print(count)

#assume all jobs finished correctly. In the future need to verify this

# get lowest energy scoring pdb
pdb_files = os.listdir(cwd+"/pdb/")
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
if not os.path.isdir("measure_rms"):
	os.makedirs("measure_rms")

for filename in glob.glob(os.path.join("/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/measure_rms_scripts", '*')):
	shutil.copy(filename, cwd+"/measure_rms/")
os.chdir("measure_rms")

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
	print(count)
import time
time.sleep(60) 
comm = 'sh post_run_command'
com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)

os.chdir("../")

if not os.path.isdir("top_models"):
	os.makedirs("top_models")

for filename in glob.glob(os.path.join("/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/top_models", '*')):
	shutil.copy(filename, cwd+"/top_models/")


models_sorted_path=cwd+"/measure_rms/cdr_co_rms"

os.chdir(cwd+"/top_models/")

comm = './creat_run.sh '+models_sorted_path
com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)

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

f= open("top_models", 'rb')
i=0
for line in f:
	i=i+1
	for file in glob.glob("../pdb/"+str(line,"utf-8").rstrip()+"*"):
		shutil.copy(file, ".")
		file=file.split("/")[-1]
		print (file)
		os.rename(file, str(i)+".pdb")
	

os.chdir(cwd)

if not os.path.isdir("plots"):
	os.makedirs("plots")

for filename in glob.glob(os.path.join("/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/plot", '*')):
	shutil.copy(filename, cwd+"/plots/")

os.chdir("plots/")

comm = '/apps/RH7U2/gnu/python/2.7.12/bin/python2 plot_rms_plots.py ../measure_rms/cdr_co_rms ../top_models/top_models'
com_proc = subprocess.run([comm], shell=True, check=True,stdout=subprocess.PIPE)













	#print(count)
		
