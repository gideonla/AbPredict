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
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
VH_starting_AAs="VQLE"
VL_starting_AAs="VLTQ"


cwd = os.getcwd()

for filename in glob.glob(os.path.join("/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/required_eden_Rosetta_files/", '*')):
	shutil.copy(filename, os.getcwd())
#	print (filename)

# setting up AbChains object
chains = ChainExtract()
num = chains.parse_fasta(argv[1])

# defining locations
chains.vsets_hmm = "/home/labs/fleishman/gideonla/ABpredict/test_cut_length/ChainExtract/hmms/all_cys.hmm"
#chains.vsets_hmm ="/home/labs/fleishman/gideonla/ABpredict/test_cut_length/ChainExtract/database/igg-chains_cdhit.hmm"

# doing the actual work
# First find sequence between cysteines 
dtbl = chains.find_domains(chains.full_seq, heavy=True, light=True)
chains.parse_hmmscan_dtblout(dtbl, heavy=True, light=True)
chains.classify_and_extract_chains(chains.full_seq, heavy=True, light=True)
print('Heavy chain sequence: ')
print (chains.heavy)
print(chains.heavy_chain_seq)
cys_heavy_seq = chains.heavy_chain_seq
print ('\n')
print ('Light chain sequence: ')
print (chains.light)
print(chains.light_chain_seq)
cys_light_seq = chains.light_chain_seq


# Getting sequencees using Jake's hmms
#chains.vsets_hmm = "/home/labs/fleishman/gideonla/ABpredict/test_cut_length/ChainExtract/hmms/all_cys.hmm"
#chains.vsets_hmm ="/home/labs/fleishman/gideonla/ABpredict/test_cut_length/ChainExtract/database/igg-chains_cdhit.hmm"
chains.vsets_hmm ="/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/hmmscan/kappa_C_F_heavy_smart_domains_aln_W.hmm"

dtbl = chains.find_domains(chains.full_seq, heavy=True, light=True)
chains.parse_hmmscan_dtblout(dtbl, heavy=True, light=True)
chains.classify_and_extract_chains(chains.full_seq, heavy=True, light=True)
print('**********Heavy chain sequence: ')
print (chains.heavy)
print(chains.heavy_chain_seq)
heavy_seq = chains.heavy_chain_seq
print ('\n')
print ('Light chain sequence: ')
print (chains.light)
print(chains.light_chain_seq)
light_seq = chains.light_chain_seq

#Get HMMs of VH from first cys to W at end of H3
fh = open("./sequences/query_heavy.fa","w")
fh.write(">heavy\n")
fh.write(str(heavy_seq))
fh.close()
vh_c_W = ChainExtract()
vh_c_W.parse_fasta("./sequences/query_heavy.fa")
vh_c_W.vsets_hmm ="/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/hmmscan/heavy_C_W_hmm/heavy_smart_domains_aln_W.hmm"
record = SeqRecord(heavy_seq,id="heavy", name="HokC")
dtbl = vh_c_W.find_domains(heavy_seq, heavy=True, light=0)
vh_c_W.parse_hmmscan_dtblout(dtbl, heavy=True, light=0)
#record = SeqRecord(heavy_seq,IUPAC.protein, id="heavy")
record = SeqRecord(heavy_seq,id="heavy", name="HokC",description="toxic membrane protein")
#pdb.set_trace()
vh_c_W.classify_and_extract_chains(record, heavy=True, light=0)
vh_c_W_seq=str(vh_c_W.heavy_chain_seq)


#Get HMMs of VL from first cys to F at end of L3
fl = open("./sequences/query_light.fa","w")
fl.write(">light\n")
fl.write(str(light_seq))
fl.close()
vl_c_W = ChainExtract()
vl_c_W.parse_fasta("./sequences/query_light.fa")
vl_c_W.vsets_hmm ="/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/hmmscan/light_C_W_hmm/kappa_C_F.hmm"
record = SeqRecord(light_seq,id="light", name="HokC")
dtbl = vl_c_W.find_domains(light_seq, heavy=0, light=1)
vl_c_W.parse_hmmscan_dtblout(dtbl, heavy=0, light=1)
#record = SeqRecord(heavy_seq,IUPAC.protein, id="heavy")
#pdb.set_trace()
vl_c_W.classify_and_extract_chains(record, heavy=0, light=1)
vl_c_W_seq=str(vl_c_W.light_chain_seq)


pdb.set_trace()



#print (chains.input_seq)
print ("\n")
H1_H2_len =  (len(cys_heavy_seq)*4)-4 #The minus 4 and times 4 is to calculate the dihedral db length
L1_L2_len =  (len(cys_light_seq)*4)-4 #The minus 4 and times 4 is to calculate the dihedral db length

#If the input sequence is too short in the tail regions then we need to pad it 
first_heavy_cys_position=(heavy_seq.find(cys_heavy_seq))
if (first_heavy_cys_position<21):
	temp_seq=VH_starting_AAs[:20-first_heavy_cys_position]
	temp_seq=temp_seq+heavy_seq
	heavy_seq=temp_seq
	
first_light_cys_position=(light_seq.find(cys_light_seq))
if (first_light_cys_position<22):
	temp_seq=VL_starting_AAs[:21-first_light_cys_position]
	temp_seq=temp_seq+light_seq
	light_seq=temp_seq


H1_H2_seq=heavy_seq[0:20+len(cys_heavy_seq)]
L1_L2_seq=light_seq[0:21+len(cys_light_seq)]

print(H1_H2_seq)
print("*****")
print(L1_L2_seq)

#print (H1_H2_seq.tostring().replace(H1_H2_seq.tostring(),''))
second_heavy_cys=str(heavy_seq).index(str(cys_heavy_seq))+len(cys_heavy_seq)
second_light_cys=str(light_seq).index(str(cys_light_seq))+len(cys_light_seq)
H3_seq=(str(heavy_seq)[second_heavy_cys:])
L3_seq=(str(light_seq)[second_light_cys:])
print (H3_seq)
H3_len = (len(H3_seq)*4)+4
L3_len = (len(L3_seq)*4)+4 
print (H3_len)
print (L3_len)

f = open('seq_flag','w')
f.write("-parser:script_vars sequence="+str(L1_L2_seq)+str(L3_seq)+str(H1_H2_seq)+str(H3_seq))
f.close()
pdb.set_trace()
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
		
