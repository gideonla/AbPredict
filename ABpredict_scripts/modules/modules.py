from ChainExtract import ChainExtract
import sys
import warnings
import Bio
import pdb
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO
    from Bio import SeqIO
import logging


VH_starting_AAs="VQL"
VH_ending_AAs="VTVS"
VL_starting_AAs="IEMT"
VL_ending_AAs="GTKL"

logger = logging.getLogger(__name__)

def is_fasta(filename):
	pdb.set_trace()
	with open(filename, "r") as handle:
		fasta = SeqIO.parse(handle, "fasta")
		pdb.set_trace()
		try:
			return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file
		except:
			logger.error("Input file is not a text file. Exiting")
			sys.exit("Input file is not a text file. Exiting")


def extract_seqeunces(chains:ChainExtract,hmm_file:str):
	chains.vsets_hmm = hmm_file

	dtbl = chains.find_domains(chains.full_seq, heavy=True, light=True)
	chains.parse_hmmscan_dtblout(dtbl, heavy=True, light=True)
	chains.classify_and_extract_chains(chains.full_seq, heavy=True, light=True)
	print('Heavy chain sequence: ')
	print (chains.heavy)
	print(chains.heavy_chain_seq)
	heavy_seq = chains.heavy_chain_seq
	print ('\n')
	print ('Light chain sequence: ')
	print (chains.light)
	print(chains.light_chain_seq)
	light_seq = chains.light_chain_seq
	return heavy_seq,light_seq


def check_variable_sequecne(variable_seq:Bio.Seq.Seq, variable_seq_CC:Bio.Seq.Seq, variable_seq_CX:Bio.Seq.Seq, cys_number,tail_length_const,starting_AAs,ending_AAs):
	first_cys_position=(variable_seq.find(variable_seq_CC))
	CDR3_heavy_position=(variable_seq.find(variable_seq_CX)+len(variable_seq_CX)-1)
	tail_length=len(variable_seq)-CDR3_heavy_position-1
	
	if (first_cys_position<cys_number):
		if (first_cys_position>(cys_number-6)): #If the first cys position is too short then padding is not enough and run fails
			temp_seq=starting_AAs[:(cys_number-1)-first_cys_position]
			temp_seq=temp_seq+variable_seq
			variable_seq=temp_seq
		else:
			logger.error("Variable domain sequence is too short in the N-terminl domain. Exiting")
			sys.exit("Variable domain sequence is too short in the N-terminl domain. Exiting")

	if (tail_length<tail_length_const):
		if (tail_length>(tail_length_const-4)):
			variable_seq+=ending_AAs[tail_length-tail_length_const:]
		else:
			logger.error("Variable domain sequence is too short in the N-terminl domain. Exiting")
			sys.exit("Variable domain sequence is too short in the N-terminl domain. Exiting")
	return variable_seq
		


