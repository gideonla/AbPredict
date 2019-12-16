from __future__ import print_function
import sys
import os
import shutil
import subprocess
from Bio import SeqIO, SearchIO
import logging
logger = logging.getLogger(__name__)

class ChainExtract:
    """
    A class for the processing of an input sequence for the antibody
    re-engineering pipeline.
    Dependencies: hmmscan3 (http://hmmer.org/)
    """

    # paths to chain HMMs
    vsets_hmm = ""
    heavy_hmm = ""
    kappa_hmm = ""
    lambda_hmm = ""
    input_seq = ""

    def __init__(self):
        ### TODO ###
        # set up dir structure for input/output files
        # logging

        # working directory
        self.home = os.getcwd()
        self.name = "query" # input name / job name
        self.full_seq = "" # biopython Seq object of input seq
        self.heavy_full_seq = ""
        self.light_full_seq = ""
        self.heavy = {'score': 0,'coords':(0,0)}
        self.light = {'score': 0,'coords':(0,0),
                      'chain_type': 'light'}
        self.heavy_chain_seq = ""
        self.light_chain_seq = ""


    def parse_fasta(self, fasta_path):
        """
        Takes the path to a fasta file and parses the contents.
        File can contain a single sequence with both heavy and light chains
        or two sequences containing only a heavy or light chain and labeled
        with 'heavy' or 'light' respectively.
        """
        if not os.path.isdir('./sequences'):
            os.mkdir('./sequences')
        shutil.copyfile(fasta_path, './sequences/input.fa')
        records = [rec for rec in SeqIO.parse(fasta_path, 'fasta')]
        self.input_seq=(records[0].seq)
        #print(records[0].seq)
        len_warning = "fasta should contain either one sequence with both "
        "chains or two sequences, each with one chain, labeled with 'light' "
        "and 'heavy' for each chain type."
        assert len(records) == 1 or len(records) == 2, len_warning

        if len(records) == 1:
            records[0].id = 'query_heavy_light'
            self.full_seq = records[0]
            SeqIO.write(records[0], './sequences/query.fa', 'fasta')

        elif len(records) == 2:
            for rec in records:
                if "heavy" in rec.description:
                    self.heavy_full_seq = rec.seq
                    SeqIO.write(rec, './sequences/query_heavy.fa', 'fasta')
                elif "light" in rec.description:
                    self.light_full_seq = rec.seq
                    SeqIO.write(rec, './sequences/query_light.fa', 'fasta')
                else:
                    exception_message = "Neither 'heavy' nor 'light' can be "
                    "found in the sequence description of {}.".format(rec)
                    raise Exception(exception_message)

        return len(records)


    def find_domains(self, seq, heavy=True, light=True):
        """
        Takes input sequence and classifies the chain as either heavy/light
        and light- kappa/lambda. Default tries to classify both chain types,
        for a heavy containing seq, pass light=Flase and visa versa.
        Returns a sequence of the correct length for AbPredict.
        """
        if not os.path.isdir('./hmmscan'):
            os.mkdir('./hmmscan')
        if heavy and light:
            infile = 'query.fa'
            outname = self.name
        elif heavy:
            infile = 'query_heavy.fa'
            outname = self.name + '_heavy'
        elif light:
            infile = 'query_light.fa'
            outname = self.name + '_light'

        hmmscan_cmd = 'hmmscan --domtblout ./hmmscan/{}.dout '.format(outname)
        hmmscan_cmd += '-o ./hmmscan/{}.out '.format(outname)
        hmmscan_cmd += '{} ./sequences/{}'.format(self.vsets_hmm, infile)
        hmmscan_proc = subprocess.run([hmmscan_cmd], shell=True, check=True)
        try:
            dtbl = SearchIO.read('./hmmscan/{}.dout'.format(outname),'hmmscan3-domtab')
        except:
            logger.error("Input file is not an antibody sequence. Exiting")
            sys.exit("Input file is not an antibody sequence. Exiting")
        return dtbl


    def parse_hmmscan_dtblout(self, dtbl, heavy, light):
        """
        Takes the dtbl output from an hmmscan run and returns the highest
        scoring domain for each HMM in the database.
        """

        for hmm in dtbl:
            dom = hmm.id.split("_")[0]
            for hit in hmm:
                score = hit.bitscore
                h = heavy and (dom == 'heavy') and (score > 80)
                l = (light and dom in ['kappa' , 'lambda']) and (score > 80)
                if h and score > self.heavy['score']:
                    self.heavy['score'] = score
                    self.heavy['coords'] = (hit.query_start,
                                            hit.query_end)
                if l and score > self.light['score']:
                    self.light['chain_type'] = dom
                    self.light['score'] = score
                    self.light['coords'] = (hit.query_start,
                                            hit.query_end)


    def classify_and_extract_chains(self, seq, heavy=True, light=True):
        """
        Chooses which light chain to use and extracts the sequences based
        on the coordinates from the hmmscan results.
        """

        if heavy:
            h_start = self.heavy['coords'][0]
            h_end = self.heavy['coords'][1]
            self.heavy_chain_seq = seq.seq[h_start:h_end]
        if light:
            l_start = self.light['coords'][0]
            l_end = self.light['coords'][1]
            self.light_chain_seq = seq.seq[l_start:l_end]
        if heavy and light:
            # do the domains overlap?
            overlap = ((h_start < l_start) and (h_end > l_start)) or \
                      ((h_start > l_start) and (h_end < l_end))
            if overlap:
                print("Warning: Highest scoring domains overlap.")
                print ("Heavy chain:\n", self.heavy)
                print ("Light chain:\n", self.light)
                exit(1)
