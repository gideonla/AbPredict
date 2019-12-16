#!/usr/bin/env python3
import logging
from collections import OrderedDict 
from Bio.PDB.PDBParser import PDBParser
from flab.pdb_.util import get_sequence
from Bio.SeqUtils import seq1, seq3
from flab.pssm.pssm import PSSM 
# from pross2.pssm.PSSM import PSSM
# from flab.pdb_.MyPDB_funcs import parse_PDB
import re
import numpy as np
from copy import deepcopy
import os


logger = logging.getLogger(__name__)


class ResFile:
    """Currently works only on resfiles from type PIKAA"""
    def __init__(self, path: str='', pdb='', pssm=''):
        """
        :param pdb: path to pdb the resfile belonges to (position nubering is
        by it)
        """
        self.path = path
        self.resfile = (self.parse_res_file(self.path) 
                        if self.path else OrderedDict())
        self.rearranged = False
        self.set_pdb(pdb) 
        self.set_pssm(pssm)  

    @staticmethod
    def parse_res_file(resfile_path):                                                
        """                                       
        :param resfile_path:
        :return: list of tuples (residueCHAIN, the allowed residues as a string)
        """
        resfile = open(resfile_path, 'r').readlines()
        resfile = OrderedDict((''.join(line.split()[:2]),
                               line.split()[-1].upper())
                              for line in resfile if 'PIKAA' in line)
        # sorts by residue number to mutate
        resfile = OrderedDict(sorted(resfile.items(), 
                                     key=lambda x: (x[0][-1], int(x[0][:-1]))))
        logger.debug(resfile)
        return resfile
    
    def init_from_dict(self, res_dict):
        """
        :param res_dict: [position] = string of AAs
        """
        for pos, allowed in res_dict.items():
            if pos in self.resfile.keys():
                logger.warning('replacing position {} from {} to {}'.format(
                    pos, self.resfile[pos], allowed))
            self.resfile[pos] = allowed
        resfile = OrderedDict(sorted(self.resfile.items(),                             
                                   key=lambda x: (x[0][-1], int(x[0][:-1]))))
        self.resfile = resfile

    def set_pdb(self, pdb_path):
        self.pdb = PDBParser().get_structure('', pdb_path) if pdb_path else ''
        # self.pdb = parse_PDB(pdb_path) if pdb_path else None 
    
    def set_pssm(self, pssm_path):
        # self.pssm = PSSM(name='', pssm_file=pssm_path) if pssm_path else None
        self.pssm = PSSM(pssm_path) if pssm_path else None
    
    def __len__(self):
        return len(self.resfile) 

    def __str__(self):
        return '\n'.join(['{:<6} {}'.format(k, v) for k, v in
                          self.resfile.items()])
    
    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        new_rf = deepcopy(self) 
        if self.pdb and other.pdb:
            if get_sequence(self.pdb) != get_sequence(other.pdb): 
                raise ValueError(('concatenated resfiles should come from the '
                                  'same protein:\n {}\n{}').format(self.pdb,
                                                                   other.pdb))
        for residue, allowed in other.resfile.items():
            if residue in new_rf.resfile.keys():
                raise ValueError('Residue %s appears in both resfiles' % residue)
            new_rf.resfile[residue] = allowed
        return new_rf
    
    def __iter__(self):
        for residue, allowed in self.resfile.items():
            yield residue, allowed

    def write_resfile(self, path):
        """Writes the resfile to a new file located in path"""
        sorted_keys = sorted(self.resfile.keys(), 
                             key=lambda x: int(x[:-1]))
        lines = ['nataa\nstart\n']
        # lines += ['{:<8}{}       PIKAA   {}\n'.format(k[:-1], k[-1],
        #                                               self.resfile[k])
        #           for k in sorted_keys]
        lines += ['{}\t{}\tPIKAA\t{}\n'.format(k[:-1], k[-1], self.resfile[k])
                  for k in sorted_keys]
        open(path, 'w').writelines(lines)

    def remove_min_pssm_probability(self, min_aa_probablity=-2):            
        """
        Assumes the pssm has only one chain !!!!
        :param min_aa_probablity: equal or bigger stay, smaller probabilities
        are removed
        """
        assert self.pssm, 'PSSM was not set to resfile'
        assert self.pdb, "pdb file was not set"
        result = OrderedDict()
        for residue, allowed in self.resfile.items():
            res = int(residue[:-1])
            chain = residue[-1]
            
            org_res = self.pdb[0][chain][res]
            line = list(self.pdb[0][chain]).index(org_res)
            if self.rearranged:
                # varifies that it found the right line
                assert self.pssm.get_wt(line) == allowed[0]
            
            new_allowed = [aa for aa in allowed
                           if self.pssm.get_score(line=line,
                                             aa=aa) >= min_aa_probablity]
            result[residue] = ''.join(new_allowed)
        self.resfile = result
    
    def rearrange_resfile(self, remove_only_wt=True):
        """
        Rearranges the allowed residues to place thw WT residues first. 
        :param resfile: a resfile parsed to tuples (from parse_res_file).
        :param pdb_path: path to a PDB file corresponding to the resfile
        :param remove_only_wt: if True, removes records with only the WT
        allowed 
        """
        assert self.pdb, "pdb file was not set"
        resfile_new = OrderedDict()
        for residue, allowed_res in self.resfile.items():
            res = int(residue[:-1])
            chain = residue[-1]
            wt = seq1(self.pdb[0][chain][res].resname)
            if 1 == len(allowed_res) and allowed_res == wt and remove_only_wt:
                continue
            try:
                wt_i = allowed_res.index(wt)
                allowed_res = wt + allowed_res[:wt_i] + allowed_res[wt_i + 1:]
            except ValueError: # if wt is not in the residues
                logger.warning(('WT residue {} was not in position {} in the res'
                             ' file').format(wt, residue)) 
                allowed_res = wt + allowed_res
            resfile_new[residue] = allowed_res 
        self.resfile  = resfile_new           
        self.rearranged = True

    def get_allowed(self, position):
        return self.resfile[position] if position in self.resfile else ''

    def get_positions(self):
        return list(self.resfile.keys())
    
    def serial2seq(self, serialnum):
        """given a serial number in the length of 2*(len of res file) returns
        the sequnce it represents according to the order of the resfile
        """
        return OrderedDict([(k, v[int(n) - 1])
                    for (k, v), n in zip(self.resfile.items(), 
                                         re.findall('.{2}', serialnum))])

    def possible_permutations(self):
        return np.prod([len(v) for k, v in self.resfile.items()])


def get_allowed_residues(resfile_path,                              
                         pdb_path,
                         pssm_path='',
                         min_aa_probablity=-2,
                         path='',
                         print_ordered=True):
    """
    Parses the resfile, removes non allowed res by PSSM score and rearrages the
    resfile to have the WT residue to be the first residue
    """
    resfile = ResFile(path=resfile_path, pdb=pdb_path, pssm=pssm_path)
    if pssm_path:
        resfile.remove_min_pssm_probability(min_aa_probablity=min_aa_probablity) 
    resfile.rearrange_resfile()
    if print_ordered:
        path = path if path else resfile_path + '.ordered'
        resfile.write_resfile(path)
    return resfile


def unify_res_file(paths, pdb=None, out_path=os.getcwd()):
    """Unifies resfiles, each threshold into a different file.
    The files shoukd obviously come from the same protein
    :param paths: paths to resfiles. should not have points in the file name
    except from the ones indicating the threshold
    :param pdb: path to a pdb file. will change the first aa in each row to the
    wt
    :param out_path: path to directory where to print the unified resfiles
    :return: a dictionary [threshold] = unified resfile
            and the paths to all the new unified resfilles
    """
    d = dict()
    for f in paths:
        threshold = float(f[f.find('.') + 1:])
        rf = ResFile(f)
        if pdb:
            rf.set_pdb(pdb)
            rf.rearrange_resfile()
        if threshold in d.keys():
            d[threshold] += rf 
        else:
            d[threshold] = rf
    out_paths = list()
    for threshold, resfile in d.items():
        path = os.path.join(out_path, 'designable_aa_resfile.' + str(threshold))
        out_paths.append(path)
        resfile.write_resfile(path)
    return d, out_paths
