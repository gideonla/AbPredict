from ...pdb_.util import *
from ..rosetta_job import *
import os
from collections import OrderedDict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


def get_blade_limits_in_chimera(template: str, start_res: int,
                                end_res: int, blade: str) -> dict:
    """
    Assumes that template starts from 1 !!!
    This function is used for splice out only (assumes only one blade inserted)
    :param template: path to template pdb
    :param start_res: real index of the start of the blade
    :param end_res: real index of the end of the blade
    :param blade: path to blade pdb (find best match) or the sequence of the blade
    :return: dictionary of limits of before, after and within the blade
    """
    numres_template = get_length(template)
    len_template_after_blade = numres_template - end_res

    if os.path.isfile(blade):
        blade_len = get_length(blade) - 2 # -2 because the find best match blade is longer from the spliced blade in both ends
    else: # if it is the sequence
        blade_len = len(blade)

    limits = dict()
    limits['before'] = (1, start_res - 1)
    limits['blade'] = (start_res, start_res + blade_len - 1)
    limits['after'] = (start_res + blade_len,
                       start_res + blade_len + len_template_after_blade - 1)
    return limits


def splice_out_create_run(template: str, start_res: int, end_res: int,
                          blade: str, dir_name=''):
    """

    :param dir_name:
    :param template: path to template pdb
    :param start_res: index of the start of the blade (one into the blade)
    :param end_res: index of the end of the blade (one into the blade)
    :param blade: path to blade pdb
    :return:
    """
    if not dir_name:
        dir_name = os.path.basename(os.getcwd())

    name = os.path.basename(blade).replace('.pdb', '')
    limits = get_blade_limits_in_chimera(template, start_res, end_res, blade)

    args = ['@flags',
            '-parser:script_vars source=' + os.path.abspath(blade),
            'db_name={}_{}'.format(dir_name, name),
            'start_res=' + str(start_res),
            'end_res=' + str(end_res),
            'before=1-' + str(limits['before'][1]),
            'after={}-{}'.format(limits['after'][0], limits['after'][1]),
            '-out:prefix ' + name]

    create_job(name, args=args, dirs=['pdb', 'db'])


def get_chimera_limits(blades: list) -> dict:
    """
    Returns the start and end residue of each inserted blade in the context of the chimeric protein
    :param blades: list of tuples of (start_res, end_res, blade_length). Assumes
    each end res is bigger from its start_Res and smaller from the next start_res
    :return:
    """
    blades.sort(key=lambda x: x[0]) # sort by starting res
    limits = list()

    first_start = blades[0][0]
    first_end = blades[0][0] + blades[0][-1] - 1
    limits.append((first_start, first_end))
    prev_end = first_end
    for i in range(1, len(blades)):
        start, end, blade_len = blades[i]
        fixed_template = start - blades[i-1][1] - 1
        new_start = prev_end + 1 + fixed_template
        new_end = new_start + blade_len - 1
        limits.append((new_start, new_end))
        prev_end = new_end


    return limits


def save_designs_fasta(designs: list, path='splice_designs.fa') -> str:
    """
    Saves a fasta of splice designs
    :param designs: list of SpliceDesigns
    :return: path where fasta file was saved
    """
    records = [SeqRecord(Seq(des.sequence), id=des.path, description='') 
               for des in designs]
    SeqIO.write(records, path, 'fasta')
    return path
