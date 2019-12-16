import subprocess as sp
import os
from time import sleep
import logging
from flab.rosetta.RosettaPaths import RosettaPaths


logger = logging.getLogger(__name__)
CHAIN = 'X'


def create_sdf_from_ligand(pdb_path, ligand_name, sdf_output):
    """

    :param pdb_path:
    :param ligand_name:
    :param ligand_name::param sdf_output:
    :param ligand_name::return:
    :param ligand_name::return:"""
    import pymol                                    
    pymol.finish_launching(['pymol', '-cq'])
    pymol.cmd.delete('all')

    pymol.cmd.load(pdb_path)
    pymol.cmd.create('lig', selection='(resn {})'.format(ligand_name))
    pymol.cmd.save(sdf_output, selection='lig', format='sdf')
    pymol.cmd.delete('all')
    if not open(sdf_output).read(): 
        raise ValueError('No ligand named {} in file {}'.format(ligand_name, 
                                                                pdb_path)) 


def mol_to_param(sdf_path, ligand_name, path, r_paths=None):
    """

    :param sdf_path: path to sdf file containing the ligand
    :param ligand_name: wanted resname in the pdb created for the ligand
    :param path: directory where to place the output
    :return: a tuple of (path to ligand's pdb, path to params file)
    """
    r_paths = r_paths if r_paths else RosettaPaths() 
    path_start = os.getcwd()
    os.chdir(path)
    script = os.path.join(r_paths.source,
                          'scripts/python/public/molfile_to_params.py')
    logger.info('mol_to_param command: ' + ' '.join(['python2', script,
                                                     '-n', ligand_name, '-p', 
                                                     ligand_name, '--chain', 
                                                     CHAIN, sdf_path]))
    process = sp.Popen(['python2', script, '-n', ligand_name, '-p', ligand_name, '--chain', CHAIN, sdf_path],
                       stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = process.communicate()
    process.wait()
    lig_pdb = path + '/' + ligand_name + '_0001.pdb'
    params = path + '/' + ligand_name + '.params'
    os.chdir(path_start)
    return lig_pdb, params


def unite_pdbs(pdb_protein, pdb_ligand, ligand_name,
               combined_path='protein_n_ligand.pdb'):
    """
    Removes the original ligand lines from the pdb, then concatenates the new
    ligand lines where the original file had a TER line
    :param combined_path: path to new pdb file
    :param pdb_protein: path to pdb of the whole enzyme
    :param pdb_ligand:
    :param ligand_name: resname of the ligand
    :return:
    """
    protein = [line for line in open(pdb_protein, 'r').readlines()
               if ligand_name not in line]
    ligand = open(pdb_ligand, 'r').readlines()
    index = [i for i in range(len(protein)) if protein[i].startswith('TER')]
    index = -1 if 0 == len(index) else index[0]
    combined = protein[0:index] + ligand + protein[index:]
    open(combined_path, 'w').writelines(combined)


def create_params_file(pdb_path, ligand_name, output_path='', r_paths=None,
                       unite_pdb=None):
    """
    Creates params file for ligand in the pdb file and updates this ligand
    according to the new pdb generated from mol_to_params
    :param pdb_path:
    :param ligand_name: resname of ligand
    :param output_path: a directory to place the output (sdf, new pdb, params
                        etc.). Default: pwd
    :param unite_pdb: If specified, will add the correct ligand to this instead
    to the original pdb
    :return: path to new pdb file and param file
    """
    logger.debug('pdb path ' + pdb_path)
    logger.debug('ligand_name ' + ligand_name)
    if not os.path.isdir(output_path):
        output_path = os.getcwd()
    logger.debug('output_path ' + output_path)
    output_path = os.path.abspath(output_path)
    sdf_output = os.path.join(output_path, 'lig_tmp.sdf')
    create_sdf_from_ligand(pdb_path, ligand_name, sdf_output)
    lig_pdb, params = mol_to_param(sdf_output, ligand_name, output_path,
                                   r_paths)
    
    unite_pdb = unite_pdb if unite_pdb else pdb_path
    combined_path = '{}/{}_{}.pdb'.format(
        output_path,
        os.path.basename(unite_pdb).replace('.pdb', ''),
        ligand_name)
    unite_pdbs(unite_pdb, lig_pdb, ligand_name, combined_path)
    return combined_path, params


if __name__ == '__main__':
    create_sdf_from_ligand('4pud.pdb', 'XYP', 'lala.sdf')



