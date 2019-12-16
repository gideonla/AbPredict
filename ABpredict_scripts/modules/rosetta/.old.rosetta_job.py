#!/bin/env python3


import os
from time import sleep
from functools import reduce


ROSETTA = '/home/labs/fleishman/rosaliel/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease'

def create_dir(path):
    """Creates the directory in path
    :param path: path to the new folder to create
    """
    if not os.path.exists(path):
        os.mkdir(path)


def create_file(path):
    if not os.path.isfile(path):
        open(path, 'a').close()


def create_dirs(path, dir_names):
    """Creates directories in path
    :param path: path to directory where to locate all newly created directories
    :param dir_names: directories names to create
    """
    for dir in dir_names:
        create_dir(os.path.join(path, dir))



def create_job(name, path='', queue='fleishman', rusage='2048', args=list(),
               dirs=list()):
    """
    Creates a job file and add command to the command file
    :param name: the name of the job
    :param path:
    :param queue:
    :param rusage:
    :param args:
    :return:
    """
    if path == '':
        path = os.getcwd()
    job_path = '{}/jobs/job.{}'.format(path, name)
    command_path = path + '/command'
    out_path = '{}/out/out.{}'.format(path, name)
    err_path = '{}/err/err.{}'.format(path, name)
    create_dirs(path, ['err', 'out', 'jobs'] + dirs)
    create_file(command_path)
    os.chmod(command_path, 0o774)

    # Creates the job
    job = list()
    job.append('#!/bin/bash')
    job.append('cd ' + path)
    last_line = ROSETTA
    if args: 
        vars = ''    
        for var in args:
            vars += str(var) + ' '
        job.append('{} {}'.format(last_line, vars))
    else:
        job.append(last_line)    
    open(job_path, 'w').writelines([line + '\n' for line in job])
    os.chmod(job_path, 0o774)
    
    # adds the line to the command file
    command = 'bsub  -C 1024 -u /dev/null '
    command += '-R rusage[mem={}] '.format(rusage)
    command += '-L /bin/bash -G fleishman-wx-grp-lsf '
    command += '-q {} '.format(queue) 
    command += '-o {} '.format(out_path)
    command += '-e {} '.format(err_path)
    command += job_path
    open(command_path, 'a').write(command + '\n')


def prepare_rosetta_flag(path, flags):
    """Creates a flag file from a flags dictionary"""
    flags_file = open(path, 'w')
    for k in list(flags.keys()):
        if k.endswith('='):
            flags_file.write(k + flags[k] + '\n')
        else:
            flags_file.write(k + ' ' + flags[k] + '\n')


def parse_rosetta_flags(path):
    """Parses a flags file to a dictionary. Ignores lines starting with #."""
    flags = dict()
    for line in open(path):
        if line.startswith('#') or line.strip() == '': continue
        if line.startswith('-parser:script_vars'):
            var = line.split()[1].split('=')[0]
            flags['-parser:script_vars ' + var + '='] = line.split('=')[1]
        else:
            k = line.split()[0]
            value = line.split()[1:]
            if 0 == len(value):
                value = ''
            elif 1 == len(value):
                value = value[0]
            else:
                value = reduce(lambda x, y: '{} {}'.format(x, y), value)
            flags[k] = value

    return flags
    
    
