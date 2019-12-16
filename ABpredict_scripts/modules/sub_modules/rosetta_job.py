#!/bin/env python3

import os
from functools import reduce
from collections import OrderedDict
from utils.file_system import create_dirs, create_file, copyanything
import logging
from rosetta.RosettaPaths import RosettaPaths

logger = logging.getLogger(__name__)


class RosettaJob:
    """Handeling a single Rosetta job"""

    def __init__(self, name: str, path: str, flags, r_paths: RosettaPaths=None):
        """
        :param path: path to directory containing all files to run the job
        :param flags: an instance of RosettaFlags:
            """
        self.r_paths = r_paths if r_paths else RosettaPaths('~/Rosetta')
        self.name = str(name)
        self.path = path
        self.flags = flags

    def create_job(self,
                   command_path='',
                   queue='fleishman',
                   rusage='2048',
                   dirs=[],
                   dir_count=''):
        """
        Creates a job file and add command to the command file
        :param path:
        :param queue:
        :param rusage:
        :param args:
        :param dir_count: the number of err/ out/ etc. dir to send th output to      
        :return:
        """
        if not os.path.isdir("jobs"):
            os.makedirs("jobs")
        if not os.path.isdir("out"):
            os.makedirs("out")
        if not os.path.isdir("err"):
            os.makedirs("err")
        if not os.path.isdir("pdb"):
            os.makedirs("pdb")

        dir_count = '_' + str(dir_count) if dir_count else ''
        self.job_path = '{}/jobs{}/job.{}'.format(self.path, dir_count, 
                                                  self.name)
        if not command_path:
            command_path = os.path.join(self.path, 'command')
        self.command_path = command_path
        out_path = os.path.join(self.path, 'out' + dir_count, 
                                'out.' + self.name)
        err_path = os.path.join(self.path, 'err' + dir_count,
                                'err.' + self.name)
        # dirs = list(set(['err', 'out', 'jobs', 'pdbs', 'bin', 'data'] + dirs))
        # create_dirs(self.path, dirs)
        # create_file(command_path)
        # os.chmod(command_path, 0o774)

        # Creates the job
        job = list()
        job.append('#!/bin/bash')
        job.append('cd ' + self.path)
        job.append(self.r_paths.r_scripts() + ' ' + str(self.flags))
        print (self.r_paths.r_scripts() + ' ' + str(self.flags))
        open(self.job_path, 'w').writelines([line + '\n' for line in job])
        os.chmod(self.job_path, 0o774)

        # adds the line to the command file
        command = 'bsub '
        command += '-R rusage[mem={}] '.format(rusage)
        command += '-G fleishman-wx-grp-lsf '
        command += '-q {} '.format(queue)
        command += '-o {} '.format(out_path)
        command += '-e {} '.format(err_path)
        command += self.job_path
        open(command_path, 'a').write(command + '\n')


class RosettaFlags:
    # TODO  parse flags file
    """Storing all types of flags and command line options for a rosetta job"""

    def __init__(self,
                 flags_files_path: list=None,
                 options: OrderedDict=None,
                 scripts_vars: OrderedDict=None):
        """
        """
        if not flags_files_path:
            self.flags_files = list()
        else:
            self.flags_files = (flags_files_path if list ==
                                type(flags_files_path) else [flags_files_path])
        self.options = options if options else OrderedDict()
        self.scripts_vars = scripts_vars if scripts_vars else OrderedDict()

    def add_option(self, option, value=''):
        """Replaces option's value if already exists"""
        if option in self.options.keys():
            logger.debug('changes {}: {} -> {}'.format(option, self.options[
                option], value))
        self.options[option] = value

    def add_script_var(self, var, value):
        if var in self.scripts_vars.keys():
            logger.debug('changes {}: {} -> {}'.format(var, self.scripts_vars[
                var], value))
        self.scripts_vars[var] = value

    def add_flags_file(self, path):
        self.flags_files.append(path)

    def __str__(self):
        result = ''
        if self.flags_files:
            result += self.__flags_str()
        result += ' ' + self.__options_str() + ' '
        if self.scripts_vars:
            result += self.__scripts_vars_str()
        return result

    def __add__(self, other):
        """
        Unites two flags. it they have common options/script vars it will take
        the value of the first one.
        """
        flags_files = self.flags_files + other.flags_files
        options = other.options
        msg = '{} is in both flags. Was:{} Set to:{}'
        for k, v in self.options.items():
            if k in options.keys():
                logger.warning(msg.format(k, options[k], v))
            options[k] = v
        scripts_vars = other.scripts_vars
        for k, v in self.scripts_vars.items():
            if k in scripts_vars.keys():
                logger.warning(msg.format(k, scripts_vars[k], v))
            scripts_vars[k] = v
        return RosettaFlags(flags_files, options, scripts_vars)

    def __flags_str(self, oneline=True):
        concat = ' ' if oneline else '\n'
        return concat.join(['@' + f for f in self.flags_files])

    def __options_str(self, oneline=True):
        """Returns a string representing all options.
        :param oneline: if false, each option is in a separate line
        """
        concat = ' ' if oneline else '\n'
        return concat.join(
            ['-{} {}'.format(k, v) for k, v in self.options.items()])

    def __scripts_vars_str(self, oneline=True):
        if oneline:
            return ('-parser:script_vars ' + ' '.join(
                ['{}={}'.format(k, v) for k, v in self.scripts_vars.items()]))
        return '\n'.join([
            '-parser:script_vars {}={}'.format(k, v)
            for k, v in self.scripts_vars.items()
        ])

    def write_flags_file(self, path='flags'):
        """Writes all flags and options in a flags file.
        :param path: a path to a flags file to write in to
        :return: path of written flags file
        """
        create_file(path)
        lines = ''
        if self.flags_files:
            lines += self.__flags_str(oneline=False) + '\n'
        lines += self.__options_str(oneline=False) + '\n'
        if self.scripts_vars:
            lines += self.__scripts_vars_str(oneline=False) + '\n'
        open(path, 'a').write(lines)
        return path


def parse_flags(path):
    flags = []
    options = OrderedDict()
    scripts_vars = OrderedDict()

    for line in open(path):
        line = line.strip()
        line = line.split('#')[0]
        # flags
        if line.startswith('@'):
            flags.append(line[1:])
        # empty or comments line
        elif not line.startswith('-'):
            continue
            # script vars
        elif line.startswith('-parser:script_vars'):
            for i in line.split()[1:]:
                k, v = i.split('=')
                scripts_vars[k] = v
        # options
        else:
            sp = line.split()
            options[sp[0][1:]] = ' '.join(sp[1:])
    return RosettaFlags(
        flags_files_path=flags, options=options, scripts_vars=scripts_vars)
