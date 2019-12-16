#!/bin/env python3


from pross2.utils import file_system
from .rosetta_job import RosettaFlags, RosettaJob, parse_flags
from .RosettaPaths import RosettaPaths
import os
import copy
from flab.lsf.job_executer import JobsExecuter
import datetime as dt
import logging


logger = logging.getLogger(__name__)


class Rosetter:
    """Creates all directories and files necessary to run a rosetta
    application
    """
    def __init__(self, name: str,
                 r_paths: RosettaPaths=None,
                 flags: RosettaFlags=None,
                 data=[],
                 dirs=[],
                 path=os.getcwd(),
                 max_dir=5000,
                 dir_count=False) -> None:
        """
        :param data: list of path of data files to copy to data directory
        :param dir_count: if true, will separate the output of each max_dir
        jobs
        """
        self.name = name
        self.r_paths = r_paths if r_paths else RosettaPaths()
        # directories
        self.dir_count = 1 if dir_count else False
        self.updated = True
        self.max_dir = max_dir
        self.path = os.path.join(path, name)
        file_system.create_dir(self.path)
        self.__create_dirs(dirs)

        # flags
        flags = flags if flags else RosettaFlags()
        flags.add_option('database', self.r_paths.db())
        # flags.add_option('out:path:pdb', self.path + '/pdbs/')
        self.flags = RosettaFlags()
        self.flags_path = os.path.join(self.path, 'bin', 'flags')
        self.set_flags(flags)
        # command file
        self.command = os.path.join(self.path, 'bin', 'command')
        file_system.create_file(self.command)
        os.chmod(self.command, 0o774)

        for d in data:
            file_system.copyanything(d, os.path.join(self.path, 'data'))

        self.jobs = list()

    def set_flags(self, flags_: RosettaFlags, overwrite=False):
        """set_flags
        Unites the new flags with the current flags of the rosetter.
        :param flags_:
        :type flags_: RosettaFlags
        :param overwrite: if true, in a conflict will keep the new flags values
        """
        self.flags = (self.flags + flags_
                      if not overwrite else flags_ + self.flags)
        self.write_flags()

    def write_flags(self):
        self.flags.write_flags_file(self.flags_path)

    def __create_dirs(self, dirs): 
        """
        This function is used ONLY IN THE INIT METHOD
        :param dirs: additional dirs to create
        """
        self.dirs = list(set(dirs + ['err', 'out', 'jobs']))#, 'pdbs']))
        if self.dir_count:
            self.dirs = [d + '_' + str(self.dir_count) for d in self.dirs]
        
        dirs = list(set(self.dirs + ['bin', 'data']))
        file_system.create_dirs(self.path, dirs)

    def __update_dir_count(self):
        if not self.dir_count:
            return
        self.dir_count += 1
        self.dirs = [d[:d.rfind('_')] + '_' + str(self.dir_count)
                     for d in self.dirs]
        file_system.create_dirs(self.path, self.dirs)
        self.updated = True
        logger.warning('updating the dir count to ' + str(self.dir_count))
        return self.dir_count

    def create_job(self, name, flags: RosettaFlags=None,
                   queue='fleishman', rusage='2048'):
        if not flags:
            flags = RosettaFlags()
        flags.add_flags_file(self.flags_path)
        job = RosettaJob(name, path=self.path,
                         r_paths=self.r_paths,
                         flags=flags)
        job.create_job(self.command, queue=queue, rusage=rusage,
                       dir_count=self.dir_count)
        self.jobs.append(job)
        if len(self.jobs) % self.max_dir == 0:
            self.__update_dir_count()

    def run(self, N_in_parallel=3000, wait_interval=2):
        # TODO ROSALIE - add something that checks if there are any jobs that
        # won't finish for soemreason, kill and resubmit them
        start = dt.datetime.now()
        logger.info('starts running %s', self.name)
        je = JobsExecuter(command_path=self.command,
                          N_in_parallel=N_in_parallel,
                          wait_interval=wait_interval)
        je.run_jobs()
        time = (dt.datetime.now() - start).total_seconds() / 60
        logger.info('%s finished running. took %i', self.name, time)
        return time
