#!/bin/env python3
import sys
sys.path.append('/home/labs/fleishman/WebServer/ABPREDICT/ABpredict_scripts/modules/sub_modules/lsfe')
from . import lsf_util
from time import sleep
from getpass import getuser
import subprocess as sp
import logging
from datetime import datetime

logger = logging.getLogger(__name__)
USER = getuser()


class Job(object):
    def __init__(self, command=''):
        """

        :param command: a string of a bsubs
        """
        self.command = command.strip()
        self.jobID = -1
        self.queue = ''
        self.resubmitID = -1
        self.rqueue = ''

    def submit(self, resubmit=False):
        """

        :return: job ID if the job was executed successfully, false otherwise
        """
        if not self.command:
            logger.warning('no command assigned to job')
            return False

        process = sp.Popen(
            self.command.split(),
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True)
        output, error = process.communicate()
        if not output:  # job failed
            logger.warning('in submit: no output for {}'.format(self.command))
            return False
        id_ = output[output.index('<') + 1:output.index('>')]
        q = output[output.rfind('<') + 1:output.rfind('>')]
        if resubmit:
            self.resubmitID = id_
            self.rqueue = q
            return self.resubmitID
        else:
            self.jobID = id_
            self.queue = q
            return self.jobID

    def get_status(self):
        if -1 == self.jobID:
            return -1
        jobs = lsf_util.bjobs(jobID=self.jobID)
        if 0 == len(jobs):
            return 'DONE'
        else:
            self.status = jobs[jobs.JOBID == self.jobID].iloc[0]['STAT']
        return self.status

    def switch_queue(self, queue):
        switched = lsf_util.bswitch(queue=queue, jobs=[self.jobID])
        return self.jobID in switched


class JobsExecuter(object):
    def __init__(self,
                 command_path='',
                 N_in_parallel=3000,
                 wait_interval=2,
                 resubmit_min=60 * 24 * 7,
                 kill_min=60 * 24 * 7):
        """

        :param command_path: path to command file
        :param N_in_parallel: How many jobs to submit in parallel
        :param wait_interval: time (in minutes) to wait between peeking the jobs
        :param resubmit_min: how long before killing a job. default: one week
        :param resubmit_times: times to resubmit job
        """
        self.user = getuser()
        self.jobs = list()
        if command_path:
            self.command_path = command_path
            self.read_command(command_path)
        self.submitted_jobs = list()
        self.N_in_parallel = N_in_parallel
        self.wait_interval = wait_interval

        self.resubmit_seconds = resubmit_min * 60
        self.kill_seconds = kill_min * 60
        self.failed = list()

    def read_command(self, path_command):
        """
        Initializes jobs from the command file
        :param path_command:
        :return:
        """
        self.jobs = [
            Job(command.strip())
            for command in open(path_command, 'r').readlines()
        ]

    def add_jobs(self, jobs: list):
        """

        :param jobs: list of Job objects
        :return:
        """
        self.jobs = self.jobs + jobs

    def __run_jobs(self, jobs):
        """
        Submits the jobs. moves the submitted ones from self.jobs to
        self.submitted_jobs
        :param jobs: a list of jobIDs
        :return:
        """
        ids = list()
        for job in jobs:
            id = job.submit()
            if id:
                ids.append(id)
        self.submitted_jobs += [job for job in self.jobs if job.jobID in ids]
        self.jobs = [job for job in self.jobs if job.jobID not in ids]
        logger.info('submitted {} jobs'.format(len(ids)))

    def is_done(self):
        """Returns True if all jobs are done"""
        return 0 == len(self.jobs) and 0 == len(self.submitted_jobs)

    def __get_job_by_id(self, ids: list):
        return [job for job in self.submitted_jobs if job.jobID in ids]

    def long_jobs(self):
        ids = [job.jobID for job in self.submitted_jobs]
        bjobs = lsf_util.bjobs(jobID=ids)
        now = datetime.now()
        too_long = [
            job.JOBID for i, job in bjobs.iterrows()
            if (lsf_util.seconds_since_submit(job.SUBMIT_TIME, now) >
                self.resubmit_seconds)
        ]
        too_long = self.__get_job_by_id(too_long)
        for job in too_long:
            # still not resubmitted
            if -1 == job.resubmitID:
                id_ = job.submit(resubmit=True)
                if not id_:
                    logger.error('cant resubmit a long job. will b killed\n' +
                                 job.command)
                    self.__finish_job(job)
                else:
                    logger.warning('the following command was resubmitted:\n' +
                                   job.command)

            else:  # was resubmitted already
                bjob = lsf_util.bjobs(jobID=job.resubmitID)
                stat = bjob.STAT[0]
                if stat in ['DONE', 'EXIT']:
                    self.__finish_job(job, fail=False)
                    logger.info('resubmission done. killing original\n' +
                                job.command)
                else:  # resubmitted not done
                    submit_time = bjob.SUBMIT_TIME[0]
                    resubited_time = lsf_util.seconds_since_submit(
                        submit_time, now)
                    if resubited_time > self.kill_seconds:
                        logger.error('resubmitted job too long. both killed\n'
                                     + job.command)
                        self.__finish_job(job, resubmit=True)

    def __finish_job(self, job: Job, fail=True, resubmit=False):
        """Kills the job, add to fails and removes from
        submitted list
        """
        ids = [job.jobID]
        if resubmit:
            ids.append(job.resubmitID)
        lsf_util.bkill(ids)
        if fail:
            self.failed.append(job.command)
        self.__remove_from_submitted([job.jobID])

    def remove_from_submitted(self, kill_resubmit=True):
        """Removes from submitted all jobs with a status DONE or EXIT.
        :param kill_resubmit: if the jobs have a resubmitted job, it will be
        killed
        """
        ids = [job.jobID for job in self.submitted_jobs]
        print (ids)
        bjobs = lsf_util.bjobs(jobID=ids)
        ids = bjobs[(bjobs.STAT == 'DONE') | (
            bjobs.STAT == 'EXIT')].JOBID.tolist()
        if kill_resubmit:
            jobs_to_end = self.__get_job_by_id(ids)
            ids += [
                job.resubmitID for job in jobs_to_end if job.resubmitID != -1
            ]
        self.__remove_from_submitted(ids)

    def __remove_from_submitted(self, ids, kill_resubmit=True):
        """removes from the submitted job list all jobs with id in ids"""
        logger.debug('removes: %i', len(ids))
        logger.debug('submitted_jobs before: ' + str(len(self.submitted_jobs)))
        self.submitted_jobs = [
            job for job in self.submitted_jobs if job.jobID not in ids
        ]
        logger.debug('submitted_jobs after: ' + str(len(self.submitted_jobs)))

    def run_jobs(self):
        while 0 != len(self.jobs):
            self.remove_from_submitted()
            self.long_jobs()
            #lsf_util.average_queue_pends()

            N = len(self.jobs)
            to_index = self.N_in_parallel - len(self.submitted_jobs)
            to_index = to_index if to_index < N else N

            self.__run_jobs(self.jobs[:to_index])
            sleep(self.wait_interval * 60)

        logger.info('All jobs were submitted')

        while not self.is_done():
            self.remove_from_submitted()
            #lsf_util.average_queue_pends()
            self.long_jobs()
            print (self.wait_interval)
            sleep(self.wait_interval * 60)

        logger.info('Finished running Rosetta modeling trajectories')

    def save_failes(self, path):
        """prints to a file all failing commands"""
        open(path, 'w').writelines(self.failed)
