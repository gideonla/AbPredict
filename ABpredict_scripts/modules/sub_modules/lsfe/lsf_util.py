#!/bin/env python3
import subprocess as sp
import pandas as pd
from math import ceil
from getpass import getuser
import logging
from datetime import datetime


USER = getuser()
QUEUES = ['new-short', 'Fleishman-service', 'new-medium', 'new-all.q', 'new-long']
logger = logging.getLogger(__name__)


def bjobs(jobID='', queue='', user=USER, w=False) -> pd.DataFrame:
    """
    :param jobID: a string or a list of strings of jobIDs
    :return: The bjobs output as a DataFrame with the following columns:
    JOBID   USER    STAT  QUEUE   FROM_HOST   EXEC_HOST  JOB_NAME SUBMIT_TIME
    """
    command = ['bjobs', '-u', user]
    if w:
        command.append('-w')
    if queue:
        command.append('-q')
        command.append(queue)
    if jobID:
        if str == type(jobID):
            command.append(jobID)
        else:
            command += jobID 
    process = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE,
                       universal_newlines=True)
    # error = process.stderr.read()
    output = process.stdout.readlines()
    if not output:
        error = process.stderr.read()
        logger.debug(error)
        bjobs_o = pd.DataFrame(columns=['JOBID', 'USER', 'STAT', 'QUEUE',
                                        'FROM_HOST', 'EXEC_HOST', 'JOB_NAME',
                                        'SUBMIT_TIME'])
    else:
        # jobs = [job.split()[:7] + [' '.join(job.split()[7:])]
        #         for job in output[1:]]
        jobs = list()
        for job in output[1:]:
            jobs.append((job[:8].strip(),    # jobID
                        job[8:16].strip(),  # user
                        job[16:22].strip(), # stat
                        job[22:33].strip(), # queue
                        job[33:45].strip(), # from_host
                        job[45:57].strip(), # exec_host
                        job[57:68].strip(), # job_name
                        job[68:].strip()))  # submit_time
        bjobs_o = pd.DataFrame(data=jobs, columns=output[0].split())
    return bjobs_o


def bqueues(queue='') -> pd.DataFrame:
    """

    :param queue:
    :return: the bqueues output as a DataFrame with the columns:
    QUEUE_NAME  PRIO STATUS  MAX JL/U JL/P JL/H NJOBS  PEND   RUN  SUSP
    """
    command = ['bqueues', '-u', USER]
    if queue:
        command.append(queue)
    process = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE,
                       universal_newlines=True)
    output = process.stdout.readlines()
    error = process.stderr.read()
    assert 'No such queue' not in error, \
        'The queue {} doesn\' exists in bqueues'.format(queue)
    queues = pd.DataFrame(data=[q.split() for q in output[1:]],
                          columns=output[0].split())
    return queues


def bswitch(queue: str, jobs: list) -> list:
    """

    :param queue:
    :param jobs: list od jobIDs
    :return:
    """
    switched = list()
    for jobID in jobs:
        command = ['bswitch', queue, jobID]
        process = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE,
                           universal_newlines=True)
        output, error = process.communicate()
        if 'is switched' in output:
            switched.append(jobID)
    return switched


def bkill(jobID='', user=USER):
    """TODO: not tested for failing to kill"""
    command = ['bkill', '-u', USER]    
    if jobID:
        if str == type(jobID):
            command.append(jobID)
        else:
            command += jobID 
    process = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE,
                       universal_newlines=True)


def average_queue_pends():
    jobs = bjobs()
    if 0 == len(jobs):
        return
    jobs_no_run = jobs[jobs.STAT != 'RUN']
    pend_average = ceil(len(jobs_no_run) / len(QUEUES))
    if 0 == pend_average:
        logger.info('all jobs are running!')
        return
    free_q = jobs_no_run.groupby(['QUEUE']).size() - pend_average
    
    # adding all queues
    for q in QUEUES:
        if q not in free_q.keys():
            free_q = free_q.set_value(q, -pend_average)    
    
    free_q.sort_values(inplace=True, ascending=False)


    logger.info('pend average ' + str(pend_average))

    to_switch = list()
    for q, free in free_q.iteritems():
        logger.debug('{} {}'.format(q, free))
        logger.debug('to switch len before: {}'.format(len(to_switch)))
        if free >= 0:
            free = max(free, 1)
            to_switch += jobs_no_run[jobs_no_run.QUEUE == q].JOBID.tolist()[:free]
        else: # queue to switch to
            switch_to_q = to_switch[:-free]
            switched = bswitch(queue=q, jobs=switch_to_q)
            to_switch = list(set(to_switch) - set(switched))
        logger.debug('to switch len after: {}'.format(len(to_switch)))


def seconds_since_submit(submit_time: str, now=None):
    """
    :param submit_time: the string the the submit time column from bjobs
    :param now: a datetime object if want to check time between different time
    than the current one 
    :returns: the number of seconds since the job was submitted
    """
    if not now:
        now = datetime.now()
    submit_time = datetime.strptime(submit_time, '%b %d %H:%M')
    
    # update year of submit_time
    year = now.year if now.month >= submit_time.month else now.year-1
    submit_time = submit_time.replace(year=year)
    
    elapsed = now - submit_time 
    return elapsed.total_seconds()
