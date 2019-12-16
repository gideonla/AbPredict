#!/bin/env python3


import os
import shutil, errno


def create_dir(path):
    """Creates the directory in path
    :param path: path to the new folder to create
    :return: True if created a new directory
    """
    if not os.path.exists(path):
        os.mkdir(path)
        return True
    else:
        return False


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


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise


def removeanything(path):
    try:
        shutil.rmtree(path)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            os.remove('a')
        else:
            raise


