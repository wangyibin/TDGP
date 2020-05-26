#!/usr/bin env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import os
import os.path as op
import logging
import sys

from multiprocessing import Pool, Process, cpu_count
from TDGP.apps.base import debug,listify


debug()

def main():

    actions = (
            ("parafly", ""),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


class Jobs(object):
    def __init__(self):
        pass

class Parallel(object):
    """
    Run commands in parallel.
    """
    def __init__(self, target, args, threads=cpu_count()):
        self.target = target
        self.args = args
        self.threads = min(len(args), threads)

    def run(self):
        p = Pool(self.threads)
        res = p.map(self.target, self.args)
        return res


def parallel(target, args, threads=cpu_count):
    p = Pool(min(len(args), threads))
    res = p.map(target, args)
    return res



class CMD(object):
    """
    Linux command execute object

    Params:
    -------
    """

    def __init__(self, cmds, threads=4):

        self.cmds = listify(cmds)
        self.threads = threads
        self.run()
    
    def run(self):
        p = Parallel(os.system, self.cmds, self.threads)
        p.run()



PBS_HEADER = """#!/bin/bash
#PBS -m ae
#PBS -j eo
#PBS -N {}
#PBS -q {}
#PBS -V 
#PBS -l nodes=1:ppn={}
cd $PBS_O_WORKDIR
"""

SGE_HEADER = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N {}
#$ -q {}
#$ -pe mpi {}
"""
class Cluster(object):
    """
    class of cluster operation
        in order to execute successful should set 
                the `CLUSTER` variants into ENV
    Params:
    --------
    
    Returns:
    ---------
    out: `str`: CLUSTER

    Functions:
    ---------
    get_header: get the header of cluster system with parameters
    get_raw_header: get the raw of cluster system with parameters

    """
    
    def __init__(self):
        self.CLUSTER = 'SGE'
        self.get()
        self.get_header()
        self.get_raw_header()

    def get(self):
        """
        To obtain the environment of `CLUSTER`,
            if not found will be set default `SGE`.
        """
        try:
            self.CLUSTER = os.environ['CLUSTER']
        except KeyError:
            logging.warning('There is not environment `CLUSTER` in PATH')
        return self.CLUSTER


    def get_header(self, name='jobs', queue='all.q', threads=1):
        """
        According to the environment of `CLUSTER` to 
            return a header of cluster system
        """
        if self.CLUSTER.upper() == "SGE":
            self.header = SGE_HEADER.format(name, queue, threads)
            
        elif self.CLUSTER.upper() == "PBS":
            queue = 'workq' if queue == 'all.q' else queue
            self.header = PBS_HEADER.format(name, queue, threads)
        else:
            logging.warning("there is not of header "
                            "of cluster:`{}`".format(self.cluster))
        return self.header

    def get_raw_header(self, name='jobs', queue='all.q', threads=1):
        """
        According to the environment of `CLUSTER` to 
            return a header of cluster system
        """
        if self.CLUSTER.upper() == "SGE":
            self.raw_header = SGE_HEADER
        elif self.CLUSTER.upper() == "PBS":
            self.raw_header = PBS_HEADER
        else:
            logging.warning("there is not of header "
                            "of cluster:`{}`".format(self.cluster))
        return self.raw_header


    def __str__(self):
        return self.CLUSTER

    __retr__ = __str__
        
        


    