#!/usr/bin env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import argparse
import os
import os.path as op
import logging
import sys

from multiprocessing import Pool, Process, cpu_count
from TDGP.apps.base import debug, listify, ActionDispatcher


debug()

def main():

    actions = (
            ("clusterHeader", "print header of cluster system"),
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
#PBS -j eo {}
#PBS -q {}
#PBS -V 
#PBS -l nodes=1:ppn={} {}
if [[ ! -z $PBS_O_WORKDIR ]]; then
    cd $PBS_O_WORKDIR
fi
"""

SGE_HEADER = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd {}
#$ -q {}
#$ -pe mpi {} {}
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
    
    def __init__(self, cluster=None, 
                    name=None, queue=None, 
                    threads=1, array=None):
        self.CLUSTER = cluster if cluster else None
        if not self.CLUSTER:
            self.get()
        self.get_header(name, queue, threads, array)
        self.get_raw_header()

    def get(self):
        """
        To obtain the environment of `CLUSTER`,
            if not found will be set default `SGE`.
        """
        try:
            self.CLUSTER = os.environ['CLUSTER']
        except KeyError:
            self.CLUSTER = 'SGE'
            logging.warning('There is not environment `CLUSTER` in PATH')

        return self.CLUSTER


    def get_header(self, name=None, queue=None, 
                        threads=1, array=None):
        """
        According to the environment of `CLUSTER` to 
            return a header of cluster system
        """
        if self.CLUSTER.upper() == "SGE":
            name = "\n#$ -N " + name  if name else ""
            queue = queue if queue else "all.q"
            array = "\n#$ -t " + array if array else ""
            self.header = SGE_HEADER.format(name, queue, threads, array)   
        elif self.CLUSTER.upper() == "PBS":
            name = "\n#PBS -N " + name if name else ""
            queue = queue if queue else "workq"
            array = "\n#PBS -J " + array if array else ""
            self.header = PBS_HEADER.format(name, queue, threads, array)

        else:
            logging.warning("there is not of header "
                            "of cluster:`{}`".format(self.CLUSTER))
            sys.exit()
        return self.header

    def get_raw_header(self):
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
                            "of cluster:`{}`".format(self.CLUSTER))
            sys.exit()
        return self.raw_header


    def __str__(self):
        return self.CLUSTER

    __retr__ = __str__
        
### out command ###
def clusterHeader(args):
    """
    %(prog)s 
    print the header of clustes
    """    
    p = p=argparse.ArgumentParser(prog=clusterHeader.__name__,
                        description=clusterHeader.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pOpt.add_argument('-c', '--cluster', default=None, 
            help='cluster system [default: auto]')
    pOpt.add_argument('-n', '--name', default=None,
            help='name of jobs in cluster [default: jobs name]')
    pOpt.add_argument('-q', '--queue', default=None, 
            help='queue of cluster [default: auto]')
    pOpt.add_argument('-t', '--threads', default=1, type=int,
            help='threads number of program [default: %(default)s]')
    pOpt.add_argument('-a', '--array', default=None, 
            help='array jobs [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    cluster = Cluster(args.cluster, args.name, 
                             args.queue, args.threads,
                            args.array)
    print(cluster.header, file=sys.stdout)




if __name__ == "__main__":
    main()