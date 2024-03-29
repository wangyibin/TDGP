#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Fastq format and tools.
"""

from __future__ import print_function

import argparse
import logging
import numpy as np
import multiprocess
import gzip
import os
import os.path as op
import random
import re
import shutil
import subprocess
import sys

from glob import glob

from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import check_file_exists, debug
from TDGP.apps.base import listify
from TDGP.apps.grid import parallel


def main():
    actions = (
        ('test', "test"), 
        ('splitFastq', "split fastq to several parts")
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def test():
    print('test')

def splitfastq(infile, outdir, nreads):
    """
    split fastq by reads number
    """
    nlines = int(nreads) * 4

    try:
        os.makedirs(outdir)
    except OSError:
        pass

    if infile.endswith('.gz'):
        prefix = re.sub('((.fastq)|(.fq)).gz', '_part', 
                op.join(outdir, op.basename(infile)))
        cmd = 'zcat {} | split -a 3 -l {} -d - {}'.format(
                infile, nlines, prefix)
    else:
        prefix = re.sub('((.fastq)|(.fq))', '_part',
                op.join(outdir, op.basename(infile)))
        cmd = 'split -a 3 -l {} -d {} {}'.format(
                nlines, infile, prefix)

    retcode = subprocess.call(cmd, shell=True)
    if retcode !=0:
        logging.error('Failed to split infile with return code {}'.format(retcode))
        sys.exit(1)
    
    files = glob(prefix + "*")
    files.sort()
    res = []
    for i in files:
        shutil.move(i, op.join(
                op.dirname(i),
                op.basename(i)[-3:] + "_" +
                op.basename(i)[:-8] + ".fastq"))
        res.append(op.join(
                op.dirname(i),
                op.basename(i)[-3:] + "_" +
                op.basename(i)[:-8] + ".fastq"))
    #res = glob(op.dirname(prefix) + "/*fastq")
    return res
    
## out command

def downsample(args):
    """
    %(prog)s R1 R2 [Options]
        random select fastq to downsample
    """
    p = argparse.ArgumentParser(prog=downsample.__name__,
                        description=downsample.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-1', required=True,
            help='left fastq file')
    pReq.add_argument('-2', required=True,
            help='right fastq file')
    pOpt.add_argument('-n', type=int, 
            help='number of reads each fastq to select')
    pOpt.add_argument('-s', type=float,
            help='size of data each fastq')
    pOpt.add_argument('-l', '--length', default=150, type=int,
            help='length of reads [default: %(default)s]')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    
    
def mv(x):
    shutil.move(x, "./")    

def splitFastq(args):
    """
    %(prog)s R1 R2 
    """
    p = argparse.ArgumentParser(prog=splitFastq.__name__,
                        description=splitFastq.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fastq', nargs="+", help='fastq file')
    pOpt.add_argument('-n', '--nreads', default=1e8, type=float,
            help='reads number of per part [default: %(default)s')
    pOpt.add_argument('-t', '--threads', type=int, default=12,
            help='the thread numbers of program [default: %(default)s]')
    pOpt.add_argument('--decompress', action='store_true', default=False, 
            help='decompress output result [default: %(default)s]')

    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    
    outs = list(map(lambda x: 
            re.sub('((.fastq.gz)|(.fq.gz)|(.fastq)|(.fq))', '_out', x), 
            args.fastq))
    thread = min(len(args.fastq), args.threads)
    pool = multiprocess.Pool(thread)
    logging.debug('Starting ...')
    def splitfq(arg_list):
        fq, outdir, nreads = arg_list
        files = splitfastq(fq, outdir, nreads)
        return files
    task_list = [(fq, out, args.nreads) for (fq, out) in zip(args.fastq, outs)]
    res = pool.map(splitfq, task_list)
    res = np.append([], res)
    cmd_func = lambda x: "gzip {0}".format(x)
    gzip_cmd = '{}_gzip_cmd.list'.format(os.getpid())
    with open(gzip_cmd, 'w') as fo:
        print("\n".join(map(cmd_func, res)), file=fo)
    
    if not args.decompress:
        cmd = "cat {} | parallel -j {} {{}}".format(gzip_cmd, args.threads)
        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            logging.error('Failed to excude command of `{}`'.format(cmd))
        else:
            logging.debug('Successful compress all fastq file.')
        res = list(map(lambda x: x + ".gz", res))
    
    
    parallel(mv, list(res), args.threads)
    os.system("rm -rf {}".format(" ".join(outs)))
    logging.debug('Done')



if __name__ == "__main__":
    main()