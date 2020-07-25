#!/usr/bin/env python

from __future__ import print_function


import gzip
import os.path as op
import sys
import random

from Bio import SeqIO
def downsample(args):
    """
    %(prog)s R1 R2 [Options]
        random select fastq to downsample
    """
    p = p=argparse.ArgumentParser(prog=downsample.__name__,
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
    pOpt.add_argument('--seed', type=int, default=12345,
            help='random seed [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    

    def random_selecet(a):
        fastq, number = a
        length = int(os.pipe("zcat {} | parallel --pipe wc -l | awk '{i+=$1}END{print i}'".format(fastq)))

        
        fastq_handle = SeqIO(gzip.open(fastq), 'fastq')
        for record in fastq_handle:
            num = random.randint(length)
            if num < number:
                print(record)
                break 
    
    
        

    