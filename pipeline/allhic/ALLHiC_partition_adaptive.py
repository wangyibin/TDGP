#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

"""

import argparse
import logging
import os
import os.path as op
import sys

import glob
import multiprocessing 
import tempfile

from collections import OrderedDict 


ENZYME_SITES = {
    'HINDIII': 'AAGCTT',
    'MBOI': 'GATC',
    'Arima': 'GATC,GANTC',
    'DPNI': 'GATC' 
}

def run_cmd(cmd):
    logging.info("CMD: " + cmd)
    os.system(cmd)

def extract(bam, ref, esite):
     
    extract_cmd = "allhic extract {} {} --RE {} 2>extract.log".format(bam, ref, esite)

    run_cmd(extract_cmd)

    counts_file = "{}.counts_{}.txt".format(bam.replace('.bam', ''), 
                                            esite)
    pairs_file = "{}.pairs.txt".format(bam.replace('.bam', ''))

    return counts_file, pairs_file 

def is_complete_partition(counts_file, k, targetDir='./'):
    """
    complete partition is that partition k groups, instead of k+.
    """
    prefix = counts_file.replace('.txt', '')
    
    targets = glob.glob('{}/{}.{}g*txt'.format(targetDir, prefix, k))
    print(targets)
    if len(targets) == k:
        return True
    else:
        return False



def partition(counts_file, paris_file, k, 
            minREs=25, maxLinkDensity=2):
    
    partition_cmd = 'allhic partition {} {} {} --minREs {} --maxLinkDensity{}'
    partition_cmd = partition_cmd.format(counts_file, paris_file, k, 
                                        minREs, maxLinkDensity)
    run_cmd(partition_cmd)

    
def main(args):
    
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-b', '--bamFile', 
            help='prunned bam', required=True)
    pReq.add_argument('-r', '--reference', 
            help='draft.asm.fasta', required=True)
    pReq.add_argument('-k', required=True, type=int,
            help='number of groups (user defined K value)')
    pOpt.add_argument('-e', '--enzyme', default='HindIII',
            help='enzyme sites (HindIII: AAGCTT; MboI: GATC)')
    pOpt.add_argument('--tmp_dir', default='ALLHiC_tmp', 
            help='temporary directory [default: $(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    bam = args.bamFile
    ref = args.reference
    k = args.k
    enzyme = args.enzyme
    try:
        esite = ENZYME_SITES[enzyme.upper()]
    except KeyError:
        esite = enzyme 
    
    counts_file, paris_file = extract(bam, ref, esite)
    print(is_complete_partition(counts_file, k))
    with tempfile.TemporaryDirectory(prefix=args.tmp_dir, dir='./') as tmpDir:
        
        pass


    
    

if __name__ == "__main__":
    main(sys.argv[1:])