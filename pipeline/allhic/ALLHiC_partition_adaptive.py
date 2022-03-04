#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

"""

import argparse
import logging
import os
import os.path as op
import sys
import shutil

import glob
import multiprocessing 
import tempfile

import pandas as pd

from collections import OrderedDict 
from joblib import Parallel, delayed
from rich.logging import Console, RichHandler


logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=Console(stderr=True))]
)
logger = logging.getLogger(__name__)

ENZYME_SITES = {
    'HINDIII': 'AAGCTT',
    'MBOI': 'GATC',
    'Arima': 'GATC,GANTC',
    'DPNI': 'GATC' 
}

def run_cmd(cmd):
    logger.info("Running CMD: " + cmd)
    os.system(cmd)

def cmd_exists(program):
    """
    Check program is exists.
    """
    from shutil import which
    return which(program) is not None

def extract(bam, ref, esite):
     
    extract_cmd = "allhic extract {} {} --RE {} 2>extract.log".format(
                                                        bam, ref, esite)
    run_cmd(extract_cmd)

    counts_file = "{}.counts_{}.txt".format(bam.replace('.bam', ''), 
                                            esite)
    pairs_file = "{}.pairs.txt".format(bam.replace('.bam', ''))

    return counts_file, pairs_file 

def get_partition_res(counts_file, k, targetDir):

    prefix = counts_file.replace('.txt', '')
    targets = glob.glob('{}/{}.{}g*txt'.format(targetDir, prefix, k))
    targets = sorted(targets)
    return targets

def is_complete_partition(targets, k, targetDir='./'):
    """
    complete partition is that partition k groups, instead of k+.
    """
    
    if len(targets) == k:
        return True
    else:
        return False

def import_countRE(count_RE):
    df = pd.read_csv(count_RE, sep='\t', header=0, index_col=None)

    return df 

def partition(counts_file, pairs_file, k, 
            minREs=25, maxLinkDensity=2):
    
    partition_cmd = ('allhic partition {} {} {} --minREs {} '
                '--maxLinkDensity {} 1>/dev/null 2>/dev/null')
    partition_cmd = partition_cmd.format(counts_file, pairs_file, k, 
                                        minREs, maxLinkDensity)
    run_cmd(partition_cmd)

    targets = get_partition_res(counts_file, k, "./")
    
    return targets

def calc_length(targets):
    data = list(map( lambda x: import_countRE(x)['Length'].sum(), targets))
    return data

def calc_range(data):
    range = max(data) - min(data)
    return range

def find_best_partition(counts_file, pairs_file, k, minREs, maxLinkDensity):
    workdir = "{}_{}".format(minREs, maxLinkDensity)
    os.makedirs(workdir)
    os.chdir(workdir)
    os.link("../../" + counts_file, "./" + counts_file)
    os.link("../../" + pairs_file, "./" + pairs_file)
    targets = partition(counts_file, pairs_file, k, minREs, maxLinkDensity)
    if is_complete_partition(targets, k):
        flag = "pass"
    else:
        flag = 'fail'
        os.chdir("../")
        shutil.rmtree(workdir)
        return None
    lengths = calc_length(targets)
    res = calc_range(lengths)
    result = (minREs, maxLinkDensity, ",".join(map(str, lengths)), res, 
                    res/sum(lengths))
    os.chdir("../")

    return result

def main(args):
    
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-c', '--counts', 
            help='counts txt', required=True)
    pReq.add_argument('-p', '--pairs', 
            help='prunned pairs txt', required=True)
    # pReq.add_argument('-r', '--reference', 
    #         help='draft.asm.fasta', required=True)
    pReq.add_argument('-k', required=True, type=int,
            help='number of groups (user defined K value)')
    pOpt.add_argument('--minREs', default='50,300', 
            help='comma seperated of minREs ranges [default: %(default)s]')
    pOpt.add_argument('--minREs_step', default=5, type=int,
            help='step of minREs [default: %(default)s]')
    pOpt.add_argument('--maxLinkDensity', default='2,10', 
            help='comma seperated of maxLinkDensity ranges [default: %(default)s]')
    pOpt.add_argument('--maxLinkDensity_step', default=2, type=int,
            help='step of maxLinkDensity [default: %(default)s]')
    # pOpt.add_argument('-e', '--enzyme', default='HindIII',
    #         help='enzyme sites (HindIII: AAGCTT; MboI: GATC)')
    pOpt.add_argument('--tmp_dir', default='ALLHiC_tmp', 
            help='temporary directory [default: $(default)s]')
    pOpt.add_argument('-t', '--threads', type=int, default=8,
            help='number of program threads[default:%(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    counts_file = args.counts
    pairs_file = args.pairs
    k = args.k
    minREs_tuple = list(map(int, args.minREs.split(",")))
    maxLinkDensity_tuple = list(map(int, args.maxLinkDensity.split(",")))
    minREs_step = args.minREs_step
    maxLinkDensity_step = args.maxLinkDensity_step

    if not cmd_exists('allhic'):
        logger.error('Command not found of `allhic`')
        sys.exit()

    prefix = counts_file.replace(".txt", "")
    total_length = import_countRE(counts_file)['Length'].sum()

    results = []
    with tempfile.TemporaryDirectory(prefix=args.tmp_dir, dir='./') as tmpDir:
        logger.info('Working on temporary directory: {}'.format(tmpDir))
        os.chdir(tmpDir)
        
        params = []
        for minREs in range(minREs_tuple[0], minREs_tuple[1], minREs_step):
            for maxLinkDensity in range(maxLinkDensity_tuple[0], 
                                        maxLinkDensity_tuple[1], 
                                        maxLinkDensity_step):
                params.append((counts_file, pairs_file, k, minREs, maxLinkDensity))
        
        logger.info("Finding the best partition result ...")
        results = Parallel(n_jobs=args.threads)(delayed(
                                        find_best_partition)(c, p, g, r, l)
                                            for c, p, g, r, l in params)

        results = filter(lambda x: x is not None, results)
        results = sorted(results, key=lambda x: x[4])
        try:
            best = results[0]
        except IndexError:
            logger.error("Couldn't found best results")
            return
        logger.info('Best result is [{}]'.format(best))

        best_targetDir = "./{}_{}".format(best[0], best[1]) 
        best_targets = get_partition_res(counts_file, k, best_targetDir)
        os.system("cp {}_{}/*{}g*.txt ../".format(best[0], best[1], k))
        os.chdir("../")
        logger.info("Removed temporary directory.")
        
    ## output partition result table
    output_table = "{}.partition.table".format(prefix)
    with open(output_table, 'w') as out:
        print("\t".join(['#minREs', 'maxLinkDensity', 'groups',
                            'range', 'value']), file=out)
        for res in results:
            print("\t".join(map(str, res)), file=out)
    
    logger.info("Done.")

if __name__ == "__main__":
    main(sys.argv[1:])