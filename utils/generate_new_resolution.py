#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%(prog)s 
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys


from TDGP.apps.base import debug, check_file_exists
from TDGP.apps.grid import parallel
debug()

#SCRIPT_PATH = "/public1/home/stu_wangyibin/software/HiC-Pro_2.11.1/scripts/"
#PYTHON_PATH = "/public1/home/stu_wangyibin/software/anaconda2/bin/python"

def build_matrix(allValidPairs, 
                    binsize, 
                    chrsizes, 
                    oprefix_path):
    
    oprefix = "{0}/raw/{1}/{2}_{1}".format(oprefix_path, binsize,
                                        op.basename(oprefix_path))
    

    cmd = "cat {} | build_matrix --matrix-format upper \
            --binsize {} --chrsizes {} --ifile /dev/stdin \
                --oprefix {}".format(allValidPairs, 
                    binsize, chrsizes, oprefix)
    logging.debug("Starting to execute build_matrix in `{}`".format(binsize))
    retcode =  os.system(cmd)
    return retcode

def iced_norm(oprefix_path, binsize):
    iprefix = "{0}/raw/{1}/{2}_{1}.matrix".format(oprefix_path, binsize,
                                        op.basename(oprefix_path))
    oprefix = "{0}/iced/{1}/{2}_{1}_iced.matrix".format(oprefix_path, binsize,
                                        op.basename(oprefix_path))
    cmd = "ice --results_filename {} --filter_low_counts_perc 0.02 \
        --filter_high_counts_perc 0 --max_iter 100 --eps 0.1 --remove-all-zeros-loci \
            --output-bias 1 --verbose 1 {}".format(oprefix, iprefix)
    logging.debug("Starting to execute iced_norm in `{}`".format(binsize))
    retcode = os.system(cmd)


def main(args):
    allValidPairs, binsize, chrsizes, oprefix_path = args
    if oprefix_path.endswith('/'):
        oprefix_path = op.split(oprefix_path)[0]
    try:
        os.makedirs("{}/raw/{}".format(oprefix_path, binsize))
    except OSError:
        pass
    try:
        os.makedirs("{}/iced/{}".format(oprefix_path, binsize))
    except OSError:
        pass
    
    retcode = build_matrix(allValidPairs, binsize, 
                                chrsizes, oprefix_path)
    if retcode == 1:
        logging.error('failed to execute build_matrix, please check `{}`'.format(binsize))
        sys.exit()
    else:
        logging.debug('build_matrix done `{}`'.format(binsize))
    retcode = iced_norm(oprefix_path, binsize)
    if retcode == 1:
        logging.error('failed to execute ice_norm, please check `{}`'.format(binsize))
        sys.exit()
    else:
        logging.debug('iced_norm done `{}/{}`'.format(oprefix_path, binsize))
    logging.debug('Successful to generate matrix in `{}/{}`'.format(oprefix_path, binsize))


if __name__ == "__main__":
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-i', '--input', required=True, 
            help='allValidParis file')
    pReq.add_argument('-b', '--binsize', nargs="+", required=True,
            help='binsize of matrix')
    pReq.add_argument('-c', '--chrsizes', required=True,
            help='chromsizes file')
    pReq.add_argument('-o', '--oprefix', required=True,
            help='out prefix of matrix')
    pOpt.add_argument('-t', '--thread', type=int, default=4, 
            help='the thread of program [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args()

    tasks = [(args.input, binsize, args.chrsizes, args.oprefix) 
            for binsize in args.binsize]
    
    parallel(main, tasks, args.thread)
