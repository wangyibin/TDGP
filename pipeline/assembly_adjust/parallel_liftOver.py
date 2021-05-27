#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
parallel run liftOverValidPairs by snakemake
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys


def parallel_liftOver(args):
    """
    %(prog)s 
    """
    p = argparse.ArgumentParser(prog=parallel_liftOver.__name__,
                        description=parallel_liftOver.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('validpairs', 
            help='validpairs from hicpro')
    pReq.add_argument('old_agp', help='old agp file')
    pReq.add_argument('new_agp', help='new agp file')
    
    pOpt.add_argument('--nosplit', action='store_true', default=False, 
            help='Do not split validpairs [default: %(default)s]')
    pOpt.add_argument('-t', '--threads', type=int, default=20,
            help='threads number of each script [default: %(default)s]')
    pOpt.add_argument('-j', '--jobs', type=int, default=10, 
            help='maximum jobs of snakemake submit to cluster')
    pOpt.add_argument('--cluster', default='qsub -l nodes=1:ppn={threads}'
            ' -j oe -q workq -V', help='snakemake cluster command [default: %(default)s]')
    pOpt.add_argument('-s', '--snakefile', default=None, 
            help='snakefile of parallel_liftOver.smk [default: auto find]')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    path = op.dirname(op.realpath(__file__))

    sample = op.basename(args.validpairs).split(".")[0]
     ## split validpairs
    if not args.nosplit:
        cmd = "split -n -a 4 -l {} data/{}".format(1000000, sample)
        #os.system(cmd)
        print(cmd)
    ## run snakemake
    if not args.snakefile:
        snakefile = op.join(path, "parallel_liftOver.smk")
    else:
        snakefile = args.snakefile
    
    if not op.exists(snakefile):
        logging.error("No such snakefile of `{}`".format(snakefile))
        sys.exit()
    cmd = "snakemake -s {} --config sample={} old_agp={} \
        new_agp={} ncpus={} -j {} "
    cmd = cmd.format(snakefile, sample, args.old_agp,
                    args.new_agp, args.threads, args.jobs)
    cmd += "--jobname liftOver.{rulename}.{jobid}.sh "
    if args.cluster:
        cmd += '--cluster "{}"'.format(args.cluster)
    
   
    
    print(cmd)
    os.system(cmd)


if __name__ == "__main__":
    parallel_liftOver(sys.argv[1:])