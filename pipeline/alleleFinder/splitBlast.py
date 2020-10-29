#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
parallel run blast by snakemake
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

def splitBlast(args):
    """
    %(prog)s 
    """
    p = p=argparse.ArgumentParser(prog=splitBlast.__name__,
                        description=splitBlast.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pSmk = p.add_argument_group('Snakemake arguments')
    pReq.add_argument('-d', '--database', required=True, 
            help='database fasta')
    pReq.add_argument('-q', '--query', required=True,
            help='query fasta')
    pOpt.add_argument('-dbtype', default='nucl', choices=['nucl', 'prot'],
            help='dtype of makeblastdb [default: %(default)s]')
    pOpt.add_argument('-nparts', default=10, type=int, 
            help='part number of query split [default: %(default)s]')
    pOpt.add_argument('-o', '--output', default='all.blast.out', 
            help='output of all blast results [default: %(default)s]')
    pOpt.add_argument('-evalue', default='1e-5', 
            help='evalue of blast [default: %(default)s]')
    pOpt.add_argument('-outfmt', default='6', 
            help='outfmt of blast [default: %(default)s]')
    pOpt.add_argument('-num_alignments', default='5', 
            help='num_alignments of blast [default: %(default)s]')
    pOpt.add_argument('-t', '--threads', type=int, default=4,
            help='threads number of each script [default: %(default)s]')
    pSmk.add_argument('-j', '--jobs', type=int, default=5, 
            help='maximum jobs of snakemake submit to cluster')
    pSmk.add_argument('--cluster', default='qsub -l nodes=1:ppn={threads}'
            ' -j oe -q workq -V', help='snakemake cluster command \n[default: %(default)s]')
    pSmk.add_argument('-s', '--snakefile', default=None, 
            help='snakefile of parallel_liftOver.smk [default: auto find]')
    # pSmk.add_argument('--rerun-incomplete', dest='rerun_incomplete', 
    #         action='store_true', default=False,
    #         help='rerun incomplete jobs [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    path = op.dirname(op.realpath(__file__))

    database = args.database
    query = args.query

    ## run snakemake
    if not args.snakefile:
        snakefile = op.join(path, "splitBlast.smk")
    else:
        snakefile = args.snakefile
    

    if not op.exists(snakefile):
        logging.error("No such snakefile of `{}`".format(snakefile))
        sys.exit()
    cmd = "snakemake -s {} --config fasta={} query={} dbtype={} \
        nparts={} ncpus={} out={} evalue={} outfmt={} \
            num_alignments={} -j {}  --rerun-incomplete "
    cmd = cmd.format(snakefile, database, query, args.dbtype,
                    args.nparts, args.threads, args.output, args.evalue, 
                    args.outfmt, args.num_alignments, args.jobs)
    cmd += "--jobname splitBlast.{rulename}.{jobid}.sh "
    if args.cluster:
        cmd += '--cluster "{}"'.format(args.cluster)


    print(cmd, file=sys.stderr)
    os.system(cmd)


if __name__ == "__main__":
    splitBlast(sys.argv[1:])