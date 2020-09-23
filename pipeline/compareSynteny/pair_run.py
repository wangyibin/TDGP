#!/usr/bin/env python
# -*- coding:utf-8 -*-


"""
%prog sample.list [Options]
    to pairs multi samples and run haplotype compartment compare.
"""

import argparse
import os
import os.path as op
import sys

from itertools import combinations
from TDGP.apps.grid import Cluster

def pairs(samples):
    samples = list(combinations(samples, 2))

    return list(map(lambda x:"{}-{}".format(x[0], x[1]), samples))


def run_pairs(pairs):
    sample1, sample2 = pairs.split("-")
    sample1_h = sample1.replace('H', '')
    sample2_h = sample2.replace('H', '')

    os.mkdir(pairs)
    os.chdir(pairs)
    header = Cluster().header
    os.system('ln -s ../data/{}* .'.format(sample1))
    os.system('ln -s ../data/{}* .'.format(sample2))
    os.system('ln -s ../data/{0}-{1}/{{{0},{1}}}.synteny.bed .'.format(sample1_h, sample2_h))
    script = 'run_{}.sh'.format(pairs)
    with open(script, 'w') as out:
        cmd = header + "\n"
        cmd += 'python -m TDGP.analysis.ab getSyntenyGenePca {0}_*.bg {1}_*.bg {2}.synteny.bed {3}.synteny.bed --plot\n'.format(sample1, sample2, sample1_h, sample2_h)
        cmd += 'cut -f 1-3,5 {0}.synteny.eigen1.bg > {0}.synteny.eigen1.nogene.bg\n'.format(sample1_h)
        cmd += 'cut -f 1-3,5 {0}.synteny.eigen1.bg > {0}.synteny.eigen1.nogene.bg\n'.format(sample2_h)
        cmd += 'python ../plotLRPerchrom_haplotype.py {0}.synteny.eigen1.nogene.bg {1}.synteny.eigen1.nogene.bg --xlabel {2} --ylabel {3}\n'.format(sample1_h, sample2_h, sample1, sample2)
        
        out.write(cmd)
    os.system('qsub {}'.format(script))
    os.chdir('../')

    

def main(args):
    p = p=argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-s', '--samples', nargs="+", 
            help='sample of pairs', required=True)
    pOpt.add_argument('-r', '--rerun', action='store_true', default=False,
            help='rerun all sample pairs [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    samples_pairs = pairs(args.samples)
    for pair in samples_pairs:

        if args.rerun:
            os.system('rm -rf {}'.format(pair))
        else:
            if op.exists(pair):
                continue

        run_pairs(pair)


if __name__ == "__main__":
    main(sys.argv[1:])