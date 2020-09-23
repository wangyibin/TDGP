#!/usr/bin/env python
# -*- coding:utf-8 -*-


"""
%prog sample.list [Options]
    to pairs multi samples and run jcvi mcscan.
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


def run_synteny(pairs):
    sample1, sample2 = pairs.split("-")
    os.mkdir(pairs)
    os.chdir(pairs)
    header = Cluster().header
    os.system('ln -s ../data/{}* .'.format(sample1))
    os.system('ln -s ../data/{}* .'.format(sample2))
    script = 'run_{}.sh'.format(pairs)
    with open(script, 'w') as out:
        cmd = header + "\n"
        cmd += 'python -m jcvi.compara.catalog ortholog {} {} --cscore=0.99 --no_strip_names --nochpf \n'.format(sample1, sample2)
        cmd += 'convert_anchors_to_link.py {0}.bed {1}.bed {0}.{1}.anchors {0}.{1}.bed\n'.format(sample1, sample2)
        cmd += 'cut -f 1-3,7 {0}.{1}.bed > {0}.synteny.bed\n'.format(sample1, sample2)
        cmd += 'cut -f 4-6,8 {0}.{1}.bed > {1}.synteny.bed'.format(sample1, sample2)
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

        run_synteny(pair)



if __name__ == "__main__":
    main(sys.argv[1:])
