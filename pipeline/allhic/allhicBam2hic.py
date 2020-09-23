#!/usr/bin/env python
#! -*- coding:utf-8 -*-


from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from TDGP.apps.grid import CMD
from agp2assembly import agp2assembly
from allhicBam2links import bam2links

DIR_3DDNA = "~/software/3d-dna"

def allhicBam2hic(args):
    """
    %(prog)s <sorted.bam> <groups.asm.fasta> <groups.agp> [Options]
        convert allhic bam file to juicer `.hic` file
    """
    p = p=argparse.ArgumentParser(prog=allhicBam2hic.__name__,
                        description=allhicBam2hic.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bam', 
            help='bam file from allhic pipeline')
    pReq.add_argument('fasta', help='groups.asm.fasta')
    pReq.add_argument('agp', help='groups.agp')
    pOpt.add_argument('--software', default=DIR_3DDNA,
            help='3d-dna path [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    bam = args.bam
    fasta = args.fasta
    agp = args.agp 

    out_prefix = bam.split(".")[0]
    out_links = "{}.links".format(out_prefix)
    out_sorted_links = "{}.sorted.links".format(out_prefix)
    out_assembly = "{}.assembly".format(out_prefix)
    # out_asm = "{}.asm".format(out_prefix)
    # out_cprops = "{}.cprops".format(out_prefix)
    bam2links([bam, '-o', out_links])
    agp2assembly([agp, '-o', out_assembly])

    sort_link_cmd = "LC_ALL=C;sort -k2,2n -k6,6n {} > {}".format(
                                        out_links, out_sorted_links)
    os.system(sort_link_cmd)

    # grep_cmd = "grep -v '>' {0} > {1}; \
    #     grep '>' {0} | sed 's/>//g' > {2}".format(out_assembly, out_asm, out_cprops)
    # os.system(grep_cmd)
    
    hic = "{}/visualize/run-assembly-visualizer.sh -p true \
        {} {} {} ".format(args.software,out_assembly, out_sorted_links)
    os.system(hic)
    

if __name__ == "__main__":
    allhicBam2hic(sys.argv[1:])