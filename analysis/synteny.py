#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Genome synteny analysis libraries.
"""

from __future__ import print_function

import argparse 
import logging
import os
import os.path as op
import sys

from collections import OrderedDict
from subprocess import Popen, PIPE
from TDGP.apps.base import debug, check_file_exists, listify
from TDGP.apps.base import ActionDispatcher


def main():

    actions = (
            ("getSyntenyBlock", "To get synteny block from anchor file, "
                "which is generate from jcvi"),
            ('convertAnchorsToLink', 'Convert anchors to link file'), 
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def create_bed_dict(bedfile):
    bed_dict = {}
    check_file_exists(bedfile)
    with open(bedfile) as fp:
        for line in fp:
            line_list = line.strip().split()
            chrom, start, end, gene = line_list[:4]

            bed_dict[gene] = (chrom, start, end)
    
    return bed_dict


def convertAnchorsToLink(args):
    """
    %(prog)s bed1 bed2 anchor [Options]
        To convert anchors file to link bed file, which is generate
            from jcvi
    """
    p = argparse.ArgumentParser(prog=convertAnchorsToLink.__name__,
                        description=convertAnchorsToLink.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bed1', help='gene bed file of species1')
    pReq.add_argument('bed2', help='gene bed file of species2')
    pReq.add_argument('anchor', help='anchor file of synteny gene pairs')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'), 
            default=sys.stdout, help='output file [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    bed1 = args.bed1
    bed2 = args.bed2
    anchor = args.anchor
    check_file_exists(anchor)

    bed1 = create_bed_dict(bed1)
    bed2 = create_bed_dict(bed2)
    popen = Popen(['sort', '-V'], 
                                stdout= PIPE,
                                stdin = PIPE)
    with open(anchor) as fp:
        for line in fp:
            if line.startswith("#"):
                
                continue

            line_list = line.strip().split()
            gene1, gene2 = line_list[:2]
            
            if gene1 not in bed1 or gene2 not in bed2:
                continue
            #chrom1, start1, end1 = bed1[gene]
            #chrom2, start2, end2 = bed2[gene]
            print("{}\n".format('\t'.join(['\t'.join(bed1[gene1]), 
                    '\t'.join(bed2[gene2]), gene1, gene2])), file=args.out)
    


def getSyntenyBlock(args):
    """
    %(prog)s bed1 bed2 anchor [Options]
        To get synteny block from anchor file, which is generate
            from jcvi
    
    """

    p = argparse.ArgumentParser(prog=getSyntenyBlock.__name__,
                        description=getSyntenyBlock.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bed1', help='gene bed file of species1')
    pReq.add_argument('bed2', help='gene bed file of species2')
    pReq.add_argument('anchor', help='anchor file of synteny gene pairs')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'), 
            default=sys.stdout, help='output file [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    bed1 = args.bed1
    bed2 = args.bed2
    anchor = args.anchor
    check_file_exists(anchor)

    bed1 = create_bed_dict(bed1)
    bed2 = create_bed_dict(bed2)
    
    block_num = 0
    tmp_outfile = anchor.replace('anchor', 'bed')

    tmp_out = open(tmp_outfile, 'w')
    with open(anchor) as fp:
        for line in fp:
            if line.startswith("#"):
                block_num += 1
                block = "block{}".format(block_num)
                continue

            line_list = line.strip().split()
            gene1, gene2 = line_list[:2]
            
            if gene1 not in bed1 or gene2 not in bed2:
                continue
            #chrom1, start1, end1 = bed1[gene]
            #chrom2, start2, end2 = bed2[gene]
            tmp_out.write('\t'.join(['\t'.join(bed1[gene1]), 
                    '\t'.join(bed2[gene2]), gene1, gene2, block]) + "\n")
    
    block_db = OrderedDict()
    with open(tmp_outfile) as fp:
        for line in fp:
            line_list = line.strip().split()
            chr1, start1, end1, chr2, start2, end2 = line_list[:6]
            block = line_list[8]
            if block not in block_db:
                block_db[block] = [chr1, start1, end1, chr2, start2, end2]
            
            if block_db[block][0] != chr1:
                continue
            if block_db[block][3] != chr2:
                continue

            if block_db[block][1] > start1:
                block_db[block][1] = start1
            if block_db[block][2] < end1:
                block_db[block][2] = end1
            if block_db[block][4] > start2:
                block_db[block][4] = start2
            if block_db[block][5] < end2:
                block_db[block][5] = end2
            
        for block in block_db:
            print("\t".join([block, "\t".join(block_db[block])]), file=args.out)
        
        logging.debug('Successful ... result is in `{}`'.format(args.out.name))


if __name__ == "__main__":
    main()