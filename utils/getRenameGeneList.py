#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
get rename list from gff file
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from collections import defaultdict, OrderedDict, Counter

def getGeneDb(gff):
    """
    get gene and chromosome database from gff file
    """
    db = defaultdict(list)
    with open(gff) as fp:
        for line in fp:
            if line.strip() == "" or line[0] == "#":
                continue
            line_list = line.strip().split()
            chrom = line_list[0]
            type_ = line_list[2]
            attribute = line_list[8]
            
            attribute_db = dict(
                map(lambda x: x.split("="), attribute.strip(";").split(";")))
            ID = attribute_db["ID"]
            if line_list[2] == 'gene':
                db[chrom].append(ID)

    return db

def getRenameGene(db, prefix):
    """
    generate a renamed list according with chromosome number
    """
    renamed_list = []
    tig_count = 1
    for chrom in db:
        if chrom[:3] == "Chr":
            chrom_num = int(chrom.replace("Chr", ""))
            gene_list = sorted(db[chrom], key=lambda x: int(x.replace(prefix, "")))
            count = 1
            for gene in gene_list:
                renamed_list.append((gene, "Sspon_mono.{:02d}G{:05d}0".format(chrom_num, count)))
                count += 1
        else:
            gene_list = db[chrom]
            for gene in gene_list:
                renamed_list.append((gene, "Sspon_mono.ctg{:05d}0".format(tig_count)))
                tig_count += 1
    return renamed_list


def getRenameGeneList(args):
    """
    %(prog)s.py <in.gff> [Options]
        get renamed gene list from a gff file.
        Examples:
        --------------------------------------
        %(prog)s.py in.gff > renamed.list
    """
    p = p=argparse.ArgumentParser(prog=getRenameGeneList.__name__,
                        description=getRenameGeneList.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('gff', 
            help='gff file')
    pReq.add_argument('prefix', help='Gene prefix before numbers')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    db = getGeneDb(args.gff)
    renamed_list = getRenameGene(db, args.prefix)
    for gene_pair in renamed_list:
        print("\t".join(gene_pair), file=args.output)


if __name__ == "__main__":
    getRenameGeneList(sys.argv[1:])
