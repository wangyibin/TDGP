#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import argparse
import os
import os.path as op
import sys

from Bio import SeqIO

def rename_fasta(infasta,rename_list, outfasta):
    rename_db = dict(i.strip().split() for i in open(rename_list) 
                    if i.strip())
    with open(infasta) as fp, open(outfasta, 'w') as out:
        fasta = SeqIO.parse(fp, 'fasta')
        for record in fasta:
            if record.name in rename_db:
                record.name = rename_db[record.name]
                record.id = record.name
            print(">{}\n{}".format(record.name, str(record.seq)), file=out)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: {} <in.fasta> <rename.list> <out.fasta>".format(sys.argv[0]))
        sys.exit()
    infasta, rename_list, outfasta = sys.argv[1:]
    rename_fasta(infasta, rename_list, outfasta)

