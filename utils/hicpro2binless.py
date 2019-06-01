#!/usr/bin/env python
# -*- coding:utf-8 -*-


# Modify from HiCPro scripts mapped_hic_fragment.py
"""
%prog test.bam MboI.bed out.tsv
    Script to convert hicpro result to binless.
    test.bam: hicpro bwt final result
    MboI.bed: restriction fragment bed file .
"""

from __future__ import  print_function
from bx.intervals.intersection import Intersecter, Interval
from collections import defaultdict

import os
import os.path as op
import pysam
import sys


def load_restriction_fragment(enzyme_bed):
    """
    Read enzyme bed file and store the intervals in a tree
    :param enzyme_bed: bed file path
    :return: res_frag
    """
    res_frag = defaultdict(lambda :Intersecter())
    with open(enzyme_bed) as fp:
        for nline, line in enumerate(fp, 1):
            line_list = line.strip().split()
            try:
                chrn, start, end, name = line_list[:4]
            except ValueError:
                print("[Warning] wrong input format in line :{}. "
                      "Not a bed file?!".format(nline))
                continue

            start = int(start)
            end = int(end)

            tree = res_frag[chrn]
            tree.add_interval(Interval(start, end, value={'name':name}))

    return res_frag

def get_read_start(read):
    """
    Return the 5' end position of the read
    :param read: AlignmentSeq
    :return: int
    """
    if read.is_reverse:
        pos = read.pos + read.alen - 1
    else:
        pos = read.pos

    return pos

def get_pos_range(res_frag, chrn, read):
    """
    Intersect a given pos with the set of restriction fragments,
    and return a range of fragments.
    """
    pos = get_read_start(read)
    if chrn in res_frag:
        resfrag = res_frag[chrn].find(pos, pos+1)
        if len(resfrag) == 1:
            return resfrag[0]

def get_ordered_reads(read1, read2):
    r1, r2 = read1, read2
    if read1.tid == read2.tid:
        if get_read_start(read1) > get_read_start(read2):
            r1, r2 = read2, read1
    else:
        if read1.tid > read2.tid:
            r1, r2 = read2, read1

    return r1, r2

def get_read_strand(read):
    strand = "1" if not read.is_reverse else "0"
    return strand


def hicpro2binless(bamfile, enzyme_bed, outtsv):
    """
    Convert hicpro result to binless.
    :param bamfile: the bamfile of bowtie_result/bwt/
    :param enzyme_bed: the bedfile of restriction fragments
    """
    if not op.exists(bamfile):
        print("[Error] No such file of {}".format(bamfile))
        return
    if not op.exists(enzyme_bed):
        print("[Error] No such file of {}".format(enzyme_bed))

    res_frag = load_restriction_fragment(enzyme_bed)
    samfile = pysam.Samfile(bamfile, 'rb')
    outtsv = open(outtsv, 'w')

    reads_num = 0
    for read in samfile.fetch(until_eof=True):
        reads_num += 1

        if read.is_read1:
            r1 = read
        elif read.is_read2:
            r2 = read
            if not r1.is_unmapped and not r2.is_unmapped:
                or1, or2 = get_ordered_reads(r1, r2)
                or1_chrom = or1.reference_name
                or2_chrom = or2.reference_name
                or1_pos = get_read_start(or1) + 1
                or2_pos = get_read_start(or2) + 1
                or1_resfrag = get_pos_range(res_frag, or1_chrom, or1)
                or2_resfrag = get_pos_range(res_frag, or2_chrom, or2)

                if or1_resfrag is not None:
                    rup1 = or1_resfrag.start
                    rdn1 = or1_resfrag.end
                if or2_resfrag is not None:
                    rup2 = or2_resfrag.start
                    rdn2 = or2_resfrag.end

                result = (or1.qname, or1_chrom, or1_pos, get_read_strand(or1),
                          or1.alen, rup1, rdn1, or2_chrom, or2_pos, get_read_strand(or2),
                          or2.alen, rup2, rdn2)
                outtsv.write("\t".join(map(str, result)) + "\n")
    outtsv.close()

if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)

    opts, args = p.parse_args()
    if len(args) != 3:
        sys.exit(p.print_help())

    bamfile, enzyme_bed, outtsv = args
    hicpro2binless(bamfile, enzyme_bed, outtsv)
