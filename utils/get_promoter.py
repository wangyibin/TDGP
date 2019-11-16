#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Time: 2019/7/13 16:33

"""
%prog sample.gff3 chrom.sizes
    To get gene promoter from a gfffile, the promoter is defined as TSS upstream 2kb.
"""
import gzip
import sys

def create_gene_db(gfffile, chrom_sizes):
    chrom_len_db = dict(i.strip().split() for i in open(chrom_sizes) if i.strip())
    if gfffile[-2:] == "gz":
        f = gzip.open(gfffile)
    else:
        f = open(gfffile)
    with f as fp:
        for line in fp:
            if line[0] != "#":
                line_list = line.strip().split()
                if line_list[2] == 'gene':
                    chrom = line_list[0]
                    chrom_len = chrom_len_db[chrom]
                    start = int(line_list[3])
                    end = int(line_list[4])
                    strand = line_list[6]
                    name = line_list[8].split(";")[0].split('=')[1]
                    if strand == "+":
                        promoter_start = start - 2000 - 1
                        promoter_end = start - 1
                        if promoter_start <= 0:
                            promoter_start = 0
                    else:
                        promoter_end = end + 2000 +1
                        promoter_start = end + 1
                        if promoter_end >= chrom_len:
                            promoter_end = chrom_len
                    print("\t".join(map(str,[chrom, promoter_start, promoter_end, name, ".", strand])))

if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) != 2:
        sys.exit(p.print_help())
    gfffile, chrom_sizes = args
    create_gene_db(gfffile, chrom_sizes)
