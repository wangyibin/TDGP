#!/usr/bin/env python
#! -*- coding:utf-8 -*-

# @Time: 2019/6/10 8:37

"""
%prog rice_100000_abs.bed rice [options]
    generate a header file for cworld from bed file.
"""

from __future__ import print_function

import sys

from collections import defaultdict


def bed2header(bedfile, genome, perchr=False):
    out_dict = defaultdict(lambda :[])
    with open(bedfile) as fp:
        for line in fp:
            line_list = line.strip().split()
            third = "{}:{}-{}".format(line_list[0], line_list[1], line_list[2])
            out = [ line_list[3], genome, third]
            if perchr:
                out_dict[line_list[0]].append(out)
            else:
                print("|".join(out), file=sys.stdout)

        if perchr:
            for chrn in out_dict:
                with open('{}_{}.header'.format(chrn, bedfile), 'w') as out:
                    out.write("\n".join(map(lambda x: "|".join(x),
                                                 out_dict[chrn])))


if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)
    p.add_option("-c", "--perchr", dest="perchr", action='store_true',
                 default=False, help='If specified the output are '
                                     'written per chromosome')
    opts, args = p.parse_args()
    if len(args) != 2:
        sys.exit(p.print_help())

    bedfile, genome = args
    bed2header(bedfile, genome, opts.perchr)
