#!/usr/bin/env python
# -*- coding:utf-8 -*-


"""
%prog name.list [options]
or singled-end data 
%prog name.list --suffix='_fp_.fastq.gz' --suffix2=False

    To rename the sra fastq by a list.
"""


import os
import sys


def rename_sra_fastq(namelist, suffix='_R1.fastq.gz', suffix2=None):
    if suffix2 is not False and not suffix2:
        suffix2 = suffix.replace('1','2')
    db = dict([i.strip().split() for i in open(namelist) if i.strip()])
    for i in db:
        try:
            os.rename(i + suffix, db[i] + suffix)
        except OSError:
            print("There is no file of {}".format(i + suffix))
        if suffix2 is not False:
            try:
                os.rename(i + suffix2, db[i] + suffix2)
            except OSError:
                print("There is no file of {}".format(i + suffix2))


if __name__ == "__main__":
    from optparse import OptionParser

    p = OptionParser(__doc__)
    p.add_option("--suffix", default='_R1.fastq.gz',
                help='The suffix of fastq [default: %default]')
    p.add_option("--suffix2", default="None",
            help='The suffix of fastq, if is single set it to False'
            '[default: %default]')

    opts, args = p.parse_args()
    if len(args) != 1:
        sys.exit(p.print_help())

    list_file, = args
    rename_sra_fastq(list_file, opts.suffix, eval(opts.suffix2))
