#!/usr/bin/env python
# -*- coding:utf-8 -*-


"""
%prog sample.cool chrom.list [Options]

"""

from __future__ import print_function

import os
import os.path as op
import sys



def main(coolfile, chrom_list, outprefix=None, bgfile=None, cmap='RdYlBu_r'):
    if op.exists(chrom_list):
        chrom_list = [ i.strip().split()[0] for i in open(chrom_list)
                if i.strip()]
    else:
        chrom_list = chrom_list.strip(",").split(",")
    chrom_order = " ".join(chrom_list)
    if not outprefix:
        outprefix = coolfile.split(".")[0]
    suffix = ""
    if bgfile:
        suffix = "--bigwig {}".format(bgfile)
    cmd_wg = "conda run -n hicexplorer hicPlotMatrix --matrix {} -o {}_wg_heatmap.pdf --log1p --dpi 300 --title 'Hi-C Heatmap for Whole Genome' --chromosomeOrder {}  --clearMaskedBins --colorMap {} {}".format(coolfile, outprefix, chrom_order, cmap, suffix)
    cmd_per =  "conda run -n hicexplorer hicPlotMatrix --matrix {} -o {}_per_heatmap.pdf --log1p --dpi 300 --title 'Hi-C Heatmap'  --perChromosome --chromosomeOrder {}  --clearMaskedBins --colorMap {} {}".format(coolfile, outprefix, chrom_order,cmap, suffix)
    print(cmd_wg)
    print(cmd_per)
    os.system(cmd_wg)
    os.system(cmd_per)


if __name__ == "__main__":
    from optparse import OptionParser

    p = OptionParser(__doc__)
    
    p.add_option("-o", "--outprefix", default=None, 
            help="the outprefix of file [default: coolfile_prefix]")
    p.add_option("--bigwig", default=None,
            help="the bigwig file of ab compartments")
    p.add_option("--cmap", default='RdYlBu_r', 
            help="colormap of heatmap [default: %default]")
    opts, args = p.parse_args()
    if len(args) != 2:
        sys.exit(p.print_help())
    coolfile, chrom_list = args
    main(coolfile, chrom_list, opts.outprefix, opts.bigwig, opts.cmap)
