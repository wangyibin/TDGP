#!/usr/bin/env python
# -*- coding:utf-8 -*-


"""
%prog sample.cool chrom.list [Options]

"""

from __future__ import print_function


import argparse
import os
import os.path as op
import sys

from TDGP.apps.grid import CMD


def main(args):
    p = argparse.ArgumentParser(prog=main.__name__,
                        description=main.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('coolfile', help='coolfile')
    pReq.add_argument('chromlist', help='chrom list')
    pOpt.add_argument('-b', '--bigwig', nargs="*", default=None,
            help='bigwig file ')
    pOpt.add_argument('--bgYlabel', nargs="*", default=None,
            help='ylabel of bigwig track')
    pOpt.add_argument('-o', '--outprefix', default=None,
            help='outprefix of heatmap [default: cool_outprefix')
    pOpt.add_argument('-c', '--cmap', default='RdYlBu_r',
            help='colormap of heatmap [default: %(default)s]')
    pOpt.add_argument('-t', '--threads', default=4, type=int,
            help='threads of programs')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    coolfile = args.coolfile 
    chrom_list = args.chromlist
    outprefix = args.outprefix
    bgfile = args.bigwig
    bgylabel = args.bgYlabel
    threads = args.threads
    cmap = args.cmap
    
    

    if op.exists(chrom_list):
        chrom_list = [ i.strip().split()[0] for i in open(chrom_list)
                if i.strip()]
    else:
        chrom_list = chrom_list.strip(",").split(",")
    chrom_order = " ".join(chrom_list)
    if not outprefix:
        outprefix = coolfile.split(".")[0]
    suffix = ""
    ylabel_suffix = ""
    if bgfile:
        suffix = "--bigwig {}".format(" ".join(bgfile))
        if bgylabel:
            ylabel_suffix = " --pyLabel {}".format(" ".join(bgylabel))
    
    suffix = suffix + ylabel_suffix
    cmd_wg = "hicPlotMatrix --matrix {} -o {}_wg_heatmap.pdf --log1p --dpi 300 --title 'Hi-C Heatmap for Whole Genome' --chromosomeOrder {}  --clearMaskedBins --colorMap {} {}".format(coolfile, outprefix, chrom_order, cmap, suffix)
    cmd_per = "hicPlotMatrix --matrix {} -o {}_per_heatmap.pdf --log1p --dpi 300 --perChromosome --chromosomeOrder {}  --clearMaskedBins --colorMap {} {}".format(coolfile, outprefix, chrom_order,cmap, suffix)
    cmd_wg_png = "hicPlotMatrix --matrix {} -o {}_wg_heatmap.png --log1p --dpi 300 --title 'Hi-C Heatmap for Whole Genome' --chromosomeOrder {}  --clearMaskedBins --colorMap {} {}".format(coolfile, outprefix, chrom_order, cmap, suffix)
    cmd_per_png = "hicPlotMatrix --matrix {} -o {}_per_heatmap.png --log1p --dpi 300 --perChromosome --chromosomeOrder {}  --clearMaskedBins --colorMap {} {}".format(coolfile, outprefix, chrom_order,cmap, suffix)
    cmds = []
    cmds.extend([cmd_wg, cmd_per, cmd_wg_png, cmd_per_png])
    
    for chrom in chrom_list:
        cmd_single = "hicPlotMatrix --matrix {0} -o {1}_{2}_heatmap.pdf --log1p --dpi 300 --title 'Hi-C Heatmap' --region {2} --clearMaskedBins --colorMap {3} {4}".format(coolfile, outprefix, chrom, cmap, suffix)
        cmds.append(cmd_single)
        cmd_single_png = "hicPlotMatrix --matrix {0} -o {1}_{2}_heatmap.png --log1p --dpi 300 --title 'Hi-C Heatmap' --region {2} --clearMaskedBins --colorMap {3} {4}".format(coolfile, outprefix, chrom, cmap, suffix)
        cmds.append(cmd_single_png)
    
    CMD(cmds, threads)
if __name__ == "__main__":
    main(sys.argv[1:])