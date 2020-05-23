#!/usr/bin/env python
# -*- coding:utf-8 -*-


import os
import os.path as op
import pandas as pd
import sys


def import_agp(agpfile,
               removeU=True):
    """
    import agp file

    Params:
    -------
    agpfile: `str`, agp file from allhic assembly
    removeU: `bool`:[True], if remove U 

    Returns:
    --------
    out: `pd.DataFrame`

    Examples:
    ---------
    >>> import_agp(`sample.agp`)

    """
    names = ('chrom', 'start', 'end',
             'idx', 'type', 'tig',
             'num', 'length', 'stand')
    df = pd.read_csv(agpfile,
                     sep='\t',
                     header=None,
                     index_col=0,
                     names=names)

    if removeU:
        df = df.loc[df.type != 'U']

    return df


def import_bed(bedfile):
    """
    import bed file

    Params:
    -------
    bedfile: `str`, bed file

    Returns:
    -------
    out: `pd.DataFrame`

    Examples:
    ---------
    >>> import_bed('sample.bed')
    """
    df = pd.read_csv(bedfile,
                     sep='\t',
                     header=None,
                     index_col=0)

    return df


def convert_coordinate(agpfile, bedfile, out):
    agp_df = import_agp(agpfile)
    bed_df = import_bed(bedfile)

    new_df = bed_df.copy()
    new_df.reset_index(inplace=True)
    for i, row in enumerate(bed_df.itertuples()):
        tig = row.Index
        tmp_df = agp_df[agp_df.tig == tig]
        chrom = tmp_df.index.to_list()[0]
        start = tmp_df.start + row[1] - 1
        end = tmp_df.start + row[2] - 1
        new_df.iloc[i, 0] = chrom
        new_df.iloc[i, 1] = start.values[0]
        new_df.iloc[i, 2] = end.values[0]

    new_df.to_csv(out, sep='\t', header=None, index=None)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: {} <agpfile> <bedfile> <outfile>".format(
            op.basename(sys.argv[0])))
        sys.exit()

    convert_coordinate(sys.argv[1], sys.argv[2], sys.argv[3])
