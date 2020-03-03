#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Util for deal vcf file.
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import check_file_exists, debug
from TDGP.apps.base import listify


def main():

    actions = (
            ("test", "test"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def filter_vcf(args):
    """
    filter_vcf
    """
    vcf_df, _type, minDP, remove_snp, remove_indel, out = args
    
    vcf_df.dropna(inplace=True)
    vcf_df = vcf_df[vcf_df.SAMPLE.map(lambda x: x.split(":")[0]).isin(_type)]
    if remove_snp:
        vcf_df = vcf_df[(vcf_df['REF'].map(lambda x: len(x)) > 1) | \
                        (vcf_df['ALT'].map(lambda x: len(x)) > 1)]
    if remove_indel:
        vcf_df = vcf_df[(vcf_df['REF'].map(lambda x: len(x)) == 1) & \
                        (vcf_df['ALT'].map(lambda x: len(x)) == 1)]
    def deal_info(info):
        return dict(map(lambda x: x.split("="), info.split(";")))
    if minDP:
        vcf_df = vcf_df[vcf_df.INFO.map(lambda x: int(deal_info(x)['DP']) >= minDP)]
    
    if out:
        try:
            vcf_df.to_csv(out, sep='\t', header=None)
        except:
            logging.warning('Warning: Failed to output filter vcf result...')
      
    return vcf_df


def import_vcf(vcf_path,
               _type='0/1',
               minDP=0,
               chunksize=1e6,
               threads=12, 
               remove_snp=False, 
               remove_indel=True,
               out=None):
    #allow_type = ('0/1', '1/1', '0/0', 'heter', 'homo', 'variant', 'all')
    standard_type = ('0/1', '1/1', '0/0')
    standard_type_dict = {'heter': '0/1', 'homo': '1/1'}
    if _type not in standard_type:
        if _type == 'all':
            _type = ['0/1', '1/1', '0/0', './.']
        elif _type == 'variant':
            _type = ['0/1', '1/1']
        elif _type == 'heter':
            _type = ['0/1']
        elif __type == 'homo':
            _type = ['1/1']
        else:
            logging.error('Error input the vcf type `{}`'.format(_type))
    else:
        _type = [_type] 
    
    compression = 'gzip' if vcf_path.endswith('.gz') else 'infer'
    vcf_df = pd.read_csv(vcf_path, sep='\t', header=None, 
                         index_col=0,comment="#",
                         chunksize=chunksize,
                         compression=compression,
                         names=('CHROM', 'POS', 
                              'ID', 'REF', 'ALT',
                              'QUAL', 'FILTER', 'INFO',
                              'FORMAT', 'SAMPLE'))
    
    args_list = [(i, _type, minDP,  remove_snp, 
                    remove_indel, out) for i in vcf_df]
    with multiprocessing.Pool(threads) as pool:
        res = pool.map(filter_vcf, args_list)
  
    logging.debug('Loaded vcf file of `{}`'.format(vcf_path))
    return pd.concat(res)




if __name__ == "__main__":
    main()