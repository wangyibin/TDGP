#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Util tools for gff formats.
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


def import_gff(gff_path, _type=None):
    compression = 'gzip' if gff_path.endswith('.gz') else 'infer'
    gff_df = pd.read_csv(gff_path, sep='\t', 
                         header=None, 
                         index_col=0, 
                         comment="#", 
                         compression = compression,
                         names=('chrom', 'source',
                                'type', 'start', 
                                'end', 'score', 
                                'strand','phase',
                                'attributes'))
    if _type in set(gff_df.head(100)['type']):
        if _type: 
            gff_df = gff_df[gff_df['type'] == _type]
    else:
        logging.warning('Warning: Failed to filter data by' 
              'type, input type is not correct')
    
    logging.debug('Loaded gff file of `{}`'.format(gff_path))
    return gff_df

def get_gene_id(attributes):
    db = dict(map(lambda x: x.split('='), 
                  [i for i in attributes.split(';') if i]))
    return db['ID']




if __name__ == "__main__":
    main()