#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
convert json to yaml
"""
from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import json
import yaml

def json2yaml(args):
    """
    %(prog)s <json> [Options]
        convert json file to yaml file
        example:
            %(prog)s example.json > example.yaml
    """
    p = argparse.ArgumentParser(prog=json2yaml.__name__,
                        description=json2yaml.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('json', 
            help='json file')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    with open(args.json) as fp:

        data = json.load(fp)
    
    yaml.dump(data, args.output)
    

if __name__ == "__main__":
    json2yaml(sys.argv[1:])