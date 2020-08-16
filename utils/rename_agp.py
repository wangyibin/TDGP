#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import argparse
import os
import os.path as op
import sys


def rename_agp(agp, rename_list):
    rename_db = dict(i.strip().split() 
                    for i in open(rename_list)
                    if i.strip())
    with open(agp) as fp:
        for line in fp:
            line_list = line.strip().split()
            group = line_list[0]
            if group in rename_db:
                line_list[0] = rename_db[group]
            print("\t".join(line_list), file=sys.stdout)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: {} <groups.agp> <rename.list> > out.agp ".format(sys.argv[0]))
        sys.exit()
    rename_agp(sys.argv[1], sys.argv[2])
