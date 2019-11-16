#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import os
import sys


def check_file(infile):
    if not os.path.exists(infile):
        sys.exit()

    
    with open(infile) as fp:
        for line in fp:
            if len(line.strip().split()) != 12:
                continue
            print(line, end="", file=sys.stdout)
    

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: {} sample.allVaildParis".format(sys.argv[0]))
        sys.exit()

    check_file(sys.argv[1])
