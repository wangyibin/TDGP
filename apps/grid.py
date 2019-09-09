#!/usr/bin env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import os
import os.path as op
import logging
import sys

from multiprocessing import Pool, Process, cpu_count
from TDGP.apps.base import debug


debug()


class Jobs(object):
    def __init__(self):
        pass

class Parallel(object):
    """
    Run commands in parallel.
    """
    def __init__(self, target, args, threads=cpu_count()):
        self.target = target
        self.args = args
        self.threads = min(len(args), threads)

    def run(self):
        p = Pool(self.threads)
        res = p.map(self.target, self.args)
        return res


def parallel(target, args, threads=cpu_count):
    p = Pool(min(len(args), threads))
    res = p.map(target, args)
    return res
