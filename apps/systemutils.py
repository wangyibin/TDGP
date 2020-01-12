#(c) 2012 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maksim Imakaev (imakaev@mit.edu)

"""
Some important utilities from Max. This includes:

Set exception hook to pdb
Run in separate process
fork-map
fork-map-reduce
fork-map-average
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import os
import sys
import pickle
import pdb
import traceback
import warnings
import subprocess
import numpy as np
from copy import copy
import logging
from functools import reduce
log = logging.getLogger(__name__)


def commandExists(command):
    "checks if the bash command exists"
    command = command.split()[0]
    if subprocess.call(['which', command]) != 0:
        return False
    return True

def gzipWriter(filename, pigzArguments=("-4",)):
    """
    creates a writing process with gzip or parallel gzip (pigz) attached to it
    """
    filename = os.path.abspath(filename)
    with open(filename, 'wb') as outFile:
        if commandExists("pigz"):
            writer = ["pigz", "-c"] + list(pigzArguments)
        else:
            writer = ["gzip", "-c", "-1"]
            warnings.warn("Please install 'pigz' parallel gzip for faster speed")

        pwrite = subprocess.Popen(writer, stdin=subprocess.PIPE, stdout=outFile, shell=False, bufsize=-1)
    log.info("""Writer created with command "{0}" """.format(writer))
    return pwrite

def _exceptionHook(infoType, value, tb):
    "Exception hook"
    traceback.print_exception(infoType, value, tb)
    print()
    pdb.post_mortem(tb)


def setExceptionHook():
    "sets exception hook to pdb"
    sys.excepthook = _exceptionHook


class transparentDict(dict):  # transparent dictionary, that returns the key
    def __missing__(self, key):
        return key


def run_in_separate_process(func, *args, **kwds):
    pread, pwrite = os.pipe()
    pid = os.fork()
    if pid > 0:
        os.close(pwrite)
        with os.fdopen(pread, 'rb') as f:
            status, result = pickle.load(f)
        os.waitpid(pid, 0)
        if status == 0:
            return result
        else:
            raise result
    else:
        os.close(pread)
        try:
            result = func(*args, **kwds)
            status = 0
        except Exception as exc:
            result = exc
            status = 1
        with os.fdopen(pwrite, 'wb') as f:
            try:
                pickle.dump((status, result), f, pickle.HIGHEST_PROTOCOL)
            except pickle.PicklingError as exc:
                pickle.dump((2, exc), f, pickle.HIGHEST_PROTOCOL)
        os._exit(0)


def deprecate(newFunction, oldFunctionName=None, message=None):
    """If you rename your function, you can use this to issue deprecation warning for the old name
    Juse use   newFunction = deprecate(oldFunction)"""
    try:
        newName = newFunction.__name__
    except:
        newName = "_UndeterminedName_"
    if oldFunctionName is None:
        oldFunctionName = "_UnspecifiedName_"
    if message == None:
        message = "Function %s was renamed to %s" % (
            oldFunctionName, newName)

    def oldFunction(*args, **kwargs):
        warnings.warn(message)
        return newFunction(*args, **kwargs)
    return oldFunction

def _nprocessors():
    if sys.platform == 'darwin':
        try:
            from multiprocessing import cpu_count
            return cpu_count()
        except NotImplementedError:
            pass
    else:
        # Cygwin (Windows) and Linuxes
        # Could try sysconf(_SC_NPROCESSORS_ONLN) (LSB) next.  Instead, count processors in cpuinfo.
        try:
            s = open('/proc/cpuinfo', 'r').read()
            return s.replace(' ', '').replace('\t', '').count('processor:')
        except:
            pass
    return 4

nproc = _nprocessors()

def fmap(f, *a, **kw):
    """
    forkmap.map(..., n=nprocessors), same as map(...).
    n must be a keyword arg; default n is number of physical processors.
    """
    n = max([kw.get(i, 0) for i in ['n','N', "nproc", "Nproc", "NProc"]])
    if n == 0:
        n = nproc

    if n == 1:
        return list(map(f, *a))

    L = list(zip(*a))
    n = min(n, len(L))

    ans = [None] * len(L)
    pipes = [os.pipe() for i in range(n - 1)]

    for i in range(n):
        if i < n - 1 and not os.fork():  # Child, and not last processor
            try:
                try:
                    obj = [f(*x) for x in L[i::n]]
                except Exception as obj:
                    pass
                with os.fdopen(pipes[i][1],'wb') as f:
                    pickle.dump(obj,f, protocol=pickle.HIGHEST_PROTOCOL)
            except:
                traceback.print_exc()
            finally:
                os._exit(0)
        elif i == n - 1:  # parent
            try:
                ans[i::n] = [f(*x) for x in L[i::n]]
                for k in range(n - 1):
                    with os.fdopen(pipes[k][0],'rb') as f:
                        obj = pickle.load(f)
                    if isinstance(obj, Exception):
                        raise obj
                    ans[k::n] = obj
            finally:
                for j in range(n - 1):
                    os.wait()
    return ans




def _testFmap():

    for i in range(1, 300):
        print(i)
        a = list(range(i))
        for j in range(1, 10):
            b = fmap(lambda x:x, a, n=j)
            assert (np.array(a) == np.array(b)).all()

def _fmapredcount(function, data, reduction=lambda x, y: x + y, n=4, exceptionList=[IOError]):
    """fork-map-reduce
    Performs fork-map of function on data, automatically reducing the data inside each worker.
    If evaluation throws the exception from exceptionList, this results are simply ignored
    """
    def funsum(x, y):
        """reduces two x[0],y[0], keeping track of # of
        successful evaluations that were made
        Also keeps track of None's that can occur if evaluation failed"""
        if x is None:
            if y is None:
                return None
            else:
                return y
        else:
            if y is None:
                return x
            else:
                return (reduction(x[0], y[0]), x[1] + y[1])

    def newfunction(x):
        try:
            "if function is evaluated, it was evaluated one time"
            return function(x), 1
        except tuple(exceptionList):
            return None

    if len(data) < n:
        n = len(data)
    datas = []

    for i in range(n):
        datas.append(copy(data[i::n]))  # split like that if beginning and end of the array have different evaluation time

    def worker(dataList):
        dataList[0] = newfunction(dataList[0])
        return reduce(lambda z, y: funsum(z, newfunction(y)), dataList)  # reducing newfunction with our new reduction algorithm

    reduced = fmap(worker, datas, n=n)
    return reduce(funsum, reduced)


def fmapred(function, data, reduction=lambda x, y: x + y, n=4, exceptionList=[IOError]):
    """reduces two x[0],y[0], keeping track of # of
    successful evaluations that were made
    Also ignores failed evaluations with exceptions from exceptionList.

    Parameters
    ----------
    function : function
        function to be applied to the data
    data : iterable
        input data
    reduction : function, optional
        Reduction function. By default - sum
    n : int, optional
        number of CPUs
    exceptionList : list, optional
        list of exceptions to be ignored during reduction. By default, only IOError is ignored.
    """
    return _fmapredcount(function, data, reduction=reduction, n=n, exceptionList=exceptionList)[0]


def fmapav(function, data, reduction=lambda x, y: x + y, n=4, exceptionList=[IOError]):
    """Calculates averate of [fucntion(i) for i in data]
    Also ignores failed evaluations with exceptions from exceptionList.

    Parameters
    ----------
    function : function
        function to be applied to the data
    data : iterable
        input data
    reduction : function, optional
        Reduction function. By default - sum
    n : int, optional
        number of CPUs
    exceptionList : list, optional
        list of exceptions to be ignored during reduction. By default, only IOError is ignored.
    """

    a = _fmapredcount(function, data, reduction=reduction, n=n,
                      exceptionList=exceptionList)
    return a[0] / float(a[1])
