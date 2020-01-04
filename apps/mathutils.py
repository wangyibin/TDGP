# reference https://bitbucket.org/mirnylab/mirnylib/src/default/mirnylib/numutils.py
# most copy from reference

from __future__ import print_function
import numpy as np

from . import systemutils
def logbinsnew(a, b, ratio=0, N=0):
    a = int(a)
    b = int(b)
    a10, b10 = np.log10([a, b])
    if ratio != 0:
        if N != 0:
            raise ValueError("Please specify N or ratio")
        N = np.log(b / a) / np.log(ratio)
    elif N == 0:
        raise ValueError("Please specify N or ratio")
    data10 = np.logspace(a10, b10, N)
    data10 = np.array(np.rint(data10), dtype=int)
    data10 = np.sort(np.unique(data10))
    assert data10[0] == a
    assert data10[-1] == b
    return data10


def observedOverExpected(matrix):
    "Calculates observedOverExpected of any contact map. Ignores NaNs"

    #cdef np.ndarray[np.double_t, ndim=2] 
    data = np.array(matrix, dtype=np.double, order="C")
    N = data.shape[0]
    
    _bins = logbinsnew(1, N, 1.03)
    _bins = [(0, 1)] + [(_bins[i], _bins[i+1]) for i in range(len(_bins) - 1)]
    #cdef np.ndarray[np.int64_t, ndim=2] 
    bins = np.array(_bins, dtype=np.int64, order="C")
    M = bins.shape[0]
    
    for k in range(M):
        start, end = bins[k, 0], bins[k, 1]
        ss = 0
        count = 0
        for offset in range(start, end):
            for j in range(0, N - offset):
                x = data[offset + j, j]
                if np.isnan(x):
                    continue
                ss += x
                count += 1
        
        meanss = ss / count
        if meanss != 0:
            for offset in range(start,end):
                for j in range(0,N-offset):
                    data[offset + j, j] /= meanss
                    if offset > 0: data[j, offset+j] /= meanss
    return data

def fastMatrixSTD(inMatrix):
    "Estimates variance of the matrix"
    inMatrix = np.asarray(inMatrix)

    if len(inMatrix.flat) < 5:
        return np.mean(inMatrix)

    if len(inMatrix.flat) < 10000:
        return  inMatrix.var()

    else:
        varvar = inMatrix.flat[::100].var() + inMatrix.flat[33::231].var()\
            + inMatrix.flat[24::242].var() + inMatrix[55::413].var()
    return np.sqrt(0.25 * varvar)

def isSymmetric(inMatrix):
    """
    Checks if the supplied matrix is symmetric.
    """

    M = inMatrix
    varDif = np.abs(M - M.T).max()
    return varDif < fastMatrixSTD(inMatrix) * 0.0000001 + 0.00001

def fillDiagonal(inArray, diag, offset=0):
    "Puts diag in the offset's diagonal of inArray"
    N = inArray.shape[0]
    assert inArray.shape[1] == N
    if offset >= 0:
        inArray.flat[offset:N * (N - offset):N + 1] = diag
    else:
        inArray.flat[(-offset * N)::N + 1] = diag


def removeDiagonals(inArray, m):
    """removes up to mth diagonal in array
    m = 0: main
    m = 1: 3 diagonals, etc.
    """
    for i in range(-m, m + 1):
        fillDiagonal(inArray, 0, i)

def ultracorrectSymmetricWithVector(x,v = None,M=None,diag = -1,
                                    tolerance=1e-5):
    """Main method for correcting DS and SS read data. Possibly excludes diagonal.
    By default does iterative correction, but can perform an M-time correction"""
    
    if M == None:
        M = 599
    totalBias = np.ones(len(x),float)
    if v == None: v = np.zeros(len(x),float)  #single-sided reads
    x = np.array(x,np.double,order = 'C')
    _x = x
    v = np.array(v,float,order = "C")
    N = len(x)

    for iternum in range(M):
        s0 = np.sum(_x,axis = 1)

        mask = s0 == 0

        v[mask] = 0   #no SS reads if there are no DS reads here
        nv = v / (totalBias * (totalBias[mask==False]).mean())
        s = s0 + nv
        for dd in range(diag + 1):   #excluding the diagonal
            if dd == 0:
                s -= np.diagonal(_x)
            else:
                dia = np.array(np.diagonal(_x,dd))
                s[dd:] = s[dd:] -  dia
                s[:len(s)-dd] = s[:len(s)-dd] - dia
        s = s / np.mean(s[s0!=0])
        s[s0==0] = 1
        s -= 1
        s *= 0.8
        s += 1
        totalBias *= s

        for i in range(N):
            for j in range(N):
                _x[i,j] = _x[i,j] / (s[i] * s[j])

        if M == 599:
            if np.abs(s-1).max() < tolerance:
                #print "IC used {0} iterations".format(iternum+1)
                break


    corr = totalBias[s0!=0].mean()  #mean correction factor
    x  = x * corr * corr #renormalizing everything
    totalBias /= corr
    return x, v/totalBias, totalBias

def ultracorrectAssymetric(x, M="auto", tolerance="unused"):  # @UnusedVariable
    "just iterative correction of an assymetric matrix"
    if M == "auto":
        M = 50
    x = np.array(x, float)
    if x.sum() == 0:
        return x
    print(np.mean(x))
    newx = np.array(x)
    for _ in range(M):
        correctInPlace(newx)

    print(np.mean(newx), end=' ')
    newx /= (1. * np.mean(newx) / np.mean(x))
    print(np.mean(newx))
    return newx

def iterativeCorrection(x, M="auto", tolerance=1e-6,
                        symmetric="auto",
                        skipDiags=-1):
    "A wrapper for iterative correction of any matrix"
    x = np.asarray(x, dtype=float)

    if len(x.shape) != 2:
        raise ValueError("Only 2D matrices are allowed!")
    if x.shape[0] != x.shape[1]:
        symmetric = False
    elif symmetric.lower() == "auto":
        symmetric = isSymmetric(x)

    if not symmetric:
        if M == "auto":
            M = 50
        print("Matrix is not symmetric, doing {0} iterations of IC".format(M))
        return ultracorrectAssymetric(x, M), 1
    if M == "auto":
        M = None  # default of ultracorrectSymmetricWithVector
    else:
        try:
            M = int(M)
        except:
            raise ValueError("Please provide integer for M; {0} provided".format(M))
    corrected, dummy, bias = ultracorrectSymmetricWithVector(x, M=M, tolerance=tolerance,
                                                             diag=skipDiags)
    return corrected, bias

def ultracorrect(*args, **kwargs):
    return iterativeCorrection(*args, **kwargs)[0]
ultracorrect = systemutils.deprecate(ultracorrect,
     message="Please use iterativeCorrection instead of ultracorrect")


def completeIC(hm, minimumSum=40, diagsToRemove=2, returnBias=False, minimumNumber=10, minimumPercent=.1, returnBiasMask=False):
    """Makes a safe iterative correction
    (i.e. with removing low-coverage regions and diagonals)
    for a symmetric heatmap
    Only keeps rows/columns with sum more than minimumSum,
    and with at least  minumumNumber or minimumPercent (default 20 & 1/5 length) non-zero entries
    """
    assert isSymmetric(hm)
    hm = np.asarray(hm, dtype=float)

    hmc = hm.copy()  # to remove diagonals safely
    removeDiagonals(hmc, diagsToRemove - 1)
    matsum = np.sum(hmc, axis=0)
    mask = matsum > minimumSum
    num = min(len(hmc) * minimumPercent, minimumNumber)
    mask = mask * (np.sum(hmc > 0, axis=0) > num)

    if mask.sum() + 3 < (matsum > 0).sum() * 0.5:
        warnings.warn("""Iterative correction will remove more than a half of the matrix
        Check that values in rows/columns represent actual reads,
        and the sum over rows/columns exceeds minimumSum""")

    hmc[~mask] = 0
    hmc[:, ~mask] = 0
    if hmc.sum() == 0:
        if returnBias:
            return np.zeros_like(hm), np.zeros(len(hm))
        if returnBiasMask:
            return np.zeros_like(hm), np.zeros(len(hm)), np.zeros(len(hm))
        return np.zeros_like(hm)

    hm, bias = iterativeCorrection(hmc, skipDiags=diagsToRemove)
    dmean = np.median(np.diagonal(hm, diagsToRemove))
    for t in range(-diagsToRemove + 1, diagsToRemove):
        fillDiagonal(hm, dmean, t)
    hm[~mask] = 0
    hm[:, ~mask] = 0

    if returnBias:
        return hm, bias
    if returnBiasMask:
        return hm, bias, mask
    return hm