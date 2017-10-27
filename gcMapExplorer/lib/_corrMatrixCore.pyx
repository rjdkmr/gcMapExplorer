#!/usr/bin/env python3
#
# Author: Rajendra Kumar
#
# This file is part of gcMapExplorer
# Copyright (C) 2016-2017  Rajendra Kumar, Ludvig Lizana, Per Stenberg
#
# gcMapExplorer is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gcMapExplorer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gcMapExplorer.  If not, see <http://www.gnu.org/licenses/>.
#
#=============================================================================


import numpy as np
cimport numpy as np
cimport cython
import scipy as sp
import sys
from libc.math cimport isnan, sqrt, pow

dtypeFloat32 = np.float32
dtypeFloat = np.float
dtypeInt = np.int
dtypeBool = np.bool
ctypedef np.int_t dtypeInt_t
ctypedef np.float_t dtypeFloat_t
ctypedef np.float32_t dtypeFloat32_t

# declare the interface to the C code
cdef extern from "corrMatrixCoreSRC.h" nogil:
    int _correlationMatrix(double *in_array, double *out_array, int size, double *maskvalue)
    int _covarianceMatrix(double *in_array, double *out_array, int size, double *maskvalue)
    double _calculateCorrelation(double *x, double*y, int n, double *maskvalue)
    double _calculateCovariance(double *x, double*y, int n, double *maskvalue)


@cython.boundscheck(False)
@cython.wraparound(False)
def calculateCorrMatrix(np.ndarray[double, ndim=2, mode="c"] in_array not None, np.ndarray[double, ndim=2, mode="c"] out=None, maskvalue=None):
    """Calculate correlation matrix from a 2D numpy array.
    During calculation, array values equal to maskvalue will not be considered.

    .. note: In the background it is implemented in C language and therefore,
        it is accurate and much faster than the ``numpy.ma.corrcoef``.

    .. warning: input and output arrays should be of data type ``numpy.float64``.
        Therefore, convert type to ``numpy.float64`` before using this function.

    Parameters
    ----------
    in_array : numpy.ndarray
        Input numpy array
    out : numpy.ndarray
        If it is ``None``, output array is returned.
    maskvalue : float
        If this value is given, all elements of input array with this equal
        value is masked during the calculation.

    Returns
    -------
    out : numpy.ndarray or None
        In case of ``out=None``, ouput array will be returned. Otherwise, ``None``
        is returned
    """
    cdef int size
    cdef double c_maskvalue

    size = in_array.shape[0]
    toReturn = False
    if out is None:
        out = np.zeros((size, size), dtype=in_array.dtype)
        toReturn = True

    if maskvalue is not None:
        c_maskvalue = dtypeFloat32(maskvalue)
        with nogil:
            _correlationMatrix(&in_array[0,0], &out[0,0], size, &c_maskvalue)
    else:
        with nogil:
            _correlationMatrix(&in_array[0,0], &out[0,0], size, NULL)


    if toReturn:
        return out
    else:
        return None

@cython.boundscheck(False)
@cython.wraparound(False)
def calculateCovMatrix(np.ndarray[double, ndim=2, mode="c"] in_array not None, np.ndarray[double, ndim=2, mode="c"] out=None, maskvalue=None):
    """Calculate covariance matrix from a 2D numpy array.
    During calculation, array values equal to maskvalue will not be considered.

    .. warning: input and output arrays should be of data type ``numpy.float64``.
        Therefore, convert type to ``numpy.float64`` before using this function.

    Parameters
    ----------
    in_array : numpy.ndarray
        Input numpy array
    out : numpy.ndarray
        If it is ``None``, output array is returned.
    maskvalue : float
        If this value is given, all elements of input array with this equal
        value is masked during the calculation.

    Returns
    -------
    out : numpy.ndarray or None
        In case of ``out=None``, ouput array will be returned. Otherwise, ``None``
        is returned

    """

    cdef int size
    cdef double c_maskvalue

    size = in_array.shape[0]
    toReturn = False
    if out is None:
        out = np.zeros((size, size), dtype=in_array.dtype)
        toReturn = True

    if maskvalue is not None:
        c_maskvalue = maskvalue
        _covarianceMatrix(&in_array[0,0], &out[0,0], size, &c_maskvalue)
    else:
        _covarianceMatrix(&in_array[0,0], &out[0,0], size, NULL)


    if toReturn:
        return out
    else:
        return None

@cython.boundscheck(False)
@cython.wraparound(False)
def calculateCorrelation(np.ndarray[double, ndim=1, mode="c"] x not None, np.ndarray[double, ndim=1, mode="c"] y not None, maskvalue=None):
    """Calculate correlation between two 1D numpy array.
    During calculation, array values equal to maskvalue will not be considered.

    .. warning: input array should be of data type ``numpy.float64``.
        Therefore, convert type to ``numpy.float64`` before using this function.

    Parameters
    ----------
    x : numpy.ndarray
        First input numpy array
    y : numpy.ndarray
        Second input numpy array
    maskvalue : float
        If this value is given, all elements of input array with this equal
        value is masked during the calculation.

    Returns
    -------
    result : float
        Correlation coefficient

    """
    cdef int size
    cdef double c_maskvalue
    cdef double result

    size = x.shape[0]

    if maskvalue is not None:
        c_maskvalue = maskvalue
        with nogil:
            result = _calculateCorrelation(&x[0], &y[0], size, &c_maskvalue)
    else:
        with nogil:
            result = _calculateCorrelation(&x[0], &y[0], size, NULL)

    return result

@cython.boundscheck(False)
@cython.wraparound(False)
def calculateCovariance(np.ndarray[double, ndim=1, mode="c"] x not None, np.ndarray[double, ndim=1, mode="c"] y not None, maskvalue=None):
    """Calculate covariance between two 1D numpy array.
    During calculation, array values equal to maskvalue will not be considered.

    .. warning: input array should be of data type ``numpy.float64``.
        Therefore, convert type to ``numpy.float64`` before using this function.

    Parameters
    ----------
    x : numpy.ndarray
        First input numpy array
    y : numpy.ndarray
        Second input numpy array
    maskvalue : float
        If this value is given, all elements of input array with this equal
        value is masked during the calculation.

    Returns
    -------
    result : float
        Covariance value


    """
    cdef int size
    cdef double c_maskvalue
    cdef double result

    size = x.shape[0]

    if maskvalue is not None:
        c_maskvalue = maskvalue
        with nogil:
            result = _calculateCovariance(&x[0], &y[0], size, &c_maskvalue)
    else:
        with nogil:
            result = _calculateCovariance(&x[0], &y[0], size, NULL)

    return result

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def _calculateCorrMatrixOLDSLOW(np.ndarray[dtypeFloat_t, ndim=2] matrix):
    assert matrix.dtype == dtypeFloat
    cdef:
        int length = matrix.shape[0]
        int boolCount = 0
        int start = 0
        int step = 1
        double r = 0
        int i = 0
        int j = 0
        int k = 0
        np.ndarray[dtypeFloat_t, ndim=2] rvec = np.zeros((2,2), dtype=np.float)
        np.ndarray[dtypeFloat_t, ndim=2] corrMatrix = np.zeros((length,length), dtype=np.float)
        np.ndarray[np.uint8_t, cast=True, ndim=2] boolMatrix = np.zeros((length,length), dtype=np.bool)
        np.ndarray[np.uint8_t, cast=True, ndim=1] boolOneC = np.zeros((length,), dtype=np.bool)

    for i in range(length):
        boolMatrix[i] = ( matrix[i] != 0)
        boolMatrix[:, i] = (matrix[:, i] != 0)

    for i in range(length):
        if i % 100 == 0:
            sys.stdout.write('\r {0}/{1} completed...'.format(i+1, length))
            sys.stdout.flush()

        for j from start <= j < i by step:

            # Get index of common
            idx = []
            for k from start <= k < length by step:
                if boolMatrix[i, k] != 0 and boolMatrix[j, k] != 0:
                    idx.append(k)

            if len(idx) > 3:
                rvec = np.corrcoef(matrix[i, idx], matrix[j, idx] )
                corrMatrix[i,j] = rvec[0, 1]
                corrMatrix[j,i] = rvec[0, 1]
            else:
                corrMatrix[i, j] = 0
                corrMatrix[j, i] = 0

    return corrMatrix
