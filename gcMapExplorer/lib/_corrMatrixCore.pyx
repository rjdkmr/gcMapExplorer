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

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def _calculateCorrMatrix(np.ndarray[dtypeFloat_t, ndim=2] matrix):
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
        boolMatrix[:, i] = (matrix[i] != 0)

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
