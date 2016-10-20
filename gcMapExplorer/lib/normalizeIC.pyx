#!/usr/bin/env python3
#
# Author: Rajendra Kumar
#
# This file is part of gcMapExplorer
# Copyright (C) 2016  Rajendra Kumar, Ludvig Lizana, Per Stenberg
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


from . import ccmapHelpers as cmh


def performIterativeCorrection(matrix, tol, iteration):
    count = 0
    B = np.ones(matrix.shape[0], dtype=np.float)

    while True:
        dB = np.sum(matrix, axis=0)              # Sum over rows
        dB = dB/np.mean(dB[np.nonzero(dB != 0)])      # Renormalization by mean value
        dB[np.nonzero(dB == 0)] = 1                   # Set zeros value to one

        # Divide W_ij by dB_i dB_j
        for i in range(matrix.shape[0]):
            matrix[i] = matrix[i] / dB
            matrix[:, i] = matrix[:, i] / dB

        B = B *  dB

        if np.var(dB) < tol or count > iteration:
            break

        count = count + 1
