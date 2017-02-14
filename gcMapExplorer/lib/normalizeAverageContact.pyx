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
from . import ccmapHelpers as cmh
from . import cmstats


def _normalizeByAvgContact(ccMapObj, normCCMap, stats='median', percentile_thershold_no_data=None, thershold_data_occup=None):

    normCCMap.make_editable()
    ccMapObj.make_readable()

    bNonZeros = cmh.get_nonzeros_index(ccMapObj.matrix, thershold_percentile=percentile_thershold_no_data, thershold_data_occup=thershold_data_occup)
    normCCMap.bNoData = ~bNonZeros

    ccMapObj.make_readable()
    avgContacts = cmstats.getAvgContactByDistance(ccMapObj, stats=stats)

    # Do along diagonal
    for m in range(ccMapObj.shape[0]):

        # do for postive diagonal offset
        idx = cmh.kth_diag_indices(m, ccMapObj.matrix)
        if avgContacts[m] != 0:
            normCCMap.matrix[idx] = ccMapObj.matrix[idx] / avgContacts[m]
        else:
            normCCMap.matrix[idx] = ccMapObj.matrix[idx]

        # Do for negative diagonal offset
        del idx
        if m != 0:
            idx = cmh.kth_diag_indices(m*-1, ccMapObj.matrix)
            if avgContacts[m] != 0:
                normCCMap.matrix[idx] = ccMapObj.matrix[idx] / avgContacts[m]
            else:
                normCCMap.matrix[idx] = ccMapObj.matrix[idx]
            del idx

    marray = np.ma.masked_equal(normCCMap.matrix, 0.0, copy=False)
    normCCMap.minvalue = marray.min()
    normCCMap.maxvalue = np.amax(normCCMap.matrix)

    return normCCMap
