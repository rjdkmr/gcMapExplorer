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

import os, re
import numpy as np
import scipy as sp

from . import gcmap as gmp
from .ccmap import checkCCMapObjectOrFile

from ._corrMatrixCore import _calculateCorrMatrix

def calculateCorrMatrix(inputCCMap, logspace=False, outFile=None, workDir=None):
	""" Calculate correlation matrix of a contact map.
	It calculates correlation between all rows and columns of contact map.


	Parameters
	----------
	ccMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ccmap file
		A CCMAP object containing observed contact frequency or a ccmap file.
	logspace : bool
		If its value is ``True``, at first map is converted as logarithm of map
		and subsequently correlation will be calculated.
	outFile : str
		Name of output ccmap file, to save directly the correlation matrix as a ccmap file. In case of this option, ``None`` will return.
	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.


	Returns
	-------
	ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ``None``
		Normalized Contact map. When ``outFile`` is provided, ``None`` is returned. In case of any other error, ``None`` is returned.

	"""

	# Check whether input is a file or a obejct
	ccmap, ccmapType = checkCCMapObjectOrFile(inputCCMap, workDir=workDir)

	try:
	    ccmap.make_readable()
	    matrix = (ccmap.matrix[~ccmap.bNoData,:])[:,~ccmap.bNoData]   # Selected row-column which are not all zeros
	    matrix = matrix.astype(np.float)
	    mshape = matrix.shape

	    if logspace:
	        matrix = np.log10(matrix)
	        matrix[np.isinf(matrix)] = 0

	    corrMatrix = _calculateCorrMatrix(matrix)
	    corrMatrix[np.isnan(corrMatrix)] = 0

	    corrCCMap = ccmap.copy(fill=0.0)
	    corrCCMap.make_editable()

	    dsm_i = 0
	    ox = corrCCMap.matrix.shape[0]
	    idx_fill = np.nonzero( ~(corrCCMap.bNoData) )
	    for i in range(ox):
	        if not corrCCMap.bNoData[i]:
	            corrCCMap.matrix[i, idx_fill] = corrMatrix[dsm_i]
	            corrCCMap.matrix[idx_fill, i] = corrMatrix[dsm_i]
	            dsm_i += 1

	    corrCCMap.matrix.flush()
	    corrCCMap.minvalue = -1
	    corrCCMap.maxvalue = 1

	except (KeyboardInterrupt, SystemExit) as e:
		# In case of program termination, delete the newly created ccmap and raise error
		if 'corrCCMap' in locals():	del corrCCMap
		if ccmapType == 'File' and 'ccmap' in locals():	del ccmap
		raise e

	# In case of other error, delete the newly created ccmap, print error and return None
	except Exception as e:
		if 'corrCCMap' in locals():	del corrCCMap
		if ccmapType == 'File' and 'ccmap' in locals():	del ccmap
		logger.warning(e)
		logger.warning('Error in Correlation Matrix Calculation !!!')
		return None

	# Save output ccmap file
	if outFile is not None:
		cmp.save_ccmap(corrCCMap, outFile, compress=True)
		del corrCCMap
		return
	else:
		return corrCCMap


def calculateCorrMatrixForGCMaps(gcMapInputFile, gcMapOutFile, logspace=False, replaceMatrix=False, compression='lzf', workDir=None):
	""" Calculate Correlation matrix for all maps present in input gcmap file
	It calculates correlation between all rows and columns of contact map.


	Parameters
	----------
	gcMapInputFile : str
	    Input gcmap file.
	gcMapOutFile : str
	    Output gcmap file.
	logspace : bool
		If its value is ``True``, at first map is converted as logarithm of map
		and subsequently correlation will be calculated.
	replaceMatrix : bool
		If its value is ``True``, the map will be replaced in output file.
		Otherwise, if a map is present, calculation will be skipped.
	compression : str
	    Compression method in output gcmap file. Presently allowed : ``lzf`` for LZF compression and ``gzip`` for GZIP compression.
	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.


	Returns
	-------
	None

	"""

	# Get list of maps in ascending order
	gcmap = gmp.GCMAP(gcMapInputFile)
	gcmap.loadSmallestMap()
	mapList = gcmap.mapNameList.copy()

	for mapName in mapList:
		# Iterate over available resolutions
		gcmap.changeMap(mapName)

		while True:
			# Check if map is present in output file
			# It check for all resolutions.
			tf = gmp.GCMAP(gcMapOutFile)
			doExist = False
			while True:
			    print('Calculating for {0}-{1}'.format(mapName, gcmap.resolution))
			    doExist = tf.checkMapExist(mapName=mapName, resolution=gcmap.resolution)
			    if doExist:
			        if not gcmap.toCoarserResolution():
			            break
			    else:
			        break

			# If map is present, either replace or skip according to input
			if doExist:
				if not replaceMatrix:
				    break
				else:
				    tf.hdf5[mapName].pop(gcmap.resolution)

			del tf

			ccMap = gmp.loadGCMapAsCCMap(gcmap.hdf5, mapName=mapName,
										resolution=gcmap.resolution,
		                                workDir=workDir)

			try:
				corrCCMap = calculateCorrMatrix(ccMap, logspace=logspace)
				if corrCCMap is not None:
					gmp.gcmap.addCCMap2GCMap(corrCCMap, gcMapOutFile,
												scaleoffset=4,
												compression=compression,
		                                       	generateCoarse=False,
												replaceCMap=False)
					del corrCCMap

				if not gcmap.toCoarserResolution():
					break

			# In case of program termination, delete the newly created ccmap and raise error
			except (KeyboardInterrupt, SystemExit) as e:
				if 'corrCCMap' in locals(): del corrCCMap
				if 'gcmap' in locals(): del gcmap
				if 'ccMap' in locals(): del ccMap
				raise e

			del ccMap
