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
import tempfile
import logging

import scipy.sparse as sparse
# import scipy.sparse.linalg as sparse_linalg
from scipy.linalg import get_blas_funcs
from scipy import linalg as sp_linalg
import string, random
import os
import time

from . import ccmapHelpers as cmh
from . import ccmap as cmp
from . import gcmap as gmp
from . import normalizeKnightRuiz as krnorm
from . import normalizeAverageContact as avcnorm
from . import normalizeIC as icnorm

dtype_npBINarray = 'float32'

logger = logging.getLogger('normalizer')
logger.setLevel(logging.INFO)


def NormalizeKnightRuizOriginal(ccMapObj, tol=1e-12, x0=None, delta=0.1, Delta=3, fl=0):
	'''Original Knight-Ruiz algorithm for matrix balancing

	Ported from a matlab script given in the supporting information of the following paper:
		* P.A. Knight and D. Ruiz (2013). A fast algorithm for matrix balancing (2013). IMA Journal of Numerical Analysis, 33, 1029-1047"
		* Matrix must be symmetric and non-negative
		* For input matrix A, this function find a vector X such that diag(X)*A*diag(X) is close to doubly stochastic.

	.. warning::
		* This is original ported code and kept here for comparison and testing.
		* Do not use it because for large matrix it may end up with consuming all the memory for large matrix.

	Parameters
	----------
	ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		A CCMAP object containing observed contact frequency

	Returns
	-------
	normCCMap : :class:`CCMAP`
		Normalized Contact map.

	'''

	normCCMap = ccMapObj.copy()
	normCCMap.make_editable()
	ccMapObj.make_readable()

	normCCMap.bNoData = np.all( ccMapObj.matrix == 0.0, axis=0)
	bNonZeros = ~normCCMap.bNoData
	A = (ccMapObj.matrix[bNonZeros,:])[:,bNonZeros]   # Selected row-column which are not all zeros

	BinNormMatrix = np.memmap(normCCMap.path2matrix, dtype=dtype_npBINarray, mode='w+', shape=normCCMap.shape)

	n = A.shape[0]                                                                 # n = size(A,1)
	e = np.ones((n, 1))                                                                         # e = ones(n,1)
	res = []
	if x0 is None:
		x0 = e


	g = 0.9        # Parameters used in inner stopping criterion.
	etamax = 0.1   # Parameters used in inner stopping criterion.
	eta = etamax
	stop_tol = tol * 0.5
	x = x0
	rt = tol**2                                 # rt = tol^2
	v = x * A.dot(x)                            # v = x.*(A*x)
	rk = 1 - v
	rho_km1 = np.dot( rk.conjugate().T, rk)     # rho_km1 = rk'*rk
	rout = rho_km1
	rold = rout

	# x, x0, e, v, rk, y, Z, w, p, ap :     vector shape(n, 1) : [ [value] [value] [value] [value] ... ... ... [value] ]
	# rho_km1, rout, rold, innertol, alpha :  scalar shape(1 ,1) : [[value]]

	MVP = 0                                                      # Well count matrix vector products.
	i = 0                                                        # Outer iteration count.

	if fl == 1:
		print('it in. it res')

	while rout > rt: 										     # Outer iteration
		i = i + 1
		k = 0
		y = e
		innertol = max( [ eta**2 * rout, rt ])                   # innertol = max([eta^2*rout,rt]);

		while rho_km1 > innertol:                                # Inner iteration by CG
			k = k + 1
			if k == 1:
				Z = rk / v                                       # Z = rk./v
				p = Z
				rho_km1 = np.dot( rk.conjugate().T, Z)           # rho_km1 = rk'*Z
			else:
				beta = rho_km1 / rho_km2
				p = Z + (beta * p)


			# Update search direction efficiently.
			w = x * A.dot((x*p)) + (v * p)                       # w = x.*(A*(x.*p)) + v.*p
			alpha = rho_km1 / np.dot( p.conjugate().T, w)        # alpha = rho_km1/(p'*w)
			ap = alpha * p                                       # ap = alpha*p (No dot function as alpha is scalar)

			# Test distance to boundary of cone.
			ynew = y + ap;
			#print(i, np.amin(ynew), delta, np.amin(ynew) <= delta)
			#print(i, np.amax(ynew), Delta, np.amax(ynew) >= Delta)
			if np.amin(ynew) <= delta:
				if delta == 0:
					break
				ind = np.nonzero(ap < 0) # ind = find(ap < 0)
				gamma = np.amin( (delta - y[ind]) / ap[ind] )    # gamma = min((delta - y(ind))./ap(ind))
				y = y + np.dot(gamma, ap)                        # y = y + gamma*ap
				break
			if np.amax(ynew) >= Delta:
				ind = np.nonzero( ynew > Delta )                 # ind = find(ynew > Delta);
				gamma = np.amin( (Delta-y[ind]) / ap[ind])       # gamma = min((Delta-y(ind))./ap(ind));
				y = y + np.dot(gamma, ap)                        # y = y + gamma*ap;
				break
			y = ynew
			rk = rk - alpha*w                                    # rk = rk - alpha*w
			rho_km2 = rho_km1
			Z = rk / v
			rho_km1 = np.dot( rk.conjugate().T, Z)               # rho_km1 = rk'*Z


		x = x * y                                                # x = x.*y
		v = x * A.dot(x)                                         # v = x.*(A*x)
		rk = 1 - v
		rho_km1 = np.dot( rk.conjugate().T, rk)                  # rho_km1 = rk'*rk
		rout = rho_km1
		MVP = MVP + k + 1

		# Update inner iteration stopping criterion.
		rat = rout/rold
		rold = rout
		res_norm = np.sqrt(rout)
		eta_o = eta
		eta = g*rat

		#print(i, res_norm)

		if g*eta_o**2 > 0.1:
			eta = np.amax([eta, g*eta_o**2])                    # eta = max([eta,g*eta_o^2])

		eta = np.amax([np.amin([eta, etamax]), stop_tol/res_norm]);   # eta = max([min([eta,etamax]),stop_tol/res_norm]);

		if fl == 1:
			print('%3d %6d %.3e %.3e %.3e \n' % (i,k, res_norm, np.amin(y), np.amin(x)))
			res=[res, res_norm]


	# Generation of Doubly stochastic matrix ( diag(X)*A*diag(X) )
	(fd, path2matrix) = tempfile.mkstemp(suffix='.bin', prefix='nparray_', dir=None, text=False)
	os.close(fd)     # Close file, error in windows OS
	A_DSMat = np.memmap(path2matrix, dtype=dtype_npBINarray, mode='w+', shape=A.shape)
	A_DSMat[:] = x.T * (A * x)

	# Assigning minvalue and maxvalue
	normCCMap.maxvalue = np.float(np.amax(A_DSMat))
	minvalue = np.amin(A_DSMat)
	v_steps = np.linspace(minvalue, normCCMap.maxvalue, 100)
	normCCMap.minvalue = minvalue - (v_steps[1] - v_steps[0])

	if (normCCMap.minvalue < 0.0):
		normCCMap.minvalue = 0.0

	dsm_i = -1
	dsm_j = 0
	for i in range(BinNormMatrix.shape[0]):
		if not normCCMap.bNoData[i]:
			dsm_i += 1

		dsm_j = 0
		for j in range(BinNormMatrix.shape[1]):
			if normCCMap.bNoData[i] or normCCMap.bNoData[j]:
				BinNormMatrix[i][j] = normCCMap.minvalue
			else:
				BinNormMatrix[i][j] = A_DSMat[dsm_i][dsm_j]
				dsm_j += 1


	# To check if sum of rows and columns are one
	r_sum = A_DSMat.sum(axis = 0)
	c_sum = A_DSMat.sum(axis = 1)
	#for i in range(A_DSMat.shape[0]):
		#print(i, x[i], r_sum[i], c_sum[i])


	BinNormMatrix.flush()
	del BinNormMatrix
	del A_DSMat

	if os.path.isfile(path2matrix):
		os.remove(path2matrix)

	return normCCMap

def normalizeCCMapByKR(ccMap, memory='RAM', tol=1e-12, outFile=None, vmin=None, vmax=None, percentile_thershold_no_data=None, thershold_data_occup=None, workDir=None):
	"""Normalize a ccmap using Knight-Ruiz matrix balancing method.

	.. note::
		* This function uses a modified version of orginal ported code given in :meth:`NormalizeKnightRuizOriginal`.
		* **Please refer to:** P.A. Knight and D. Ruiz (2013). A fast algorithm for matrix balancing (2013). IMA Journal of Numerical Analysis, 33, 1029-1047

	Parameters
	----------
	ccMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ccmap file
		A CCMAP object containing observed contact frequency or a ccmap file.

	memory : str
		Accepted keywords are ``RAM`` and ``HDD``:

			* ``RAM``: All intermediate arrays are generated in memory(RAM). This version is faster, however, it requires RAM depending on the input matrix size.
			* ``HDD``: All intermediate arrays are generated as memory mapped array files on hard-disk.

	tol : float
		Tolerance for matrix balancing. Smaller tolreance increases accuracy in sums of rows and columns.

	outFile : str
		Name of output ccmap file, to save directly the normalized map as a ccmap file. In case of this option, ``None`` will return.

	vmin : float
		Minimum thershold value for normalization. If contact frequency is less than or equal to this thershold value, this value is discarded during normalization.

	vmax : float
		Maximum thershold value for normalization. If contact frequency is greater than or equal to this thershold value, this value is discarded during normalization.

	percentile_thershold_no_data : int
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		``percentile_thershold_no_data`` should be between 1 and 100. This options discard the rows and columns which are above this percentile.
		For example: if this value is 99, those row or columns will be discarded which contains larger than number of zeros (missing data) at 99 percentile.

		To calculate percentile, all blank rows are removed, then in all rows, number of zeros are counted. Afterwards, number of zeros at
		`percentile_thershold_no_data` percentile is obtained. In next step, if a row contain number of zeros larger than this percentile value,
		the whole row and column is assigned to have missing data. This percentile indicates highest numbers of zeros (missing data) in given rows/columns.

	thershold_data_occup : float
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		This ratio is (number of bins with data) / (total number of bins in the given row/column).
		For example: if `thershold_data_occup = 0.8`, then all rows containing more than 20\% of missing data will be discarded.

		Note that this parameter is suitable for low resolution data because maps are likely to be much less sparse.

	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.


	Returns
	-------
	ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ``None``
		Normalized Contact map. When ``outFile`` is provided, ``None`` is returned. In case of any other error, ``None`` is returned.

	"""

	# Check whether input is a file or a obejct
	ccMapObjOrig, ccmapType = cmp.checkCCMapObjectOrFile(ccMap, workDir=workDir)

	# Make another copy here for maximum and minimum thershold value
	if vmin is not None or vmax is not None:
		ccMapObj = ccMapObjOrig.copy()
		ccMapObj.make_editable()

		if ccmapType == 'File':
			del ccMapObjOrig

		ccmapType = 'File' # This temporary file should be deleted when neccessary
		if vmin is not None:
			ccMapObj.matrix[ np.nonzero(ccMapObj.matrix <= vmin) ] = 0.0
		if vmax is not None:
			ccMapObj.matrix[ np.nonzero(ccMapObj.matrix >= vmax) ] = 0.0
		ccMapObj.matrix.flush()
		ccMapObj.make_unreadable()

	else:
		ccMapObj = ccMapObjOrig

	normCCMap = ccMapObj.copy(fill=0.0)
	normCCMap.make_editable()
	ccMapObj.make_readable()

	logger.info(' KR Normalization will be done through {0}.'.format(memory))
	if ccMapObj.xlabel is not None:
		logger.info(' KR Normalization is in process for {0} map...'.format(ccMapObj.xlabel))
	else:
		logger.info(' KR Normalization is in process for UNKNOWN map!!')

	# Default Parameters
	delta=0.1
	Delta=3
	fl=0

	# Main section
	try:

		if memory=='RAM':
			bNonZeros = cmh.get_nonzeros_index(ccMapObj.matrix, thershold_percentile=percentile_thershold_no_data, thershold_data_occup=thershold_data_occup)

			A = (ccMapObj.matrix[bNonZeros,:])[:,bNonZeros]   # Selected row-column which are not all zeros

			normCCMap.bNoData = ~bNonZeros
			del bNonZeros

			# Added one so that all value with zero will be one
			KRObj = krnorm.KnightRuizNorm(A, tol=tol, delta=delta, Delta=Delta, workDir=workDir)
			minvalue, maxvalue = KRObj.run(A, fl, normCCMap.matrix, normCCMap.bNoData)

		elif memory == 'HDD':

			# Selected row-column which are not all zeros
			A, normCCMap.bNoData  = cmh.remove_zeros(ccMapObj.matrix, thershold_percentile=percentile_thershold_no_data,
			                                            thershold_data_occup=thershold_data_occup, workDir=workDir)

			# Try for removing temporary files related to above variable A
			try:
				KRObj = krnorm.KnightRuizNorm(A.arr, memory='HDD', tol=tol, delta=delta, Delta=Delta, workDir=workDir)
				minvalue, maxvalue = KRObj.run(A.arr, fl, normCCMap.matrix, normCCMap.bNoData)
			except(KeyboardInterrupt, SystemExit) as e:
				del A
				raise e
			del A

		else:
			raise ValueError ('This memory={0} is not is not understandable.. Please use \'RAM\' or \'HDD\'.' .format(memory))

		normCCMap.maxvalue = float(maxvalue)
		normCCMap.minvalue = float(minvalue)

		normCCMap.make_unreadable()
		ccMapObj.make_unreadable()

		if ccMapObj.xlabel is not None:
			logger.info(' 	...Finished KR Normalization for {0} map...'.format(ccMapObj.xlabel))
		else:
			logger.info(' 	...Finished KR Normalization for UNKNOWN map!!')


		# When outFile is provided write a output
		if outFile is not None:
			cmp.save_ccmap(normCCMap, outFile, compress=True)

		# Delete ccmap object if input was a file
		if ccmapType == 'File':
			del ccMapObj

		# Whether outFile is given.
		if outFile is None:
			return normCCMap
		else:
			del normCCMap
			return None

	# In case of program termination, delete the newly created ccmap and raise error
	except (KeyboardInterrupt, SystemExit) as e:
		if 'normCCMap' in locals():	del normCCMap
		if ccmapType == 'File' and 'ccMapObj' in locals():	del ccMapObj
		raise e

	# In case of other error, delete the newly created ccmap, print error and return None
	except Exception as e:
		if 'normCCMap' in locals():	del normCCMap
		if ccmapType == 'File' and 'ccMapObj' in locals():	del ccMapObj
		logger.warning(e)
		logger.warning('Error in normalizing map!!!')
		raise e

def normalizeGCMapByKR(gcMapInputFile, gcMapOutFile, mapSizeCeilingForMemory=20000, vmin=None, vmax=None, tol=1e-12, percentile_thershold_no_data=None, thershold_data_occup=None, compression='lzf', workDir=None, logHandler=None):
	"""Normalize a gcmap using Knight-Ruiz matrix balancing method.


	Parameters
	----------
	gcMapInputFile : str
		Name of input gcmap file.

	gcMapOutFile : str
		Name of output gcmap file.

	mapSizeCeilingForMemory : int
		Maximum size of contact map allowed for calculation using RAM. If map size or shape is larger than this value,
		normalization will be performed using disk (HDD).

	vmin : float
		Minimum thershold value for normalization. If contact frequency is less than or equal to this thershold value, this value is discarded during normalization.

	vmax : float
		Maximum thershold value for normalization. If contact frequency is greater than or equal to this thershold value, this value is discarded during normalization.

	tol : float
		Tolerance for matrix balancing. Smaller tolreance increases accuracy in sums of rows and columns.

	percentile_thershold_no_data : int
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		``percentile_thershold_no_data`` should be between 1 and 100. This options discard the rows and columns which are above this percentile.
		For example: if this value is 99, those row or columns will be discarded which contains larger than number of zeros (missing data) at 99 percentile.

		To calculate percentile, all blank rows are removed, then in all rows, number of zeros are counted. Afterwards, number of zeros at
		`percentile_thershold_no_data` percentile is obtained. In next step, if a row contain number of zeros larger than this percentile value,
		the whole row and column is assigned to have missing data. This percentile indicates highest numbers of zeros (missing data) in given rows/columns.

	thershold_data_occup : float
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		This ratio is (number of bins with data) / (total number of bins in the given row/column).
		For example: if `thershold_data_occup = 0.8`, then all rows containing more than 20\% of missing data will be discarded.

		Note that this parameter is suitable for low resolution data because maps are likely to be much less sparse.

	compression : str
	    Compression method in output gcmap file. Presently allowed : ``lzf`` for LZF compression and ``gzip`` for GZIP compression.

	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.


	Returns
	-------
	None


	.. seealso::
		:meth:`gcMapExplorer.lib.normalizer.normalizeCCMapByKR`

	"""

	# Get list of maps in ascending order
	gcmap = gmp.GCMAP(gcMapInputFile)
	gcmap.loadSmallestMap()
	mapList = gcmap.mapNameList.copy()
	del gcmap

	for mapName in mapList:
		ccMap = gmp.loadGCMapAsCCMap(gcMapInputFile, mapName=mapName, workDir=workDir)

		# Because ccMap is already loaded here, directly edit matrix here
		# No need to pass vmin and vmx during normalization as
		# it is already taken care here
		if vmin is not None or vmax is not None:
			ccMap.make_editable()
			if vmin is not None:
				ccMap.matrix[ np.nonzero(ccMap.matrix <= vmin) ] = 0.0
			if vmax is not None:
				ccMap.matrix[ np.nonzero(ccMap.matrix >= vmax) ] = 0.0
			ccMap.matrix.flush()
			ccMap.make_unreadable()

		try:

			if ccMap.shape[0] > mapSizeCeilingForMemory:
				memory = 'HDD'
			else:
				memory = 'RAM'

			norm_ccmap = normalizeCCMapByKR(ccMap, memory=memory, tol=tol,
						percentile_thershold_no_data=percentile_thershold_no_data,
						thershold_data_occup=thershold_data_occup,
						workDir=workDir)

			if norm_ccmap is not None:
				gmp.addCCMap2GCMap(norm_ccmap, gcMapOutFile,
									compression=compression,
									generateCoarse=True, coarsingMethod='sum',
									logHandler=logHandler)

		# In case of program termination, delete the newly created ccmap and raise error
		except (KeyboardInterrupt, SystemExit) as e:
			if 'ccMap' in locals():	del ccMap
			if 'norm_ccmap' in locals():	del norm_ccmap
			if 'gcmap' in locals():	del gcmap
			raise e

		del ccMap
		del norm_ccmap

def normalizeCCMapByIC(ccMap, tol=1e-4, vmin=None, vmax=None, outFile=None, iteration=500, percentile_thershold_no_data=None, thershold_data_occup=None, workDir=None):
	""" Normalize a ccmap by Iterative correction method

	This method normalize the raw contact map by removing biases from experimental procedure.
	For more details, see `this publication <http://www.nature.com/nmeth/journal/v9/n10/full/nmeth.2148.html>`_.

	Parameters
	----------
	ccMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ccmap file.
		A CCMAP object containing observed contact frequency or a ccmap file

	tol : float
		Tolerance value. The relative increment in the results before declaring convergence.

	vmin : float
		Minimum thershold value for normalization. If contact frequency is less than or equal to this thershold value, this value is discarded during normalization.

	vmax : float
		Maximum thershold value for normalization. If contact frequency is greater than or equal to this thershold value, this value is discarded during normalization.

	outFile : str
		Name of output ccmap file, to save directly the normalized map as a ccmap file. In case of this option, ``None`` will return.

	iteration : int
		Number of iteration to stop the normalization.

	percentile_thershold_no_data : int
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		``percentile_thershold_no_data`` should be between 1 and 100. This options discard the rows and columns which are above this percentile.
		For example: if this value is 99, those row or columns will be discarded which contains larger than number of zeros (missing data) at 99 percentile.

		To calculate percentile, all blank rows are removed, then in all rows, number of zeros are counted. Afterwards, number of zeros at
		`percentile_thershold_no_data` percentile is obtained. In next step, if a row contain number of zeros larger than this percentile value,
		the whole row and column is assigned to have missing data. This percentile indicates highest numbers of zeros (missing data) in given rows/columns.

	thershold_data_occup : float
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		This ratio is (number of bins with data) / (total number of bins in the given row/column).
		For example: if `thershold_data_occup = 0.8`, then all rows containing more than 20\% of missing data will be discarded.

		Note that this parameter is suitable for low resolution data because maps are likely to be much less sparse.

	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.

	Returns
	-------
	normCCMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ``None``
		Normalized Contact map. When ``outFile`` is provided, ``None`` is returned. In case of any other error, ``None`` is returned.

	"""

	# Check whether input is a file or a obejct
	ccMapObj, ccmapType = cmp.checkCCMapObjectOrFile(ccMap, workDir=workDir)

	tmap = ccMapObj.copy()

	# In case if vmin and vmax is given
	if vmin is not None or vmax is not None:
		tmap.make_editable()
		if vmin is not None:
			tmap.matrix[ np.nonzero(tmap.matrix <= vmin) ] = 0.0
		if vmax is not None:
			tmap.matrix[ np.nonzero(tmap.matrix >= vmax) ] = 0.0
		tmap.make_unreadable()

	if ccMapObj.xlabel is not None:
		logger.info(' Iterative Correction is in process for {0} map...'.format(ccMapObj.xlabel))
	else:
		logger.info(' Iterative Correction is in process for UNKNOWN map!!')

	# Check if any filter is applied
	applyFilter = True
	if percentile_thershold_no_data is None and thershold_data_occup is None:
	    applyFilter = False

	try:
		tmap.make_editable()

		matrix = None
		bNonZeros = None
		if applyFilter:
		    bNonZeros = cmh.get_nonzeros_index(tmap.matrix, thershold_percentile=percentile_thershold_no_data, thershold_data_occup=thershold_data_occup)
		    matrix = (tmap.matrix[bNonZeros,:])[:,bNonZeros]   # Selected row-column which are not all zeros
		    icnorm.performIterativeCorrection(matrix, tol, iteration)
		else:
			matrix = tmap.matrix
			icnorm.performIterativeCorrection(matrix, tol, iteration)
			bNonZeros = ~np.all( tmap.matrix == 0.0, axis=0)

		ma = np.ma.masked_equal(matrix, 0.0, copy=False)
		tmap.minvalue = ma.min()
		tmap.maxvalue = ma.max()
		tmap.bNoData = ~bNonZeros

		if applyFilter:
		    # Fill original output matrix
		    dsm_i = 0
		    ox = tmap.matrix.shape[0]
		    idx_fill = np.nonzero( bNonZeros )
		    for i in range(ox):
		        if not tmap.bNoData[i]:
		            tmap.matrix[i, idx_fill] = matrix[dsm_i]
		            tmap.matrix[idx_fill, i] = matrix[dsm_i]
		            dsm_i += 1
		        else:
		            idx_nonzero = np.nonzero( tmap.matrix[i] > 0 )
		            tmap.matrix[i][idx_nonzero].fill(0.0)

		if ccMapObj.xlabel is not None:
			logger.info(' 	...Finished Iterative Correction for {0} map...'.format(ccMapObj.xlabel))
		else:
			logger.info(' 	...Finished Iterative Correction for UNKNOWN map!!')

		# Save output ccmap file
		if outFile is not None:
			cmp.save_ccmap(tmap, outFile, compress=True)

		# Delete ccmap object if input was a file
		if ccmapType == 'File':
			del ccMapObj

		# Whether outFile is given.
		if outFile is None:
			return tmap
		else:
			del tmap
			return None

	# In case of program termination, delete the newly created ccmap and raise error
	except (SystemExit, KeyboardInterrupt) as e:
		if 'tmap' in locals():	del tmap
		if ccmapType == 'File' and 'ccMapObj' in locals():	del ccMapObj
		raise e

	# In case of other error, delete the newly created ccmap, print error and return None
	except Exception as e:
		if 'tmap' in locals():	del tmap
		if ccmapType == 'File' and 'ccMapObj' in locals():	del ccMapObj
		logger.warning(e)
		logger.warning('Error in Iterative Correction!!!')
		return None

def normalizeGCMapByIC(gcMapInputFile, gcMapOutFile, vmin=None, vmax=None, tol=1e-12, iteration=500, percentile_thershold_no_data=None, thershold_data_occup=None, compression='lzf', workDir=None, logHandler=None):
	"""Normalize a gcmap using Iterative Correction.

	This method normalize the raw contact map by removing biases from experimental procedure.
	For more details, see `this publication <http://www.nature.com/nmeth/journal/v9/n10/full/nmeth.2148.html>`_.

	Parameters
	----------
	gcMapInputFile : str
		Name of input gcmap file.

	gcMapOutFile : str
		Name of output gcmap file.

	vmin : float
		Minimum thershold value for normalization. If contact frequency is less than or equal to this thershold value, this value is discarded during normalization.

	vmax : float
		Maximum thershold value for normalization. If contact frequency is greater than or equal to this thershold value, this value is discarded during normalization.

	tol : float
		Tolerance value. The relative increment in the results before declaring convergence.

	iteration : int
		Number of iteration to stop the normalization.

	percentile_thershold_no_data : int
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		``percentile_thershold_no_data`` should be between 1 and 100. This options discard the rows and columns which are above this percentile.
		For example: if this value is 99, those row or columns will be discarded which contains larger than number of zeros (missing data) at 99 percentile.

		To calculate percentile, all blank rows are removed, then in all rows, number of zeros are counted. Afterwards, number of zeros at
		`percentile_thershold_no_data` percentile is obtained. In next step, if a row contain number of zeros larger than this percentile value,
		the whole row and column is assigned to have missing data. This percentile indicates highest numbers of zeros (missing data) in given rows/columns.

	thershold_data_occup : float
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		This ratio is (number of bins with data) / (total number of bins in the given row/column).
		For example: if `thershold_data_occup = 0.8`, then all rows containing more than 20\% of missing data will be discarded.

		Note that this parameter is suitable for low resolution data because maps are likely to be much less sparse.

	compression : str
	    Compression method in output gcmap file. Presently allowed : ``lzf`` for LZF compression and ``gzip`` for GZIP compression.

	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.


	Returns
	-------
	None


	.. seealso::
		:meth:`gcMapExplorer.lib.normalizer.normalizeCCMapByIC`


	"""

	# Get list of maps in ascending order
	gcmap = gmp.GCMAP(gcMapInputFile)
	gcmap.loadSmallestMap()
	mapList = gcmap.mapNameList.copy()
	del gcmap

	for mapName in mapList:
		ccMap = gmp.loadGCMapAsCCMap(gcMapInputFile, mapName=mapName, workDir=workDir)

		try:
			norm_ccmap = normalizeCCMapByIC(ccMap, tol=tol, iteration=iteration,
						vmin=vmin, vmax=vmax,
						percentile_thershold_no_data=percentile_thershold_no_data,
						thershold_data_occup=thershold_data_occup,
						workDir=workDir)

			if norm_ccmap is not None:
				gmp.addCCMap2GCMap(norm_ccmap, gcMapOutFile,
									compression=compression,
									generateCoarse=True, coarsingMethod='sum',
									logHandler=logHandler)

		# In case of program termination, delete the newly created ccmap and raise error
		except (KeyboardInterrupt, SystemExit) as e:
			if 'ccMap' in locals():	del ccMap
			if 'norm_ccmap' in locals():	del norm_ccmap
			if 'gcmap' in locals():	del gcmap
			raise e

		del ccMap
		del norm_ccmap

def normalizeCCMapByMCFS(ccMap, stats='median', vmin=None, vmax=None, stype='o/e', outFile=None, percentile_thershold_no_data=None, thershold_data_occup=None, workDir=None):
	""" Scale ccmap using Median Contact Frequency

	This method can be used to normalize contact map with expected values.
	These expected values could be either Median or Average contact values
	for particular distance between two locations/coordinates. At first,
	Median/Average distance contact frequency for each distance is calculated.
	Subsequently, the observed contact frequency is either divided ('o/e') or
	substracted ('o-e') by median/average contact frequency obtained for
	distance between the two locations.

	.. note:
		In place of median, mean can be also used for normalization. See below for options.

	Parameters
	----------
	ccMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ccmap file
		A CCMAP object containing observed contact frequency or a ccmap file

	stats : str
	    Statistics to be calculated along diagonals: It may be either "mean" or "median". By default, it is "median".

	vmin : float
		Minimum thershold value for normalization. If contact frequency is less than or equal to this thershold value, this value is discarded during normalization.

	vmax : float
		Maximum thershold value for normalization. If contact frequency is greater than or equal to this thershold value, this value is discarded during normalization.

	stype : str
		Type of scaling. It may be either 'o/e' or 'o-e'. In case of 'o/e',
		Observed/Expected will be calculated while (Observed - Expected)
		will be calculated for 'o-e'.

	outFile : str
		Name of output ccmap file, to save directly the normalized map as a ccmap file. In case of this option, ``None`` will return.

	percentile_thershold_no_data : int
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		``percentile_thershold_no_data`` should be between 1 and 100. This options discard the rows and columns which are above this percentile.
		For example: if this value is 99, those row or columns will be discarded which contains larger than number of zeros (missing data) at 99 percentile.

		To calculate percentile, all blank rows are removed, then in all rows, number of zeros are counted. Afterwards, number of zeros at
		`percentile_thershold_no_data` percentile is obtained. In next step, if a row contain number of zeros larger than this percentile value,
		the whole row and column is assigned to have missing data. This percentile indicates highest numbers of zeros (missing data) in given rows/columns.

	thershold_data_occup : float
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		This ratio is (number of bins with data) / (total number of bins in the given row/column).
		For example: if `thershold_data_occup = 0.8`, then all rows containing more than 20\% of missing data will be discarded.

		Note that this parameter is suitable for low resolution data because maps are likely to be much less sparse.

	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.

	Returns
	-------
	ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ``None``
		Normalized Contact map. When ``outFile`` is provided, ``None`` is returned. In case of any other error, ``None`` is returned.

	"""

	if stype not in ['o/e', 'o-e']:
		logger.warning('Wrong value {0} for stype: {1}.'.format(stype, ['o/e', 'o-e']))
		return None

	# Check whether input is a file or a obejct
	ccMapObjOrig, ccmapType = cmp.checkCCMapObjectOrFile(ccMap, workDir=workDir)

	# Make another copy here for maximum and minimum thershold value
	if vmin is not None or vmax is not None:
		ccMapObj = ccMapObjOrig.copy()
		ccMapObj.make_editable()

		if ccmapType == 'File':
			del ccMapObjOrig

		ccmapType = 'File' # This temporary file should be deleted when neccessary
		if vmin is not None:
			ccMapObj.matrix[ np.nonzero(ccMapObj.matrix <= vmin) ] = 0.0
		if vmax is not None:
			ccMapObj.matrix[ np.nonzero(ccMapObj.matrix >= vmax) ] = 0.0
		ccMapObj.matrix.flush()
		ccMapObj.make_unreadable()

	else:
		ccMapObj = ccMapObjOrig

	normCCMap = ccMapObj.copy()

	if ccMapObj.xlabel is not None:
		logger.info(' Median Contact Frequency Scaling is in process for {0} map...'.format(ccMapObj.xlabel))
	else:
		logger.info(' Median Contact Frequency Scaling is in process for UNKNOWN map!!')

	# Main section
	try:
		if stype == 'o/e':
			avcnorm._normalizeByAvgContactByDivision(ccMapObj, normCCMap, stats=stats,
						percentile_thershold_no_data=percentile_thershold_no_data,
						thershold_data_occup=thershold_data_occup)

		if stype == 'o-e':
			avcnorm._normalizeByAvgContactBySubstraction(ccMapObj, normCCMap, stats=stats,
						percentile_thershold_no_data=percentile_thershold_no_data,
						thershold_data_occup=thershold_data_occup)

		if ccMapObj.xlabel is not None:
			logger.info(' 	...Finished Median Contact Frequency Scaling for {0} map...'.format(ccMapObj.xlabel))
		else:
			logger.info(' 	...Finished Median Contact Frequency Scaling for UNKNOWN map!!')

		# Save output ccmap file
		if outFile is not None:
			cmp.save_ccmap(normCCMap, outFile, compress=True)

		# Delete ccmap object if input was a file
		if ccmapType == 'File':
			del ccMapObj

		# Whether outFile is given.
		if outFile is None:
			return normCCMap
		else:
			del normCCMap
			return None

	# In case of program termination, delete the newly created ccmap and raise error
	except (KeyboardInterrupt, SystemExit) as e:
		if 'normCCMap' in locals():	del normCCMap
		if ccmapType == 'File' and 'ccMapObj' in locals():	del ccMapObj
		raise e

	# In case of other error, delete the newly created ccmap, print error and return None
	except Exception as e:
		if 'normCCMap' in locals():	del normCCMap
		if ccmapType == 'File' and 'ccMapObj' in locals():	del ccMapObj
		logger.warning(e)
		logger.warning('Error in Median Contact Frequency Scaling!!!')
		return None

def normalizeGCMapByMCFS(gcMapInputFile, gcMapOutFile, stats='median', vmin=None, vmax=None, stype='o/e', percentile_thershold_no_data=None, thershold_data_occup=None, compression='lzf', workDir=None, logHandler=None):
	""" Scale all maps in gcmap using Median Contact Frequency

	This method can be used to normalize contact map with expected values.
	These expected values could be either Median or Average contact values
	for particular distance between two locations/coordinates. At first,
	Median/Average distance contact frequency for each distance is calculated.
	Subsequently, the observed contact frequency is either divided ('o/e') or
	substracted ('o-e') by median/average contact frequency obtained for
	distance between the two locations.

	Parameters
	----------
	gcMapInputFile : str
		Name of input gcmap file.

	gcMapOutFile : str
		Name of output gcmap file.

	stats : str
	    Statistics to be calculated along diagonals: It may be either "mean" or "median". By default, it is "median".

	vmin : float
		Minimum thershold value for normalization. If contact frequency is less than or equal to this thershold value, this value is discarded during normalization.

	vmax : float
		Maximum thershold value for normalization. If contact frequency is greater than or equal to this thershold value, this value is discarded during normalization.

	stype : str
		Type of scaling. It may be either 'o/e' or 'o-e'. In case of 'o/e',
		Observed/Expected will be calculated while (Observed - Expected)
		will be calculated for 'o-e'.

	percentile_thershold_no_data : int
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		``percentile_thershold_no_data`` should be between 1 and 100. This options discard the rows and columns which are above this percentile.
		For example: if this value is 99, those row or columns will be discarded which contains larger than number of zeros (missing data) at 99 percentile.

		To calculate percentile, all blank rows are removed, then in all rows, number of zeros are counted. Afterwards, number of zeros at
		`percentile_thershold_no_data` percentile is obtained. In next step, if a row contain number of zeros larger than this percentile value,
		the whole row and column is assigned to have missing data. This percentile indicates highest numbers of zeros (missing data) in given rows/columns.

	thershold_data_occup : float
		It can be used to filter the map, where rows/columns with largest numbers of missing data can be discarded.
		This ratio is (number of bins with data) / (total number of bins in the given row/column).
		For example: if `thershold_data_occup = 0.8`, then all rows containing more than 20\% of missing data will be discarded.

		Note that this parameter is suitable for low resolution data because maps are likely to be much less sparse.

	compression : str
	    Compression method in output gcmap file. Presently allowed : ``lzf`` for LZF compression and ``gzip`` for GZIP compression.

	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.


	Returns
	-------
	None


	.. seealso::
		:meth:`gcMapExplorer.lib.normalizer.normalizeCCMapByMCFS`


	"""

	if stype not in ['o/e', 'o-e']:
		raise ValueError('Wrong value {0} for stype: {1}.'.format(stype, ['o/e', 'o-e']))

	# Get list of maps in ascending order
	gcmap = gmp.GCMAP(gcMapInputFile)
	gcmap.loadSmallestMap()
	mapList = gcmap.mapNameList.copy()

	for mapName in mapList:

		# Iterate over available resolutions
		gcmap.changeMap(mapName)
		previousResolution = gcmap.resolution

		while True:
			ccMap = gmp.loadGCMapAsCCMap(gcmap.hdf5, mapName=mapName, resolution=gcmap.resolution, workDir=workDir)

			# Because ccMap is already loaded here, directly edit matrix here
			# No need to pass vmin and vmx during normalization as
			# it is already taken care here
			if vmin is not None or vmax is not None:
				ccMap.make_editable()
				if vmin is not None:
					ccMap.matrix[ np.nonzero(ccMap.matrix <= vmin) ] = 0.0
				if vmax is not None:
					ccMap.matrix[ np.nonzero(ccMap.matrix >= vmax) ] = 0.0
				ccMap.matrix.flush()
				ccMap.make_unreadable()

			try:
				norm_ccmap = normalizeCCMapByMCFS(ccMap, stats=stats, stype=stype,
							percentile_thershold_no_data=percentile_thershold_no_data,
							thershold_data_occup=thershold_data_occup,
							workDir=workDir)

				if norm_ccmap is not None:
					gmp.addCCMap2GCMap(norm_ccmap, gcMapOutFile,
										compression=compression,
										generateCoarse=False, replaceCMap=False,
										logHandler=logHandler)

					del norm_ccmap

				else:
					break

				gcmap.toCoarserResolution()
				if previousResolution == gcmap.resolution:
					break
				else:
					previousResolution = gcmap.resolution

			# In case of program termination, delete the newly created ccmap and raise error
			except (KeyboardInterrupt, SystemExit) as e:
				if 'ccMap' in locals():	del ccMap
				if 'norm_ccmap' in locals():	del norm_ccmap
				if 'gcmap' in locals():	del gcmap
				raise e

			del ccMap
