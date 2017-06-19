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
# import scipy.sparse.linalg as sparse_linalg
from scipy.linalg import get_blas_funcs
from scipy import linalg as sp_linalg
import random
import os, re
import logging

from gcMapExplorer.config import getConfig
from . import ccmap as cmp
from . import gcmap as gmp


# get configuration
config = getConfig()

logger = logging.getLogger('StatDist')
logger.setLevel(logging.INFO)

# Calling blas dot function for matrix
def mdot(v, w):
	gemm = get_blas_funcs("gemm", [v, w])
	return gemm(alpha=1., a=v, b=w)

def calculateTransitionProbablityMatrix(A, Out=None, bNoData=None):
	""" Core function to calculate transition probablity matrix.

	It is a core function to calculate transition probablity matrix.

	Parameters
	----------
	A : numpy.memmap or :attr:`gcMapExplorer.lib.ccmap.CCMAP.matrix` or numpy.ndarray
		Input map or matrix
	Out : numpy.memmap or :attr:`gcMapExplorer.lib.ccmap.CCMAP.matrix` or numpy.ndarray
		Ouput transition matrix. In case if it is ``None``, ouput matrix will be returned.
	bNoData : numpy.array[bool]
		1D-array containing ``True`` and ``False`` values. Its size should be
		equal to input array row/column size. Row/Column with ``False`` value
		will be considered during the calculation.

	Returns
	-------
	Out : numpy.ndarray or ``None``.
		Ouput transition matrix. In case if ``Out`` is passed, ``None`` will be returned.

	"""
	mshape = A.shape

	toReturn = False

	if Out is None:
		Out = np.zeros(A.shape, dtype=A.dtype)
		toReturn = True

	for m in range(mshape[0]):
		if bNoData is not None:
			if bNoData[m]:
				continue
		sumr = A[m].sum()
		Out[m] = (A[m] / sumr)[:]

	ma = np.ma.masked_equal(Out, 0.0, copy=False)
	minvalue = ma.min()

	for m in range(mshape[0]):
		if bNoData[m]:
			Out[m].fill(minvalue)
			Out[:,m].fill(minvalue)
		else:
			idx = np.nonzero( Out[m] == 0 )
			Out[m][idx] = minvalue

	if toReturn:
		return Out

def transitionProbablityMatrixForCCMap(ccMap,  outFile=None, percentile_thershold_no_data=None, thershold_data_occup=None):
	""" To calculate transition probablity matrix.

	This method can be used to calculate transition probablity matrix. This is similar to markov-chain transition matrix.

	:note: This transition matrix is not symmetric, because each row represents stochastic row vector,
		which contains contact probablity of this bin with every other bins and sum of row is always equal to one. See here: |markov chain link|.

	Parameters
	----------
	ccMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ccmap file
		A CCMAP object containing observed contact frequency or a ccmap file

	outFile : str
		Name of output ccmap file, to save directly the map as a ccmap file. In case of this option, ``None`` will return.

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

	Returns
	-------
	normCCMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ``None``
		Transition matrix. When ``outFile`` is provided, ``None`` is returned. In case of any other error, ``None`` is returned.

	"""

	# Check whether input is a file or a obejct
	ccMapObj, ccmapType = cmp.checkCCMapObjectOrFile(ccMap, workDir=workDir)

	normCCMap = ccMapObj.copy(fill=0.0)
	normCCMap.make_editable()
	ccMapObj.make_readable()

	bNoData = None
	if percentile_thershold_no_data is not None or thershold_data_occup is not None or normCCMap.bNoData is None:
		bNonZeros = cmh.get_nonzeros_index(ccMapObj.matrix, thershold_percentile=percentile_thershold_no_data, thershold_data_occup=thershold_data_occup)
		bNoData = ~bNonZeros
		normCCMap.bNoData = bNoData

	# Calculate transition probablity
	calculateTransitionProbablityMatrix(ccMapObj.matrix, normCCMap.bNoData, Out=normCCMap.matrix)
	ma = np.ma.masked_equal(normCCMap.matrix, 0.0, copy=False)
	normCCMap.minvalue = ma.min()
	normCCMap.maxvalue = ma.max()

	normCCMap.make_unreadable()

	# Delete ccmap object if input was a file
	if ccmapType == 'File':
		del ccMapObj

	# Save output ccmap file
	if outFile is not None:
		cmp.save_ccmap(normCCMap, outFile, compress=True)

	# Whether outFile is given.
	if outFile is None:
		return normCCMap
	else:
		del normCCMap
		return None

def transitionProbablityMatrixForGCMap(gcMapInputFile, gcMapOutFile, resolution, compression='lzf', workDir=None, logHandler=None):
	""" To calculate transition matrices using a gcmap file.

	It can be used to calculate transition matrices (markov-chain) for all maps
	present in a gcmap file for given resolution.

	.. note:: Matrices will be calculated for only input resolution. For coarser resolutions, data will be downsampled,
			and therefore in output gcmap, only matrices corresponding to input resolution is correct. Other coarsed
			resolutionsm matrices are only for visualization purpose.

	Parameters
	----------
	gcMapInputFile : str
		Name of input gcmap file.

	gcMapOutFile : str
		Name of output gcmap file.

	resolution : str
		Input resolution at which transition matrix will be calculated.

	compression : str
	    Compression method in output gcmap file. Presently allowed : ``lzf`` for LZF compression and ``gzip`` for GZIP compression.

	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.

	"""


	# Get list of maps in ascending order
	gcmap = gmp.GCMAP(gcMapInputFile)
	gcmap.loadSmallestMap()
	mapList = gcmap.mapNameList.copy()
	del gcmap

	for mapName in mapList:
		logger.info('Calculating transition matrix for {0}...'.format(mapName))

		ccMap = gmp.loadGCMapAsCCMap(gcMapInputFile, mapName=mapName, resolution=resolution, workDir=workDir)

		try:
			trProbMap = transitionProbablityMatrixForCCMap(ccMap, thershold_percentile=percentile_thershold_no_data, thershold_data_occup=thershold_data_occup)

			gmp.addCCMap2GCMap(trProbMap, gcMapOutFile, compression=compression,
								generateCoarse=True, coarsingMethod='sum',
								logHandler=logHandler)

		# In case of program termination, delete the newly created ccmap and raise error
		except (KeyboardInterrupt, SystemExit) as e:
			if 'ccMap' in locals():	del ccMap
			if 'trProbMap' in locals():	del trProbMap
			raise e

		logger.info('       ... Finished Calculating transition matrix.')

		del ccMap
		del trProbMap

def statDistrByEigenDecompForCCMap(ccMap, chrom=None, hdf5Handle=None, compression='lzf'):
	""" Calcualte stationary distribution from a ccmap file or object.

	It uses eigendecomposition method to calculate stationary distribution from
	transition matrix.

	Parameters
	----------
	ccMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ccmap file
		A CCMAP object containing observed contact frequency or a ccmap file

	chrom : str
		Name of chromosome -- neccessary as used in output HDF5 file.

	hdf5Handle : :class:`gcMapExplorer.lib.genomicsDataHandler.HDF5Handler`
		If it is provided, stationary distribution will be directly added to this file.
		In case of ``None``, stationary distribution will be returned as a
		1D array.

	compression : str
	    Compression method in output HDF5 file. Presently allowed : ``lzf``
		for LZF compression and ``gzip`` for GZIP compression.

	Returns
	-------
	sdist : numpy.1darray or ``None``
		If ``hdf5Handle`` is provided ``None`` is returned otherwise stationary distribution as 1D array is returned.

	"""

	if hdf5Handle is not None and chrom is None:
		logger.info('No chromosome name provided. Skiping calculation...')
		return None

	# Check whether input is a file or a obejct
	ccMapObj, ccmapType = cmp.checkCCMapObjectOrFile(ccMap, workDir=workDir)
	ccMapObj.make_readable()

	sdist = stationaryDistributionByEigenDecomp(ccMapObj.matrix, ccMapObj.bNoData)
	ccMapObj.make_unreadable()

	# Remove large outliers
	minvalue = np.amin(sdist)
	mx = np.ma.masked_equal(sdist, minvalue)
	mfiltered = util.detectOutliersMasked1D(mx, thresh=8)
	sdist = mfiltered.filled(0.0)[:]

	# Renormalize
	sdist = sdist/np.sum(sdist)

	resolution = util.binsizeToResolution(ccMapObj.binsize)

	if ccmapType == 'File':
		del ccMapObj

	if hdf5Handle is not None:
		hdf5Handle.addDataByArray(chrom, resolution, 'maximum', sdist, compression=compression)
		hdf5Handle.addDataByArray(chrom, resolution, 'average', sdist, compression=compression)
		hdf5Handle.addDataByArray(chrom, resolution, 'sum', sdist, compression=compression)

		for level in [2, 4, 5, 8, 10, 16, 20, 32]:
			maxs = cmp.downSample1D(sdist, level=level, func='max')
			avgs = cmp.downSample1D(sdist, level=level, func='mean')
			sums = cmp.downSample1D(sdist, level=level, func='sum')
			new_resolution = util.binsizeToResolution( util.resolutionToBinsize(resolution) * level)

			hdf5Handle.addDataByArray(chrom, new_resolution, 'maximum', maxs, compression=compression)
			hdf5Handle.addDataByArray(chrom, new_resolution, 'average', avgs, compression=compression)
			hdf5Handle.addDataByArray(chrom, new_resolution, 'sum', sums, compression=compression)

	else:
		return sdist

def statDistrByEigenDecompForGCMap(gcMapInputFile, outFile, resolution, compression='lzf', workDir=None):
	""" Calculate stationary distribution using transition matrices from gcmap file for given resolution.

	It uses eigendecomposition method to calculate stationary distribution from
	transition matrix.

	.. note:: Use same input resolution as used during calculation of transition matrix using gcmap file.

	Parameters
	----------
	gcMapInputFile : str
		Name of input gcmap file.

	outFile : str
		Name of output HDF5 file to store calculated stationary distribution.

	resolution : str
		Input resolution at which transition matrix was calculated. And
		stationary distribution will be calculated.

	compression : str
	    Compression method in output HDF5 file. Presently allowed : ``lzf`` for LZF compression and ``gzip`` for GZIP compression.

	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.

	"""
	# Get list of maps in ascending order
	gcmap = gmp.GCMAP(gcMapInputFile)
	gcmap.loadSmallestMap()
	mapList = gcmap.mapNameList.copy()
	del gcmap

	hdf5Handle = gmlib.genomicsDataHandler.HDF5Handler(outFile, title='Stationary Distribution')
	hdf5Handle.open()

	for mapName in mapList:

		logger.info('Calculating stationary distribution for {0} ...'.format(mapName))

		ccMap = gmp.loadGCMapAsCCMap(gcMapInputFile, mapName=mapName, resolution=resolution, workDir=workDir)

		try:
			statDistrByEigenDecompForCCMap(ccMap, chrom=mapName, hdf5Handle=hdf5Handle, compression=compression)

		except (KeyboardInterrupt, SystemExit) as e:
			if 'ccMap' in locals():	del ccMap
			if 'hdf5Handle' in locals():	del hdf5Handle
			raise e

		logger.info('             ... Finished Calculating stationary distribution.')
		del ccMap

def stationaryDistributionByEigenDecomp(prob_matrix, bNoData):
	""" To calculate stationary distribution from probablity transition matrix.

	Stationary distribution is calculated using probablity transition matrix
	with eigendecomposition.

	Parameters
	----------
	prob_matrix : numpy.memmap or :attr:`gcMapExplorer.lib.ccmap.CCMAP.matrix` or numpy.ndarray
		Input probablity transition matrix
	bNoData : numpy.array[bool]
		1D-array containing ``True`` and ``False`` values. Its size should be
		equal to input matrix row/column size. Row/Column with ``False`` value
		will be considered during the calculation. Row/Column with ``True``
		value will not be considered during calculation and in these locations,
		minimum stationary distribution will be filled in the output.

	Returns
	-------
	sdist : numpy.ndarray
		Ouput stationary distribution as 1D numpy array.

	"""

	orig_shape = prob_matrix.shape
	matrix = (prob_matrix[~bNoData,:])[:,~bNoData]   # Selected row-column which are not all zeros
	mshape = matrix.shape


	A = MemoryMappedArray(mshape, dtype='float64')
	A.arr[:] = matrix[:]

	# Check if minimum is zero, replace all zeros with minimum value
	ma = np.ma.masked_equal(A.arr, 0.0, copy=False)
	idx = np.nonzero( A.arr == 0 )
	A.arr[idx] = ma.min()

	# Try using scipy sparse package because it is very fast
	#eigvalue, eigvector = sparse_linalg.eigs(A.arr.T, k=1, sigma=1.000)
	#eigvector = np.array(eigvector.flatten(), dtype=np.float64)
	#distrib = eigvector

	# Check if any egienvector value is negative, if yes, re-calculate eigenvector using scipy linear algebra package
	#idx = np.nonzero( distrib < 0 )[0]
	#if idx.any():
	eigvalue, eigvector = 	sp_linalg.eig(A.arr.T)
	distrib = eigvector[:,0]

	distrib = distrib/distrib.sum()

	prob_eigenvector = np.zeros((orig_shape[0],))
	vecmin = np.amin(distrib)

	j = 0
	for i in range(orig_shape[0]):
		if bNoData[i]:
			prob_eigenvector[i] = vecmin
		else:
			prob_eigenvector[i] = distrib[j]
			j = j + 1

	return prob_eigenvector/prob_eigenvector.sum()


def MatrixMultiplyHDD(A, B, Out=None):
	bOut = False
	if type(A) is not MemoryMappedArray:
		raise ValueError ('Not a MemoryMappedArray !!')
	if type(B) is not MemoryMappedArray:
		raise ValueError ('Not a MemoryMappedArray !!')

	if len(A.arr.shape) == 1:
		Am = 1
		An = A.arr.shape[0]
	elif len(A.arr.shape) == 2:
		Am = A.arr.shape[0]
		An = A.arr.shape[1]
	else:
		raise ValueError ('input A must be a vector or square array')

	if len(B.arr.shape) == 1:
		Bm = B.arr.shape[0]
		Bn = 1
	elif len(B.arr.shape) == 2:
		Bm = B.arr.shape[0]
		Bn = B.arr.shape[1]
	else:
		raise ValueError ('input B must be a vector or square array')

	if An != Bm:
		raise ArithmeticError ('Multiplication not possible !!!')

	dtype='float64'

	dtypes=[A.dtype, B.dtype]

	if 'int' in str(A.dtype) and 'int' in str(B.dtype):
		dtype=sorted(dtypes)[-1]

	if 'float' in str(A.dtype) or 'float' in str(B.dtype):
		if 'float' in dtypes:
			dtype = 'float'
		if 'float32' in dtypes:
			dtype = 'float32'
		if 'float64' in dtypes:
			dtype = 'float64'


	if Out is None:
		if Am == 1 and Bn == 1:
			Out = MemoryMappedArray((1, ), dtype=dtype)
		elif Am == 1 and Bn > 1:
			Out = MemoryMappedArray((Bn, ), dtype=dtype)
		elif Am > 1 and Bn == 1:
			Out = MemoryMappedArray((Am, ), dtype=dtype)
		else:
			Out = MemoryMappedArray((Am, Bn), dtype=dtype)
	else:
		bOut = True

	if Am == 1 and Bn == 1:
		Out.arr[0] = np.dot(A.arr, B.arr)
	elif Am == 1 and Bn > 1:
		for j in range(Bn):
			Out.arr[j] = np.dot(A.arr, B.arr[:,j])
	elif Am > 1 and Bn == 1:
		for i in range(Am):
			Out.arr[i] = np.dot(A.arr[i], B.arr)
	else:
		for i in range(Am):
			for j in range(Bn):
				Out.arr[i][j] = np.dot(A.arr[i], B.arr[:,j])

	if bOut:
		return 0
	else:
		if len(Out.arr.shape) == 1 and Out.arr.shape[0] == 1:
			return Out
		else:
			return Out

def MatrixPowerHDD(M, n, result=None):
	if type(M) is not MemoryMappedArray:
		raise ValueError ('Not a MemoryMappedArray !!')

	if len(M.arr.shape) != 2 or M.arr.shape[0] != M.arr.shape[1]:
		raise ValueError("input must be a square array")

	if not np.issubdtype(type(n), int):
		raise TypeError("exponent must be an integer")

	if result is None:
		result = M.copy()
	tmp_result = M.copy()

	if n <= 3:
		for _ in range(n-1):
			MatrixMultiplyHDD(result, M, Out=tmp_result)
			result.copy_from(tmp_result)
		del tmp_result
		return result

    # binary decomposition to reduce the number of Matrix
    # multiplications for n > 3.
	beta = np.binary_repr(n)
	Z = M.copy()
	Zout = M.copy()
	q = 0
	t = len(beta)
	while beta[t-q-1] == '0':
		MatrixMultiplyHDD(Z, Z, Out=Zout)
		Z.copy_from(Zout)
		q += 1

	result.copy_from(Z)
	for k in range(q+1, t):
		MatrixMultiplyHDD(Z, Z, Out=Zout)
		Z.copy_from(Zout)
		if beta[t-k-1] == '1':
			MatrixMultiplyHDD(result, Z, Out=tmp_result)
			result.copy_from(tmp_result)

	del tmp_result
	del Z
	del Zout

	return result

def MatrixPowerRAM(M, n):
	"""Raise a square matrix to the (integer) power `n`.

	.. note:: Ported from numpy/matrixlib/defmatrix.py. This function was rewritten to speed-up the calculation. In place of numpy
			dot function, a mdot function is used, which directly uses blas dgemm function.


	For positive integers `n`, the power is computed by repeated matrix
	squarings and matrix multiplications. If ``n == 0``, the identity matrix
	of the same shape as M is returned. If ``n < 0``, the inverse
	is computed and then raised to the ``abs(n)``.


	Parameters
	----------
	M : ndarray or matrix object
	    Matrix to be "powered."  Must be square, i.e. ``M.shape == (m, m)``,
	    with `m` a positive integer.
	n : int
	    The exponent can be any integer or long integer, positive,
	    negative, or zero.

	Returns
	-------
	M**n : ndarray or matrix object
	    The return value is the same shape and type as `M`;
	    if the exponent is positive or zero then the type of the
	    elements is the same as those of `M`. If the exponent is
	    negative the elements are floating-point.

	Raises
	------
	LinAlgError
	    If the matrix is not numerically invertible.

	"""

	M = np.asanyarray(M)

	if len(M.shape) != 2 or M.shape[0] != M.shape[1]:
		raise ValueError("input must be a square array")
	if not np.issubdtype(type(n), int):
		raise TypeError("exponent must be an integer")

	from numpy.linalg import inv

	if n==0:
		M = M.copy()
		M[:] = np.identity(M.shape[0])
		return M
	elif n<0:
		M = inv(M)
		n *= -1

	result = M
	if n <= 3:
		for _ in range(n-1):
			result = mdot(result, M)
		return result

	# binary decomposition to reduce the number of Matrix
	# multiplications for n > 3.
	beta = np.binary_repr(n)
	Z, q, t = M, 0, len(beta)
	while beta[t-q-1] == '0':
		Z = mdot(Z, Z)
		q += 1
	result = Z
	for k in range(q+1, t):
		Z = mdot(Z, Z)
		if beta[t-k-1] == '1':
			result = mdot(result, Z)
	return result

def CalcStationaryDistributionRAM(orig_matrix, shape, bNonZeros, stop_tol=1e-12, iteration=None):
	matrix = (orig_matrix[bNonZeros,:])[:,bNonZeros]   # Selected row-column which are not all zeros
	mshape = matrix.shape

	A = MemoryMappedArray(mshape, dtype='float64')
	B = MemoryMappedArray(mshape, dtype='float64')
	C = MemoryMappedArray(mshape, dtype='float64')

	A.arr[:] = matrix[:]
	B.arr[:] = matrix[:]
	C.arr[:] = matrix[:]

	if iteration is None:
		iteration = int(shape[0]/2)

	tol = 100
	prev_tol = 100
	c = 1
	while (tol >= stop_tol):
		print('Before Iteration: {0}, Tol: {1}' .format(c, tol))
		tol = LoopStaionaryDistributionRAM(matrix, A, B, C, mshape, c, 60)
		c += 60
		print('After Iteration: {0}, Tol: {1}' .format(c, tol))

		if prev_tol == tol:
			break

		if c >= iteration:
			break

		prev_tol = tol

	vector = np.zeros((shape[0],))
	v_j = 0
	vecmin = np.amin(B.arr)
	for i in range(shape[0]):
		if bNonZeros[i]:
			vector[i] = B.arr[0][v_j]
			v_j = v_j + 1
		else:
			vector[i] = vecmin

	del A
	del B
	del C

	return vector

def CalcStationaryDistributionHDD(orig_matrix, shape, bNonZeros, stop_tol=1e-12, iteration=None):

	matrix = (orig_matrix[bNonZeros,:])[:,bNonZeros]   # Selected row-column which are not all zeros
	mshape = matrix.shape

	MA = MemoryMappedArray(mshape, dtype='float64')
	A = MemoryMappedArray(mshape, dtype='float64')
	B = MemoryMappedArray(mshape, dtype='float64')
	C = MemoryMappedArray(mshape, dtype='float64')

	MA.arr[:] = matrix[:]
	A.arr[:] = matrix[:]
	B.arr[:] = matrix[:]
	C.arr[:] = matrix[:]

	if iteration is None:
		iteration = int(shape[0]/2)

	tol = 100
	prev_tol = 100
	c = 1
	while (tol >= stop_tol):
		print('Before Iteration: {0}, Tol: {1}' .format(c, tol))
		tol = LoopStaionaryDistributionHDD(MA, A, B, C, mshape, c, 10)
		c += 10
		print('After Iteration: {0}, Tol: {1}' .format(c, tol))
		if prev_tol == tol:
			break

		if c >= iteration:
			break

		prev_tol = tol

	vector = np.zeros((shape[0],))
	v_j = 0
	vecmin = np.amin(B.arr[0])
	for i in range(shape[0]):
		if bNonZeros[i]:
			vector[i] = B.arr[0][v_j]
			v_j = v_j + 1
		else:
			vector[i] = vecmin

	del A
	del B
	del C
	del MA

	return vector

def LoopStaionaryDistributionRAM(matrix, A, B, C, mshape, c, nc):
	if c==1:
		A.arr[:] = MatrixPowerRAM(matrix, nc-1)[:]
	else:
		C.arr[:] = MatrixPowerRAM(matrix, nc-1)[:]
		A.arr[:] = mdot(B.arr, C.arr)
	B.arr[:] = mdot(A.arr, matrix)[:]

	#B.arr[:] = np.linalg.matrix_power(matrix, nc)[:]
	#A.arr[:] = np.dot(B.arr, matrix)[:]
	tidx = np.random.randint(mshape[0]-1, size=1)
	tol = np.amax((B.arr[tidx] - A.arr[tidx])**2)

	return tol

def LoopStaionaryDistributionHDD(MA, A, B, C, mshape, c, nc):
	if c==1:
		MatrixPowerHDD(MA, nc-1, result=A)
	else:
		MatrixPowerHDD(MA, nc-1, result=C)
		MatrixMultiplyHDD(B, C, Out=A)
	MatrixMultiplyHDD(A, MA, Out=B)

	tidx = np.random.randint(mshape[0]-1, size=1)
	tol = np.amax((B.arr[tidx] - A.arr[tidx])**2)

	return tol


def StationaryDistributionByMultiply(ccMapObj, stop_tol=1e-12, vector=None, iteration=None):

	ccMapObj.make_readable()
	bNonZeros = ~ccMapObj.bNoData
	mem = psutil.virtual_memory()

	if ccMapObj.matrix.nbytes < mem.available:
		print("Using RAM")
		vector = CalcStationaryDistributionRAM(ccMapObj.matrix, ccMapObj.shape, bNonZeros, stop_tol=stop_tol, iteration=None)
	else:
		print("Using HDD")
		vector = CalcStationaryDistributionHDD(ccMapObj.matrix, ccMapObj.shape, bNonZeros, stop_tol=stop_tol, iteration=None)

	return vector
