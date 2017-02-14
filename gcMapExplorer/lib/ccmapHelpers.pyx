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
import scipy.sparse as sparse
# import scipy.sparse.linalg as sparse_linalg
from scipy.linalg import get_blas_funcs
from scipy import linalg as sp_linalg
import scipy.stats as stats
import string, random
import os, re

from scipy.stats import mstats
import numpy.ma as ma

dtype_npBINarray = 'float32'


def sorted_nicely( l ):
	""" Sorts the given iterable in the way that is expected.

	Required arguments:
	l -- The iterable to be sorted.

	"""
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	return sorted(l, key = alphanum_key)

def locate_significant_digit_after_decimal(value):
	""" Get location at which significant digit start after decimal
	"""
	location = 0
	base = 0
	while( base <=0 ):
		base = int( value * 10 ** location  )
		location = location + 1

	return location

def kth_diag_indices(k, a):
	""" Get diagonal indices of array 'a' offset by 'k'
	"""
	rows, cols = np.diag_indices_from(a)
	if k < 0:
		return rows[:k], cols[-k:]
	elif k > 0:
		return rows[k:], cols[:-k]
	else:
		return rows, cols

# Calling blas dot function for matrix
def mdot(v, w):
	gemm = get_blas_funcs("gemm", [v, w])
	return gemm(alpha=1., a=v, b=w)

class MemoryMappedArray:
	"""Convenient wrapper for numpy memory mapped array file

	For more details, see here: (See: |numpy memmap|).

	Attributes
	----------
	path2matrix : str
		Path to numpy memory mapped array file
	arr : numpy.memmap
		Pointer to memory mapped numpy array
	workDir : str
		Path to the directory where temporary intermediate files are generated. If ``None``, files are generated in the temporary directory according to the OS type.
	dtype : str
		Data type of array

	Parameters
	----------
	shape : tuple
		Shape of array
	fill : int or float (Optional)
		Fill array with this value
	dtype : str
		Data type of array

	"""
	def __init__(self, shape, fill=None, workDir=None, dtype='float32'):
		self.workDir = workDir

		# Generating a file with random name
		(fd, fname) = tempfile.mkstemp(suffix='.tmp', prefix='nparray_', dir=self.workDir, text=False)
		os.close(fd)     # Close file, error in windows OS

		# Generating a memmap array
		self.path2matrix = fname
		self.arr = np.memmap(self.path2matrix, dtype=dtype, mode='w+', shape=shape)
		self.dtype = dtype
		if fill is not None:
			self.arr.fill(fill)

	def copy(self):
		"""Copy this numpy memory mapped array and generate new

		Returns
		-------
		out : :class:`MemoryMappedArray`
			A new :class:`MemoryMappedArray` instance with copied arrays

		"""
		out = MemoryMappedArray(self.arr.shape, workDir=self.workdir, dtype=self.dtype)
		out.arr[:] = self.arr[:]
		return out

	def copy_from(self, src):
		"""Copy values from source :class:`MemoryMappedArray`

		Parameters
		----------
		src : :class:`MemoryMappedArray`
			Source memory mapped arrays for new values

		Returns
		-------
		None


		Raises
		------
		ValueError
			if src is not of :class:`MemoryMappedArray` instance

		"""
		if not isinstance(src, MemoryMappedArray):
			raise ValueError ('Not a MemoryMappedArray !!')

		self.arr[:] = src.arr[:]
		self.arr.flush()

	def copy_to(self, dest):
		"""Copy values to destination :class:`MemoryMappedArray`

		Parameters
		----------
		dest : :class:`MemoryMappedArray`
			Destination memory mapped arrays

		Returns
		-------
		None


		Raises
		------
		ValueError
			if dest is not of :class:`MemoryMappedArray` instance

		"""

		if type(dest) is not MemoryMappedArray:
			raise ValueError ('Not a MemoryMappedArray !!')
		dest.arr[:] = self.arr[:]
		#dest.arr.flush()

	def __del__(self):
		del self.arr
		try:
			os.remove(self.path2matrix)
		except:
			pass

def get_nonzeros_index(matrix, thershold_percentile=None, thershold_data_occup=None):
	"""To get a numpy array of bool values for all rows/columns which have **NO** missing data

	Parameters
	----------
	matrix : numpy.memmap or :attr:`gcMapExplorer.lib.ccmap.CCMAP.matrix`
		Input matrix
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
	bData : numpy.array[bool]
		1D-array containing ``True`` and ``False`` values.
		* If ``True``: row/column has data above the thershold
		* If ``False``: row/column has no data under the thershold
	"""
	bData = np.empty(shape=matrix.shape[0], dtype='bool')
	mx = matrix.shape[0]
	#my = matrix.shape[1]

	zero_count = []
	for i in range(mx):
		if np.sum( matrix[i] ) == 0.0:
			bData[i] = False
			# Although whole column/row is zero, count of zero taken as zero
			zero_count.append(0)
		else:
			bData[i] = True
			zero_count.append(np.nonzero( matrix[i] == 0 )[0].shape[0])


	zero_count = np.asarray(zero_count)

	if thershold_percentile is not None and thershold_data_occup is not None:
		raise AssertionError("Both 'thershold_percentile' and 'thershold_count_ratio' cannot be used simultaneously!")

	# Make false if number of zeros is less than thershold
	if thershold_percentile is not None:
		percentile = np.percentile(zero_count[np.nonzero(zero_count != 0)], thershold_percentile)
		for i in range(mx):
			if bData[i]:
				if zero_count[i] >= percentile:
					bData[i] = False
					#print(i, zero_count[i], percentile)

	# Make false if number of zeros is less than thershold
	if thershold_data_occup is not None:
		for i in range(mx):
			if bData[i]:
				if zero_count[i]/mx >=  (1.0 - thershold_data_occup):
					bData[i] = False
					#print(i, count[i], real_data, thershold_percentile)

	bData = np.asarray(bData)

	length = np.nonzero( bData == True )[0].shape[0]
	if length == 0:
		if thershold_percentile is not None:
			raise ValueError("All rows contain zero values!!! Try increasing 'thershold_percentile' value or do not use it.")
		if thershold_data_occup is not None:
			raise ValueError("All rows contain zero values!!! Try decreasing 'thershold_data_occup' value or do not use it.")

	return bData

def remove_zeros(matrix, thershold_percentile=None, thershold_data_occup=None, workDir=None):
	"""To remove rows/columns with missing data (zero values)

	Parameters
	----------
	matrix : numpy.memmap or :attr:`gcMapExplorer.lib.ccmap.CCMAP.matrix`
		Input matrix
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
		Path to the directory where temporary intermediate files are generated. If ``None``, files are generated in the temporary directory according to the OS type.


	Returns
	-------
	A : :class:`MemoryMappedArray`
		:class:`MemoryMappedArray` instance containing new truncated array as memory mapped file
	bNoData : numpy.array[bool]
		1D-array containing ``True`` and ``False`` values.
		* If ``True``: row/column has no data under the thershold
		* If ``False``: row/column has data above the thershold

	"""

	bData = np.empty(shape=matrix.shape[0], dtype='bool')
	mx = matrix.shape[0]
	my = matrix.shape[1]

	bData = get_nonzeros_index(matrix, thershold_percentile=thershold_percentile, thershold_data_occup=thershold_data_occup)

	length = np.nonzero( bData == True )[0].shape[0]
	A = MemoryMappedArray((length, length),  workDir=workDir, dtype='float64')

	ci = -1
	cj = 0
	for i in range(mx):
		if not bData[i]:
			continue
		else:
			ci += 1
		cj = 0
		for j in range(my):
			if bData[j]:
				A.arr[ci][cj] = matrix[i][j]
				cj += 1

	return A, ~bData

def genFullMapFromTruncMap(inMatrix, outMatrix, bNoData, minvalue=None):

	if minvalue is None:
		minvalue = np.amin(inMatrix)

	inI = -1
	inJ = 0
	ox = outMatrix.shape[0]
	oy = outMatrix.shape[1]

	for outI in range(ox):
		if not bNoData[outI]:
			inI += 1

		inJ = 0
		for outJ in range(oy):
			if bNoData[outI] or bNoData[outJ]:
				outMatrix[outI][outJ] = minvalue
			else:
				if inMatrix[inI][inJ] == 0:
					outMatrix[outI][outJ] = minvalue
				else:
					outMatrix[outI][outJ] = inMatrix[inI][inJ]
				inJ += 1



def CalculateTransitionProbablity(A, bNoData, Out):
	mshape = A.shape

	minval = 9e10
	maxval = -9e10
	for m in range(mshape[0]):
		if bNoData[m]:
			continue
		sumr = A[m].sum()
		Out[m] = (A[m] / sumr)[:]

		x = np.sort(Out[m,np.nonzero(Out[m])])[0]
		minval = min(minval, x[0])
		maxval = max(maxval, x[-1])
		x=None

	for m in range(mshape[0]):
		if bNoData[m]:
			Out[m].fill(minval)
			Out[:,m].fill(minval)
		else:
			idx = np.nonzero( Out[m] == 0 )
			Out[m][idx] = minval

	return float(minval), float(maxval)

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

	.. note::	Ported from numpy/matrixlib/defmatrix.py. This function was rewritten to speed-up the calculation. In pace of numpy
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
	See Also
	--------
	matrix
	    Provides an equivalent function as the exponentiation operator
	    (``**``, not ``^``).
	Examples
	--------
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
	""" Calculate stationary distibution
	Hellloooo

	"""
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

def CalcStationaryDistributionByEigen(orig_matrix, shape, bNoData):
	matrix = (orig_matrix[~bNoData,:])[:,~bNoData]   # Selected row-column which are not all zeros
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

	prob_eigenvector = np.zeros((shape[0],))
	vecmin = np.amin(distrib)

	j = 0
	for i in range(shape[0]):
		if bNoData[i]:
			prob_eigenvector[i] = vecmin
		else:
			prob_eigenvector[i] = distrib[j]
			j = j + 1

	return prob_eigenvector/prob_eigenvector.sum()


def make_smooth_map(matrix, shape, filter_pass=1):

	def calc_mean(A, B, i, j, symmetric):
		istart = i-1
		iend = i+2
		jstart = j-1
		jend = j+2

		if istart < 0:
			istart = 1
		if jstart < 0:
			jstart = 1
		if iend > shape[0]:
			iend = shape[0]
		if jend > shape[1]:
			jend = shape[1]

		value = A.arr[istart:iend, jstart:jend].mean()
		B.arr[i][j] = value
		if symmetric:
			B.arr[j][i] = value


	A = MemoryMappedArray(shape, dtype='float32')
	A.arr[:] = matrix[:]
	B = MemoryMappedArray(shape, dtype='float32')

	symmetric = False
	if np.allclose(matrix.T, matrix):
		symmetric = True

	pass_num = 1
	while(pass_num <= filter_pass):
		for i in range(shape[0]):
			if symmetric:
				for j in range(i+1):
					calc_mean(A, B, i, j, symmetric)
			else:
				for j in range(shape[1]):
					calc_mean(A, B, i, j, symmetric)

			A.arr.flush()
			B.arr.flush()

		pass_num += 1
		A.arr[:] = B.arr[:]

	matrix[:] = A.arr[:]
	del A
	del B

def get_min_max_values(matrix, get_bNoData=False):
	bData = np.empty(shape=matrix.shape[0], dtype='bool')

	count = 0
	mx = matrix.shape[0]
	my = matrix.shape[1]

	for i in range(mx):
		if np.sum(matrix[i]) == 0.0:
			bData[i] = False
		else:
			count += 1
			bData[i] = True

	bData = np.asarray(bData)

	ci = -1
	cj = 0
	maxvalue = -9e10
	minvalue = 9e10
	for i in range(mx):
		if not bData[i]:
			continue
		for j in range(my):
			s = np.sort(matrix[i, np.nonzero(matrix[i])])[0]
			minvalue = min(minvalue, s[0])
			maxvalue = max(maxvalue, s[-1])

	if get_bNoData:
		return minvalue, maxvalue, ~bData
	else:
		return minvalue, maxvalue

def name_generator(size=10, chars=string.ascii_letters + string.digits ):
	s = ''.join(random.choice(chars) for _ in range(size))
	return s
