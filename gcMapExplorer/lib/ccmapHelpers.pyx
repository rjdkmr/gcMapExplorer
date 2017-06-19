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

from gcMapExplorer.config import getConfig

# get configuration
config = getConfig()

dtype_npBINarray = 'float32'

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

		# Working and output directory
		if workDir is None:
			self.workDir = config['Dirs']['WorkingDirectory']
		else:
			self.workDir = workDir

		# Generating a file with random name
		(fd, fname) = tempfile.mkstemp(suffix='.tmp', prefix='gcx_nparray_', dir=self.workDir, text=False)
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
		if os.path.isfile(self.path2matrix):
			os.remove(self.path2matrix)

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
