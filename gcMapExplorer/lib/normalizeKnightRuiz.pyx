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
import time

import string, random
import os

from . import ccmapHelpers as cmh

dtype_npBINarray = 'float32'


class KnightRuizNorm:
	"""A modified Knight-Ruiz algorithm for matrix balancing

	The original ported Knight-Ruiz algorithm is modifed to implement the normalization using both memory/RAM and disk.
	It allows the normalization of small Hi-C maps to huge maps that could not be accomodated in RAM.

	Parameters
	----------
	A : numpy.ndarray or :class:`MemoryMappedArray`
		Input matrix.

		.. note::

			* Matrix should not contain any row or column with all zero values (missing data for row/column). This type of matrix can be obtained from :meth:`remove_zeros`.
			* If ``memory='HDD'``, ``A`` should be :class:`MemoryMappedArray`

	memory : str
		Accepted keywords are ``RAM`` and ``HDD``:

			* ``RAM``: All intermediate arrays are generated in memory(RAM). This version is faster, however, it requires RAM depending on the input matrix size.
			* ``HDD``: All intermediate arrays are generated as memory mapped array files on hard-disk.

	workDir : str
		Path to the directory where temporary intermediate files are generated. If ``None``, files are generated in the temporary directory according to the OS type.


	"""
	def __init__(self, A, memory='RAM', tol=1e-12, delta=0.1, Delta=3, workDir=None):

		self.workDir = workDir

		if memory == 'RAM':
			self.init_for_RAM(A, tol=tol, delta=delta, Delta=Delta)
		elif memory == 'HDD':
			self.init_for_HDD(A, tol=tol, delta=delta, Delta=Delta)
		else:
			raise ValueError ('This memory={0} is not is not understandable.. Please use \'RAM\' or \'HDD\'.' .format(memory))
		self.memory = memory

	def init_for_RAM(self, A, tol=1e-12, delta=0.1, Delta=3):
		# x, x0, e, v, rk, y, Z, w, p, ap :     vector shape(n, 1) : [ [value] [value] [value] [value] ... ... ... [value] ]
		# rho_km1, rout, rold, innertol, alpha :  scalar shape(1 ,1) : [[value]]

		# Dummy Initialazation of x, x0, e, v, rk, y, Z, w, p, ap

		self.x = None
		self.x0 = None
		self.v = None
		self.w = None
		self.rk = None
		self.y = None
		self.ynew = None
		self.Z = None
		self.p = None
		self.ap = None

		# Dummy initialization rho_km1, rout, rold, innertol, alpha, k
		self.rho_km1 = None
		self.rout = None
		self.rold = None
		self.innertol = None
		self.alpha = None
		self.k = None

		# Top level initialization
		self.tol = tol
		self.delta = delta
		self.Delta = Delta
		self.MVP = 0                                  # Well count matrix vector products.
		self.i = 0                                      # Outer iteration count.

		# Real initialization
		self.n = A.shape[0]                                                                 # n = size(A,1)
		self.e = np.ones((self.n, 1))                                                            # e = ones(n,1)
		self.res = []
		self.x0 = self.e


		self.g = 0.9        # Parameters used in inner stopping criterion.
		self.etamax = 0.1   # Parameters used in inner stopping criterion.
		self.eta = self.etamax
		self.stop_tol = self.tol * 0.5

		self.x = self.x0
		self.rt = self.tol**2                                 # rt = tol^2
		self.v = self.x * np.dot(A, self.x)                            # v = x.*(A*x)
		self.rk = 1 - self.v
		self.rho_km1 = np.dot( self.rk.conjugate().T, self.rk)     # rho_km1 = rk'*rk
		self.rout = self.rho_km1
		self.rold = self.rout

	def doInnerLoopRAM(self, A):
		self.k = self.k + 1

		if self.k == 1:
			self.Z = self.rk / self.v                                       # Z = rk./v
			self.p = self.Z
			self.rho_km1 = np.dot( self.rk.conjugate().T, self.Z)           # rho_km1 = rk'*Z
		else:
			beta = self.rho_km1 / self.rho_km2
			self.p = self.Z + (beta * self.p)


		# Update search direction efficiently.
		self.w = self.x * np.dot(A, (self.x * self.p)) + (self.v * self.p)                       # w = x.*(A*(x.*p)) + v.*p
		self.alpha = self.rho_km1 / np.dot( self.p.conjugate().T, self.w)        # alpha = rho_km1/(p'*w)
		self.ap = self.alpha * self.p                                       # ap = alpha*p (No dot function as alpha is scalar)

		# Test distance to boundary of cone.
		self.ynew = self.y + self.ap;
		#print(i, np.amin(ynew), delta, np.amin(ynew) <= delta)
		#print(i, np.amax(ynew), Delta, np.amax(ynew) >= Delta)
		if np.amin(self.ynew) <= self.delta:
			if self.delta == 0:
				return 'break'
			ind = np.nonzero(self.ap < 0) # ind = find(ap < 0)
			gamma = np.amin( (self.delta - self.y[ind]) / self.ap[ind] )    # gamma = min((delta - y(ind))./ap(ind))
			self.y = self.y + np.dot(gamma, self.ap)                        # y = y + gamma*ap
			return 'break'
		if np.amax(self.ynew) >= self.Delta:
			ind = np.nonzero( self.ynew > self.Delta )                 # ind = find(ynew > Delta);
			gamma = np.amin( (self.Delta - self.y[ind]) / self.ap[ind])       # gamma = min((Delta-y(ind))./ap(ind));
			self.y = self.y + np.dot(gamma, self.ap)                        # y = y + gamma*ap;
			return 'break'
		self.y = self.ynew
		self.rk = self.rk - self.alpha * self.w                                    # rk = rk - alpha*w
		self.rho_km2 = self.rho_km1
		self.Z = self.rk / self.v
		self.rho_km1 = np.dot( self.rk.conjugate().T, self.Z)               # rho_km1 = rk'*Z

		return 'complete'

	def doOuterLoopRAM(self, A, fl):
		self.x = self.x * self.y                                                # x = x.*y
		self.v = self.x * np.dot(A, self.x)                                         # v = x.*(A*x)
		self.rk = 1 - self.v
		self.rho_km1 = np.dot( self.rk.conjugate().T, self.rk)                  # rho_km1 = rk'*rk
		self.rout = self.rho_km1
		self.MVP =self.MVP + self.k + 1

		# Update inner iteration stopping criterion.
		rat = self.rout/self.rold
		self.rold = self.rout
		res_norm = np.sqrt(self.rout)
		eta_o = self.eta
		self.eta = self.g*rat

		#print(i, res_norm)

		if self.g*eta_o**2 > 0.1:
			self.eta = np.amax([self.eta, self.g*eta_o**2])                    # eta = max([eta,g*eta_o^2])

		self.eta = np.amax([np.amin([self.eta, self.etamax]), self.stop_tol/res_norm]);   # eta = max([min([eta,etamax]),stop_tol/res_norm]);

		if fl == 1:
			print('%3d %6d %.3e %.3e %.3e \n' % (self.i, self.k, res_norm, np.amin(self.y), np.amin(self.x)))
			self.res=[self.res, res_norm]

	def GenNormMatrixRAM(self, A, OutMatrix, bNoData):
		# Generation of Doubly stochastic matrix ( diag(X)*A*diag(X) )
		(fd, path2matrix) = tempfile.mkstemp(suffix='.bin', prefix='nparray_', dir=self.workDir, text=False)
		os.close(fd)     # Close file, error in windows OS

		A_DSMat = np.memmap(path2matrix, dtype=dtype_npBINarray, mode='w+', shape=A.shape)

		try:
			# To perform (self.x.T * (A * self.x)) in steps to remove load on RAM
			minvalue = 1e10
			maxvalue = -1e10
			mx = self.x.shape[0]
			x_d = cmh.MemoryMappedArray((mx, mx), workDir=self.workDir, dtype='float64')

			for i in range(mx):
				x_d.arr[i] = self.x[:,0] * A[i]
			x_d.arr.flush()

			for i in range(mx):
				A_DSMat[i] = x_d.arr[i] * self.x[i,:]

			ma = np.ma.masked_equal(A_DSMat, 0.0, copy=False)
			minvalue = ma.min()
			maxvalue = ma.max()

			del x_d

			dsm_i = 0
			ox = OutMatrix.shape[0]
			idx_fill = np.nonzero( ~bNoData )
			for i in range(ox):
				if not bNoData[i]:
					OutMatrix[i, idx_fill] = A_DSMat[dsm_i]
					OutMatrix[idx_fill, i] = A_DSMat[dsm_i]
					dsm_i += 1

			''' Old One, extremely slow, kept here for reference
			dsm_i = -1
			dsm_j = 0
			ox = OutMatrix.shape[0]
			oy = OutMatrix.shape[1]

			time1 = time.time()
			for i in range(ox):
				if not bNoData[i]:
					dsm_i += 1

				dsm_j = 0
				for j in range(oy):
					if bNoData[i] or bNoData[j]:
						OutMatrix[i][j] = 0.0
					else:
						OutMatrix[i][j] = A_DSMat[dsm_i][dsm_j]
						dsm_j += 1

			time2 = time.time()
			print('Re-assign:', time2-time1)
			'''

			#r_sum = A_DSMat.sum(axis = 0)
			#c_sum = A_DSMat.sum(axis = 1)
			#for i in range(A_DSMat.shape[0]):
				#print(i, self.x[i], r_sum[i], c_sum[i])

			OutMatrix.flush()

		except (KeyboardInterrupt, SystemExit) as e:
			del A_DSMat
			os.remove(path2matrix)
			raise e

		# Removing temporary memory-mapped array file
		del A_DSMat
		if os.path.isfile(path2matrix):
			os.remove(path2matrix)

		return minvalue, maxvalue

	def run_for_RAM(self, A, fl, OutMatrix, bNoData):
		if fl == 1:
			print('it in. it res')

		while self.rout > self.rt: 										     # Outer iteration
			self.i = self.i + 1
			self.k = 0
			self.y = self.e
			self.innertol = max( [ self.eta**2 * self.rout, self.rt ])                   # innertol = max([eta^2*rout,rt]);

			while self.rho_km1 > self.innertol:                                # Inner iteration by CG
				if self.doInnerLoopRAM(A) == 'break':
					break

			self.doOuterLoopRAM(A, fl)

		minvalue, maxvalue = self.GenNormMatrixRAM(A, OutMatrix, bNoData)
		return minvalue, maxvalue

	def init_for_HDD(self, A, tol=1e-12, delta=0.1, Delta=3):
		# x, x0, e, v, rk, y, Z, w, p, ap :     vector shape(n, 1) : [ [value] [value] [value] [value] ... ... ... [value] ]
		# rho_km1, rout, rold, innertol, alpha :  scalar shape(1 ,1) : [[value]]

		self.n = A.shape[0]                                                             # n = size(A,1)
		self.e = cmh.MemoryMappedArray((self.n, 1), fill=1.0,  workDir=self.workDir,  dtype='float64')                            # e = ones(n,1)

		# Dummy Initialazation of x, x0, e, v, rk, y, Z, w, p, ap
		self.x = cmh.MemoryMappedArray((self.n, 1), workDir=self.workDir,  dtype='float64')
		self.x0 = cmh.MemoryMappedArray((self.n, 1),  workDir=self.workDir, dtype='float64')
		self.v = cmh.MemoryMappedArray((self.n, 1),  workDir=self.workDir, dtype='float64')
		self.w = cmh.MemoryMappedArray((self.n, 1),  workDir=self.workDir, dtype='float64')
		self.rk = cmh.MemoryMappedArray((self.n, 1),  workDir=self.workDir, dtype='float64')
		self.y = cmh.MemoryMappedArray((self.n, 1),  workDir=self.workDir, dtype='float64')
		self.ynew = cmh.MemoryMappedArray((self.n, 1),  workDir=self.workDir, dtype='float64')
		self.Z = cmh.MemoryMappedArray((self.n, 1),  workDir=self.workDir, dtype='float64')
		self.p = cmh.MemoryMappedArray((self.n, 1), workDir=self.workDir,  dtype='float64')
		self.ap = cmh.MemoryMappedArray((self.n, 1),  workDir=self.workDir, dtype='float64')


		# Dummy initialization rho_km1, rout, rold, innertol, alpha, k
		self.rho_km1 = None
		self.rout = None
		self.rold = None
		self.innertol = None
		self.alpha = None
		self.k = None

		# Top level initialization
		self.tol = tol
		self.delta = delta
		self.Delta = Delta
		self.MVP = 0                                  # Well count matrix vector products.
		self.i = 0                                      # Outer iteration count.

		# Real initialization
		self.res = []
		self.x0.copy_from(self.e)


		self.g = 0.9        # Parameters used in inner stopping criterion.
		self.etamax = 0.1   # Parameters used in inner stopping criterion.
		self.eta = self.etamax
		self.stop_tol = self.tol * 0.5

		self.x.copy_from(self.x0)
		self.rt = self.tol**2                                 # rt = tol^2
		self.v.arr[:] = (self.x.arr * np.dot(A, self.x.arr))[:]                            # v = x.*(A*x)
		self.rk.arr[:] = (1 - self.v.arr)[:]
		self.rho_km1 = np.dot( self.rk.arr.conjugate().T, self.rk.arr)     # rho_km1 = rk'*rk
		self.rout = self.rho_km1
		self.rold = self.rout

	def doInnerLoopHDD(self, A):

		self.k = self.k + 1

		if self.k == 1:
			self.Z.arr[:] = self.rk.arr / self.v.arr                                       # Z = rk./v
			self.p.copy_from(self.Z)
			self.rho_km1 = np.dot( self.rk.arr.conjugate().T, self.Z.arr)           # rho_km1 = rk'*Z
		else:
			beta = self.rho_km1 / self.rho_km2
			self.p.arr[:] = (self.Z.arr + (beta * self.p.arr))[:]


		# Update search direction efficiently.
		self.w.arr[:] = self.x.arr * np.dot(A, (self.x.arr * self.p.arr)) + (self.v.arr * self.p.arr)                       # w = x.*(A*(x.*p)) + v.*p
		self.alpha = self.rho_km1 / np.dot( self.p.arr.conjugate().T, self.w.arr)        # alpha = rho_km1/(p'*w)
		self.ap.arr[:] = self.alpha * self.p.arr                                       # ap = alpha*p (No dot function as alpha is scalar)

		# Test distance to boundary of cone.
		self.ynew.arr[:] = self.y.arr + self.ap.arr
		#print(i, np.amin(ynew), delta, np.amin(ynew) <= delta)
		#print(i, np.amax(ynew), Delta, np.amax(ynew) >= Delta)
		if np.amin(self.ynew.arr) <= self.delta:
			if self.delta == 0:
				return 'break'
			ind = np.nonzero(self.ap.arr < 0) # ind = find(ap < 0)
			gamma = np.amin( (self.delta - self.y.arr[ind]) / self.ap.arr[ind] )    # gamma = min((delta - y(ind))./ap(ind))
			self.y.arr[:] = self.y.arr + np.dot(gamma, self.ap.arr)                        # y = y + gamma*ap
			return 'break'
		if np.amax(self.ynew.arr) >= self.Delta:
			ind = np.nonzero( self.ynew.arr > self.Delta )                 # ind = find(ynew > Delta);
			gamma = np.amin( (self.Delta - self.y.arr[ind]) / self.ap.arr[ind])       # gamma = min((Delta-y(ind))./ap(ind));
			self.y.arr[:] = self.y.arr + np.dot(gamma, self.ap.arr)                        # y = y + gamma*ap;
			return 'break'
		self.y.copy_from(self.ynew)
		self.rk.arr[:] = self.rk.arr - self.alpha * self.w.arr                                    # rk = rk - alpha*w
		self.rho_km2 = self.rho_km1
		self.Z.arr[:] = self.rk.arr / self.v.arr
		self.rho_km1 = np.dot( self.rk.arr.conjugate().T, self.Z.arr)               # rho_km1 = rk'*Z

		return 'complete'

	def doOuterLoopHDD(self, A, fl):
		self.x.arr[:] = self.x.arr * self.y.arr                                                # x = x.*y
		self.v.arr[:] = self.x.arr * np.dot(A, self.x.arr)                                         # v = x.*(A*x)
		self.rk.arr[:] = (1 - self.v.arr)[:]
		self.rho_km1 = np.dot( self.rk.arr.conjugate().T, self.rk.arr)                  # rho_km1 = rk'*rk
		self.rout = self.rho_km1
		self.MVP =self.MVP + self.k + 1

		# Update inner iteration stopping criterion.
		rat = self.rout/self.rold
		self.rold = self.rout
		res_norm = np.sqrt(self.rout)
		eta_o = self.eta
		self.eta = self.g*rat

		#print(i, res_norm)

		if self.g*eta_o**2 > 0.1:
			self.eta = np.amax([self.eta, self.g*eta_o**2])                    # eta = max([eta,g*eta_o^2])

		self.eta = np.amax([np.amin([self.eta, self.etamax]), self.stop_tol/res_norm]);   # eta = max([min([eta,etamax]),stop_tol/res_norm]);

		if fl == 1:
			print('%3d %6d %.3e %.3e %.3e \n' % (self.i, self.k, res_norm, np.amin(self.y.arr), np.amin(self.x.arr)))
			self.res=[self.res, res_norm]

	def GenNormMatrixHDD(self, A, OutMatrix, bNoData):
		# Generation of Doubly stochastic matrix ( diag(X)*A*diag(X) )
		(fd, path2matrix) = tempfile.mkstemp(suffix='.bin', prefix='nparray_', dir=self.workDir, text=False)
		os.close(fd)     # Close file, error in windows OS
		A_DSMat = np.memmap(path2matrix, dtype=dtype_npBINarray, mode='w+', shape=A.shape)

		try:
			minvalue = 1e10
			maxvalue = -1e10
			mx = self.x.arr.shape[0]
			x_d = cmh.MemoryMappedArray((mx, mx),  workDir=self.workDir, dtype='float64')
			for i in range(mx):
				x_d.arr[i] = self.x.arr[:,0] * A[i]
			x_d.arr.flush()

			for i in range(mx):
				A_DSMat[i] = x_d.arr[i] * self.x.arr[i,:]

			ma = np.ma.masked_equal(A_DSMat, 0.0, copy=False)
			minvalue = ma.min()
			maxvalue = ma.max()

			del x_d

			dsm_i = 0
			ox = OutMatrix.shape[0]
			idx_fill = np.nonzero( ~bNoData )
			for i in range(ox):
				if not bNoData[i]:
					OutMatrix[i, idx_fill] = A_DSMat[dsm_i]
					OutMatrix[idx_fill, i] = A_DSMat[dsm_i]
					dsm_i += 1

			''' OLD ONE
			dsm_i = -1
			dsm_j = 0

			ox = OutMatrix.shape[0]
			oy = OutMatrix.shape[1]

			for i in range(ox):
				if not bNoData[i]:
					dsm_i += 1

				dsm_j = 0
				for j in range(oy):
					if bNoData[i] or bNoData[j]:
						OutMatrix[i][j] = 0.0
					else:
						OutMatrix[i][j] = A_DSMat[dsm_i][dsm_j]
						dsm_j += 1
			'''

			#r_sum = A_DSMat.sum(axis = 0)
			#c_sum = A_DSMat.sum(axis = 1)
			#for i in range(A_DSMat.shape[0]):
				#print(i, x[i], r_sum[i], c_sum[i])

			OutMatrix.flush()

		except (KeyboardInterrupt, SystemExit) as e:
			del A_DSMat
			os.remove(path2matrix)
			raise e


		del A_DSMat
		if os.path.isfile(path2matrix):
			os.remove(path2matrix)

		return minvalue, maxvalue

	def del_intermediate_files(self):
		if isinstance(self.e, cmh.MemoryMappedArray):
			del self.e
		if isinstance(self.x, cmh.MemoryMappedArray):
			del self.x
		if isinstance(self.x0, cmh.MemoryMappedArray):
			del self.x0
		if isinstance(self.v, cmh.MemoryMappedArray):
			del self.v
		if isinstance(self.w, cmh.MemoryMappedArray):
			del self.w
		if isinstance(self.rk, cmh.MemoryMappedArray):
			del self.rk
		if isinstance(self.y, cmh.MemoryMappedArray):
			del self.y
		if isinstance(self.ynew, cmh.MemoryMappedArray):
			del self.ynew
		if isinstance(self.Z, cmh.MemoryMappedArray):
			del self.Z
		if isinstance(self.p, cmh.MemoryMappedArray):
			del self.p
		if isinstance(self.ap, cmh.MemoryMappedArray):
			del self.ap

	def run_for_HDD(self, A, fl, OutMatrix, bNoData):
		if fl == 1:
			print('it in. it res')

		while self.rout > self.rt: 										     # Outer iteration
			self.i = self.i + 1
			self.k = 0
			self.y.copy_from(self.e)
			self.innertol = np.amax( [ self.eta**2 * self.rout, self.rt ])                   # innertol = max([eta^2*rout,rt]);

			while self.rho_km1 > self.innertol:                                # Inner iteration by CG
				if self.doInnerLoopHDD(A) == 'break':
					break

			self.doOuterLoopHDD(A, fl)

		minvalue, maxvalue = self.GenNormMatrixHDD(A, OutMatrix, bNoData)

		return minvalue, maxvalue

	def run(self, A, fl, OutMatrix, bNoData):
		"""Perform Knight-Ruiz normalization

		Parameters
		----------
		A : numpy.ndarray or :attr:`MemoryMappedArray.arr`
			Input matrix.

			.. note::
				* Matrix should not contain any row or column with all zero values (missing data for row/column). This type of matrix can be obtained from :meth:`remove_zeros`.

			.. warning::
				If A was :class:`MemoryMappedArray` in :class:`KnightRuizNorm`. Here ``A`` should be :attr:`MemoryMappedArray.arr` instead of :class:`MemoryMappedArray`.

		fl : int
			Its value should be zero

		OutMatrix : :attr:`gcMapExplorer.lib.ccmap.CCMAP.matrix`
			Output matrix of Hi-C map to which normalized matrix is returned.

		bNoData : numpy.ndarray[bool]
			A numpy.array containing bool to show if rows/columns have missing data. It can be obtained from :meth:`remove_zeros`.

		"""
		minvalue = 0.0
		maxvalue = 0.0

		if self.memory == 'RAM':
			minvalue, maxvalue = self.run_for_RAM(A, fl, OutMatrix, bNoData)
		if self.memory == 'HDD':
			try:
				minvalue, maxvalue = self.run_for_HDD(A, fl, OutMatrix, bNoData)
			except (KeyboardInterrupt, SystemExit) as e:
				self.del_intermediate_files()
				raise e

		return minvalue, maxvalue
