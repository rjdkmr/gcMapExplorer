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

import os
import logging
import numpy as np
from scipy import stats
import numpy.ma as ma
import h5py

from . import ccmap as cmp
from . import gcmap as gmp
from . import ccmapHelpers as cmh
from . import util
from . import genomicsDataHandler as gdh

from gcMapExplorer.config import getConfig

# get configuration
config = getConfig()

def correlateCMaps(ccMapObjOne, ccMapObjTwo, ignore_triangular=True, diagonal_offset=1, corrType='pearson', blockSize=None, slideStepSize=1, cutoffPercentile=None, workDir=None, outFile=None, logHandler=None):
	"""To calculate correlation between two Hi-C maps

	This function can be used to calculate either Pearson or Spearman rank-order correlation
	between two Hi-C maps. It also ignore lower-trangular matrix with diagnonal offset to avoid duplicate and large values.

	Parameters
	----------
	ccMapObjOne : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		First :class:`gcMapExplorer.lib.ccmap.CCMAP` instance containing Hi-C data

	ccMapObjTwo : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		Second :class:`gcMapExplorer.lib.ccmap.CCMAP` instance containing Hi-C data

	ignore_triangular : bool
		Whether entire matrix is considered or only one half triangular region of matrixis considered.

	diagonal_offset : int
		If ``ignore_triangular=True``, it is used to determine how much bins are ignored from the diagonal in one half triangular region of matrix.
		``diagonal_offset = 0`` is the main diagonal, ``diagonal_offset > 0`` means ignore this many bins from the diagonal.

	corrType : str
		Correlation type. For Pearson and Spearman rank-order correlation, use ``pearson`` and ``spearman``, respectively.

	blockSize : str
		To calculate block-wise correlations by sliding block of given size along diagonals. It should be in resolution.
		For example, ``1mb``, ``500kb``, ``5mb``, ``2.5mb`` etc.
		If ``None``, correlation of whole map is calculated. Sliding step of block depends on ``slideStepSize``.

	slideStepSize : int
		Step-size in bins by which blocks will be shifted for block-wise correlation. If slideStepSize is large then blocks might not be overlapped.

	cutoffPercentile : float
		Cutoff percentile to discard values during correlation calculation. If a cutoff percentile is given, values less than this
		percentile value will not be considered during correlation calculation.

	workDir : str
		Name of working directory, where temporary files will be kept.If ``workDir = None``, file will be generated in OS based temporary directory.

	outFile : str
		Name of output file. Only written for block-wise correlation.


	Returns
	-------
	corr : float or list
		Correlation coefficient

	pvalue/centers : float or list
		If ``blockSize=None`` 2-tailed p-value is returned. For block-wise correlation, list of block-center is returned.


	.. seealso::
		* `scipy.stats.pearsonr <http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html#scipy.stats.pearsonr>`_ for Pearson correlation.
		* `scipy.stats.spearmanr <http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html#scipy.stats.spearmanr>`_ for Spearman rank-order correlation.


	"""
	# logger
	logger = logging.getLogger('correlateCMaps')
	if logHandler is not None:
		logger.propagate = False
		logger.addHandler(logHandler)
	logger.setLevel(logging.INFO)

	# To check, user have given correct keywords
	if not (corrType == 'pearson' or corrType == 'spearman' or corrType == 'kendell-tau'):
		raise NotImplementedError ('{0} not implemented. Use "pearson" or "spearman".' .format(corrType))

	# To check if resolution of maps are same
	if ccMapObjOne.binsize != ccMapObjTwo.binsize:
		raise AssertionError ("Resolution of first Hi-C map does not match with Second Hi-C map.")

	# Determine smallest shape between two maps
	if ccMapObjOne.shape[0] >= ccMapObjTwo.shape[0]:
		smallest_shape = ccMapObjTwo.shape[0]
	else:
		smallest_shape = ccMapObjOne.shape[0]

	# Calculation for whole map
	if blockSize is not None:
		logger.info(' Block-wise correlation with [{0}] block-size'.format(blockSize))

		# checking if block-size is larger than smallest shape
		size = int(util.resolutionToBinsize(blockSize)/ccMapObjOne.binsize)
		if size >= smallest_shape:
			raise AssertionError ("Size of input block [{0}] is larger than smallest size [{1}] of Hi-C map".format(size, smallest_shape))
		else:
			blockSize = size

		# Print some information about blocks
		logger.info(' Number of Blocks: {0} '.format( int(smallest_shape/blockSize) ) )
		logger.info(' Size of each Block in bins: {0} '.format( blockSize ) )
		if size-slideStepSize < 0:
			logger.info(' Blocks are not overlapping. Distance between two adjacent blocks in bins : {0}'.format(slideStepSize-size))
		else:
			logger.info(' Number of Overlapping bins between sliding blocks:  {0}' .format( size-slideStepSize ) )


	corr, pvalue, centers = None, None, None

	try:

		# generate boolean arrays to store mask
		m1 = cmh.MemoryMappedArray(shape=(smallest_shape, smallest_shape), dtype=np.bool, workDir=workDir)
		m2 = cmh.MemoryMappedArray(shape=(smallest_shape, smallest_shape), dtype=np.bool, workDir=workDir)
		mask = cmh.MemoryMappedArray(shape=(smallest_shape, smallest_shape), dtype=np.bool, workDir=workDir)

		# Determine masks for two maps and combine it
		m1.arr[:] = ccMapObjOne.matrix[:smallest_shape, :smallest_shape] < ccMapObjOne.minvalue
		m2.arr[:] = ccMapObjTwo.matrix[:smallest_shape, :smallest_shape] < ccMapObjTwo.minvalue

		if cutoffPercentile is not None:
			percentileOne = np.percentile(ma.array(ccMapObjOne.matrix[:smallest_shape, :smallest_shape], mask=m1.arr).compressed(), cutoffPercentile)
			percentileTwo = np.percentile(ma.array(ccMapObjTwo.matrix[:smallest_shape, :smallest_shape], mask=m2.arr).compressed(), cutoffPercentile)
			m1.arr[:] = ccMapObjOne.matrix[:smallest_shape, :smallest_shape] > percentileOne
			m2.arr[:] = ccMapObjTwo.matrix[:smallest_shape, :smallest_shape] > percentileTwo

			mask.arr[:] = ~( m1.arr | m2.arr )
			mask.arr[ np.nonzero( ccMapObjOne.matrix[:smallest_shape, :smallest_shape] == 0.0) ] = True
			mask.arr[ np.nonzero( ccMapObjTwo.matrix[:smallest_shape, :smallest_shape] == 0.0) ] = True
		else:
			mask.arr[:] = ( m1.arr & m2.arr )

		# Mask lower diagonal with diagonal_offset
		if ignore_triangular:
			mask.arr[np.tril_indices_from(mask.arr, k=diagonal_offset)] = True

		if blockSize is None:
			# This section for correlation between whole maps
			ccMapObjOneMa = ma.array(ccMapObjOne.matrix[:smallest_shape, :smallest_shape], mask=mask.arr)
			ccMapObjTwoMa = ma.array(ccMapObjTwo.matrix[:smallest_shape, :smallest_shape], mask=mask.arr)

			if ccMapObjOneMa[~ccMapObjOneMa.mask].shape[0] > 10:
				if corrType == 'pearson':
					corr, pvalue = stats.pearsonr(ccMapObjOneMa.compressed(), ccMapObjTwoMa.compressed())
				elif corrType == 'kendall-tau':
					corr, pvalue = stats.kendelltau(ccMapObjOneMa.compressed(), ccMapObjTwoMa.compressed())
				else:
					corr, pvalue = stats.spearmanr(ccMapObjOneMa.compressed(), ccMapObjTwoMa.compressed())

			else:
				corr = 0
				pvalue = 0

		else:
			# This section for block-wise correlation between maps
			corr, centers = [], []
			csp = 0           # Current start position
			cep = blockSize   # Current end position
			while(cep < smallest_shape):

				ccMapObjOneMa = ma.array(ccMapObjOne.matrix[csp:cep, csp:cep], mask=mask.arr[csp:cep, csp:cep])
				ccMapObjTwoMa = ma.array(ccMapObjTwo.matrix[csp:cep, csp:cep], mask=mask.arr[csp:cep, csp:cep])

				# At least 10 valid values are present
				if ccMapObjOneMa[~ccMapObjOneMa.mask].shape[0] > 10:
					if corrType == 'pearson':
						tcorr, tpvalue = stats.pearsonr(ccMapObjOneMa.compressed(), ccMapObjTwoMa.compressed())
					elif corrType == 'kendall-tau':
						corr, pvalue = stats.kendelltau(ccMapObjOneMa.compressed(), ccMapObjTwoMa.compressed())
					else:
						tcorr, tpvalue = stats.spearmanr(ccMapObjOneMa.compressed(), ccMapObjTwoMa.compressed())
				else:
					tcorr = 0

				if tcorr is ma.masked:
					corr.append(0.0)
				else:
					corr.append(float(tcorr))

				centers.append( (csp + (blockSize/2) ) * ccMapObjOne.binsize )

				csp = csp + slideStepSize
				cep = csp + blockSize

		del m1
		del m2
		del mask
	except (KeyboardInterrupt, SystemExit) as e:
		del m1
		del m2
		del mask
		raise e

	if logHandler is not None:
		logger.removeHandler( logHandler )

	if outFile is not None and blockSize is not None:
		fout = open(outFile, 'w')
		for i in range(len(corr)):
			if corr[i] != 0:
				fout.write('{0}\t{1}\n'.format(int(centers[i]), corr[i]))
		fout.close()


	if blockSize is None:
		return corr, pvalue
	else:
		return corr, centers


def correlateGCMaps(gcmapOne, gcmapTwo, outFile=None, blockSize=None, slideStepSize=1, name=None, cutoffPercentile=None, ignore_triangular=True, diagonal_offset=1, corrType='pearson', workDir=None, logHandler=None):
	"""To calculate correlation between common Hi-C maps from two gcmap files

	This function can be used to calculate either Pearson or Spearman rank-order correlation between common maps present in two gcmap files.
	It also ignore lower-trangular matrix with diagnonal offset to avoid duplicate and large values.

	.. note:: If block-wise correlation calculation will be initiated by ``blockSize`` option, a ``outFile`` and ``name`` is also
			required for further processing. The block-wise correlation will stored in output HDF5 format file.


	Parameters
	----------
	gcmapOne : str
		First gcmap file.

	gcmapTwo : str
		Second gcmap file

	outFile : str
		Name of output file. Only written for block-wise correlation.

	blockSize : str
		To calculate block-wise correlations by sliding block of given size along diagonals. It should be in resolution.
		For example, ``1mb``, ``500kb``, ``5mb``, ``2.5mb`` etc.
		If ``None``, correlation of whole map is calculated. Sliding step of block depends on ``slideStepSize``.

	slideStepSize : int
		Step-size in bins by which blocks will be shifted for block-wise correlation. If slideStepSize is large then blocks might not be overlapped.

	name : str
		Title of dataset in HDF5 output file. If ``blockSize`` option is used, ``name`` is an essential argument.

	cutoffPercentile : float
		Cutoff percentile to discard values during correlation calculation. If a cutoff percentile is given, values less than this
		percentile value will not be considered during correlation calculation.

	ignore_triangular : bool
		Whether entire matrix is considered or only one half triangular region of matrixis considered.

	diagonal_offset : int
		If ``ignore_triangular=True``, it is used to determine how much bins are ignored from the diagonal in one half triangular region of matrix.
		``diagonal_offset = 0`` is the main diagonal, ``diagonal_offset > 0`` means ignore this many bins from the diagonal.

	corrType : str
		Correlation type. For Pearson and Spearman rank-order correlation, use ``pearson`` and ``spearman``, respectively.

	workDir : str
		Name of working directory, where temporary files will be kept.If ``workDir = None``, file will be generated in OS based temporary directory.

	Returns
	-------
	mapList : list
		List of chromosomes

	corrs : list
		Correlation coefficient of each chromosome

	pvalue : list
		2-tailed p-value for correlation coefficient of each chromosome.

	"""

	# logger
	logger = logging.getLogger('correlateGCMaps')
	if logHandler is not None:
		logger.propagate = False
		logger.addHandler(logHandler)
	logger.setLevel(logging.INFO)

	if blockSize is not None:
		if outFile is None:
			raise ValueError('No output file provided.')
		if name is None:
			raise ValueError('No data name provided.')

	mapListOne = []
	hdf5 = h5py.File(gcmapTwo)
	for key in hdf5:
		if hdf5[key].attrs['xlabel'] == hdf5[key].attrs['ylabel']:
			mapListOne.append(key)
	hdf5.close()

	mapListTwo = []
	hdf5 = h5py.File(gcmapOne)
	for key in hdf5:
		if hdf5[key].attrs['xlabel'] == hdf5[key].attrs['ylabel']:
			mapListTwo.append(key)
	hdf5.close()

	mapList = util.sorted_nicely( list(set(mapListOne).intersection(mapListTwo)) )

	if blockSize is not None:
		h5Handle = gdh.HDF5Handler(outFile)
	else:
		corrs, pvalues = [], []

	for key in mapList:

		logger.info(' Performing calculation for {0} ...'.format(key))


		ccMapObjOne = gmp.loadGCMapAsCCMap(gcmapOne, mapName=key, workDir=workDir)
		ccMapObjOne.make_readable()

		ccMapObjTwo = gmp.loadGCMapAsCCMap(gcmapTwo, mapName=key, workDir=workDir)
		ccMapObjTwo.make_readable()


		#ccMapObjOne = gmp.GCMAP(gcmapOne, mapName=key)
		#ccMapObjTwo = gmp.GCMAP(gcmapTwo, mapName=key)

		binsize = ccMapObjOne.binsize

		dataArray = np.zeros(shape=(max(ccMapObjOne.shape[0], ccMapObjTwo.shape[0]),) )

		try:
			results = correlateCMaps(ccMapObjOne, ccMapObjTwo, ignore_triangular=ignore_triangular, diagonal_offset=diagonal_offset, corrType=corrType,
		                blockSize=blockSize, cutoffPercentile=cutoffPercentile, slideStepSize=slideStepSize, workDir=workDir, outFile=None, logHandler=logHandler)
		except	(KeyboardInterrupt, SystemExit) as e:
			del ccMapObjOne
			del ccMapObjTwo
			raise e

		del ccMapObjOne
		del ccMapObjTwo

		corr, centers = None, None
		if results is not None:
			corr, centers = results

		if blockSize is not None:
			centers = np.asarray(centers) / binsize
			centers = centers.astype(int)
			dataArray[(centers,)] = corr
			h5Handle.addDataByArray(key, util.binsizeToResolution(binsize), name, dataArray)
		else:
			corrs.append(corr)
			pvalues.append(centers)

		logger.info('     Finished calculation for {0} ...\n'.format(key))

	if blockSize is not None:
		h5Handle.close()

	if blockSize is None:
		return mapList, corrs, pvalues


def correlateCMapsBinWise(ccMapObjOne, ccMapObjTwo, corrType='pearson', cutoffPercentile=None, workDir=None, outFile=None, logHandler=None):
	"""To calculate correlation between two Hi-C maps

	This function can be used to calculate either Pearson or Spearman rank-order correlation
	between two Hi-C maps. It also ignore lower-trangular matrix with diagnonal offset to avoid duplicate and large values.


	Parameters
	----------
	ccMapObjOne : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		First :class:`gcMapExplorer.lib.ccmap.CCMAP` instance containing Hi-C data
	ccMapObjTwo : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		Second :class:`gcMapExplorer.lib.ccmap.CCMAP` instance containing Hi-C data
	ignore_triangular : bool
		Whether entire matrix is considered or only one half triangular region of matrixis considered.
	diagonal_offset : int
		If ``ignore_triangular=True``, it is used to determine how much bins are ignored from the diagonal in one half triangular region of matrix.

		``diagonal_offset = 0`` is the main diagonal, ``diagonal_offset > 0`` means ignore this many bins from the diagonal.

	corrType : str
		Correlation type. For Pearson and Spearman rank-order correlation, use ``pearson`` and ``spearman``, respectively.
	blockSize : str
		To calculate block-wise correlations by sliding block of given size along diagonals. It should be in resolution.
		For example, ``1mb``, ``500kb``, ``5mb``, ``2.5mb`` etc.
		If ``None``, correlation of whole map is calculated. Sliding step of block depends on ``slideStepSize``.
	slideStepSize : int
	  	Step-size in bins by which blocks will be shifted for block-wise correlation. If slideStepSize is large then blocks might not be overlapped.

	workDir : str
		Name of working directory, where temporary files will be kept.If ``workDir = None``, file will be generated in OS based temporary directory.

	outFile : str
		Name of output file. Only written for block-wise correlation.

	Returns
	-------
	corr : float or list
		Correlation coefficient
	pvalue/centers : float or list
		If ``blockSize=None`` 2-tailed p-value is returned. For block-wise correlation, list of block-center is returned.


	.. seealso::
		* `scipy.stats.stats.pearsonr <http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html#scipy.stats.pearsonr>`_ for Pearson correlation.
		* `scipy.stats.stats.spearmanr <http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html#scipy.stats.spearmanr>`_ for Spearman rank-order correlation.


	"""
	# logger
	logger = logging.getLogger('correlateCMaps')
	if logHandler is not None:
		logger.propagate = False
		logger.addHandler(logHandler)
	logger.setLevel(logging.INFO)

	# To check, user have given correct keywords
	if not (corrType == 'pearson' or corrType == 'spearman' or corrType == 'kendell-tau'):
		raise NotImplementedError ('{0} not implemented. Use "pearson" or "spearman".' .format(corrType))

	# To check if resolution of maps are same
	if ccMapObjOne.binsize != ccMapObjTwo.binsize:
		raise AssertionError ("Resolution of first Hi-C map does not match with Second Hi-C map.")

	# Determine smallest shape between two maps
	if ccMapObjOne.shape[0] >= ccMapObjTwo.shape[0]:
		smallest_shape = ccMapObjTwo.shape[0]
	else:
		smallest_shape = ccMapObjOne.shape[0]

	corr, centers = [], []
	try:

		# generate boolean arrays to store mask
		m1 = cmh.MemoryMappedArray(shape=(smallest_shape, smallest_shape), dtype=np.bool, workDir=workDir)
		m2 = cmh.MemoryMappedArray(shape=(smallest_shape, smallest_shape), dtype=np.bool, workDir=workDir)
		mask = cmh.MemoryMappedArray(shape=(smallest_shape, smallest_shape), dtype=np.bool, workDir=workDir)

		# Determine masks for two maps and combine it
		m1.arr[:] = ccMapObjOne.matrix[:smallest_shape, :smallest_shape] == ccMapObjOne.minvalue
		m2.arr[:] = ccMapObjTwo.matrix[:smallest_shape, :smallest_shape] == ccMapObjTwo.minvalue

		if cutoffPercentile is not None:
			percentileOne = np.percentile(ma.array(ccMapObjOne.matrix[:smallest_shape, :smallest_shape], mask=m1.arr).compressed(), cutoffPercentile)
			percentileTwo = np.percentile(ma.array(ccMapObjTwo.matrix[:smallest_shape, :smallest_shape], mask=m2.arr).compressed(), cutoffPercentile)
			m1.arr[:] = ccMapObjOne.matrix[:smallest_shape, :smallest_shape] <= percentileOne
			m2.arr[:] = ccMapObjTwo.matrix[:smallest_shape, :smallest_shape] <= percentileTwo

		mask.arr[:] = ( m1.arr | m2.arr )

		for b in range(smallest_shape):
			ccMapObjOneMa = ma.array(ccMapObjOne.matrix[b], mask=mask.arr[b])
			ccMapObjTwoMa = ma.array(ccMapObjTwo.matrix[b], mask=mask.arr[b])

			tcorr = 0
			if ccMapObjOneMa[~ccMapObjOneMa.mask].shape[0] > 10:
				if corrType == 'pearson':
					tcorr, _ = stats.pearsonr(ccMapObjOneMa.compressed(), ccMapObjTwoMa.compressed())
				elif corrType == 'kendall-tau':
					tcorr, _ = stats.kendelltau(ccMapObjOneMa.compressed(), ccMapObjTwoMa.compressed())
				else:
					tcorr, _ = stats.spearmanr(ccMapObjOneMa.compressed(), ccMapObjTwoMa.compressed())

			if tcorr is ma.masked:
				corr.append(0.0)
			else:
				corr.append(float(tcorr))
			centers.append( b * ccMapObjOne.binsize )

		del m1
		del m2
		del mask
	except (KeyboardInterrupt, SystemExit) as e:
		del m1
		del m2
		del mask
		raise e

	if logHandler is not None:
		logger.removeHandler( logHandler )

	if outFile is not None:
		fout = open(outFile, 'w')
		for i in range(len(corr)):
			if corr[i] != 0:
				fout.write('{0}\t{1}\n'.format(int(centers[i]), corr[i]))
		fout.close()

	return corr, centers


def getAvgContactByDistance(ccmaps, stats='median', removeOutliers=False, outliersThershold=3.5):
	"""To calcualte average contact as a function of distance

	Parameters
	----------
	ccmaps : :class:`gcMapExplorer.lib.ccmap.CCMAP` or list[:class:`gcMapExplorer.lib.ccmap.CCMAP`]
		Input contact maps

	stats : str
		Statistics for scaling. Accepted methods are ``mean``
		and ``median``.

	removeOutliers : bool
		If ``True``, outliers will be removed before calculating input
		statistics.

	outliersThershold : float
		The modified z-score to use as a threshold. Observations with
		a modified z-score (based on the median absolute deviation) greater
		than this value will be classified as outliers.

	Returns
	-------
		avg_contacts : numpy.array
		A one-dimensional numpy array containing average contacts, where index is distance between two locations for given resolution/binsize.
		For example, if ``ccmap.binsize=100000`` and ``avg_contacts[4]=1234.56``, then at distance of 400000 b, average contact is ``1234.56``.

	"""
	ccmapList = []
	if isinstance(ccmaps, cmp.CCMAP) or isinstance(ccmaps, gmp.GCMAP):
		ccmapList = [ccmaps]
	elif isinstance(ccmaps, list):
		for ccmap in ccmaps:
			if not isinstance(ccmap, cmp.CCMAP):
				raise TypeError ('Not a list of CCMAP instances')

		ccmapList = ccmaps
	else:
		raise TypeError ('Not a CCMAP instance or a list of CCMAP instances')

	# Determine max shape
	largest_shape = 0
	for ccmap in ccmapList:
		largest_shape = max(ccmap.shape[0], largest_shape)

	# Start calculation stats
	avg_contacts = []
	for i in range(largest_shape):

		data = np.empty(0)

		for t in range(len(ccmapList)):
		    idx = util.kth_diag_indices(i, ccmapList[t].matrix)

		    # Mask all zero values, then calculate stats
		    marray = ma.masked_equal(ccmapList[t].matrix[idx], 0.0, copy=False)
		    data = np.hstack((marray.compressed(), data))

		if data.shape[0] >= 1:

			notOutlierIndex = None
			outliersCount = 0
			if removeOutliers:
				outliers = util.detectOutliers1D(data, thresh=outliersThershold)
				outliersCount = np.sum(outliers)
				notOutlierIndex = np.nonzero(~outliers)

			if stats == 'median':
				if removeOutliers and outliersCount>0:
					t_avg = np.median(data[notOutlierIndex])
				else:
					t_avg = np.median(data)
			else:
				if removeOutliers and outliersCount>0:
					t_avg = np.mean(data[notOutlierIndex])
				else:
					t_avg = np.mean(data)

			avg_contacts.append( t_avg )

		else:
			avg_contacts.append( 0.0 )

		del data

	return np.asarray(avg_contacts)
