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
import re, sys, os, shutil
import tarfile
import subprocess as sub
import string, random, tempfile
import json
import logging
import gzip
import h5py

from .util import *

from gcMapExplorer.config import getConfig

# get configuration
config = getConfig()

from . import ccmapHelpers as cmh

dtype_npBINarray = 'float32'

logger = logging.getLogger('ccmap')
logger.setLevel(logging.INFO)


class CCMAP:
	"""This class contains variables to store Hi-C Data.

	The class is instantiated by two methods:
		>>> ccMapObj = gcMapExplorer.lib.ccmap.CCMAP()
		>>> ccMapObj = gcMapExplorer.lib.ccmap.CCMAP(dtype='float32')

	Parameters
	----------
	dtype : str, Optional
			Data type for matrix. [``Default='float32'``]


	Attributes
	----------

	path2matrix : str
		Path to numpy array binary file on local file system
	yticks		:	list
		Minimum and maximum locations along Y-axis. e.g. ``yticks=[0, 400000]``
	xticks		:	list
		Minimum and maximum locations along X-axis. e.g. ``xticks=[0, 400000]``
	binsize		:  int
		Resolution of data. In case of 10kb resolution, binsize is 10000.
	title		:	str
		Title of the data
	xlabel		:	str
		Title for X-axis
	ylabel		:	str
		Title for Y-axis
	shape		:	tuple
		Overall shape of matrix
	minvalue	:	float
		Minimum value in matrix
	maxvalue	:	float
		Maximum value in matrix
	matrix		: 	numpy.memmap
		A memmap object pointing to matrix.

		HiC map data is saved as a numpy array binary file on local file system. This file can be only read after mapping to a numpy memmap object. After mapping, ``matrix``
		can be used as a numpy array.
		The file name is randomly generated as ``npBinary_XXXXXXXXXX.tmp``, where ``X`` can be a alphanumeric charecter. Please see details in |numpy memmap|.

		When ccmap is saved, this file is renamed with '.npbin' extension.

	bNoData		:	numpy.ndarray
		A boolean numpy array of matrix shape
	bLog		:	bool
		If values in matrix are in log
	state       :	str
		State of CCMAP object

		This keyword stores the state of the object. The state ensures when the numpy array binary file should be deleted from the local file system.

		Three keywords are used:
			* ``temporary`` :
				When object is created, it is in temporary state. After executing the script, numpy array binary file is automatically deleted from the local file system.

			* ``saved``:
				When a temporary or new object is saved, the numpy array binary file is copied to the destination directory and state is changed to saved.
				However, after saving, state become temporary. This method ensures that the saved copy is not deleted and only temporary copy is deleted after executing the script.

				When a already saved object is loaded, it is in saved state. The numpy array binary file is read from the original location and remaines saved
				at the original location after executing the script.

			* ``compressed``:
				When object is saved and numpy array binary file is simultaneously compressed, it is saved as compressed state. When this object is loaded,
				the numpy array binary file is decompressed into the working directory. This decompressed file is automatically deleted after
				execution of script while compressed file remains saved at the original location.

	dtype       :	str
		Data type of matrix

	"""

	def __init__(self, dtype=dtype_npBINarray):
		self.path2matrix = ''
		self.xticks = None
		self.yticks = None
		self.binsize = 0
		self.title = None
		self.xlabel = None
		self.ylabel = None
		self.shape = None
		self.minvalue = 9e10
		self.maxvalue = -9e10
		self.bNoData = None
		self.bLog = False
		self.dtype = dtype
		self.matrix = None

		self.state = 'temporary'

	def __del__(self):
		if self.matrix is not None:
			self.matrix._mmap.close()
		if self.matrix is not None:
			self.matrix
		if self.state == 'temporary' or self.state == 'compressed':
			if os.path.exists(self.path2matrix):
				os.remove(self.path2matrix)
		del self.xticks
		del self.yticks
		del self.shape

	def get_ticks(self, binsize=None):
		"""To get xticks and yticks for the matrix

		Parameters
		----------
		binsize : int
			Number of base in each bin or pixel or box of contact map.

		Returns
		-------
		xticks	:	numpy.array
			1D array containing positions along X-axis
		yticks	:	numpy.array
			1D array containing positions along X-axis

		"""
		if binsize is None:
			xticks = np.arange(self.xticks[0], self.xticks[1], self.binsize)
			yticks = np.arange(self.yticks[0], self.yticks[1], self.binsize)
		else:
			yticks = np.arange(self.yticks[0], self.yticks[1], binsize)
			xticks = np.arange(self.xticks[0], self.xticks[1], binsize)
		return xticks, yticks

	def gen_matrix_file(self, workDir=None):
		# Working and output directory
		if workDir is None:
			workDir = config['Dirs']['WorkingDirectory']

		# Only filename is generated
		fd, self.path2matrix = tempfile.mkstemp(prefix='gcx_npBinary_', suffix='.tmp', dir=workDir)
		os.close(fd)
		if os.path.isfile(self.path2matrix):
			os.remove(self.path2matrix)

	def make_readable(self):
		"""Enable reading the numpy array binary file.

		Matrix file is saved on local file system. This file can be only read after mapping to a memmap object. This method maps the numpy memmap object to self.matrix variable.
		After using this method, :attr:`gcMapExplorer.lib.ccmap.CCMAP.matrix` can be used directly as similar to numpy array. Please see details in |numpy memmap|

		"""
		if self.matrix is not None:
			del self.matrix
		self.matrix = np.memmap(self.path2matrix, dtype=self.dtype, mode='r', shape=self.shape)

	def make_unreadable(self):
		"""Disable reading the numpy array binary file from local file system

		"""
		if self.matrix is not None:
			del self.matrix
			self.matrix = None

	def make_writable(self):
		"""Create new numpy array binary file on local file system and enable reading/writing to this file

		.. note::
		 	If a matrix file with similar name is already present, old file will be backed up.

		"""
		if self.matrix is not None:
			del self.matrix
		if os.path.exists(self.path2matrix):
			gen_backup(self.path2matrix)
		self.matrix = np.memmap(self.path2matrix, dtype=self.dtype, mode='w+', shape=self.shape)

	def make_editable(self):
		"""Enable editing numpy array binary file

		"""
		if self.matrix is not None:
			del self.matrix
		self.matrix = np.memmap(self.path2matrix, dtype=self.dtype, mode='r+', shape=self.shape)

	def copy(self, fill=None):
		"""To create a new copy of CCMAP object

		This method can be used to create a new copy of :class:`gcMapExplorer.lib.ccmap.CCMAP`.
		A new numpy array binary file will be created and all values from old file will be copied.

		Parameters
		----------
		fill : float
			Fill map with the value. If not given, map values will be copied.

		"""

		ccMapObj = CCMAP(dtype=self.dtype)

		for key in self.__dict__:
			ccMapObj.__dict__[key] = self.__dict__[key]

		# Generaing temporary numpy array file
		dirname = os.path.dirname( self.path2matrix )
		ccMapObj.gen_matrix_file(workDir=dirname)

		#ccMapObj.path2matrix = os.path.join(os.getcwd(), 'nparray_' + getRandomName() + '.bin')

		self.make_readable()
		ccMapObj.make_writable()

		if fill is None:
			ccMapObj.matrix[:] = self.matrix[:]
		else:
			ccMapObj.matrix.fill(fill)

		ccMapObj.matrix.flush()
		self.make_unreadable()
		ccMapObj.make_unreadable()

		return ccMapObj

def jsonify(ccMapObj):
	"""Changes data type of attributes in CCMAP object for |json link|.

	Before saving the CCMAP object, its attributes data types are neccessary to change because few data types are not supported by json.

	Therefore, it is converted into other data types which are supported by json.
	These are the following attributes which are changed:

	+-------------+----------------------------------+----------------------------+
	| Attributes  | Original                         | modified                   |
	+-------------+----------------------------------+----------------------------+
	| bNoData     | numpy boolean array              | string of 0 and 1          |
	+-------------+----------------------------------+----------------------------+
	| xticks      | list of integer                  | list of string             |
	+-------------+----------------------------------+----------------------------+
	| yticks      | list of integer                  | list of string             |
	+-------------+----------------------------------+----------------------------+
	| minvalue    | float                            | string                     |
	+-------------+----------------------------------+----------------------------+
	| maxvalue    | float                            | string                     |
	+-------------+----------------------------------+----------------------------+
	| binsize     | integer                          | string                     |
	+-------------+----------------------------------+----------------------------+
	| shape       | tuple of integer                 | list of string             |
	+-------------+----------------------------------+----------------------------+
	| dtype       | Numpy dtype                      | string                     |
	+-------------+----------------------------------+----------------------------+


	..	warning::
		If a object is passed through this method, it should be again passed through :meth:`gcMapExplorer.lib.ccmap.dejsonify` for any further use.
		Otherwise, this object cannot be used in any other methods because of the attributes data type modifications.


	Parameters
	----------
	ccMapObj	:	:class:`gcMapExplorer.lib.ccmap.CCMAP`
		A CCMAP object

	Returns
	-------
	None



	"""
	if ccMapObj.bNoData is not None:
		ccMapObj.bNoData = ''.join(list(map(str, list(np.array(ccMapObj.bNoData).astype(int)) )))

	if ccMapObj.xticks is not None:
		ccMapObj.xticks = list(map(str, ccMapObj.xticks))

	if ccMapObj.yticks is not None:
		ccMapObj.yticks = list(map(str, ccMapObj.yticks))

	ccMapObj.minvalue = str(float(ccMapObj.minvalue))
	ccMapObj.maxvalue = str(float(ccMapObj.maxvalue))
	ccMapObj.shape = list(map(str, ccMapObj.shape))

	if hasattr(ccMapObj, 'binsize'):
		ccMapObj.binsize = str(ccMapObj.binsize)

	if isinstance(ccMapObj.dtype, np.dtype):
		ccMapObj.dtype = ccMapObj.dtype.name

def dejsonify(ccMapObj, json_dict=None):
	"""Change back the data type of attributes in CCMAP object.

	Before loading the CCMAP object, its attributes data types are neccessary to change back.

	Therefore, it is converted into original data types as shown in a table (see :meth:`gcMapExplorer.lib.ccmap.jsonify`)

	Parameters
	----------
	ccMapObj	:	:class:`gcMapExplorer.lib.ccmap.CCMAP`
		A CCMAP object

	json_dict	:	dict, Optional
		A directory obtained after loading the saved file through json.

	"""
	if json_dict is not None:
		for key in json_dict:
			ccMapObj.__dict__[key] = json_dict[key]

	if ccMapObj.bNoData is not None:
		ccMapObj.bNoData = np.array( list(map(int, list(ccMapObj.bNoData))), dtype=bool)

	if ccMapObj.xticks is not None:
		ccMapObj.xticks = list(map(int, ccMapObj.xticks))

	if ccMapObj.yticks is not None:
		ccMapObj.yticks = list(map(int, ccMapObj.yticks))

	ccMapObj.minvalue = float(ccMapObj.minvalue)
	ccMapObj.maxvalue = float(ccMapObj.maxvalue)
	ccMapObj.shape = tuple(map(int, ccMapObj.shape))

	if hasattr(ccMapObj, 'binsize'):
		ccMapObj.binsize = int(ccMapObj.binsize)


def save_ccmap(ccMapObj, outfile,  compress=False, logHandler=None):
	""" Save CCMAP object on file

	CCMAP object can be saved as file for easy use. |json link| is used to save the object.
	The binary numpy array file is copied in the destination directory. If ``compress=True``, the array file will be compressed in gzip format.

	.. note::
		* Compression significantly reduces the array file size. However, its loading is slow during initiation when file is decompressed.

		* After loading, the decompressed binary numpy array file takes additional memory on local file system.




	Parameters
	----------
	ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		A CCMAP object, which has to be saved
	outfile : str
		Name of output file including path to the directory/folder where file should be saved.
	compress	:	bool
		If ``True``, numpy array file will be compressed.

	Returns
	-------
	None

	"""

	# logger
	logHandler = logHandler
	logger = logging.getLogger('save_ccmap')
	if logHandler is not None:
		logger.propagate = False
		logger.addHandler(logHandler)
	logger.setLevel(logging.INFO)

	ccMapObj.make_unreadable()

	# Get file extension of output file name. If output name does not contain .ccmap extension it will be added
	file_extension = os.path.splitext(outfile)[1]
	if file_extension != '.ccmap':
		outfile = outfile + '.ccmap'

	# Getting output directory and generating path name for numpy array file
	outdir = os.path.dirname( os.path.abspath( os.path.expanduser(outfile)) )
	basename = os.path.basename( outfile )

	nparrayFileName = os.path.splitext(basename)[0] + '.npbin'
	nparrayfile = os.path.join(outdir,  nparrayFileName)

	logger.info(" Saving ccmap to file [{0}] and [{1}] ..." .format(outfile, nparrayfile))
	# Copying numpy array file to the output directory
	if ccMapObj.path2matrix == nparrayfile:
		logging.warning('File {0} already exists !!! If you have used make_writable or make_editable, this file is already modified and saved...' .format(ccMapObj.path2matrix))
	else:
		# Generating backup of old numpy array file
		if os.path.exists(nparrayfile):
			gen_backup(nparrayfile)
		shutil.copy(ccMapObj.path2matrix, nparrayfile)

	# Generating .ccmap backup file if it is already present in the output directory
	if os.path.exists(outfile):
		gen_backup(outfile)

	# Only store numpy array file in object because this file and .ccmap file should be present in same directory
	t_path2matrix = ccMapObj.path2matrix
	ccMapObj.path2matrix = nparrayFileName

	# converting data type before saving .ccmap
	jsonify(ccMapObj)
	if compress:
		ccMapObj.state = 'compressed'
	else:
		ccMapObj.state = 'saved'

	# saving .ccmap
	fout =  open( os.path.abspath( os.path.expanduser(outfile)), "w" )
	json.dump(ccMapObj.__dict__, fout, indent=4, separators=(',', ':'))
	fout.close()

	logger.info(" Compressing [{0}] ..." .format(nparrayfile))

	# Compressing numpy array file and removing original file
	if compress:

		if os.path.exists(nparrayfile+'.gz'):
			gen_backup(nparrayfile+'.gz')

		with open(nparrayfile, 'rb') as f_in:
			with gzip.open(nparrayfile + '.gz', 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)

		if os.path.exists(nparrayfile):
			os.remove(nparrayfile)

	logger.info("       Finished!!!\n")

	# Reverting the path of matrix file so that this object can be still used
	ccMapObj.path2matrix = t_path2matrix
	ccMapObj.state = 'temporary'
	dejsonify(ccMapObj)

def load_ccmap(infile, workDir=None):
	""" Load CCMAP object from an input file

	CCMAP object can be created from the input file, which was earlier saved using :meth:`gcMapExplorer.lib.ccmap.save_ccmap`.
	If the binary numpy array is compressed, this file is automatically extracted in the current working directory. After completion of the execution,
	this decompressed file will be automatically deleted. The compressed saved file will be remained unchanged.

	Parameters
	----------
	infile : str
		Name of the inputp file including path to the directory/folder where file is saved.

	workDir : str
		Name of working directory, where temporary files will be kept.If ``workDir = None``, file will be generated in OS based temporary directory.

	Returns
	-------
	ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		A CCMAP object

	"""

	# Working directory
	if workDir is None:
		workDir = config['Dirs']['WorkingDirectory']

	# Open .ccmap file and load all keywords
	fin = open( infile, "r" )
	json_dict = json.load( fin )
	fin.close()

	# New object
	ccMapObj = CCMAP(dtype=np.dtype(json_dict['dtype']))

	# Generating the object from the json dictionary
	dejsonify(ccMapObj, json_dict)

	# Generating full path to numpy array file
	indir = os.path.dirname( os.path.abspath( os.path.expanduser(infile) ) )
	nparrayInFile = os.path.join( indir, ccMapObj.path2matrix)

	# if numpy array file is compressed, extract it and move to current working directory
	# else use this array file directly from original location
	if os.path.exists(nparrayInFile + '.gz'):
		ccMapObj.gen_matrix_file(workDir=workDir)      # Generate matrix file

		# try to catch system-exit and keyboard intrrupt from user. Have to delete manually created file
		try:
			with gzip.open(nparrayInFile + '.gz', 'rb') as f_in:
				with open(ccMapObj.path2matrix, 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
		except (KeyboardInterrupt, SystemExit) as e:
			os.remove(ccMapObj.path2matrix)
			raise e

		# Making this obeject as temorary
		ccMapObj.state = 'temporary'
	else:
		ccMapObj.path2matrix = nparrayInFile

	return ccMapObj

def export_cmap(ccmap, outfile, doNotWriteZeros=True):
	"""To export ``.ccmap`` as text file

	This function export ``.ccmap`` as coordinate list (COO) format sparse matrix file.
	In COO format, lists of (row, column, value) as three tab seprated columns are written in output file.

	Parameters
	----------
	ccmap : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		An instance of :class:`gcMapExplorer.lib.ccmap.CCMAP`, which is need to be exported.
	outfile : str
		Output file name.
	doNotWriteZeros : bool
		Do not write Zero values. It reduces memory of file.

	"""

	fout = open(outfile, 'w')

	if hasattr(ccmap, 'make_readable'):
		ccmap.make_readable()

	try:
		for i in range(ccmap.shape[0]):
			for j in range(i+1):
				if doNotWriteZeros:
					if ccmap.matrix[i][j] != 0:
						fout.write('{0}\t{1}\t{2}\n'.format(i, j, ccmap.matrix[i][j] ))
					else:
						fout.write('{0}\t{1}\t{2}\n'.format(i, j, ccmap.matrix[i][j] ))

	except (SystemExit, KeyboardInterrupt) as e:
		del ccmap
		raise e

def checkCCMapObjectOrFile(ccMap, workDir=None):
	"""Check whether ccmap is a obejct or file

	It can be used to check whether input is a :class:`gcMapExplorer.lib.ccmap.CCMAP`
	or a ccmap file.

	It returns the :class:`gcMapExplorer.lib.ccmap.CCMAP` and input type name:
	i.e. ``File`` or ``Object`` as an identification keyword for the input.

	In case if ``ccMap`` argument is a filename, this file
	will be opened as a :class:`gcMapExplorer.lib.ccmap.CCMAP` object and
	will be returned with ``ccmapType`` as ``File``.

	In case if ``ccMap`` argument is a :class:`gcMapExplorer.lib.ccmap.CCMAP`
	object, this file, same object will be returned with ``ccmapType`` as
	``Object``.

	Parameters
	----------
	ccMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or str
		CCMAP object or ccmap file.

	workDir : str
		Path to the directory where temporary intermediate files are generated.
		If ``None``, files are generated in the temporary directory according to the main configuration.

	Returns
	-------
	ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		CCMAP object

	ccmapType : str
		'File' or 'Object'

	"""
	# Check whether input is a file or a obejct
	ccmapType = 'Object'
	ccMapObj = None
	if not isinstance(ccMap, CCMAP):
		ccmapType = 'File'
		if not os.path.isfile(ccMap):
			logger.info(' {0} file not found.'.format(ccMap))
			logger.info(' Not able to normalize.')
			return None
		else:
			ccMapObj = load_ccmap(ccMap, workDir=workDir)
	else:
		ccMapObj = ccMap

	return ccMapObj, ccmapType

def downSampleCCMap(cmap, level=2, workDir=None):

	level = int(level)
	# Working and output directory
	if workDir is None:
		workDir = config['Dirs']['WorkingDirectory']

	# To determine output shape but it does not contain remainder regions
	# subtracted one as real data start from first index
	outShapeX = int( np.ceil(float(cmap.shape[1] - 1 )/level) )
	outShapeY = int( np.ceil(float(cmap.shape[0] - 1 )/level) )

	# To determine remainder regions(sourceMasked.max() - sourceMasked.min())
	pad_size_outer = int( ( outShapeY * level) - (cmap.shape[0] - 1) )
	pad_size_inner = int( ( outShapeX * level) - (cmap.shape[1] - 1) )

	# Otput shape
	outputShape = (outShapeX + 1 , outShapeY + 1)

	#rint(pad_size_outer, pad_size_inner, outputShape)

	outCmap = cmp.CCMAP(dtype=np.float)
	outCmap.shape = outputShape
	outCmap.gen_matrix_file(workDir=workDir)
	outCmap.make_writable()

	# Ignore first row and column --- no data
	count = 1
	i = 1
	while i < cmap.shape[0]:

	    # Here handles inner
	    if pad_size_inner != 0:
	        zero_part = []
	        for n in range(level):
	            zero_part.append( np.zeros(pad_size_inner) )

	        # Here handles outer
	        if i+level < cmap.shape[0]:
	            x_padded = np.hstack( (cmap.matrix[i:i + level, 1:], zero_part) )
	        else:
	            x_padded = np.hstack( (cmap.matrix[i:, 1:], zero_part[:cmap.matrix[i:].shape[0]] ) )

	        # print(x_padded.shape)

	        r = x_padded.sum(axis=0).reshape(-1, level).sum(axis=1)

	        #rint(r.shape)

	    else:

	        # Here handles outer
	        if i+level < cmap.shape[0]:
	            r = cmap.matrix[i:i+level, 1:].sum(axis=0).reshape(-1, level).sum(axis=1)
	        else:
	            remain = cmap.shape[0]-i
	            r = cmap.matrix[i:, 1:].sum(axis=0).reshape(-1, remain).sum(axis=1)

	    i = i + level

	    outCmap.matrix[count, 1:] = r
	    count = count + 1

	outCmap.matrix.flush()

	# Assign minimum value
	ma = np.ma.masked_equal(outCmap.matrix, 0.0, copy=False)
	outCmap.minvalue = ma.min()
	del ma

	# Assign Maximum value
	outCmap.maxvalue = np.amax(outCmap.matrix)

	# New binsize
	outCmap.binsize = cmap.binsize * level

	# xlabel and ylabel
	outCmap.xlabel = cmap.xlabel
	outCmap.ylabel = cmap.ylabel

	return outCmap

def downSample1D(array, level=2, func='max'):
	""" Downsample or coarse a one-dimensional array
	"""

	level = int(level)
	outShape = int( np.ceil(float(array.shape[0] - 1 )/level) )
	pad_size = int( ( outShape * level) - (array.shape[0] - 1) )

	outputShape = outShape + 1

	outputArray = np.zeros(outputShape, dtype=array.dtype)
	if pad_size != 0:
		x_padded = np.hstack( (array[1:], np.zeros(pad_size)) )
	else:
		x_padded = array[1:]


	if func =='sum':
		outputArray[1:] = x_padded.reshape(-1, level).sum(axis=1)
	elif func =='mean':
		outputArray[1:] = x_padded.reshape(-1, level).mean(axis=1)
	else:
		outputArray[1:] = x_padded.reshape(-1, level).max(axis=1)

	return outputArray


def smoothen_map(ccMapObj, filter_pass=1):
	""" Smoothen map using unweighted sliding-average smooth method.

	This function can be used to smooth out the map. It implements unweighted sliding-average smooth method.

	Parameters
	----------
	ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		A CCMAP object, which has to be smoothen.
	filter_pass : int
		Number of pass to smooth. Larger the number of pass, more smoother is the map. However, pattern in maps may be difficult to identify when
		pass number is high.

	Returns
	-------
	ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP`
		The smoothed CCMAP.

	"""

	smooth_map = ccMapObj.copy()
	smooth_map.make_editable()
	cmh.make_smooth_map(smooth_map.matrix, smooth_map.shape, filter_pass=filter_pass)
	smooth_map.make_unreadable()
	return smooth_map


def LogOfMatrix(ccMapObj):

	ccMapObj.make_readable()

	LogHiCmap = CCMAP()
	LogHiCmap.path2matrix = os.getcwd() + '/nparray_' + getRandomName() + '.bin'

	LogHiCmap.shape = ccMapObj.shape
	LogHiCmap.xticks = ccMapObj.xticks
	LogHiCmap.yticks = ccMapObj.yticks
	LogHiCmap.binsize = ccMapObj.binsize
	LogHiCmap.bLog = True

	bNonZeros = None
	#if ccMapObj.bNoData is not None:
	#	LogHiCmap.bNoData = ccMapObj.bNoData
	#	bNonZeros = ~LogHiCmap.bNoData
	#else:
	LogHiCmap.bNoData = np.all( ccMapObj.matrix == 0.0, axis=0)
	bNonZeros = ~LogHiCmap.bNoData

	# Log of part of matrix containing data
	path2matrixA = os.getcwd() + '/nparray_' + getRandomName() + '.bin'
	A = (ccMapObj.matrix[bNonZeros,:])[:,bNonZeros]   # Selected row-column which are not all zeros
	BinMatrixA = np.memmap(path2matrixA, dtype=dtype_npBINarray, mode='w+', shape=A.shape)
	BinMatrixA[:] = np.log10(A)[:]
	BinMatrixA.flush()

	# Assigning minvalue and maxvalue
	LogHiCmap.maxvalue = float(np.amax(BinMatrixA))
	minvalue = np.amin(BinMatrixA)
	v_steps = np.linspace(minvalue, LogHiCmap.maxvalue, 100)
	LogHiCmap.minvalue = minvalue - (v_steps[1] - v_steps[0])

	# Making full matrix
	BinLogMatrix = np.memmap(LogHiCmap.path2matrix, dtype=dtype_npBINarray, mode='w+', shape=LogHiCmap.shape)
	A_i = -1
	A_j = 0
	for i in range(BinLogMatrix.shape[0]):
		if not LogHiCmap.bNoData[i]:
			A_i += 1

		A_j = 0
		for j in range(BinLogMatrix.shape[1]):
			if LogHiCmap.bNoData[i] or LogHiCmap.bNoData[j]:
				BinLogMatrix[i][j] = LogHiCmap.minvalue
			else:
				BinLogMatrix[i][j] = BinMatrixA[A_i][A_j]
				A_j += 1
	BinLogMatrix.flush()

	del BinLogMatrix
	del BinMatrixA

	try:
		os.remove(path2matrixA)
	except:
		pass

	return LogHiCmap


def gen_HiC_svg_plot(ccMapObj, outfile):

	import matplotlib
	matplotlib.use('SVG')
	import matplotlib.pyplot as plt

	# File namings
	outdir = os.path.dirname( os.path.abspath(outfile) )
	outfilebase = os.path.splitext(os.path.basename(os.path.abspath(outfile)))[0]
	####

	def make_plot(idx, jdx, mat, i, j, minvalue, maxvalue, idiv=2000, jdiv=10000):
		fig = plt.figure(111)

		plt.imshow(mat[i:i+idiv, j:j+jdiv], origin='lower', aspect='equal', cmap='Greys', interpolation='none', vmax=maxvalue, vmin=minvalue)

		ax = plt.gca()

		plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)

		name = os.path.join(outdir, '{0}.{1}.{2}.{3}.{4}.svg' .format(outfilebase, idx, jdx, idiv, jdiv))

		plt.savefig(name)
		plt.clf()
		del fig

	def combine_svg(data):
		fout = open(outfile, 'w')
		fout.write('<?xml version="1.0" encoding="utf-8" standalone="no"?>\n')
		fout.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')
		fout.write('<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"> \n')
		fout.write('<g id="figure_1">\n')

		y = 0
		idx = 1
		for i in range(len(data)-1, -1, -1):
			x = 0
			tmp_y = 0
			for j in range(len(data[i])):
				filename = os.path.join( outdir, '{0}.{1}.{2}.{3}.{4}.svg' .format(outfilebase, i, j, data[i][j][0], data[i][j][1]))
				tmp_y = data[i][j][0]

				fin = open(filename, 'r')
				read_png = False
				png_data = ''
				for line in fin:
					if re.search('iVBORw0', line):
						png_data = line.split('"')[0]
						break
				fin.close()

				fout.write('<image x="{0}" xlink:href="data:image{1}/png;base64, \n' .format(x, idx))
				fout.write(png_data)
				fout.write('" y="{0}" transform="scale(0.1)" />\n' .format(y))

				x += data[i][j][1]

				os.remove(filename)

			y += tmp_y

		fout.write('</g>\n')
		fout.write('</svg>\n')
		fout.close()

	# Main part of function start here
	ccMapObj.make_readable()

	# Minimum and maximum value
	minvalue = ccMapObj.minvalue
	maxvalue = ccMapObj.maxvalue

	print (minvalue, maxvalue)

	idiv = 1000
	jdiv = 10000
	i = 0
	idx = 0
	file_data = []

	print("Generating plot segements as svg files ...\n")
	while(1):
		outer_loop_finished = False
		j = 0
		jdx = 0
		tmp_file_data = []
		while(1):
			if i+idiv >= ccMapObj.shape[0] and j+jdiv >=ccMapObj.shape[0] :
				make_plot(idx, jdx, ccMapObj.matrix, i, j, minvalue, maxvalue, idiv=ccMapObj.shape[0]-i, jdiv=ccMapObj.shape[0]-j)
				outer_loop_finished = True
				tmp_file_data.append([ccMapObj.shape[0]-i, ccMapObj.shape[0]-j])
				break

			elif i+idiv >= ccMapObj.shape[0] and j+jdiv < ccMapObj.shape[0] :
				make_plot(idx, jdx, ccMapObj.matrix, i, j, minvalue, maxvalue, idiv=ccMapObj.shape[0]-i, jdiv=jdiv)
				tmp_file_data.append([ccMapObj.shape[0]-i, jdiv])

			elif i+idiv < ccMapObj.shape[0] and j+jdiv >= ccMapObj.shape[0] :
				make_plot(idx, jdx, ccMapObj.matrix, i, j, minvalue, maxvalue, idiv=idiv, jdiv=ccMapObj.shape[0]-j)
				tmp_file_data.append([idiv, ccMapObj.shape[0]-j])
				break

			else:
				make_plot(idx, jdx, ccMapObj.matrix, i, j, minvalue, maxvalue, idiv=idiv, jdiv=jdiv)
				tmp_file_data.append([idiv, jdiv])

			j += jdiv
			jdx += 1

		file_data.append(tmp_file_data)

		if outer_loop_finished:
			break

		i += idiv
		idx += 1

	print("Merging plot segments in one svg file...\n")
	combine_svg(file_data)

def gen_backup(filename):
	fc = 1
	full_backup_path = ''
	indir = os.path.dirname( filename )
	basename = os.path.basename( filename )
	while 1:
		backup_filename = '#' + basename + '.{0}#' .format(fc)
		full_backup_path = os.path.join(indir,  backup_filename)
		if not os.path.exists(full_backup_path):
			break
		fc += 1

	shutil.move(filename, full_backup_path)
	logger.info('File {0} already exists !!! Backing up as {1} .' .format(filename, full_backup_path))

	return
