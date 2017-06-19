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
from scipy import stats as spstats
import re
import os
import copy
import shlex
import shutil
import subprocess
import tempfile
import logging
import h5py
import json
import csv
from urllib.request import urlopen, Request
from urllib.error import URLError

from gcMapExplorer.config import getConfig, updateConfig

from . import util

# get configuration
config = getConfig()

def downloadFile(url, output):
    """ Download a file from url and save to a file

    Parameters
    ----------
    ouput : str
        Output file name.

    Returns
    -------
    sucess : bool
        If downloaded successfully ``True`` otherwise ``False``
    """
    try:
        response = urlopen(url)
    except (ValueError, URLError):
        return False

    fout = open(output, 'wb')
    shutil.copyfileobj(response, fout)
    return True


def check_resolution_list(resolutions):
    new_resolutions = [ '1kb', '2kb', '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb', '100kb', '160kb', '200kb', '320kb', '500kb', '640kb', '1mb']

    if resolutions is not None:
        if isinstance(resolutions, list):
            new_resolutions = list( set(new_resolutions + resolutions) )
        else:
            raise TypeError(' Resolutions should be of list type.')

    return new_resolutions

def check_coarsening_method(methods):
    acepted_methods = ['min', 'max', 'amean', 'hmean', 'gmean', 'median']
    if methods is not None:
        for method in methods:
            if method not in acepted_methods:
                raise ValueError( ' Coarsening method {0} is not implemented..\
                \n Use these: {1}'.format(method, acepted_methods) )

        return methods
    else:
        return acepted_methods

class TempNumpyArrayFiles:
    """To handle temporary numpy array files

    To convert a Wig file to hdf5 file, data are parsed, and further stored temporarily in these memory-mapped numpy array files.
    Use of numpy arrays avoid dependency from storing chromosome location/coordinates because array index is used as the location/coordinates.
    Additionaly, chromosome could be very large and to store these arrays could be memory expensive, these arrays are stored as binary files on the disk.

    .. note::
        These generated files are either automatically deleted after execution of script or by deleting [``del``] :class:`TempNumpyArrayFiles` instance.

    Attributes
    ----------
    chromSizeInfo : dict
        Dictionary contains chromosome size. Numpy array files will be generated on the basis of these sizes.
    files : dict
        Dictionary for names of temporary numpy array files, where keys are chromosomes and values are respective file names.
    arrays : dict
        Dictionary of memory mapped numpy array (`numpy.memmap <http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.memmap.html>`_), where
        keys are chromosomes and values are respective `numpy.memmap <http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.memmap.html>`_ arrays.
    workDir : str
        Working directory where temporary numpy array files will be generated.



    """
    def __init__(self, workDir=None):
        self.chromSizeInfo = None
        self.files = None
        self.arrays = None
        self.printInfo = dict()

        # Working and output directory
        if workDir is None:
            self.workDir = config['Dirs']['WorkingDirectory']
        else:
            self.workDir = workDir

        self.logger = logging.getLogger(__name__+'.TempNumpyArrayFiles')
        self.logger.setLevel(logging.INFO)

    def __del__(self):
        self._removeTempNumpyFiles()

    def _removeTempNumpyFiles(self):
        if self.arrays is not None:
            for key in self.arrays:
                self._removeTempNumpyFile(key)

    def _removeTempNumpyFile(self, key):
        self.logger.info(' Closing and removing temporary numpy array file [{0}] .' .format(self.files[key]))
        self.arrays[key]._mmap.close()

        if os.path.isfile(self.files[key]):
            os.remove(self.files[key])

    def updateArraysByBigWig(self, bigWigFileName):
        """Update/resize all array files using given bigWig file

        Parameters
        ----------
        bigWigFileName : str
            Name of input bigWig file


        """
        chromSizeInfo = self._getBigWigInfo(bigWigFileName)
        for key in chromSizeInfo:
            if key in self.chromSizeInfo:
                if chromSizeInfo[key] > self.chromSizeInfo[key]:
                    self.chromSizeInfo[key] = chromSizeInfo[key]
                    self.generateTempNumpyFile(key, regenerate=True)
            else:
                self.chromSizeInfo[key] = chromSizeInfo[key]
                self.generateTempNumpyFile(key)

    def updateArraysByChromSize(self, chrom, size):
        """Update/resize an array file using given chromosome and its size

        Parameters
        ----------
        chrom : str
            Chromosome name
        size : int
            Total size of chromosome

        """
        if chrom in self.chromSizeInfo:
            if size > self.chromSizeInfo[chrom]:
                self.chromSizeInfo[chrom] = size
                self.generateTempNumpyFile(chrom, regenerate=True)
        else:
            self.chromSizeInfo[chrom] = size
            self.generateTempNumpyFile(chrom)

    def addChromSizeInfo(self, bigWigFileName):
        """ Update chromosome sizes using new bigWig file

        To update :attr:`TempNumpyArrayFiles.chromSizeInfo` for new bigWig files, this method can be used.

        .. note::
            This method only updates the :attr:`TempNumpyArrayFiles.chromSizeInfo` dictionary. It does not resize numpy array files.

        Parameters
        ----------
        bigWigFileName : str
            Name of input bigWig file

        """
        chromSizeInfo = self._getBigWigInfo(bigWigFileName)

        for key in chromSizeInfo:
            if key in self.chromSizeInfo:
                self.chromSizeInfo[key] = max( chromSizeInfo[key], self.chromSizeInfo[key] )
            else:
                self.chromSizeInfo[key] = chromSizeInfo[key]

    def _getBigWigInfo(self, filename):
        """Base method to Retrieve chromosome names and their sizes

        * Use :meth:`TempNumpyArrayFiles.addChromSizeInfo` to automatically retrieve chromosome size information from bigWig file and to store in :attr:`TempNumpyArrayFiles.chromSizeInfo`.

        .. warning::
            **Private method.** Use it at your own risk. It is used internally in :meth:`TempNumpyArrayFiles.addChromSizeInfo`

        Parameters
        ----------
        filename : str
            Input bigWig file

        """

        cmd = '{0} -chroms {1}' .format(self.pathTobigWigInfo, filename)

        args = shlex.split(cmd)
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        self.logger.info(' Started running command: [{0}] ...' .format(cmd))
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            self.logger.info(stdout)
            self.logger.info(stderr)
            self.logger.error(' Not able to run command: [ {0} ]. See details above.' .format(cmd))
        self.logger.info('         ... Finished!!!' .format(cmd))

        chromSizeInfo = dict()

        chroms = None
        read = False
        counter = 0
        for line in re.split('\n', stdout.decode("utf-8")):
            if re.search('chromCount', line) != None:
                temp = re.split('\s+', line)
                chroms = int(temp[1])
                read = True
                continue

            if read:
                line = line.lstrip().rstrip()
                temp = re.split('\s+', line)
                chromSizeInfo[temp[0]] = int(temp[2])
                counter += 1

                if counter >= chroms:
                    read = False

        return chromSizeInfo

    def genrateAllTempNumpyFiles(self):
        """Generate all memory mapped numpy array files

        It is used to generate all memory mapped numpy array files using the :attr:`TempNumpyArrayFiles.chromSizeInfo` dictionary.

        """

        if self.self.chromSizeInfo is None:
            raise ValueError(" Chromosome sizes are not available. ")
        if self.files is None:
            self.files = dict()

        self.logger.info(' Generating temporary numpy array files...')
        for key in self.chromSizeInfo:
            self.generateTempNumpyFile(key)

    def generateTempNumpyFile(self, key, regenerate=False):
        """Generate a memory mapped numpy array file

        It is used to generate a memory mapped numpy array for given chromosome for which size is already the :attr:`TempNumpyArrayFiles.chromSizeInfo` dictionary.

        Parameters
        ----------
        key : str
            chromosome name. Should be present as key in :attr:`TempNumpyArrayFiles.chromSizeInfo`.
        regenerate : bool
            Replace or regenerate new memory mapped array for given chromosome

        """

        # Initialize print information
        for chrom in self.chromSizeInfo:
            if chrom not in self.printInfo:
                self.printInfo[chrom] = True

        try:
            self._generateTempNumpyFile(key, regenerate=regenerate)
        except (SystemExit, KeyboardInterrupt) as e:
            self._removeTempNumpyFiles()
            raise e

    def _generateTempNumpyFile(self, key, regenerate=False):
        """enerate a memory mapped numpy array file

        It is used to generate a memory mapped numpy array file for given chromosome for which size is already in the :attr:`TempNumpyArrayFiles.chromSizeInfo` dictionary.

        .. warning::
            **Private method.** Use it at your own risk. It is used internally in :meth:`TempNumpyArrayFiles.generateTempNumpyFile`.

        Parameters
        ----------
        key : str
            chromosome name. Should be present as key in :attr:`TempNumpyArrayFiles.chromSizeInfo`.
        regenerate : bool
            Replace or regenerate new memory mapped array for given chromosome

        """
        if self.files is None:
            self.files = dict()

        if self.arrays is None:
            self.arrays = dict()

        if (key in self.files) and not regenerate:
            if self.printInfo[key]:
                self.logger.info(' {0} is already present in temporary numpy file list, skipping ...' .format(key))
                self.printInfo[key] = False
        else:
            if regenerate:
                self._removeTempNumpyFile(key)

            (fd, fname) = tempfile.mkstemp(suffix='.npy', prefix='gcx_'+key+'_', dir=self.workDir, text=False)
            os.close(fd)     # Close file, error in windows OS
            self.logger.info(' Generating temporary numpy array file [{0}] for {1} ...' .format(fname, key))

            size = self.chromSizeInfo[key] + 1        # Be careful: added one to easiliy handle real locations, zeroth index is dummy, dont use zeroth location
            np.save(fname, np.zeros(size, dtype=np.float))
            npy = np.load(fname, mmap_mode='r+')

            self.arrays[key] = npy
            self.files[key] = fname

            # print(npy.shape)
            self.printInfo[key] = True
            self.logger.info(' \t... Finished generating temporary numpy array file.')

    def fillAllArraysWithZeros(self):
        """Fill all arrays with zeros.

        To fill all memory mapped array with zero. It is used in :class:`WigHandler` so that new data extracted from Wig files can be stored in these array files.

        """
        if self.arrays is not None:
            for key in self.chromSizeInfo:
                self.arrays[key].fill(0.0)

class HDF5Handler:
    """Handler for genomic data HDF5 file.

    This class acts like a handler and can be used to read, write, and modify genomic data file in `HDF5 format <https://www.hdfgroup.org/why_hdf/>`_.
    This is a binary file and is compressed using ``zlib`` method to reduce the storage memory.


    Structure of HDF5 file: ``/<Chromosome>/<Resolution>/<1D Numpy Array>``

        ::

          HDF5 ──────────────────────────> title
            ├──────── chr1
            │           ├───── 1kb
            │           │        ├──────── amean  ( Arithmatic mean) (type: 1D Numpy Array)
            │           │        ├──────── median ( Median value   ) (type: 1D Numpy Array)
            │           │        ├──────── hmean  ( Harmonic mean  ) (type: 1D Numpy Array)
            │           │        ├──────── gmean  ( Geometric mean ) (type: 1D Numpy Array)
            │           │        ├──────── min    ( Minimum value  ) (type: 1D Numpy Array)
            │           │        └──────── max    ( Maximum value  ) (type: 1D Numpy Array)
            │           │
            │           ├────  5kb
            │           │        ├──────── amean  ( Arithmatic mean) (type: 1D Numpy Array)
            │           │        ├──────── median ( Median value   ) (type: 1D Numpy Array)
            │           │        ├──────── hmean  ( Harmonic mean  ) (type: 1D Numpy Array)
            │           │        ├──────── gmean  ( Geometric mean ) (type: 1D Numpy Array)
            │           │        ├──────── min    ( Minimum value  ) (type: 1D Numpy Array)
            │           │        └──────── max    ( Maximum value  ) (type: 1D Numpy Array)
            │           │
            │           └────  ...
            │
            ├──────── chr2
            │           ├───── 1kb
            │           │        ├──────── amean  ( Arithmatic mean) (type: 1D Numpy Array)
            │           │        ├──────── median ( Median value   ) (type: 1D Numpy Array)
            │           │        ├──────── hmean  ( Harmonic mean  ) (type: 1D Numpy Array)
            │           │        ├──────── gmean  ( Geometric mean ) (type: 1D Numpy Array)
            │           │        ├──────── min    ( Minimum value  ) (type: 1D Numpy Array)
            │           │        └──────── max    ( Maximum value  ) (type: 1D Numpy Array)
            │           └────  ..
            :
            :
            :
            └───── ...


    Attributes
    ----------
    filename : str
        HDF5 file name

    title : str
        Title of the data

    hdf5 : h5py.File
        input/output stream to HDF5 file

    data : dict
        This dictionary is generated by :meth:`HDF5Handler.buildDataTree()`. This dictionary gives access to all data arrays.



    Parameters
    ----------
    filename : str
        HDF5 file name. e.g.: ``abcxyz.h5``
    title : str
        title or name of the data



    Examples
    --------
    .. code-block:: python

        from hiCMapAnalyze import genomicsDataHandler as gdh
        import numpy as np

        # Load available file
        hdf5Hand = gdh.HDF5Handler('test.h5')

        # Build the data structure, this is an essential step to retrieve the data easiliy as shown below
        hdf5Hand.buildDataTree()

        # Print shape and maximum value of chr1->1kb->mean array
        print(hdf5Hand.data['chr1']['1kb']['mean'].shape, np.amax(hdf5Hand.data['chr1']['1kb']['mean']))


    """

    def __init__(self, filename, title=None):
        self.logger = logging.getLogger(__name__+'.HDF5Handler')
        self.logger.setLevel(logging.INFO)

        self.filename = filename
        self.hdf5 = None

        self.title = None
        if title is not None:
            self.title = title

        self.data = None

    def __del__(self):
        if self.hdf5 is not None:
            self.close()

    def close(self):
        """ close hdf5 file
        """
        self.hdf5.close()
        self.logger.info(' Closed {0} ...' .format(self.filename))
        self.hdf5 = None

    def open(self):
        """ open hdf5 file
        """
        self.hdf5 = h5py.File(self.filename)

        # Fail-safe mechanism for title in case of either new file or already opened file
        if 'title' not in self.hdf5.attrs:
            if self.title is None:
                self.hdf5.attrs['title'] = os.path.splitext(os.path.basename(self.filename))[0]
                self.title = os.path.splitext(os.path.basename(self.filename))[0]
            else:
                self.hdf5.attrs['title'] = self.title
        else:
            if self.title is None:
                self.title = self.hdf5.attrs['title']

        self.logger.info(' Opened {0} ...' .format(self.filename))

    def setTitle(self, title):
        """ Set title of the dataset

        It can be used to set or replace the title of the dataset.
        If file is not yet opened, title will be stored to file when file
        will be opened.

        Parameters
        ----------
        title : str
            The title of dataset.

        """
        self.title = title
        if self.hdf5 is not None:
            self.hdf5.attrs['title'] = os.path.splitext(os.path.basename(self.filename))[0]
        else:
            self.logger.info(" File is not opened!!! Title will be stored to file on file opening.")

    def getChromList(self):
        """To get list of all chromosomes present in hdf5 file

        Returns
        -------
        chroms : list
            List of all chromosomes present in hdf5 file

        """
        if self.hdf5 is None:
            self.open()

        chroms = []
        for g in self.hdf5.keys():
            chroms.append(g)

        return chroms

    def hasChromosome(self, chrom):
        """To get list of all chromosomes present in hdf5 file

        Parameters
        ----------
        chrom : str
            Chromosome name to be look up in file.

        Returns
        -------
        gotChromsome : bool
            If queried chromosome present in file ``True`` otherwise ``False``.

        """

        if self.hdf5 is None:
            self.open()

        gotChromsome = False
        if chrom in self.hdf5:
            gotChromsome = True

        return gotChromsome


    def getResolutionList(self, chrom, dataName=None):
        """ To get all resolutions for given chromosome from hdf5 file

        Parameters
        ----------
        chrom : str
            chromosome name
        dataName : str
            Options to get list of all resolution list for given data name

        Returns
        -------
        resolutionList : list[str]
            A list of all available resolutions for the given chromosome

        Raises
        ------
        KeyError
            If chromosome not found in hdf5 file

        """
        chromList = self.getChromList()

        resolutionList = []
        if chrom in chromList:
            for r in self.hdf5[chrom].keys():
                if dataName is not None and self.hasDataName(chrom, r, dataName):
                    resolutionList.append(r)
                else:
                    resolutionList.append(r)
        else:
            raise KeyError(' Chromosome [{0}] not found in [{1}] file...' .format(chrom, self.filename))

        return resolutionList

    def hasResolution(self, chrom, resolution, dataName=None):
        """To check if given resolution for given chromosome is present

        Parameters
        ----------
        chrom : str
            Chromosome name to be look up in file.
        resolution : str
            Data Resolution for queried Chromosome
        dataName : str
            Options to check if resolution for given data name is present

        Returns
        -------
        gotResolution : bool
            If queried resolution of given chromsome present in file ``True``
            otherwise ``False``.

        """

        if self.hdf5 is None:
            self.open()

        gotResolution = False
        if self.hasChromosome(chrom):
            if resolution in self.getResolutionList(chrom, dataName=dataName):
                gotResolution = True

        return gotResolution

    def getDataNameList(self, chrom, resolution):
        """ List of all available arrays by respecitve coarse method name for given chromosome and resolution

        Parameters
        ----------
        chrom : str
            chromosome name
        resolution : str
            resolution

        Returns
        -------
        nameList : list[str]
            List of arrays by name of dataset

        Raises
        ------
        KeyError
            If chromosome not found in hdf5 file. If input resolution keyword is not found for input chromosome.

        """

        nameList = []

        if chrom in self.getChromList():
            if resolution in self.getResolutionList(chrom):
                for m in self.hdf5[chrom][resolution].keys():
                    nameList.append(m)
            else:
                raise KeyError(' Resolution [{0}] for Chromosome [{1}] not found in [{2}] file...' .format(resolution, chrom, self.filename))

        else:
            raise KeyError(' Chromosome [{0}] not found in [{1}] file...' .format(chrom, self.filename))

        return nameList

    def hasDataName(self, chrom, resolution, dataName):
        """To check if given data for given resolution for given chromosome is present

        Parameters
        ----------
        chrom : str
            Chromosome name to be look up in file.
        resolution : str
            Data Resolution for queried Chromosome
        dataName : str
            Name of data to be queried in given Chromosome.

        Returns
        -------
        gotDataName : bool
            If queried data in given chromsome at given resolution is
            present in file ``True`` otherwise ``False``.

        """

        if self.hdf5 is None:
            self.open()

        gotDataName = False
        if self.hasResolution(chrom, resolution):
            if dataName in self.hdf5[chrom][resolution]:
                gotDataName = True

        return gotDataName

    def buildDataTree(self):
        """ Build data dictionary from the input hdf5 file

        To retrieve the data from hdf5 file, this function should be used to built the dictionary :attr:`HDF5Handler.data`.
        This dictionary gives access directly to data of any chromosome with specific resolution.

        """

        if self.data is None:
            self.data = dict()
        for chrom in self.getChromList():
            self.data[chrom] = dict()
            for resolution in self.getResolutionList(chrom):
                self.data[chrom][resolution] = dict()
                for cm in self.getDataNameList(chrom, resolution):
                    self.data[chrom][resolution][cm] = self.hdf5[chrom][resolution][cm]

    def addDataByArray(self, Chrom, resolution, data_name, value_array, compression='lzf'):
        """ Add array to the hdf5 file for given chromosome, resolution and data name.
            It can be used either to add new data array or to replace existing data.

        Parameters
        ----------
        Chrom : str
            Chromosome Name
        resolution : str
            Reslution of data
        data_name : str
            Name of data.
        value_array : numpy.ndarray
            An array containing values.

        """

        if self.hdf5 is None:
            self.open()

        firstLevel = '/' + Chrom
        secondLevel = '/' + Chrom + '/' + resolution
        thirdLevel = '/' + Chrom + '/' + resolution + '/' + data_name

        if Chrom not in self.hdf5:
            self.hdf5.create_group(Chrom)

        if resolution not in self.hdf5[Chrom]:
            self.hdf5[Chrom].create_group(resolution)

        array = np.asarray(value_array)
        if data_name not in self.hdf5[Chrom][resolution]:
            if compression == 'gzip':
                self.hdf5[Chrom][resolution].create_dataset(data_name, array.shape, dtype=array.dtype, data=array, compression="gzip", shuffle=True, compression_opts=4)
            else:
                self.hdf5[Chrom][resolution].create_dataset(data_name, array.shape, dtype=array.dtype, data=array, compression="lzf", shuffle=True)
        else:
            self.hdf5[Chrom][resolution].pop(data_name)
            if compression == 'gzip':
                self.hdf5[Chrom][resolution].create_dataset(data_name, array.shape, dtype=array.dtype, data=array, compression="gzip", shuffle=True, compression_opts=4)
            else:
                self.hdf5[Chrom][resolution].create_dataset(data_name, array.shape, dtype=array.dtype, data=array, compression="lzf", shuffle=True)

        self.hdf5[Chrom][resolution][data_name].attrs['compression'] = compression


class BigWigHandler:
    """ To handle bigWig files and to convert it to h5 file

    This class can be used to convert bigWig file to h5 file.
    It can also be used to combine several bigWig files that are originated from replicated experiments.

    .. warning:: Presently ``bigWigToWig`` and ``bigWigInfo`` is not available
                 for Windows OS. Therefore, this class will fail in this OS.

    Attributes
    ----------
    bigWigFileNames : str or list[str]
        List of bigWig file names including path
    pathTobigWigToWig : str
        Path to ``bigWigToWig`` program. It can be downloaded from
        http://hgdownload.cse.ucsc.edu/admin/exe/ for MacOSX and Linux.
        If path to program is already present in configuration file, it will be
        taken from the configuration.

        If it is not present in configuration file, the input path **should**
        be provided. It will be stored in configuration file for later use.
    pathTobigWigInfo : str
        Path to ``bigWigInfo`` program. It can be downloaded from
        http://hgdownload.cse.ucsc.edu/admin/exe/ for MacOSX and Linux.
        If path to program is already present in configuration file, it will be
        taken from the configuration.

        If it is not present in configuration file, the input path **should**
        be provided. It will be stored in configuration file for later use.
    WigFileNames : str
        List of Wig file names, either autmoatically generated or given by user
    chromName : str
        Name of input target chromosome. If this is provided, only this chromosome
        data is extracted and stored in h5 file.
    wigHandle : WigHandler
        WigHandler instance to parse Wig file and save data as hdf5 file
    chromSizeInfo : dict
        A dictionary containing chromosome size information
    methodToCombine : str
        method to combine bigWig/Wig files, Presently, accepted keywords are: ``mean``, ``min`` and ``max``
    maxEntryWrite : int
        Number of lines read from Wig file at an instant, after this, data is dumped in temporary numpy array file

    Parameters
    ----------
    filenames : str or list[str]
        A bigWig file or list of bigWig files including path
    pathTobigWigToWig: str
        Path to ``bigWigToWig`` program. It can be downloaded from
        http://hgdownload.cse.ucsc.edu/admin/exe/ for MacOSX and Linux.
        If path to program is already present in configuration file, it will be
        taken from the configuration.

        If it is not present in configuration file, the input path **should**
        be provided. It will be stored in configuration file for later use.
    pathTobigWigInfo : str
        Path to ``bigWigInfo`` program. It can be downloaded from
        http://hgdownload.cse.ucsc.edu/admin/exe/ for MacOSX and Linux.
        If path to program is already present in configuration file, it will be
        taken from the configuration.

        If it is not present in configuration file, the input path **should**
        be provided. It will be stored in configuration file for later use.
    chromName : str
        Name of input target chromosome. If this is provided, only this chromosome
        data is extracted and stored in h5 file.
    methodToCombine : str
        method to combine bigWig/Wig files, Presently, accepted keywords are: ``mean``, ``min`` and ``max``
    maxEntryWrite : int
        Number of lines read from Wig file at an instant, after this, data is dumped in temporary numpy array file. To reduce memory (RAM) occupancy,
        reduce this number because large numbers need large RAM.

    """

    def __init__(self, filenames, pathTobigWigToWig=None, pathTobigWigInfo=None, chromName=None, methodToCombine='mean', workDir=None, maxEntryWrite=10000000):
        self.bigWigFileNames = filenames
        self.pathTobigWigToWig = pathTobigWigToWig
        self.pathTobigWigInfo = pathTobigWigInfo
        self.chromName = chromName            # If only an input chromosome is required

        self.WigFileNames = None           # Wig file name
        self.wigToDel = False             # Flag to delete or retain wig file
        self.wigHandle = None             # Used for WigHandler object

        self.chromSizeInfo = None         # dictionary for chromosome sizes, extracted from bigWigInfo program

        # Working and output directory
        if workDir is None:
            self.workDir = config['Dirs']['WorkingDirectory']
        else:
            self.workDir = workDir

        self.logger = logging.getLogger(__name__+'.BigWigHandler')
        self.logger.setLevel(logging.INFO)

        self.maxEntryWrite = maxEntryWrite # maximum count before writing to output file will be attempted

        # Check for programs
        self._checkBigWigInfoProgram(pathTobigWigInfo)
        self._checkBigWigToWigProgram(pathTobigWigToWig)

        # Convert to list
        if not isinstance(filenames, list):
            self.bigWigFileNames = [ filenames ]

        if methodToCombine not in ['mean', 'max', 'min']:
            raise NotImplementedError(' Method [{0}] to combine bigwig file not implemented. Use: \'mean\', \'max\' or \'min\' ' .format(methodToCombine))
        self.methodToCombine = methodToCombine

    def __del__(self):
        self._removeTempWigFiles()

    def _checkBigWigInfoProgram(self, pathTobigWigInfo):
        """ Check if bigWigInfo program is available or accessible.

        If program is not available in configuration file, the given
        path will be stored in the file after checking its accessibility.

        The path is stored in :attr:`gcMapExplorer.lib.genomicsDataHandler.BigWigHandler.pathTobigWigInfo`

        Parameters
        ----------
        pathTobigWigInfo : str
            Path to bigWigInfo program

        """
        if pathTobigWigInfo is None or pathTobigWigInfo == config['Programs']['bigWigInfo']:
            self.pathTobigWigInfo = config['Programs']['bigWigInfo']
            if self.pathTobigWigInfo == 'None':
                raise ValueError('Path to bigwiginfo is not set/provided...')

            if not os.path.isfile(self.pathTobigWigInfo):
                raise IOError('bigwiginfo is not found at {0} . Set or provide different path.'.format(self.pathTobigWigInfo))
        else:
            if not os.path.isfile(pathTobigWigInfo):
                raise IOError('bigwiginfo is not found at {0} . Set or provide different path.'.format(pathTobigWigInfo))

            # Add to config file
            if pathTobigWigInfo != config['Programs']['bigWigInfo']:
                self.logger.info('Adding bigwiginfo path in configuration...')
                updateConfig('Programs', 'bigWigInfo', self.pathTobigWigInfo)
                config['Programs']['bigWigInfo'] = self.pathTobigWigInfo

            self.pathTobigWigInfo = pathTobigWigInfo

    def _checkBigWigToWigProgram(self, pathTobigWigToWig):
        """ Check if bigWigToWig program is available or accessible.

        If program is not available in configuration file, the given
        path will be stored in the file after checking its accessibility.

        The path is stored in :attr:`gcMapExplorer.lib.genomicsDataHandler.BigWigHandler.pathTobigWigToWig`

        Parameters
        ----------
        pathTobigWigToWig : str
            Path to bigWigToWig program

        """
        if pathTobigWigToWig is None or pathTobigWigToWig == config['Programs']['bigWigToWig']:
            self.pathTobigWigToWig = config['Programs']['bigWigToWig']
            if self.pathTobigWigToWig == 'None':
                raise ValueError('Path to bigWigToWig is not set/provided...')

            if not os.path.isfile(self.pathTobigWigToWig):
                raise IOError('bigWigToWig is not found at {0} . Set or provide different path.'.format(self.pathTobigWigToWig))
        else:
            if not os.path.isfile(pathTobigWigToWig):
                raise IOError('bigWigToWig is not found at {0} . Set or provide different path.'.format(pathTobigWigToWig))

            # Add to config file
            if pathTobigWigToWig != config['Programs']['bigWigToWig']:
                self.logger.info('Adding bigWigToWig path in configuration...')
                updateConfig('Programs', 'bigWigToWig', self.pathTobigWigToWig)
                config['Programs']['bigWigToWig'] = self.pathTobigWigToWig

            self.pathTobigWigToWig = pathTobigWigToWig

    def _removeTempWigFiles(self):
        if self.wigToDel:
            for f in self.WigFileNames:
                self.logger.info(' Removing temporary Wig file [{0}] .' .format(f))
                if os.path.isfile(f):
                    os.remove(f)

    def _getBigWigInfo(self, filename):
        """Base method to Retrieve chromosome names and their sizes

        * Chromosome size information is stored for a given bigWig file.
          If size of chromosome is already present in dictionary, largest size
          is stored in dictionary.

        * Use :meth:`BigWigHandler.getBigWigInfo` to automatically retrieve
          chromosome size information from all bigWig files.

        .. warning::
            **Private method.** Use it at your own risk. It is used internally
            in :meth:`BigWigHandler.getBigWigInfo`

        Parameters
        ----------
        filename : str
            Input bigWig file

        """
        cmd = '{0} -chroms {1}' .format(self.pathTobigWigInfo, filename)

        args = shlex.split(cmd)
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        self.logger.info(' Started running command: [{0}] ...' .format(cmd))
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            self.logger.info(stdout)
            self.logger.info(stderr)
            self.logger.error(' Not able to run command: [ {0} ]. See details above.' .format(cmd))
        self.logger.info('         ... Finished!!!' .format(cmd))

        if self.chromSizeInfo is None:
            self.chromSizeInfo = dict()

        chroms = None
        read = False
        counter = 0
        for line in re.split('\n', stdout.decode("utf-8")):
            if re.search('chromCount', line) != None:
                temp = re.split('\s+', line)
                chroms = int(temp[1])
                read = True
                continue

            if read:
                line = line.lstrip().rstrip()
                temp = re.split('\s+', line)

                if temp[0] in self.chromSizeInfo:
                    self.chromSizeInfo[temp[0]] = max( int(temp[2]), self.chromSizeInfo[temp[0]] )
                else:
                    self.chromSizeInfo[temp[0]] = int(temp[2])
                counter += 1

                if counter >= chroms:
                    read = False

    def getBigWigInfo(self):
        """Retrieve chromosome names and their sizes

        BigWigInfo progrom is executed on all listed bigWig files and
        chromosomes name with respecitve size is stored in
        :attr:`BigWigHandler.chromSizeInfo` variable.
        From the several listed bigWig files, only largest size of chromosomes
        are considered.

        If :attr:`BigWigHandler.chromName` is provided, only target chromsome
        information is kept in :attr:`BigWigHandler.chromSizeInfo` dictionary.

        """
        for f in self.bigWigFileNames:
            self._getBigWigInfo(f)

        # In case if a single chromosome is required to be extracted
        if self.chromName is not None:
            chromSizeInfo = dict()
            if self.chromName not in self.chromSizeInfo:
                self.logger.error(' The requested chromosome {0} is not found in bigWig files!!.'.format(self.chromName))

            chromSizeInfo[self.chromName] = self.chromSizeInfo[self.chromName]
            self.chromSizeInfo = chromSizeInfo

        # To print chromosome information
        output_line = ''
        chroms = 0
        for chrom in self.chromSizeInfo:
            output_line = output_line + '\n\t\t\t{0}\t\t{1}' .format(chrom, self.chromSizeInfo[chrom])
            chroms += 1
        self.logger.info(' bigWig files {0} contains {1} chromosomes. \n\t\t\tChromosome\tSize {2}' .format(self.bigWigFileNames, chroms, output_line))

    def _bigWigtoWig(self, bigWigFileName, outfilename):
        """Base method to generate Wig file from a bigWig file

        Use :meth:`BigWigHandler.bigWigtoWig` to automatically convert all
        bigWig files to Wig files.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally
            in :meth:`BigWigHandler.bigWigtoWig`


        Parameters
        ----------
        bigWigFileName : str
            Input bigWig file names.

        outfilename : str
            Name of output Wig file.

        """
        cmd = None
        if self.chromName is None:
            cmd = '{0} {1} {2}' .format(self.pathTobigWigToWig, bigWigFileName, outfilename)
        else:
            cmd = '{0} -chrom={1} {2} {3}' .format(self.pathTobigWigToWig, self.chromName, bigWigFileName, outfilename)

        args = shlex.split(cmd)
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        self.logger.info(' Started running command: [{0}] ...' .format(cmd))
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            self.logger.info(stdout)
            self.logger.info(stderr)
            self.logger.error(' Not able to convert {0} to bigWig format by command: [ {1} ]. See details above.' .format(self.bigWigFileNames, cmd))
        self.logger.info('         ... Finished!!!' .format(cmd))

    def bigWigtoWig(self, outfilenames=None):
        """ To generate Wig file

        It uses bigWigToWig program to convert bigWig to Wig file. It uses
        :attr:`BigWigHandler.chromSizeInfo` to extract the listed chromosome
        data.

        If outfilenames are provided, wig files are generated with these names.
        Otherwise, Wig file names are generated randomly and listed in
        :attr:`BigWigHandler.WigFileNames`. If these files are generated with
        random names, these will be deleted after execution.

        Parameters
        ----------
        outfilenames : str or list of strip
            List of Wig file names. If ``None``, names are automatically
            generated, files are temporarily created and after execution,
            all files are deleted.

        """
        if self.WigFileNames is None:
            self.WigFileNames = []

        # Convert to list
        if outfilenames is not None:
            if not isinstance(outfilenames, list):
                if len(self.bigWigFileNames) == 1:
                    self.WigFileNames.append(outfilenames)
                else:
                    dirname = os.path.dirname(outfilenames)
                    basename = os.path.basename(outfilenames)
                    fname, ext = os.path.splitext(basename)
                    for i in range(len(self.bigWigFileNames)):
                        outfile = os.path.join(dirname, '{0}_{1}{2}' .format(fname, i+1, ext))
                        self.WigFileNames.append(outfile)
            else:
                if len(self.bigWigFileNames) != len(outfilenames):
                    raise ValueError(' Number of input bigwig files [{0}] does not match with number of output wig files [{1}].' .format(len(self.bigWigFileNames), len(self.WigFileNames)))
                self.WigFileNames = outfilenames
        else:
            for i in range(len(self.bigWigFileNames)):
                tName = os.path.join(self.workDir, '{0}_{1}.wig.tmp' .format(util.getRandomName(), i+1) )
                self.WigFileNames.append(tName)
                self.wigToDel = True

        for i in range(len(self.bigWigFileNames)):
            self._bigWigtoWig(self.bigWigFileNames[i], self.WigFileNames[i])

    def saveAsH5(self, filename, tmpNumpyArrayFiles=None, title=None, resolutions=None, coarsening_methods=None, compression='lzf', keep_original=False):
        """Save data to h5 file.

        Parameters
        ----------
        filename : str
            Output hdf5 file name with h5 extension.
        tmpNumpyArrayFiles : :class:`TempNumpyArrayFiles` (optional)
             Usually not required. This :class:`TempNumpyArrayFiles` instance stores the temporary numpy array files information.
             To convert large number of bigWig files, its use increases the conversion speed significantly because new temporary array files takes time to generate
             and frequent generation of these files can be avoided.
        title : str (optional)
            Title of the data
        resolutions : list of str
            Additional input resolutions other than these default resolutions:
            1kb', '2kb', '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb',
            '100kb', '160kb','200kb', '320kb', '500kb', '640kb',  and '1mb'.

            For Example: use ``resolutions=['25kb', '50kb', '75kb']`` to add
            additional 25kb, 50kb and 75kb resolution data.
        coarsening_methods : list of str
            Methods to coarse or downsample the data for converting from 1-base
            to coarser resolutions. Presently, five methods are implemented.

            * ``'min'``    -> Minimum value
            * ``'max'``    -> Maximum value
            * ``'amean'``  -> Arithmatic mean or average
            * ``'hmean'``  -> Harmonic mean
            * ``'gmean'``  -> Geometric mean
            * ``'median'`` -> Median

            In case of ``None``, all five methods will be considered. User may
            use only subset of these methods. For example:
            ``coarse_method=['max', 'amean']`` can be used for downsampling by
            only these two methods.
        compression : str
            data compression method in HDF5 file : ``lzf`` or ``gzip`` method.
        keep_original : bool
            Whether original data present in bigwig file should be incorporated in HDF5 file. This will significantly increase size of HDF5 file.

        Examples
        --------
        .. code-block:: python

            from gcMapExplorer.lib import genomicsDataHandler as gdh

            # start BigWigHandler to combine and convert two bigWig files
            bigwig = gdh.BigWigHandler(['first.bigWig', 'second.bigWig'], './bigWigToWig', './bigWigInfo')

            # Save hdf5 file with two additional resolutions
            # and only two downsampling method.
            bigwig.saveAsH5('converted.h5', resolutions=['25kb', '50kb'], coarsening_methods=['max', 'amean'])


        """
        if self.chromSizeInfo is None:
            self.logger.info(' No information of chromsomes and their sizes from BigWig File!!! Need to get information...')
            self.getBigWigInfo()

        if self.WigFileNames is None:
            self.logger.info(' No Wig File!!! Need to generate Wig File...')
            self.bigWigtoWig()

        self.wigHandle = WigHandler(self.WigFileNames, self.chromSizeInfo, tmpNumpyArrayFiles=tmpNumpyArrayFiles,
                                     methodToCombine=self.methodToCombine, workDir=self.workDir, maxEntryWrite=self.maxEntryWrite)
        self.wigHandle.saveAsH5(filename, title=title, resolutions=resolutions,
                                coarsening_methods=coarsening_methods,
                                compression=compression,
                                keep_original=keep_original)


class WigHandler:
    """To convert Wig files to hdf5 file

    It parses wig files and save all data to a hdf5 file for given resolutions.

    Attributes
    ----------
    WigFileNames : list[str]
        List of input  Wig files.

        .. note:: In case if :attr:`WigHandler.chromName` is provided, only
                  one wig file is accepted.
    chromName : str
        Name of target chromosome name need to be extracted from wig file.
    chromSizeInfo : dict
        A dictionary containing chromosome size information.
    _chromPointerInFile : dict
        A dictionary containing position index of each chromosome in wig file.
    indexFile : str
        A file in json format containing indices (position in wig file) and
        sizes of chromosomes. If this file is not present and given as input,
        a new file will be generated. If this file is present, indices and
        sizes will be taken from this file. If index and size of input
        chromosome is not present in json file, these will be determined from
        wig file and stored in same json file. This file could be very helpful
        in case when same wig file has to be read many times because step to
        determine index and size of chromosome is skipped.
    methodToCombine : str
        method to combine bigWig/Wig files, Presently, accepted keywords are: ``mean``, ``min`` and ``max``
    tmpNumpyArrayFiles : :class:`TempNumpyArrayFiles`
        This :class:`TempNumpyArrayFiles` instance stores the temporary numpy array files information.
    isWigParsed : bool
        Whether Wig files are already parsed.
    maxEntryWrite : int
        Number of lines read from Wig file at an instant, after this, data is dumped in temporary numpy array file

    Parameters
    ----------
    filenames : str or list(str)
        List of input  Wig files.

        .. note:: In case if :attr:`WigHandler.chromName` is provided, only
                  one wig file is accepted.
    chromName : str
        Name of target chromosome name need to be extracted from wig file.
    chromSizeInfo : dict
        A dictionary containing chromosome size information. Generated by :meth:`BigWigHandler.getBigWigInfo`.
    indexFile : str
        A file in json format containing indices (position in wig file) and
        sizes of chromosomes. If this file is not present and given as input,
        a new file will be generated. If this file is present, indices and
        sizes will be taken from this file. If index and size of input
        chromosome is not present in json file, these will be determined from
        wig file and stored in same json file. This file could be very helpful
        in case when same wig file has to be read many times because step to
        determine index and size of chromosome is skipped.
    tmpNumpyArrayFiles : :class:`TempNumpyArrayFiles`
        This :class:`TempNumpyArrayFiles` instance stores the temporary numpy array files information.
    methodToCombine : str
        method to combine bigWig/Wig files, Presently, accepted keywords are: ``mean``, ``min`` and ``max``
    maxEntryWrite : int
        Number of lines read from Wig file at an instant, after this, data is dumped in temporary numpy array file. To reduce memory (RAM) occupancy,
        reduce this number because large numbers need large RAM.


    """
    def __init__(self, filenames, chromSizeInfo=None, chromName=None, indexFile=None, tmpNumpyArrayFiles=None, methodToCombine='mean',
                 workDir=None, maxEntryWrite=10000000):
        self.WigFileNames = filenames

        self.chromSizeInfo = chromSizeInfo         # dictionary for chromosome sizes, extracted from bigWigInfo program
        self.tmpNumpyArrayFiles = None             # List of temporary file to store array
        self.isWigParsed  = False

        # Working and output directory
        if workDir is None:
            self.workDir = config['Dirs']['WorkingDirectory']
        else:
            self.workDir = workDir

        self.logger = logging.getLogger(__name__+'.WigHandler')
        self.logger.setLevel(logging.INFO)

        self.maxEntryWrite = maxEntryWrite # maximum count before writing to output file will be attempted

        # For an input chromosome
        self.chromName = chromName
        self._chromPointerInFile = None # Dictionary for pointers in files
        self.indexFile = indexFile

        # Convert to list
        if not isinstance(filenames, list):
            self.WigFileNames = [ filenames ]

        if methodToCombine not in ['mean', 'max', 'min']:
            raise NotImplementedError(' Method [{0}] to combine bigwig file not implemented. Use: \'mean\', \'max\' or \'min\' ' .format(methodToCombine))
        self.methodToCombine = methodToCombine

        if self.chromSizeInfo is None:
            if len(self.WigFileNames) > 1:
                raise NotImplementedError('Presently only coversion for single wig file is implemented!!')

            if self.indexFile is not None:
                self._loadChromSizeAndIndex()
            self._getChromSizeInfo(self.WigFileNames[0], inputChrom=self.chromName)

        if tmpNumpyArrayFiles is None:
            self.tmpNumpyArrayFiles = TempNumpyArrayFiles(workDir=self.workDir)
            self.tmpNumpyArrayFiles.chromSizeInfo = copy.deepcopy(self.chromSizeInfo)
        else:
            self._checkExternTempNumpArrayFiles(tmpNumpyArrayFiles)
            self.tmpNumpyArrayFiles = tmpNumpyArrayFiles
            self.tmpNumpyArrayFiles.fillAllArraysWithZeros()

    def _checkExternTempNumpArrayFiles(self, tmpNumpyArrayFiles):
        if not isinstance(tmpNumpyArrayFiles, TempNumpyArrayFiles):
            raise TypeError(' Provided tmpNumpyArrayFiles is not instance of TempNumpyArrayFiles.')

        # check if instance already have information about chromosome size
        if tmpNumpyArrayFiles.chromSizeInfo is None:
            tmpNumpyArrayFiles.chromSizeInfo = copy.deepcopy(self.chromSizeInfo)
        else:
            for key in self.chromSizeInfo:
                tmpNumpyArrayFiles.updateArraysByChromSize(key, self.chromSizeInfo[key])

    def _PerformDataCoarsening(self, Chrom, resolution, coarsening_method):
        """Base method to perform Data coarsening.

        This method read temporary Numpy array files and perform data coarsening using the given input method.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`WigHandler._StoreInHdf5File`.

        Parameters
        ----------
        Chrom : str
            Chromosome name

        resolution : str
            resolution in word.

        coarsening_method : str
            Name of method to use for data coarsening. Accepted keywords: min, max, median, amean, gmean and hmean.


        """

        output = []
        binsize = util.resolutionToBinsize(resolution)
        size = self.chromSizeInfo[Chrom] + 1

        for i in range(1, size, binsize):
            tmpx = None
            if i+binsize >= size:
                tmpx = self.tmpNumpyArrayFiles.arrays[Chrom][i : size]
            else:
                tmpx = self.tmpNumpyArrayFiles.arrays[Chrom][i : i+binsize]

            int_idx = np.nonzero(tmpx > 0)

            if int_idx[0].shape[0] == 0:
                output.append(0.0)
                continue

            #print(Chrom, tmpx.shape, i, i+binsize, tmpx)
            if coarsening_method == 'max':
                output.append(np.amax(tmpx[int_idx]))

            if coarsening_method == 'min':
                output.append(np.amin(tmpx[int_idx]))

            if coarsening_method == 'amean':
                output.append(np.mean(tmpx[int_idx]))

            if coarsening_method == 'hmean':
                output.append(spstats.hmean(tmpx[int_idx]))

            if coarsening_method == 'gmean':
                output.append(spstats.gmean(tmpx[int_idx]))

            if coarsening_method == 'median':
                output.append(np.median(tmpx[int_idx]))

        # print(Chrom, resolution, coarse_method, size, binsize, size/binsize, len(output), np.amax(output))

        return np.asarray(output)

    def _getChromTitle_parseWig(self, line):
        """To parse chromosome title from the format line of fixedStep and variableStep format Wig file.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`WigHandler._parseWig`.

        Parameters
        ----------
        line : str
            The line containing chromosome information from Wig file

        """
        s = re.search('(?<=chrom=)\w+', line)
        if s:
            return s.group(0)

    def _getChromTitleBedgraph_parseWig(self, line):
        """To parse chromosome title from the format line of bedGraph format Wig file.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`WigHandler._parseWig`.

        Parameters
        ----------
        line : str
            The line containing chromosome information from Wig file

        """
        s = re.split('\s+', line)
        chrm = re.split(':', s[2])[0]
        return chrm

    def _getStartStepFixedStep_parseWig(self, line):
        """To parse start and step values the format line of fixedStep format Wig file.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`WigHandler._parseWig`.

        Parameters
        ----------
        line : str
            The line containing chromosome information from Wig file

        """
        s1 = re.search('(?<=start=)\d+', line)
        if s1:
            start = int(s1.group(0))

        s2 = re.search('(?<=step=)\d+', line)
        if s2:
            step = int(s2.group(0))
        else:
            step = 1

        return int(start), int(step)

    def _getSpan_parseWig(self, line):
        """To parse span value the format line of fixedStep format Wig file.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`WigHandler._parseWig`.

        Parameters
        ----------
        line : str
            The line containing chromosome information from Wig file

        """
        s = re.search('(?<=span=)\d+', line)
        if s:
            return int(s.group(0))
        else:
            return 1

    def _FillDataInNumpyArrayFile(self, ChromTitle, location_list, value_list):
        """Fill the extracted data from Wig file to temporary numpy array file

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`WigHandler._parseWig`.

        Parameters
        ----------
        ChromTitle : str
            Name of chromosome
        location_list : list of int
            List of locations for given chromosome
        value_list : list of float
            List of values for respecitve chromosome location

        """

        # Write numpy array data to numpy file
        self.tmpNumpyArrayFiles.generateTempNumpyFile(ChromTitle)

        valueInNumpyArray = self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list]
        boolNonZero = (valueInNumpyArray != 0)
        noZeroIdx = np.nonzero( boolNonZero  )
        zeroIdx = np.nonzero( ~boolNonZero  )
        value_list = np.asarray( value_list )

        if self.methodToCombine == 'mean':
            valueInNumpyArray[noZeroIdx] = (valueInNumpyArray[noZeroIdx] + value_list[noZeroIdx] ) / 2
            valueInNumpyArray[zeroIdx] = value_list[zeroIdx]
        elif self.methodToCombine == 'max':
            valueInNumpyArray = np.maximum( valueInNumpyArray, value_list )
        elif self.methodToCombine == 'min':
            valueInNumpyArray[noZeroIdx] = np.minimum( valueInNumpyArray[noZeroIdx], value_list[noZeroIdx] )
        else:
            valueInNumpyArray = value_list

        self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list] = valueInNumpyArray

        self.tmpNumpyArrayFiles.arrays[ChromTitle].flush()

    def _StoreInHdf5File(self, hdf5Out, title, resolutions=None, coarsening_methods=None, compression='lzf', keep_original=False):
        """ Base method to store coarsed data in hdf5 file.

        At first data is coarsened and subsequently stored in h5 file.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`WigHandler.saveAsH5`.

        Parameters
        ----------
        hdf5Out : str or :class:`HDF5Handler`
            Name of output hdf5 file or instance of :class:`HDF5Handler`
        title : str
            Title of data
        resolutions : list of str
            Additional input resolutions other than these default resolutions:
            1kb', '2kb', '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb',
            '100kb', '160kb','200kb', '320kb', '500kb', '640kb',  and '1mb'.

            For Example: use ``resolutions=['25kb', '50kb', '75kb']`` to add
            additional 25kb, 50kb and 75kb resolution data.
        coarsening_methods : list of str
            Methods to coarse or downsample the data for converting from 1-base
            to coarser resolutions. Presently, five methods are implemented.
            * ``'min'``    -> Minimum value
            * ``'max'``    -> Maximum value
            * ``'amean'``  -> Arithmatic mean or average
            * ``'hmean'``  -> Harmonic mean
            * ``'gmean'``  -> Geometric mean
            * ``'median'`` -> Median
            In case of ``None``, all five methods will be considered. User may
            use only subset of these methods. For example:
            ``coarse_method=['max', 'amean']`` can be used for downsampling by
            only these two methods.
        compression : str
            data compression method in HDF5 file : ``lzf`` or ``gzip`` method.
        keep_original : bool
            Whether original data present in wig file should be incorporated in HDF5 file. This will significantly increase size of HDF5 file.

        """
        # In cases when hdf5 filename or hdf5handler object is given as input
        # Initialize new hdf5handler in case file name is given
        toCloseHdfOut = False
        if not isinstance(hdf5Out, HDF5Handler):
            hdf5Out = HDF5Handler(hdf5Out, title=title)
            toCloseHdfOut = True

        resolutions = check_resolution_list(resolutions)

        # Writing to hdf5 aur appending to dictionary to return
        self.logger.info(' Writing output file [{0}] ...' .format(hdf5Out.filename))

        for Chrom in self.chromSizeInfo:
            if self.chromName is not None and self.chromName != Chrom:
                continue
            for res in resolutions:
                if keep_original:
                    hdf5Out.addDataByArray(Chrom, '1b', 'Original', self.tmpNumpyArrayFiles.arrays[Chrom][:], compression='gzip')
                for coarsening_method in check_coarsening_method(coarsening_methods):
                    hdf5Out.addDataByArray(Chrom, res, coarsening_method, self._PerformDataCoarsening(Chrom, res, coarsening_method), compression=compression )

        self.logger.info(' \t\t ... Finished writing output file [{0}] ' .format(hdf5Out.filename))
        if toCloseHdfOut:
            hdf5Out.close()

    def parseWig(self):
        """ To parse Wig files

        This method parses all Wig files listed in :attr:`WigHandler.WigFileNames`.
        The extracted data is further stored in temporary numpy array files of respective chromosome.
        These numpy array files can be used either for data coarsening or for further analysis.


        * **To save as h5:** Use :meth:`WigHandler.saveAsH5`.
        * **To perform analysis:** Use :meth:`WigHandler.getRawWigDataAsDictionary` to get a dictionary of numpy arrays.

        """
        # Parsing each file
        try:
            for f in self.WigFileNames:
                self._parseWig(f)
        except (SystemExit, KeyboardInterrupt) as e:
            del self.tmpNumpyArrayFiles
            raise e

        self.isWigParsed = True

    def _parseWig(self, wigFileName):
        """ Base method to parse a Wig file.

        This method parses a Wig file and extracted data are copied in temporary numpy array files.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`WigHandler.parseWig`.

        Parameters
        ----------
        wigFileName : str
            Name of Wig File

        """

        # Different format of output wig file
        bedGraph = False
        variableStep = False
        fixedStep = False
        ftype = None

        PreviousChromTitle, ChromTitle, location, value = 'dummy', 'dummy', None, None
        start, step, incr = None, None, None     # for fixedstep type
        span = 1 # default for both fixed and variable step
        dictOut = None # For default return
        printed = False

        fin = open(wigFileName, 'r')

        # If location in file is known previously, jump there
        if self.chromName is not None and self._chromPointerInFile is not None:
            if self.chromName in self._chromPointerInFile:
                fin.seek(self._chromPointerInFile[self.chromName])


        location_list, value_list = [], []
        ChromList = []
        counter = 0
        for line in fin:

            # To store data in temporary numpy array
            if (ChromTitle != PreviousChromTitle or counter >= self.maxEntryWrite) and not (len(location_list) < 1) :

                self._FillDataInNumpyArrayFile(PreviousChromTitle, location_list, value_list)

                del location_list
                del value_list
                location_list, value_list = [], []
                counter = 0


            # In case if input chromosome is already read, break here
            # Also Logging the information
            if PreviousChromTitle != ChromTitle:
                if PreviousChromTitle != 'dummy':
                    self.logger.info('         ... Finished reading and processing for {0} ' .format(PreviousChromTitle))
                    if PreviousChromTitle == self.chromName:
                        break

                if not printed:
                    self.logger.info(' #### Identified \'{0}\' format in WIG file... ##### ' .format(ftype))
                    printed = True

                self.logger.info(' Started reading and processing for {0} ' .format(ChromTitle))

            PreviousChromTitle = ChromTitle
            line = line.rstrip().lstrip()

            # Skip blank line
            if not line.strip():
                continue

            if re.search('variableStep', line) is not None:
                variableStep = True
                bedGraph = False
                fixedStep = False
                ChromTitle = self._getChromTitle_parseWig(line)
                ftype = 'variableStep'
                continue

            if re.search('#bedGraph', line) is not None:
                bedGraph = True
                variableStep = False
                fixedStep = False
                ChromTitle = self._getChromTitleBedgraph_parseWig(line)
                ftype = 'bedGraph'
                continue

            if re.search('fixedStep', line) is not None:
                fixedStep = True
                bedGraph = False
                variableStep = False
                ftype = 'fixedStep'
                ChromTitle = self._getChromTitle_parseWig(line)
                start, step = self._getStartStepFixedStep_parseWig(line)
                span = self._getSpan_parseWig(line)
                incr = 0
                continue

            temp = re.split('\s+', line)

            if variableStep:
                location = int(temp[0])
                value = float(temp[1])

            if fixedStep:
                location = start + incr
                value = float(temp[0])
                incr += 1

            if bedGraph:
                value = float(temp[3])

                #Zero based - start from zero - add one

                loc_start = int(temp[1]) + 1
                loc_end = int(temp[2]) + 1

                for l in range(loc_start, loc_end):
                    location_list.append(l)
                    value_list.append(value)
                    counter += 1
            else:
                for g in range(0, span):
                    location = location + g
                    location_list.append(location)
                    value_list.append(value)
                    counter += 1


            location = None
            value = None

        fin.close()

        # Last segment of lines are added in temporary numpy array file
        if self.chromName is None or self.chromName == ChromTitle:
            self._FillDataInNumpyArrayFile(ChromTitle, location_list, value_list)
            self.logger.info('         ... Finished reading and processing for {0} ' .format(ChromTitle))

    def saveAsH5(self, hdf5Out, title=None, resolutions=None, coarsening_methods=None, compression='lzf', keep_original=False):
        """To convert Wig files to hdf5 file

        Parameters
        ----------
        hdf5Out : :class:`HDF5Handler` or str
            Output hdf5 file name or :class:`HDF5Handler` instance
        title : str
            Title of the data
        resolutions : list of str
            Additional input resolutions other than these default resolutions:
            1kb', '2kb', '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb',
            '100kb', '160kb','200kb', '320kb', '500kb', '640kb',  and '1mb'.

            For Example: use ``resolutions=['25kb', '50kb', '75kb']`` to add
            additional 25kb, 50kb and 75kb resolution data.
        coarsening_methods : list of str
            Methods to coarse or downsample the data for converting from 1-base
            to coarser resolutions. Presently, five methods are implemented.

            * ``'min'``    -> Minimum value
            * ``'max'``    -> Maximum value
            * ``'amean'``  -> Arithmatic mean or average
            * ``'hmean'``  -> Harmonic mean
            * ``'gmean'``  -> Geometric mean
            * ``'median'`` -> Median

            In case of ``None``, all five methods will be considered. User may
            use only subset of these methods. For example:
            ``coarse_method=['max', 'amean']`` can be used for downsampling by
            only these two methods.
        compression : str
            data compression method in HDF5 file : ``lzf`` or ``gzip`` method.
        keep_original : bool
            Whether original data present in bigwig file should be incorporated in HDF5 file. This will significantly increase size of HDF5 file.

        """
        if not self.isWigParsed:
            self.parseWig()

        # Storing data in hdf5 file
        self._StoreInHdf5File(hdf5Out, title, compression=compression, coarsening_methods=coarsening_methods, resolutions=resolutions, keep_original=keep_original)

    def getRawWigDataAsDictionary(self, dicOut=None):
        """To get a entire dictionary of data from Wig file

        It generates a dictionary of numpy arrays for each chromosome. These arrays are stored in temporary numpy array files of :class:`TempNumpyArrayFiles`.

        Parameters
        ----------
        dicOut : dict
            The output dictionary to which data will be added or replaced.

        Return
        ------
        dicOut : dict
            The output dictionary.

        """

        noReturn = True
        if dicOut is None:
            dictOut = dict()
            noReturn = False

        if not self.isWigParsed:
            self.parseWig()

        for key in self.chromSizeInfo:
            dictOut[key] = self.tmpNumpyArrayFiles.arrays[key]

        if not noReturn:
            return dicOut

    def setChromosome(self, chromName):
        """ Set the target chromosome for reading and extracting from wig file

        To read and convert data of another chromsome from a wig file, it
        can be set here. After this, directly use :meth:`WigHandler.saveAsH5`
        to save data in H5 file.

        Parameters
        ----------
        chromName : str
            Name of new target chromosome

        """
        if len(self.WigFileNames) > 1:
            raise NotImplementedError('Presently only coversion for single wig file is implemented!!')

        self.chromName = chromName
        self.isWigParsed = False           # Reset Wig Parsing flag

        if self.chromName in self.chromSizeInfo and self.chromName in self._chromPointerInFile:
            return
        else:
            # Get size of new chromosome
            self._getChromSizeInfo(self.WigFileNames[0], inputChrom=self.chromName)

            # Also update to chromSizeInfo in TempNumpyArrayFiles
            if self.chromName in self.tmpNumpyArrayFiles.chromSizeInfo:
                 self.tmpNumpyArrayFiles.chromSizeInfo[self.chromName] = max(self.chromSizeInfo[self.chromName], self.tmpNumpyArrayFiles.chromSizeInfo[self.chromName])
            else:
                self.tmpNumpyArrayFiles.chromSizeInfo[self.chromName] = self.chromSizeInfo[self.chromName]

    def _getChromSizeInfo(self, wigFileName, inputChrom=None):
        """ Get chromosome size and index wig file

        This method parses a Wig file, extracts chromosome size and index it
        for each chromosome.

        It sets :attr:`WigHandler._chromPointerInFile` and
        :attr:`WigHandler.chromSizeInfo`.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally
            durng initialization and in :meth:`WigHandler.setChromosome`.



        Parameters
        ----------
        wigFileName : str
            Name of Wig File

        inputChrom : str
            Name of target chromosome

        """

        # Initialize dictionary
        if self.chromSizeInfo is None:
            self.chromSizeInfo = dict()

        # Different format of output wig file
        bedGraph = False
        variableStep = False
        fixedStep = False
        ftype = None

        PreviousChromTitle, ChromTitle, location, PreviousLine = 'dummy', 'dummy', None, None
        start, step, incr = None, None, None     # for fixedstep type
        span = 1 # default for both fixed and variable step
        dictOut = None # For default return
        printed = False

        fin = open(wigFileName, 'r')

        # In case if file is already read previously, continue from there
        pointer = 0
        if self._chromPointerInFile is not None:
            pointer = max(self._chromPointerInFile.values())
            fin.seek(pointer)

        # Initialize pointer for chromosome in wig file
        if self._chromPointerInFile is None:
            self._chromPointerInFile = dict()

        # if chromName is already in both dictionary, just return from here
        if self.chromName in self.chromSizeInfo and self.chromName in self._chromPointerInFile:
            return

        currentLocation = None
        ChromList = []
        counter = 0
        for line in fin:

            # To store data in temporary numpy array
            if (ChromTitle != PreviousChromTitle) and currentLocation is not None :
                self.chromSizeInfo[PreviousChromTitle] = currentLocation
                currentLocation = None

            if PreviousChromTitle != ChromTitle:

                # Store location of new chromosome
                self._chromPointerInFile[ChromTitle] = pointer - len(PreviousLine)

                # If index file is provided, update the index file
                if self.indexFile is not None:
                    self._saveChromSizeAndIndex()

                # Logging the information
                if PreviousChromTitle != 'dummy':
                    self.logger.info('         ... Got maximum size of {0} for {1}' .format(self.chromSizeInfo[PreviousChromTitle], PreviousChromTitle))

                    # If got input chromosome, abort further processing
                    if inputChrom == PreviousChromTitle:
                        break


                if not printed:
                    self.logger.info(' #### Identified \'{0}\' format in WIG file... ##### ' .format(ftype))
                    printed = True

                self.logger.info(' Searching and Indexing for chromosome {0} ... ' .format(ChromTitle))


            pointer += len(line)
            PreviousChromTitle = ChromTitle
            PreviousLine = line

            line = line.rstrip().lstrip()

            # Skip blank line
            if not line.strip():
                continue

            if re.search('variableStep', line) is not None:
                variableStep = True
                bedGraph = False
                fixedStep = False
                ChromTitle = self._getChromTitle_parseWig(line)
                ftype = 'variableStep'
                continue

            if re.search('#bedGraph', line) is not None:
                bedGraph = True
                variableStep = False
                fixedStep = False
                ChromTitle = self._getChromTitleBedgraph_parseWig(line)
                ftype = 'bedGraph'
                continue

            if re.search('fixedStep', line) is not None:
                fixedStep = True
                bedGraph = False
                variableStep = False
                ftype = 'fixedStep'
                ChromTitle = self._getChromTitle_parseWig(line)
                start, step = self._getStartStepFixedStep_parseWig(line)
                span = self._getSpan_parseWig(line)
                incr = 0
                continue

            temp = re.split('\s+', line)

            if variableStep:
                location = int(temp[0])

            if fixedStep:
                location = start + incr
                incr += 1

            if bedGraph:
                location = int(temp[2]) + 1
            else:
                location = location + span

            if location is not None:
                if currentLocation is None:
                    currentLocation = location
                else:
                    currentLocation = max(location, currentLocation)
            location = None

        fin.close()

        # Last segment of lines
        if inputChrom is None or inputChrom == ChromTitle:
            self.chromSizeInfo[ChromTitle] = currentLocation
            self.logger.info('         ... Got maximum size of {0} for {1}' .format(self.chromSizeInfo[PreviousChromTitle], PreviousChromTitle))


        # If index file is provided, update the index file
        if self.indexFile is not None:
            self._saveChromSizeAndIndex()

    def _saveChromSizeAndIndex(self):
        """ Save chromosomes sizes and indices dictionary to a json file
        """
        data = dict()
        if self.chromSizeInfo is not None:
            data['size-info'] = self.chromSizeInfo

        if self._chromPointerInFile is not None:
            data['index'] = self._chromPointerInFile

        fout =  open( os.path.abspath( os.path.expanduser(self.indexFile)), "w" )
        json.dump(data, fout, indent=4, separators=(',', ':'))
        fout.close()

    def _loadChromSizeAndIndex(self):
        """ Load chromosome sizes and indices from a json file
        """
        try:
            fin = open( self.indexFile, "r" )
        except:
            return

        data = json.load( fin )
        fin.close()

        if not data:
            return

        if 'size-info' in data:
            self.chromSizeInfo = data['size-info']

        if 'index' in data:
            self._chromPointerInFile = data['index']


class BEDHandler:
    """To convert BED files to hdf5/h5 file

    It parses bed files and save all data to a hdf5/h5 file for
    given resolutions.

    Attributes
    ----------
    BedFileNames : list[str]
        List of input  bed files.

        .. note:: In case if :attr:`BEDHandler.chromName` is provided, only
                  one wig file is accepted.

    column : int
        The column number, which is considered as data column. Column number
        could vary and depends on BED format. For example:

        * ENCODE broadPeak format (BED 6+3): 7th column
        * ENCODE gappedPeak format (BED 12+3): 13th column
        * ENCODE narrowPeak format (BED 6+4): 7th column
        * ENCODE RNA elements format (BED 6+3): 7th column

    chromName : str
        Name of target chromosome name need to be extracted from bed file.
    chromSizeInfo : dict
        A dictionary containing chromosome size information.
    _chromPointerInFile : dict
        A dictionary containing position index of each chromosome in bed file.
    indexFile : str
        A file in json format containing indices (position in bed file) and
        sizes of chromosomes. If this file is not present and given as input,
        a new file will be generated. If this file is present, indices and
        sizes will be taken from this file. If index and size of input
        chromosome is not present in json file, these will be determined from
        bed file and stored in same json file. This file could be very helpful
        in case when same wig file has to be read many times because step to
        determine index and size of chromosome is skipped.
    methodToCombine : str
        method to combine bed files, Presently, accepted keywords are: ``mean``, ``min`` and ``max``
    tmpNumpyArrayFiles : :class:`TempNumpyArrayFiles`
        This :class:`TempNumpyArrayFiles` instance stores the temporary numpy array files information.
    isBedParsed : bool
        Whether bed files are already parsed.
    maxEntryWrite : int
        Number of lines read from bed file at an instant, after this, data is dumped in temporary numpy array file

    Parameters
    ----------
    filenames : str or list(str)
        List of input  bed files.

        .. note:: In case if :attr:`BEDHandler.chromName` is provided, only
                  one bed file is accepted.

    column : int
        The column number, which is considered as data column. Column number
        could vary and depends on BED format. For example:

        * ENCODE broadPeak format (BED 6+3): 7th column
        * ENCODE gappedPeak format (BED 12+3): 13th column
        * ENCODE narrowPeak format (BED 6+4): 7th column
        * ENCODE RNA elements format (BED 6+3): 7th column

    chromName : str
        Name of target chromosome name need to be extracted from bed file.
    indexFile : str
        A file in json format containing indices (position in bed file) and
        sizes of chromosomes. If this file is not present and given as input,
        a new file will be generated. If this file is present, indices and
        sizes will be taken from this file. If index and size of input
        chromosome is not present in json file, these will be determined from
        bed file and stored in same json file. This file could be very helpful
        in case when same bed file has to be read many times because step to
        determine index and size of chromosome is skipped.
    tmpNumpyArrayFiles : :class:`TempNumpyArrayFiles`
        This :class:`TempNumpyArrayFiles` instance stores the temporary numpy array files information.
    methodToCombine : str
        method to combine bed files, Presently, accepted keywords are: ``mean``, ``min`` and ``max``
    maxEntryWrite : int
        Number of lines read from bed file at an instant, after this, data is dumped in temporary numpy array file. To reduce memory (RAM) occupancy,
        reduce this number because large numbers need large RAM.


    """
    def __init__(self, filenames, column=7, chromName=None, indexFile=None, tmpNumpyArrayFiles=None, methodToCombine='mean',
                 workDir=None, maxEntryWrite=10000000):
        self.bedFileNames = filenames
        self.column = column

        self.chromSizeInfo = None         # dictionary for chromosome sizes, extracted from bigWigInfo program
        self.tmpNumpyArrayFiles = None             # List of temporary file to store array
        self.isBedParsed  = False

        # Working and output directory
        if workDir is None:
            self.workDir = config['Dirs']['WorkingDirectory']
        else:
            self.workDir = workDir

        self.logger = logging.getLogger(__name__+'.WigHandler')
        self.logger.setLevel(logging.INFO)

        self.maxEntryWrite = maxEntryWrite # maximum count before writing to output file will be attempted

        # For an input chromosome
        self.chromName = chromName
        self._chromPointerInFile = None # Dictionary for pointers in files
        self.indexFile = indexFile

        # Convert to list
        if not isinstance(filenames, list):
            self.bedFileNames = [ filenames ]

        if methodToCombine not in ['mean', 'max', 'min']:
            raise NotImplementedError(' Method [{0}] to combine bed file not implemented. Use: \'mean\', \'max\' or \'min\' ' .format(methodToCombine))
        self.methodToCombine = methodToCombine

        if self.chromSizeInfo is None:
            if len(self.bedFileNames) > 1:
                raise NotImplementedError('Presently only coversion for single bed file is implemented!!')

            if self.indexFile is not None:
                self._loadChromSizeAndIndex()
            self._getChromSizeInfo(self.bedFileNames[0], inputChrom=self.chromName)

        if tmpNumpyArrayFiles is None:
            self.tmpNumpyArrayFiles = TempNumpyArrayFiles(workDir=self.workDir)
            self.tmpNumpyArrayFiles.chromSizeInfo = copy.deepcopy(self.chromSizeInfo)
        else:
            self._checkExternTempNumpArrayFiles(tmpNumpyArrayFiles)
            self.tmpNumpyArrayFiles = tmpNumpyArrayFiles
            self.tmpNumpyArrayFiles.fillAllArraysWithZeros()

    def _checkExternTempNumpArrayFiles(self, tmpNumpyArrayFiles):
        if not isinstance(tmpNumpyArrayFiles, TempNumpyArrayFiles):
            raise TypeError(' Provided tmpNumpyArrayFiles is not instance of TempNumpyArrayFiles.')

        # check if instance already have information about chromosome size
        if tmpNumpyArrayFiles.chromSizeInfo is None:
            tmpNumpyArrayFiles.chromSizeInfo = copy.deepcopy(self.chromSizeInfo)
        else:
            for key in self.chromSizeInfo:
                tmpNumpyArrayFiles.updateArraysByChromSize(key, self.chromSizeInfo[key])

    def _PerformDataCoarsening(self, Chrom, resolution, coarse_method):
        """Base method to perform Data coarsening.

        This method read temporary Numpy array files and perform data coarsening using the given input method.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`BEDHandler._StoreInHdf5File`.

        Parameters
        ----------
        Chrom : str
            Chromosome name

        resolution : str
            resolution in word.

        coarse_method : str
            Name of method to use for data coarsening. Accepted keywords: min, max, median, amean, gmean and hmean.


        """

        output = []
        binsize = util.resolutionToBinsize(resolution)
        size = self.chromSizeInfo[Chrom] + 1

        for i in range(1, size, binsize):
            tmpx = None
            if i+binsize >= size:
                tmpx = self.tmpNumpyArrayFiles.arrays[Chrom][i : size]
            else:
                tmpx = self.tmpNumpyArrayFiles.arrays[Chrom][i : i+binsize]

            int_idx = np.nonzero(tmpx > 0)

            if int_idx[0].shape[0] == 0:
                output.append(0.0)
                continue

            #print(Chrom, tmpx.shape, i, i+binsize, tmpx)
            if coarse_method == 'max':
                output.append(np.amax(tmpx[int_idx]))

            if coarse_method == 'min':
                output.append(np.amin(tmpx[int_idx]))

            if coarse_method == 'amean':
                output.append(np.mean(tmpx[int_idx]))

            if coarse_method == 'hmean':
                output.append(spstats.hmean(tmpx[int_idx]))

            if coarse_method == 'gmean':
                output.append(spstats.gmean(tmpx[int_idx]))

            if coarse_method == 'median':
                output.append(np.median(tmpx[int_idx]))

        # print(Chrom, resolution, coarse_method, size, binsize, size/binsize, len(output), np.amax(output))

        return np.asarray(output)

    def _FillDataInNumpyArrayFile(self, ChromTitle, location_list, value_list):
        """Fill the extracted data from bed file to temporary numpy array file

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`BEDHandler._parseBed`.

        Parameters
        ----------
        ChromTitle : str
            Name of chromosome
        location_list : list of int
            List of locations for given chromosome
        value_list : list of float
            List of values for respecitve chromosome location

        """

        # Write numpy array data to numpy file
        self.tmpNumpyArrayFiles.generateTempNumpyFile(ChromTitle)

        valueInNumpyArray = self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list]
        boolNonZero = (valueInNumpyArray != 0)
        noZeroIdx = np.nonzero( boolNonZero  )
        zeroIdx = np.nonzero( ~boolNonZero  )
        value_list = np.asarray( value_list )

        if self.methodToCombine == 'mean':
            valueInNumpyArray[noZeroIdx] = (valueInNumpyArray[noZeroIdx] + value_list[noZeroIdx] ) / 2
            valueInNumpyArray[zeroIdx] = value_list[zeroIdx]
        elif self.methodToCombine == 'max':
            valueInNumpyArray = np.maximum( valueInNumpyArray, value_list )
        elif self.methodToCombine == 'min':
            valueInNumpyArray[noZeroIdx] = np.minimum( valueInNumpyArray[noZeroIdx], value_list[noZeroIdx] )
        else:
            valueInNumpyArray = value_list

        self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list] = valueInNumpyArray

        self.tmpNumpyArrayFiles.arrays[ChromTitle].flush()

    def _StoreInHdf5File(self, hdf5Out, title, resolutions=None, coarsening_methods=None, compression='lzf', keep_original=False):
        """ Base method to store coarsened data in hdf5/h5 file.

        At first data is coarsened and subsequently stored in h5 file.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`BEDHandler.saveAsH5`.

        Parameters
        ----------
        hdf5Out : str or :class:`HDF5Handler`
            Name of output hdf5 file or instance of :class:`HDF5Handler`
        title : str
            Title of data
        resolutions : list of str
            Additional input resolutions other than these default resolutions:
            1kb', '2kb', '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb',
            '100kb', '160kb','200kb', '320kb', '500kb', '640kb',  and '1mb'.

            For Example: use ``resolutions=['25kb', '50kb', '75kb']`` to add
            additional 25kb, 50kb and 75kb resolution data.
        coarsening_methods : list of str
            Methods to coarse or downsample the data for converting from 1-base
            to coarser resolutions. Presently, five methods are implemented.

            * ``'min'``    -> Minimum value
            * ``'max'``    -> Maximum value
            * ``'amean'``  -> Arithmatic mean or average
            * ``'hmean'``  -> Harmonic mean
            * ``'gmean'``  -> Geometric mean
            * ``'median'`` -> Median

            In case of ``None``, all five methods will be considered. User may
            use only subset of these methods. For example:
            ``coarse_method=['max', 'amean']`` can be used for downsampling by
            only these two methods.
        compression : str
            data compression method in HDF5 file : ``lzf`` or ``gzip`` method.
        keep_original : bool
            Whether original data present in wig file should be incorporated in HDF5 file. This will significantly increase size of HDF5 file.

        """
        # In cases when hdf5 filename or hdf5handler object is given as input
        # Initialize new hdf5handler in case file name is given
        toCloseHdfOut = False
        if not isinstance(hdf5Out, HDF5Handler):
            hdf5Out = HDF5Handler(hdf5Out, title=title)
            toCloseHdfOut = True


        resolutions = check_resolution_list(resolutions)

        # Writing to hdf5 aur appending to dictionary to return
        self.logger.info(' Writing output file [{0}] ...' .format(hdf5Out.filename))

        for Chrom in self.chromSizeInfo:
            if self.chromName is not None and self.chromName != Chrom:
                continue
            for res in resolutions:
                if keep_original:
                    hdf5Out.addDataByArray(Chrom, '1b', 'Original', self.tmpNumpyArrayFiles.arrays[Chrom][:], compression='gzip')
                for coarsening_method in check_coarsening_method(coarsening_methods):
                    hdf5Out.addDataByArray(Chrom, res, coarsening_method, self._PerformDataCoarsening(Chrom, res, coarsening_method), compression=compression )

        self.logger.info(' \t\t ... Finished writing output file [{0}] ' .format(hdf5Out.filename))
        if toCloseHdfOut:
            hdf5Out.close()

    def parseBed(self):
        """ To parse bed files

        This method parses all bed files listed in :attr:`BEDHandler.bedFileNames`.
        The extracted data is further stored in temporary numpy array files of respective chromosome.
        These numpy array files can be used either for data coarsening or for further analysis.


        * **To save as h5:** Use :meth:`BEDHandler.saveAsH5`.
        * **To perform analysis:** Use :meth:`BEDHandler.getRawWigDataAsDictionary` to get a dictionary of numpy arrays.

        """
        # Parsing each file
        try:
            for f in self.bedFileNames:
                self._parseBed(f)
        except (SystemExit, KeyboardInterrupt) as e:
            del self.tmpNumpyArrayFiles
            raise e

        self.isBedParsed = True

    def _parseBed(self, bedFileName):
        """ Base method to parse a bed file.

        This method parses a bed file and extracted data are copied in temporary numpy array files.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`BEDHandler.parseBed`.

        Parameters
        ----------
        bedFileName : str
            Name of bed File

        """

        PreviousChromTitle, ChromTitle, location, value = 'dummy', 'dummy', None, None
        dictOut = None # For default return

        fin = open(bedFileName, 'r')

        # If location in file is known previously, jump there
        if self.chromName is not None and self._chromPointerInFile is not None:
            if self.chromName in self._chromPointerInFile:
                fin.seek(self._chromPointerInFile[self.chromName])


        location_list, value_list = [], []
        ChromList = []
        counter = 0
        for line in fin:

            temp = re.split('\s+', line.rstrip().lstrip())

            # Skip lines that start with track or browser, but count the words
            if re.search('track', line) != None or re.search('browser', line) != None:
                continue

            ChromTitle = temp[0]

            # To store data in temporary numpy array
            if (ChromTitle != PreviousChromTitle or counter >= self.maxEntryWrite) and not (len(location_list) < 1) :

                self._FillDataInNumpyArrayFile(PreviousChromTitle, location_list, value_list)

                del location_list
                del value_list
                location_list, value_list = [], []
                counter = 0


            # In case if input chromosome is already read, break here
            # Also Logging the information
            if PreviousChromTitle != ChromTitle:

                # If index file is provided, update the index file
                if self.indexFile is not None:
                    self._saveChromSizeAndIndex()

                if PreviousChromTitle != 'dummy':
                    self.logger.info('         ... Finished reading and processing for {0} ' .format(PreviousChromTitle))
                    if PreviousChromTitle == self.chromName:
                        break

            PreviousChromTitle = ChromTitle
            loc_start = int(temp[1]) + 1
            loc_end = int(temp[2]) + 1
            value = float( temp[self.column-1] )

            for l in range(loc_start, loc_end):
                location_list.append(l)
                value_list.append(value)
                counter += 1

        fin.close()

        # Last segment of lines are added in temporary numpy array file
        if self.chromName is None or self.chromName == ChromTitle:
            self._FillDataInNumpyArrayFile(ChromTitle, location_list, value_list)
            self.logger.info('         ... Finished reading and processing for {0} ' .format(ChromTitle))

        # If index file is provided, update the index file
        if self.indexFile is not None:
            self._saveChromSizeAndIndex()


    def saveAsH5(self, hdf5Out, title=None, resolutions=None, coarsening_methods=None, compression='lzf', keep_original=False):
        """To convert bed files to hdf5 file

        It parses bed files, coarsened the data and store in an input hdf5/h5
        file.

        Parameters
        ----------
        hdf5Out : :class:`HDF5Handler` or str
            Output hdf5 file name or :class:`HDF5Handler` instance
        title : str
            Title of the data
        resolutions : list of str
            Additional input resolutions other than these default resolutions:
            1kb', '2kb', '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb',
            '100kb', '160kb','200kb', '320kb', '500kb', '640kb',  and '1mb'.

            For Example: use ``resolutions=['25kb', '50kb', '75kb']`` to add
            additional 25kb, 50kb and 75kb resolution data.
        coarsening_methods : list of str
            Methods to coarse or downsample the data for converting from 1-base
            to coarser resolutions. Presently, five methods are implemented.

            * ``'min'``    -> Minimum value
            * ``'max'``    -> Maximum value
            * ``'amean'``  -> Arithmatic mean or average
            * ``'hmean'``  -> Harmonic mean
            * ``'gmean'``  -> Geometric mean
            * ``'median'`` -> Median

            In case of ``None``, all five methods will be considered. User may
            use only subset of these methods. For example:
            ``coarse_method=['max', 'amean']`` can be used for downsampling by
            only these two methods.
        compression : str
            data compression method in HDF5 file : ``lzf`` or ``gzip`` method.
        keep_original : bool
            Whether original data present in bigwig file should be incorporated in HDF5 file. This will significantly increase size of HDF5 file.

        """
        if not self.isBedParsed:
            self.parseBed()

        # Storing data in hdf5 file
        self._StoreInHdf5File(hdf5Out, title, resolutions=resolutions, coarsening_methods=coarsening_methods, compression=compression, keep_original=keep_original)

    def getRawWigDataAsDictionary(self, dicOut=None):
        """To get a entire dictionary of data from bed file

        It generates a dictionary of numpy arrays for each chromosome.
        These arrays are stored in temporary numpy array files of
        :class:`TempNumpyArrayFiles`.

        Parameters
        ----------
        dicOut : dict
            The output dictionary to which data will be added or replaced.

        Return
        ------
        dicOut : dict
            The output dictionary.

        """

        noReturn = True
        if dicOut is None:
            dictOut = dict()
            noReturn = False

        if not self.isBedParsed:
            self.parseBed()

        for key in self.chromSizeInfo:
            dictOut[key] = self.tmpNumpyArrayFiles.arrays[key]

        if not noReturn:
            return dicOut

    def setChromosome(self, chromName):
        """ Set the target chromosome for reading and extracting from bed file

        To read and convert data of another chromsome from a bed file, it
        can be set here. After this, directly use :meth:`BEDHandler.saveAsH5`
        to save data in H5 file.

        Parameters
        ----------
        chromName : str
            Name of new target chromosome

        """
        if len(self.bedFileNames) > 1:
            raise NotImplementedError('Presently only coversion for single wig file is implemented!!')

        self.chromName = chromName
        self.isBedParsed = False           # Reset Wig Parsing flag

        if self.chromName in self.chromSizeInfo and self.chromName in self._chromPointerInFile:
            return
        else:
            # Get size of new chromosome
            self._getChromSizeInfo(self.bedFileNames[0], inputChrom=self.chromName)

            # Also update to chromSizeInfo in TempNumpyArrayFiles
            if self.chromName in self.tmpNumpyArrayFiles.chromSizeInfo:
                 self.tmpNumpyArrayFiles.chromSizeInfo[self.chromName] = max(self.chromSizeInfo[self.chromName], self.tmpNumpyArrayFiles.chromSizeInfo[self.chromName])
            else:
                self.tmpNumpyArrayFiles.chromSizeInfo[self.chromName] = self.chromSizeInfo[self.chromName]

    def _getChromSizeInfo(self, bedFileName, inputChrom=None):
        """ Get chromosome size and index bed file

        This method parses a bed file, extracts chromosome size and index it
        for each chromosome.

        It sets :attr:`BEDHandler._chromPointerInFile` and
        :attr:`BEDHandler.chromSizeInfo`.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally
            durng initialization and in :meth:`BEDHandler.setChromosome`.



        Parameters
        ----------
        bedFileName : str
            Name of Wig File

        inputChrom : str
            Name of target chromosome

        """

        # Initialize dictionary
        if self.chromSizeInfo is None:
            self.chromSizeInfo = dict()

        PreviousChromTitle, ChromTitle, location, PreviousLine = 'dummy', 'dummy', None, None
        start, end = None, None
        dictOut = None # For default return

        fin = open(bedFileName, 'r')

        # In case if file is already read previously, continue from there
        pointer = 0
        if self._chromPointerInFile is not None:
            pointer = max(self._chromPointerInFile.values())
            fin.seek(pointer)

        # Initialize pointer for chromosome in bed file
        if self._chromPointerInFile is None:
            self._chromPointerInFile = dict()

        # if chromName is already in both dictionary, just return from here
        if self.chromName in self.chromSizeInfo and self.chromName in self._chromPointerInFile:
            return

        currentLocation = None
        ChromList = []
        counter = 0
        for line in fin:

            temp = re.split('\s+', line.rstrip().lstrip())
            if len(temp) < self.column:
                raise AssertionError(' Number of columns [{0}] is less than the requested column [{1}] .'.format(len(temp), self.column))

            # Skip lines that start with track or browser, but count the words
            if re.search('track', line) != None or re.search('browser', line) != None:
                pointer += len(line)
                continue

            ChromTitle = temp[0]

            if (ChromTitle != PreviousChromTitle) and currentLocation is not None :
                self.chromSizeInfo[PreviousChromTitle] = currentLocation
                currentLocation = None

            if (ChromTitle != PreviousChromTitle):
                self._chromPointerInFile[ChromTitle] = pointer
                # Logging the information
                if PreviousChromTitle != 'dummy':
                    self.logger.info('         ... Got maximum size of {0} for {1}' .format(self.chromSizeInfo[PreviousChromTitle], PreviousChromTitle))

                    # If got input chromosome, abort further processing
                    if inputChrom == PreviousChromTitle:
                        break

                self.logger.info(' Searching and Indexing for chromosome {0} ... ' .format(ChromTitle))

            pointer += len(line)
            PreviousChromTitle = ChromTitle

            location = int(temp[2])
            if currentLocation is None:
                currentLocation = location
            else:
                currentLocation = max(location, currentLocation)

        fin.close()

        # Last segment of lines
        if inputChrom is None or inputChrom == ChromTitle:
            self.chromSizeInfo[ChromTitle] = currentLocation
            self.logger.info('         ... Got maximum size of {0} for {1}' .format(self.chromSizeInfo[PreviousChromTitle], PreviousChromTitle))

        # If index file is provided, update the index file
        if self.indexFile is not None:
            self._saveChromSizeAndIndex()

    def _saveChromSizeAndIndex(self):
        """ Save chromosomes sizes and indices dictionary to a json file
        """
        data = dict()
        if self.chromSizeInfo is not None:
            data['size-info'] = self.chromSizeInfo

        if self._chromPointerInFile is not None:
            data['index'] = self._chromPointerInFile

        fout =  open( os.path.abspath( os.path.expanduser(self.indexFile)), "w" )
        json.dump(data, fout, indent=4, separators=(',', ':'))
        fout.close()

    def _loadChromSizeAndIndex(self):
        """ Load chromosome sizes and indices from a json file
        """
        try:
            fin = open( self.indexFile, "r" )
        except:
            return

        data = json.load( fin )
        fin.close()

        if not data:
            return

        if 'size-info' in data:
            self.chromSizeInfo = data['size-info']

        if 'index' in data:
            self._chromPointerInFile = data['index']


class EncodeDatasetsConverter:
    """ Download and convert datasets from ENCODE Experiments matrix

    It can be used to download and convert multiple datasets from ENCODE Experiment
    matrix (https://www.encodeproject.org/matrix/?type=Experiment).
    Presently, only bigWig files are downloaded and then converted.

    At first search the datasets on https://www.encodeproject.org/matrix/?type=Experiment .
    Then click on download button on top of the page. A text file will be downloaded.
    This text file can be used as input in this program. All bigWig files will be
    downloaded and converted to gcMapExplorer compatible hdf5 format.

    .. note:: At first a metafile is automatically downloaded and then files
              are filtered according to bigWig format and Assembly. Subsequently,
              if several replicates are present, only datasets with combined
              replicates are considered. In case if two replicates are present
              and combined replicates are not present, first replicate will be
              considered. Combining replicates are not yet implemented

    .. warning:: Presently ``bigWigToWig`` and ``bigWigInfo`` is not available
                 for Windows OS. Therefore, this class will fail in this OS.

    Attributes
    ----------
    inputFile : str
        Name of input file downloaded from ENCODE Experiments matrix website.
    assembly : str
        Name of reference genome. Example: hg19, GRCh38 etc.
    pathTobigWigToWig : str
        Path to ``bigWigToWig`` program. It can be downloaded from
        http://hgdownload.cse.ucsc.edu/admin/exe/ for MacOSX and Linux.
        If path to program is already present in configuration file, it will be
        taken from the configuration.

        If it is not present in configuration file, the input path **should**
        be provided. It will be stored in configuration file for later use.
    pathTobigWigInfo : str
        Path to ``bigWigInfo`` program. It can be downloaded from
        http://hgdownload.cse.ucsc.edu/admin/exe/ for MacOSX and Linux.
        If path to program is already present in configuration file, it will be
        taken from the configuration.

        If it is not present in configuration file, the input path **should**
        be provided. It will be stored in configuration file for later use.
    metafile : str
        Name of metafile downloaded from ENCODE website. It is automatically
        downlaoded from input file. It contains all the meta-data required
        for processing.

    metaData = list of dictionary
        A list of dictionary read from metafile. It is already filtered
        according to the different criteria such as file-format, assembly,
        replicates.

    Parameters
    ----------
    inputFile : str
        Name of input file downloaded from ENCODE Experiments matrix website.
    assembly : str
        Name of reference genome. Example: hg19, GRCh38 etc.
    pathTobigWigToWig : str
        Path to ``bigWigToWig`` program. It can be downloaded from
        http://hgdownload.cse.ucsc.edu/admin/exe/ for MacOSX and Linux.
        If path to program is already present in configuration file, it will be
        taken from the configuration.

        If it is not present in configuration file, the input path **should**
        be provided. It will be stored in configuration file for later use.
    pathTobigWigInfo : str
        Path to ``bigWigInfo`` program. It can be downloaded from
        http://hgdownload.cse.ucsc.edu/admin/exe/ for MacOSX and Linux.
        If path to program is already present in configuration file, it will be
        taken from the configuration.

        If it is not present in configuration file, the input path **should**
        be provided. It will be stored in configuration file for later use.

    """

    def __init__(self, inputFile, assembly, methodToCombine='mean', pathTobigWigToWig=None, pathTobigWigInfo=None, workDir=None):
        self.inputFile = inputFile
        self.assembely = assembly
        self.metafile = None
        self.metaData = None
        self.indexes = None
        self._bigWigFiles = None
        self._checkPointFile = 'checkpoint.txt'
        self._checkPoint = []
        self._implementedAssays = ['ChIP-seq', 'siRNA knockdown followed by RNA-seq', 'RNA-seq', 'DNase-seq', 'FAIRE-seq']

        if methodToCombine not in ['mean', 'max', 'min']:
            raise NotImplementedError(' Method [{0}] to combine bigwig file not implemented. Use: \'mean\', \'max\' or \'min\' ' .format(methodToCombine))
        self.methodToCombine = methodToCombine

        # Working and output directory
        if workDir is None:
            self.workDir = config['Dirs']['WorkingDirectory']
        else:
            self.workDir = workDir

        self.logger = logging.getLogger(__name__+'.DownloadEncodeDatasets')
        self.logger.setLevel(logging.INFO)

        self._checkBigWigInfoProgram(pathTobigWigInfo)
        self._checkBigWigToWigProgram(pathTobigWigToWig)
        self.downloadMetaData()

    def __del__(self):
        if os.path.isfile(self.metafile):
            os.remove(self.metafile)

        if self._bigWigFiles is not None:
            self._removeBigWigFiles()

    def _checkBigWigInfoProgram(self, pathTobigWigInfo):
        """ Check if bigWigInfo program is available or accessible.

        If program is not available in configuration file, the given
        path will be stored in the file after checking its accessibility.

        The path is stored in :attr:`gcMapExplorer.lib.genomicsDataHandler.EncodeDatasetsConverter.pathTobigWigInfo`

        Parameters
        ----------
        pathTobigWigInfo : str
            Path to bigWigInfo program

        """
        if pathTobigWigInfo is None or pathTobigWigInfo == config['Programs']['bigWigInfo']:
            self.pathTobigWigInfo = config['Programs']['bigWigInfo']
            if self.pathTobigWigInfo == 'None':
                raise ValueError('Path to bigwiginfo is not set/provided...')

            if not os.path.isfile(self.pathTobigWigInfo):
                raise IOError('bigwiginfo is not found at {0} . Set or provide different path.'.format(self.pathTobigWigInfo))
        else:
            if not os.path.isfile(pathTobigWigInfo):
                raise IOError('bigwiginfo is not found at {0} . Set or provide different path.'.format(pathTobigWigInfo))

            # Add to config file
            if pathTobigWigInfo != config['Programs']['bigWigInfo']:
                self.logger.info('Adding bigwiginfo path in configuration...')
                updateConfig('Programs', 'bigWigInfo', self.pathTobigWigInfo)
                config['Programs']['bigWigInfo'] = self.pathTobigWigInfo

            self.pathTobigWigInfo = pathTobigWigInfo

    def _checkBigWigToWigProgram(self, pathTobigWigToWig):
        """ Check if bigWigToWig program is available or accessible.

        If program is not available in configuration file, the given
        path will be stored in the file after checking its accessibility.

        The path is stored in :attr:`gcMapExplorer.lib.genomicsDataHandler.EncodeDatasetsConverter.pathTobigWigToWig`

        Parameters
        ----------
        pathTobigWigToWig : str
            Path to bigWigToWig program

        """
        if pathTobigWigToWig is None or pathTobigWigToWig == config['Programs']['bigWigToWig']:
            self.pathTobigWigToWig = config['Programs']['bigWigToWig']
            if self.pathTobigWigToWig == 'None':
                raise ValueError('Path to bigWigToWig is not set/provided...')

            if not os.path.isfile(self.pathTobigWigToWig):
                raise IOError('bigWigToWig is not found at {0} . Set or provide different path.'.format(self.pathTobigWigToWig))
        else:
            if not os.path.isfile(pathTobigWigToWig):
                raise IOError('bigWigToWig is not found at {0} . Set or provide different path.'.format(pathTobigWigToWig))

            # Add to config file
            if pathTobigWigToWig != config['Programs']['bigWigToWig']:
                self.logger.info('Adding bigWigToWig path in configuration...')
                updateConfig('Programs', 'bigWigToWig', self.pathTobigWigToWig)
                config['Programs']['bigWigToWig'] = self.pathTobigWigToWig

            self.pathTobigWigToWig = pathTobigWigToWig

    def _getCreationDateOfFile(self, fileId):
        # This URL locates the ENCODE biosample with accession number ENCBS000AAA
        URL = r'https://www.encodeproject.org/file/{0}/?format=json'.format(fileId)

        # GET the object
        response = urlopen(URL)

        # Decode to string and save as json
        readableResponse = response.read().decode('utf8')
        response_json_dict = json.loads(readableResponse)

        return re.split('T', response_json_dict['date_created'])[0]

    def _readMetaDataChIPseq(self, row):
        data = dict()
        name = re.split('-', row['Experiment target'])[0]

        if 'signal p-value' in row['Output type']:
            data['type'] = 'signal-' + name
        elif 'fold change' in row['Output type']:
            data['type'] = 'fold-' + name
        else:
            return

        return data

    def _readMetaData_siRNA_RNAseq(self, row):
        data = dict()
        name = re.split('-', row['Experiment target'])[0]

        if re.search('signal', row['Output type']) is not None:

            if re.search('unique reads', row['Output type']) is not None:
                data['type'] = 'uniq-reads-signal' + '-' + name
            if re.search('all reads', row['Output type']) is not None:
                data['type'] = 'all-reads-signal' + '-' + name

            if 'type' not in data:
                data['type'] = 'signal' + '-' + name

        else:
            return

        return data

    def _readMetaDataRNAseq(self, row):
        data = dict()

        if re.search('signal', row['Output type']) is not None:

            date = self._getCreationDateOfFile(row['File accession'])
            self.logger.info(' Fetching creation date for file {0} : {1} '.format(row['File accession'], date))
            data['date'] = date

            if re.search('unique reads', row['Output type']) is not None:
                if re.search('plus strand', row['Output type']) is not None:
                    data['type'] = 'plus-uniq-reads' + '-' + date
                elif re.search('minus strand', row['Output type']) is not None:
                    data['type'] = 'minus-uniq-reads' + '-' + date
                else:
                    data['type'] = 'uniq-reads' + '-' + date

            if re.search('all reads', row['Output type']) is not None:
                if re.search('plus strand', row['Output type']) is not None:
                    data['type'] = 'plus-all-reads' + '-' + date
                elif re.search('minus strand', row['Output type']) is not None:
                    data['type'] = 'minus-all-reads' + '-' + date
                else:
                    data['type'] = 'all-reads' + '-' + date

            if 'type' not in data:
                data['type'] = 'signal' + '-' + date

        else:
            return

        return data

    def _readMetaDataDNaseSeq(self, row):
        data = dict()

        if re.search('signal', row['Output type']) is not None:

            date = self._getCreationDateOfFile(row['File accession'])
            self.logger.info(' Fetching creation date for file {0} : {1} '.format(row['File accession'], date))
            data['date'] = date

            if re.search('unique reads', row['Output type']) is not None:
                data['type'] = 'uniq-reads-signal' + '-' + date
            if re.search('raw', row['Output type']) is not None:
                data['type'] = 'raw-signal' + '-' + date
            if re.search('all reads', row['Output type']) is not None:
                data['type'] = 'all-reads-signal' + '-' + date

            if 'type' not in data:
                data['type'] = 'signal' + '-' + date

        else:
            return

        return data

    def downloadMetaData(self):
        """ Download the metadata file

        It downloads the metadata file and stored at
        :attr:`gcMapExplorer.lib.genomicsDataHandler.EncodeDatasetsConverter.metafile`

        """
        # Read metadata file url
        fin = open(self.inputFile, 'r')
        url = fin.readline()
        url = url.rstrip()
        fin.close()

        # Generate metadata file name
        name = 'gcx_' + util.getRandomName()+ '.tsv'
        self.metafile = os.path.join(self.workDir, name)

        # Download the file
        self.logger.info(' Downloading metadata...')
        downloadFile(url, self.metafile)
        self.logger.info('                        ... Finished.')

    def readMetaData(self):
        """ Read the metafile and extract the information

        It reads the metafile, filter the datasets according to assembly
        and file-format and make a list as dictionary.

        Each dictionary contains following field:
            * title : ``Experiment target``-``Experiment accession``-``File accession``-[``fold``/``signal``]
            * type : ``fold`` for fold change over control or ``signal`` for signal p-value.
            * url : URL to file
            * ``Experiment accession``
            * ``File accession``
            * ``Biological replicate(s)``

        """

        metaData = []

        with open(self.metafile) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                if row['File format'] == 'bigWig' and row['Assembly'] == self.assembely:

                    # Skipping if assay is not implemented
                    if row['Assay'] not in self._implementedAssays:
                        self.logger.info(' Assay [{0}] not implemented. Skipping it...'.format(row['Assay']))
                        continue

                    data = None
                    if row['Assay'] == 'ChIP-seq':
                        data = self._readMetaDataChIPseq(row)
                    if row['Assay'] == 'siRNA knockdown followed by RNA-seq':
                        data = self._readMetaData_siRNA_RNAseq(row)
                    if row['Assay'] == 'RNA-seq':
                        data = self._readMetaDataRNAseq(row)
                    if row['Assay'] in ['DNase-seq', 'FAIRE-seq'] :
                        data = self._readMetaDataDNaseSeq(row)

                    if data is not None:
                        data['title'] = data['type'] + '-' + row['Experiment accession']  + '-' +  row['File accession']
                        data['url'] = row['File download URL']
                        data['Experiment accession'] = row['Experiment accession']
                        data['File accession'] = row['File accession']
                        data['Assay'] = row['Assay']

                        # Sometime biological replicates are not provided
                        if not row['Biological replicate(s)']:
                            data['Biological replicate(s)'] = '1'
                        else:
                            data['Biological replicate(s)'] = row['Biological replicate(s)']

                        metaData.append(data)

        self.metaData = metaData

    def filterReplicates(self):
        """ It filters the metaData according to the replicates.
        If several replicates for a dataset are present, it only reads the
        dataset with combined replicates.

        In case if two replicates are present and combined replicates are not
        present, first replicate will be considered. Combining replicates are
        not yet implemented
        """
        experiments = dict()
        for i in range(len(self.metaData)):
            accession = self.metaData[i]['Experiment accession']
            etype = self.metaData[i]['type']
            if accession not in experiments:
                experiments[accession] = {etype : [i]}
            else:
                if etype not in experiments[accession]:
                    experiments[accession][etype] =  [i]
                else:
                    experiments[accession][etype] += [i]

        index = []
        for expID in experiments:
            for etype in experiments[expID]:
                count = 0
                currIdxes = []
                maxIdx = None
                length = 0
                for idx in experiments[expID][etype]:
                    currIdxes.append(idx)
                    count += 1
                    temp = re.split(',\s+', self.metaData[idx]['Biological replicate(s)'])
                    # In case if combined data is already present
                    if length < len(temp):
                        length = len(temp)
                        maxIdx = [ idx ]

                # In case if combined data is not present
                if length == 1 and count > 1:
                    maxIdx = currIdxes

                #print(maxIdx, self.metaData[maxIdx]['Biological replicate(s)'], self.metaData[maxIdx]['title'])
                index.append(maxIdx)

        self.indexes = index

    def _readFromCheckPoint(self):
        """ Read the titles from checkpoint file.
        """
        try:
            fin = open(self._checkPointFile, 'r')
        except:
            return

        json_dict = json.load( fin )
        if json_dict is not None:
            if 'titles' in json_dict:
                self._checkPoint = json_dict['titles']
        fin.close()

        return

    def _writeToCheckPoint(self, idx):
        """ Write to done titles to checkpoint file
        """

        haveInCP = False
        for midx in idx:
            if self.metaData[midx]['title'] in self._checkPoint:
                haveInCP = True

        if not haveInCP:
            self._checkPoint.append(self.metaData[idx[0]]['title'])
            try:
                fin = open(self._checkPointFile, 'w')
                json.dump({'titles' :  self._checkPoint }, fin, indent=4, separators=(',', ':') )
                fin.close()
            except IOError:
                pass
            except:
                raise

        return

    def saveAsH5(self, outDir, resolutions=None, coarsening_methods=None, compression='lzf', keep_original=False):
        """Download the files and convert to gcMapExplorer compatible hdf5 file.

        Name of output files:

            (1) For ChIP-seq assay:
                a. signal-<Experiment target>-<Experiment accession>-<File accessions.h
                b. fold-<Experiment target>-<Experiment accession>-<File accessions>.h5


            (2) For RNA-seq:
                a. uniq-reads-<date>-<Experiment accession>-<File accessions>.h
                b. plus-uniq-reads-<date>-<Experiment accession>-<File accessions>.h
                c. minus-uniq-reads-<date>-<Experiment accession>-<File accessions>.h
                d. all-reads-<date>-<Experiment accession>-<File accessions>.h5
                e. plus-all-reads-<date>-<Experiment accession>-<File accessions>.h5
                f. minus-all-reads-<date>-<Experiment accession>-<File accessions>.h5
                g. signal-<date>-<Experiment accession>-<File accession>.h5


            (3) For DNase-seq:
                a. uniq-reads-signal-<date>-<Experiment accession>-<File accessions>.h
                b. raw-signal-<date>-<Experiment accession>-<File accessions>.h
                c. all-reads-signal-<date>-<Experiment accession>-<File accessions>.h
                d. signal-<date>-<Experiment accession>-<File accessions>.h5


            (4) For siRNA + RNA-seq:
                a. uniq-reads-signal-<Experiment target>-<Experiment accession>-<File accessions>.h
                b. all-reads-signal-<Experiment target>-<Experiment accession>-<File accessions>.h
                c. signal-<Experiment target>-<Experiment accession>-<File accessions>.h5


        .. note:: Because downloading and conversion might take very long time,
                  it also generates a checkpoint file in the output directory.
                  Therefore, in case of crash or abrupt exit, the process can
                  be continued from the last file.


        Parameters
        ----------
        outDir : str
            Output directory where all files will be saved. Checkpoint file
            will be stored in same directory.
        resolutions : list of str
            Additional input resolutions other than these default resolutions:
            1kb', '2kb', '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb',
            '100kb', '160kb','200kb', '320kb', '500kb', '640kb',  and '1mb'.

            For Example: use ``resolutions=['25kb', '50kb', '75kb']`` to add
            additional 25kb, 50kb and 75kb resolution data.
        coarsening_methods : list of str
            Methods to coarse or downsample the data for converting from 1-base
            to coarser resolutions. Presently, five methods are implemented.

            * ``'min'``    -> Minimum value
            * ``'max'``    -> Maximum value
            * ``'amean'``  -> Arithmatic mean or average
            * ``'hmean'``  -> Harmonic mean
            * ``'gmean'``  -> Geometric mean
            * ``'median'`` -> Median

            In case of ``None``, all five methods will be considered. User may
            use only subset of these methods. For example:
            ``coarse_method=['max', 'amean']`` can be used for downsampling by
            only these two methods.
        compression : str
            data compression method in HDF5 file : ``lzf`` or ``gzip`` method.
        keep_original : bool
            Whether original data present in bigwig file should be incorporated in HDF5 file. This will significantly increase size of HDF5 file.

        """
        self.readMetaData()
        self.filterReplicates()

        # get working directory
        if self.workDir is None or self.workDir == config['Dirs']['workingdirectory']:
            tmpNumpArrayFiles = TempNumpyArrayFiles()
        else:
            tmpNumpArrayFiles = TempNumpyArrayFiles(workDir=self.workDir)

        self._checkPointFile = os.path.join(outDir, self._checkPointFile)
        self._readFromCheckPoint()


        for idx in self.indexes:

            # Skip those that are already converted
            haveInCP = False
            for midx in idx:
                if self.metaData[midx]['title'] in self._checkPoint:
                    haveInCP = True
            if haveInCP:
                self.logger.info(' {0} is already converted... Skipping...'.format(self.metaData[idx[0]]['title']))
                continue

            self._removeBigWigFiles()
            self._bigWigFiles = dict()
            for midx in idx:
                output = '{0}.bigWig' .format(self.metaData[midx]['File accession'])
                self._bigWigFiles[midx] = os.path.join(self.workDir, output)

            # Download the file
            downloadSuccess = True
            for midx in idx:
                self.logger.info(' Downloading file from  {0}...'.format(self.metaData[midx]['url']))
                if not downloadFile(self.metaData[midx]['url'], self._bigWigFiles[midx]):
                    self.logger.warning('Not able to download: {0} '.format(self.metaData[midx]['url']))
                    downloadSuccess = False
                    break
            if not downloadSuccess:
                continue
            self.logger.info('                ... Finished Downloading.')

            # Convert the file
            bigwigHandle = BigWigHandler(list(self._bigWigFiles.values()),
                                pathTobigWigToWig=self.pathTobigWigToWig,
                                pathTobigWigInfo=self.pathTobigWigInfo,
                                methodToCombine=self.methodToCombine,
                                workDir=self.workDir, maxEntryWrite=10000000)

            # Name of output file
            name = self.metaData[idx[0]]['type'] + '-' + self.metaData[idx[0]]['Experiment accession']
            for midx in idx:
                name += '-' +  self.metaData[midx]['File accession']
            outputH5 = os.path.join(outDir, name +'.h5')

            try:
                bigwigHandle.saveAsH5(outputH5, tmpNumpyArrayFiles=tmpNumpArrayFiles,
                                title=name,
                                resolutions=resolutions,
                                coarsening_methods=coarsening_methods,
                                compression=compression,
                                keep_original=keep_original)


            except (SystemExit, KeyboardInterrupt) as e:
                del tmpNumpArrayFiles
                raise e
            except Exception as e:
                self.logger.info(' Not able to convert {0} files.'.format(self._bigWigFiles.values()))
            finally:
                del bigwigHandle
                self._removeBigWigFiles()
                self._writeToCheckPoint(idx)

        del tmpNumpArrayFiles

    def _removeBigWigFiles(self):
        """ Remove downloaded bigwig file
        """
        # Remove bigwig file
        if self._bigWigFiles is None:
            return

        for key in self._bigWigFiles:
            self.logger.info(' Removing temporary bigwig file [{0}] '.format(self._bigWigFiles[key]))
            if os.path.isfile(self._bigWigFiles[key]):
                os.remove(self._bigWigFiles[key])

        self._bigWigFiles = None

class TextFileHandler:
    """To import a genomic data from column text file format

    It reads text file, make an full array for given shape and fills the missing place with zeros.
    These zeros could be later masked to perform any analysis.

    Example file format:

        ::

            15670000	0.2917373776435852
            15680000	0.2292359322309494
            15690000	0.023434270173311234
            15700000	0.06813383102416992
            15710000	0.13660947978496552
            15720000	0.17478400468826294
            15730000	0.20540907979011536
            .
            .
            .


    This class is used in browser to visualize the genomic data, which are directly imported from text file.

    Attributes
    ----------
    filename : str
        Input text file
    shape : int
        Size of array required to built
    binsize : int
        Size of bins expected in input file. If ``binsize = None``, binsize will be determined from the files, however, it is good to give
        expected binsize to check whether expected binsize match with binsize present in input text file.
    title : str
        Title of the input data
    workDir : str
        Directory where temporary files will be generated. If ``None``, defualt temporary directory of the respective OS will be used.
    data : numpy.ndarray or numpy.memmap
        One-dimensional array containing the data
    tmpNumpyFileName : str
        Name of temporary numpy memory-mapped file

    Parameters
    ----------
    filename : str
        Input text file
    shape : int
        Size of array required to built
    binsize : int
        Size of bins expected in input file. If ``binsize = None``, binsize will be determined from the files, however, it is good to give
        expected binsize to check whther expected binsize match with binsize present in input text file.
    title : str
        Title of the input data
    workDir : str
        Directory where temporary files will be generated. If ``None``, defualt temporary directory of the respective OS will be used.


    """

    def __init__(self, filename, shape, binsize=None, title=None, workDir=None):
        self.filename = filename
        self.shape = shape

        # Working and output directory
        if workDir is None:
            self.workDir = config['Dirs']['WorkingDirectory']
        else:
            self.workDir = workDir

        self.data = None
        self.binsize = None
        self.title = title

        self.tmpNumpyFileName = None

        self.logger = logging.getLogger('TextFileHandler')
        self.logger.setLevel(logging.INFO)

        binsizeInFile = self._getBinSize()
        if binsize is not None:
            if binsizeInFile != binsize:
                raise AssertionError ('Input binsize {0} does not match with the binsize {1} present in file.'.format(binsize, binsizeInFile))
        self.binsize = binsizeInFile

    def __del__(self):
        self._removeTempNumpyFile()

    def _removeTempNumpyFile(self):
        """Remove temporary numpy memory-mapped file
        """
        if self.tmpNumpyFileName is not None:
            del self.data
            if os.path.isfile(self.tmpNumpyFileName):
                os.remove(self.tmpNumpyFileName)

    def _generateTempNumpyFile(self):
        """Generate temporary numpy memory-mapped file
        """
        (fd, fname) = tempfile.mkstemp(suffix='.npy', prefix='gcx_np_', dir=self.workDir, text=False)
        os.close(fd)     # Close file, error in windows OS
        self.logger.info(' Generating temporary numpy array file [{0}] ...' .format(fname))

        size = self.shape + 1        # Be careful: added one to easiliy handle real locations, zeroth index is dummy, dont use zeroth location
        np.save(fname, np.zeros(size, dtype=np.float))
        self.data = np.load(fname, mmap_mode='r+')
        self.tmpNumpyFileName = fname

    def _getBinSize(self):
        """Determine binsize from input file
        """
        fin = open(self.filename, 'r')

        listBinSize = []
        count = 0
        prevLocation = None
        for line in fin:
            line = line.lstrip().rstrip()

            if not line.strip():    continue

            temp = re.split('\s+', line)

            newLocation = int(temp[0])
            if prevLocation is not None:
                listBinSize.append( abs(newLocation - prevLocation) )
            prevLocation = newLocation

            if count > 100:
                break

            count = count + 1

        fin.close()

        return np.amin( listBinSize )

    def readData(self):
        """Read data from input file

        Read data from input file and store in :attr:`TextFileHandler.data` as one-dimensional array
        """
        if self.shape > 100000:
            self._generateTempNumpyFile()
        else:
            size = self.shape + 1
            self.data = np.zeros(size, dtype=np.float)

        fin = open(self.filename, 'r')

        try:
            for line in fin:
                line = line.lstrip().rstrip()

                if not line.strip():    continue

                temp = re.split('\s+', line)

                idx = int( int(temp[0])/self.binsize )

                try:
                    self.data[idx] = float(temp[1])
                except IndexError:
                    break
                    pass

        except (SystemExit, KeyboardInterrupt) as e:
            self._removeTempNumpyFile()
            raise e

        fin.close()
