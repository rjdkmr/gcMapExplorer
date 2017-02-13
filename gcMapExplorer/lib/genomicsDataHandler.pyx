#!/usr/bin/env python3
#
# Author: Rajendra Kumar
#
# This file is part of gcMapExplorer
# Copyright (C) 2016  Rajendra Kumar, Ludvig Lizana, Per Stenberg
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
import re, os, sys
import copy
import shlex, subprocess
import tempfile
import logging
import urllib
import shutil
import h5py

from gcMapExplorer.config import getConfig

# get configuration
config = getConfig()

from . import ccmap as cmp


def check_resolution_list(resolutions):
    new_resolutions = [ '1kb', '2kb', '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb', '100kb', '160kb', '200kb', '320kb', '500kb', '640kb', '1mb']

    if resolutions is not None:
        if isinstance(resolutions, list):
            new_resolutions = list( set(new_resolutions + resolutions) )
        else:
            raise ValueError(' Resolutions should be of list type.')

    return new_resolutions

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
        try:
            self.arrays[key]._mmap.close()
        except:
            pass

        try:
            os.remove(self.files[key])
        except:
            pass

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
            self.logger.info('{0} is already present in temporary numpy file list, skipping ...' .format(key))
        else:
            if regenerate:
                self._removeTempNumpyFile(key)

            (fd, fname) = tempfile.mkstemp(suffix='.npy', prefix=key+'_', dir=self.workDir, text=False)
            os.close(fd)     # Close file, error in windows OS
            self.logger.info(' Generating temporary numpy array file [{0}] for {1} ...' .format(fname, key))

            size = self.chromSizeInfo[key] + 1        # Be careful: added one to easiliy handle real locations, zeroth index is dummy, dont use zeroth location
            np.save(fname, np.zeros(size, dtype=np.float))
            npy = np.load(fname, mmap_mode='r+')

            self.arrays[key] = npy
            self.files[key] = fname

            # print(npy.shape)
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
        try:
            self.hdf5.close()
            self.logger.info(' Closed {0} ...' .format(self.filename))
            self.hdf5 = None
        except:
            pass

    def open(self):
        """ open hdf5 file
        """
        self.hdf5 = h5py.File(self.filename)

        # Fail-safe mechanism for title in case of either new file or already opened file
        if 'title' not in self.hdf5.attrs:
            if self.title is None:
                self.hdf5.attrs['title'] = "Genome Data"
                self.title = "Genome Data"
            else:
                self.hdf5.attrs['title'] = self.title
        else:
            if self.title is None:
                self.title = self.hdf5.attrs['title']

        self.logger.info(' Opened {0} ...' .format(self.filename))

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

    def getResolutionList(self, chrom):
        """ To get all resolutions for given chromosome from hdf5 file

        Parameters
        ----------
        chrom : str
            chromosome name

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
                resolutionList.append(r)
        else:
            raise KeyError(' Chromosome [{0}] not found in [{1}] file...' .format(chrom, self.filename))

        return resolutionList

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
    """ To handle bigWig files and to convert it to hdf5 file

    This class can be used to convert bigWig file to hdf5 file.
    It can also be used to combine several bigWig files that are originated from replicated experiments.

    Attributes
    ----------
    bigWigFileNames : str or list[str]
        List of bigWig file names including path
    pathTobigWigToWig : str
        Path to bigWigToWig program
    pathTobigWigInfo : str
        Path to bigWigInfo program
    WigFileNames : str
        List of Wig file names, either autmoatically generated or given by user
    wigHandle : WigHandler
        WigHandler instance to parse Wig file and save data as hdf5 file
    chromSizeInfo : dict
        A dictionary containing chromosome size information
    methodToCombine : str
        method to combine bigWig/Wig files, Presently, accepted keywords are: ``mean``, ``min`` and ``max``
    chromName : str
        Name of chromosome, which has to be extracted from bigWig file.
    resolutions : list[str]
        List of input resolutions for which data will be stored in the hdf5 file
    maxEntryWrite : int
        Number of lines read from Wig file at an instant, after this, data is dumped in temporary numpy array file

    Parameters
    ----------
    filenames : str or list[str]
        A bigWig file or list of bigWig files including path
    pathTobigWigToWig: str
        Path to bigWigToWig program
    pathTobigWigInfo : str
        Path to bigWigInfo program
    resolutions : list[str] (optional)
        List of input resolutions for which data will be stored inside hdf5 file. If its value is ``None``, a new default list is generated by :meth:`check_resolution_list`
    methodToCombine : str
        method to combine bigWig/Wig files, Presently, accepted keywords are: ``mean``, ``min`` and ``max``
    chromName : str
        Name of chromosome, which has to be extracted from bigWig file.
    maxEntryWrite : int
        Number of lines read from Wig file at an instant, after this, data is dumped in temporary numpy array file. To reduce memory (RAM) occupancy,
        reduce this number because large numbers need large RAM.

    """

    def __init__(self, filenames, pathTobigWigToWig, pathTobigWigInfo, resolutions=None, chromName=None, methodToCombine='mean', workDir=None, maxEntryWrite=10000000):
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

        # Convert to list
        if not isinstance(filenames, list):
            self.bigWigFileNames = [ filenames ]

        if not (methodToCombine == 'mean' or methodToCombine == 'max' or methodToCombine == 'min'):
            raise NotImplementedError(' Method [{0}] to combine bigwig file not implemented. Use: \'mean\', \'max\' or \'min\' ' .format(methodToCombine))
        self.methodToCombine = methodToCombine

        self.resolutions = check_resolution_list(resolutions)

    def __del__(self):
        self._removeTempWigFiles()

    def _removeTempWigFiles(self):
        if self.wigToDel:
            for f in self.WigFileNames:
                self.logger.info(' Removing temporary Wig file [{0}] .' .format(f))
                try:
                    os.remove(f)
                except:
                    pass

    def _getBigWigInfo(self, filename):
        """Base method to Retrieve chromosome names and their sizes

        * Chromosome size information is stored for a given bigWig file. If size of chromosome is already present in dictionary, largest size is stored in dictionary.

        * Use :meth:`BigWigHandler.getBigWigInfo` to automatically retrieve chromosome size information from all bigWig files.

        .. warning::
            **Private method.** Use it at your own risk. It is used internally in :meth:`BigWigHandler.getBigWigInfo`

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

        BigWigInfo progrom is executed on all listed bigWig files and chromosomes name with respecitve size is stored in :attr:`BigWigHandler.chromSizeInfo` variable.
        From the several listed bigWig files, only largest size of chromosomes are considered.

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

        Use :meth:`BigWigHandler.bigWigtoWig` to automatically convert all bigWig files to Wig files.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`BigWigHandler.bigWigtoWig`


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
            self.logger.error(' Not able to convert {0} to bigWig format by command: [ {1} ]. See details above.' .format(self.bigWigFileNames), cmd)
        self.logger.info('         ... Finished!!!' .format(cmd))

    def bigWigtoWig(self, outfilenames=None):
        """ To generate Wig file

        It uses bigWigToWig program to convert bigWig to Wig file.
        Newly generated Wig files name are listed in :attr:`BigWigHandler.WigFileNames`.

        Parameters
        ----------
        outfilenames : str or list of strip
            List of Wig file names. If ``None``, names are automatically generated, files are temporarily created and after execution, all files are deleted.

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
                tName = os.path.join(self.workDir, '{0}_{1}.wig.tmp' .format(cmp.name_generator(), i+1) )
                self.WigFileNames.append(tName)
                self.wigToDel = True

        for i in range(len(self.bigWigFileNames)):
            self._bigWigtoWig(self.bigWigFileNames[i], self.WigFileNames[i])

    def saveToH5(self, filename, tmpNumpyArrayFiles=None, title=None, compression='lzf', keep_original=False):
        """Save data to hdf5 file.

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
        compression : str
            data compression method in HDF5 file : ``lzf`` or ``gzip`` method.
        keep_original : bool
            Whether original data present in bigwig file should be incorporated in HDF5 file. This will significantly increase size of HDF5 file.

        Examples
        --------
        .. code-block:: python

            from hiCMapAnalyze import genomicsDataHandler as gdh

            # start BigWigHandler to combine and convert two bigWig files
            bigwig = gdh.BigWigHandler(['first.bigWig', 'second.bigWig'], './bigWigToWig', './bigWigInfo', maxEntryWrite=1000000)

            # Save hdf5 file
            bigwig.saveToH5('converted.h5')


        """
        if self.chromSizeInfo is None:
            self.logger.info(' No information of chromsomes and their sizes from BigWig File!!! Need to get information...')
            self.getBigWigInfo()

        if self.WigFileNames is None:
            self.logger.info(' No Wig File!!! Need to generate Wig File...')
            self.bigWigtoWig()

        self.wigHandle = WigHandler(self.WigFileNames, self.chromSizeInfo, resolutions=self.resolutions, tmpNumpyArrayFiles=tmpNumpyArrayFiles,
                                     methodToCombine=self.methodToCombine, workDir=self.workDir, maxEntryWrite=self.maxEntryWrite)
        self.wigHandle.convertWigToH5(filename, title=title, compression=compression, keep_original=keep_original)


class WigHandler:
    """To convert Wig files to hdf5 file

    It parses wig files and save all data to a hdf5 file for given resolutions.

    Attributes
    ----------
    WigFileNames : list[str]
        List of Wig files
    chromSizeInfo : dict
        A dictionary containing chromosome size information.
    methodToCombine : str
        method to combine bigWig/Wig files, Presently, accepted keywords are: ``mean``, ``min`` and ``max``
    resolutions : list of str
        List of input resolutions for which data will be stored inside hdf5 file
    tmpNumpyArrayFiles : :class:`TempNumpyArrayFiles`
        This :class:`TempNumpyArrayFiles` instance stores the temporary numpy array files information.
    isWigParsed : bool
        Whether Wig files are already parsed.
    maxEntryWrite : int
        Number of lines read from Wig file at an instant, after this, data is dumped in temporary numpy array file

    Parameters
    ----------
    filenames : str or list(str)
        A Wig file or list of Wig files including path.
    chromSizeInfo : dict
        A dictionary containing chromosome size information. Generated by :meth:`BigWigHandler.getBigWigInfo`.
    resolutions : list[str]
        List of input resolutions for which data will be stored inside hdf5 file
    tmpNumpyArrayFiles : :class:`TempNumpyArrayFiles`
        This :class:`TempNumpyArrayFiles` instance stores the temporary numpy array files information.
    methodToCombine : str
        method to combine bigWig/Wig files, Presently, accepted keywords are: ``mean``, ``min`` and ``max``
    maxEntryWrite : int
        Number of lines read from Wig file at an instant, after this, data is dumped in temporary numpy array file. To reduce memory (RAM) occupancy,
        reduce this number because large numbers need large RAM.


    """
    def __init__(self, filenames, chromSizeInfo=None, resolutions=None, tmpNumpyArrayFiles=None, methodToCombine='mean', workDir=None, maxEntryWrite=10000000):
        self.WigFileNames = filenames

        self.chromSizeInfo = chromSizeInfo                  # dictionary for chromosome sizes, extracted from bigWigInfo program
        self.tmpNumpyArrayFiles = None                        # List of temporary file to store array
        self.isWigParsed  = False

        # Working and output directory
        if workDir is None:
            self.workDir = config['Dirs']['WorkingDirectory']
        else:
            self.workDir = workDir

        self.resolutions = check_resolution_list(resolutions)

        self.logger = logging.getLogger(__name__+'.WigHandler')
        self.logger.setLevel(logging.INFO)

        self.maxEntryWrite = maxEntryWrite # maximum count before writing to output file will be attempted

        # Convert to list
        if not isinstance(filenames, list):
            self.WigFileNames = [ filenames ]

        if not (methodToCombine == 'mean' or methodToCombine == 'max' or methodToCombine == 'min'):
            raise NotImplementedError(' Method [{0}] to combine wig file not implemented. Use: \'mean\', \'max\' or \'min\' ' .format(methodToCombine))
        self.methodToCombine = methodToCombine

        if tmpNumpyArrayFiles is None:
            self.tmpNumpyArrayFiles = TempNumpyArrayFiles(workDir=self.workDir)
            self.tmpNumpyArrayFiles.chromSizeInfo = copy.deepcopy(chromSizeInfo)
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

    def _PerformDataCoarsing(self, Chrom, resolution, coarse_method):
        """Base method to perform Data coarsing.

        This method read temporary Numpy array files and perform data coarsing using the given input method.

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`WigHandler._StoreInHdf5File`.

        Parameters
        ----------
        Chrom : str
            Chromosome name

        resolution : str
            resolution in word.

        coarse_method : str
            Name of method to use for data coarsing. Accepted keywords: min, max, median, amean, gmean and hmean.


        """

        output = []
        binsize = cmp.resolutionToBinsize(resolution)
        size = self.chromSizeInfo[Chrom] + 1

        for i in range(1, size, binsize):
            tmpx = None
            if i+binsize >= size:
                tmpx = self.tmpNumpyArrayFiles.arrays[Chrom][i : size]
            else:
                tmpx =  self.tmpNumpyArrayFiles.arrays[Chrom][i : i+binsize]

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

        for n in range(len(location_list)):
            if self.methodToCombine == 'mean':
                self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list[n]] = (value_list[n]  + self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list[n]]) / 2
            elif self.methodToCombine == 'max':
                self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list[n]] = max(value_list[n], self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list[n]])
            elif self.methodToCombine == 'min':
                self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list[n]] = min(value_list[n], self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list[n]])
            else:
                self.tmpNumpyArrayFiles.arrays[ChromTitle][location_list[n]] = value_list[n]

        self.tmpNumpyArrayFiles.arrays[ChromTitle].flush()


    def _StoreInHdf5File(self, hdf5Out, title, compression='lzf', keep_original=False):
        """ Base method to perform data coarsing and to store coarsed data in hdf5 file

        .. warning::
            **Private method**. Use it at your own risk. It is used internally in :meth:`WigHandler._parseWig`.

        Parameters
        ----------
        hdf5Out : str or :class:`HDF5Handler`
            Name of output hdf5 file or instance of :class:`HDF5Handler`
        title : str
            Title of data
        compression : str
            data compression method in HDF5 file : ``lzf`` or ``gzip`` method.
        keep_original : bool
            Whether original data present in bigwig file should be incorporated in HDF5 file. This will significantly increase size of HDF5 file.

        """
        # In cases when hdf5 filename or hdf5handler object is given as input
        # Initialize new hdf5handler in case file name is given
        if not isinstance(hdf5Out, HDF5Handler):
            hdf5Out = HDF5Handler(hdf5Out, title=title)

        # Writing to hdf5 aur appending to dictionary to return
        self.logger.info(' Writing output file [{0}] ...' .format(hdf5Out.filename))

        for Chrom in self.chromSizeInfo:
            for res in self.resolutions:
                if keep_original:
                    hdf5Out.addDataByArray(Chrom, '1b', 'Original', self.tmpNumpyArrayFiles.arrays[Chrom][:], compression='gzip')
                for coarsing_method in ['min', 'max', 'amean', 'hmean', 'gmean', 'median']:
                    hdf5Out.addDataByArray(Chrom, res, coarsing_method, self._PerformDataCoarsing(Chrom, res, coarsing_method), compression=compression )

        self.logger.info(' \t\t ... Finished writing output file [{0}] ' .format(hdf5Out.filename))
        hdf5Out.close()


    def parseWig(self):
        """ To parse Wig files

        This method parses all Wig files listed in :attr:`WigHandler.WigFileNames`.
        The extracted data is further stored in temporary numpy array files of respective chromosome.
        These numpy array files can be used either for data coarsing or for further analysis.


        * **To save as hdf5:** Use :meth:`WigHandler.convertWigToH5`.
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

        This method parses a Wig file and extrated data are copied in temporary numpy array files.

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

        fin = open(wigFileName, 'r')

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

            # Logging the information
            if PreviousChromTitle != ChromTitle:
                if PreviousChromTitle != 'dummy':
                    self.logger.info('         ... Finished reading and processing for {0} ' .format(PreviousChromTitle))
                self.logger.info(' Identified \'{0}\' format in WIG file; Reading and processing for {1} ... ' .format(ftype, ChromTitle))

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
        self._FillDataInNumpyArrayFile(ChromTitle, location_list, value_list)
        self.logger.info('         ... Finished reading and processing for {0} ' .format(ChromTitle))


    def convertWigToH5(self, hdf5Out, title=None, compression='lzf', keep_original=False):
        """To convert Wig files to hdf5 file

        Parameters
        ----------
        hdf5Out : class:`HDF5Handler` or str
            Output hdf5 file name or class:`HDF5Handler` instance
        title : str
            Title of the data
        compression : str
            data compression method in HDF5 file : ``lzf`` or ``gzip`` method.
        keep_original : bool
            Whether original data present in bigwig file should be incorporated in HDF5 file. This will significantly increase size of HDF5 file.

        """
        if not self.isWigParsed:
            self.parseWig()

        # Storing data in hdf5 file
        self._StoreInHdf5File(hdf5Out, title, compression=compression, keep_original=keep_original)

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
            try:
                os.remove(self.tmpNumpyFileName)
            except:
                pass

    def _generateTempNumpyFile(self):
        """Generate temporary numpy memory-mapped file
        """
        (fd, fname) = tempfile.mkstemp(suffix='.npy', prefix='np_', dir=self.workDir, text=False)
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
