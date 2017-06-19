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

import re, os
import logging
import tempfile
import numpy as np
import tarfile
import zipfile

from gcMapExplorer.config import getConfig

# get configuration
config = getConfig()

from . import ccmap as cmp
from . import gcmap as gmp
from . import util

class CooMatrixHandler:
    """To import ccmap from files similar to sparse matrix in Coordinate (COO) format

    Two types of coordinates are accepted:
        * with ``coordinate='real'`` pair of absolute binned locations on chromosome
        * with ``coordinate='index'`` row and column index of matrix. index **should** always start from zero for absoulte beginning of chromosome.
          e.g. for 10kb, 0-10000 should have index of zero, 10000-20000 have index of one. If this is file format, resolution should be provided explicitly.

    .. warning::
        Input file should contain matrix for **only one** chromosome.

    Following file format can be read as a text file, where first and second column is location on chromosome and third column is the value:
    ::

        20000000	20000000	2692.0
        20000000	20100000	885.0
        20100000	20100000	6493.0
        20000000	20200000	15.0
        20100000	20200000	52.0
        20200000	20200000	2.0
        20000000	20300000	18.0
        20100000	20300000	40.0
        .
        .
        .
        .
        .
        .


    To Instantiate this class, three scenarios are possible:
        * ``Text in Archive``: Both a list of input files and a compressed file is given. It means look for input files in compressed file.
        * ``Text``: Only a list of input files is given. It means, read data directly from input file.
        * ``Archive``: Only a compressed file is given. It means, read all files present in compressed files.

    Parameters
    ----------
    inputFiles : str or list
        name of a input file or list of input files
    inputCompressedFile : str
        name of input tar archive
    mapType : str
    	Type of HiC map

    		* ``intra``: for intra-chromosomal map
    		* ``inter``: for inter-chromosomal map
    resolution : str
        resolution of HiC map. Example: '100kb', '10kb' or '25kb' etc. If ``None``, resolution will be determined from the input file.
    workDir : str
        Directory where temporary files will be stored.


    Attributes
    ----------
    inputFileList : list
        List of input files. It could be ``None`` when not provided.
    inputType : str
        Type of input. Three types: ``Text``, ``Text in Archive`` and ``Archive`` are determined from user input.
    inputCompressedFile : str
        Name of input compressed file. It could be ``None`` when not provided.
    outputFileList : str
        List of Output files
    compressType : str
        Format of compressed file. It could be etiher ``.tar`` or ``.zip`` or ``None``
    compressHandle : ZipFile or TarFile
        An object to handle compressed file. It could be `ZipFile <https://docs.python.org/3/library/zipfile.html#zipfile-objects>`_ or
        `TarFile <https://docs.python.org/3/library/tarfile.html#tarfile-objects>`_ instance depending on compressed format.
    mapType : str
    	Types of HiC map
    	* ``intra``: for intra-chromosomal map
    	* ``inter``: for inter-chromosomal map
    coordinate : str
        Coordinate type in input text file. It could be either ``real`` for real locations or ``index`` for rows and column indices.
    resolution : str
        resolution of HiC map. Example: '100kb', '10kb' or '25kb' etc. If ``None``, resolution will be determined from the input file.

        If ``coordinate='index'``, resolution is essential for further processing.
    workDir : str
        Directory where temporary files will be stored. If not provided, directory name will be taken from configuration file.

    """

    def __init__(self, inputFiles=None, inputCompressedFile=None, mapType='intra', resolution=None, coordinate='real', workDir=None, logHandler=None):

        # Checking both text and compressed input files
        self.inputType = self._checkInputType(inputFiles, inputCompressedFile)

        # Initial processing for archive file
        self.inputCompressedFile = inputCompressedFile
        self.compressType = self._checkCompressedFile()

        # Make a input file list if user has provided one
        self.inputFileList = None
        if inputFiles is not None:
            self.inputFileList = []
            self._addInputTextFiles(inputFiles)

        # If user has not provided input files, read file list from compressed file and return handle object of compressed file
        self.compressHandle = None
        self._extractCompressedFile()

        # Working and output directory
        if workDir is not None:
            self.workDir = workDir
        else:
            self.workDir = config['Dirs']['WorkingDirectory']

        # map type
        self.mapType = mapType

        # Checking corrdinate type: index or real location
        if coordinate == 'real':
            pass
        elif coordinate == 'index':
            if resolution is None:
                raise ValueError ("Please provide resolution because coordinate values are indexes of rows and columns")
        else:
            raise ValueError (" Keyword [{0}] is not allowed. Use 'real' or 'index'." .format(coordinate))
        self.resolution = resolution
        self.coordinate = coordinate

        # Output file list
        self.outputFileList = None

        # Labels
        self.xlabels = None
        self.ylabels = None

        # logger
        self.logHandler = logHandler
        self.logger = logging.getLogger('CooMatrixHandler')
        if logHandler is not None:
            self.logger.propagate = False
            self.logger.addHandler(logHandler)
        self.logger.setLevel(logging.INFO)

    def __del__(self):
        if self.logHandler is not None:
            self.logger.removeHandler( self.logHandler )

    def _addInputTextFiles(self, inputFiles):
        """ Add Input Files

        It add input files and also convert to list

        Parameters
        ----------
        inputFiles : str or list
            name of a input file or list of input files

        """
        if isinstance(inputFiles, list):
            self.inputFileList = self.inputFileList + inputFiles
        else:
            self.inputFileList.append(inputFiles)

    def _checkInputType(self, inputFiles, compressedFile):
        """To check input type given by user

        There are three possible scenarios:
            * ``Text in Archive``: Both a list of input files and a compressed file is given. It means look for input files in compressed file.
            * ``Text``: Only a list of input files is given. It means, read data directly from input file.
            * ``Archive``: Only a compressed file is given. It means, read all files present in compressed files.

        Parameters
        ----------
        inputFiles : str or list
            name of a input file or list of input files.
        compressedFile : str
            Name of compressed file

        Raises
        ------
        ValueError
            If both input files and compressedFile is ``None``

        """
        inputType = None
        # Three porabable type of inputs 'Text' 'Text in Archive' 'Archive'
        if inputFiles is not None and compressedFile is not None:
            inputType = 'Text in Archive'
        elif inputFiles is not None and compressedFile is None:
            inputType = 'Text'
        elif inputFiles is None and compressedFile is not None:
            inputType = 'Archive'
        else:
            raise ValueError ("Both inputFiles and inputCompressedFile cannot be None.")

        return inputType

    def _checkCompressedFile(self):
        """Determine format of compressed file.

        Returns
        -------
        compressType : str
            One of the two possible keywords (``tar`` and ``zip``) or ``None``.

        Raises
        ------
        NotImplementedError
            If format is not ``zip``, ``tar.gz`` or ``tar.bz2``.

        """
        if self.inputCompressedFile is None:  return
        compressType = None
        file_extension = os.path.splitext(self.inputCompressedFile)[1]
        if file_extension == '.zip':
            compressType = 'zip'
        elif file_extension == '.gz' or file_extension == '.tar' or file_extension == '.bz2':
            compressType = 'tar'
        else:
            raise NotImplementedError (" [{0}] compression format is not implemented. Please use '.zip' or '.tar.gz' or '.tar.bz2' compression format." .format(file_extension))

        return compressType

    def _extractCompressedFile(self):
        """Generate input file list from compressed file
        """
        if self.inputCompressedFile is None:    return

        if self.compressType == 'zip':
            self.compressHandle = zipfile.ZipFile(self.inputCompressedFile, 'r')
            if self.inputFileList is None:
                self.inputFileList = self.compressHandle.namelist()
        else:
            self.compressHandle = tarfile.open(self.inputCompressedFile, 'r')
            if self.inputFileList is None:
                self.inputFileList = self.compressHandle.getnames()

    def _read_a_map_text_file(self, filename):
        """ Generate CCMAP object from a text file with column format --- i  j  value

        To generate new CCMAP oject from a input text file. The file should contain three column where first and second column indicate
        location on the chromosome and third column indicate the contact frequency

        Parameters
        ----------
        filename : str
        	Name of input text file.

        Returns
        -------
        ccMapObj : :class:`hiCMapAnalyze.HiCMapMain.CCMAP`
        	A CCMAP object

        """

        fin = open(filename, 'r')

        i, j, value = [], [], []

        print('  ')
        self.logger.info(' Reading file: [{0}]... ' .format(filename))

        binsize = None
        if self.coordinate != 'real':
            binsize = util.resolutionToBinsize(self.resolution)

        for line in fin:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            temp = line.split()
            if self.coordinate == 'real':
                i.append(int(temp[0]))
                j.append(int(temp[1]))
            else:
                i.append(int(temp[0])*binsize)
                j.append(int(temp[1])*binsize)

            value.append(float(temp[2]))

        self.logger.info('     ... Finished reading file: [{0}] ' .format(filename))

        ccMapObj = gen_map_from_locations_value(i, j, value, resolution=self.resolution, workDir=self.workDir, mapType=self.mapType)
        self.resolution = util.binsizeToResolution(ccMapObj.binsize)

        del i
        del j
        del value

        return ccMapObj

    def _read_map_from_compressed_file(self, mapfile):
        """ Generate CCMAP object from a file, which is inside a tar archive

        If HiC map file is inside a tar archive, this method generate CCMAP object by directly reading the map file without extracting the tar archive.

        Parameters
        ----------
        mapfile : str
        	Name of input map file with full relative path inside the tar archive

        Returns
        -------
        ccMapObj : ``None`` or :class:`hiCMapAnalyze.HiCMapMain.CCMAP`
        	A CCMAP object is retured. However, if input file is not found inside tar archive, ``None`` is returned.

        """

        # Opening tar file
        # tar = tarfile.open(self.inputCompressedFile, "r")

        print('  ')
        self.logger.info(' Extracting-Reading [{0}] from [{1}]... ' .format(mapfile, self.inputCompressedFile))

        try:
            # extract the specific file to file stream
            if self.compressType == 'tar':
                fin = self.compressHandle.extractfile(mapfile)
            else:
                fin = self.compressHandle.open(mapfile)
        except KeyError:
            self.logger.warning(" [{0}] not found in [{1}]... \n\t\t...Skipping this file..." .format(mapfile, self.inputCompressedFile))
            return


        i, j, value = [], [], []

        binsize = None
        if self.coordinate != 'real':
            binsize = util.resolutionToBinsize(self.resolution)

        while (1):

            line = fin.readline().decode()

            if not line :
            	break

            line = line.rstrip('\n')

            if not line.strip():
            	continue

            temp = line.split()
            if self.coordinate == 'real':
                i.append(int(temp[0]))
                j.append(int(temp[1]))
            else:
                i.append(int(temp[0])*binsize)
                j.append(int(temp[1])*binsize)
            value.append(float(temp[2]))

        self.logger.info('   ...Finished extracting and reading [{0}]' .format(mapfile))

        ccMapObj = gen_map_from_locations_value(i, j, value, resolution=self.resolution, workDir=self.workDir, mapType=self.mapType)
        self.resolution = util.binsizeToResolution(ccMapObj.binsize)

        del i
        del j
        del value

        return ccMapObj

    def setOutputFileList(self, outputFiles):
        """To set list of output files

        Parameters
        ----------
        outputFiles : str or list
            Name of a output file or list of output files. For each input file,
            a respective output file will be generated, therefore, number of
            input and output files should match.

        """
        if not isinstance(outputFiles, list):
            self.outputFileList = [ outputFiles ]
        else:
            self.outputFileList = outputFiles

        if len(self.outputFileList) != len(self.inputFileList):
            raise AssertionError (" Number of input files [{0}] does not match"
                    + "with output files [{1}]. " .format(
                        len(self.inputFileList), len(self.outputFileList)) )

    def setLabels(self, xlabels, ylabels):
        """To set xlabels and ylabels for contact maps

        xlabel and ylabel act as a title of data along X-axis
        and Y-axis respectively.

        Parameters
        ----------
        xlabels : str or list
            Name of the data along X-axis or list of names of the data for respective input files.
        ylabels : str or list
            Name of the data along y-axis or list of names of the data for respective input files.

        """
        if not isinstance(xlabels, list):
            self.xlabels = [ xlabels ]
        else:
            self.xlabels = xlabels

        if not isinstance(xlabels, list):
            self.ylabels = [ ylabels ]
        else:
            self.ylabels = ylabels

        if len(self.xlabels) != len(self.inputFileList):
            raise AssertionError (" Number of input files [{0}] does not match"
                                    + "with xlabels [{1}]. " .format(
                                    len(self.inputFileList),
                                    len(self.xlabels)) )

        if len(self.ylabels) != len(self.inputFileList):
            raise AssertionError (" Number of input files [{0}] does not match"
                                    + "with ylabels [{1}]. " .format(
                                    len(self.inputFileList),
                                    len(self.ylabels)) )

    def save_ccmaps(self, outputFiles=None, xlabels=None, ylabels=None,
                    compress=True):
        """To Save all Hi-C maps

        This function reads input files one by one and save it as a ``.ccmap``
        file.

        Parameters
        ----------
        outputFiles : str or list
            Name of a output file or list of output files. For each input file,
            a respective output file will be generated, therefore, number of
            input and output files should match.
        xlabels : str or list
            Name of the data along X-axis or list of names of the data for
            respective input files.
        ylabels : str or list
            Name of the data along y-axis or list of names of the data for
            respective input files. if it is ``None``, ylabels will be same as
            xlabels.
        compress : bool
        	If ``True``, numpy array (matrix) file will be compressed to reduce
            storage memory.

        """

        if outputFiles is not None:
            self.setOutputFileList(outputFiles)

        if outputFiles is None and self.outputFileList is None:
            raise ValueError('Output file names not provided.')

        # xlabels and ylabels handling
        if xlabels is not None and ylabels is not None:
            self.setLabels(xlabels, ylabels)

        if xlabels is not None and ylabels is None:
            self.setLabels(xlabels, xlabels)

        if xlabels is None and self.xlabels is None:
            raise ValueError('xlabels not provided.')

        for i in range(len(self.inputFileList)):
            ccmap = None

            try:
                if self.inputType == 'Text in Archive' \
                        or self.inputType == 'Archive':
                    ccmap = self._read_map_from_compressed_file(
                                    self.inputFileList[i])
                else:
                    ccmap = self._read_a_map_text_file(self.inputFileList[i])

                if ccmap is not None:
                    ccmap.xlabel = self.xlabels[i]
                    ccmap.ylabel = self.ylabels[i]

                    cmp.save_ccmap(ccmap, self.outputFileList[i],
                                    compress=compress)
                    del ccmap
                else:
                    self.logger.warning("   Not able to read [{0}] ... "
                                + "Skipped!!" .format(self.inputFileList[i]))

            except (KeyboardInterrupt, SystemExit) as e:
                if ccmap is not None:
                    del ccmap
                raise e

    def save_gcmap(self, outputFile, xlabels=None, ylabels=None,
                    coarsingMethod='sum', compression='lzf'):
        """To Save all Hi-C maps as a gcmap file

        This function reads input files one by one and save it as a ``.gcmap``
        file.

        Parameters
        ----------
        outputFile : str
            Name of a output gcmap file.
        xlabels : str or list
            Name of the data along X-axis or list of names of the data for
            respective input files.
        ylabels : str or list
            Name of the data along y-axis or list of names of the data for
            respective input files. if it is ``None``, ylabels will be same as
            xlabels.
        coarsingMethod : str
            Method of downsampling. Three accepted methods are ``sum``: sum all
            values, ``mean``: Average of all values and ``max``: Maximum of
            all values.
        compression : str
            Compression method. Presently allowed : ``lzf`` for LZF compression
            and ``gzip`` for GZIP compression.

        """

        # xlabels and ylabels handling
        if xlabels is not None and ylabels is not None:
            self.setLabels(xlabels, ylabels)

        if xlabels is not None and ylabels is None:
            self.setLabels(xlabels, xlabels)

        if xlabels is None and self.xlabels is None:
            raise ValueError('xlabels are not provided.')

        for i in range(len(self.inputFileList)):
            ccmap = None

            try:
                if self.inputType == 'Text in Archive' \
                        or self.inputType == 'Archive':
                    ccmap = self._read_map_from_compressed_file(
                                    self.inputFileList[i])
                else:
                    ccmap = self._read_a_map_text_file(self.inputFileList[i])

                if ccmap is not None:
                    ccmap.xlabel = self.xlabels[i]
                    ccmap.ylabel = self.ylabels[i]
                    gmp.addCCMap2GCMap(ccmap, outputFile,
                                       compression = compression,
                                       generateCoarse = True,
                                       coarsingMethod = coarsingMethod)

                    del ccmap
                else:
                    self.logger.warning("   Not able to read [{0}] ... "
                                + "Skipped!!" .format(self.inputFileList[i]))

            except (KeyboardInterrupt, SystemExit) as e:
                if ccmap is not None:
                    del ccmap
                raise e

class PairCooMatrixHandler:
    """To import ccmap from files similar to paired sparse matrix Coordinate (COO) format

    This format is very similar to COO format with addiotional infromation of chromosome. Therefore,
    maps for all chromosome could be contained in a single file.

    This type of format appeared with following publication:
        * http://dx.doi.org/10.1016/j.cell.2015.10.026 --- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72512

    Following file format can be read as a text file, where first and second column is location on chromosome and third column is the value:
    ::

        chr4     60000   75000   chr4    60000   75000   0.1163470887070292
        chr4     60000   75000   chr4    105000  120000  0.01292745430078102
        chr4     60000   75000   chr4    435000  450000  0.01292745430078102
        chr4     75000   90000   chr4    75000   90000   0.05170981720312409
        chr4     75000   90000   chr4    345000  360000  0.01292745430078102
        chr4     90000   105000  chr4    90000   105000  0.01292745430078102
        .
        .
        .
        .
        .
        .


    Parameters
    ----------
    inputFile : str
        name of a input file
    ccmapOutDir : str
        name of directory where all ccmap file will be stored.
    ccmapSuffix : str
    	Suffix for ccmap file name.
    gcmapOut : str
        Name of output gcmap file.
    workDir : str
        Directory where temporary files will be stored.


    Attributes
    ----------
    inputFile : str
        name of a input file
    ccmapOutDir : str
        name of directory where all ccmap file will be stored.
    ccmapSuffix : str
    	Suffix for ccmap file name.
    gcmapOut : str
        Name of output gcmap file.
    gcmapOutOptions : dict
        Dictionary for gcmap output options.
    workDir : str
        Directory where temporary files will be stored.


    Examples
    --------
        >>> pair_map_handle = PairCooMatrixHandler('GSM1863750_tethered_rep1_contacts.txt', gcmapOut='GSM1863750_tethered_rep1_contacts.gcmap')
        >>> pair_map_handle.setGCMapOptions()
        >>> pair_map_handle.runConversion()


    """

    def __init__(self, inputFile, ccmapOutDir=None, ccmapSuffix=None, gcmapOut=None, workDir=None, logHandler=None):
        self.inputFile = inputFile
        self.ccmapOutDir = ccmapOutDir
        self.ccmapSuffix = ccmapSuffix
        self.gcmapOut = gcmapOut
        self.gcmapOutOptions = None
        self._fin = None
        self.workDir = workDir
        self.logHandler = logHandler

        # Working and output directory
        if workDir is not None:
            self.workDir = workDir
        else:
            self.workDir = config['Dirs']['WorkingDirectory']
        # logger
        self.logHandler = logHandler
        self.logger = logging.getLogger('PairCooMatrixHandler')
        if logHandler is not None:
            self.logger.propagate = False
            self.logger.addHandler(logHandler)
        self.logger.setLevel(logging.INFO)

    def setGCMapOptions(self, compression='lzf', generateCoarse=True, coarsingMethod='sum', replaceCMap=True):
        """ Set options for output gcmap file

        Parameters
        ----------
        compression : str
            Compression method. Presently allowed : ``lzf`` for LZF compression and ``gzip`` for GZIP compression.
        generateCoarse : bool
            Also generates all coarser maps where resolutions will be coarsed by a factor of two, consequetively.
            e.g.: In case of 10 kb input resolution, downsampled maps of ``20kb``, ``40kb``, ``80kb``, ``160kb``, ``320kb`` etc.
            will be generated untill, map size is less than 500.
        coarsingMethod : str
            Method of downsampling. Three accepted methods are ``sum``: sum all values, ``mean``: Average of all values
            and ``max``: Maximum of all values.
        replaceCMap : bool
            Replace entire old ccmap data including resolutions and coarsed data.

        """
        self.gcmapOutOptions = dict()
        self.gcmapOutOptions['compression'] = compression
        self.gcmapOutOptions['generateCoarse'] = generateCoarse
        self.gcmapOutOptions['coarsingMethod'] = coarsingMethod
        self.gcmapOutOptions['replaceCMap'] = replaceCMap

    def _convertProcessedValues(self, i, j, value, chrom):
        self.logger.info('##### CONVERTING FOR {0} #######'.format(chrom))
        ccmap = gen_map_from_locations_value(i, j, value, mapType='intra', workDir=self.workDir, logHandler=self.logHandler)
        ccmap.xlabel = chrom
        ccmap.ylabel = chrom

        if self.ccmapOutDir is not None or self.ccmapSuffix is not None:
            self._save_ccmap(ccmap, chrom, util.binsizeToResolution(ccmap.binsize))

        if self.gcmapOut is not None:
            self._save_gcmap(ccmap)

        del ccmap
        self.logger.info('##### ################# #######\n'.format(chrom))

    def _save_ccmap(self, ccmap, title, resolution):
        if self.ccmapOutDir is None:
            self.ccmapOutDir = os.getcwd()

        if self.ccmapSuffix is None:
            outFileName = title + '_' + resolution + '.ccmap'
        else:
            outFileName = title + '_' + resolution + '_' + self.ccmapSuffix + '.ccmap'

        fullOutPath = os.path.join(self.ccmapOutDir, outFileName)
        cmp.save_ccmap(ccmap, fullOutPath, compress=True)


    def _save_gcmap(self, ccmap):
        if self.gcmapOutOptions is None:
            raise ValueError('No options set for output gcmap file!!. Use PairCooMatrix.setGCMapOptions() to set the options.')

        compression = self.gcmapOutOptions['compression']
        generateCoarse = self.gcmapOutOptions['generateCoarse']
        coarsingMethod = self.gcmapOutOptions['coarsingMethod']
        replaceCMap = self.gcmapOutOptions['replaceCMap']
        gmp.addCCMap2GCMap(ccmap, self.gcmapOut,
                                compression=compression,
                                generateCoarse=generateCoarse,
                                coarsingMethod=coarsingMethod,
                                replaceCMap=replaceCMap,
                                logHandler=self.logHandler)

    def runConversion(self):
        """ Perform conversion and save to ccmap and/or gcmap file.

        Read the input file, process the data, and convert it to ccmap
        or gcmap file. For output gcmap, :meth:`PairCooMatrixHandler.setGCMapOptions`
        should be called to set the neccessary options.
        """
        # Check if gcmap output is given and gcmap options are set
        if self.gcmapOutOptions is None and self.gcmapOut is not None:
            raise ValueError('No options set for output gcmap file!!. Use PairCooMatrix.setGCMapOptions() to set the options.')

        # At the start of file
        if self._fin is None:
            self._fin = open(self.inputFile, 'r')
        else:
            self._fin.seek(0)

        prevChrom = 'dummy'
        i, j, value = [], [], []
        for line in self._fin:

            # Process and break line
            line = line.rstrip().lstrip()
            temp = re.split('\s+', line)

            # Check for intra-chromosomal contact map
            if temp[0] == temp[3] and not temp[0] == prevChrom and value:
                self._convertProcessedValues(i, j, value, prevChrom)
                self.logger.info('##### Started Reading For {0} #######'.format(temp[0]))

                i, j, value = [], [], []

            # assign chromosome name
            prevChrom = temp[0]

            # Append only in case of intra-chromosomal map
            if temp[0] == temp[3]:
                i.append(int(temp[2]))
                j.append(int(temp[5]))
                value.append(float(temp[6]))

        # Convert for last map --- intra-chromosomal
        if temp[0] == temp[3]:
            self._convertProcessedValues(i, j, value, temp[0])
            i, j, value = [], [], []

class HomerInputHandler:
    """To import ccmap from Hi-C maps generated by HOMER

    HOMER package generates the `Hi-C interaction matrices in text file <http://homer.salk.edu/homer/interactions/HiCmatrices.html>`_. This Hi-C interaction file can be imported using this class.
    HOMER format interaction matrix file may contain data for all the chromsomes while this class separatly read and save matrix of each chromosome.

    To Instantiate this class, three scenarios are possible:
        * ``Text in Archive``: Both a list of input files and a compressed file is given. It means look for input files in compressed file.
        * ``Text``: Only a list of input files is given. It means, read data directly from input file.
        * ``Archive``: Only a compressed file is given. It means, read all files present in compressed files.

    Parameters
    ----------
    inputFiles : str or list
        Name of a input file or list of input files. If ``None``, all files from compressed files are used as input files.
    inputCompressedFile : str
        Name of input compressed file. Accepted formats: ``tar.gz``, ``tar.bz2`` and ``zip``.
    workDir : str
        Directory where temporary files will be stored. If it is not provided, this value is taken from configuration file.

    Attributes
    ----------
    inputFileList : list
        List of input files. It could be ``None`` when not provided.
    inputType : str
        Type of input. Three types: ``Text``, ``Text in Archive`` and ``Archive`` are determined from user input.
    inputCompressedFile : str
        Name of input compressed file. It could be ``None`` when not provided.
    compressType : str
        Format of compressed file. It could be etiher ``.tar`` or ``.zip`` or ``None``
    compressHandle : ZipFile or TarFile
        An object to handle compressed file. It could be `ZipFile <https://docs.python.org/3/library/zipfile.html#zipfile-objects>`_ or
        `TarFile <https://docs.python.org/3/library/tarfile.html#tarfile-objects>`_ instance depending on compressed format.
    workDir : str
        Directory where temporary files will be stored.
    fIns : list[output file stream]
        Input file stream for each input files
    chromList : list[str]
        List of chromosome found in input files
    resolution : str
        Resolution of map
    fTmpOutNames : list[str]
        List of temporary output files where data for each chromosomes are extracted separately
    fTmpOut : list[output file stream]
        List of output file streams for respective temporary output files

    """

    def __init__(self, inputFiles=None, inputCompressedFile=None, workDir=None, logHandler=None):

        # Checking both text and compressed input files
        self.inputType = self._checkInputType(inputFiles, inputCompressedFile)  # raise error when both inputs are None

        # Initial processing for archive file
        self.inputCompressedFile = inputCompressedFile
        self.compressType = self._checkCompressedFile()

        # Make a input file list if user has provided one
        self.inputFileList = None
        if inputFiles is not None:
            self.inputFileList = []
            self._addInputTextFiles(inputFiles)

        # If user has not provided input files, read file list from compressed file and return handle object of compressed file
        self.compressHandle = None
        self._extractCompressedFile()

        # Set working directory
        if workDir is not None:
            self.workDir = workDir
        else:
            self.workDir = config['Dirs']['WorkingDirectory']

        self.chromList = None
        self.resolution = None
        self.headers = None
        self.fTmpOut = None
        self.fTmpOutNames = None
        self.fIns = None

        # logger
        # logger
        self.logHandler = logHandler
        self.logger = logging.getLogger('HomerInputHandler')
        if logHandler is not None:
            self.logger.propagate = False
            self.logger.addHandler(logHandler)
        self.logger.setLevel(logging.INFO)

    def __del__(self):
        self._removeTemporaryOutputFiles()
        if self.logHandler is not None:
            self.logger.removeHandler( self.logHandler )

    def _checkInputType(self, inputFiles, compressedFile):
        """To check input type given by user

        There are three possible scenarios:
            * ``Text in Archive``: Both a list of input files and a compressed file is given. It means look for input files in compressed file.
            * ``Text``: Only a list of input files is given. It means, read data directly from input file.
            * ``Archive``: Only a compressed file is given. It means, read all files present in compressed files.

        Parameters
        ----------
        inputFiles : str or list
            name of a input file or list of input files.
        compressedFile : str
            Name of compressed file

        Raises
        ------
        ValueError
            If both input files and compressedFile is ``None``

        """
        inputType = None
        # Three porabable type of inputs 'Text' 'Text in Archive' 'Archive'
        if inputFiles is not None and compressedFile is not None:
            inputType = 'Text in Archive'
        elif inputFiles is not None and compressedFile is None:
            inputType = 'Text'
        elif inputFiles is None and compressedFile is not None:
            inputType = 'Archive'
        else:
            raise ValueError ("Both inputFiles and inputCompressedFile cannot be None.")

        return inputType

    def _checkCompressedFile(self):
        """Determine format of compressed file.

        Returns
        -------
        compressType : str
            One of the two possible keywords (``.tar`` and ``.zip``) or ``None``.

        Raises
        ------
        NotImplementedError
            If format is not ``zip``, ``tar.gz`` or ``tar.bz2``.

        """
        if self.inputCompressedFile is None:  return
        compressType = None
        file_extension = os.path.splitext(self.inputCompressedFile)[1]
        if file_extension == '.zip':
            compressType = 'zip'
        elif file_extension == '.gz' or file_extension == '.tar' or file_extension == '.bz2':
            compressType = 'tar'
        else:
            raise NotImplementedError (" [{0}] compression format is not implemented. Please use '.zip' or '.tar.gz' or '.tar.bz2' compression format." .format(file_extension))

        return compressType

    def _addInputTextFiles(self, inputFiles):
        """ Add Input Files

        It add input files and also convert to list

        Parameters
        ----------
        inputFiles : str or list
            name of a input file or list of input files

        """
        if isinstance(inputFiles, list):
            self.inputFileList = self.inputFileList + inputFiles
        else:
            self.inputFileList.append(inputFiles)

    def _extractCompressedFile(self):
        """Generate input file list from compressed file
        """
        if self.inputCompressedFile is None:    return

        if self.compressType == 'zip':
            self.compressHandle = zipfile.ZipFile(self.inputCompressedFile, 'r')
            if self.inputFileList is None:
                self.inputFileList = self.compressHandle.namelist()
        else:
            self.compressHandle = tarfile.open(self.inputCompressedFile, 'r')
            if self.inputFileList is None:
                self.inputFileList = self.compressHandle.getnames()

    def _openInputFiles(self):
        """Open input files
        Open all input files and store input file stream in :attr:`HOMERInputHandler.fIns` list

        """
        if self.fIns is None:
            self.fIns = []
        else:
            self._closeInputFiles()
            self.fIns = []

        for f in self.inputFileList:
            # First text file
            if self.inputType == 'Text':
                fin = open(f, 'r')
                self.fIns.append(fin)
            else:
                fin = None
                # Zip file
                if self.compressType == 'zip':
                    try:
                        fin = self.compressHandle.open(f)
                    except KeyError:
                        self.logger.warning(' {0} not found in [{1}]... \n\t\t...Skipping this file...'.format(f, self.inputCompressedFile))
                else:
                    # Tar file
                    try:
                        fin = self.compressHandle.extractfile(f)
                    except KeyError:
                        self.logger.warning(' [{0}] not found in [{1}]... \n\t\t...Skipping this file...'.format(f, self.inputCompressedFile))
                self.fIns.append(fin)

    def _closeInputFiles(self):
        """Close all input files
        Close all input file stream present in :attr:`HOMERInputHandler.fIns` list
        """
        if self.fIns is not None:
            for fin in self.fIns:
                if fin is not None:
                    # close function is not available when tar file object is used
                    if self.compressType != 'tar':
                        fin.close()
            del self.fIns
            self.fIns = None

    def _generateTemporaryOutputFiles(self):
        """Generate temporary text files

        Generate temporary text file for each chromosome and matrix is stored as --- i j value --- format.

        """
        if self.fTmpOut is None:
            self.fTmpOut = dict()
        else:
            self._closeTemporaryOutputFiles()
            self.fTmpOut = dict()

        if self.fTmpOutNames is None:
            self.fTmpOutNames = dict()
        else:
            self._removeTemporaryOutputFiles()
            self.fTmpOutNames = dict()

        self._getChromListAndResolution()

        for chrom in self.chromList:

            # Only filename is generated
            fd, fname = tempfile.mkstemp(prefix='gcx_{0}_' .format(chrom), suffix='.tmp', dir=self.workDir)
            os.close(fd)
            if os.path.isfile(fname):
                os.remove(fname)

            fout = open(fname, 'w')
            self.fTmpOut[chrom] = fout
            self.fTmpOutNames[chrom] = fname

    def _closeTemporaryOutputFiles(self):
        """Close all temporary files
        Close all output file stream present in :attr:`HOMERInputHandler.fTmpOut` list
        """
        if self.fTmpOut is not None:
            for key in self.fTmpOut:
                self.fTmpOut[key].close()
            self.fTmpOut = None

    def _openTemporaryOutputFiles(self):
        """Open all temporary files
        Open all output files present in :attr:`HOMERInputHandler.fTmpOutNames` list
        """
        if self.fTmpOut is not None:
            for key in self.fTmpOut:
                self.fTmpOut[key].close()
            self.fTmpOut = None

        self.fTmpOut = dict()
        for key in self.fTmpOutNames:
            self.fTmpOut[key] = open(self.fTmpOutNames[key], 'r')

    def _removeTemporaryOutputFiles(self):
        """Remove all temporary files
        Delete all temporary files present in :attr:`HOMERInputHandler.fTmpOutNames` list from external storage
        """
        self._closeTemporaryOutputFiles()
        if self.fTmpOutNames is not None:
            for key in self.fTmpOutNames:
                if os.path.isfile(self.fTmpOutNames[key]):
                    os.remove(self.fTmpOutNames[key])


    def _getChromListAndResolution(self):
        """ Get chromosome list and resolution

        Read input files headers and determine chromosome names and resolution. After call this function, chromosome list can be found in :attr:`HOMERInputHandler.chromList`
        and Resolution can be found in :attr:`HOMERInputHandler.resolution`

        """

        self._openInputFiles()

        # Removing previous header list
        if self.headers is None:
            self.headers = []
        else:
            del self.headers
            self.headers = []

        chroms = []
        self.logger.info(" Getting chromosome list and resolution from Input Files ...")
        for fin in self.fIns:
            if fin is None: continue
            binsize = None

            # Only read first line
            if self.compressType is None:
                line = fin.readline()
            else:
                line = fin.readline().decode()

            line = line.lstrip().rstrip()
            match = re.search('Regions', line)
            header = re.split('\s+', line[match.start():])

            count = 0
            prev_j = None
            step = []
            for m in range(len(header[1:])):
                chrm, j = re.split('-', header[m+1])
                chroms.append( self._getChromName(chrm) )

                if count < 100 and prev_j is not None:
                    ds = abs( int(j) - int(prev_j) )
                    if ds > 0:
                        step.append(ds)
                        count += 1
                prev_j = j

            if binsize is None:
                self.resolution = util.binsizeToResolution( int( np.amin(step) ) )
            else:
                if binsize != int( np.amin(step) ):
                    raise AssertionError ('Data Resolution does not match in input files.')

            # Appending header of each file
            self.headers.append(header)

        self.chromList = util.sorted_nicely(list(set(chroms)))

        self.logger.info(" Resolution: {0}" .format(self.resolution))
        outline = ' '
        for chrm in self.chromList:
            outline = outline + "\n                              " + chrm
        self.logger.info(" Following chromsomes found in input files: {0}" .format(outline))

    def _getChromName(self, chrm):
        """Generate chromosome name
        Sometimes, input file only contains chromosome names. This function change number to name by prefixing it with ``chr``

        """
        if re.search('(chr)|(Chr)', chrm) == None:
            return 'chr' + chrm
        else:
            return chrm

    def _ProcessInputs(self):
        """Read input files and save map to temporary text file

        This function reads input files and save data temporarily in output
        text files. Name of these temporary files are listed
        in :attr:`HOMERInputHandler.fTmpOutNames` .

        """

        self._generateTemporaryOutputFiles()

        try:
            for n in range(len(self.fIns)):

                # When a file is not found in compressed file, file object could be none
                if self.fIns[n] is None:    continue

                self.logger.info(" Reading [{0}] file ..." .format(self.inputFileList[n]))

                while(1):
                    if self.compressType is None:
                        line = self.fIns[n].readline()
                    else:
                        line = self.fIns[n].readline().decode()

                    if not line:    break

                    line = line.lstrip().rstrip()
                    temp = re.split('\s+', line)

                    # Header is taken from available header list. Also, first line is already read from each input file
                    chrm1, i = re.split('-', temp[1])
                    #print(len(temp[2:]), chrm1, i, len(header[1:]))
                    for m in range(len(self.headers[n][1:])):
                        chrm2, j = re.split('-', self.headers[n][m+1])
                        if chrm1 == chrm2:
                            val = float(temp[m+2])
                            if val != 0:
                                self.fTmpOut[self._getChromName(chrm1)].write('{0}\t{1}\t{2}\n' .format(i, j, val))

                self.logger.info("           ... Finished reading [{0}] file ..." .format(self.inputFileList[n]))

        except (KeyboardInterrupt, SystemExit) as e:
            self._removeTemporaryOutputFiles()
            raise e

        self._closeInputFiles()
        self._closeTemporaryOutputFiles()

    def save_ccmaps(self, outdir, suffix=None):
        """Import and save ccmap file

        Read input files, save data temporarily in text file for each chromosome and import these data to native ccmap format using
        :class:`CooMatrixHandler` class.

        .. note::
            * Output file names will be automatically generated as ``<chromosome>_<resolution>.ccmap`` format. e.g. ``chr12_10kb.ccmap``.
            * A suffix can be added to all files as ``<chromosome>_<resolution>_<suffix>.ccmap`` format. e.g. if ``suffix='_RawObserved'``, file name is ``chr12_10kb_RawObserved.ccmap``.

        Parameters
        ----------
        outdir : str
            Path to directory where all ccmaps have to be saved
        suffix : str
            Any suffix to file name
        """
        self._ProcessInputs()

        outputFiles = []
        inputFiles = []
        try:
            for i in range(len(self.chromList)):

                if suffix is None:
                    outFileName = self.chromList[i] + '_' + self.resolution + '.ccmap'
                else:
                    outFileName = self.chromList[i] + '_' + self.resolution + '_' + suffix + '.ccmap'

                outFile = os.path.join(outdir, outFileName)
                outputFiles.append(outFile)

                inputFiles.append(self.fTmpOutNames[self.chromList[i]])

            reader = CooMatrixHandler(inputFiles, resolution=self.resolution)
            reader.save_ccmaps(outputFiles,
                                xlabels=self.chromList,
                                ylabels=self.chromList)
            self.logger.info(" Saved {0} files." .format(outputFiles))

        except (KeyboardInterrupt, SystemExit) as e:
            self._removeTemporaryOutputFiles()
            raise e

    def save_gcmap(self, outputFile, coarsingMethod='sum', compression='lzf'):
        """To Save all Hi-C maps as a gcmap file

        This function reads input files one by one and save it as a ``.gcmap``
        file.

        Parameters
        ----------
        outputFile : str
            Name of a output gcmap file.
        coarsingMethod : str
            Method of downsampling. Three accepted methods are ``sum``: sum all
            values, ``mean``: Average of all values and ``max``: Maximum of
            all values.
        compression : str
            Compression method. Presently allowed : ``lzf`` for LZF compression
            and ``gzip`` for GZIP compression.

        """

        self._ProcessInputs()

        inputFiles = []
        try:
            for i in range(len(self.chromList)):
                inputFiles.append(self.fTmpOutNames[self.chromList[i]])

            reader = CooMatrixHandler(inputFiles, resolution=self.resolution)
            reader.save_gcmap(outputFile,
                                xlabels=self.chromList,
                                ylabels=self.chromList,
                                coarsingMethod=coarsingMethod,
                                compression=compression)

        except (KeyboardInterrupt, SystemExit) as e:
            self._removeTemporaryOutputFiles()
            raise e

class BinsNContactFilesHandler:
    """To import Hi-C map from bin and contact file in list format

    These types of files are appeared in following GEO data:
        * http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61471
        * http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34453

    Parameters
    ----------
    binFile : str
        Bin file
    contactFile : str
        Contact file in list format
    workDir : str
        Directory where temporary files will be stored. If it is not provided, this value is taken from configuration file.

    Attributes
    ----------
    binFile : str
        Bin file
    contactFile : str
        Contact file in list format
    ChromSize : dict
        Dictionary of Chromosome size
    ChromBinsInfo : dict
        Dictionary of min (start) and max (end) bin number for each Chromosome
    npyBinFileList : dict
        Dictionary containing tuple (memmap stream, Temporary numpy array file name)
    binsize : int
        Bin size of Hi-C map
    ccmaps : dict
        Dictionary of CCMAP instances for each Chromosome


    """
    def __init__(self, binFile, contactFile, workDir=None, logHandler=None):
        self.binFile = binFile
        self.contactFile = contactFile

        self.ChromSize = None           # Store Chromosome size
        self.ChromBinsInfo = None       # Store max and min Bins for each Chromosome
        self.npyBinFileList = None      # Temporary numpy array dict containing tuple (npy, filename)

        self.binsize = None

        self.ccmaps = None            # ccmap list for chromosomes

        if workDir is not None:
            self.workDir = workDir
        else:
            self.workDir = config['Dirs']['WorkingDirectory']

        # logger
        self.logHandler = logHandler
        self.logger = logging.getLogger('BinsNContactFilesHandler')
        if logHandler is not None:
            self.logger.propagate = False
            self.logger.addHandler(logHandler)
        self.logger.setLevel(logging.INFO)

    def __del__(self):
        self._removeTempFiles()
        if self.logHandler is not None:
            self.logger.removeHandler( self.logHandler )

    def _removeTempFiles(self):
        """Remove temporary numpy array files
        """
        if self.npyBinFileList is not None:
            for key in self.npyBinFileList:
                self.npyBinFileList[key][0].close()

                self.logger.info(' Removing temporary numpy array file [{0}] for {1} ...' .format(self.npyBinFileList[key][1], key))
                if os.path.isfile(self.npyBinFileList[key][1]):
                    os.remove(self.npyBinFileList[key][1])

    def _getChromName(self, chrm):
        """Generate chromosome name
        Sometimes, input file only contains chromosome names. This function change number to name by prefixing it with ``chr``

        """
        if re.search('(chr)|(Chr)', chrm) == None:
            return 'chr' + chrm
        else:
            return chrm

    def _printChromosomeInfo(self):
        """Print information about chromosomes
        """

        outline = '  Chromosome Size: \n'
        for key in self.ChromSize:
            outline = outline + '{0:>50} : {1}\n' .format(key, self.ChromSize[key])
        self.logger.info(outline)

        if self.ChromBinsInfo is not None:
            outline = '  Chromosome Bins info:\n'
            for key in self.ChromBinsInfo:
                outline = outline +  '{0:>55}: {1}\n' .format(key, self.ChromBinsInfo[key])
            self.logger.info(outline)

        if self.npyBinFileList is not None:
            outline = '   Temporary array size info:\n'
            for key in self.npyBinFileList:
                outline = outline + '{0:>70}: {1}\n' .format( key, self.npyBinFileList[key][0].shape[0] )
            self.logger.info(outline)

    def _printHiCmapInfo(self):
        """Print information about Hi-C maps
        """
        if self.ccmaps is not None:
            outline = '   Hi-C Maps Summary:\n'
            outline = outline + '{0:>60}\tSize\t\tMax. \tMin. \n'.format('Chromosome')
            for key in self.ccmaps:
                self.ccmaps[key].make_readable()
                outline = outline + '{0:>60}\t{1}\t{2}\t{3}\n' .format(key, self.ccmaps[key].shape, self.ccmaps[key].maxvalue, self.ccmaps[key].minvalue)
            self.logger.info(outline)

    def _readBinFile(self):
        """Read bin file

        Reads bin file and generates ``ChromSize`` and ``ChromBinsInfo`` dictionary
        """

        fin = open(self.binFile, 'r')

        ChromTitle = 'dummy'
        PrevChromTitle = 'dummy'

        for line in fin:
            line = line.lstrip().rstrip()

            # Do not read a line that start with word
            if re.match('\d+', line) == None:
                continue

            temp = re.split('\s+', line)
            s0 = int(temp[2])
            s1 = int(temp[3])

            if self.binsize is not None:
                if self.binsize != s1-s0:
                    raise ValueError('Length of bins are not constant in files. Exiting.')
            else:
                self.binsize = s1-s0

            PrevChromTitle = ChromTitle
            ChromTitle = temp[1]

            if ChromTitle !=  PrevChromTitle:
                if self.ChromSize is None:
                    self.ChromSize = dict()

                if self.ChromBinsInfo is None:
                    self.ChromBinsInfo = dict()

                if ChromTitle not in self.ChromBinsInfo:
                    self.ChromBinsInfo[ChromTitle] = { 'min': int(temp[0]), 'max':  None }

                self.ChromSize[ChromTitle] = s1
                self.ChromBinsInfo[ChromTitle]['max'] = int(temp[0])

            else:
                self.ChromSize[ChromTitle] = s1
                self.ChromBinsInfo[ChromTitle]['max'] = int(temp[0])

        # Print information about chromosome
        self._printChromosomeInfo()

        # Generating temporary numpy array files
        for key in self.ChromSize:
            if self.npyBinFileList is None:
                self.npyBinFileList = dict()

            (fout, fname) = tempfile.mkstemp(suffix='.tmp', prefix='gcx_nparray_'+key+'_', dir=self.workDir, text=False)
            os.close(fout)
            self.logger.info(' Generating temporary numpy array file [{0}] for {1} ...' .format(fname, key))

            size = int(self.ChromSize[key]/self.binsize)+1
            np.save(fname, np.zeros(size, dtype=np.int))
            npy = np.load(fname, mmap_mode='r+')
            self.npyBinFileList[key] = (npy, fname)

            self.logger.info(' Finished. \n' .format(fname, key))


        # Re-read file to store bin index
        fin.seek(0)
        for line in fin:
            line = line.lstrip().rstrip()

            # Do not read a lin that start with word
            if re.match('\d+', line) == None:
                continue

            temp = re.split('\s+', line)
            ChromTitle = temp[1]

            idx = int(int(temp[3])/self.binsize)
            self.npyBinFileList[ChromTitle][0][idx] = int(temp[0])

        fin.close()

    def _readContactFile(self):
        """Read contact file and generate CCMAP instances
        """
        fin = open(self.contactFile, 'r')

        ChromTitle = 'dummy'
        PrevChromTitle = 'dummy'

        self.logger.info(' Reading contact file ...\n')
        i, j, value = [], [], []
        for line in fin:
            line = line.lstrip().rstrip()

            # Do not read a lin that start with word
            if re.match('\d+', line) == None:
                continue

            temp = re.split('\s+', line)
            iBin = int(temp[0])
            jBin = int(temp[1])
            c = float(temp[3])

            for key in self.ChromBinsInfo:
                if (min(iBin, jBin) >= self.ChromBinsInfo[key]['min']) and (max(iBin, jBin) <=  self.ChromBinsInfo[key]['max']):
                    idx = np.nonzero( self.npyBinFileList[key][0] == iBin )
                    jdx = np.nonzero( self.npyBinFileList[key][0] == jBin )
                    i.append(idx[0][0]*self.binsize)
                    j.append(jdx[0][0]*self.binsize)
                    value.append(c)

                    ChromTitle = key

            if ChromTitle != PrevChromTitle and PrevChromTitle != 'dummy':
                self.logger.info(' \tGenerating Hi-C Map for [{0}] ... \n' .format(PrevChromTitle))
                ccmap = gen_map_from_locations_value(i, j, value, resolution=util.binsizeToResolution(self.binsize), workDir=self.workDir, mapType='intra')
                ccmap.xlabel = self._getChromName(PrevChromTitle)
                ccmap.ylabel = self._getChromName(PrevChromTitle)

                if self.ccmaps is None:
                    self.ccmaps = dict()

                self.ccmaps[PrevChromTitle] = ccmap
                i, j, value = [], [], []
                self.logger.info(' Finished\n' .format(PrevChromTitle))


            PrevChromTitle = ChromTitle


        # Last One
        self.logger.info(' \tGenerating Hi-C Map for [{0}] ... ' .format(ChromTitle))
        ccmap = gen_map_from_locations_value(i, j, value, resolution=util.binsizeToResolution(self.binsize),  workDir=self.workDir, mapType='intra')
        ccmap.xlabel = self._getChromName(PrevChromTitle)
        ccmap.ylabel = self._getChromName(PrevChromTitle)

        self.ccmaps[ChromTitle] = ccmap
        i, j, value = [], [], []
        self.logger.info(' Finished\n' .format(ChromTitle))

        self.logger.info(' Finished reading contact file.\n')

        # Print Hi-C Map information summary
        self._printHiCmapInfo()

    def _genHiCMaps(self):
        """Read both bin and contact file and generated CCMAP
        """
        try:
            self._readBinFile()
            self._readContactFile()
        except (SystemExit, KeyboardInterrupt) as e:
            self._removeTempFiles()
            if self.ccmaps is not None:
                for key in self.ccmaps:
                    del self.ccmaps[key]
            raise e

    def save_ccmaps(self, outdir, suffix=None):
        """Import and save ccmap file

        Read input files, save data temporarily for each chromosome and import these data to native ccmap format.

        ..note::
            * Output file names will be automatically generated as ``<chromosome>_<resolution>.ccmap`` format. e.g. ``chr12_10kb.ccmap``.
            * A suffix can be added to all files as ``<chromosome>_<resolution><suffix>.ccmap`` format. e.g. if ``suffix='_RawObserved'``, file name is ``chr12_10kb_RawObserved.ccmap``.

        Parameters
        ----------
        outdir : str
            Path to directory where all ccmaps have to be saved
        suffix : str
            Any suffix to file name
        """
        if self.ccmaps is None:
            self._genHiCMaps()

        resolution = util.binsizeToResolution(self.binsize)
        for key in self.ccmaps:
            if suffix is None:
                outFileName = self._getChromName(key) + '_' + resolution + '.ccmap'
            else:
                outFileName = self._getChromName(key) + '_' + resolution + suffix + '.ccmap'

            fullOutPath = os.path.join(outdir, outFileName)

            cmp.save_ccmap(self.ccmaps[key], fullOutPath, compress=True)

    def save_gcmap(self, outputFile, coarsingMethod='sum', compression='lzf'):
        """To Save all Hi-C maps as a gcmap file

        This function reads input files one by one and save it as a ``.gcmap``
        file.

        Parameters
        ----------
        outputFile : str
            Name of a output gcmap file.
        coarsingMethod : str
            Method of downsampling. Three accepted methods are ``sum``: sum all
            values, ``mean``: Average of all values and ``max``: Maximum of
            all values.
        compression : str
            Compression method. Presently allowed : ``lzf`` for LZF compression
            and ``gzip`` for GZIP compression.

        """

        if self.ccmaps is None:
            self._genHiCMaps()

        resolution = util.binsizeToResolution(self.binsize)
        for key in self.ccmaps:
            gmp.addCCMap2GCMap(self.ccmaps[key], outputFile,
                                compression=compression,
                                coarsingMethod=coarsingMethod)


def gen_map_from_locations_value(i, j, value, resolution=None, mapType='intra', workDir=None, logHandler=None):
    """To generate CCMAP object from three lists -- i, j, value

    Parameters
    ----------
    i : list[int]
        List of first location from each pair
    j : list[int]
        List of second location from each pair
    resolution : str
        Resolution of Hi-C map
    mapType : str
        Hi-C map type: ``intra`` or ``inter`` chromosomal map
    value : list[float]
        List of values for respective location
    workDir : str
        Directory where temporary files will be stored. If it is not provided, this value is taken from configuration file.

    Returns
    -------
    ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP`
    	A CCMAP object

    """

    # logger
    logHandler = logHandler
    logger = logging.getLogger('genMapFromLists')
    if logHandler is not None:
        logger.propagate = False
        logger.addHandler(logHandler)
    logger.setLevel(logging.INFO)

    # set working directory
    if workDir is None:
        workDir = config['Dirs']['WorkingDirectory']

    ccMapObj = cmp.CCMAP()

    maxLi = np.amax(i)
    maxLj = np.amax(j)

    maxL = np.amax([i, j])

    if resolution is not None:
        ccMapObj.binsize = util.resolutionToBinsize(resolution)
    else:
        step = []
        max_count = len(i)
        for n in range(1, 100):
            ds = abs(i[n] - i[n-1])
            if ds > 0:
                step.append(ds)
            ds = abs(j[n] - j[n-1])
            if ds > 0:
                step.append(ds)

            if n > max_count:
                break

        ccMapObj.binsize = int( np.amin(step) )

    logger.info(" Total number of data in input file: {0}" .format(len(i)))
    if mapType == 'intra':
    	minL = np.amin([i, j])
    	logger.info("Minimum base-pair: {0} and Maximum base-pair: {1} are present in input data" .format(minL, maxL))

    	xticks = np.arange(0, maxL + ccMapObj.binsize, ccMapObj.binsize)
    	yticks = np.arange(0, maxL + ccMapObj.binsize, ccMapObj.binsize)
    	ccMapObj.xticks = [0, maxL + ccMapObj.binsize]
    	ccMapObj.yticks = [0, maxL + ccMapObj.binsize]
    else:
    	xticks = np.arange(0, maxLi + ccMapObj.binsize, ccMapObj.binsize)
    	yticks = np.arange(0, maxLj + ccMapObj.binsize, ccMapObj.binsize)
    	ccMapObj.xticks = [0, maxLi + ccMapObj.binsize]
    	ccMapObj.yticks = [0, maxLj + ccMapObj.binsize]

    xdlen = len(xticks)
    ydlen = len(yticks)
    ccMapObj.shape = (xdlen, ydlen)

    # Creating a file name for binary numpy array on disk
    ccMapObj.gen_matrix_file(workDir=workDir)

    # Creating a file for binary numpy array on disk
    ccMapObj.make_writable()

    logger.info("Shape of overall map: ({0}, {0})\n" .format(xdlen, ydlen))

    # print(len(i), len(value), ccMapObj.binsize, ccMapObj.shape)
    try:
        for n in range(len(i)):
            ni = int( i[n]/ccMapObj.binsize )
            nj = int( j[n]/ccMapObj.binsize )

            # In case if two observation is present for same pair of locations. take average of them
            if ccMapObj.matrix[ni][nj] != 0:
                ccMapObj.matrix[ni][nj] = (ccMapObj.matrix[ni][nj] + value[n])/2
            else:
                ccMapObj.matrix[ni][nj] = value[n]

            if mapType == 'intra':
                ccMapObj.matrix[nj][ni] = ccMapObj.matrix[ni][nj]

            if (ccMapObj.maxvalue < value[n]):
            	ccMapObj.maxvalue = value[n]
    except (KeyboardInterrupt, SystemExit) as e:
    	del ccMapObj
    	raise e

    # Get minimum value after zero
    ma = np.ma.masked_equal(ccMapObj.matrix, 0.0, copy=False)
    ccMapObj.minvalue = ma.min()
    del ma

    ccMapObj.make_unreadable()

    if logHandler is not None:
        logger.removeHandler( logHandler )

    return ccMapObj
