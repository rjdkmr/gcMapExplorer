gcMapExplorer.genomicsDataHandler
=================================

This module is developed to visualize and analyze Genomics data with respect to Hi-C maps. This module contains method to convert bigWig and Wig file to hdf5 file.
The hdf5 file gives us flexibility to access the data for given range of location of a specific chromosome at particular resolution.


.. currentmodule:: gcMapExplorer.lib.genomicsDataHandler


HDF5Handler
-----------

.. autoclass:: HDF5Handler
  :members: getChromList, getResolutionList, getDataNameList, buildDataTree, _addDataByArray


BigWigHandler
-------------

.. autoclass:: BigWigHandler
  :members: _getBigWigInfo, _bigWigtoWig, getBigWigInfo, bigWigtoWig, saveToH5



TextFileHandler
---------------

.. autoclass:: TextFileHandler
  :members: _getBinSize, readData


WigHandler
----------

.. autoclass:: WigHandler
  :members: _PerformDataCoarsing, _getChromTitle_parseWig, _getChromTitleBedgraph_parseWig, _getStartStepFixedStep_parseWig, _getSpan_parseWig, _FillDataInNumpyArrayFile, _StoreInHdf5File, parseWig, _parseWig, convertWigToH5, getRawWigDataAsDictionary


TempNumpyArrayFiles
-------------------

.. autoclass:: TempNumpyArrayFiles
  :members: _getBigWigInfo, _generateTempNumpyFile, updateArraysByBigWig, updateArraysByChromSize, addChromSizeInfo, genrateAllTempNumpyFiles, generateTempNumpyFile, fillAllArraysWithZeros
