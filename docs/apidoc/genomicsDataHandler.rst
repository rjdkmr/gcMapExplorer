
genomicsDataHandler module
==========================

This module is developed to visualize and analyze Genomics data with respect to Hi-C maps. This module contains method to convert bigWig and Wig file to hdf5 file.
The hdf5 file gives us flexibility to access the data for given range of location of a specific chromosome at particular resolution.


.. currentmodule:: gcMapExplorer.lib

.. autosummary::
		genomicsDataHandler.HDF5Handler
		genomicsDataHandler.HDF5Handler.getChromList
		genomicsDataHandler.HDF5Handler.getResolutionList
		genomicsDataHandler.HDF5Handler.getDataNameList
		genomicsDataHandler.HDF5Handler.buildDataTree
		genomicsDataHandler.BigWigHandler
		genomicsDataHandler.BigWigHandler.getBigWigInfo
		genomicsDataHandler.BigWigHandler.bigWigtoWig
		genomicsDataHandler.BigWigHandler.saveToH5
		genomicsDataHandler.WigHandler
		genomicsDataHandler.WigHandler.parseWig
		genomicsDataHandler.WigHandler.convertWigToH5
		genomicsDataHandler.WigHandler.getRawWigDataAsDictionary
		genomicsDataHandler.TextFileHandler
		genomicsDataHandler.TextFileHandler.readData
		genomicsDataHandler.TempNumpyArrayFiles
		genomicsDataHandler.TempNumpyArrayFiles.updateArraysByBigWig
		genomicsDataHandler.TempNumpyArrayFiles.updateArraysByChromSize
		genomicsDataHandler.TempNumpyArrayFiles.addChromSizeInfo
		genomicsDataHandler.TempNumpyArrayFiles.genrateAllTempNumpyFiles
		genomicsDataHandler.TempNumpyArrayFiles.generateTempNumpyFile
		genomicsDataHandler.TempNumpyArrayFiles.fillAllArraysWithZeros


HDF5Handler class
-----------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.HDF5Handler
  :members: getChromList, getResolutionList, getDataNameList, buildDataTree, _addDataByArray


BigWigHandler class
-------------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.BigWigHandler
  :members: _getBigWigInfo, _bigWigtoWig, getBigWigInfo, bigWigtoWig, saveToH5



TextFileHandler class
---------------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.TextFileHandler
  :members: _getBinSize, readData


WigHandler class
----------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.WigHandler
  :members: _PerformDataCoarsing, _getChromTitle_parseWig, _getChromTitleBedgraph_parseWig, _getStartStepFixedStep_parseWig, _getSpan_parseWig, _FillDataInNumpyArrayFile, _StoreInHdf5File, parseWig, _parseWig, convertWigToH5, getRawWigDataAsDictionary


TempNumpyArrayFiles class
-------------------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.TempNumpyArrayFiles
  :members: _getBigWigInfo, _generateTempNumpyFile, updateArraysByBigWig, updateArraysByChromSize, addChromSizeInfo, genrateAllTempNumpyFiles, generateTempNumpyFile, fillAllArraysWithZeros
