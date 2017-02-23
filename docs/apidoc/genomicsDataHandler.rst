
genomicsDataHandler module
==========================

This module is developed to visualize and analyze Genomics data with respect to Hi-C maps. This module contains method to convert bigWig and Wig file to hdf5 file.
The hdf5 file gives us flexibility to access the data for given range of location of a specific chromosome at particular resolution.


.. currentmodule:: gcMapExplorer.lib.genomicsDataHandler

.. autosummary::
		HDF5Handler
		HDF5Handler.getChromList
		HDF5Handler.getResolutionList
		HDF5Handler.getDataNameList
		HDF5Handler.buildDataTree
		BigWigHandler
		BigWigHandler.getBigWigInfo
		BigWigHandler.bigWigtoWig
		BigWigHandler.saveAsH5
		WigHandler
		WigHandler.parseWig
		WigHandler.setChromosome
		WigHandler.saveAsH5
		WigHandler.getRawWigDataAsDictionary
		BEDHandler
		BEDHandler.parseBed
		BEDHandler.setChromosome
		BEDHandler.saveAsH5
		TextFileHandler
		TextFileHandler.readData
		TempNumpyArrayFiles
		TempNumpyArrayFiles.updateArraysByBigWig
		TempNumpyArrayFiles.updateArraysByChromSize
		TempNumpyArrayFiles.addChromSizeInfo
		TempNumpyArrayFiles.genrateAllTempNumpyFiles
		TempNumpyArrayFiles.generateTempNumpyFile
		TempNumpyArrayFiles.fillAllArraysWithZeros


HDF5Handler class
-----------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.HDF5Handler
		:members:
		:private-members:

BigWigHandler class
-------------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.BigWigHandler
		:members:
		:private-members:

WigHandler class
----------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.WigHandler
		:members:
		:private-members:


BEDHandler class
----------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.BEDHandler
		:members:
		:private-members:

TempNumpyArrayFiles class
-------------------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.TempNumpyArrayFiles
  :members: _getBigWigInfo, _generateTempNumpyFile, updateArraysByBigWig, updateArraysByChromSize, addChromSizeInfo, genrateAllTempNumpyFiles, generateTempNumpyFile, fillAllArraysWithZeros


TextFileHandler class
---------------------

.. autoclass:: gcMapExplorer.lib.genomicsDataHandler.TextFileHandler
  :members: _getBinSize, readData
