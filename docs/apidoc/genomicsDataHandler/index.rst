genomicsDataHandler module
==========================

This module is developed to visualize and analyze Genomics data with respect
to Hi-C maps. This module contains method to convert bigWig and Wig file to
hdf5 file.

The hdf5 file gives us flexibility to access the data for given range of
location of a specific chromosome at particular resolution.


List of class
-------------

.. toctree::
   :maxdepth: 2

   HDF5Handler <hdf5handler>
   BigWigHandler <bigwighandler>
   WigHandler <wighandler>
   BEDHandler <bedhandler>
   EncodeDatasetsConverter <EncodeDatasetsConverter>
   TextFileHandler <textHandler>
   TempNumpyArrayFiles <tempnumpyarrayfiles>
