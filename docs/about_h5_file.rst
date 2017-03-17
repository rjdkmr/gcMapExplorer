.. |hdf5| raw:: html

    <a href="https://www.hdfgroup.org/HDF5/" target="_blank"> HDF5 </a>

.. |lzf| raw:: html
    
    <a href="http://www.h5py.org/lzf/" target="_blank"> LZF </a>


Genomic track HDF5 (.h5) file
=============================
To enable rapid visualization of genomic track datasets along with contact map,
we used |hdf5| format file to store these datasets at various resolutions. HDF5
is a binary indexed file and therefore, entire datasets or a portion of dataset
can be accessed rapidly.

HDF5 library is available for C, C++, R, Java and Python programming language,
and therefore these files can be directly read through these languages.


Downsampling or Coarsening of datasets
--------------------------------------
Genomic contact map can be of different resolutions, and therefore, resolution
of corresponding genomic track datasets should match during the visualization/analysis.
Therefore, genomic dataset need to be downsampled or coarsened. However,
there are several possible methods for downsampling that can be suitable for
different purposes. Therefore, we have implmented six different methods as
follows,

* Arithmatic mean
* Geometric mean
* Harmonic mean
* Median
* Maximum
* Minimum

When a file is opened in the ``browser``, user receives a prompt for selection of
downsampling method, and subsequently, the selected data is loaded into the
``browser``.

Structure of genomic track (.h5) file
-------------------------------------

Format: ``/<Chromosome>/<Resolution>/<1D Numpy Array>``

::

  HDF5 ──────────────────────────> title
    ├──────── chr1
    │           ├───── 1kb
    │           │        ├──────── amean  ( Arithmatic mean) (type: 1D Array)
    │           │        ├──────── median ( Median value   ) (type: 1D Array)
    │           │        ├──────── hmean  ( Harmonic mean  ) (type: 1D Array)
    │           │        ├──────── gmean  ( Geometric mean ) (type: 1D Array)
    │           │        ├──────── min    ( Minimum value  ) (type: 1D Array)
    │           │        └──────── max    ( Maximum value  ) (type: 1D Array)
    │           │
    │           ├────  5kb
    │           │        ├──────── amean  ( Arithmatic mean) (type: 1D Array)
    │           │        ├──────── median ( Median value   ) (type: 1D Array)
    │           │        ├──────── hmean  ( Harmonic mean  ) (type: 1D Array)
    │           │        ├──────── gmean  ( Geometric mean ) (type: 1D Array)
    │           │        ├──────── min    ( Minimum value  ) (type: 1D Array)
    │           │        └──────── max    ( Maximum value  ) (type: 1D Array)
    │           │
    │           └────  ...
    │
    ├──────── chr2
    │           ├───── 1kb
    │           │        ├──────── amean  ( Arithmatic mean) (type: 1D Array)
    │           │        ├──────── median ( Median value   ) (type: 1D Array)
    │           │        ├──────── hmean  ( Harmonic mean  ) (type: 1D Array)
    │           │        ├──────── gmean  ( Geometric mean ) (type: 1D Array)
    │           │        ├──────── min    ( Minimum value  ) (type: 1D Array)
    │           │        └──────── max    ( Maximum value  ) (type: 1D Array)
    │           └────  ..
    :
    :
    :
    └───── ...



Compression
-----------
In h5 file, dataset is stored as an 1D array. Presently, two compression methods 
are allowed in the h5 file:

* |lzf|
* GZIP

By default, |lzf| is used to compress arrays. This method is very fast, and allow 
the reading.

.. Warning::
    |lzf| method is only avaiable through **Python h5py** module, and
    therefore, this file cannot be read by another programming language through
    standard library.
    
    For portablity, use GZIP compression method, which is available in standard 
    HDF5 library.


Convert bigWig/wig/bed to genomic track h5 file
-----------------------------------------------
To convert bigWig/wig/bed files to genomic track files a GUI application and
several commands are available.

.. toctree::
    h5Converter : A GUI application to convert bigWig/wig/bed <commands/h5Converter>
    bigwig2h5 : convert bigWig to h5 <commands/bigWig2h5>
    wig2h5 : convert wig to h5 <commands/wig2h5>
    bed2h5 : convert bed to h5 <commands/bed2h5>


Convert using ``gcMapExplorer`` Python modules:
    * bigWig file: :class:`gcMapExplorer.lib.genomicsDataHandler.BigWigHandler`
    * wig file: :class:`gcMapExplorer.lib.genomicsDataHandler.WigHandler`
    * bed file: :class:`gcMapExplorer.lib.genomicsDataHandler.BEDHandler`


