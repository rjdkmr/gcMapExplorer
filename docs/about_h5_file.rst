Genomic track HDF5 (.h5) file
=============================
To enable rapid visualization of genomic track datasets along with contact map,
we used HDF5 format file to store these datasets at various resolutions. HDF5
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


Convert bigWig/wig/bed to genomic track h5 file
-----------------------------------------------
To convert bigWig/wig/bed files to genomic track files a GUI application and
several commands are available.

* `h5Converter`_: A GUI application to convert all the command.
* `bigwig2h5`_: A command to convert bigWig file
* `wig2h5`_: A command to convert wig file
* `bed2h5`_: A command to convert bed file



h5Converter
~~~~~~~~~~~
This is a GUI application, which can be used to convert bigWig/wig/bed file into
``gcMapExplorer browser`` compatible h5 file. To open this interface, use following
command:

::

    gcMapExplorer h5Converter


.. figure:: images/h5Converter.png
      :scale: 90%
      :alt: Screen snapshot of genomic track dataset converter

      h5Converter Interface

This interface contains in-built
help (**Click on** ``??`` **button**) to understand the functionality of interfaces.


bigwig2h5
~~~~~~~~~

.. include:: commands/bigWig2h5.rst

wig2h5
~~~~~~

.. include:: commands/wig2h5.rst

bed2h5
~~~~~~

.. include:: commands/bed2h5.rst
