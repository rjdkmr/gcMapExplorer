
.. _browser: http://gcmapexplorer.readthedocs.io/en/latest/mapBrowser.html
.. _genome contact map: http://gcmapexplorer.readthedocs.io/en/latest/about_gcmap_file.html
.. _genomic track datasets: http://gcmapexplorer.readthedocs.io/en/latest/about_h5_file.html
.. _Normalization of contact maps: http://gcmapexplorer.readthedocs.io/en/latest/cmapNormalization.html
.. _convert bigWig/wig/bed file to genomic track dataset h5 file: http://gcmapexplorer.readthedocs.io/en/latest/about_h5_file.html#convert-bigwig-wig-bed-to-genomic-track-h5-file

.. image:: https://travis-ci.org/rjdkmr/gcMapExplorer.svg?branch=master
    :target: https://travis-ci.org/rjdkmr/gcMapExplorer


Genome Contact Map Explorer - gcMapExplorer
===========================================

It is a platform to visualize and analyze the contact maps that are generated from Hi-C experiments. This package is developed by considering the huge size of contact maps at very fine resolution. It contains

* Graphical user interface - Several windows like applications to perform tasks (See below tables).
* Command Line Interface - Several commands to perform tasks (See below tables).
* `Application Programming Interface <http://gcmapexplorer.readthedocs.io/en/latest/apidoc/summary.html>`_
  - It can be used to perform analysis by any mathematical operations through programming.


**For more details, visit:** `gcMapExplorer Homepage <http://gcmapexplorer.readthedocs.io/>`_

**For Discussion and Questions, visit** `this forum <https://groups.google.com/forum/#!forum/gcmapexplorer>`_

Features
--------

* Support for **huge contact maps** - Use of Disk instead of RAM - Matrices/arrays are stored in Disks -
  mathematical operations by directly reading/writing from/to Disks, **without loading them into RAM**
* A browser_ with rich interfaces
  for **Comparative** and **Interactive** visualization of **two dimensional contact maps** along
  with **genomic datasets** such as produced by DNase-seq, ChIP-seq, RNA-seq etc.
* Contact maps can be **zoomed in/out** from finest resolution to whole chromosome level.
* Rich customizations of **color scale for contact maps** visualization
* Rich customizations of **X- and Y- axis properties**.

* `Normalization of contact maps`_ by
    * **Iterative Correction** (IC)
    * **Knight-Ruiz Matrix Balancing** (KR)
    * **Distance-Frequency**

* A **new file format** based on HDF5 for `genome contact map`_ and `genomic track datasets`_.
    * **Portable**, **platform independent** and can be read through C/C++, JAVA, Python and R programming language.
    * **Very fast to read** - fast browsing of contact maps and genomic datasets

* Another file format for `chormosomal contact map <http://gcmapexplorer.readthedocs.io/en/latest/about_ccmap_file.html>`_
  - much faster than above format to read/write but not compact. Suitable for performing calculations.
* `A GUI interface and commands <http://gcmapexplorer.readthedocs.io/en/latest/about_gcmap_file.html#convert-hi-c-data-to-gcmap>`_
  to convert Coordinate Sparse, Pair Coordinate Sparse, HOMER Interaction matrix, Bin-Contact formats into the new gcmap and ccmap formats.
* Interface and commands to `convert bigWig/wig/bed file to genomic track dataset h5 file`_.
* Interface and commands for `Normalization of contact maps`_.
* Publication ready images at one click.


----


Interfaces and Commands
-----------------------

Usage
-----

Run ``gcMapExplorer`` command on terminal to get list of all sub-commands.

**Following sub-commands are available:**

.. list-table:: Graphical User Interface Applications
    :widths: 1, 4
    :header-rows: 1
    :name: gui-table

    * - Command
      - Function

    * - browser_
      - Interactive Browser for genomic contact maps

    * - `cmapImporter <http://gcmapexplorer.readthedocs.io/en/latest/commands/cmapImporter.html>`_
      - Interface to import contact maps and datasets

    * - `cmapNormalizer <http://gcmapexplorer.readthedocs.io/en/latest/commands/cmapNormalizer.html>`_
      - Interface to normalize contact maps

    * - `h5Converter <http://gcmapexplorer.readthedocs.io/en/latest/commands/h5Converter.html>`_
      - Interface to convert bigWig/wig/bed file to h5 file


.. list-table::  Commands to import Hi-C data
    :widths: 1, 4
    :header-rows: 1
    :name: import-hic-command-table

    * - Command
      - Function

    * - `coo2cmap <http://gcmapexplorer.readthedocs.io/en/latest/commands/coo2cmap.html>`_
      - Import COO sparse matrix format to ccmap or gcmap

    * - `pairCoo2cmap <http://gcmapexplorer.readthedocs.io/en/latest/commands/pairCoo2cmap.html>`_
      - Import map from files similar to paired COO format

    * - `homer2cmap <http://gcmapexplorer.readthedocs.io/en/latest/commands/homer2cmap.html>`_
      - Import HOMER Hi-C interaction matrix to ccmap or gcmap

    * - `bc2cmap <http://gcmapexplorer.readthedocs.io/en/latest/commands/bc2cmap.html>`_
      - Import Bin-Contact format files to ccmap or gcmap


.. list-table:: Commands to convert bigWig/wig/bed to h5
    :widths: 1, 4
    :header-rows: 1
    :name: convert-to-h5-file-table

    * - Command
      - Function

    * - `bigwig2h5 <http://gcmapexplorer.readthedocs.io/en/latest/commands/bigWig2h5.html>`_
      - Convert a bigWig file to HDF5 format h5 file

    * - `wig2h5 <http://gcmapexplorer.readthedocs.io/en/latest/commands/wig2h5.html>`_
      - Convert a wig file to HDF5 format h5 file

    * - `bed2h5 <http://gcmapexplorer.readthedocs.io/en/latest/commands/bed2h5.html>`_
      - Convert a bed file to HDF5 format h5 file

    * - `encode2h5 <http://gcmapexplorer.readthedocs.io/en/latest/commands/encode2h5.html>`_
      - Download and convert ENCODE datasets to HDF5 format h5 files


.. list-table:: Commands to normalize Hi-C map
    :widths: 1, 4
    :header-rows: 1
    :name: normalize-maps-table

    * - Command
      - Function

    * - `normKR <http://gcmapexplorer.readthedocs.io/en/latest/commands/normKR.html>`_
      - Normalization using Knight-Ruiz matrix balancing

    * - `normIC <http://gcmapexplorer.readthedocs.io/en/latest/commands/normIC.html>`_
      - Normalization using Iterative Correction

    * - `normMCFS <http://gcmapexplorer.readthedocs.io/en/latest/commands/normMCFS.html>`_
      - Scale maps using Median/Mean Contact Frequency


.. list-table:: Commands for Analysis
    :widths: 1, 4
    :header-rows: 1

    * - Command
      - Function

    * - corrBWcmaps
      - Calculate correlation between contact maps


Command help
------------
Run ``gcMapExplorer <sub-commands> -h`` command.

For example:
	* ``gcMapExplorer normKR -h``
	* ``gcMapExplorer coo2cmap -h``
