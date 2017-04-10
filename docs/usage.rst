

How to use gcMapExplorer?
=========================

Several interfaces are available as following.

* :ref:`gui-table` - Several windows like applications to perform tasks.
* Command Line Interface - Several commands to perform tasks.
    * :ref:`import-hic-command-table`
    * :ref:`convert-to-h5-file-table`
    * :ref:`normalize-maps-table`
* `Application Programming Interface <http://gcmapexplorer.readthedocs.io/en/latest/apidoc/summary.html>`_
  - It can be used to perform analysis by any mathematical operations through programming.

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

    * - `browser <http://gcmapexplorer.readthedocs.io/en/latest/mapBrowser.html>`_
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

    * - `encodeToH5 <http://gcmapexplorer.readthedocs.io/en/latest/commands/encodeToH5.html>`_
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
