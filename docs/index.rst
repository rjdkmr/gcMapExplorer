.. hiCMapAnalyze documentation master file, created by
   sphinx-quickstart on Wed Sep 30 22:10:12 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Genome Contact Map Explorer - gcMapExplorer
===========================================

It is a platform to visualize and analyze the contact maps that are generated from Hi-C experiments. This package is developed by considering the huge size of contact maps at very fine resolution. It contains

  * Graphical User Interface - Several windows like applications to perform tasks.
  * Command Line Interface - Several commands to perform tasks.
  * Application Programming Interface - It can be used to perform analysis by any mathematical operations through programming.

Features:
---------

  * Support for **huge contact maps** - Use of Disk instead of RAM - Matrices/arrays are stored in Disks - mathematical operations by directly reading/writing from/to Disks, **without loading them into RAM**
  * A browser with rich interfaces for **Comparative** and **Interactive** visualization of **two dimensional contact maps** along with **genomic datasets** such as produced by DNase-seq, ChIP-seq, RNA-seq etc.
  * Contact maps can be **zoomed in/out** from finest resolution to whole chromosome level.
  * Rich customizations of **color scale for contact maps** visualization
  * Rich customizations of **X- and Y- axis properties**.

  * Normalization of contact maps by

    * **Iterative Correction** (IC)
    * **Knight-Ruiz Matrix Balancing** (KR)
    * **Distance-Frequency**

  * A **new file format** for contact map  and genomic datasets:

    * **Portable**, **platform independent** and can be read through C/C++, JAVA, Python and R programming language.
    * **Very fast to read** - fast browsing of contact maps and genomic datasets

  * Another file format for chormosomal contact map - much faster than above format to read/write but not compact
  * Easy import of Coordinate Sparse, HOMER Interaction matrix and Bin-Contact formats to the new formats.
  * Interface for file format conversions
  * Interface for Normalization
  * Publication ready images at one click.



Screen-shots
------------

.. figure:: images/browser.png
      :scale: 40%
      :alt: Screen snapshot of Genome Contact Map browser

      Genome Contact Map browser

.. figure:: images/axis.png
      :scale: 35%
      :alt: Screen snapshot of Axis Properties interface in Browser

      Axis Properties interface in Browser


.. figure:: images/importer.png
      :scale: 35%
      :alt: Screen snapshot of Importer Interface

      gcmap Importer Interface


.. figure:: images/h5Converter.png
      :scale: 35%
      :alt: Screen snapshot of genomic track dataset converter

      genomic track dataset converter Interface


.. figure:: images/normalizer.png
      :scale: 35%
      :alt: Screen snapshot of normalizer interface

      Contact map normalization Interface

****


Contents
========

.. py:module:: gcMapExplorer

.. toctree::
   :maxdepth: 2

   Requirements and Installation <install>
   How to use gcMapExplorer? <usage.rst>
   About gcmap file <about_gcmap_file>
   About ccmap and npbin files <about_ccmap_file>
   About Genomic track h5 file <about_h5_file>
   Download example datasets <dLcmaps>
   Summary of Python Modules <apidoc/summary>
   Examples using Python Modules <modules_examples/index>
   Python Modules documentation <apidoc/index>


Indices
=======

* :ref:`genindex`
* :ref:`modindex`
