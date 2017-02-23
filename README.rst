Genome Contact Map Explorer - gcMapExplorer
===========================================

It is a platform to visualize and analyze the contact maps that are generated from Hi-C experiments. This package is developed by considering the huge size of contact maps at very fine resolution. It contains

  * Graphical User Interface - Several windows like applications to perform tasks.
  * Command Line Interface - Several commands to perform tasks.
  * Application Programming Interface - It can be used to perform analysis by any mathematical operations through programming.


**For more details, visit:** `gcMapExplorer Homepage <http://gcmapexplorer.readthedocs.io/>`_

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
* Interface for data conversion
* Interface for Normalization
* Publication ready images at one click.


----


Interfaces and Commands
-----------------------

Usage
~~~~~

``gcMapExplorer [Command]``

Run ``gcMapExplorer`` command on terminal to get list of all sub-commands.

Following sub-commands are available:

Graphical User Interface
~~~~~~~~~~~~~~~~~~~~~~~~
* **browser** : Interactive Browser for genomic contact maps
* **importer** : Interface to import contact maps and datasets
* **normalizer** : Interface to normalize contact maps

Commands to convert or import data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **coo2cmap** : Import COO sparse matrix format to ccmap or gcmap
* **pairCoo2cmap** : Import map from files similar to paired COO format
* **homer2cmap** : Import HOMER Hi-C interaction matrix to ccmap or gcmap
* **bc2cmap** : Import Bin-Contact format files to ccmap or gcmap
* **bigwig2h5** : Import a bigWig file to HDF5 format h5 file

Commands to normalize contact map
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **normKR** : Normalization using Knight-Ruiz matrix balancing
* **normIC** : Normalization using Iterative Correction
* **normMCFS** : Normalization by Median Contact Frequency Scaling

Commands for Analysis
~~~~~~~~~~~~~~~~~~~~~
* **corrBWcmaps** : Calculate correlation between contact maps
