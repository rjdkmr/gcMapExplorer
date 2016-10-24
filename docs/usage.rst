

How to use gcMapExplorer?
=========================

Several interfaces are available as following.

* Graphical User Interface - Several windows like applications to perform tasks.
* Command Line Interface - Several commands to perform tasks.
* Application Programming Interface - It can be used to perform analysis by any mathematical operations through programming.

Usage
-----

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

To Launch Browser
-----------------
Run ``gcMapExplorer browser`` command.

To Launch Normalizer
--------------------
Run ``gcMapExplorer normalizer`` command.

To Launch Importer
------------------
Run ``gcMapExplorer importer`` command.

To get help for a command
-------------------------
Run ``gcMapExplorer <sub-commands> -h`` command. For examples:

* ``gcMapExplorer normKR -h``
* ``gcMapExplorer coo2cmap -h``
