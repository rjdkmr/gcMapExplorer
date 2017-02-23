.. |numpy memmap| raw:: html

	<a href="http://docs.scipy.org/doc/numpy/reference/generated/numpy.memmap.html" target="_blank">numpy memmap</a>


.. |numpy array| raw:: html

  	<a href="http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.array.html" target="_blank">numpy array</a>


.. |numpy array concept| raw:: html

  	<a href="http://www.scipy-lectures.org/intro/numpy/array_object.html" target="_blank">here</a>

.. |indexing and slicing| raw:: html

  	<a href="http://www.scipy-lectures.org/intro/numpy/array_object.html#indexing-and-slicing" target="_blank">indexing and slicing</a>

.. |numpy routines| raw:: html

  <a href="http://docs.scipy.org/doc/numpy/reference/routines.html" target="_blank">numpy</a>

.. |scipy routines| raw:: html

  <a href="http://docs.scipy.org/doc/scipy/reference/" target="_blank">scipy</a>



``*.ccmap`` and ``*.npbin`` files
---------------------------------

This package implements the Chromosome Contact Map (ccmap) data in a specific format, which contains two inter-related files,
  * ``*.ccmap``: It is a text file and contains meta-data for the Chromosome Contact Map.
  * ``*.npbin`` or ``*.npbin.gz``: This file contains the Chromosome Contact Map data, which is a memory mapped 2D matrix.


In contrast to ``gcmap`` file, the ``*.ccmap`` and ``*.npbin`` paired files contain only one contact map. To perform mathematical
operations directly using ``gcmap`` is slow
(`see here <modules_examples/access_gcmap_data.html#Execution-Time-Comparison-between-np.ndarray,-ccmap.matrix-and-gcmap.matrix>`_),
therefore, gcMapExplorer uses ``*.ccmap`` and ``*.npbin`` paired files during mathematical calculations.


Why two files?
==============
Contact map is a two-dimensional matrix. Size of matrix can be very huge for large chromosome at high resolutions.
Memory required to handle such huge matrices could be very large and sometimes beyond the available hardware.
Therefore, we used a matrix mapped to the file ``*.npbin`` (compressed ``*.npbin.gz``), which is stored in external disk.
For each Contact map, we also need to store some properties like its title/name, size, minimum and maximum values, columns/rows with missing data, and
path to memory mapped matrix file. These properties are stored in ``*.ccmap`` file.

Advantages of memory mapped matrix file
=======================================
* It is a binary indexed file and any particular region of the matrix can be rapidly accessed.
* This file is generated using |numpy memmap|, and therefore, it can be used as |numpy array|. See also |numpy array concept| for more about numpy array.
* Because, it can be used as a numpy array, |indexing and slicing| operations can be performed to access the data.
* All mathematical operations available in |numpy routines| and |scipy routines| modules can be directly performed.


Contents of ``*.ccmap`` file
=============================

``*.ccmap`` is a text file and its content is shown as example:

::

  {
      "title":null,
      "path2matrix":"chr22_100kb_normKR.npbin",
      "xlabel":null,
      "minvalue":"7.207987891888479e-06",
      "bLog":false,
      "state":"saved",
      "maxvalue":"0.28213343024253845",
      "binsize":100000,
      "shape":[
          "513",
          "513"
      ],
      "matrix":null,
      "yticks":[
          "0",
          "51300000"
      ],
      "bNoData":"111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111100000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
      "xticks":[
          "0",
          "51300000"
      ],
      "ylabel":null,
      "dtype":"float32"
  }


Following properties are included in Contact Map metadata file:
  * ``title``: Title of the data. Used to display in browser.
  * ``path2matrix``: Path to ``*.npbin`` or ``*.npbin.gz`` file
  * ``bLog``: Whether values in matrix is Logarithm values
  * ``maxvalue``: Maximum value
  * ``minvalue``: Minimum value
  * ``binsize``: Resolution of data.
  * ``shape``: Shape of matrix along X and Y axis
  * ``xticks``: Upper and lower limits of X-axis
  * ``yticks``: Upper and lower limits of Y-axis
  * ``xlabel``: Label for x-axis
  * ``ylabel``: Label for y-axis
  * ``matrix``: See :attr:`gcMapExplorer.lib.CCMAP.matrix`
  * ``state``: See :attr:`gcMapExplorer.lib.CCMAP.state`
  * ``dtype``: Data type for memory mapped matrix file. e.g. float, float32, float64 etc.
  * ``bNoData``: Whether data is missing for entire row/column
