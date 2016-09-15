.. |hdf5| raw:: html

  <a href="https://www.hdfgroup.org/HDF5/" target="_blank"> HDF5 </a>

.. |lzf| raw:: html

    <a href="http://www.h5py.org/lzf/" target="_blank"> LZF </a>

.. |h5| raw:: html

    <a href="https://github.com/mannau/h5" target="_blank"> h5 </a>

.. |rhdf5| raw:: html

    <a href="http://bioconductor.org/packages/release/bioc/html/rhdf5.html" target="_blank"> rhdf5 </a>

``gcmap`` file
==============

The ``gcmap`` file is in |hdf5| format to store Genome Contact Map (gcmap). It is implemented by considering both portability and readability. A single file may contains maps of
various resolutions of all chromosomes. It also contains properties, which are attributes, of each map.


Structure of ``gcmap`` file
---------------------------
As |hdf5| format supports Hierarchical Data Model, therefore we implemented the contact maps in format way.
Overall structure format is as follows:

::

    HDF5
    │
    ├──────── chr1 ──── Attributes : ['xlabel', 'ylabel', 'compression']
    │           │
    │           ├────── 10kb  ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 20kb  ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 40kb  ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 60kb  ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 80kb  ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 160kb ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 320kb ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 640kb ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           │
    │           ├────── 10kb-bNoData  ( 1D Numpy Array )
    │           ├────── 20kb-bNoData  ( 1D Numpy Array )
    │           ├────── 40kb-bNoData  ( 1D Numpy Array )
    │           ├────── 60kb-bNoData  ( 1D Numpy Array )
    │           ├────── 80kb-bNoData  ( 1D Numpy Array )
    │           ├────── 160kb-bNoData ( 1D Numpy Array )
    │           ├────── 320kb-bNoData ( 1D Numpy Array )
    │           └────── 640kb-bNoData ( 1D Numpy Array )
    │
    ├──────── chr2 ──── Attributes : ['xlabel', 'ylabel', 'compression']
    │           │
    │           ├────── 10kb  ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 20kb  ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 40kb  ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 60kb  ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 80kb  ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 160kb ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 320kb ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           ├────── 640kb ( 2D Numpy Array ) ─── Attributes : ['minvalue', 'maxvalue', 'xshape', 'yshape', 'binsize']
    │           │
    │           ├────── 10kb-bNoData  ( 1D Numpy Array )
    │           ├────── 20kb-bNoData  ( 1D Numpy Array )
    │           ├────── 40kb-bNoData  ( 1D Numpy Array )
    │           ├────── 60kb-bNoData  ( 1D Numpy Array )
    │           ├────── 80kb-bNoData  ( 1D Numpy Array )
    │           ├────── 160kb-bNoData ( 1D Numpy Array )
    │           ├────── 320kb-bNoData ( 1D Numpy Array )
    │           └────── 640kb-bNoData ( 1D Numpy Array )
    :
    :
    :
    └───── ...


Compression
-----------
In gcmap file, contact map is stored as compressed 2D matrix. Presently, two compression method are allowed in the gcmap file:

* |lzf|
* GZIP

By default, |lzf| is used to compress arrays. This method is very fast, and allow the rapid contact map reading.
However, the size reduction is moderate in comparison with GZIP compression method.

.. Warning::
  |lzf| method is only avaiable through **Python h5py** module, and therefore, this file cannot be read by another programming language through standard library.
  For portablity, use GZIP compression method, which is available in standard HDF5 library.


Portability and Readability
---------------------------
The ``gcmap`` file with **GZIP** compressed arrays can be read and write from any programming language. For C/C++/Java, a standard HDF5 library is available from |hdf5| group.
For R programming language, |h5| and |rhdf5| are available.

Both GZIP and |lzf| compression reduces the file size significantly as compare to respective flat text file. Therefore, this file is also
suitable for storage and transfer.
