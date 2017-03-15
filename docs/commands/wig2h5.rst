
**Description:**

.. code-block:: bash

  Import a wig file to HDF5 format h5 file
  ============================================

  wig file can be converted into gcMapExplorer compaitable HDF5 file using
  this tool. This HDF5 file can be loaded into gcMapExplorer browser for
  interactive visualization.

  This tool does not require any external program.

  Resolutions
  ===========
  By default, original data are downsampled to following resolutions:  '1kb',
  '2kb', '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb', '100kb', '160kb',
  '200kb', '320kb', '500kb', '640kb',  and '1mb'.

  The data are downsampled at this stage only to speed up the visualization
  process as downsampling might slow down the interactive visualization.

  Downsampling/Coarsening method
  ==============================
  Presently, six methods are implemented:
  1) min    -> Minimum value
  2) max    -> Maximum value
  3) amean  -> Arithmatic mean or average
  4) hmean  -> Harmonic mean
  5) gmean  -> Geometric mean
  6) median -> Median

  All these methods are used by default.
  See below help for "-dm/--downsample-method" option to change the methods.

  To keep original 1 base resolution data
  =======================================
  By default, the output h5 file does not contain original 1-base resolution
  data to reduce the file size. To keep the original data in h5 file, used
  -ko/--keep-original flag.



**Usage:**

.. code-block:: bash

  usage: gcMapExplorer wig2h5 [-h] [-i input.wig] [-t "Genomic Dataset"]
                              [-r "List of Resolutions"]
                              [-dm "List of downsampling method"]
                              [-icn CHROMNAME] [-cmeth lzf] [-o out.h5] [-ow]
                              [-ko] [-idf index.json]
                              [-wd /home/rajendra/deskForWork/scratch]



**Optional arguments:**

.. code-block:: bash

  -h, --help            show this help message and exit
  -i input.wig, --input input.wig
                        Input wig file.

  -t "Genomic Dataset", --title "Genomic Dataset"
                        Title of the dataset.
  -r "List of Resolutions", --resolutions "List of Resolutions"
                        Additional input resolutions other than these resolutions: 1kb', '2kb',
                        '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb', '100kb', '160kb','200kb',
                        '320kb', '500kb', '640kb',  and '1mb'.

                        Resolutions should be provided in comma seprated values. For Example:
                        -r "25kb, 50kb, 75kb"

  -dm "List of downsampling method", --downsample-method "List of downsampling method"
                        Methods to coarse or downsample the data for converting from 1-base
                        to coarser resolutions. If this option is not provided, all six methods (see
                        above) will be considered. User may use only subset of these methods.
                        For example: -dm "max, amean" can be used for downsampling by only these
                        two methods.

  -icn CHROMNAME, --input-chromosome CHROMNAME
                        Input Chromosome Name.
                        If this is provided, only this chromosome data is extracted and stored in h5
                        file.

  -cmeth lzf, --compression-method lzf
                        Data compression method in h5 file.
  -o out.h5, --out out.h5
                        Output h5 file.

                        If file is already present, it will replace the data. Therefore, be careful
                        if a file with same name is present.

  -ow, --overwrite      If a output file is already present, overwrite the datasets in the output
                        file.

  -ko, --keep-original  To copy original 1-base resolution data in h5 file. This will increase the
                        file size significantly.

  -idf index.json, --index-file index.json
                        Index file in json format.
                        A file in json format containing indices (position in wig file) and sizes of
                        chromosomes. If this file is not present and given as input, a new file will be
                        generated. If this file is present, indices andsizes will be taken from this
                        file. If index and size of input chromosome is not present in json file, these
                        will be determined from wig file and stored in same json file. This file could
                        be very helpful in case when same wig file has to be read many times because
                        step to determine index and size of chromosome is skipped.

  -wd /home/rajendra/deskForWork/scratch, --work-dir /home/rajendra/deskForWork/scratch
                        Directory where temporary files will be stored.
