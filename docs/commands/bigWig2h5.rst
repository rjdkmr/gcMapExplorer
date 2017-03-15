
**Description:**

.. code-block:: bash

  Import a bigWig file to HDF5 format h5 file
  ============================================

  bigWig file can be converted into gcMapExplorer compatible HDF5 file using
  this tool. This HDF5 file can be loaded into gcMapExplorer browser for
  interactive visualization.

  Requirements
  ============
  1) bigWigToWig : It converts binary bigWig file to ascii Wig file.
  2) bigWigInfo : It fetches the information about chromosomes from bigWig file.

  Both tools can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/
  for linux and Mac platform. However, these tools are not yet available for
  Windows OS.

  Path to these tools can be set using gcMapExplorer configure utility or can be
  given with the command.

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

  usage: gcMapExplorer bigwig2h5 [-h] [-i input.bigWig] [-t "Genomic Dataset"]
                                 [-b2w bigWigToWig] [-binfo bigWigInfo]
                                 [-r "List of Resolutions"] [-icn CHROMNAME]
                                 [-dm "List of downsampling method"]
                                 [-cmeth lzf] [-o out.h5] [-ow] [-ko]
                                 [-wd /home/rajendra/deskForWork/scratch]

**Optional arguments:**

.. code-block:: bash

  -h, --help            show this help message and exit
  -i input.bigWig, --input input.bigWig
                        Input bigWig file.

  -t "Genomic Dataset", --title "Genomic Dataset"
                        Title of the dataset.
  -b2w bigWigToWig, --bigWigToWig bigWigToWig
                        Path to bigWigToWig tool.

                        This is not neccessary when bigWigToWig path is already set using gcMapExplorer
                        configure utility.

                        It can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/
                        for linux and Mac platform.

  -binfo bigWigInfo, --bigWigInfo bigWigInfo
                         Path to bigWigInfo tool.

                        This is not neccessary when bigWigInfo path is already set using gcMapExplorer
                        configure utility.

                        It can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/
                        for linux and Mac platform.

  -r "List of Resolutions", --resolutions "List of Resolutions"
                        Additional input resolutions other than these resolutions: 1kb', '2kb',
                        '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb', '100kb', '160kb','200kb',
                        '320kb', '500kb', '640kb',  and '1mb'.

                        Resolutions should be provided in comma seprated values. For Example:
                        -r "25kb, 50kb, 75kb"

  -icn CHROMNAME, --input-chromosome CHROMNAME
                        Input Chromosome Name.
                        If this is provided, only this chromosome data is extracted and stored in h5
                        file.

  -dm "List of downsampling method", --downsample-method "List of downsampling method"
                        Methods to coarse or downsample the data for converting from 1-base
                        to coarser resolutions. If this option is not provided, all six methods (see
                        above) will be considered. User may use only subset of these methods.
                        For example: -dm "max, amean" can be used for downsampling by only these
                        two methods.

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

  -wd /home/rajendra/deskForWork/scratch, --work-dir /home/rajendra/deskForWork/scratch
                        Directory where temporary files will be stored.
