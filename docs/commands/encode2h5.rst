encode2h5
~~~~~~~~~~

**Description:**

.. code-block:: bash

  Download and Convert ENCODE datasets to h5 files
  =================================================

  It can be used to download and convert multiple datasets from ENCODE Experiment
  matrix (https://www.encodeproject.org/matrix/?type=Experiment).
  Presently, only bigWig files are downloaded and then converted.

  At first search the datasets on https://www.encodeproject.org/matrix/?type=Experiment .
  Then click on download button on top of the page. A text file will be downloaded.
  This text file can be used as input in this program. All bigWig files will be
  downloaded and converted to gcMapExplorer compatible hdf5 format.

  NOTE: At first a metafile is automatically downloaded and then files
        are filtered according to bigWig format and Assembly. Subsequently,
        if several replicates are present, only datasets with combined
        replicates are considered. In case if two replicates are present
        and combined replicates are not present, replicates will be combined with
        '-mtc/--method-to-combine' option.

  NOTE: Because downloading and conversion might take very long time, it also
        generates a checkpoint file in the output directory. Therefore,
        in case of crash or abrupt exit, the process can be continued from the
        last file.

  Name of output files:

      (1) For ChIP-seq assay:
          a. signal-<Experiment target>-<Experiment accession>-<File accessions.h
          b. fold-<Experiment target>-<Experiment accession>-<File accessions>.h5

      (2) For RNA-seq:
          a. uniq-reads-<date>-<Experiment accession>-<File accessions>.h
          b. plus-uniq-reads-<date>-<Experiment accession>-<File accessions>.h
          c. minus-uniq-reads-<date>-<Experiment accession>-<File accessions>.h
          d. all-reads-<date>-<Experiment accession>-<File accessions>.h5
          e. plus-all-reads-<date>-<Experiment accession>-<File accessions>.h5
          f. minus-all-reads-<date>-<Experiment accession>-<File accessions>.h5
          g. signal-<date>-<Experiment accession>-<File accession>.h5

      (2) For DNase-seq:
          a. uniq-reads-signal-<date>-<Experiment accession>-<File accessions>.h
          b. raw-signal-<date>-<Experiment accession>-<File accessions>.h
          c. all-reads-signal-<date>-<Experiment accession>-<File accessions>.h
          d. signal-<date>-<Experiment accession>-<File accessions>.h5

      Note that name of cell-line is not included here. Therefore, use the
      directory name as a identfiers for cell-lines or species. The Experiment
      and File accession can be used to back-track about the dataset on ENCODE
      website.

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

    usage: gcMapExplorer encode2h5 [-h] [-i input.txt] [-amb hg19] [-asy ChIP-seq]
                                   [-b2w bigWigToWig] [-binfo bigWigInfo]
                                   [-r "List of Resolutions"]
                                   [-dm "List of downsampling method"]
                                   [-cmeth lzf] [-mtc mean] [-od outDir] [-ko]
                                   [-wd /home/rajendra/deskForWork/scratch]


**Optional arguments:**

.. code-block:: bash

  -h, --help            show this help message and exit
  -i input.txt, --input input.txt
                        Input text file.
                        At first search the datasets on https://www.encodeproject.org/matrix/?type=Experiment.
                        Then click on download button on top of the page. A text file will be downloaded.
                        This text file can be used as input in this program.

  -amb hg19, --assembly hg19
                         Name of reference genome.
                        Example: hg19, GRCh38 etc.

  -asy ChIP-seq, --assay ChIP-seq
                         Name of assay.
                        Presently, four assays are implemented:
                        'ChIP-seq', 'RNA-seq', 'DNase-seq' and 'FAIRE-seq'.

  -b2w bigWigToWig, --bigWigToWig bigWigToWig
                        Path to bigWigToWig tool.

                        This is not neccessary when bigWigToWig path is already set using gcMapExplorer
                        configure utility.

                        It can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/
                        for linux and Mac platform.

                        If it is not present in configuration file, the input path should
                        be provided. It will be stored in configuration file for later use.

  -binfo bigWigInfo, --bigWigInfo bigWigInfo
                         Path to bigWigInfo tool.

                        This is not neccessary when bigWigInfo path is already set using gcMapExplorer
                        configure utility.

                        It can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/
                        for linux and Mac platform.

                        If it is not present in configuration file, the input path should
                        be provided. It will be stored in configuration file for later use.

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

  -cmeth lzf, --compression-method lzf
                        Data compression method in h5 file.
  -mtc mean, --method-to-combine mean
                        Methods to combine data from more than two input file. Presently, three
                        methods can be used: 'mean', 'max' and 'min' for average, maximum and minimum
                        value, respectively.

  -od outDir, --outDir outDir
                         Directory to save all h5 files. It is an essential input.

  -ko, --keep-original  To copy original 1-base resolution data in h5 file. This will increase the
                        file size significantly.

  -wd /home/rajendra/deskForWork/scratch, --work-dir /home/rajendra/deskForWork/scratch
                        Directory where temporary files will be stored.
