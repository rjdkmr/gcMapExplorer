coo2cmap - convert COO sparse matrix format to ccmap or gcmap
-------------------------------------------------------------

As shown below in example, in this format, first and second column is location
on chromosome and third column is the respective value:

::

    20000000        20000000        2692.0
    20000000        20100000        885.0
    20100000        20100000        6493.0
    20000000        20200000        15.0
    20100000        20200000        52.0
    20200000        20200000        2.0
    20000000        20300000        18.0
    20100000        20300000        40.0

**NOTE** that, above location is real value. However, with ``-idx/--index``
option, these two same column willbe considered as index value. index should 
always start from zero for absoulte beginning of chromosome.e.g. for 10kb, 
0-10000 should have index of zero, 10000-20000 have index of one. If this 
is file format,resolution should be provided with ``-r/--resolution`` option.

Usage:
    .. code-block:: bash

        usage: gcMapExplorer coo2cmap [-h] [-i input.txt] [-ic input.tar.gz]
                                      [-mt intra] [-r 10kb] [-idx]
                                      [-ccm 10kb_RawObserved] [-od OUTDIR]
                                      [-gcm inOut.gcmap] [-cmeth lzf] [-dmeth sum]
                                      [-wd /home/rajendra/deskForWork/scratch]


**Optional arguments:**

.. code-block:: bash

  -h, --help            show this help message and exit
  -i input.txt, --input input.txt
                        Meta input file containing input contact map files list with respective
                        xlabel and ylabel. xlabel should be always provided. In case of intra-
                        chromosomal map, only xlabel is sufficient because both x and y axis are of
                        same chromosome. However for inter-chromosomal map, both xlabel and ylabel
                        should be provided. Example format:
                        
                        100kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_100kb.RAWobserved    chr1
                        100kb_resolution_intrachromosomal/chr5/MAPQGE30/chr5_100kb.RAWobserved    chr5
                        100kb_resolution_intrachromosomal/chr15/MAPQGE30/chr15_100kb.RAWobserved  chr15
                        100kb_resolution_intrachromosomal/chr20/MAPQGE30/chr20_100kb.RAWobserved  chr20
                        100kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_100kb.RAWobserved  chr21
                        100kb_resolution_intrachromosomal/chr22/MAPQGE30/chr22_100kb.RAWobserved  chr22
                        
  -ic input.tar.gz, --input-compressed input.tar.gz
                        Input compressed archive file containing all the listed contact maps.
                        Presently, only "tar.gz" and "zip" compressed files are supported.
                        
                        If -i/--input is not provided, all files from compressed file will be tried for
                        processing.
                        
  -mt intra, --mapType intra
                         Type of listed contact maps: "intra" or "inter" chromosomal map.
                        
  -r 10kb, --resolution 10kb
                        Resolution of all maps. It is an optional argument. Note that, if this
                        option is not provided, resolution will be automatically determined from the
                        contact map file. However, in case of -idx/--index option, resolution
                        should be provided as resolution cannot be determined from input contact map
                        file.
                        
  -idx, --index         It determines whether contact map files have real coordinate of chromosome
                        or index number. If this option is enabled, -r/--resolution option should be
                        provided.
                        
  -ccm 10kb_RawObserved, --ccmap 10kb_RawObserved
                         Use this to convert all contact maps to ccmap format files. Provide suffix
                        of ccmap file names with this option and it will enable the conversion.
                        
                        Ouput ccmap file name is generated outmatically as follows;
                        if xlabel is not equal to ylabel: <xlabel>_<ylabel>_<suffix>.ccmap
                        else: <xlabel>_<suffix>.ccmap
                        
                        Note that -od/--out-dir option is also required because all ccmaps will be
                        saved in this directory.
                        
  -od OUTDIR, --out-dir OUTDIR
                        Directory where all ccmap files will be saved.
  -gcm inOut.gcmap, --gcmap inOut.gcmap
                        Provide gcmap file to convert all contact maps into one gcmap file.
                        File name should contain full path because -od/--out-dir is not considered
                        for this conversion.
                        
  -cmeth lzf, --compression-method lzf
                        Data compression method in gcmap file.
  -dmeth sum, --downsample-method sum
                        Downsampling method to coarsen the resolution in gcmap file. The option is
                        intended to use with -gcm/--gcmap option. Three accepted methods are
                                sum  : sum of values,
                                mean : Average of values and
                                max  : Maximum of values.
                        
                        This option generates all coarser maps where resolutions will be coarsened by
                        a factor of two, consequetively. e.g.: In case of 10 kb input resolution,
                        downsampled maps of "20kb", "40kb", "80kb", "160kb", "320kb" etc. will be
                        generated until, map size is less than 500.
                        
  -wd /home/rajendra/deskForWork/scratch, --work-dir /home/rajendra/deskForWork/scratch
                        Directory where temporary files will be stored.

