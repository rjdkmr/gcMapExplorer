bc2cmap -Bin-Contact files pair to ccmap or gcmap
-------------------------------------------------

In this format, two separate files are available. One file contains bins
information and other contains contact frequency.

These types of files are present in following GEO data:
    * http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61471
    * http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34453

This format contains a pair of file:
    BIN file:
        ::
    
            cbin	chr	from.coord	to.coord	count
            1	2L	0	160000	747
            2	2L	160000	320000	893
            3	2L	320000	480000	1056
            4	2L	480000	640000	1060
            5	2L	640000	800000	978
            6	2L	800000	960000	926
            .
            .
            .



    CONTACT file in list format:
        ::
    
            cbin1	cbin2	expected_count	observed_count
            1	1	40.245201	21339
            1	2	83.747499	5661
            1	3	92.12501	1546
            1	4	93.401273	864
            1	5	87.265472	442
            .
            .
            .



Both BIN and CONTACT files are **neccessary** for the conversion.


Usage:
    .. code-block:: bash

        usage: gcMapExplorer bc2cmap [-h] [-ib nm_none_160000.bins]
                                     [-ic nm_none_160000.n_contact] [-ccm RawObserved]
                                     [-od OUTDIR] [-gcm inOut.gcmap] [-cmeth lzf]
                                     [-dmeth sum]
                                     [-wd /home/rajendra/deskForWork/scratch]



**Optional arguments:**

.. code-block:: bash

  -h, --help            show this hel
  p message and exit
  -ib nm_none_160000.bins, --input-bin nm_none_160000.bins
                         Input BIN file as shown above.
                        
  -ic nm_none_160000.n_contact, --input-contact nm_none_160000.n_contact
                         Input CONTACT file as shown above.
                        
  -ccm RawObserved, --ccmap RawObserved
                         Use this to convert all contact maps to ccmaps file. Provide suffix of
                        ccmap file names with this option and it will enable the conversion.
                        
                        Ouput ccmap file name is generated outmatically as follows;
                        <chromosome>_<resolution>_<suffix>.ccmap
                        
                        Note that -od/--out-dir option is also required because all ccmaps will be
                        saved in this directory.
                        
  -od OUTDIR, --out-dir OUTDIR
                        Directory where all ccmap files will be saved.
  -gcm inOut.gcmap, --gcmap inOut.gcmap
                        Provide gcmap file to convert all contact maps into one gcmap file.
                        File name should contain full path because -od/--out-dir is not considered
                        for thi conversion.
                        
  -cmeth lzf, --compression-method lzf
                        Data compression method for gcmap file.
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


