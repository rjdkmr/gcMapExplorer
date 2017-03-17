homer2cmap - HOMER Hi-C matrix to ccmap or gcmap
------------------------------------------------

`HOMER package <http://homer.salk.edu/homer/interactions/>`_ contains modules to
analyze genome wide interaction data. It creates Hi-C matrix in a specific
format as shown as shown `here <http://homer.salk.edu/homer/interactions/HiCmatrices.html>`_.

This format contains contact map in a matrix format.

Usage:
    .. code-block:: bash

        usage: gcMapExplorer homer2cmap [-h] [-i matrix.txt] [-ccm RawObserved]
                                        [-od OUTDIR] [-gcm inOut.gcmap] [-cmeth lzf]
                                        [-dmeth sum]
                                        [-wd /home/rajendra/deskForWork/scratch]



**Optional arguments:**

.. code-block:: bash

  -h, --help            show this help message and exit
  -i matrix.txt, --input matrix.txt
                         File containing HOMER Hi-C ineraction matrix format contact map
                        
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

