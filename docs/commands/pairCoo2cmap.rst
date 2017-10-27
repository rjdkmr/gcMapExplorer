pairCoo2cmap - paired COO sparse matrix to ccmap or gcmap
---------------------------------------------------------

This format is very similar to COO format with an additional information of
chromosome. Therefore, maps for all chromosome could be contained in a single
file.

This type of format appeared with this `publication
<http://dx.doi.org/10.1016/j.cell.2015.10.026>`_ 
(`GEO datasets <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72512>`_).

Following file format can be read as a text file, where first and second
column is location on chromosome and third column is the value:

::

    chr4     60000   75000   chr4    60000   75000   0.1163470887070292
    chr4     60000   75000   chr4    105000  120000  0.01292745430078102
    chr4     60000   75000   chr4    435000  450000  0.01292745430078102
    chr4     75000   90000   chr4    75000   90000   0.05170981720312409
    chr4     75000   90000   chr4    345000  360000  0.01292745430078102
    chr4     90000   105000  chr4    90000   105000  0.01292745430078102
    .
    .
    .
    .
    .
    .


Usage:
    .. code-block:: bash
        
        usage: gcMapExplorer pairCoo2cmap [-h] [-i maps.txt] [-ccm RawObserved]
                                          [-od OUTDIR] [-gcm inOut.gcmap] [-cmeth lzf]
                                          [-dmeth sum]
                                          [-wd /home/rajendra/deskForWork/scratch]

                                          
**Optional arguments:**

.. code-block:: bash

  -h, --help            show this help message and exit
  -i maps.txt, --input maps.txt
                         Input file name.
                        
  -ccm RawObserved, --ccmap RawObserved
                         Use this to convert all contact maps to ccmaps file. Provide suffix of
                        ccmap file names with this option and it will enable the conversion.
                        
                        Output ccmap file name is generated automatically as follows;
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
                        a factor of two, consecutively. e.g.: In case of 10 kb input resolution,
                        downsampled maps of "20kb", "40kb", "80kb", "160kb", "320kb" etc. will be
                        generated until, map size is less than 500.
                        
  -wd /home/rajendra/deskForWork/scratch, --work-dir /home/rajendra/deskForWork/scratch
                        Directory where temporary files will be stored.

