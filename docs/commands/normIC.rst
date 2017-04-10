normIC - Normalization by Iterative Correction
----------------------------------------------

This method normalize the raw contact map by removing biases from
experimental procedure. This is an method of matrix balancing, however, in the
normalized, sum of rows and columns are **not** equal to one.

Usage:

  .. code-block:: bash

    usage: gcMapExplorer normIC [-h] [-i input.gcmap] [-fi gcmap]
                                [-o output.gcmap] [-fo gcmap] [-t 0.0001]
                                [-vmax VMAX] [-vmin VMIN] [-c 500] [-ptnd 99]
                                [-tdo 0.8] [-cmeth lzf]
                                [-wd /home/rajendra/deskForWork/scratch]



**Optional arguments:**

.. code-block:: bash

  -h, --help            show this help message and exit
  -i input.gcmap, --input input.gcmap
                         Input ccmap or gcmap file.

  -fi gcmap, --format-input gcmap
                         Input format: 'ccmap' or 'gcmap'.

  -o output.gcmap, --output output.gcmap
                         Output ccmap or gcmap file.

                        When input file is ccmap, ouput file can be gcmap. However, when a input file
                        is gcmap, output file will be only in gcmap.

  -fo gcmap, --format-output gcmap
                         Input format: 'ccmap' or 'gcmap'.

                        When input file is ccmap, ouput file can be gcmap. However, when a input file
                        is gcmap, output file will be only in gcmap.

  -t 0.0001, --tolerance 0.0001
                         Tolerance for matrix balancing.
                        Smaller tolreance increases accuracy in sums of rows and columns.

  -vmax VMAX, --maximum-value VMAX
                         Minimum thershold value for normalization.
                        If contact frequency is less than or equal to this thershold value,
                        this value is discarded during normalization.

  -vmin VMIN, --minimum-value VMIN
                         Maximum thershold value for normalization.
                        If contact frequency is greater than or equal to this thershold value,
                        this value is discarded during normalization.

  -c 500, --iteration 500
                         Number of iteration to stop the normalization.

  -ptnd 99, --percentile-thershold-no-data 99
                         It can be used to filter the map, where rows/columns with largest numbers
                        of missing data can be discarded. Its value should be between 1 and 100.
                        This options discard the rows and columns which are above this percentile.
                        For example: if this value is 99, those rows or columns will be discarded which
                        contains larger than number of zeros (missing data) at 99 percentile.

                        To calculate percentile, all blank rows are removed, then in all rows, number
                        of zeros are counted. Afterwards, number of zeros at input percentile is
                        obtained. In next step, if a row contain number of zeros larger than this
                        percentile value, the whole row and column is assigned to have missing data.
                        This percentile indicates highest numbers of zeros (missing data) in given
                        rows/columns.

  -tdo 0.8, --thershold-data-occupancy 0.8
                         It can be used to filter the map, where rows/columns with largest numbers
                        of missing data can be discarded.This ratio is:
                          (number of bins with data) / (total number of bins in the given row/column)

                        For example: if -tdo = 0.8, then all rows containing more than 20% of
                        missing data will be discarded.

  -cmeth lzf, --compression-method lzf
                        Data compression method for output gcmap file.
  -wd /home/rajendra/deskForWork/scratch, --work-dir /home/rajendra/deskForWork/scratch
                        Directory where temporary files will be stored.
