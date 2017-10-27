normKR - Normalization by Vanilla-Coverage
------------------------------------------

This method was first used for inter-chromosomal map. Later it was used for
intra-chromosomal map by Rao et al., 2014 (http://dx.doi.org/10.1016/j.cell.2014.11.021).
This is a simple method where at first each element is divided by sum of
respective row and subsequently divided by sum of respective column.

For more details, see publications:
  1) https://doi.org/10.1126/science.1181369
  2) http://dx.doi.org/10.1016/j.cell.2014.11.021

Usage:

  .. code-block:: bash

      usage: gcMapExplorer normVC [-h] [-i input.gcmap] [-fi gcmap]
                                  [-o output.gcmap] [-fo gcmap] [-sq] [-vmax VMAX]
                                  [-vmin VMIN] [-ptnd 99] [-tdo 0.8] [-cmeth lzf]
                                  [-wd /home/*****/scratch]

**Optional arguments:**

.. code-block:: bash

    -h, --help            show this help message and exit
    -i input.gcmap, --input input.gcmap
                         Input ccmap or gcmap file.

    -fi gcmap, --format-input gcmap
                         Input format: 'ccmap' or 'gcmap'.

    -o output.gcmap, --output output.gcmap
                         Output ccmap or gcmap file.

                        When input file is ccmap, output file can be gcmap. However, when a input file
                        is gcmap, output file will be only in gcmap.

    -fo gcmap, --format-output gcmap
                         Output format: 'ccmap' or 'gcmap'.

                        When input file is ccmap, output file can be gcmap. However, when a input file
                        is gcmap, output file will be only in gcmap.

    -sq, --sqroot          Square-root of normalized map

    -vmax VMAX, --maximum-value VMAX
                         Minimum threshold value for normalization.
                        If contact frequency is less than or equal to this threshold value,
                        this value is discarded during normalization.

    -vmin VMIN, --minimum-value VMIN
                         Maximum threshold value for normalization.
                        If contact frequency is greater than or equal to this threshold value,
                        this value is discarded during normalization.

    -ptnd 99, --percentile-threshold-no-data 99
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

    -tdo 0.8, --threshold-data-occupancy 0.8
                         It can be used to filter the map, where rows/columns with largest numbers
                        of missing data can be discarded.This ratio is:
                          (number of bins with data) / (total number of bins in the given row/column)

                        For example: if -tdo = 0.8, then all rows containing more than 20% of
                        missing data will be discarded.

    -cmeth lzf, --compression-method lzf
                        Data compression method for output gcmap file.
    -wd /home/*****/scratch, --work-dir /home/*****/scratch
                        Directory where temporary files will be stored.
