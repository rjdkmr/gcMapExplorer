hic2gcmap - convert hic file to gcmap
-------------------------------------

Hic files as described in [Juicer]_ can be converted to gcmap.

Usage:
    .. code-block:: bash

        usage: gcMapExplorer hic2gcmap [-h] (-a | -c A B | -l) [-cmp {lzf,gzip}]
                                       [-r R] [-n {VC,VC_SQRT,KR,none}]
                                       [-co {sum,mean,max,none}]
                                       input_file [output]

*Help:*

    .. code-block:: bash

        positional arguments:
          input_file            hic input file
          output                output file or directory (default: .)

        optional arguments:
          -h, --help            show this help message and exit
          -a, --all             all chromosome pairs (default: False)
          -c A B, --chromosomes A B
                                a pair of chromosomes A B (default: None)
          -l, --list            list all available chromosomes (default: False)
          -cmp {lzf,gzip}, --compression {lzf,gzip}
                                compression type (default: lzf)
          -r R, --resolution R  the resolution R, as an integer or as kb (default:
                                finest)
          -n {VC,VC_SQRT,KR,none}, --norm {VC,VC_SQRT,KR,none}
                                the type of norm to use (default: none)
          -co {sum,mean,max,none}, --coarsening {sum,mean,max,none}
                                the coarsening method (default: sum)


.. [Juicer] Juicer Provides a One-Click System for Analyzing Loop-Resolution Hi-C Experiments, Durand, Neva C. et al. Cell Systems, Volume 3, Issue 1, p. 95-98.

