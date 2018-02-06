#!/usr/bin/env python3
#
# Author: Rajendra Kumar
#
# This file is part of gcMapExplorer
# Copyright (C) 2016-2017  Rajendra Kumar, Ludvig Lizana, Per Stenberg
#
# gcMapExplorer is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gcMapExplorer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gcMapExplorer.  If not, see <http://www.gnu.org/licenses/>.
#
#=============================================================================

import os
import sys
import re
import argparse

import gcMapExplorer.lib as gmlib

from gcMapExplorer.config import getConfig

# get configuration
config = getConfig()

description = \
"""Import Bin-Contact format files
=================================================
In this format, two separate files are available. One file contains bins
information and other contains contact frequency.

These types of files are present in following GEO data:
* http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61471
* http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34453

This format contains a pair of file:
BIN file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CONTACT file in list format:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cbin1	cbin2	expected_count	observed_count
1	1	40.245201	21339
1	2	83.747499	5661
1	3	92.12501	1546
1	4	93.401273	864
1	5	87.265472	442
.
.
.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Both BIN and CONTACT files are necessary for the conversion.

=================================================
"""

inputBinFileHelp = \
""" Input BIN file as shown above.

"""

inputContactFileHelp = \
""" Input CONTACT file as shown above.

"""

ccmapSuffixHelp = \
""" Use this to convert all contact maps to ccmaps file. Provide suffix of
ccmap file names with this option and it will enable the conversion.

Output ccmap file name is generated automatically as follows;
<chromosome>_<resolution>_<suffix>.ccmap

Note that -od/--out-dir option is also required because all ccmaps will be
saved in this directory.

"""

outDirHelp = \
""" Directory to save all ccmap files. It should be provided when -ccmap/--ccmap
option is used.

"""

fileGCMapHelp = \
"""Provide gcmap file to convert all contact maps into one gcmap file.
File name should contain full path because -od/--out-dir is not considered
for thi conversion.

"""

downsampleMethodHelp = \
"""Downsampling method to coarsen the resolution in gcmap file. The option is
intended to use with -gcm/--gcmap option. Three accepted methods are
        sum  : sum of values,
        mean : Average of values and
        max  : Maximum of values.

This option generates all coarser maps where resolutions will be coarsened by
a factor of two, consecutively. e.g.: In case of 10 kb input resolution,
downsampled maps of "20kb", "40kb", "80kb", "160kb", "320kb" etc. will be
generated until, map size is less than 500.

"""

def main():

    # Construct command line arguments and parsed it
    parser, args = parseArguments()

    # Input Bin file
    inputBinFile = None
    if args.inputBinFile is not None:
        inputBinFile = args.inputBinFile.name
        args.inputBinFile.close()
    else:
        showErrorAndExit(parser, "No BIN Input File!!!\n")

    # Input contact file
    inputContactFile = None
    if args.inputContactFile is not None:
        inputContactFile = args.inputContactFile.name
        args.inputContactFile.close()
    else:
        showErrorAndExit(parser, "No CONTACT Input File!!!\n")

    if args.ccmapSuffix is not None and args.outDir is None:
        msg = "No output directory is given for ccmap files!!!\n"
        showErrorAndExit(parser, msg)

    if args.ccmapSuffix is None and args.fileGCMap is None:
        showErrorAndExit(parser, "No output format directed!!!\n")

    # Check for scratch directory
    if not os.path.isdir(args.workDir):
        showErrorAndExit(parser, '\nScratch Directory "{0}" not found !!!\n'.format(args.workDir))

    # Initialize
    binContactReader = gmlib.importer.BinsNContactFilesHandler(
                                        inputBinFile,
                                        inputContactFile,
                                        workDir=os.path.normpath(args.workDir))

    # Convert ccmaps
    if args.ccmapSuffix is not None:
        # Save ccmaps
        binContactReader.save_ccmaps(os.path.normpath(args.outDir), args.ccmapSuffix)

    # Convert gcmap
    if args.fileGCMap is not None:
        binContactReader.save_gcmap(os.path.normpath(args.fileGCMap),
                                    coarseningMethod=args.coarseningMethod,
                                    compression=args.compression)

def parseArguments():
    parser = argparse.ArgumentParser(
                prog='gcMapExplorer bc2cmap',
                description=description,
                formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-ib', '--input-bin', action='store',
                        type=argparse.FileType('r'),
                        metavar='nm_none_160000.bins',
                        dest='inputBinFile', help=inputBinFileHelp)

    parser.add_argument('-ic', '--input-contact', action='store',
                        type=argparse.FileType('r'),
                        metavar='nm_none_160000.n_contact',
                        dest='inputContactFile', help=inputContactFileHelp)

    parser.add_argument('-ccm', '--ccmap', action='store', dest='ccmapSuffix',
                        metavar='RawObserved', help=ccmapSuffixHelp)

    parser.add_argument('-od', '--out-dir', action='store', dest='outDir',
                        help='Directory where all ccmap files will be saved.\n')

    parser.add_argument('-gcm', '--gcmap', action='store', dest='fileGCMap',
                        metavar='inOut.gcmap', help=fileGCMapHelp)

    parser.add_argument('-cmeth', '--compression-method', action='store',
                        dest='compression', metavar='lzf',
                        choices=['lzf', 'gzip'], default='lzf',
                        help='Data compression method for gcmap file.\n')

    parser.add_argument('-dmeth', '--downsample-method', action='store',
                        dest='coarseningMethod', metavar='sum',
                        choices=['max', 'mean', 'sum'], default='sum',
                        help=downsampleMethodHelp)

    parser.add_argument('-wd', '--work-dir', action='store', dest='workDir',
                        default=config['Dirs']['WorkingDirectory'],
                        metavar=config['Dirs']['WorkingDirectory'],
                        help='Directory where temporary files will be stored.')

    idx = sys.argv.index("bc2cmap")+1
    args = parser.parse_args(args=sys.argv[idx:])

    return parser, args


def checkFileExist(filename, parser):
    if not os.path.isfile(filename):
        parser.print_help()
        print("\n===== ERROR =======")
        print("File: {0} does not exist !!!\n".format(filename))
        print("See Usage Above!!!")
        sys.exit(-1)


def showErrorAndExit(parser, message):
    parser.print_help()
    print("\n===== ERROR =======")
    print(message)
    print("See Usage Above!!!")
    sys.exit(-1)


if __name__ == "__main__":
    main()
