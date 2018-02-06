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
"""Import map from files similar to paired sparse matrix Coordinate (COO) format
===============================================================================

This format is very similar to COO format with additional information of
chromosome. Therefore, maps for all chromosome could be contained in a single
file.

This type of format appeared with following publication:
            http://dx.doi.org/10.1016/j.cell.2015.10.026
       GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72512

Following file format can be read as a text file, where first and second
column is location on chromosome and third column is the value:

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

===============================================================================
"""

inputFileHelp = \
""" Input file name.

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

    # Input File
    inputFile = None
    if args.inputFile is not None:
        inputFile = args.inputFile.name
        args.inputFile.close()
    else:
        showErrorAndExit(parser, "No Input File!!!\n")

    if args.ccmapSuffix is not None and args.outDir is None:
        msg = "No output directory is given for ccmap files!!!\n"
        showErrorAndExit(parser, msg)

    if args.ccmapSuffix is None and args.fileGCMap is None:
        showErrorAndExit(parser, "No output format directed!!!\n")

    # Check for scratch directory
    if not os.path.isdir(args.workDir):
        showErrorAndExit(parser, '\nScratch Directory "{0}" not found !!!\n'.format(args.workDir))

    if args.outDir is not None:
        args.outDir = os.path.normpath(args.outDir)

    if args.fileGCMap is not None:
        args.fileGCMap = os.path.normpath(args.fileGCMap)

    # Main region here
    pair_coo_reader = gmlib.importer.PairCooMatrixHandler(
                                    inputFile,
                                    ccmapOutDir = args.outDir,
                                    ccmapSuffix=args.ccmapSuffix,
                                    gcmapOut=args.fileGCMap)

    pair_coo_reader.setGCMapOptions(
            coarseningMethod=args.coarseningMethod,
            compression=args.compression
        )

    pair_coo_reader.runConversion()

def parseArguments():
    parser = argparse.ArgumentParser(
                prog='gcMapExplorer pairCoo2cmap',
                description=description,
                formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('r'), metavar='maps.txt',
                        dest='inputFile', help=inputFileHelp)

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

    idx = sys.argv.index("pairCoo2cmap")+1
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
