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
"""Scale maps using Median/Mean Contact Frequency
=================================================

This method can be used to normalize contact map with expected values.
These expected values could be either Median or Average contact values
for particular distance between two locations/coordinates. At first,
Median/Average distance contact frequency for each distance is calculated.
Subsequently, the observed contact frequency is either divided ('o/e') or
substracted ('o-e') by median/average contact frequency obtained for
distance between the two locations.

=================================================

"""

inputFileHelp = \
""" Input ccmap or gcmap file.

"""

inputFileFormatHelp = \
""" Input format: 'ccmap' or 'gcmap'.

"""


outputFileHelp = \
""" Output ccmap or gcmap file.

When input file is ccmap, ouput file can be gcmap. However, when a input file
is gcmap, output file will be only in gcmap.

"""

outputFileFormatHelp = \
""" Input format: 'ccmap' or 'gcmap'.

When input file is ccmap, ouput file can be gcmap. However, when a input file
is gcmap, output file will be only in gcmap.

"""

statsHelp = \
""" Statistics to be considered for scaling.
It may be either “mean” or “median”. By default, it is “median”.

"""

vminHelp = \
""" Minimum thershold value for normalization.
If contact frequency is less than or equal to this thershold value,
this value is discarded during normalization.

"""

vmaxHelp = \
""" Maximum thershold value for normalization.
If contact frequency is greater than or equal to this thershold value,
this value is discarded during normalization.

"""

stypeHelp = \
""" Type of scaling.
It may be either 'o/e' or 'o-e'. In case of 'o/e',
Observed/Expected will be calculated while (Observed - Expected)
will be calculated for 'o-e'.

"""

percentile_thershold_no_data_help = \
""" It can be used to filter the map, where rows/columns with largest numbers
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

"""

thershold_data_occup_help = \
""" It can be used to filter the map, where rows/columns with largest numbers
of missing data can be discarded.This ratio is:
  (number of bins with data) / (total number of bins in the given row/column)

For example: if -tdo = 0.8, then all rows containing more than 20%% of
missing data will be discarded.

"""

def main():

    # Construct command line arguments and parsed it
    parser, args = parseArguments()

    # check input file exist
    if args.inputFile is not None:
        checkFileExist(args.inputFile, parser)
    else:
        msg = 'No input file!!!'
        showErrorAndExit(parser, msg)

    # check input file format provided
    if args.inputFileFormat is None:
        msg = 'No input file format given !!!'
        showErrorAndExit(parser, msg)

    # check output file provided
    if args.outputFile is None:
        msg = 'No output file name given !!!'
        showErrorAndExit(parser, msg)

    # check output file format provided
    if args.outputFileFormat is None:
        msg = 'No output file format given !!!'
        showErrorAndExit(parser, msg)

    if args.outputFileFormat == 'ccmap' and args.inputFileFormat == 'gcmap':
        msg = 'With gcmap input format, output format should be gcmap format.'
        showErrorAndExit(parser, msg)

    if args.percentile_thershold_no_data is not None and args.thershold_data_occup is not None:
        msg = 'Both Percentile and Fraction filters cannot be used simultaneously!!'
        showErrorAndExit(parser, msg)

    if args.outputFileFormat == 'ccmap' and args.inputFileFormat == 'ccmap':
        gmlib.normalizer.normalizeCCMapByMCFS(args.inputFile, stats=args.stats, outFile=args.outputFile,
                                            vmin=args.vmin, vmax=args.vmax, stype=args.stype,
                                            percentile_thershold_no_data=args.percentile_thershold_no_data,
                                            thershold_data_occup=args.thershold_data_occup,
                                            workDir=args.workDir)

    if args.outputFileFormat == 'gcmap' and args.inputFileFormat == 'ccmap':
        norm_ccmap = gmlib.normalizer.normalizeCCMapByMCFS(args.inputFile, stats=args.stats,
                                            vmin=args.vmin, vmax=args.vmax, stype=args.stype,
                                            percentile_thershold_no_data=args.percentile_thershold_no_data,
                                            thershold_data_occup=args.thershold_data_occup,
                                            workDir=args.workDir)
        gmlib.gcmap.addCCMap2GCMap(norm_ccmap, args.outputFile, compression=args.compression)


    if args.outputFileFormat == 'gcmap' and args.inputFileFormat == 'gcmap':
        gmlib.normalizer.normalizeGCMapByMCFS(args.inputFile, args.outputFile,
                                            stats=args.stats, vmin=args.vmin, vmax=args.vmax, stype=args.stype,
                                            percentile_thershold_no_data=args.percentile_thershold_no_data,
                                            thershold_data_occup=args.thershold_data_occup,
                                            compression=args.compression,
                                            workDir=args.workDir)

def parseArguments():
    parser = argparse.ArgumentParser(
                prog='gcMapExplorer normMCFS',
                description=description,
                formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        metavar='input.gcmap',
                        dest='inputFile', help=inputFileHelp)

    parser.add_argument('-fi', '--format-input', action='store',
                        metavar='gcmap',
                        choices=['gcmap', 'ccmap', 'hicmap'],
                        dest='inputFileFormat', help=inputFileFormatHelp)

    parser.add_argument('-o', '--output', action='store',
                        metavar='output.gcmap',
                        dest='outputFile', help=outputFileHelp)

    parser.add_argument('-fo', '--format-output', action='store',
                        metavar='gcmap',
                        choices=['gcmap', 'ccmap'],
                        dest='outputFileFormat', help=outputFileFormatHelp)

    parser.add_argument('-s', '--stats', action='store',
                        dest='stats', metavar='median',
                        choices=['median', 'mean'], default='median',
                        help=statsHelp)
    parser.add_argument('-vmax', '--maximum-value', action='store',
                        dest='vmax', type=float, metavar=None,
                        help=vminHelp)

    parser.add_argument('-vmin', '--minimum-value', action='store',
                        dest='vmin', metavar=None, type=float,
                        help=vmaxHelp)

    parser.add_argument('-st', '--stype', action='store',
                        dest='stype', metavar='o/e', type=str,
                        choices=['o/e', 'o-e'], default='o/e',
                        help=stypeHelp)

    parser.add_argument('-ptnd', '--percentile-thershold-no-data', action='store',
                        dest='percentile_thershold_no_data', metavar=99,
                        type = float,
                        help=percentile_thershold_no_data_help)

    parser.add_argument('-tdo', '--thershold-data-occupancy', action='store',
                        dest='thershold_data_occup', metavar=0.8,
                        type = float,
                        help=thershold_data_occup_help)

    parser.add_argument('-cmeth', '--compression-method', action='store',
                        dest='compression', metavar='lzf',
                        choices=['lzf', 'gzip'], default='lzf',
                        help='Data compression method for output gcmap file.\n')

    parser.add_argument('-wd', '--work-dir', action='store', dest='workDir',
                        default=config['Dirs']['WorkingDirectory'],
                        metavar=config['Dirs']['WorkingDirectory'],
                        help='Directory where temporary files will be stored.')

    idx = sys.argv.index("normMCFS")+1
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
