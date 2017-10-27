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
"""To calculate overall/block-wise correlation between contact maps
===================================================================

This tool can be used to calculate correlation between two contact maps.

With option "-bs/--block-size", it calculates correlation along the blocks of
input size. Without this option, it calculates correlation between entire maps.

NOTE: Both input files should be in SAME format.
===================================================================

"""

inputOneHelp = \
""" First input file in either ccmap or gcmap format

"""

inputTwoHelp = \
""" Second input file in either ccmap or gcmap format

"""

blockSizeHelp = \
""" Size of the sliding block in words (kb, mb). If this option is not provided
Correlation between whole map will be computed.

"""

slideStepSizeHelp = \
""" Step size to slide the block along diagonal.

"""

nameHelp = \
""" Name of output data
When results are written in h5 file, a name is necessary to label the data.
The correlation results will be stored with this name in h5 file given with
"-o/--output" option.

"""

outFileHelp = \
""" Name of output h5 file. The obtained correlation will be written in this
output h5 file.

"""


ignore_triangular_help = \
"""To NOT ignore upper triangle of the map.
Intra-chromosomal contact maps are symmetric, therefore, during correlation
calculation, by default, only one side of the diagonal is considered.
To consider the whole map use this option.

"""

cutoffPercentileHelp = \
"""To ignore bins that are below this percentile.
The maps contains large numbers of low contact frequency. This frequencies
might not be significant and inclusion in correlation calculation may affect the
correlation value. Therefore, this option can be used discard these lower
frequencies using percentile cutoff.

"""

diagonal_offset_help = \
""" To ignore diagonal including adjacent diagonal up to the given offset

"""

corrTypeHelp = \
""" Method of correlation: Pearson or Spearman.

"""


def main():
    # Construct command line arguments and parsed it
    parser, args = parseArguments()

    # Check if no input files are provided
    if args.inputOne is None or args.inputTwo is None:
        showErrorAndExit(parser, "No input files.\n")

    # First input file
    inputOne = None
    inputOne = args.inputOne.name
    args.inputOne.close()

    # Second input file
    inputTwo = None
    inputTwo = args.inputTwo.name
    args.inputTwo.close()

    # Compare extension of input files to check file format
    extOne = os.path.splitext(os.path.basename(inputOne))[-1]
    extTwo = os.path.splitext(os.path.basename(inputTwo))[-1]
    if extOne != extTwo:
        msg = 'Input file format appears to be different. Use input files with same format.\n'
        showErrorAndExit(parser, msg)

    # gcmap or ccmap
    gcmap = False
    if extOne == '.gcmap':
        gcmap = True

    # Check for block-size
    blockSize = None
    if args.blockSize is not None:
        try:
            gmlib.util.resolutionToBinsize(args.blockSize)
        except ValueError:
            msg = '{0} is not an accepted resolution unit!!\n'.format(args.blockSize)
            showErrorAndExit(parser, msg)

        blockSize = args.blockSize

    # Check for data label
    name = None
    if args.blockSize is not None and args.name is None:
        showErrorAndExit(parser, 'No name is given to label the data!!!\n')
    else:
        name = args.name

    # Check if outfile is given
    if args.outFile is None and args.blockSize is not None:
        showErrorAndExit(parser, 'No output file is given!!!\n')

    # If out file given, store name
    outFile = None
    if args.outFile is not None:
        outFile = args.outFile.name
        args.outFile.close()

    # Check for cutoff percentile option
    cutoffPercentile = None
    if args.cutoffPercentile is not None:
        if args.cutoffPercentile <0 or args.cutoffPercentile > 100:
            msg = 'Percentile cutoff should be in range of 0 to 100.\n'
            showErrorAndExit(parser, msg)
        else:
            cutoffPercentile = args.cutoffPercentile

    if gcmap:
        results = gmlib.cmstats.correlateGCMaps(inputOne, inputTwo,
                                        outFile=outFile,
                                        blockSize=blockSize, name=name,
                                        slideStepSize=args.slideStepSize,
                                        cutoffPercentile=cutoffPercentile,
                                        ignore_triangular=args.ignore_triangular,
                                        diagonal_offset=args.diagonal_offset,
                                        corrType=args.corrType,
                                        workDir=args.workDir)

        if blockSize is None and results is not None:
            printResultMany(results[0], results[1], results[2])

    else:
        ccMapObjOne = gmlib.ccmap.load_ccmap(inputOne, workDir=workDir)
        ccMapObjOne.make_readable()

        ccMapObjTwo = gmlib.ccmap.load_ccmap(inputTwo, workDir=workDir)
        ccMapObjTwo.make_readable()

        dataArray = np.zeros(shape=(max(ccMapObjOne.shape[0], ccMapObjTwo.shape[0]),) )

        if blockSize is not None:
            h5Handle = gdh.HDF5Handler(outFile)

        try:
        	results = gmlib.cmstats.correlateCMaps(ccMapObjOne, ccMapObjTwo,
                                ignore_triangular=args.ignore_triangular,
                                diagonal_offset=args.diagonal_offset,
                                corrType=args.corrType,
                                blockSize=blockSize,
                                cutoffPercentile=cutoffPercentile,
                                slideStepSize=args.slideStepSize,
                                workDir=workDir, outFile=None)

        except	(KeyboardInterrupt, SystemExit) as e:
        	del ccMapObjOne
        	del ccMapObjTwo
        	raise e

        corr, centers = None, None
        if results is not None:
        	corr, centers = results

        if blockSize is not None:
        	centers = np.asarray(centers) / binsize
        	centers = centers.astype(int)
        	dataArray[(centers,)] = corr
        	h5Handle._addDataByArray(key, cmp.binsizeToResolution(binsize), name, dataArray)
        else:
            printResultOne(corr, centers)


def printResultOne(corr, pvalue):
    line = "\n ===================== \n"
    line += " \n Correlation = {0} \n".format(corr)
    line += " \n P-Value = {1} \n".format(pvalue)
    line += "\n ===================== \n"
    sys.stdout.write(line)

def printResultMany(mapList, corrs, pvalues):
    sys.stdout.write("\n ===================== \n")
    sys.stdout.write("\n # Name \t Correlation \t P-Value \n")
    for i in range(len(mapList)):
        sys.stdout.write("\n{0}\t{1}\t{2}\n".format(mapList[i], corrs[i], pvalues[i]))
    sys.stdout.write("\n ===================== \n")


def parseArguments():
    parser = argparse.ArgumentParser(
                prog='gcMapExplorer corrBWcmaps',
                description=description,
                formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i1', '--input-one', action='store',
                        type=argparse.FileType('r'), metavar='input1.gcmap',
                        dest='inputOne', help=inputOneHelp)

    parser.add_argument('-i2', '--input-two', action='store',
                        type=argparse.FileType('r'), metavar='input2.gcmap',
                        dest='inputTwo', help=inputTwoHelp)

    parser.add_argument('-bs', '--block-size', action='store', type=str,
                        metavar='1mb', dest='blockSize',
                        help=blockSizeHelp)

    parser.add_argument('-sss', '--slide-step-size', action='store', type=int,
                        metavar=1, dest='slideStepSize', default=1,
                        help=slideStepSizeHelp)

    parser.add_argument('-o', '--output', action='store',
                        type=argparse.FileType('a'), metavar='output.h5',
                        dest='outFile', help=outFileHelp)

    parser.add_argument('-n', '--name', action='store', type=str,
                        metavar='corr', dest='name', help=nameHelp)

    parser.add_argument('-nit', '--no-ignore-triangular', action='store_false',
                        dest='ignore_triangular', help=ignore_triangular_help)

    parser.add_argument('-do', '--diagonal-offset', action='store', type=int,
                        metavar='1', dest='diagonal_offset', default=1,
                        help=diagonal_offset_help)

    parser.add_argument('-rmeth', '--correlation-method', action='store',
                        dest='corrType', metavar='spearman',
                        choices=['pearson', 'spearman'], default='spearman',
                        help=corrTypeHelp)

    parser.add_argument('-ct', '--cutoff-percentile', action='store',
                        dest='cutoffPercentile', metavar='95', type=float,
                        help=cutoffPercentileHelp)

    parser.add_argument('-wd', '--work-dir', action='store', dest='workDir',
                        default=config['Dirs']['WorkingDirectory'],
                        metavar=config['Dirs']['WorkingDirectory'],
                        help='Directory where temporary files will be stored.')

    idx = sys.argv.index("corrBWcmaps")+1
    args = parser.parse_args(args=sys.argv[idx:])

    return parser, args


def showErrorAndExit(parser, message):
    parser.print_help()
    print("\n===== ERROR =======")
    print(message)
    print("See Usage Above!!!")
    sys.exit(-1)

if __name__ == "__main__":
    main()
