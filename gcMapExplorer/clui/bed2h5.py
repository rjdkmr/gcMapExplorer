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
"""Import a bed file to HDF5 format h5 file
============================================

bed file can be converted into gcMapExplorer compatible HDF5 file using
this tool. This HDF5 file can be loaded into gcMapExplorer browser for
interactive visualization.

This tool does not require any external program.

Resolutions
===========
By default, original data are downsampled to following resolutions:  '1kb',
'2kb', '4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb', '100kb', '160kb',
'200kb', '320kb', '500kb', '640kb',  and '1mb'.

The data are downsampled at this stage only to speed up the visualization
process as downsampling might slow down the interactive visualization.

Downsampling/Coarsening method
==============================
Presently, six methods are implemented:
1) min    -> Minimum value
2) max    -> Maximum value
3) amean  -> Arithmatic mean or average
4) hmean  -> Harmonic mean
5) gmean  -> Geometric mean
6) median -> Median

All these methods are used by default.
See below help for "-dm/--downsample-method" option to change the methods.

To keep original 1 base resolution data
=======================================
By default, the output h5 file does not contain original 1-base resolution
data to reduce the file size. To keep the original data in h5 file, used
-ko/--keep-original flag.

"""

inputFileHelp = \
"""Input wig file.

"""
dataColumnHelp = \
"""The column number, which is considered as data column. Column number
could vary and depends on BED format. For example:
1) ENCODE broadPeak format (BED 6+3): 7th column
2) ENCODE gappedPeak format (BED 12+3): 13th column
3) ENCODE narrowPeak format (BED 6+4): 7th column
4) ENCODE RNA elements format (BED 6+3): 7th column

"""

resolutionHelp = \
"""Additional input resolutions other than these resolutions: 1kb', '2kb',
'4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb', '100kb', '160kb','200kb',
'320kb', '500kb', '640kb',  and '1mb'.

Resolutions should be provided in comma seprated values. For Example:
-r "25kb, 50kb, 75kb"

"""

coarseningMethodHelp = \
"""Methods to coarse or downsample the data for converting from 1-base
to coarser resolutions. If this option is not provided, all six methods (see
above) will be considered. User may use only subset of these methods.
For example: -dm "max, amean" can be used for downsampling by only these
two methods.

"""

inputChromosomeHelp = \
"""Input Chromosome Name.
If this is provided, only this chromosome data is extracted and stored in h5
file.

"""

outHelp = \
"""Output h5 file.

If file is already present, it will replace the data. Therefore, be careful
if a file with same name is present.

"""

overwriteHelp = \
"""If a output file is already present, overwrite the datasets in the output
file.

"""

keepOriginalHelp = \
"""To copy original 1-base resolution data in h5 file. This will increase the
file size significantly.

"""

indexFileHelp = \
"""Index file in json format.
A file in json format containing indices (position in bed file) and sizes of
chromosomes. If this file is not present and given as input, a new file will be
generated. If this file is present, indices andsizes will be taken from this
file. If index and size of input chromosome is not present in json file, these
will be determined from bed file and stored in same json file. This file could
be very helpful in case when same bed file has to be read many times because
step to determine index and size of chromosome is skipped.

"""

def main():
    # Construct command line arguments and parsed it
    parser, args = parseArguments()

    # Check for input bigWig File
    if args.inputBedFile is None:
        showErrorAndExit(parser, '\nInput BED file is not given!!!\n')
    inputBedFile = args.inputBedFile.name
    args.inputBedFile.close()

    # Check for additional resolutions
    resolutions = args.resolutions
    if resolutions is not None:
        resolutions = get_resolution_list(resolutions, parser)

    # Check for input coarsening methods
    coarsening_methods = args.coarsening_methods
    if coarsening_methods is not None:
        coarsening_methods = get_coarsening_methods_list( coarsening_methods, parser )

    # Check for output file
    outFile = args.h5Name
    if outFile is not None:
        check_overwrite_status(outFile, args.overwrite, parser)
    else:
        showErrorAndExit(parser, '\nNo output file!!!\n')

    # Check for scratch directory
    if not os.path.isdir(args.workDir):
        showErrorAndExit(parser, '\nScratch Dirctory "{0}" not found !!!\n'.format(args.workDir))

    # Main conversion start here
    bed = gmlib.genomicsDataHandler.BEDHandler(inputBedFile,
                             indexFile=args.indexFile, column=args.column,
                            chromName=args.chromName, workDir=args.workDir)

    bed.saveAsH5(outFile, title=args.title, compression=args.compression,
                resolutions=resolutions, coarsening_methods=coarsening_methods,
                keep_original=args.keep_original)
    del bed

def get_resolution_list(resolution, parser):
    rlist = []
    temp = re.split(',', resolution)
    for r in temp:
        r = r.rstrip().lstrip()
        if not r.strip():
            continue
        try:
            b = gmlib.util.resolutionToBinsize(r)
            rlist.append( gmlib.util.binsizeToResolution(b) )
        except ValueError:
            showErrorAndExit(parser, '\n"{0}" contains resolution "{1}", which is not an acceptable resolution.\n'.format(resolution, r))

    return rlist

def get_coarsening_methods_list(coarsening_methods, parser):
    cmlist = []
    temp = re.split(',', coarsening_methods)
    for cm in temp:
        cm = cm.rstrip().lstrip()
        cmlist.append(cm)

    # Kept here to check if user has given correct input
    dummy = gmlib.genomicsDataHandler.check_coarsening_method(cmlist)

    return cmlist

def check_overwrite_status(outFile, overwrite, parser):

    # Nothing to do if overwrite
    if overwrite:
        return

    # Else do check other stuffs
    if os.path.isfile(outFile):
        print('\nWARNING: Output File "{0}" already exist.!!!'.format(outFile))

        s = input('\n Are you sure to overwrite datasets [y/n] (default: n): ')
        while(True):
            if not s:
                s = 'n'
                break

            s = s.lstrip().rstrip()

            if s != 'n' and s != 'y':
                s = input('\n Are you sure to overwrite datasets [y/n] (default: n): ')
            else:
                break

        if str(s) == 'n':
            showErrorAndExit(parser, '\nPlease provide another output file name.\n')

def parseArguments():
    parser = argparse.ArgumentParser(
                prog='gcMapExplorer bed2h5',
                description=description,
                formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('rb'), metavar='input.bed',
                        dest='inputBedFile', help=inputFileHelp)
    parser.add_argument('-t', '--title', action='store',
                        dest='title', metavar='"Genomic Dataset"',
                        default='"Genomic Dataset"',
                        help='Title of the dataset.\n')
    parser.add_argument('-dtc', '--data-column', action='store', type=int,
                        metavar=7, dest='column', default=7,
                        help=dataColumnHelp)
    parser.add_argument('-r', '--resolutions', action='store', type=str,
                        metavar='"List of Resolutions"', dest='resolutions',
                        help=resolutionHelp)
    parser.add_argument('-dm', '--downsample-method', action='store', type=str,
                        metavar='"List of downsampling method"',
                        dest='coarsening_methods',
                        help=coarseningMethodHelp)
    parser.add_argument('-icn', '--input-chromosome', action='store',
                        dest='chromName',
                        help=inputChromosomeHelp)
    parser.add_argument('-cmeth', '--compression-method', action='store',
                        dest='compression', metavar='lzf',
                        choices=['lzf', 'gzip'], default='lzf',
                        help='Data compression method in h5 file.\n')

    parser.add_argument('-o', '--out', action='store', dest='h5Name',
                        metavar='out.h5', help=outHelp)

    parser.add_argument('-ow', '--overwrite', action='store_true',
                        default = False,
                        dest='overwrite', help=overwriteHelp)

    parser.add_argument('-ko', '--keep-original', action='store_true',
                        default = False,
                        dest='keep_original', help=keepOriginalHelp)

    parser.add_argument('-idf', '--index-file', action='store',
                        metavar='index.json', dest='indexFile',
                        help=indexFileHelp)

    parser.add_argument('-wd', '--work-dir', action='store', dest='workDir',
                        default=config['Dirs']['WorkingDirectory'],
                        metavar=config['Dirs']['WorkingDirectory'],
                        help='Directory where temporary files will be stored.')

    idx = sys.argv.index("bed2h5")+1
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
