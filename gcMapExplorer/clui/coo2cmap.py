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
"""Import COO sparse matrix format to ccmap or gcmap
=================================================

As shown below in example, in this format, first and second column is location
on chromosome and third column is the respective value:

20000000        20000000        2692.0
20000000        20100000        885.0
20100000        20100000        6493.0
20000000        20200000        15.0
20100000        20200000        52.0
20200000        20200000        2.0
20000000        20300000        18.0
20100000        20300000        40.0

NOTE that, above location is real value. However, with -idx/--index option,
these two same column willbe considered as index value. index should always
start from zero for absoulte beginning of chromosome.e.g. for 10kb, 0-10000
should have index of zero, 10000-20000 have index of one. If this is file
format,resolution should be provided with -r/--resolution option.

=================================================
"""

inputMetaFileHelp = \
"""Meta input file containing input contact map files list with respective
xlabel and ylabel. xlabel should be always provided. In case of intra-
chromosomal map, only xlabel is sufficient because both x and y axis are of
same chromosome. However for inter-chromosomal map, both xlabel and ylabel
should be provided. Example format:

100kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_100kb.RAWobserved    chr1
100kb_resolution_intrachromosomal/chr5/MAPQGE30/chr5_100kb.RAWobserved    chr5
100kb_resolution_intrachromosomal/chr15/MAPQGE30/chr15_100kb.RAWobserved  chr15
100kb_resolution_intrachromosomal/chr20/MAPQGE30/chr20_100kb.RAWobserved  chr20
100kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_100kb.RAWobserved  chr21
100kb_resolution_intrachromosomal/chr22/MAPQGE30/chr22_100kb.RAWobserved  chr22

"""

inputCompressedFileHelp = \
"""Input compressed archive file containing all the listed contact maps.
Presently, only "tar.gz" and "zip" compressed files are supported.

If -i/--input is not provided, all files from compressed file will be tried for
processing.

"""

mapTypeHelp = \
""" Type of listed contact maps: "intra" or "inter" chromosomal map.

"""

resolutionHelp = \
"""Resolution of all maps. It is an optional argument. Note that, if this
option is not provided, resolution will be automatically determined from the
contact map file. However, in case of -idx/--index option, resolution
should be provided as resolution cannot be determined from input contact map
file.

"""

indexHelp = \
"""It determines whether contact map files have real coordinate of chromosome
or index number. If this option is enabled, -r/--resolution option should be
provided.

"""

ccmapSuffixHelp = \
""" Use this to convert all contact maps to ccmap format files. Provide suffix
of ccmap file names with this option and it will enable the conversion.

Ouput ccmap file name is generated outmatically as follows;
if xlabel is not equal to ylabel: <xlabel>_<ylabel>_<suffix>.ccmap
else: <xlabel>_<suffix>.ccmap

Note that -od/--out-dir option is also required because all ccmaps will be
saved in this directory.

"""

outDirHelp = \
""" Directory to save all ccmap files. It should be provided when -ccma/--ccmap
option is used.

"""

fileGCMapHelp = \
"""Provide gcmap file to convert all contact maps into one gcmap file.
File name should contain full path because -od/--out-dir is not considered
for this conversion.

"""

downsampleMethodHelp = \
"""Downsampling method to coarsen the resolution in gcmap file. The option is
intended to use with -gcm/--gcmap option. Three accepted methods are
        sum  : sum of values,
        mean : Average of values and
        max  : Maximum of values.

This option generates all coarser maps where resolutions will be coarsened by
a factor of two, consequetively. e.g.: In case of 10 kb input resolution,
downsampled maps of "20kb", "40kb", "80kb", "160kb", "320kb" etc. will be
generated until, map size is less than 500.

"""

def main():
    # Construct command line arguments and parsed it
    parser, args = parseArguments()

    if args.inputMetaFile is None and args.inputCompressed is None:
        msg = "Neither input meta file nor input" \
                + " compressed file was provided.\n"
        showErrorAndExit(parser, msg)

    # args.inputcompressed is a IO object, so close it but retains file name
    inputCompressed = None
    if args.inputCompressed is not None:
        inputCompressed = args.inputCompressed.name
        args.inputCompressed.close()

    # Make a list of input contact files, xlabel and ylabel
    inputFiles, xlabels, ylabels = None, None, None
    if args.inputMetaFile is not None:
        inputFiles, xlabels, ylabels = makeInputFileList(
                                        args.inputMetaFile,
                                        parser,
                                        inputCompressed=inputCompressed)

    if args.ccmapSuffix is not None and args.outDir is None:
        msg = "No output directory is given for ccmap files!!!\n"
        showErrorAndExit(parser, msg)

    if args.ccmapSuffix is None and args.fileGCMap is None:
        showErrorAndExit(parser, "No output format directed!!!\n")

    if args.index:
        coordinate = 'index'
    else:
        coordinate = 'real'

    # Check for scratch directory
    if not os.path.isdir(args.workDir):
        showErrorAndExit(parser, '\nScratch Dirctory "{0}" not found !!!\n'.format(args.workDir))

    cooReader = gmlib.importer.CooMatrixHandler(inputFiles,
                                                inputCompressed,
                                                mapType=args.mapType,
                                                resolution=args.resolution,
                                                coordinate=coordinate,
                                                workDir=args.workDir)

    outputFileList = []
    if  args.ccmapSuffix is not None:
        for i in range(len(inputFiles)):
            if xlabels[i] == ylabels[i]:
                name = xlabels[i] + '_' + args.ccmapSuffix + '.ccmap'
            else:
                name = xlabels[i] + '_' + ylabels[i] + '_' \
                        + args.ccmapSuffix + '.ccmap'

            outputFileList.append(os.path.join(args.outDir, name))


    if  args.ccmapSuffix is not None and args.fileGCMap is not None:
        cooReader.save_ccmaps(outputFileList, xlabels=xlabels, ylabels=ylabels)
        del cooReader

        for inFile in outputFileList:
            ccmap = gmlib.ccmap.load_ccmap(infile=inFile, workDir=args.workDir)
            gmlib.gcmap.addCCMap2GCMap(ccmap, args.fileGCMap,
                                        compression=args.compression,
                                        coarsingMethod=args.coarsingMethod)
            del ccmap

    elif args.ccmapSuffix is not None:
        cooReader.save_ccmaps(outputFileList, xlabels=xlabels, ylabels=ylabels)
    elif args.fileGCMap is not None:
        cooReader.save_gcmap(args.fileGCMap, xlabels=xlabels, ylabels=ylabels,
                                coarsingMethod=args.coarsingMethod,
                                compression=args.compression)
    else:
        pass


def parseArguments():
    parser = argparse.ArgumentParser(
                prog='gcMapExplorer coo2cmap',
                description=description,
                formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('r'), metavar='input.txt',
                        dest='inputMetaFile', help=inputMetaFileHelp)

    parser.add_argument('-ic', '--input-compressed', action='store',
                        type=argparse.FileType('r'), metavar='input.tar.gz',
                        dest='inputCompressed', help=inputCompressedFileHelp)

    parser.add_argument('-mt', '--mapType', action='store', type=str,
                        metavar='intra', default='intra', dest='mapType',
                        help=mapTypeHelp, choices=['intra', 'inter'])

    parser.add_argument('-r', '--resolution', action='store', type=str,
                        metavar='10kb', dest='resolution', help=resolutionHelp)

    parser.add_argument('-idx', '--index', action='store_true', dest='index',
                        help=indexHelp)

    parser.add_argument('-ccm', '--ccmap', action='store', dest='ccmapSuffix',
                        metavar='10kb_RawObserved', help=ccmapSuffixHelp)

    parser.add_argument('-od', '--out-dir', action='store', dest='outDir',
                        help='Directory where all ccmap files will be saved.\n')

    parser.add_argument('-gcm', '--gcmap', action='store', dest='fileGCMap',
                        metavar='inOut.gcmap', help=fileGCMapHelp)

    parser.add_argument('-cmeth', '--compression-method', action='store',
                        dest='compression', metavar='lzf',
                        choices=['lzf', 'gzip'], default='lzf',
                        help='Data compression method in gcmap file.\n')

    parser.add_argument('-dmeth', '--downsample-method', action='store',
                        dest='coarsingMethod', metavar='sum',
                        choices=['max', 'mean', 'sum'], default='sum',
                        help=downsampleMethodHelp)


    parser.add_argument('-wd', '--work-dir', action='store', dest='workDir',
                        default=config['Dirs']['WorkingDirectory'],
                        metavar=config['Dirs']['WorkingDirectory'],
                        help='Directory where temporary files will be stored.')

    idx = sys.argv.index("coo2cmap")+1
    args = parser.parse_args(args=sys.argv[idx:])

    return parser, args


def makeInputFileList(fileObj, parser, inputCompressed=None):
    """ Make list of input contact map files
    """

    mapfiles, xlabels, ylabels = [], [], []
    for line in fileObj:

        line = line.lstrip().rstrip()  # Leading and trailing white spaces
        if not line.strip():           # Skip blank line
            continue

        # If files are not inside compressed file, check file exist
        temp = re.split('\s+', line)
        if inputCompressed is None:
            checkFileExist(temp[0], parser)

        # Check if xlabel is missing
        if len(temp) < 2:
            msg = "ERROR: {0} contains only map contact file for [{1}]. \
                    Need at least xlabel !!!\n".format(fileObj.name, temp[0])
            showErrorAndExit(parser, msg)

        # Append contact map file and xlabel
        mapfiles.append(temp[0])
        xlabels.append(temp[1])

        # Determine whether ylabel is present
        if len(temp) >= 3:
            ylabels.append(temp[2])
        else:
            ylabels.append(temp[1])

    return mapfiles, xlabels, ylabels


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
