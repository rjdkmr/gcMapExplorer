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
"""Download and Convert ENCODE datasets to h5 files
=================================================

It can be used to download and convert multiple datasets from ENCODE Experiment
matrix (https://www.encodeproject.org/matrix/?type=Experiment).
Presently, only bigWig files are downloaded and then converted.

At first search the datasets on https://www.encodeproject.org/matrix/?type=Experiment .
Then click on download button on top of the page. A text file will be downloaded.
This text file can be used as input in this program. All bigWig files will be
downloaded and converted to gcMapExplorer compatible hdf5 format.

NOTE: At first a metafile is automatically downloaded and then files
      are filtered according to bigWig format and Assembly. Subsequently,
      if several replicates are present, only datasets with combined
      replicates are considered. In case if two replicates are present
      and combined replicates are not present, replicates will be combined with
      '-mtc/--method-to-combine' option.

NOTE: Because downloading and conversion might take very long time, it also
      generates a checkpoint file in the output directory. Therefore,
      in case of crash or abrupt exit, the process can be continued from the
      last file.

Name of output files:

    (1) For ChIP-seq assay:
        a. signal-<Experiment target>-<Experiment accession>-<File accessions.h
        b. fold-<Experiment target>-<Experiment accession>-<File accessions>.h5

    (2) For RNA-seq:
        a. uniq-reads-<date>-<Experiment accession>-<File accessions>.h
        b. plus-uniq-reads-<date>-<Experiment accession>-<File accessions>.h
        c. minus-uniq-reads-<date>-<Experiment accession>-<File accessions>.h
        d. all-reads-<date>-<Experiment accession>-<File accessions>.h5
        e. plus-all-reads-<date>-<Experiment accession>-<File accessions>.h5
        f. minus-all-reads-<date>-<Experiment accession>-<File accessions>.h5
        g. signal-<date>-<Experiment accession>-<File accession>.h5

    (3) For DNase-seq:
        a. uniq-reads-signal-<date>-<Experiment accession>-<File accessions>.h
        b. raw-signal-<date>-<Experiment accession>-<File accessions>.h
        c. all-reads-signal-<date>-<Experiment accession>-<File accessions>.h
        d. signal-<date>-<Experiment accession>-<File accessions>.h5

    (4) For siRNA + RNA-seq:
        a. uniq-reads-signal-<Experiment target>-<Experiment accession>-<File accessions>.h
        b. all-reads-signal-<Experiment target>-<Experiment accession>-<File accessions>.h
        c. signal-<Experiment target>-<Experiment accession>-<File accessions>.h5

    Note that name of cell-line is not included here. Therefore, use the
    directory name as a identifiers for cell-lines or species. The Experiment
    and File accession can be used to back-track about the dataset on ENCODE
    website.

Requirements
============
1) bigWigToWig : It converts binary bigWig file to ascii Wig file.
2) bigWigInfo : It fetches the information about chromosomes from bigWig file.

Both tools can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/
for linux and Mac platform. However, these tools are not yet available for
Windows OS.

Path to these tools can be set using gcMapExplorer configure utility or can be
given with the command.

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
3) amean  -> Arithmetic mean or average
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
"""Input text file.
At first search the datasets on https://www.encodeproject.org/matrix/?type=Experiment.
Then click on download button on top of the page. A text file will be downloaded.
This text file can be used as input in this program.

"""

assemblyHelp = \
""" Name of reference genome.
Example: hg19, GRCh38 etc.

"""

bigWigToWigHelp = \
"""Path to bigWigToWig tool.

This is not necessary when bigWigToWig path is already set using gcMapExplorer
configure utility.

It can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/
for linux and Mac platform.

If it is not present in configuration file, the input path should
be provided. It will be stored in configuration file for later use.

"""

bigWigInfoHelp = \
""" Path to bigWigInfo tool.

This is not necessary when bigWigInfo path is already set using gcMapExplorer
configure utility.

It can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/
for linux and Mac platform.

If it is not present in configuration file, the input path should
be provided. It will be stored in configuration file for later use.

"""

resolutionHelp = \
"""Additional input resolutions other than these resolutions: 1kb', '2kb',
'4kb', '5kb', '8kb', '10kb', '20kb', '40kb', '80kb', '100kb', '160kb','200kb',
'320kb', '500kb', '640kb',  and '1mb'.

Resolutions should be provided in comma separated values. For Example:
-r "25kb, 50kb, 75kb"

"""

coarseningMethodHelp = \
"""Methods to coarse or downsample the data for converting from 1-base
to coarser resolutions. If this option is not provided, all six methods (see
above) will be considered. User may use only subset of these methods.
For example: -dm "max, amean" can be used for downsampling by only these
two methods.

"""

methodToCombineHelp = \
"""Methods to combine data from more than two input file. Presently, three
methods can be used: 'mean', 'max' and 'min' for average, maximum and minimum
value, respectively.

"""
outDirHelp = \
""" Directory to save all h5 files. It is an essential input.

"""

keepOriginalHelp = \
"""To copy original 1-base resolution data in h5 file. This will increase the
file size significantly.

"""

def main():
    # Construct command line arguments and parsed it
    parser, args = parseArguments()

    # Check for bigWigToWig program
    bigWigToWig = args.bigWigToWig
    if bigWigToWig == 'None':
        showErrorAndExit(parser, '\nPath to bigWigToWig not provided!!!\n')
    bigWigToWig = os.path.normpath(bigWigToWig)
    checkFileExist(bigWigToWig, parser)

    # Check for bigWigInfo program
    bigWigInfo = args.bigWigInfo
    if bigWigInfo == 'None':
        showErrorAndExit(parser, '\nPath to bigWigInfo not provided!!!\n')
    bigWigInfo = os.path.normpath(bigWigInfo)
    checkFileExist(bigWigInfo, parser)

    # Check for input bigWig File
    if args.inputFile is None:
        showErrorAndExit(parser, '\nInput file is not given!!!\n')
    inputFile = args.inputFile.name
    args.inputFile.close()

    # Check for additional resolutions
    resolutions = args.resolutions
    if resolutions is not None:
        resolutions = get_resolution_list(resolutions, parser)

    # Check for input coarsening methods
    coarsening_methods = args.coarsening_methods
    if coarsening_methods is not None:
        coarsening_methods = get_coarsening_methods_list( coarsening_methods, parser )

    # Check for output directory
    outDir = args.outDir
    if outDir is not None:
        outDir = os.path.normpath(outDir)
        if not os.path.isdir(outDir):
            showErrorAndExit(parser, '\nOutput directory not found/accessible.\n')
    else:
        outDir = os.getcwd()
        print(' WARNING: No directory is provided for output file.\n \
        All files will be stored here: [{0}]'.format(outDir))

    # Check for scratch directory
    workDir = args.workDir
    if not os.path.isdir(workDir):
        showErrorAndExit(parser, '\nScratch Directory "{0}" not found !!!\n'.format(workDir))

    encodeDatasets = gmlib.genomicsDataHandler.EncodeDatasetsConverter(inputFile, args.assembly, pathTobigWigToWig=bigWigToWig,
                                                            methodToCombine=args.methodToCombine, pathTobigWigInfo=bigWigInfo, workDir=workDir)
    encodeDatasets.saveAsH5(outDir, resolutions=resolutions, coarsening_methods=coarsening_methods, compression=args.compression, keep_original=args.keep_original)
    del encodeDatasets


def get_resolution_list(resolution, parser):
    rlist = []
    temp = re.split(',', resolution)
    for r in temp:
        r = r.rstrip().lstrip()
        if not r.strip():
            continue
        try:
            b = gmlib.ccmap.resolutionToBinsize(r)
            rlist.append( gmlib.ccmap.binsizeToResolution(b) )
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

def parseArguments():
    parser = argparse.ArgumentParser(
                prog='gcMapExplorer encode2h5',
                description=description,
                formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('r'), metavar='input.txt',
                        dest='inputFile', help=inputFileHelp)

    parser.add_argument('-amb', '--assembly', action='store',
                        metavar = 'hg19', default='hg19',
                        dest='assembly',help=assemblyHelp)

    parser.add_argument('-b2w', '--bigWigToWig', action='store',
                        type=str, metavar='bigWigToWig',
                        default=config['Programs']['bigWigToWig'],
                        dest='bigWigToWig', help=bigWigToWigHelp)

    parser.add_argument('-binfo', '--bigWigInfo', action='store', type=str,
                        metavar='bigWigInfo', dest='bigWigInfo',
                        default=config['Programs']['bigWigInfo'],
                        help=bigWigInfoHelp)

    parser.add_argument('-r', '--resolutions', action='store', type=str,
                        metavar='"List of Resolutions"', dest='resolutions',
                        help=resolutionHelp)

    parser.add_argument('-dm', '--downsample-method', action='store', type=str,
                        metavar='"List of downsampling method"',
                        dest='coarsening_methods',
                        help=coarseningMethodHelp)

    parser.add_argument('-cmeth', '--compression-method', action='store',
                        dest='compression', metavar='lzf',
                        choices=['lzf', 'gzip'], default='lzf',
                        help='Data compression method in h5 file.\n')

    parser.add_argument('-mtc', '--method-to-combine', action='store',
                        dest='methodToCombine', metavar='mean',
                        choices=['mean', 'max', 'min'], default='mean',
                        help=methodToCombineHelp)

    parser.add_argument('-od', '--outDir', action='store', dest='outDir',
                        metavar='outDir', help=outDirHelp)

    parser.add_argument('-ko', '--keep-original', action='store_true',
                        default = False,
                        dest='keep_original', help=keepOriginalHelp)

    parser.add_argument('-wd', '--work-dir', action='store', dest='workDir',
                        default=config['Dirs']['WorkingDirectory'],
                        metavar=config['Dirs']['WorkingDirectory'],
                        help='Directory where temporary files will be stored.')

    idx = sys.argv.index("encode2h5")+1
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
