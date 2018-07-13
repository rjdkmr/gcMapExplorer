#!/usr/bin/env python3
#
# This file is part of gcMapExplorer
# Copyright (C) 2016-2018  Rajendra Kumar, Ludvig Lizana, Per Stenberg
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
# =============================================================================
import logging
import sys
from argparse import ArgumentParser, ArgumentTypeError, FileType
from math import inf
from os import path

import numpy as np

from gcMapExplorer.config import getConfig
from gcMapExplorer.lib import util, ccmap, gcmap
from gcMapExplorer.lib.hicparser import HicParser, Unit, NormType

config = getConfig()


def gen_map_from_generator(values, c1_norm=None, c2_norm=None, resolution=None, map_type='intra', work_dir=None,
                           log_handler=None):
    """
    Create CCMAP from generator

    Adapted from work by Rajendra Kumar
    See gen_map_from_locations_value in gcMapExplorer.lib.importer

    :param values: generator yielding (x, y, count)
    :param c1_norm: norm vector for x
    :param c2_norm: norm vector for y
    :param resolution: Resolution of Hi-C map
    :param map_type: Hi-C map type: "intra" or "inter" chromosomal map
    :param work_dir: Directory where temporary files will be stored.
        If it is not provided, this value is taken from configuration file.
    :param log_handler: Optional log handler
    :return: CCMAP
    """

    logger = logging.getLogger('genMapFromGenerator')
    if log_handler is not None:
        logger.propagate = False
        logger.addHandler(log_handler)
    logger.setLevel(logging.INFO)

    if work_dir is None:
        work_dir = config['Dirs']['WorkingDirectory']

    ccmap_obj = ccmap.CCMAP()

    maxLi = maxLj = -inf
    minLi = minLj = inf
    i, j, value = [], [], []

    if c1_norm is not None:
        logger.info("Normalizing and finding min and max values ...")
    else:
        logger.info("Finding min and max values ...")

    for bin_x, bin_y, count in values:
        x = bin_x * values.bin_size
        y = bin_y * values.bin_size
        maxLi = max(maxLi, x)
        minLi = min(minLi, x)
        maxLj = max(maxLj, y)
        minLj = min(minLj, y)
        i.append(x)
        j.append(y)
        if c1_norm is not None and c2_norm is not None:
            norm = c1_norm[bin_x] * c2_norm[bin_y]
            if norm != 0.0:
                value.append(count / norm)
            else:
                value.append(inf)
        else:
            value.append(count)

    maxL = max(maxLi, maxLj)
    minL = min(minLi, minLj)

    if resolution is not None:
        ccmap_obj.binsize = util.resolutionToBinsize(resolution)
    else:
        step = []
        for n in range(1, min(100, len(i))):
            di = abs(i[n] - i[n - 1])
            if di > 0:
                step.append(di)
            dj = abs(j[n] - j[n - 1])
            if dj > 0:
                step.append(dj)

        ccmap_obj.binsize = int(np.amin(step))

    logger.info("Total number of data in input file: {}".format(len(i)))

    if map_type == 'intra':
        logger.info("Minimum base-pair: {} and Maximum base-pair: {} are present in input data".format(minL, maxL))

        xticks = np.arange(0, maxL + ccmap_obj.binsize, ccmap_obj.binsize)
        yticks = np.arange(0, maxL + ccmap_obj.binsize, ccmap_obj.binsize)
        ccmap_obj.xticks = [0, maxL + ccmap_obj.binsize]
        ccmap_obj.yticks = [0, maxL + ccmap_obj.binsize]
    else:
        xticks = np.arange(0, maxLi + ccmap_obj.binsize, ccmap_obj.binsize)
        yticks = np.arange(0, maxLj + ccmap_obj.binsize, ccmap_obj.binsize)
        ccmap_obj.xticks = [0, maxLi + ccmap_obj.binsize]
        ccmap_obj.yticks = [0, maxLj + ccmap_obj.binsize]

    ccmap_obj.shape = len(xticks), len(yticks)

    # Creating a file name for binary numpy array on disk
    ccmap_obj.gen_matrix_file(workDir=work_dir)

    # Creating a file for binary numpy array on disk
    ccmap_obj.make_writable()

    logger.info("Shape of overall map: {}".format(ccmap_obj.shape))

    logger.info("Copying data ...")

    ccmap_obj.maxvalue = -inf
    ccmap_obj.minvalue = inf

    for n in range(len(i)):
        ni = int(i[n] / ccmap_obj.binsize)
        nj = int(j[n] / ccmap_obj.binsize)

        # In case if two observation is present for same pair of locations. take average of them
        if ccmap_obj.matrix[ni][nj] != 0:
            ccmap_obj.matrix[ni][nj] = (ccmap_obj.matrix[ni][nj] + value[n]) / 2
        else:
            ccmap_obj.matrix[ni][nj] = value[n]

        if map_type == 'intra':
            ccmap_obj.matrix[nj][ni] = ccmap_obj.matrix[ni][nj]

        ccmap_obj.maxvalue = max(ccmap_obj.maxvalue, value[n])

        if ccmap_obj.matrix[ni][nj] > 0:
            ccmap_obj.minvalue = min(ccmap_obj.minvalue, ccmap_obj.matrix[ni][nj])

    ccmap_obj.make_unreadable()

    if log_handler is not None:
        logger.removeHandler(log_handler)

    return ccmap_obj


def hic2gcmap(hic, chr1, chr2, output_file, resolution="finest", norm_type=None, compression="lzf", coarsening="sum",
              log_handler=None):
    """
    Convert hic file to gcmap file.

    :param hic: the hic file object
    :param chr1: chromosome name
    :param chr2: chromosome name
    :param output_file: output file
    :param resolution: the resolution
    :param norm_type: the norm type (see hicparser.NormType)
    :param compression: compression type (lzf or gzip)
    :param coarsening: coarsening method (sum, mean, max or none)
    :param log_handler: the log handler
    """
    logger = logging.getLogger('hic2gcmap')
    if log_handler is not None:
        logger.propagate = False
        logger.addHandler(log_handler)
    logger.setLevel(logging.INFO)

    logger.info("Converting map for pair {} {} ...".format(chr1, chr2))

    if log_handler is not None:
        logger.removeHandler(log_handler)

    res = resolution
    if resolution == "finest":
        res = min(hic.bp_resolutions)  # gcmap downsamples from this, so choose the finest

    try:
        values = hic.blocks(chr1, chr2, Unit.BP, res)
        if norm_type is not None:
            c1_norm = hic.norm_vector(chr1, norm_type, Unit.BP, res)
            c2_norm = hic.norm_vector(chr2, norm_type, Unit.BP, res)
        else:
            c1_norm, c2_norm = None, None
    except LookupError as e:
        print("Error: {}".format(e))
        return

    if chr1 == chr2:
        map_type = 'intra'
    else:
        map_type = 'inter'

    ccmap_obj = gen_map_from_generator(values, c1_norm=c1_norm, c2_norm=c2_norm, map_type=map_type)
    ccmap_obj.xlabel = chr1
    ccmap_obj.ylabel = chr2

    gcmap.addCCMap2GCMap(ccmap_obj, output_file, compression=compression, generateCoarse=coarsening != "none",
                         coarseningMethod=coarsening)


def resolution(string):
    if string == "finest":
        return string

    if string.endswith("kb"):
        string = string[:-2] + "000"

    msg = "Resulution must be a positive integer or integer ending with 'kb'"
    try:
        value = int(string)
        if value <= 0:
            raise ArgumentTypeError(msg)
        return value
    except ValueError:
        raise ArgumentTypeError(msg)


def norm(string):
    if string == "none":
        return None
    return NormType(string)


def main():
    parser = ArgumentParser(prog="gcMapExplorer hic2gcmap", description="Convert hic files to gcmap",
                            allow_abbrev=False)
    parser.add_argument("input", type=FileType("rb"), help="hic input file")
    parser.add_argument("output", type=str, nargs="?", help="output file or directory", default=".")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-c", "--chromosomes", type=str, nargs=2, metavar=("A", "B"), help="a pair of chromosomes A B")
    group.add_argument("-l", "--list", action="store_true", help="list all available chromosomes")
    parser.add_argument("--compression", type=str, choices=["lzf", "gzip"], default="lzf", metavar="C",
                        help="compression type, choose between lzf, gzip (default: lzf)")
    parser.add_argument("-r", "--resolution", type=resolution, default="finest", metavar="R",
                        help="the resolution, as an integer or as kb (default: finest)")
    parser.add_argument("-n", "--norm", type=str, choices=["VC", "VC_SQRT", "KR", "none"], metavar="N",
                        default="none",
                        help="the type of norm to use, choose between VC, VC_SQRT, KR, none (default: none)")
    parser.add_argument("--coarsening", type=str, choices=["sum", "mean", "max", "none"], default="sum", metavar="C",
                        help="the coarsening method to use, choose between sum, mean, max, none (default: sum)")

    args = parser.parse_args(args=sys.argv[sys.argv.index("hic2gcmap") + 1:])

    hic = HicParser(args.input)

    if path.isdir(args.output):
        filename = path.splitext(path.basename(args.input.name))[0]
        output_file = path.join(args.output, filename)
        if args.chromosomes:
            output_file += "_" + "_".join(args.chromosomes)
        if args.coarsening == "none" and args.resolution != "finest":
            output_file += "_" + str(args.resolution)[:-3] + "kb"
        if args.norm != "none":
            output_file += "_" + args.norm
        output_file += "_" + args.compression
        output_file += ".gcmap"
    else:
        output_file = args.output

    if args.list:
        if len(hic.chromosomes):
            print(", ".join(hic.chromosomes))

    elif args.chromosomes:
        chr1, chr2 = args.chromosomes
        hic2gcmap(hic, chr1, chr2, output_file, resolution=args.resolution, norm_type=norm(args.norm),
                  compression=args.compression,
                  coarsening=args.coarsening)
    else:
        for record in hic.records:
            chr1, chr2 = hic.chromosome_names(record)
            if "All" in (chr1, chr2):
                continue
            hic2gcmap(hic, chr1, chr2, output_file, resolution=args.resolution, norm_type=norm(args.norm),
                      compression=args.compression,
                      coarsening=args.coarsening)


if __name__ == "__main__":
    main()
