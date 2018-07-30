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
from math import inf

import numpy as np

from gcMapExplorer.config import getConfig
from gcMapExplorer.lib import util, ccmap, gcmap

config = getConfig()


def gen_map_from_generator(values, c1_norm=None, c2_norm=None, resolution=None, map_type='intra', work_dir=None,
                           log_handler=None):
    """Generate ccmap from generator

    Adapted from work by Rajendra Kumar.

    Parameters
    ----------
    values : BlockReader
        generator yielding (x, y, count)
    c1_norm : array of float, optional
        norm vector for x values
    c2_norm : array of float, optional
        norm vector for y values
    resolution : int, optional
        the resolution
    map_type : str, optional
        the map type: "intra" or "inter" chromosomal map
    work_dir : str, optional
        directory where temporary files will be stored.
    log_handler
        the log handler

    Returns
    -------
    ccmap.CCMAP

    See Also
    --------
    gen_map_from_locations_value in gcMapExplorer.lib.importer
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


def hic2ccmap(hic, chr1, chr2, resolution="finest", norm_type=None, log_handler=None):
    """Generate ccmap from hic

    Parameters
    ----------
    hic : HicParser
        the hic parser object
    chr1 : str
        chromosome name
    chr2 : str
        chromosome name
    resolution : int or str, optional
        the resolution (default "finest" uses the finest available)
    norm_type : NormType, optional
        the norm type
    log_handler : optional
        the log handler

    Returns
    -------
    ccmap.CCMAP

    See Also
    --------
    hic2gcmap
    """
    logger = logging.getLogger('hic2gcmap')
    if log_handler is not None:
        logger.propagate = False
        logger.addHandler(log_handler)
    logger.setLevel(logging.INFO)
    logger.info("Converting map for pair {} {} ...".format(chr1, chr2))
    if log_handler is not None:
        logger.removeHandler(log_handler)

    res = min(hic.bp_resolutions) if resolution == "finest" else resolution

    record = hic.record(chr1, chr2)
    values = hic.blocks(record, res)

    if norm_type is not None:
        c1_norm = hic.norm_vector(chr1, norm_type, res)
        c2_norm = hic.norm_vector(chr2, norm_type, res)
    else:
        c1_norm, c2_norm = None, None

    if chr1 == chr2:
        map_type = 'intra'
    else:
        map_type = 'inter'

    ccmap_obj = gen_map_from_generator(values, c1_norm=c1_norm, c2_norm=c2_norm, map_type=map_type,
                                       log_handler=log_handler)
    ccmap_obj.xlabel = chr1
    ccmap_obj.ylabel = chr2

    return ccmap_obj


def hic2gcmap(hic, chr1, chr2, output_file, resolution="finest", norm_type=None, compression="lzf", downsampling="sum",
              log_handler=None):
    """Generate gcmap from hic and write to file

    Parameters
    ----------
    hic : HicParser
        the hic parser object
    chr1 : str
        chromosome name
    chr2 : str
        chromosome name
    output_file : str
        the output filename
    resolution : int or str, optional
        the resolution (default "finest" uses the finest available)
    norm_type : NormType, optional
        the norm type
    compression : str, optional
        compression type (lzf or gzip)
    downsampling : str, optional
        downsampling method (sum, mean, max, or none)
    log_handler : optional
        the log handler

    Returns
    -------
    bool
        True or False depending on success

    See Also
    --------
    hic2ccmap
    """
    ccmap_obj = hic2ccmap(hic, chr1, chr2, resolution, norm_type, log_handler)

    return gcmap.addCCMap2GCMap(ccmap_obj, output_file, compression=compression, generateCoarse=downsampling != "none",
                                coarseningMethod=downsampling, logHandler=log_handler)
