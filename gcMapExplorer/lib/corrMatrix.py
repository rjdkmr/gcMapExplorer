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

import os, re
import numpy as np
import scipy as sp
import logging

from . import gcmap as gmp
from .ccmap import checkCCMapObjectOrFile
from . import util

from ._corrMatrixCore import calculateCorrMatrix, calculateCovMatrix, _calculateCorrMatrixOLDSLOW
from ._corrMatrixCore import calculateCorrelation, calculateCovariance


__all__ = [ 'calculateCorrMatrix',
            'calculateCovMatrix',
            'calculateCorrelation',
            'calculateCovariance',
            'calculateCorrMatrixForCCMap',
            'calculateCorrMatrixForGCMaps'
          ]

logger = logging.getLogger('CorrMatrix')
logger.setLevel(logging.INFO)

def calculateCorrMatrixForCCMap(inputCCMap, logspace=False, maskvalue=0.0, vmin=None, vmax=None, outFile=None, workDir=None):
    """ Calculate correlation matrix of a contact map.
    It calculates correlation between all rows and columns of contact map.


    Parameters
    ----------
    ccMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ccmap file
        A CCMAP object containing observed contact frequency or a ccmap file.
    logspace : bool
        If its value is ``True``, at first map is converted as logarithm of map
        and subsequently correlation will be calculated.
    maskvalue : float
        Do not consider bins with this value during calculation. By default here
        it is zero because bins with zero is considered to be have missing data.
    vmin : float
        Minimum threshold value for normalization. If contact frequency is less than
        or equal to this threshold value, this value is discarded during normalization.
    vmax : float
        Maximum threshold value for normalization. If contact frequency is greater than
        or equal to this threshold value, this value is discarded during normalization.
    outFile : str
        Name of output ccmap file, to save directly the correlation matrix as a ccmap file. In case of this option, ``None`` will return.
    workDir : str
        Path to the directory where temporary intermediate files are generated.
        If ``None``, files are generated in the temporary directory according to the main configuration.


    Returns
    -------
    ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ``None``
        Normalized Contact map. When ``outFile`` is provided, ``None`` is returned. In case of any other error, ``None`` is returned.

    """

    # Check whether input is a file or a object
    ccMapObjOrig, ccmapType = checkCCMapObjectOrFile(inputCCMap, workDir=workDir)

    # Make another copy here for maximum and minimum threshold value
    # Also mask all values that are outside of cutoff
    if vmin is not None or vmax is not None:
        ccMapObj = ccMapObjOrig.copy()
        ccMapObj.make_editable()

        if ccmapType == 'File':
            del ccMapObjOrig

        ccmapType = 'File' # This temporary file should be deleted when necessary
        if vmin is not None:
            ccMapObj.matrix[ np.nonzero(ccMapObj.matrix <= vmin) ] = maskvalue
        if vmax is not None:
            ccMapObj.matrix[ np.nonzero(ccMapObj.matrix >= vmax) ] = maskvalue
        ccMapObj.matrix.flush()
        ccMapObj.make_unreadable()

    else:
        ccMapObj = ccMapObjOrig

    try:
        ccMapObj.make_readable()
        matrix = (ccMapObj.matrix[~ccMapObj.bNoData,:])[:,~ccMapObj.bNoData]   # Selected row-column which are not all zeros
        matrix = matrix.T.astype(np.float64, order='C')

        if logspace:
            matrix = np.log10(matrix)
            matrix[np.isinf(matrix)] = 0

        #diagIdx = util.kth_diag_indices(0, matrix)
        #matrix[diagIdx] = 0

        corrMatrix = calculateCorrMatrix(matrix, maskvalue=maskvalue)
        corrMatrix[np.isnan(corrMatrix)] = 0

        corrCCMap = ccMapObj.copy(fill=0.0)
        corrCCMap.make_editable()

        dsm_i = 0
        ox = corrCCMap.matrix.shape[0]
        idx_fill = np.nonzero( ~(corrCCMap.bNoData) )
        for i in range(ox):
            if not corrCCMap.bNoData[i]:
                corrCCMap.matrix[i, idx_fill] = corrMatrix[dsm_i]
                corrCCMap.matrix[idx_fill, i] = corrMatrix[dsm_i]
                dsm_i += 1

        corrCCMap.matrix.flush()
        corrCCMap.minvalue = -1
        corrCCMap.maxvalue = 1

    except (KeyboardInterrupt, SystemExit) as e:
        # In case of program termination, delete the newly created ccmap and raise error
        if 'corrCCMap' in locals():	del corrCCMap
        if ccmapType == 'File' and 'ccMapObj' in locals():	del ccmap
        raise e

    # In case of other error, delete the newly created ccmap, print error and return None
    except Exception as e:
        if 'corrCCMap' in locals():	del corrCCMap
        if ccmapType == 'File' and 'ccMapObj' in locals():	del ccmap
        logger.warning(e)
        logger.warning('Error in Correlation Matrix Calculation !!!')
        return None

    # Save output ccmap file
    if outFile is not None:
        cmp.save_ccmap(corrCCMap, outFile, compress=True)
        del corrCCMap
        return
    else:
        return corrCCMap

def calculateCorrMatrixForGCMaps(gcMapInputFile, gcMapOutFile, logspace=False, maskvalue=0.0, vmin=None, vmax=None, replaceMatrix=False, compression='lzf', workDir=None):
    """ Calculate Correlation matrix for all maps present in input gcmap file
    It calculates correlation between all rows and columns of contact map.


    Parameters
    ----------
    gcMapInputFile : str
        Input gcmap file.
    gcMapOutFile : str
        Output gcmap file.
    logspace : bool
        If its value is ``True``, at first map is converted as logarithm of map
        and subsequently correlation will be calculated.
    maskvalue : float
        Do not consider bins with this value during calculation. By default here
        it is zero because bins with zero is considered to be have missing data.
    vmin : float
        Minimum threshold value for normalization. If contact frequency is less than
        or equal to this threshold value, this value is discarded during normalization.
    vmax : float
        Maximum threshold value for normalization. If contact frequency is greater than
        or equal to this threshold value, this value is discarded during normalization.
    replaceMatrix : bool
        If its value is ``True``, the map will be replaced in output file.
        Otherwise, if a map is present, calculation will be skipped.
    compression : str
        Compression method in output gcmap file. Presently allowed : ``lzf`` for LZF compression and ``gzip`` for GZIP compression.
    workDir : str
        Path to the directory where temporary intermediate files are generated.
        If ``None``, files are generated in the temporary directory according to the main configuration.


    Returns
    -------
    None

    """

    # Get list of maps in ascending order
    gcmap = gmp.GCMAP(gcMapInputFile)
    gcmap.loadSmallestMap()
    mapList = gcmap.mapNameList.copy()

    for mapName in mapList:
        # Iterate over available resolutions
        gcmap.changeMap(mapName)


        while True:
            # Check if map is present in output file
            # It check for all resolutions.
            tf = gmp.GCMAP(gcMapOutFile)
            doExist = False
            while True:
                logger.info('Calculating for {0}-{1}'.format(mapName, gcmap.resolution))
                doExist = tf.checkMapExist(mapName=mapName, resolution=gcmap.resolution)
                if doExist:
                    if not gcmap.toCoarserResolution():
                        break
                else:
                    break

            # If map is present, either replace or skip according to input
            if doExist:
                if not replaceMatrix:
                    break
                else:
                    tf.hdf5[mapName].pop(gcmap.resolution)

            del tf

            ccMap = gmp.loadGCMapAsCCMap(gcmap.hdf5, mapName=mapName,
                                        resolution=gcmap.resolution,
                                        workDir=workDir)

            # Because ccMap is already loaded here, directly edit matrix here
            # No need to pass vmin and vmx during next step as
            # it is already taken care here
            if vmin is not None or vmax is not None:
                ccMap.make_editable()
                if vmin is not None:
                    ccMap.matrix[np.nonzero(ccMap.matrix <= vmin)] = mask
                if vmax is not None:
                    ccMap.matrix[np.nonzero(ccMap.matrix >= vmax)] = mask
                ccMap.matrix.flush()
                ccMap.make_unreadable()

            try:
                corrCCMap = calculateCorrMatrixForCCMap(ccMap, logspace=logspace, maskvalue=maskvalue)
                if corrCCMap is not None:
                    gmp.addCCMap2GCMap(corrCCMap, gcMapOutFile,
                                            scaleoffset=4,
                                            compression=compression,
                                            generateCoarse=False,
                                            replaceCMap=False)
                    del corrCCMap

                if not gcmap.toCoarserResolution():
                    break

            # In case of program termination, delete the newly created ccmap and raise error
            except (KeyboardInterrupt, SystemExit) as e:
                if 'corrCCMap' in locals(): del corrCCMap
                if 'gcmap' in locals(): del gcmap
                if 'ccMap' in locals(): del ccMap
                raise e

            del ccMap


def calculateCorrMatrixOLDSLOW(inputCCMap, logspace=False, outFile=None, workDir=None):
    """ Calculate correlation matrix of a contact map.
    It calculates correlation between all rows and columns of contact map.


    Parameters
    ----------
    ccMap : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ccmap file
        A CCMAP object containing observed contact frequency or a ccmap file.
    logspace : bool
        If its value is ``True``, at first map is converted as logarithm of map
        and subsequently correlation will be calculated.
    outFile : str
        Name of output ccmap file, to save directly the correlation matrix as a ccmap file. In case of this option, ``None`` will return.
    workDir : str
        Path to the directory where temporary intermediate files are generated.
        If ``None``, files are generated in the temporary directory according to the main configuration.


    Returns
    -------
    ccMapObj : :class:`gcMapExplorer.lib.ccmap.CCMAP` or ``None``
        Normalized Contact map. When ``outFile`` is provided, ``None`` is returned. In case of any other error, ``None`` is returned.

    """

    # Check whether input is a file or a object
    ccmap, ccmapType = checkCCMapObjectOrFile(inputCCMap, workDir=workDir)

    try:
        ccmap.make_readable()
        matrix = (ccmap.matrix[~ccmap.bNoData,:])[:,~ccmap.bNoData]   # Selected row-column which are not all zeros
        matrix = matrix.astype(np.float)
        mshape = matrix.shape

        if logspace:
            matrix = np.log10(matrix)
            matrix[np.isinf(matrix)] = 0

        corrMatrix = _calculateCorrMatrixOLDSLOW(matrix)
        corrMatrix[np.isnan(corrMatrix)] = 0

        corrCCMap = ccmap.copy(fill=0.0)
        corrCCMap.make_editable()

        dsm_i = 0
        ox = corrCCMap.matrix.shape[0]
        idx_fill = np.nonzero( ~(corrCCMap.bNoData) )
        for i in range(ox):
            if not corrCCMap.bNoData[i]:
                corrCCMap.matrix[i, idx_fill] = corrMatrix[dsm_i]
                corrCCMap.matrix[idx_fill, i] = corrMatrix[dsm_i]
                dsm_i += 1

        corrCCMap.matrix.flush()
        corrCCMap.minvalue = -1
        corrCCMap.maxvalue = 1

    except (KeyboardInterrupt, SystemExit) as e:
        # In case of program termination, delete the newly created ccmap and raise error
        if 'corrCCMap' in locals():	del corrCCMap
        if ccmapType == 'File' and 'ccmap' in locals():	del ccmap
        raise e

    # In case of other error, delete the newly created ccmap, print error and return None
    except Exception as e:
        if 'corrCCMap' in locals():	del corrCCMap
        if ccmapType == 'File' and 'ccmap' in locals():	del ccmap
        logger.warning(e)
        logger.warning('Error in Correlation Matrix Calculation !!!')
        return None

    # Save output ccmap file
    if outFile is not None:
        cmp.save_ccmap(corrCCMap, outFile, compress=True)
        del corrCCMap
        return
    else:
        return corrCCMap
