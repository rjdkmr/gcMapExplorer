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

import gcMapExplorer.lib as gmlib
import numpy as np
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.decomposition import PCA

def correctForPeaks(matrix, pc, correctionFactor):
    negIdx = np.nonzero(pc < 0)
    posIdx = np.nonzero(pc > 0)

    for i in range(matrix.shape[0]):
        posSum = np.median(matrix[i][posIdx])
        negSum = np.median(matrix[i][negIdx])

        # print(posSum, negSum, pc[i])
        factor = 1
        if posSum > negSum * correctionFactor and pc[i] < 0:
            factor = -1
            # print(i, posSum, negSum)

        if negSum > posSum * correctionFactor and pc[i] > 0:
            factor = -1
            # print(i, posSum, negSum)

        pc[i] = pc[i] * factor

    return pc


def calculateCompartments(ccmap, method='k-means', peakCorrectionFactor=None):
    if method not in ['k-means', 'hierarchical', 'PCA']:
        raise ValueError('Method should be : {0}'.format(['k-means', 'hierarchical']))

    ccmap.make_readable()
    matrix = (ccmap.matrix[~ccmap.bNoData, :])[:, ~ccmap.bNoData]  # Selected row-column which are not all zeros
    matrix[np.isnan(matrix)] = 0

    if method == 'k-means':
        cls = KMeans(n_clusters=2).fit(matrix)

    if method == 'hierarchical':
        cls = AgglomerativeClustering(n_clusters=2).fit(matrix)

    if method in ['k-means', 'hierarchical']:
        labels = cls.labels_[:]
        labels[np.nonzero(labels == 0)] = -1

    if method == 'PCA':
        pca = PCA(n_components=2)
        pca.fit(matrix)
        t_pc = pca.transform(matrix)
        labels = t_pc.T[0][:]

    # Correction for wrong peaks
    if peakCorrectionFactor is not None:
        labels = correctForPeaks(matrix, labels, peakCorrectionFactor)

    j = 0
    compartment = np.zeros((ccmap.matrix.shape[0],))
    for i in range(ccmap.matrix.shape[0]):
        if not ccmap.bNoData[i]:
            compartment[i] = labels[j]
            j = j + 1

    return compartment


def determineABforGCMap(gcMapInputFile, h5OutFile, minResolution, method='k-means', peakCorrectionFactor=None,
                        compression='lzf', workDir=None, logHandler=None):
    if method not in ['k-means', 'hierarchical', 'PCA']:
        raise ValueError('Method should be : {0}'.format(['k-means', 'hierarchical']))

    # Get list of maps in ascending order
    gcmap = gmlib.gcmap.GCMAP(gcMapInputFile)
    gcmap.loadSmallestMap()
    mapList = gcmap.mapNameList.copy()

    hdf5Handle = gmlib.genomicsDataHandler.HDF5Handler(h5OutFile)
    hdf5Handle.open()

    minBinsize = gmlib.ccmap.resolutionToBinsize(minResolution)

    for mapName in mapList:

        print(mapName)

        # Iterate over available resolutions
        gcmap.changeMap(mapName)

        while minBinsize > gmlib.ccmap.resolutionToBinsize(gcmap.resolution):
            gcmap.toCoarserResolution()

        previousResolution = gcmap.resolution

        while True:
            ccMap = gmlib.gcmap.loadGCMapAsCCMap(gcmap.hdf5, mapName=mapName, resolution=gcmap.resolution,
                                                 workDir=workDir)
            data = calculateCompartments(ccMap, method=method, peakCorrectionFactor=peakCorrectionFactor)
            del ccMap

            if peakCorrectionFactor is not None:
                name = 'Compartment by ' + str(method) + ' ' + str(peakCorrectionFactor)
            else:
                name = 'Compartment by ' + str(method)

            hdf5Handle.addDataByArray(mapName, gcmap.resolution, name, data, compression=compression)

            gcmap.toCoarserResolution()
            if previousResolution == gcmap.resolution:
                break
            else:
                previousResolution = gcmap.resolution

    del gcmap
    del hdf5Handle