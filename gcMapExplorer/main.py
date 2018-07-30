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

import sys, os, random
import numpy as np
import logging

from .config import cleanScratch

# Set logging level
logging.basicConfig(level=logging.INFO)

gui_tools =      ['browser', 'cmapImporter', 'cmapNormalizer', 'h5Converter']
convert_tools =  ['coo2cmap', 'pairCoo2cmap', 'homer2cmap', 'bc2cmap', 'hic2gcmap', 'bigwig2h5', 'wig2h5', 'bed2h5', 'encode2h5']
norm_tools =     ['normKR', 'normVC', 'normIC', 'normMCFS']
analysis_tools = ['corrBWcmaps']
utility_tools = ['config']

def main():
    options = {'browser':'Interactive Browser for genomic contact maps',
               'cmapImporter':'Interface to import contact maps and datasets',
               'cmapNormalizer' : 'Interface to normalize contact maps',
               'h5Converter' : 'Convert bigWig/wig/bed file to browser compatible h5 format',
               'coo2cmap':'Import COO sparse matrix format to ccmap or gcmap',
               'pairCoo2cmap': 'Import map from files similar to paired sparse matrix Coordinate (COO) format',
               'homer2cmap':'Import HOMER Hi-C interaction matrix to ccmap or gcmap',
               'bc2cmap': 'Import Bin-Contact format files to ccmap or gcmap',
               'hic2gcmap': 'Import hic file to gcmap',
               'bigwig2h5': 'Import a bigWig file to browser compatible h5 format',
               'wig2h5': 'Import a wig file to browser compatible h5 format',
               'bed2h5' : 'Import a bed file to browser compatible h5 format',
               'encode2h5' : 'Download and convert datasets from ENCODE project',
               'normKR' : 'Normalization using Knight-Ruiz matrix balancing',
               'normVC' : 'Normalization using Vanilla-Coverage method',
               'normIC' : 'Normalization using Iterative Correction',
               'normMCFS' : 'Normalization by Median Contact Frequency Scaling',
               'corrBWcmaps' : 'Calculate correlation between contact maps',
               'config': 'To print configuration file and clean scratch directory'}

    if len(sys.argv)<=1:
        show_help(options)
        sys.exit(-1)
    if sys.argv[1] not in options:
        print(' ERROR: "{0}" is not an accepted option.\n' .format(sys.argv[1]))
        show_help(options)
        sys.exit(-1)

    if sys.argv[1] == 'config':
        from .clui import config
        config.main()

    if sys.argv[1] == 'browser':
        from .gui import browser
        browser.main()

    if sys.argv[1] == 'cmapImporter':
        from .gui import importer_ui as impui
        impui.main()

    if sys.argv[1] == 'cmapNormalizer':
        from .gui import normalizer_ui
        normalizer_ui.main()

    if sys.argv[1] == 'h5Converter':
        from .gui import h5Converter
        h5Converter.main()

    if sys.argv[1] == 'coo2cmap':
        from .clui import coo2cmap
        coo2cmap.main()

    if sys.argv[1] == 'pairCoo2cmap':
        from .clui import pairCoo2cmap
        pairCoo2cmap.main()

    if sys.argv[1] == 'homer2cmap':
        from .clui import homer2cmap
        homer2cmap.main()

    if sys.argv[1] == 'bc2cmap':
        from .clui import bc2cmap
        bc2cmap.main()

    if sys.argv[1] == 'hic2gcmap':
        from .clui import hic2gcmap
        hic2gcmap.main()

    if sys.argv[1] == 'bigwig2h5':
        from .clui import bigwig2h5
        bigwig2h5.main()

    if sys.argv[1] == 'wig2h5':
        from .clui import wig2h5
        wig2h5.main()

    if sys.argv[1] == 'bed2h5':
        from .clui import bed2h5
        bed2h5.main()

    if sys.argv[1] == 'encode2h5':
        from .clui import encode2h5
        encode2h5.main()

    if sys.argv[1] == 'normKR':
        from .clui import normKR
        normKR.main()

    if sys.argv[1] == 'normVC':
        from .clui import normVC
        normVC.main()

    if sys.argv[1] == 'normIC':
        from .clui import normIC
        normIC.main()

    if sys.argv[1] == 'normMCFS':
        from .clui import normMCFS
        normMCFS.main()

    if sys.argv[1] == 'corrBWcmaps':
        from .clui import corrBWcmaps
        corrBWcmaps.main()

def show_help(options):
    print(' ==============================')
    print(' Usage:')
    print(' gcMapExplorer <Option>')
    print(' Use following options:\n')

    print('\n-------------------------')
    print(' Graphical User Interface')
    print('-------------------------')
    for tool in gui_tools:
        print('\t{0} : {1}\n'.format(tool, options[tool]))

    print('\n--------------------------')
    print(' To convert or import data')
    print('--------------------------')
    for tool in convert_tools:
        print('\t{0} : {1}\n'.format(tool, options[tool]))

    print('\n-------------------------')
    print(' To normalize contact map')
    print('-------------------------')
    for tool in norm_tools:
        print('\t{0} : {1}\n'.format(tool, options[tool]))

    print('\n-------------------------')
    print(' Analysis')
    print('-------------------------')
    for tool in analysis_tools:
        print('\t{0} : {1}\n'.format(tool, options[tool]))

    print('\n-------------------------')
    print(' Utility')
    print('-------------------------')
    for tool in utility_tools:
        print('\t{0} : {1}\n'.format(tool, options[tool]))

    print(' ==============================')

if __name__=="__main__":
    main()
