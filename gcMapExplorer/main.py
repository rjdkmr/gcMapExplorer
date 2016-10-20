#!/usr/bin/env python3
#
# Author: Rajendra Kumar
#
# This file is part of gcMapExplorer
# Copyright (C) 2016  Rajendra Kumar, Ludvig Lizana, Per Stenberg
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
logging.basicConfig(level=logging.INFO)

gui_tools =      ['browser', 'importer', 'normalizer']
convert_tools =  ['coo2cmap', 'homer2cmap', 'bc2cmap', 'bigwig2h5']
norm_tools =     ['normKR', 'normIC', 'normMCFS']
analysis_tools = ['corrBWcmaps']

def main():
    options = {'browser':'Interactive Browser for genomic contact maps',
               'importer':'Interface to import contact maps and datasets',
               'normalizer' : 'Interface to normalize contact maps',
               'coo2cmap':'Import COO sparse matrix format to ccmap or gcmap',
               'homer2cmap':'Import HOMER Hi-C interaction matrix to ccmap or gcmap',
               'bc2cmap': 'Import Bin-Contact format files to ccmap or gcmap',
               'bigwig2h5': 'Import a bigWig file to HDF5 format h5 file',
               'normKR' : 'Normalization using Knight-Ruiz matrix balancing',
               'normIC' : 'Normalization using Iterative Correction',
               'normMCFS' : 'Normalization by Median Contact Frequency Scaling',
               'corrBWcmaps' : 'Calculate correlation between contact maps' }

    if len(sys.argv)<=1:
        show_help(options)
        sys.exit(-1)
    if sys.argv[1] not in options:
        print(' ERROR: "{0}" is not an accepted option.\n' .format(sys.argv[1]))
        show_help(options)
        sys.exit(-1)

    if sys.argv[1] == 'browser':
        from .gui import browser
        browser.main()

    if sys.argv[1] == 'importer':
        from .gui import importer_ui as impui
        impui.main()

    if sys.argv[1] == 'normalizer':
        from .gui import normalizer_ui
        normalizer_ui.main()

    if sys.argv[1] == 'coo2cmap':
        from .clui import coo2cmap
        coo2cmap.main()

    if sys.argv[1] == 'homer2cmap':
        from .clui import homer2cmap
        homer2cmap.main()

    if sys.argv[1] == 'bc2cmap':
        from .clui import bc2cmap
        bc2cmap.main()

    if sys.argv[1] == 'bigwig2h5':
        from .clui import bigwig2h5
        bigwig2h5.main()

    if sys.argv[1] == 'normKR':
        from .clui import normKR
        normKR.main()

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

    print(' ==============================')



if __name__=="__main__":
    main()
