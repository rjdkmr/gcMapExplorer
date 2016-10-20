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

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from setuptools.extension import Extension
# To use a consistent encoding
from codecs import open
from os import path
from Cython.Build import cythonize

here = path.abspath(path.dirname(__file__))


ext_modules = [
    Extension("gcMapExplorer.lib.ccmapHelpers",            ["gcMapExplorer/lib/ccmapHelpers.pyx"],            compiler_directives={'embedsignature': True}),
    Extension("gcMapExplorer.lib.genomicsDataHandler",     ["gcMapExplorer/lib/genomicsDataHandler.pyx"],     compiler_directives={'embedsignature': True}),
    Extension("gcMapExplorer.lib.normalizeAverageContact", ["gcMapExplorer/lib/normalizeAverageContact.pyx"], compiler_directives={'embedsignature': True}),
    Extension("gcMapExplorer.lib.normalizeKnightRuiz",     ["gcMapExplorer/lib/normalizeKnightRuiz.pyx"],     compiler_directives={'embedsignature': True}),
    Extension("gcMapExplorer.lib.normalizeIC",             ["gcMapExplorer/lib/normalizeIC.pyx"],       compiler_directives={'embedsignature': True}),
    Extension("gcMapExplorer.lib.TadFinder",               ["gcMapExplorer/lib/TadFinder.pyx"],               compiler_directives={'embedsignature': True}),
    ]

setup(
    name = 'gcMapExplorer',
    version = '1.0a',

    # Required packages
    install_requires = [ 'appdirs>=1.4', 'numpy>=1.6',  'scipy>=0.9', 'matplotlib>=1.1.0', 'dask>=0.7.3', 'toolz>=0.7.4', 'h5py>=2.2.1' ],
    #ext_modules = cythonize("gcMapExplorer/lib/*.pyx", compiler_directives={'embedsignature': True}),
    ext_modules = ext_modules,
    packages=find_packages(),

    package_data = { '': ['*.ico', '*.png', '*.ui'] },

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={'console_scripts': [ 'gcMapExplorer=gcMapExplorer:main.main',], },

    # metadata for upload to pypi
    author = "Rajendra Kumar",
    author_email = "rjdkmr@gmail.com",
    description = "This package is used for analysis of Hi-C Map",
    keywords = "Hi-C genomics genome",
    include_package_data=True,

)
