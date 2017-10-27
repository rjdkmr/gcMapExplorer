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
import numpy

here = path.abspath(path.dirname(__file__))

def read(fname):
    return open(path.join(path.dirname(__file__), fname)).read()

exec(open('gcMapExplorer/_version.py').read())

ext_modules = [
    Extension("gcMapExplorer.lib.ccmapHelpers",            ["gcMapExplorer/lib/ccmapHelpers.pyx"]             ),
    Extension("gcMapExplorer.lib.normalizeKnightRuiz",     ["gcMapExplorer/lib/normalizeKnightRuiz.pyx"]      ),
    Extension("gcMapExplorer.lib.normalizeCore",           ["gcMapExplorer/lib/normalizeCore.pyx"]            ),
    Extension("gcMapExplorer.lib._corrMatrixCore",         ["gcMapExplorer/lib/_corrMatrixCore.pyx", "gcMapExplorer/lib/corrMatrixCoreSRC.c"] , include_dirs=[numpy.get_include()],
                extra_compile_args=['-fopenmp'], extra_link_args=['-fopenmp']),
    Extension("gcMapExplorer.lib.TadFinder",               ["gcMapExplorer/lib/TadFinder.pyx"]                ),
    ]

setup(
    name = 'gcMapExplorer',
    version = __version__,

    # Required packages
    install_requires = [ 'appdirs>=1.4', 'numpy>=1.6',  'scipy>=0.9', 'matplotlib>=1.1.0', 'dask>=0.7.3', 'toolz>=0.7.4', 'psutil>=5.2.0', 'h5py>=2.2.1', 'Cython>=0.23.0' ],
    #ext_modules = cythonize("gcMapExplorer/lib/*.pyx", compiler_directives={'embedsignature': True}),
    ext_modules = cythonize(ext_modules, compiler_directives={'embedsignature': True}),
    include_dirs=[numpy.get_include()],
    packages=find_packages(),

    package_data = { '': ['*.ico', '*.png', '*.ui'] },

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={'console_scripts': [ 'gcMapExplorer=gcMapExplorer:main.main',], },

    include_package_data=True,

    # metadata for upload to pypi
    author = "Rajendra Kumar",
    author_email = "rjdkmr@gmail.com",
    url = 'https://github.com/rjdkmr/gcMapExplorer',
    description = "A platform to visualize and analyze genome contact maps",
    long_description = read('README.rst'),
    keywords = ["Hi-C", "Genome Contact Map Explorer", "Contact Map Explorer", "3D Genome Organization"],
    license = 'GNU General Public License v3 (GPLv3)',
    classifiers = [
        'Environment :: Console',
        'Environment :: X11 Applications :: Qt',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
