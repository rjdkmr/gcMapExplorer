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

import os, re, sys
import tempfile
import configparser
from appdirs import system
from appdirs import AppDirs
import psutil


defaultConfigText="""
[Dirs]
WorkingDirectory = {0}

[Programs]
bigWigInfo = None
bigWigToWig = None

""".format(tempfile.gettempdir())

configAppDirs = AppDirs('gcMapExplorer')

if system == 'win32':
    configFileName = 'gcMapExplorer.ini'
else:
    configFileName = 'gcMapExplorer.conf'

configFile = os.path.join(configAppDirs.user_config_dir, configFileName)


def _defaultConfiguration():
    config = configparser.ConfigParser()
    config.read_string(defaultConfigText)
    return config

def updateConfig(section, option, value):
    """ Update configuration file

    Parameters
    ----------
    section : str
        Section of the configuration files. It could be ``Dirs`` or ``Programs``.
    option : str
        Input option, for which value is to be changed.
    value : str or int
        New value of the input option.
    """

    config = getConfig()

    if section not in config:
        return
    if option not in config[section]:
        return

    config[section][option] = value

    with open(configFile, 'w') as fout:
        config.write(fout)

def getConfig():
    """ To get the present configuration.

    Configuration file has the following organization.

    ::

        Configuration
            ├─────────── Dirs
            │             └──────── WorkingDirectory
            │
            └─────────── Programs
                          ├──────── bigWigInfo
                          └──────── bigWigToWig


    In case no configuration file is found, a new file is generated and default
    value is assigned to each option.

    Returns
    -------
    config : dict
        Dictionary of Dictionaries with option name and value pair.
        For example, config['Dirs']['WorkingDirectory'] contains path to
        scratch directory. Similarly, config['Programs']['bigWigInfo']
        contains path to bigWigInfo program.


    """
    if not os.path.exists(configFile):

        if not os.path.exists(configAppDirs.user_config_dir):
            try:
                os.makedirs(configAppDirs.user_config_dir)
            except OSError:
                pass


        config = _defaultConfiguration()
        print(" No configuration file found... Generating a new configuration file with default values as follows:")
        print("==============================")
        print("Configuration file: {0}".format(configFile))
        print("==============================")
        config.write(sys.stdout)
        print("==============================")

        with open(configFile, 'w') as fout:
            config.write(fout)

    else:
        config = configparser.ConfigParser()
        config.read(configFile)

    return config

def printConfig():
    """Print configuration file

    It can be used to print the configuration file. It shows the current
    configuration of gcMapExplorer.

    """
    config = getConfig()
    print("==============================")
    print("Configuration file: {0}".format(configFile))
    print("==============================")
    config.write(sys.stdout)
    print("==============================")


def cleanScratch():
    """ Clean scratch directory.

    It checks whether any other gcMapExplorer process is running. In case,
    when only one process (i.e. current) is running, all files with "gcx"
    prefix will be deleted from default scratch directory.

    """
    config = getConfig()
    count = 0
    for pid in psutil.pids():
        p = None
        try:
            p = psutil.Process(pid)
        except psutil.NoSuchProcess:
            pass

        if p is not None:
            if 'gcMapExplorer' in p.name():
                count += 1

    # If only one gcMapExplorer is running, it is the current one
    if count == 1:
        for f in os.listdir(config['Dirs']['WorkingDirectory']):
            if not os.path.isfile(f):
                continue
            basename = os.path.basename(f)
            base = os.path.splitext(basename)[0]
            if base in ["gcx"]:
                try:
                    os.remove(f)
                except IOError:
                    pass

        print('     ... Finished Cleaning')
