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

import os, re, sys
import tempfile
import configparser
from appdirs import system
from appdirs import AppDirs


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


def defaultConfiguration():
    config = configparser.ConfigParser()
    config.read_string(defaultConfigText)
    return config

def updateConfig(section, option, value):
    if section not in config:
        return
    if option not in config[section]:
        return

    config[section][option] = value

    with open(configFile, 'w') as fout:
        config.write(fout)

def getConfig():
    if not os.path.exists(configFile):

        if not os.path.exists(configAppDirs.user_config_dir):
            try:
                os.makedirs(configAppDirs.user_config_dir)
            except OSError:
                pass


        config = defaultConfiguration()
        with open(configFile, 'w') as fout:
            config.write(fout)

    else:
        config = configparser.ConfigParser()
        config.read(configFile)

    return config


def initConfig():
    print('Configuration file not found...')
    print('Starting configuration:')

    while True:
        s = input('--> Write working directory name [current: {0}]'.format(config['Dirs']['WorkingDirectory']))
        if s:
            if not os.path.exists(s):
                print('Path {0} does not exist... Please retry...'.format(s))
            else:
                config['Dirs']['WorkingDirectory'] = s
                break

    with open(configFileName, 'w') as fout:
        config.write(fout)

    return config
