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
# =============================================================================

import sys
import argparse

from gcMapExplorer.config import cleanScratch, printConfig

description = \
"""To print configuration file and clean scratch directory

This can be used to print configuration file and its content.

=================================================
"""

cleanScratchHelp=\
"""Clean scratch directory
gcMapexplorer try to remove all temporary files present in scracth 
directory. However, due to crash and errors, sometimes these files cannot be removed.
Therefore, by this simple command, all temporary files from the scracth directory 
is removed.

Note: Temporary files could be huge and therefore it is advised to run this command periodically.

Note: Do not use this command when any of the gcMapExplorer tools or module are in execution process.

"""

printConfigHelp=\
"""Print configuration file
It will print location configuration file and its content.

"""


def main():
    # Construct command line arguments and parsed it
    parser, args = parseArguments()

    if args.cleanScratch:
        cleanScratch()

    if args.printConfig:
        printConfig()

def parseArguments():
    parser = argparse.ArgumentParser(
                prog='gcMapExplorer config',
                description=description,
                formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-cs', '--clean-scratch', action='store_true', dest='cleanScratch',
                        default=False, help=cleanScratchHelp)

    parser.add_argument('-pc', '--print-config', action='store_true', dest='printConfig',
                        default=False, help=printConfigHelp)

    idx = sys.argv.index("config")+1
    args = parser.parse_args(args=sys.argv[idx:])

    return parser, args

if __name__ == "__main__":
    main()
