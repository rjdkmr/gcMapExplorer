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

import os, sys, re
import tempfile
import shlex, subprocess


from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.uic import loadUiType

from gcMapExplorer import importer
from gcMapExplorer.config import getConfig

# get configuration
config = getConfig()

from . import guiHelpers

# Determine absolute path to UIs directory. Relative path from this directory does not work.
DirToThisScript = os.path.dirname(os.path.abspath(__file__))
PathToUIs = os.path.join(DirToThisScript, 'UIs')

def main():
    app = QApplication(sys.argv)
    importer_interface = ImporterWindow()
    importer_interface.show()
    app.exec_()
    app.exit()

class cooMatFormatTabWidgetHelper:
    """ Helper class containing all member-functions related to COO matrix format tab-widget
    """

    def initCooMatFormatTabWidget(self):
        self.cooMatCommand = None

        self.cooMatIOTable.setColumnWidth(0, 350)
        self.cooMatIOTable.setColumnWidth(1, 100)
        self.cooMatIOTable.setColumnWidth(2, 100)

        self.cooMatArchivFileBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.cooMatIOAddFileButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.cooMatOutDirBrowseButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.cooMatScratchDirBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.cooMatGCMapOutSelectButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )

        self.cooMatTabRunButton.setIcon( self.style().standardIcon(QStyle.SP_MediaPlay) )
        self.cooMatTabStopButton.setIcon( self.style().standardIcon(QStyle.SP_MediaStop) )

        self.connectCooMatrixTabWidgets()

    def connectCooMatrixTabWidgets(self):
        self.cooMatIOTable.cellChanged.connect( self.cooMatInputFileCellChanged )
        self.cooMatIOAddFromFileButton.clicked.connect( self.cooMatAddInputFilesFromMetaFile )
        self.cooMatIOAddFileButton.clicked.connect( self.cooMatAddFileToSelectedRowCell )
        self.cooMatIOTableAddRowButton.clicked.connect( self.cooMatAddRowToTable )
        self.cooMatIOTableRemovRowButton.clicked.connect( self.cooMatRemoveRowFromTable )
        self.cooMatArchivFileBrowsButton.clicked.connect( self.cooMatBrowseCompressedFile )
        self.cooMatOutDirBrowseButton.clicked.connect( self.cooMatBrowseOutputDir )
        self.cooMatGCMapOutSelectButton.clicked.connect( self.cooMatOpenGCMapFile )
        self.cooMatScratchDirBrowsButton.clicked.connect( self.cooMatBrowseScratchDir )
        self.cooMatTabRunButton.clicked.connect( self.runCooMatrixCommand )
        self.cooMatTabStopButton.clicked.connect( self.terminateProcessing )

        self.cooMatArchivFileLineEdit.editingFinished.connect( lambda: checkFileExist(self.cooMatArchivFileLineEdit, self) )
        self.cooMatScratchDirLineEdit.editingFinished.connect( lambda: checkDirExist(self.cooMatScratchDirLineEdit, self) )
        self.cooMatOutDIrLineEdit.editingFinished.connect( lambda: checkDirExist(self.cooMatOutDIrLineEdit, self) )

    def cooMatSetLabelsToCells(self, row, xlabel=None, ylabel=None):
        """Convenient function to show xlabel and ylabel to table
        It is called several times below.
        """

        # Try to determine xlabel and ylabel from file name
        if xlabel is None:
            basename = os.path.basename( str( self.cooMatIOTable.item(row, 0).text() ) )
            g = re.match("Chr|chr", basename)
            if g:
                xlabel = re.split('_|-', basename[g.span()[0]:])[0]
                ylabel = xlabel
            else:
                xlabel = basename
                ylabel = xlabel

        if self.cooMatIOTable.item(row, 1) is None:
            self.cooMatIOTable.setItem( row, 1, QTableWidgetItem(0) )
        if self.cooMatIOTable.item(row, 2) is None:
            self.cooMatIOTable.setItem( row, 2, QTableWidgetItem(0) )

        self.cooMatIOTable.item(row, 1).setText(xlabel)
        self.cooMatIOTable.item(row, 2).setText(ylabel)

    def cooMatInputFileCellChanged(self, row, col):
        """When user write a input file name, automatically try to choose Labels
        """
        # Nothing to do when entered on second column
        if col != 0:    return
        self.cooMatSetLabelsToCells(row)

    def cooMatAddInputFilesFromMetaFile(self):
        """To built input file list and labels from a meta input file
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " Text file (*.txt *.*)"
        path = QFileDialog.getOpenFileName(self, 'Load File', '', file_choices)

        # Make list of input files
        inputFiles, xlabels, ylabels = [], [], []
        if path[0]:
            fin = open(path[0], 'r')

            # Reading meta input file
            for line in fin:
                line = line.rstrip().lstrip()
                if not line.strip():    continue
                temp = re.split('\s+', line)

                # Check if meta input file has only one column
                if len(temp) < 2:
                    msg = "{0} contains only one column [{1}].\n \
                        Need at least two column with contact map file \
                        and xlabel !!!\n".format(path[0], temp[0])
                    msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg, QMessageBox.Ok, self)
                    msgBox.exec_()
                    msgBox.close()

                # Check is meta input file contains a file path
                if not self.cooMatArchivFileLineEdit.text():
                    if not os.path.exists( temp[0] ):
                        msg = 'File "{0}" not found. \n\nIf it is inside the compressed file, first select it.'.format(temp[0])
                        msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg, QMessageBox.Ok, self)
                        msgBox.exec_()
                        msgBox.close()
                        return

                inputFiles.append(temp[0])
                xlabels.append(temp[1])

                # Determine whether ylabel is present
                if len(temp) >= 3:
                    ylabels.append(temp[2])
                else:
                    ylabels.append(temp[1])

        # If there is list, add this list to table
        if inputFiles:
            for fidx in range(len(inputFiles)):

                r = self.cooMatIOTable.rowCount()
                if fidx > r-1 :
                    self.cooMatIOTable.insertRow(fidx)

                if self.cooMatIOTable.item(fidx, 0) is None:
                    self.cooMatIOTable.setItem( fidx, 0, QTableWidgetItem(0) )

                self.cooMatIOTable.item(fidx, 0).setText(inputFiles[fidx])
                self.cooMatSetLabelsToCells(fidx, xlabels[fidx], ylabels[fidx])

    def cooMatAddFileToSelectedRowCell(self):
        """Add a input file to selected Input Files cell
        """

        row, col = getSelectedRowColumnFromTable(self.cooMatIOTable)

        # If no cell is selected in "Input Files row", raise a message
        if row is None or col != 0:
            msgBox = QMessageBox(self)
            msgBox.setWindowModality(Qt.WindowModal)
            msgBox.setWindowTitle('Information')
            msgBox.setIcon(QMessageBox.Information)
            msgBox.setText('         No Input Files cell selected. \n\nSelect a Input Files cell on table to add file.')
            msgBox.setStandardButtons(QMessageBox.Cancel)
            msgBox.exec_()
            msgBox.close()
            return

        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " Text file (*.txt *.*)"
        path = QFileDialog.getOpenFileName(self, 'Load File', '', file_choices)
        if path[0]:
            # Add input file to first column
            self.cooMatIOTable.item(row, 0).setText(path[0])
            self.cooMatSetLabelsToCells(row)

    def cooMatAddRowToTable(self):
        """Add a new row at the end of cooMatIOTable
        """
        self.cooMatIOTable.insertRow(self.cooMatIOTable.rowCount())

    def cooMatRemoveRowFromTable(self):
        """Remove a selected row from cooMatIOTable
        """
        # Total number of row
        row, col = getSelectedRowColumnFromTable(self.cooMatIOTable)

        # If no cell is selected in "Input Files row", raise a message
        if row is None:
            msgBox = QMessageBox(self)
            msgBox.setWindowModality(Qt.WindowModal)
            msgBox.setWindowTitle('Information')
            msgBox.setIcon(QMessageBox.Information)
            msgBox.setText('       A row is not selected \n\nSelect a Row to remove from table')
            msgBox.setStandardButtons(QMessageBox.Cancel)
            msgBox.exec_()
            msgBox.close()
            return

        self.cooMatIOTable.removeRow(row)

    def cooMatBrowseCompressedFile(self):
        """To get compressed file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " Archive file (*.zip *.tar.gz *.tar.bz2)"
        path = QFileDialog.getOpenFileName(self, 'Open File', '', file_choices)
        if path[0]:
            self.cooMatArchivFileLineEdit.setText(path[0])

    def cooMatBrowseOutputDir(self):
        """Browse and choose output directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Output Directory')
        if path:
            self.cooMatOutDIrLineEdit.setText(path)

    def cooMatBrowseScratchDir(self):
        """Browse and choose scratch directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Scratch Directory')
        if path:
            self.cooMatScratchDirLineEdit.setText(path)

    def cooMatOpenGCMapFile(self):
        """To open gcmap file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " gcmap file (*.gcmap);;All files(*.*)"
        path = QFileDialog.getSaveFileName(self, 'Select or Create File', '', file_choices, options=QFileDialog.DontConfirmOverwrite)
        if path[0]:
            self.cooMatGCMapOutLineEdit.setText(path[0])

    def readAndConstructCooMatCommand(self):
        """Read and Construct the command line
        """
        self.cooMatCommand = None
        options = dict()

        if self.cooMatArchivFileLineEdit.text():
            inputCompressedFile = str( self.cooMatArchivFileLineEdit.text() )
        else:
            inputCompressedFile = 'None'

        options['inputCompressedFile'] = inputCompressedFile

        inputFiles, xlabels, ylabels = [], [], []
        rows = self.cooMatIOTable.rowCount()
        for row in range(rows):
            if self.cooMatIOTable.item(row, 0) is not None and self.cooMatIOTable.item(row, 0).text() :
                inputFiles.append( str( self.cooMatIOTable.item(row, 0).text() ) )
                xlabels.append( str( self.cooMatIOTable.item(row, 1).text() ) )
                ylabels.append( str( self.cooMatIOTable.item(row, 2).text() ) )

        if inputCompressedFile == 'None' and not inputFiles:
            msg = "Neither input compressed file, nor Input Files are given!!!"
            showWarningMessageBox(msg, self)
            return False

        options['inputFiles'] = inputFiles
        options['xlabels'] = xlabels
        options['ylabels'] = ylabels

        options['mapType'] = str( self.cooMatMapTypeCBox.currentText() )
        if self.cooMatResolutionLineEdit.text():
            options['resolution'] = self.cooMatResolutionLineEdit.text()
        else:
            options['resolution'] = 'None'

        options['coordinate'] = str( self.cooMatCoordTypeCBox.currentText() )
        options['workDir'] = os.path.normpath(self.cooMatScratchDirLineEdit.text())

        if not self.cooMatCCMapGroupBox.isChecked() \
                        and not self.cooMatGCMapBoxGroupBox.isChecked():
            showWarningMessageBox("No ccmap or gcmap output !!!", self)
            return False

        options['ccmap'] = False
        if self.cooMatCCMapGroupBox.isChecked():
            ccmapSuffix = str( self.cooMatCCMapOutSuffixLineEdit.text() )
            if not ccmapSuffix:
                msg = "No suffix provided for ccmap files. \n" \
                        + "Please provide a suffix for proper name."
                showWarningMessageBox(msg, self)
                self.cooMatCCMapOutSuffixLineEdit.setFocus()
                return False

            outDir = str( self.cooMatOutDIrLineEdit.text() )
            if not outDir:
                msg = "No Output Directory provided for ccmap files \n" \
                        + "Please select a directory to save ccmap files."
                showWarningMessageBox(msg, self)
                self.cooMatOutDIrLineEdit.setFocus()
                return False

            options['ccmap'] = True
            options['ccmapSuffix'] = ccmapSuffix
            options['outDir'] = '"{0}"'.format(outDir)

        options['gcmap'] = False
        if self.cooMatGCMapBoxGroupBox.isChecked():
            fileGCMap = str( self.cooMatGCMapOutLineEdit.text() )
            if not fileGCMap:
                msg = "No Output gcmap file is created or selected \n" \
                        + "Please select or create a gcmap file."
                showWarningMessageBox(msg, self)
                self.cooMatGCMapOutLineEdit.setFocus()
                return False

            options['gcmap'] = True
            options['fileGCMap'] = '"{0}"'.format(fileGCMap)
            options['compression'] = str( self.cooMatGCMapCompressCBox.currentText() ).lower()
            options['coarseningMethod'] = str( self.cooMatGCMapDownsampleCBox.currentText() ).lower()

        self.cooMatrixConstructCommand(options)

    def cooMatrixConstructCommand(self, opts):
        """Construct the command line
        """

        command = ' coo2cmap '

        if opts['inputFiles']:
            of, inputMetaFileName = tempfile.mkstemp(suffix='.temp', prefix='gcx_', dir=opts['workDir'], text=True)
            os.close(of)

            fout = open(inputMetaFileName, 'w')
            for i in range(len(opts['inputFiles'])):
                fout.write('{0}\t{1}\t{2}\n'.format(opts['inputFiles'][i],
                                                    opts['xlabels'][i],
                                                    opts['ylabels'][i]))
            fout.close()

            command += ' -i ' + '"{0}"'.format(os.path.normpath(inputMetaFileName))
            self.temporaryFiles.append(inputMetaFileName)

        if opts['inputCompressedFile'] != 'None':
            command += ' -ic ' +  '"{0}"'.format(os.path.normpath(opts['inputCompressedFile']))

        if opts['resolution'] != 'None':
            command += ' -r ' + opts['resolution']

        if opts['coordinate'] == 'index':
            command += ' -idx '

        if opts['ccmap']:
            command += ' -ccm ' + opts['ccmapSuffix']
            command += ' -od ' + opts['outDir']

        if opts['gcmap']:
            command += ' -gcm ' + opts['fileGCMap']
            command += ' -cmeth ' + opts['compression']
            command += ' -dmeth ' + opts['coarseningMethod']

        command += ' -mt ' + opts['mapType']
        command += ' -wd ' + '"{0}"'.format(opts['workDir'])

        self.cooMatCommand = command

class homerFormatTabWidgetHelper:
    """ Helper class containing all member-functions related to HOMER matrix format tab-widget
    """

    def initHomerFormatTabWidget(self):
        self.homerCommand = None

        self.homerInputFileBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.homerOutDirBrowseButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.homerScratchDirBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.homerGCMapOutSelectButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )

        self.homerTabRunButton.setIcon( self.style().standardIcon(QStyle.SP_MediaPlay) )
        self.homerTabStopButton.setIcon( self.style().standardIcon(QStyle.SP_MediaStop) )

        self.connectHomerTabWidgets()

    def connectHomerTabWidgets(self):
        self.homerInputFileBrowsButton.clicked.connect( self.homerBrowseInputFile )
        self.homerScratchDirBrowsButton.clicked.connect( self.homerBrowseScratchDir )
        self.homerOutDirBrowseButton.clicked.connect( self.homerBrowseOutputDir )
        self.homerGCMapOutSelectButton.clicked.connect( self.homerOpenGCMapFile )

        self.homerTabRunButton.clicked.connect( self.runHomerCommand )
        self.homerTabStopButton.clicked.connect( self.terminateProcessing )

        self.homerInputFileLineEdit.editingFinished.connect( lambda: checkFileExist(self.homerInputFileLineEdit, self) )
        self.homerScratchDirLineEdit.editingFinished.connect( lambda: checkDirExist(self.homerScratchDirLineEdit, self) )
        self.homerOutDIrLineEdit.editingFinished.connect( lambda: checkDirExist(self.homerOutDIrLineEdit, self) )


    def homerBrowseInputFile(self):
        """To get compressed file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " Text file (*.txt *.dat);; All file (*.*)"
        path = QFileDialog.getOpenFileName(self, 'Open File', '', file_choices)
        if path[0]:
            self.homerInputFileLineEdit.setText(path[0])

    def homerBrowseOutputDir(self):
        """Browse and choose output directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Output Directory')
        if path:
            self.homerOutDIrLineEdit.setText(path)

    def homerBrowseScratchDir(self):
        """Browse and choose scratch directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Scratch Directory')
        if path:
            self.homerScratchDirLineEdit.setText(path)

    def homerOpenGCMapFile(self):
        """To open gcmap file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " gcmap file (*.gcmap);;All files(*.*)"
        path = QFileDialog.getSaveFileName(self, 'Select or Create File', '', file_choices, options=QFileDialog.DontConfirmOverwrite)
        if path[0]:
            self.homerGCMapOutLineEdit.setText(path[0])

    def readAndConstructHomerCommand(self):
        """Read and construct the command line
        """
        self.homerCommand = None
        options = dict()

        inputFile = None
        if self.homerInputFileLineEdit.text():
            inputFile = str( self.homerInputFileLineEdit.text() )
        else:
            self.homerInputFileLineEdit.setFocus()
            showWarningMessageBox("No input file given !!!", self)
            return False

        options['inputFile'] = os.path.normpath(inputFile)

        options['workDir'] = os.path.normpath(self.homerScratchDirLineEdit.text())

        if not self.homerCCMapGroupBox.isChecked() \
                        and not self.homerGCMapGroupBox.isChecked():
            showWarningMessageBox("No ccmap or gcmap output !!!", self)
            return False

        options['ccmap'] = False
        if self.homerCCMapGroupBox.isChecked():
            ccmapSuffix = str( self.homerOutSuffixLineEdit.text() )
            if not ccmapSuffix:
                msg = "No suffix provided for ccmap files. \n" \
                        + "Please provide a suffix for proper name."
                showWarningMessageBox(msg, self)
                self.homerOutSuffixLineEdit.setFocus()
                return False

            outDir = str( self.homerOutDIrLineEdit.text() )
            if not outDir:
                msg = "No Output Directory provided for ccmap files \n" \
                        + "Please select a directory to save ccmap files."
                showWarningMessageBox(msg, self)
                self.homerOutDIrLineEdit.setFocus()
                return False

            options['ccmap'] = True
            options['ccmapSuffix'] = ccmapSuffix
            options['outDir'] = os.path.normpath(outDir)

        options['gcmap'] = False
        if self.homerGCMapGroupBox.isChecked():
            fileGCMap = str( self.homerGCMapOutLineEdit.text() )
            if not fileGCMap:
                msg = "No Output gcmap file is created or selected \n" \
                        + "Please select or create a gcmap file."
                showWarningMessageBox(msg, self)
                self.homerGCMapOutLineEdit.setFocus()
                return False

            options['gcmap'] = True
            options['fileGCMap'] = os.path.normpath(fileGCMap)
            options['compression'] = str( self.homerGCMapCompressCBox.currentText() ).lower()
            options['coarseningMethod'] = str( self.homerGCMapDownsampleCBox.currentText() ).lower()

        self.homerConstructCommand(options)

    def homerConstructCommand(self, opts):
        """Construct the command line
        """
        command = ' homer2cmap '

        command += ' -i ' +  '"{0}"'.format(opts['inputFile'])

        if opts['ccmap']:
            command += ' -ccm ' + opts['ccmapSuffix']
            command += ' -od ' + '"{0}"'.format(opts['outDir'])

        if opts['gcmap']:
            command += ' -gcm ' + '"{0}"'.format(opts['fileGCMap'])
            command += ' -cmeth ' + opts['compression']
            command += ' -dmeth ' + opts['coarseningMethod']

        command += ' -wd ' + '"{0}"'.format(opts['workDir'])

        self.homerCommand = command

class binContactFormatTabWidgetHelper:
    """ Helper class containing all member-functions related to Bin-Contact files pair tab-widget
    """

    def initBinContactFormatTabWidget(self):
        self.binContactCommand = None

        self.binContactInputBinFileBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.binContactInputContactFileBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.binContactOutDirBrowseButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.binContactScratchDirBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.binContactGCMapOutSelectButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )

        self.binContactTabRunButton.setIcon( self.style().standardIcon(QStyle.SP_MediaPlay) )
        self.binContactTabStopButton.setIcon( self.style().standardIcon(QStyle.SP_MediaStop) )

        self.connectBinContactTabWidgets()

    def connectBinContactTabWidgets(self):
        self.binContactInputBinFileBrowsButton.clicked.connect( self.binContactBrowseInputBinFile )
        self.binContactInputContactFileBrowsButton.clicked.connect( self.binContactBrowseInputContactFile )
        self.binContactScratchDirBrowsButton.clicked.connect( self.binContactBrowseScratchDir )
        self.binContactOutDirBrowseButton.clicked.connect( self.binContactBrowseOutputDir )
        self.binContactGCMapOutSelectButton.clicked.connect( self.binContactOpenGCMapFile )

        self.binContactTabRunButton.clicked.connect( self.runBinContactCommand )
        self.binContactTabStopButton.clicked.connect( self.terminateProcessing )

        self.binContactInputBinFileLineEdit.editingFinished.connect( lambda: checkFileExist(self.binContactInputBinFileLineEdit, self) )
        self.binContactInputContactFileLineEdit.editingFinished.connect( lambda: checkFileExist(self.binContactInputContactFileLineEdit, self) )
        self.binContactScratchDirLineEdit.editingFinished.connect( lambda: checkDirExist(self.binContactScratchDirLineEdit, self) )
        self.binContactOutDIrLineEdit.editingFinished.connect( lambda: checkDirExist(self.binContactOutDIrLineEdit, self) )


    def binContactBrowseInputBinFile(self):
        """To get compressed file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " Text file (*.txt *.dat);; All file (*.*)"
        path = QFileDialog.getOpenFileName(self, 'Open File', '', file_choices)
        if path[0]:
            self.binContactInputBinFileLineEdit.setText(path[0])

    def binContactBrowseInputContactFile(self):
        """To get compressed file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " Text file (*.txt *.dat);; All file (*.*)"
        path = QFileDialog.getOpenFileName(self, 'Open File', '', file_choices)
        if path[0]:
            self.binContactInputContactFileLineEdit.setText(path[0])

    def binContactBrowseOutputDir(self):
        """Browse and choose output directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Output Directory')
        if path:
            self.binContactOutDIrLineEdit.setText(path)

    def binContactBrowseScratchDir(self):
        """Browse and choose scratch directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Scratch Directory')
        if path:
            self.binContactScratchDirLineEdit.setText(path)

    def binContactOpenGCMapFile(self):
        """To open gcmap file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " gcmap file (*.gcmap);;All files(*.*)"
        path = QFileDialog.getSaveFileName(self, 'Select or Create File', '', file_choices, options=QFileDialog.DontConfirmOverwrite)
        if path[0]:
            self.binContactGCMapOutLineEdit.setText(path[0])

    def readAndConstructBinContactCommand(self):
        """Read and construct the command line
        """
        self.binContactCommand = None
        options = dict()

        inputBinFile = None
        if self.binContactInputBinFileLineEdit.text():
            inputFile = str( self.binContactInputBinFileLineEdit.text() )
        else:
            self.binContactInputBinFileLineEdit.setFocus()
            showWarningMessageBox("No input file given !!!", self)
            return False
        options['inputBinFile'] = os.path.normpath(inputFile)

        inputContactFile = None
        if self.binContactInputContactFileLineEdit.text():
            inputFile = str( self.binContactInputContactFileLineEdit.text() )
        else:
            self.binContactInputContactFileLineEdit.setFocus()
            showWarningMessageBox("No input file given !!!", self)
            return False
        options['inputContactFile'] = os.path.normpath(inputFile)

        options['workDir'] = os.path.normpath(self.binContactScratchDirLineEdit.text())

        if not self.binContactCCMapGroupBox.isChecked() \
                        and not self.binContactGCMapGroupBox.isChecked():
            showWarningMessageBox("No ccmap or gcmap output !!!", self)
            return False

        options['ccmap'] = False
        if self.binContactCCMapGroupBox.isChecked():
            ccmapSuffix = str( self.binContactOutSuffixLineEdit.text() )
            if not ccmapSuffix:
                msg = "No suffix provided for ccmap files. \n" \
                        + "Please provide a suffix for proper name."
                showWarningMessageBox(msg, self)
                self.binContactOutSuffixLineEdit.setFocus()
                return False

            outDir = str( self.binContactOutDIrLineEdit.text() )
            if not outDir:
                msg = "No Output Directory provided for ccmap files \n" \
                        + "Please select a directory to save ccmap files."
                showWarningMessageBox(msg, self)
                self.binContactOutDIrLineEdit.setFocus()
                return False

            options['ccmap'] = True
            options['ccmapSuffix'] = ccmapSuffix
            options['outDir'] = os.path.normpath(outDir)

        options['gcmap'] = False
        if self.binContactGCMapGroupBox.isChecked():
            fileGCMap = str( self.binContactGCMapOutLineEdit.text() )
            if not fileGCMap:
                msg = "No Output gcmap file is created or selected \n" \
                        + "Please select or create a gcmap file."
                showWarningMessageBox(msg, self)
                self.binContactGCMapOutLineEdit.setFocus()
                return False

            options['gcmap'] = True
            options['fileGCMap'] = os.path.normpath(fileGCMap)
            options['compression'] = str( self.binContactGCMapCompressCBox.currentText() ).lower()
            options['coarseningMethod'] = str( self.binContactGCMapDownsampleCBox.currentText() ).lower()

        self.binContactConstructCommand(options)

    def binContactConstructCommand(self, opts):
        """Construct the command line
        """
        command = ' bc2cmap '

        command += ' -ib ' +  '"{0}"'.format(opts['inputBinFile'])
        command += ' -ic ' +  '"{0}"'.format(opts['inputContactFile'])

        if opts['ccmap']:
            command += ' -ccm ' + opts['ccmapSuffix']
            command += ' -od ' + '"{0}"'.format(opts['outDir'])

        if opts['gcmap']:
            command += ' -gcm ' + '"{0}"'.format(opts['fileGCMap'])
            command += ' -cmeth ' + opts['compression']
            command += ' -dmeth ' + opts['coarseningMethod']

        command += ' -wd ' + '"{0}"'.format(opts['workDir'])

        self.binContactCommand = command

class pairCooMatFormatTabWidgetHelper:
    """ Helper class containing all member-functions related to paired-COO matrix file tab-widget
    """

    def initPairCooMatFormatTabWidget(self):
        self.pairCooMatCommand = None

        self.pairCooMatInputFileBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.pairCooMatOutDirBrowseButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.pairCooMatScratchDirBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.pairCooMatGCMapOutSelectButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )

        self.pairCooMatTabRunButton.setIcon( self.style().standardIcon(QStyle.SP_MediaPlay) )
        self.pairCooMatTabStopButton.setIcon( self.style().standardIcon(QStyle.SP_MediaStop) )

        self.connectPairCooMatTabWidgets()

    def connectPairCooMatTabWidgets(self):
        self.pairCooMatInputFileBrowsButton.clicked.connect( self.pairCooMatBrowseInputFile )
        self.pairCooMatScratchDirBrowsButton.clicked.connect( self.pairCooMatBrowseScratchDir )
        self.pairCooMatOutDirBrowseButton.clicked.connect( self.pairCooMatBrowseOutputDir )
        self.pairCooMatGCMapOutSelectButton.clicked.connect( self.pairCooMatOpenGCMapFile )

        self.pairCooMatTabRunButton.clicked.connect( self.runPairCooMatCommand )
        self.pairCooMatTabStopButton.clicked.connect( self.terminateProcessing )

        self.pairCooMatInputFileLineEdit.editingFinished.connect( lambda: checkFileExist(self.pairCooMatInputFileLineEdit, self) )
        self.pairCooMatScratchDirLineEdit.editingFinished.connect( lambda: checkDirExist(self.pairCooMatScratchDirLineEdit, self) )
        self.pairCooMatOutDIrLineEdit.editingFinished.connect( lambda: checkDirExist(self.pairCooMatOutDIrLineEdit, self) )


    def pairCooMatBrowseInputFile(self):
        """To get compressed file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " Text file (*.txt *.dat);; All file (*.*)"
        path = QFileDialog.getOpenFileName(self, 'Open File', '', file_choices)
        if path[0]:
            self.pairCooMatInputFileLineEdit.setText(path[0])

    def pairCooMatBrowseOutputDir(self):
        """Browse and choose output directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Output Directory')
        if path:
            self.pairCooMatOutDIrLineEdit.setText(path)

    def pairCooMatBrowseScratchDir(self):
        """Browse and choose scratch directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Scratch Directory')
        if path:
            self.pairCooMatScratchDirLineEdit.setText(path)

    def pairCooMatOpenGCMapFile(self):
        """To open gcmap file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " gcmap file (*.gcmap);;All files(*.*)"
        path = QFileDialog.getSaveFileName(self, 'Select or Create File', '', file_choices, options=QFileDialog.DontConfirmOverwrite)
        if path[0]:
            self.pairCooMatGCMapOutLineEdit.setText(path[0])

    def readAndConstructPairCooMatCommand(self):
        """Read and construct the command line
        """
        self.pairCooMatCommand = None
        options = dict()

        inputFile = None
        if self.pairCooMatInputFileLineEdit.text():
            inputFile = str( self.pairCooMatInputFileLineEdit.text() )
        else:
            self.pairCooMatInputFileLineEdit.setFocus()
            showWarningMessageBox("No input file given !!!", self)
            return False
        options['-i'] = os.path.normpath(inputFile)

        options['-wd'] = os.path.normpath(self.pairCooMatScratchDirLineEdit.text())

        if not self.pairCooMatCCMapGroupBox.isChecked() \
                        and not self.pairCooMatGCMapGroupBox.isChecked():
            showWarningMessageBox("No ccmap or gcmap output !!!", self)
            return False

        ccmapSuffix = None
        if self.pairCooMatCCMapGroupBox.isChecked():
            ccmapSuffix = str( self.pairCooMatOutSuffixLineEdit.text() )
            if not ccmapSuffix:
                msg = "No suffix provided for ccmap files. \n" \
                        + "Please provide a suffix for proper name."
                showWarningMessageBox(msg, self)
                self.pairCooMatOutSuffixLineEdit.setFocus()
                return False

            outDir = str( self.pairCooMatOutDIrLineEdit.text() )
            if not outDir:
                msg = "No Output Directory provided for ccmap files \n" \
                        + "Please select a directory to save ccmap files."
                showWarningMessageBox(msg, self)
                self.pairCooMatOutDIrLineEdit.setFocus()
                return False

            options['-ccm'] = ccmapSuffix
            options['-od'] = os.path.normpath(outDir)

        if self.pairCooMatGCMapGroupBox.isChecked():
            fileGCMap = str( self.pairCooMatGCMapOutLineEdit.text() )
            if not fileGCMap:
                msg = "No Output gcmap file is created or selected \n" \
                        + "Please select or create a gcmap file."
                showWarningMessageBox(msg, self)
                self.pairCooMatGCMapOutLineEdit.setFocus()
                return False

            options['-gcm'] = os.path.normpath(fileGCMap)
            options['-cmeth'] = str( self.pairCooMatGCMapCompressCBox.currentText() ).lower()
            options['-dmeth'] = str( self.pairCooMatGCMapDownsampleCBox.currentText() ).lower()

        self.pairCooMatConstructCommand(options)

    def pairCooMatConstructCommand(self, opts):
        """Construct the command line
        """
        command = ' pairCoo2cmap '
        command += ' -i ' +  '"{0}"'.format(opts['-i'])

        if '-ccm' in opts:
            command += ' -ccm ' + opts['-ccm']
            command += ' -od ' + '"{0}"'.format(opts['-od'])

        if '-gcm' in opts:
            command += ' -gcm ' + '"{0}"'.format(opts['-gcm'])
            command += ' -cmeth ' + opts['-cmeth']
            command += ' -dmeth ' + opts['-dmeth']

        command += ' -wd ' + '"{0}"'.format(opts['-wd'])

        self.pairCooMatCommand = command

# Main Window Of Importer
pathToThisUI = os.path.join(PathToUIs, 'importer.ui')
Ui_ImporterWindow, ImporterWindowBase = loadUiType(pathToThisUI)
class ImporterWindow(ImporterWindowBase, Ui_ImporterWindow, cooMatFormatTabWidgetHelper,
                     homerFormatTabWidgetHelper, binContactFormatTabWidgetHelper,
                     pairCooMatFormatTabWidgetHelper):
    def __init__(self):
        super(ImporterWindow, self).__init__()
        self.setupUi(self)

        # Hide tabbars
        tabbars = self.findChildren(QTabBar)
        for tabbar in tabbars:
            tabbar.hide()

        # Resize height and reduce size of log text box
        self.resize(self.width(), 680)
        self.splitter.setSizes([500, 180])

        self.temporaryFiles = []
        self.process = None

        self.initCooMatFormatTabWidget()
        self.initHomerFormatTabWidget()
        self.initBinContactFormatTabWidget()
        self.initPairCooMatFormatTabWidget()

        self.setDefaultScratchDirs()
        self.connectMainButtons()

    def closeEvent(self, event):

        if self.process is not None:
            if self.process.state() == QProcess.Running:
                msg = " A process is still Running... \n" \
                        + "Are you sure to close?"
                msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg, QMessageBox.Yes | QMessageBox.No, self)
                msgBox.exec_()

                if msgBox.result() == QMessageBox.No:
                    close = False
                    event.ignore()
                else:
                    self.terminateProcessing()

                msgBox.close()

        if self.temporaryFiles:
            for i in range(len(self.temporaryFiles)):
                if os.path.isfile( self.temporaryFiles[i] ):
                    os.remove( self.temporaryFiles[i] )

    def connectMainButtons(self):
        self.inputSelectorQCBox.currentIndexChanged.connect( self.InputsTabWidget.setCurrentIndex )
        self.whatsThisButton.clicked.connect( QWhatsThis.enterWhatsThisMode )
        self.logOutputClearButton.clicked.connect( self.logOutputPlainTextEdit.clear )

    def setDefaultScratchDirs(self):
        defaultDir = config['Dirs']['WorkingDirectory']
        self.cooMatScratchDirLineEdit.setText(defaultDir)
        self.homerScratchDirLineEdit.setText(defaultDir)
        self.binContactScratchDirLineEdit.setText(defaultDir)
        self.pairCooMatScratchDirLineEdit.setText(defaultDir)

    def runCooMatrixCommand(self):
        self.readAndConstructCooMatCommand()
        if self.cooMatCommand is None:  return
        self.startProcess(self.cooMatCommand, self.cooMatTabRunButton)

    def runHomerCommand(self):
        self.readAndConstructHomerCommand()
        if self.homerCommand is None:  return
        self.startProcess(self.homerCommand, self.homerTabRunButton)

    def runBinContactCommand(self):
        self.readAndConstructBinContactCommand()
        if self.binContactCommand is None:  return
        self.startProcess(self.binContactCommand, self.binContactTabRunButton)

    def runPairCooMatCommand(self):
        self.readAndConstructPairCooMatCommand()
        if self.pairCooMatCommand is None:  return
        self.startProcess(self.pairCooMatCommand, self.pairCooMatTabRunButton)

    def startProcess(self, command, button):
        self.process = QProcess(self)
        self.process.start('gcMapExplorer', shlex.split(command))
        self.process.setProcessChannelMode( QProcess.MergedChannels)
        # self.process.setReadChannel( QProcess.StandardOutput )
        self.process.waitForStarted()

        self.process.readyReadStandardOutput.connect( self.writeLogOutputFromProcess )
        self.process.readyReadStandardError.connect( self.writeLogOutputFromProcess )
        self.process.finished.connect( lambda: self.finishedProcessing( button ) )

        self.inputSelectorQCBox.setEnabled(False)
        button.setEnabled(False)
        self.logOutputPlainTextEdit.clear()
        self.logOutputPlainTextEdit.appendPlainText(
            '##### Started Running #### \n gcMapExplorer ' + command + "\n ##### ####### ####")


    def writeLogOutputFromProcess(self):
        out = bytes(self.process.readAllStandardOutput()).decode("utf-8")
        if out:
            out = out.rstrip()
            self.logOutputPlainTextEdit.appendPlainText( out  )
        out = bytes(self.process.readAllStandardError()).decode("utf-8")
        if out:
            out = out.rstrip()
            self.logOutputPlainTextEdit.appendPlainText( out  )

    def terminateProcessing(self):
        if self.process is None: return
        self.process.terminate()

    def finishedProcessing(self, button):
        button.setEnabled(True)
        self.inputSelectorQCBox.setEnabled(True)
        self.process = None


def getSelectedRowColumnFromTable(table):
    # Total number of row
    r = table.rowCount()
    c = table.columnCount()

    # Determine which cell of row is selected
    selectedRow = None
    selectedCol = None
    for i in range(r):
        for j in range(c):
            if table.item(i, j) is None:
                table.setItem( i, j, QTableWidgetItem(0) )

            if table.item(i, j).isSelected():
                selectedRow = i
                selectedCol = j
                break

    return selectedRow, selectedCol


def checkFileExist(lineEdit, qwidget):
    if not lineEdit.text(): return
    filename = str( lineEdit.text() )
    if not os.path.isfile( filename ):
        msg = "[ {0} ] \n Not found !!!".format(filename)
        msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg,
                                                    QMessageBox.Ok, qwidget)
        msgBox.exec_()
        msgBox.close()

        lineEdit.selectAll()
        lineEdit.setFocus()

def checkDirExist(lineEdit, qwidget):
    if not lineEdit.text(): return
    dirname = str( lineEdit.text() )
    if not os.path.isdir( dirname ):
        msg = "[ {0} ] \n Not found !!!".format(dirname)
        msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg,
                                                    QMessageBox.Ok, qwidget)
        msgBox.exec_()
        msgBox.close()

        lineEdit.selectAll()
        lineEdit.setFocus()

def showWarningMessageBox(msg, qwidget):
    msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg, QMessageBox.Ok, qwidget)
    msgBox.exec_()
    msgBox.close()
