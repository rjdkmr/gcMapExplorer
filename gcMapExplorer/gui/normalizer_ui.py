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

# Main Window Of Importer
pathToThisUI = os.path.join(PathToUIs, 'normalizer.ui')
Ui_ImporterWindow, ImporterWindowBase = loadUiType(pathToThisUI)
class ImporterWindow(ImporterWindowBase, Ui_ImporterWindow):
    def __init__(self):
        super(ImporterWindow, self).__init__()
        self.setupUi(self)

        # Hide tabbars from tab widget of specific options
        tabbars = self.specOptsTabWidget.findChildren(QTabBar)
        for tabbar in tabbars:
            tabbar.hide()

        # Resize hight and reduce size of log text box
        self.resize(self.width(), 680)
        self.splitter.setSizes([520, 160])

        # Remove maximize window buttons
        self.setWindowFlags( (self.windowFlags() | Qt.CustomizeWindowHint) & ~Qt.WindowMaximizeButtonHint)

        self.temporaryFiles = []
        self.command = None
        self.process = None

        self.setDefaultScratchDirs()
        self.connectButtons()
        self.connectLineEdits()

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

    def setDefaultScratchDirs(self):
        defaultDir = config['Dirs']['WorkingDirectory']
        self.genOptsScratchDirLineEdit.setText(defaultDir)

    def connectButtons(self):
        self.methodCBox.currentIndexChanged.connect( self.specOptsTabWidget.setCurrentIndex )
        self.whatsThisButton.clicked.connect( QWhatsThis.enterWhatsThisMode )

        self.inputFileButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.inputFileButton.clicked.connect( self.openInputFile )

        self.outputFileButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.outputFileButton.clicked.connect( self.browseOutputFile )

        self.genOptsScratchDirButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.genOptsScratchDirButton.clicked.connect( self.browseScratchDir)

        # Set icon on start and stop button
        self.executeStartButton.setIcon( self.style().standardIcon(QStyle.SP_MediaPlay) )
        self.executeStopButton.setIcon( self.style().standardIcon(QStyle.SP_MediaStop) )

        # Connect start and stop button
        self.executeStartButton.clicked.connect( self.runCommand )
        self.executeStopButton.clicked.connect( self.terminateProcessing )

        # Connect clear button
        self.logOutputClearButton.clicked.connect( self.logOutputPlainTextEdit.clear )

    def connectLineEdits(self):
        # Filters Line Edit - Check only double is accepted
        self.genOptsPtndLineEdit.setValidator(QDoubleValidator())
        self.genOptsTdoLineEdit.setValidator(QDoubleValidator())

        # KR tolerance and MSCM -- only numbers
        self.specOptsKrTolLineEdit.setValidator(QDoubleValidator())
        self.specOptsKrMscmLineEdit.setValidator(QIntValidator())

        # IC tolerance and Iteration -- only numbers
        self.specOptsIcTolLineEdit.setValidator(QDoubleValidator())
        self.specOptsIcIterLineEdit.setValidator(QIntValidator())

        # vmin and vmax
        self.vminLineEdit.setValidator(QDoubleValidator())
        self.vmaxLineEdit.setValidator(QDoubleValidator())

        # Check for input file - file exist; determine file format
        self.inputFileLineEdit.editingFinished.connect( lambda: guiHelpers.checkFileExist(self.inputFileLineEdit, self) )
        self.inputFileLineEdit.editingFinished.connect( lambda: self.setFileFormat(self.inputFileLineEdit, self.inputFormatCBox) )

        # Check for output file determine file format
        self.outFileLineEdit.editingFinished.connect( lambda: self.setFileFormat(self.outFileLineEdit, self.outFormatCBox) )

        # Filters line edit - limits the value that should be entered in these line edits
        self.genOptsPtndLineEdit.editingFinished.connect( lambda: guiHelpers.constrainValueInLineEdit(self, self.genOptsPtndLineEdit, 0, 100,))
        self.genOptsTdoLineEdit.editingFinished.connect( lambda: guiHelpers.constrainValueInLineEdit(self, self.genOptsTdoLineEdit, 0, 1,))

        # Check if scratch directory exist when user try to write own directory
        self.genOptsScratchDirLineEdit.editingFinished.connect( lambda: guiHelpers.checkDirExist(self.genOptsScratchDirLineEdit, self) )


    def openInputFile(self):
        """To get compressed file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " map files (*.ccmap *.gcmap *.hicmap);; All files (*.*)"
        path = QFileDialog.getOpenFileName(self, 'Open File', '', file_choices)
        if path[0]:
            self.inputFileLineEdit.setText(path[0])
            self.setFileFormat(self.inputFileLineEdit, self.inputFormatCBox)   # Set input file format

    def setFileFormat(self, lineEdit, comboBox):
        """ Set file format in input and output combo box

        It try to get file extension from input file name and set this
        format to combo box.

        """
        filename = lineEdit.text()
        if filename:
            ext = guiHelpers.getFileExtension(filename)
            if ext == '.ccmap' or ext == '.hicmap':
                comboBox.setCurrentIndex(1)
            elif ext == '.gcmap':
                comboBox.setCurrentIndex(2)
            else:
                comboBox.setCurrentIndex(0)
        else:
            comboBox.setCurrentIndex(0) # Set input file format to None

    def browseScratchDir(self):
        """Browse and choose scratch directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Scratch Directory')
        if path:
            self.genOptsScratchDirLineEdit.setText(path)

    def browseOutputFile(self):
        """To open gcmap or ccmap file with full path
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " gcmap or ccmap file (*.gcmap *.ccmap);;All files(*.*)"
        path = QFileDialog.getSaveFileName(self, 'Select or Create File', '', file_choices, options=QFileDialog.DontConfirmOverwrite)
        if path[0]:
            self.outFileLineEdit.setText(path[0])
            self.setFileFormat(self.outFileLineEdit, self.outFormatCBox)   # Set input file format

    def isInputOutputFormatValid(self):
        """ Check Input and Output formats selected by user.
        """
        valid = True
        if self.inputFormatCBox.currentText() == 'gcmap' and self.outFormatCBox.currentText() == 'ccmap':
            msg = 'When input format is gcmap, output format should be only gcmap.'
            guiHelpers.showWarningMessageBox(msg, self)
            self.outFormatCBox.setFocus()
            valid = False

        if self.inputFormatCBox.currentIndex() == 0:
            msg = 'Input file format not known!!'
            guiHelpers.showWarningMessageBox(msg, self)
            self.inputFormatCBox.setFocus()
            valid = False

        if self.outFormatCBox.currentIndex() == 0:
            msg = 'Output file format not known!!'
            guiHelpers.showWarningMessageBox(msg, self)
            self.outFormatCBox.setFocus()
            valid = False

        return valid

    def isPTNDandTDOvalid(self):
        """ Check if both filters are applied
        """
        ptndText = self.genOptsPtndLineEdit.text()
        tdoText = self.genOptsTdoLineEdit.text()

        if ptndText and tdoText:
            msg = 'Both Percentile and Fraction filters cannot be used simultaneously!!'
            guiHelpers.showWarningMessageBox(msg, self)
            self.genOptsPtndLineEdit.setFocus()
            return False

        return True

    def isKrToleranceValid(self):
        """ Check if KR tolerance is given properly.

        * if tolerance < 1e-12 : Value is too small
        * if no tolerance : Need it

        """
        value = self.specOptsKrTolLineEdit.text()

        valid = True
        if value:
            value = float(value)
            if value < 1e-12:
                msg = 'Tolerance value too small!!\nThis precision is difficult to achieve.'
                guiHelpers.showWarningMessageBox(msg, self)
                self.specOptsKrTolLineEdit.selectAll()
                self.specOptsKrTolLineEdit.setFocus()
                valid = False
        else:
            msg = 'Provide the Tolerance value.'
            guiHelpers.showWarningMessageBox(msg, self)
            self.specOptsKrTolLineEdit.setFocus()
            valid = False

        return valid

    def isKrMscmValueValid(self):
        """ Check Map Size Celing value for memory is given
        It is required when input file is in gcmap format otherwise ignore it.
        """
        value = self.specOptsKrMscmLineEdit.text()

        valid = True
        if self.inputFormatCBox.currentText() == 'gcmap' and not value:
            msg = 'Provide Map Size Celing value for memory!'
            guiHelpers.showWarningMessageBox(msg, self)
            self.specOptsKrMscmLineEdit.setFocus()
            valid = False

        return valid

    def isIcTolAndIterationValid(self):
        Tol = self.specOptsIcTolLineEdit.text()
        Iteration = self.specOptsIcIterLineEdit.text()

        valid = True
        if Tol:
            value = float(Tol)
            if value < 1e-12:
                msg = 'Tolerance value too small!!\nThis precision is difficult to achieve.'
                guiHelpers.showWarningMessageBox(msg, self)
                self.specOptsIcTolLineEdit.setFocus()
                valid = False

        if not Tol and not Iteration:
            msg = 'Neither Tolerance nor Iteration value found !!! \n'
            guiHelpers.showWarningMessageBox(msg, self)
            self.specOptsIcTolLineEdit.setFocus()
            valid = False

        if not Iteration:
            msg = 'Iteration value not found.\nSetting it to default value.'
            guiHelpers.showWarningMessageBox(msg, self)
            self.specOptsIcIterLineEdit.setText('500')

        if not Tol:
            msg = 'Tolerance value not found.\nSetting it to default value.'
            guiHelpers.showWarningMessageBox(msg, self)
            self.specOptsIcTolLineEdit.setText('1e-4')

        return valid

    def readAndConstructCommand(self):
        """ Read all options from interface and construct the command
        """

        # Name of command
        programs = ['normKR', 'normIC', 'normMCFS', 'normVC']
        program = programs[self.methodCBox.currentIndex()]

        cmdDict = dict()

        # Input File
        inputFile = self.inputFileLineEdit.text()
        if not inputFile:
            msg = 'No input file!!'
            guiHelpers.showWarningMessageBox(msg, self)
            self.inputFileLineEdit.setFocus()
            return False
        cmdDict['-i'] = os.path.normpath(inputFile)

        # Output File
        outputFile = self.outFileLineEdit.text()
        if not outputFile:
            msg = 'No output file!!'
            guiHelpers.showWarningMessageBox(msg, self)
            self.outFileLineEdit.setFocus()
            return False
        cmdDict['-o'] = os.path.normpath(outputFile)

        # Input and output format
        inputFormat = None
        outputFormat = None
        if self.isInputOutputFormatValid():
            inputFormat = self.inputFormatCBox.currentText()
            outputFormat = self.outFormatCBox.currentText()
        else:
            return False
        cmdDict['-fi'] = inputFormat
        cmdDict['-fo'] = outputFormat

        # Vmin and vmax
        if self.vminLineEdit.text():
            cmdDict['-vmin'] = self.vminLineEdit.text()
        if self.vmaxLineEdit.text():
            cmdDict['-vmax'] = self.vmaxLineEdit.text()

        # Working or scratch directory
        workDir = self.genOptsScratchDirLineEdit.text()
        if not workDir:
            msg = 'No output file!!'
            guiHelpers.showWarningMessageBox(msg, self)
            self.genOptsScratchDirLineEdit.setFocus()
            return False
        cmdDict['-wd'] = os.path.normpath(workDir)

        # Compression method for gcmap output file
        if self.outFormatCBox.currentText() == 'gcmap':
            cmdDict['-cmeth'] = str(self.outFileCompressionCBox.currentText()).lower()

        # Filters
        if self.isPTNDandTDOvalid():
            if self.genOptsPtndLineEdit.text():
                cmdDict['-ptnd'] = float(self.genOptsPtndLineEdit.text())

            if self.genOptsTdoLineEdit.text():
                cmdDict['-tdo'] = float(self.genOptsTdoLineEdit.text())
        else:
            return False

        # Options related to KR normalization
        if program == 'normKR':

            # Tolerance
            krTol = None
            if self.isKrToleranceValid():
                krTol = float(self.specOptsKrTolLineEdit.text())
            else:
                return False
            cmdDict['-t'] = krTol

            # map size ceiling for memory - only need when input file is gcmap
            krMscm = None
            if self.isKrMscmValueValid():
                if self.specOptsKrMscmLineEdit.text():
                    krMscm = int( self.specOptsKrMscmLineEdit.text() )
                    cmdDict['-mscm'] = krMscm
            else:
                return False

            # Memory - only need when input file is ccmap
            if self.specOptsKrMemCBox.currentIndex() == 0:
                krMemory = 'RAM'
            else:
                krMemory = 'HDD'
            cmdDict['-m'] = krMemory

            # Construct the command
            self.constructKrCommand(cmdDict)

        # Options related to Iterative Correction
        if program == 'normIC':

            icTol = None
            icIter = None
            if self.isIcTolAndIterationValid():
                icTol = float( self.specOptsIcTolLineEdit.text() )
                icIter = int( self.specOptsIcIterLineEdit.text() )
            else:
                return False

            cmdDict['-t'] = icTol
            cmdDict['-c'] = icIter

            # Construct the command
            self.constructIcCommand(cmdDict)

        # Options related to MCFS
        if program == 'normMCFS':
            cmdDict['-s'] = str( self.specOptsMcfsStatsCBox.currentText() ).lower()
            cmdDict['-st'] = str( self.specOptsMcfsSTypeCBox.currentText() ).lower()

            cidx = self.specOptsMcfsScaleInputCBox.currentIndex()
            if cidx == 1:
                cmdDict['-sui'] = '-sui'

            # Construct the command
            self.constructMcfsCommand(cmdDict)

        # Options related to normVC
        if program == 'normVC':
            cidx = self.specOptsVcovTabSqrtCBox.currentIndex()
            if cidx == 1:
                cmdDict['-sq'] = '-sq'

            # Construct the command
            self.constructVanillaCovCommand(cmdDict)

    def constructKrCommand(self, cmdDict):
        """ Construct normKI command
        """
        command = ' normKR '
        command += ' -i ' +  '"{0}"'.format(cmdDict['-i'])
        command += ' -fi ' + cmdDict['-fi']
        command += ' -o ' +  '"{0}"'.format(cmdDict['-o'])
        command += ' -fo ' + cmdDict['-fo']
        command += ' -t ' + str(cmdDict['-t'])
        command += ' -m ' + cmdDict['-m']

        if '-vmin' in cmdDict:
            command += ' -vmin ' + str(cmdDict['-vmin'])
        if '-vmax' in cmdDict:
            command += ' -vmax ' + str(cmdDict['-vmax'])
        if '-mscm' in cmdDict:
            command += ' -mscm ' + str(cmdDict['-mscm'])
        if '-cmeth' in cmdDict:
            command += ' -cmeth ' + cmdDict['-cmeth']
        command += ' -wd ' +  '"{0}"'.format(cmdDict['-wd'])
        if '-ptnd' in cmdDict:
            command += ' -ptnd ' + str(cmdDict['-ptnd'])
        if '-tdo' in cmdDict:
            command += ' -tdo ' + str(cmdDict['-tdo'])

        self.command = command

    def constructIcCommand(self, cmdDict):
        """ Construct normIC command
        """
        command = ' normIC '
        command += ' -i ' +  '"{0}"'.format(cmdDict['-i'])
        command += ' -fi ' + cmdDict['-fi']
        command += ' -o ' +  '"{0}"'.format(cmdDict['-o'])
        command += ' -fo ' + cmdDict['-fo']
        command += ' -t ' + str(cmdDict['-t'])
        command += ' -c ' + str(cmdDict['-c'])
        if '-vmin' in cmdDict:
            command += ' -vmin ' + str(cmdDict['-vmin'])
        if '-vmax' in cmdDict:
            command += ' -vmax ' + str(cmdDict['-vmax'])
        if '-cmeth' in cmdDict:
            command += ' -cmeth ' + cmdDict['-cmeth']
        command += ' -wd ' +  '"{0}"'.format(cmdDict['-wd'])
        if '-ptnd' in cmdDict:
            command += ' -ptnd ' + str(cmdDict['-ptnd'])
        if '-tdo' in cmdDict:
            command += ' -tdo ' + str(cmdDict['-tdo'])

        self.command = command

    def constructMcfsCommand(self, cmdDict):
        """ Construct normMCFS command
        """
        command = ' normMCFS '
        command += ' -i ' +  '"{0}"'.format(cmdDict['-i'])
        command += ' -fi ' + cmdDict['-fi']
        command += ' -o ' +  '"{0}"'.format(cmdDict['-o'])
        command += ' -fo ' + cmdDict['-fo']

        if '-vmin' in cmdDict:
            command += ' -vmin ' + str(cmdDict['-vmin'])
        if '-vmax' in cmdDict:
            command += ' -vmax ' + str(cmdDict['-vmax'])
        if '-cmeth' in cmdDict:
            command += ' -cmeth ' + cmdDict['-cmeth']
        command += ' -wd ' +  '"{0}"'.format(cmdDict['-wd'])
        if '-ptnd' in cmdDict:
            command += ' -ptnd ' + str(cmdDict['-ptnd'])
        if '-tdo' in cmdDict:
            command += ' -tdo ' + str(cmdDict['-tdo'])

        command += ' -s ' + str(cmdDict['-s'])
        command += ' -st ' + str(cmdDict['-st'])
        if '-su' in cmdDict:
            command += ' -su '

        self.command = command

    def constructVanillaCovCommand(self, cmdDict):
        """ Construct normVC command
        """
        command = ' normVC '
        command += ' -i ' +  '"{0}"'.format(cmdDict['-i'])
        command += ' -fi ' + cmdDict['-fi']
        command += ' -o ' +  '"{0}"'.format(cmdDict['-o'])
        command += ' -fo ' + cmdDict['-fo']

        if '-vmin' in cmdDict:
            command += ' -vmin ' + str(cmdDict['-vmin'])
        if '-vmax' in cmdDict:
            command += ' -vmax ' + str(cmdDict['-vmax'])
        if '-cmeth' in cmdDict:
            command += ' -cmeth ' + cmdDict['-cmeth']
        command += ' -wd ' +  '"{0}"'.format(cmdDict['-wd'])
        if '-ptnd' in cmdDict:
            command += ' -ptnd ' + str(cmdDict['-ptnd'])
        if '-tdo' in cmdDict:
            command += ' -tdo ' + str(cmdDict['-tdo'])

        if '-sq' in cmdDict:
            command += ' -sq '

        self.command = command

    def runCommand(self):
        self.command = None
        self.readAndConstructCommand()
        if self.command is None:
            return
        self.startProcess(self.command, self.executeStartButton)

    def startProcess(self, command, button):
        self.process = QProcess(self)
        self.process.start('gcMapExplorer', shlex.split(command))
        self.process.setProcessChannelMode( QProcess.MergedChannels )
        #self.process.setReadChannel( QProcess.StandardOutput )
        self.process.waitForStarted()

        self.process.readyReadStandardOutput.connect( self.writeLogOutputFromProcess )
        self.process.readyReadStandardError.connect( self.writeLogOutputFromProcess )
        self.process.finished.connect( lambda: self.finishedProcessing( button ) )

        self.methodCBox.setEnabled(False)
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
        self.methodCBox.setEnabled(True)
        self.process = None
