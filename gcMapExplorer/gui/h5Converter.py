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

import os, re
import shlex
import sys
import string, random

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtPrintSupport import *
from PyQt5.uic import loadUiType

from gcMapExplorer import lib as gmlib

from gcMapExplorer.config import getConfig

# get configuration
config = getConfig()

# Determine absolute path to UIs directory. Relative path from this directory does not work.
DirToThisScript = os.path.dirname(os.path.abspath(__file__))
PathToUIs = os.path.join(DirToThisScript, 'UIs')

def main():
    app = QApplication(sys.argv)
    converter_interface = DialogOther1DFormatLoader()
    converter_interface.setAttribute(Qt.WA_DeleteOnClose)
    converter_interface.show()
    converter_interface.removeH5RemoveRadioButton()
    app.exec_()
    app.exit()

# Dialog box to change page size by custom length
pathToThisUI = os.path.join(PathToUIs, 'other1DFormatFileConversion.ui')
Ui_DialogOther1DFormatLoader, DialogOther1DFormatLoaderBase = loadUiType(pathToThisUI)
class DialogOther1DFormatLoader(DialogOther1DFormatLoaderBase, Ui_DialogOther1DFormatLoader):
    def __init__(self, parent=None):
        super(DialogOther1DFormatLoader, self).__init__(parent=parent)
        self.setupUi(self)

        # Resize height and reduce size of log text box
        self.resize(self.width(), 600)
        self.splitter.setSizes([460, 140])
        self.command = None
        self.process = None
        self.hideMode = False
        self.indexFile = None
        self.hasConvertedSuccessfully = False


        # Hide tabbars
        tabbars = self.fileFormatSpecificTabWidget.findChildren(QTabBar)
        for tabbar in tabbars:
            tabbar.hide()

        # Set default scratch directory, bigWigInfo and bigWigToWig paths
        self.setDefaultScratchDirs()
        self.setBigWigInfoPathConfig()
        self.setBigWigToWigPathConfig()


        self.inFileButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.bigWigToWigPathButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.bigWigInfoButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.h5FileButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.scratchDirButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )

        self.ameanCheckBox.clicked.connect( self.checkDownsamplingMethods )
        self.hmeanCheckBox.clicked.connect( self.checkDownsamplingMethods )
        self.gmeanCheckBox.clicked.connect( self.checkDownsamplingMethods )
        self.minCheckBox.clicked.connect( self.checkDownsamplingMethods )
        self.maxCheckBox.clicked.connect( self.checkDownsamplingMethods )
        self.medianCheckBox.clicked.connect( self.checkDownsamplingMethods )

        self.inFileButton.clicked.connect( self.browseInputFile )
        self.fileFormatCBox.currentIndexChanged.connect( self.fileFormatSpecificTabWidget.setCurrentIndex )
        self.scratchDirButton.clicked.connect( self.browseScratchDir )
        self.h5FileButton.clicked.connect( self.browseH5File )
        if hasattr(self, 'h5RemoveRadioButton'):
            self.h5RemoveRadioButton.clicked.connect( self.isValidH5OutFileFromWidget )
        self.bigWigInfoButton.clicked.connect( self.browseBigWigInfoProgram )
        self.bigWigToWigPathButton.clicked.connect( self.browseBigWigToWigProgram )
        self.logClearButton.clicked.connect( self.logOutputPlainTextEdit.clear )

        self.inFileLineEdit.editingFinished.connect( self.checkInputFileFormat )
        self.h5FileLineEdit.editingFinished.connect( self.isValidH5OutFileFromWidget )
        self.scratchDirLineEdit.editingFinished.connect( self.checkScratchDirectory )
        self.resolutionLineEdit.editingFinished.connect( self.checkResolutions )
        self.bigWigInfoLineEdit.editingFinished.connect( self.checkBigWigInfoProgram )
        self.bigWigToWigPathLineEdit.editingFinished.connect( self.checkBigWigToWigProgram )

        self.whatIsThisButton.clicked.connect( QWhatsThis.enterWhatsThisMode )
        self.convertButton.clicked.connect( self.startProcess )
        self.stopButton.clicked.connect( self.terminateProcessing )
        self.closeButton.clicked.connect( self.close )

    def closeEvent(self, event):

        # if process is running
        if self.process is not None:
            msg = 'Conversion is still running...\n First stop the process.'
            showWarningMessageBox(msg, self)
            event.ignore()
            return

        # In case when hide mode is on, do not close, but hide
        if self.hideMode:
            event.ignore()
            self.hide()
            return

        # Remove temporary h5 file
        if hasattr(self, 'h5RemoveRadioButton'):
            if self.h5RemoveRadioButton.isChecked() and \
                                            self.h5FileLineEdit.text():
                try:
                    os.remove( self.h5FileLineEdit.text() )
                except:
                    pass

        # Remove index file
        if self.indexFile is not None:
            try:
                os.remove(self.indexFile)
            except:
                pass

    def connect2Hide(self):
        """ When opened from browser, hide dialog on close
        """
        try:
            self.closeButton.clicked.disconnect()
        except:
            pass
        self.closeButton.clicked.connect( self.hideDialog )
        self.hideMode = True

    def hideDialog(self):
        if self.process is not None:
            msg = 'Conversion is still running...\n First stop the process.'
            showWarningMessageBox(msg, self)
        else:
            self.hide()

    def removeH5RemoveRadioButton(self):
        """ Remove h5RemoveRadioButton from form.
        It can be useful if this GUI is executed standalone instead from the
        browser
        """
        self.h5RemoveRadioButton.setParent(None)
        del self.h5RemoveRadioButton
        self.h5Line.setParent( None )
        del self.h5Line

    def setH5Name(self, name=None):
        """Set h5 file name.
        In case if h5 name is not provided a random name is generated.
        Also, return the name.
        """
        if name is None:
            self.genRandomH5Name()
        else:
            self.h5FileLineEdit.setText( name )

        return self.h5FileLineEdit.text()

    def setChromName(self, chromName=None):
        """ Set the input chromosome name.
        """
        if chromName is not None:
            self.chromNameLineEdit.setText(chromName)
            self.hasConvertedSuccessfully = False
            self.command = None

    def setInputFile(self, inputFile):
        """Set the input file through code
        """
        self.inFileLineEdit.setText(inputFile)
        self.checkInputFileFormat()

    def enableFileIndexing(self):
        """Generate index file
        """
        chars = string.ascii_letters + string.digits
        name = 	''.join(random.choice(chars) for _ in range(8))
        name = name + '.json'
        self.indexFile = os.path.join( self.scratchDirLineEdit.text(), name )

    def setDefaultScratchDirs(self):
        """ Set default scratch directory
        """
        defaultDir = config['Dirs']['WorkingDirectory']
        self.scratchDirLineEdit.setText(defaultDir)

    def setBigWigInfoPathConfig(self):
        """ Set path of bigWigInfo program from configuration
        """
        path = config['Programs']['bigWigInfo']
        if path != 'None':
            self.bigWigInfoLineEdit.setText( path )

    def setBigWigToWigPathConfig(self):
        """ Set path of bigWigToWig program from configuration
        """
        path = config['Programs']['bigWigToWig']
        if path != 'None':
            self.bigWigToWigPathLineEdit.setText( path )

    def browseScratchDir(self):
        """Browse and choose scratch directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Scratch Directory')
        if path:
            self.scratchDirLineEdit.setText(path)
        self.checkScratchDirectory(fromButton=True)

    def checkScratchDirectory(self, fromButton=False):
        """ Check whether scratch directory is provided correctly in the
        line edit.
        Return: None or scratch directory name
        """
        if (self.scrollArea.focusWidget() is  not self.scratchDirButton) or (fromButton):
            checkDirExist(self.scratchDirLineEdit, self)
        scratchDir = self.scratchDirLineEdit.text()
        if scratchDir:
            return scratchDir
        else:
            if fromButton:
                msg = 'No scratch directory'
                showWarningMessageBox(msg, self)
                self.scratchDirLineEdit.selectAll()
                self.scratchDirLineEdit.setFocus()

    def genRandomH5Name(self):
        """ Generate random name of h5 file in scratch directory
        """
        chars = string.ascii_letters + string.digits
        name = 	''.join(random.choice(chars) for _ in range(8))
        name = name + '.h5'
        h5 = os.path.join( self.scratchDirLineEdit.text(), name )
        self.h5FileLineEdit.setText( h5 )

    def browseInputFile(self):
        """To get input file with full path
        """
        # A dialog box will be displayed to select a file and path will be stored in the cell
        file_choices = " dataset files (*.bigWig *.wig *.bed);;All Files(*.*)"
        path = QFileDialog.getOpenFileName(self, 'Open File', '', file_choices)
        if path[0]:
            self.inFileLineEdit.setText(path[0])
        self.checkInputFileFormat(fromButton=True)

    def checkInputFileFormat(self, fromButton=False):
        """ Check file extension and accessibility of input file.
        """
        if (self.scrollArea.focusWidget() is  not self.inFileButton) or (fromButton):
            fileExist = checkFileExist(self.inFileLineEdit, self)
            if not fileExist:
                return

        fileName = self.inFileLineEdit.text()
        name = os.path.basename( fileName )
        if name and ((self.scrollArea.focusWidget() is  not self.inFileButton) \
                        or (fromButton)):
            extension = os.path.splitext( name )[1]
            fileOK = True
            if extension.lower() == '.bigwig':
                self.fileFormatCBox.setCurrentIndex( 0 )
            elif extension.lower() == '.wig':
                self.fileFormatCBox.setCurrentIndex( 1 )
            elif extension.lower() == '.bed':
                self.fileFormatCBox.setCurrentIndex( 2 )
            else:
                msg = 'File extension [{0}] is not acceptable.'.format(extension)
                showWarningMessageBox(msg, self)
                self.inFileLineEdit.selectAll()
                self.inFileLineEdit.setFocus()
                fileOK = False

            if fileOK:
                return os.path.normpath(fileName)
            else:
                return

        else:
            return

    def browseH5File(self):
        """To get output file with full path
        """
        # A dialog box will be displayed to select a file and path will be stored in the cell
        file_choices = " h5/hdf5 files (*.h5 *.hdf5);;All Files(*.*)"
        path = QFileDialog.getSaveFileName(self, 'Select or Create File', '', file_choices, options=QFileDialog.DontConfirmOverwrite)
        if path[0]:
            self.h5FileLineEdit.setText(os.path.normpath(path[0]))
        self.validateH5OutFile()

    def validateH5OutFile(self):
        """ This can be used to validate internally.
        Return True or None
        """
        h5Out = self.h5FileLineEdit.text()
        if not h5Out:
            return

        h5OutDir = os.path.dirname( self.h5FileLineEdit.text() )
        h5BaseName = os.path.basename( self.h5FileLineEdit.text() )
        h5Extension = os.path.splitext( h5BaseName )[1]
        scratchDir = self.scratchDirLineEdit.text()

        if hasattr(self, 'h5RemoveRadioButton'):
            if (scratchDir == h5OutDir) and not self.h5RemoveRadioButton.isChecked():
                msg = 'You are trying to save h5 file in scratch directory.\n\
                       Please, use another directory.'
                showWarningMessageBox(msg, self)
                self.h5FileLineEdit.selectAll()
                self.h5FileLineEdit.setFocus()
                return
        else:
            if (scratchDir == h5OutDir):
                msg = 'You are trying to save h5 file in scratch directory.\n\
                       Please, use another directory.'
                showWarningMessageBox(msg, self)
                self.h5FileLineEdit.selectAll()
                self.h5FileLineEdit.setFocus()
                return

        # In extension is different, change it to correct one
        if not (h5Extension == '.h5' or h5Extension == '.hdf5'):
            h5Out = os.path.join(h5OutDir, h5BaseName+'.h5')
            self.h5FileLineEdit.setText(h5Out)

        return os.path.normpath(self.h5FileLineEdit.text())

    def isValidH5OutFileFromWidget(self, checked=None):
        """ Check whether scratch directory and h5 directory is same.
        Also, check that path given in h5 box is valid.

        Return:
            If everything is fine: True otherwise: None
        """

        h5Out = self.h5FileLineEdit.text()
        if not h5Out:
            return

        h5OutDir = os.path.dirname( self.h5FileLineEdit.text() )
        h5BaseName = os.path.basename( self.h5FileLineEdit.text() )
        h5Extension = os.path.splitext( h5BaseName )[1]
        scratchDir = self.scratchDirLineEdit.text()

        # Check if h5outdir is available
        if not os.path.isdir( h5OutDir ):
            msg = 'Directory {0} not found.'.format(h5OutDir)
            showWarningMessageBox(msg, self)
            self.h5FileLineEdit.selectAll()
            self.h5FileLineEdit.setFocus()
            return


        # In case if there is no remove button make checked False.
        # Do not save file in scratch directory
        if not hasattr(self, 'h5RemoveRadioButton'):
            checked = False

        if (scratchDir == h5OutDir) and not checked and checked is not None\
            and (self.scrollArea.focusWidget() is  not self.h5FileButton):
            msg = 'You are trying to save h5 file in scratch directory.\n\
                   Please, use another directory.'
            showWarningMessageBox(msg, self)
            self.h5FileLineEdit.selectAll()
            self.h5FileLineEdit.setFocus()
            return

        # If Remove button is checked, then check if scratch directory
        # is not same to output directory
        # This is done only in case remove button is present
        if hasattr(self, 'h5RemoveRadioButton'):
            if self.h5RemoveRadioButton.isChecked():
                checked = True
            else:
                checked = False

            ## show only when user changes focus from line-edit to any
            ## other widget except remove button.
            ## User should be able to click on remove button
            if (scratchDir == h5OutDir) and (not checked) and \
            (self.scrollArea.focusWidget() is  not self.h5RemoveRadioButton) \
            and (self.scrollArea.focusWidget() is  not self.h5FileButton):
                msg = 'You are trying to save h5 file in scratch directory.\n\
                       Please, use another directory.'
                showWarningMessageBox(msg, self)
                self.h5FileLineEdit.selectAll()
                self.h5FileLineEdit.setFocus()
                return

        # In extension is different, change it to correct one
        if not (h5Extension == '.h5' or h5Extension == '.hdf5'):
            h5Out = os.path.join(h5OutDir, h5BaseName+'.h5')
            self.h5FileLineEdit.setText(h5Out)

        return True

    def checkResolutions(self):
        """ Check resolutions inputs.
        Return resolutions list or None
        """
        rInList = self.resolutionLineEdit.text()
        rlist = []
        if rInList:
            temp = re.split(',', rInList)
            for r in temp:
                r = r.rstrip().lstrip()
                if not r.strip():
                    continue
                try:
                    b = gmlib.util.resolutionToBinsize(r)
                    rlist.append( gmlib.util.binsizeToResolution(b) )
                except ValueError:
                    msg = '"{0}" contains resolution "{1}"...\n which is not an acceptable resolution.'.format(rInList, r)
                    showWarningMessageBox(msg, self )

            return rlist

    def checkDownsamplingMethods(self, checked=None):
        """Check downsampling method.
        Return method list or None
        """
        methods = []
        if self.ameanCheckBox.isChecked():
            methods.append('amean')
        if self.hmeanCheckBox.isChecked():
            methods.append('hmean')
        if self.gmeanCheckBox.isChecked():
            methods.append('gmean')
        if self.maxCheckBox.isChecked():
            methods.append('max')
        if self.minCheckBox.isChecked():
            methods.append('min')
        if self.medianCheckBox.isChecked():
            methods.append('median')

        if not methods:
            msg = 'Please select at least one downsampling method.'
            showWarningMessageBox(msg, self )
            self.coarseGroupBox.setFocus()
            return None
        else:
            return methods

    def browseBigWigInfoProgram(self):
        """To get bigWigInfo path
        """
        # A dialog box will be displayed to select a file and path will be stored in the cell
        file_choices = " Binary Executable (* *.exe)"
        path = QFileDialog.getOpenFileName(self, 'Open File', '', file_choices)
        if path[0]:
            self.bigWigInfoLineEdit.setText(path[0])
        self.checkBigWigInfoProgram(fromButton=True)

    def checkBigWigInfoProgram(self, fromButton=False):
        """Check if path to bigWigInfo is valid
        """
        if (self.scrollArea.focusWidget() is  not self.bigWigInfoButton) or (fromButton):
            checkFileExist(self.bigWigInfoLineEdit, self)
        bigWigInfoProgram = self.bigWigInfoLineEdit.text()
        if bigWigInfoProgram:
            return bigWigInfoProgram
        else:
            if fromButton:
                msg = 'No bigWigInfo !!!'
                showWarningMessageBox(msg, self)
                self.bigWigInfoLineEdit.selectAll()
                self.bigWigInfoLineEdit.setFocus()

    def browseBigWigToWigProgram(self):
        """To get bigWigToWig path
        """
        # A dialog box will be displayed to select a file and path will be stored in the cell
        file_choices = " Binary Executable (*)"
        path = QFileDialog.getOpenFileName(self, 'Open File', '', file_choices)
        if path[0]:
            self.bigWigToWigPathLineEdit.setText(path[0])
        self.checkBigWigToWigProgram(fromButton=True)

    def checkBigWigToWigProgram(self, fromButton=False):
        """ Check path to bigWigToWig is valid
        """
        if (self.scrollArea.focusWidget() is  not self.bigWigToWigPathButton) or (fromButton):
            checkFileExist(self.bigWigToWigPathLineEdit, self)
        bigWigToWigProgram = self.bigWigToWigPathLineEdit.text()
        if bigWigToWigProgram:
            return bigWigToWigProgram
        else:
            if fromButton:
                msg = 'No bigWigToWig !!!'
                showWarningMessageBox(msg, self)
                self.bigWigToWigPathLineEdit.selectAll()
                self.bigWigToWigPathLineEdit.setFocus()

    def buildCommand(self):
        """ Build the command
        """

        self.command = None
        command = dict()

        # Input file
        command['-i'] = self.checkInputFileFormat(fromButton=True)
        if command['-i'] is None:
            showWarningMessageBox('No input file.', self)
            self.inFileLineEdit.selectAll()
            self.inFileLineEdit.setFocus()
            return False
        command['-i'] = '"{0}"'.format(command['-i'])

        # Input file format
        inFormat = self.fileFormatCBox.currentText()

        # Chromosome name
        if self.chromNameLineEdit.text():
            command['-icn'] = self.chromNameLineEdit.text()

        # Dataset Title
        command['-t'] = '"' + self.titleLineEdit.text() + '"'

        # Resolutions list
        resolutions = self.checkResolutions()
        if resolutions is not None:
            command['-r'] = '"' + ', '.join(resolutions) + '"'

        # List of downsampling methods
        downsampling_methods= self.checkDownsamplingMethods()
        if downsampling_methods is None:
            return False
        else:
            command['-dm'] = '"' + ', '.join(downsampling_methods) + '"'

        # Scratch directory
        workDir = self.checkScratchDirectory(fromButton=True)
        if workDir is None:
            return False
        else:
            if workDir != config['Dirs']['WorkingDirectory']:
                command['-wd'] = '"{0}"'.format(workDir)

        # Output file
        command['-o'] = self.validateH5OutFile()
        if command['-o'] is None:
            return False
        command['-o'] = '"{0}"'.format(command['-o'])

        # Compression method in output file
        command['-cmeth'] = self.compressionCBox.currentText().lower()

        # Path to bigwig programs
        if inFormat == 'bigWig':
            bigWigInfo = self.checkBigWigInfoProgram(fromButton=True)
            if bigWigInfo is None:
                return False
            else:
                if bigWigInfo != config['Programs']['bigWigInfo']:
                    command['-binfo'] = '"{0}"'.format(bigWigInfo)

            bigWigToWig = self.checkBigWigToWigProgram(fromButton=True)
            if bigWigToWig is None:
                return False
            else:
                if bigWigToWig != config['Programs']['bigWigToWig']:
                    command['-b2w'] = '"{0}"'.format(bigWigToWig)

        # Data column in bed file
        if inFormat == 'bed':
            command['-dtc'] = str( self.bedDataColumnSpinBox.value() )

        # Keep original data
        if self.keepOriginalButton.isChecked():
            command['-ko'] = ' '

        # Overwrite h5 dataset
        command['-ow'] = ' '

        # json Index file
        if self.indexFile is not None:
            command['-idf'] = '"{0}"'.format(self.indexFile)

        # Construct command
        commandText = ' '
        if inFormat == 'bigWig':
            commandText += 'bigwig2h5 '

        if inFormat == 'wig':
            commandText += 'wig2h5 '

        if inFormat == 'bed':
            commandText += 'bed2h5 '

        for key in command:
            commandText += key + ' ' + command[key] + ' '

        self.command = commandText

        return True

    def startProcess(self):
        if not self.buildCommand():
            return

        if self.command is None:
            return

        self.process = QProcess(self)
        self.hasConvertedSuccessfully = False

        self.process.start('gcMapExplorer', shlex.split(self.command))
        self.process.setProcessChannelMode( QProcess.MergedChannels )
        #self.process.setReadChannel( QProcess.StandardOutput )
        self.process.waitForStarted()

        self.process.readyReadStandardOutput.connect( self.writeLogOutputFromProcess )
        self.process.readyReadStandardError.connect( self.writeLogOutputFromProcess )
        self.process.finished.connect( self.finishedProcessing )

        self.convertButton.setEnabled(False)
        self.stopButton.setEnabled(True)
        self.logOutputPlainTextEdit.clear()

    def writeLogOutputFromProcess(self):
        """Converting process output to log box
        """
        if self.process is None:
            return

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

    def finishedProcessing(self):
        self.convertButton.setEnabled(True)
        self.stopButton.setEnabled(False)
        exitStatus = self.process.exitStatus()
        self.process = None
        # Store result if process was normal exited
        if exitStatus == QProcess.NormalExit:
            self.hasConvertedSuccessfully = True
            # Hide automatically when process exited normally
            if self.hideMode:
                self.accept()
                self.disconnectFromAccepted()

    def connectToBrowserDataLoader(self, browser, gpa):
        self.accepted.connect( lambda : self.loadDataToBrowser(browser, gpa) )

    def disconnectFromAccepted(self):
        """ disconnect accepted
        """
        try:
            self.accepted.disconnect()
        except:
            pass

    def loadDataToBrowser(self, browser, gpa):
        """ This is used in the browser.
        Do not use in standalone mode.
        """

        if self.result() == QDialog.Accepted:
            # Open dialogbox so user can select a dataset
            gpa.selectGenomicDataHdf5ByDialogBox(self.h5FileLineEdit.text(), browser.filesOpened)
            self.centralOperationWidget.setEnabled( False )

        browser.loadDataToPlot(gpa)

    def connectToBrowseMapChanger(self, browser, newMapName):
        """ Connect accepted to browser map changer
        """
        self.accepted.connect( lambda : browser.changeMapNames(newMapName) )

def showWarningMessageBox(msg, qwidget):
    msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg, QMessageBox.Ok, qwidget)
    msgBox.exec_()
    msgBox.close()

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

def checkFileExist(lineEdit, qwidget):
    if not lineEdit.text(): return False
    filename = str( lineEdit.text() )
    if not os.path.isfile( filename ):
        msg = "[ {0} ] \n Not found !!!".format(filename)
        msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg,
                                                    QMessageBox.Ok, qwidget)
        msgBox.exec_()
        msgBox.close()

        lineEdit.selectAll()
        lineEdit.setFocus()

        return False
    else:
        return True

if __name__ == '__main__':
    main()
