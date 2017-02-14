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

import logging
import os

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.uic import loadUiType


# Determine absolute path to UIs directory. Relative path from this directory does not work.
DirToThisScript = os.path.dirname(os.path.abspath(__file__))
PathToUIs = os.path.join(DirToThisScript, 'UIs')

## To divert log output into QPlainTextWidget
class pyQtLogHandler(QObject, logging.Handler):
    messageEmitted = pyqtSignal( str )
    def __init__(self):
        super().__init__()
        self.setFormatter(logging.Formatter('%(levelname)s:%(name)s:%(message)s'))

    def emit(self, record):
        msg = self.format(record)
        self.messageEmitted.emit(msg)

class qtThread(QThread):
    resultReady = pyqtSignal()

    def __init__(self, parent=None, target=None, args=(), kwargs=None):
        super(qtThread, self).__init__(parent)
        self.results = None
        self._target = target
        self._args = args
        self._kwargs = kwargs


    def callRun(self):
        try:
            if self._target:
                if self._kwargs is not None:
                    self.results = self._target(*self._args, **self._kwargs)
                else:
                    self.results = self._target(*self._args)
        finally:
            # Avoid a refcycle if the thread is running a function with
            # an argument that has a member that points to the thread.
            del self._target, self._args, self._kwargs

        self.resultReady.emit()

    def run(self):
        self.callRun()


def qThreadTerminated():
    currentThread = qtThread.currentThread()
    if currentThread is not None:
        if currentThread.isInterruptionRequested():
            return True
        else:
            return False
    else:
        return False


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

        lineEdit.undo()
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

        lineEdit.undo()
        lineEdit.setFocus()

def showWarningMessageBox(msg, qwidget):
    msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg, QMessageBox.Ok, qwidget)
    msgBox.exec_()
    msgBox.close()

def constrainValueInLineEdit(qwidget, lineEdit, minvalue, maxvalue,):
    """ Whether value in a lineEdit is within limit.

    It checks whether value in a lineEdit is between the limit, otherwise
    throw a warning to user.

    """
    value = float(lineEdit.text())

    if not value:   return

    if value > maxvalue or value < minvalue:
        msg = "Value should be between {0} and {1} !!!".format(minvalue, maxvalue)
        msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg,
                                                    QMessageBox.Ok, qwidget)
        msgBox.exec_()
        msgBox.close()

        lineEdit.selectAll()
        lineEdit.setFocus()

def getFileExtension(name):
    return os.path.splitext(name)[1]
