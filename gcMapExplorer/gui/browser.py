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

import sys, os, random
from enum import Enum
import pdb
import numpy as np
import math
from numpy import ma
import h5py
np.seterr(all='ignore')

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtPrintSupport import *
from PyQt5.uic import loadUiType

import matplotlib as mpl
from matplotlib import figure as mplFigure
from matplotlib import widgets as malWidgets
mpl.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import gcMapExplorer.lib as gmlib

from . import browserHelpers
from . import guiHelpers
from . import h5Converter

# Determine absolute path to UIs directory. Relative path from this directory does not work.
DirToThisScript = os.path.dirname(os.path.abspath(__file__))
PathToUIs = os.path.join(DirToThisScript, 'UIs')

#import mainwindow
pathToThisUI = os.path.join(PathToUIs, 'mainWindow.ui')
Ui_MainWindow, QMainWindow = loadUiType(pathToThisUI)

class ePageType(Enum):  A5=1; A4=2; A3=3; A2=4; A1=5; A0=6; Letter=7; Custom=8
PageSize = {  ePageType.A4 : (8.27, 11.7),
              ePageType.Letter : (8.5, 11),
              ePageType.A0 : (33.1, 46.8),
              ePageType.A1 : (23.4, 33.1),
              ePageType.A2 : (16.5, 23.4),
              ePageType.A3 : (11.7, 16.5),
              ePageType.A5 : (5.83, 8.27),
              ePageType.Custom : (8.27, 11.7),
            }


def main():
    app = QApplication(sys.argv)
    FrontGUI = Main()
    FrontGUI.show()
    app.exec_()
    app.exit()

class GenomicDataSetSubPlotHelper:
    ''' This is a helper class for genomic dataset subplots.
        Do not use it separately.
    '''
    def connectGenomicDatasetOptions(self):
        # Connect y-scalar for genomic plots
        self.genomicDataYScalingSlider.valueChanged.connect(self.changeYScaleGenomicSubplotSlider)
        self.genomicDataYScalingSpinBox.valueChanged.connect(self.changeYScaleGenomicSubplotSpinbox)
        self.genomicDataYScalingMaxLineEdit.setValidator(QDoubleValidator())
        self.genomicDataYScalingMaxLineEdit.returnPressed.connect( self.updateGenomicDataYScaleLimits )
        self.genomicDataYScalingMinLineEdit.setValidator(QDoubleValidator())
        self.genomicDataYScalingMinLineEdit.returnPressed.connect( self.updateGenomicDataYScaleLimits )
        self.genomicDataYScalingGoButton.clicked.connect( self.updateGenomicDataYScaleLimits )
        self.genomicDataYScalingResetButton.clicked.connect( self.resetGenomicDataYScaleLimits )

        self.genomicDataColorButton.clicked.connect( self.changeGenomicDatasetPlotColor )
        self.genomicDataLineWidthSpinBox.valueChanged.connect( self.changeGenomicDatasetPlotLineWidth )
        self.genomicDataPercentileSpinBox.valueChanged.connect(self.displayGenomicDataPercentile)


    def openGenomicDatasetFile(self, hiCmapAxis):
        """Open genomic dataset file, load it and plot it.
        It is a front level method, which uses several lower level methods.

            * Make a new axis, and plot the data
            * Make a new tree widget and select it
            * enables all connected GUI options.
        """
        idx = hiCmapAxis.index

        # Open file selector dialog box
        file_choices = " Compatible files (*.h5 *.hdf *.bigWig *.wig *.bed)"
        file_choices += ";;HDF5 file (*.h5 *.hdf)"
        file_choices += ";;IGB file  (*.bigWig *.wig *.bed)"
        file_choices += ";;Text file (*.txt *.dat)"
        file_choices += ";;All files (*.*)"

        path = QFileDialog.getOpenFileName(self, 'Load File', '/home', file_choices)

        # If a file is selected by user
        if path[0]:
            file_extension = os.path.splitext(path[0])[1]

            if self.hiCmapAxes[idx].genomicPlotAxes is None:
                self.hiCmapAxes[idx].genomicPlotAxes = []

            # Instantiate GenomicDataPlotAxis
            gpa = GenomicDataPlotAxis( len(self.hiCmapAxes[idx].genomicPlotAxes),  self.hiCmapAxes[idx])

            if file_extension in ['.h5', '.hdf5']:
                # Open dialogbox so user can select a dataset
                gpa.selectGenomicDataHdf5ByDialogBox(path[0], self.filesOpened)
            if file_extension in ['.bigWig', '.wig', '.bed']:
                # Open dialogbox for conversion
                gpa.selectGenomicDataByH5Conversion(self, path[0], self.hiCmapAxes[idx].ccmap.xlabel, self.filesOpened)
            elif file_extension == '.txt' or file_extension == '.dat':
                gpa.readDataFromTextFile(path[0])

            # plot genomic dataset on the browser
            self.loadDataToPlot(gpa)

    def loadDataToPlot(self, gpa):
        """ plot genomic dataset on the browser
        """
        idx = gpa.hiCmapAxis.index

        # if a Dataset is selected by user, then process further
        if gpa.dataArray is not None:
            # Generate treeWidget item and add it as a child to ccmap axis treeWidget
            gpa.set_tree_widget_item(self.axisTreeWidget)


            # Append GenomicDataPlotAxis to list
            self.hiCmapAxes[idx].genomicPlotAxes.append(gpa)

            # Getting upper-most and lower-most plot
            if gpa.plotLocation == 'top':
                self.hiCmapAxes[idx].upperMostGenomicPlotAxes = gpa.index
            if gpa.plotLocation == 'bottom':
                self.hiCmapAxes[idx].lowerMostGenomicPlotAxes = gpa.index

            # Tight layout at start and then turn off it
            self.figure.set_tight_layout(True)

            # Plot the dataset
            self.hiCmapAxes[idx].updateGenomicPlotAxes()

            # set the treewidget to current widget
            self.axisTreeWidget.setCurrentItem(gpa.treeWidgetItem)

            self.canvas.draw()
            self.figure.set_tight_layout(False)

            # make active the GUI options
            self.makeGenomicSupPlotOptionsActive(idx, self.hiCmapAxes[idx].genomicPlotAxes[-1].index)

            # Get horizontal and vertical space between subplots
            self.get_horizontal_vertical_space_from_figure()

            # Make list of common data name
            self.hiCmapAxes[idx].tryEnableInterchangeCMapName()

    def makeGenomicSupPlotOptionsActive(self, hidx, gidx):
        """Set and Make genomic subplot y-scalar active
        To enable all GUI options related to genomic subplot
        """
        # At first entire group box is enabled for interaction with user
        self.genomicDataYScalingGroupBox.setEnabled(True)
        self.genomicDataSubPlotOptionsGBox.setEnabled(True)

        # If genomic dataset is plotted for first time, initialize spinbox and slider value
        if self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yscaleSpinbox is None:
            self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yscaleSpinbox = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].ylimit[1]
        if self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yscaleSlider is None:
            self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yscaleSlider = 99

        # Because, all GUI widget are linked to plotting, internal change in their states trigger all the connected methods.
        # All plottings are already done, so turn off the plotting.
        self.hiCmapAxes[hidx].doNotPlot = True

        minvalue = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[0]
        maxvalue = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[-1]
        dstep = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[1] - self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[0]

        # Change maximum and minimum value
        self.genomicDataYScalingMinLineEdit.setText(str(minvalue))
        self.genomicDataYScalingMaxLineEdit.setText(str(maxvalue))

        # VERY IMPORTANT for SpinBox: always call setRange, then setSingleStep and then setValue. NEVER EVER call setValue before.
        self.genomicDataYScalingSpinBox.setRange(minvalue, maxvalue)
        self.genomicDataYScalingSpinBox.setSingleStep(dstep)

        # change spinbox and slider value
        self.genomicDataYScalingSlider.setValue(self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yscaleSlider)
        self.genomicDataYScalingSpinBox.setValue(self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yscaleSpinbox)

        self.hiCmapAxes[hidx].doNotPlot = False

        # Change color and linewidth
        color = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].plotColor
        style = 'background-color: rgb({0}, {1}, {2});'.format(int(color[0]*255), int(color[1]*255), int(color[2]*255))
        self.genomicDataColorButton.setStyleSheet(style)

    def makeGenomicSupPlotOptionsInactive(self):
        # At first entire group box is enabled for interaction with user
        self.genomicDataYScalingGroupBox.setEnabled(False)
        self.genomicDataSubPlotOptionsGBox.setEnabled(False)

    def getCurrentHicmapAndGenomicSubplotIndex(self):
        """Method to get indexes of currently active ccmap axes and respective genomic subplot
        """
        hidx = None
        gidx = None

        if self.hiCmapAxes is None: return hidx, gidx

        hidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[hidx].image is None: return hidx, gidx

        if self.hiCmapAxes[hidx].genomicPlotAxes is not None:
            for gax in self.hiCmapAxes[hidx].genomicPlotAxes:
                if gax.treeWidgetItem is self.axisTreeWidget.currentItem():
                    gidx = gax.index
                    break

        return hidx, gidx

    def changeYScaleGenomicSubplotSlider(self, newValue=None):
        """When either user changes slider value or value are changed internally, this function is called

        This is the only method, which actively updates the plot.
        Other GUI widget connected methods either directly call this method or invoke it by changing slider value internally.

        """
        # Get index of axes
        hidx, gidx = self.getCurrentHicmapAndGenomicSubplotIndex()

        if hidx is None and gidx is None:   return
        if hidx is None:    return
        if gidx is None:    return

        if self.hiCmapAxes[hidx].doNotPlot: return

        # Get value from slider
        sidx = self.genomicDataYScalingSlider.value()
        ylimUpper = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[sidx]
        ylimLower = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[0]

        # To prevent plotting when user slides slider at the lowest value
        if ylimLower == ylimUpper:  return

        # If function is not called via spinbox value changed
        if ylimUpper != self.genomicDataYScalingSpinBox.value():
            self.genomicDataYScalingSpinBox.setValue(ylimUpper)

        # set upper and lower ylimit for plot, here plot is also updated, see: ylimit.setter
        self.hiCmapAxes[hidx].genomicPlotAxes[gidx].ylimit = (ylimLower, ylimUpper)

        # Save current values from slider and spinbox
        self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yscaleSlider = sidx
        self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yscaleSpinbox = ylimUpper

        # Redraw everything
        self.canvas.draw()

    def changeYScaleGenomicSubplotSpinbox(self, newValue=None):
        """When user changes spinbox, this method is called

        This function is called when user changes maximum and minimum limit of y-scale.
        This function update both spinbox and slider values, which eventually update and redraw plot by method connected to slider.

        A the end, it either directly calls slider connected function or indirectly calls it by setting slider value.

        """
        # Get index of axes
        hidx, gidx = self.getCurrentHicmapAndGenomicSubplotIndex()

        if hidx is None and gidx is None:   return
        if hidx is None:    return
        if gidx is None:    return

        if self.hiCmapAxes[hidx].doNotPlot: return

        # Determine value for slider for current spinbox value
        idx = ( np.abs(self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps - self.genomicDataYScalingSpinBox.value() )).argmin()

        # If determined value is not equal to present value, set slider value.
        # By setting slider value, plot is automatically updated and redrawn.
        if self.genomicDataYScalingSlider.value() != idx:
            self.genomicDataYScalingSlider.setValue(idx)
        else:
            # However, when maximum and minimum limits are changed by user, sometimes, determined slider value,
            # and present slider value remains same. Therefore, setting the slider value does not invoke its connected function.
            # Here, this connected function is called to update and redraw the plot
            self.changeYScaleGenomicSubplotSlider()

    def updateGenomicDataYScaleLimits(self):
        """When user changes y scale maximum and minimum values on GUI.

        At the end, it either directly calls spinbox connected function or indirectly by setting spinbox value.
        """
        # Get index of axes

        hidx, gidx = self.getCurrentHicmapAndGenomicSubplotIndex()

        if hidx is None and gidx is None:   return
        if hidx is None:    return
        if gidx is None:    return

        if self.hiCmapAxes[hidx].doNotPlot: return

        # to get present lower and upper limit entered by user
        lowerLimit = float( self.genomicDataYScalingMinLineEdit.text() )
        upperLimit = float( self.genomicDataYScalingMaxLineEdit.text() )

        # To get previous lower and upper limit already used in plot
        prev_lowerLimit = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[0]
        prev_upperLimit = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[-1]

        # When user make a mistake, so resetting  the upper and lower limits
        if lowerLimit >= upperLimit:
            # Display message in message box
            msg = ''' Minimum value ≥ Maximum value !!!
                            {0}  ≥   {1}
                  ''' .format(lowerLimit, upperLimit)

            msgBox = QMessageBox()
            msgBox.setWindowTitle('Warning')
            msgBox.setText(msg)
            msgBox.setInformativeText('Click Cancel to revert the values.')
            msgBox.setStandardButtons(QMessageBox.Cancel)
            msgBox.setDefaultButton(QMessageBox.Cancel)
            msgBox.exec_()

            # Revert the values in both boxes
            self.hiCmapAxes[hidx].genomicPlotAxes[gidx].doNotPlot = True
            self.genomicDataYScalingMinLineEdit.setText(str(prev_lowerLimit))
            self.genomicDataYScalingMaxLineEdit.setText(str(prev_upperLimit))
            self.hiCmapAxes[hidx].genomicPlotAxes[gidx].doNotPlot = False
            return

        # Update y-scale steps for plots
        self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps = (lowerLimit, upperLimit)

        # New delta-step
        dstep = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[1] - self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[0]

        # New range for spin-box and delta-step
        self.genomicDataYScalingSpinBox.setRange(lowerLimit, upperLimit)
        self.genomicDataYScalingSpinBox.setSingleStep(dstep)
        self.genomicDataYScalingSpinBox.setDecimals(self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yticksDecimals+1)

        # If user changes only lower limit, spinbox connected method is called because spinbox holds only upper limit
        if prev_lowerLimit != lowerLimit:
            self.changeYScaleGenomicSubplotSpinbox()
        else:
            # If user changes upper limit, spinbox value is changed, which calls its connected method
            self.genomicDataYScalingSpinBox.setValue(upperLimit)

    def resetGenomicDataYScaleLimits(self, hidx=None, gidx=None):
        """When user clicks reset button
        """

        if gidx is None:
            hidx, gidx = self.getCurrentHicmapAndGenomicSubplotIndex()

        if hidx is None and gidx is None:   return
        if hidx is None:    return
        if gidx is None:    return

        # Remember, by setting y-scale steps to None actually changes steps to original. See @yScaleSteps.setter
        self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps = None

        # Determine new upper and lower limit, and delta-step value
        lowerLimit = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[0]
        upperLimit = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[-1]
        dstep = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[1] - self.hiCmapAxes[hidx].genomicPlotAxes[gidx].yScaleSteps[0]

        # Change range and delta-step of spinbox
        self.genomicDataYScalingSpinBox.setRange(lowerLimit, upperLimit)
        self.genomicDataYScalingSpinBox.setSingleStep(dstep)

        # Change value in maximum and minimum box
        self.genomicDataYScalingMinLineEdit.setText(str(lowerLimit))
        self.genomicDataYScalingMaxLineEdit.setText(str(upperLimit))

        # Change spinbox value, which changes slider value, which updates and redraw plots
        # In extremely rare case, when new upper limit equals to present spinbox value, call spinbox connected function
        if self.genomicDataYScalingSpinBox.value() == upperLimit:
            self.changeYScaleGenomicSubplotSpinbox()
        else:
            self.genomicDataYScalingSpinBox.setValue(upperLimit)

    def displayGenomicDataPercentile(self, value):
        """Display percentile value for a given percentile

        It only gives information about the percentile value
        """
        # Get index of axes

        hidx, gidx = self.getCurrentHicmapAndGenomicSubplotIndex()

        if hidx is None and gidx is None:   return
        if hidx is None:    return
        if gidx is None:    return

        if self.hiCmapAxes[hidx].doNotPlot: return

        if self.hiCmapAxes[hidx].genomicPlotAxes[gidx].dataArray is None:  return

        dataArray = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].dataArray

        if value not in self.hiCmapAxes[hidx].genomicPlotAxes[gidx].percentileValues:
            percentile = np.percentile(dataArray[dataArray[:] != 0 ], value)
            self.hiCmapAxes[hidx].genomicPlotAxes[gidx].percentileValues[value] = percentile
        else:
            percentile = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].percentileValues[value]

        self.percentileValueLineEdit.setText(str(percentile))


    def changeGenomicDatasetPlotColor(self):
        """ Set the color of genomic dataset plot
        """

        # Get index of axes
        hidx, gidx = self.getCurrentHicmapAndGenomicSubplotIndex()

        if hidx is None and gidx is None:   return
        if hidx is None:    return
        if gidx is None:    return

        # Convert RGB tuple to QColor, QColor range from 0 to 255 while matplotlib color range from 0 to 1.
        color = self.hiCmapAxes[hidx].genomicPlotAxes[gidx].plotColor
        qcolor = QColor.fromRgb( int(color[0]*255), int(color[1]*255), int(color[2]*255), 255 )

        # QColorDialog open here
        colorDialog = QColorDialog(qcolor, self)
        colorDialog.setWindowTitle("Choose a Color")
        colorDialog.exec_()

        if colorDialog.result() == QDialog.Accepted:
            # Get new color from the user
            pickedColor = colorDialog.selectedColor()
            newColor = (pickedColor.red()/255, pickedColor.green()/255, pickedColor.blue()/255)

            # Update color and plot
            self.hiCmapAxes[hidx].genomicPlotAxes[gidx].plotColor = newColor
            self.hiCmapAxes[hidx].genomicPlotAxes[gidx].updatePlot()
            self.canvas.draw()

            # Set background color of button
            style = 'background-color: rgb({0}, {1}, {2});'.format(pickedColor.red(), pickedColor.green(), pickedColor.blue())
            self.genomicDataColorButton.setStyleSheet(style)

        del colorDialog
        del qcolor

    def changeGenomicDatasetPlotLineWidth(self, width):
        """ Change plot line width
        """
        # Get index of axes
        hidx, gidx = self.getCurrentHicmapAndGenomicSubplotIndex()

        if hidx is None and gidx is None:   return
        if hidx is None:    return
        if gidx is None:    return

        self.hiCmapAxes[hidx].genomicPlotAxes[gidx].plotLineWidth = width
        self.hiCmapAxes[hidx].genomicPlotAxes[gidx].updatePlot()
        self.canvas.draw()

    def get_subplot_obj_axis_under_mouse_mpl(self, event):
        if self.hiCmapAxes is not None:
            for hax in self.hiCmapAxes:
                if hax.genomicPlotAxes is not None:
                    for gax in hax.genomicPlotAxes:
                        contains, attrib = gax.ax.contains(event)
                        if contains:
                            return gax

    def get_subplot_obj_axis_under_mouse_qt(self, qpoint):
        if self.hiCmapAxes is not None:
            cw, ch = self.canvas.get_width_height()
            for hax in self.hiCmapAxes:
                if hax.genomicPlotAxes is not None:
                    for gax in hax.genomicPlotAxes:
                        bbox = gax.ax.get_position()
                        if bbox.contains(qpoint.x()/cw, 1-(qpoint.y()/ch)):
                            return gax


class Main(QMainWindow, Ui_MainWindow, GenomicDataSetSubPlotHelper):

    def __init__(self, ):
        super(Main, self).__init__()
        self.setupUi(self)

        self.setWindowState(self.windowState() or Qt.WindowMaximized)
        """ Window maximized at start """

        self.menuRightClick = None
        """ Right click QMenu object """

        self.DialogAxisProps = None
        """ Axis Properties QTabWidget Object """

        ## Analysis dialogs
        self.corrDialog = None # correlation dialog box

        self.binsDisplayed = 1000
        """ Number of bins shown on plots along X and Y axis. This will control the X/Y range for all plots.
            Default number is the lowest of all maps if it is less than 1000.
        """

        self.interchangeableResolutions = None
        """List of resolution in ascending order which are common between maps
        """

        self.currentResolutionIndex = None
        """Holds the index of current resolution index of self.interchangeableResolutions
        """

        self.interchangeableCMapNames = None
        """ List of interchangeable contact map names
        """

        self.menuActionsToAddGenomicDataset = None
        """This list stores menu QAction of all CCMAPAXIS
        """


        self.filesOpened = dict()
        """ Dictionary of all opened files to its file stream object """

        self.colorMapsDictionary = None
        """ List of color maps """

        self.press_on_plot = None
        self.figure = None
        self.canvas = None
        self.hiCmapAxes = None
        """ List of HiC Map axis """
        self.ActiveHiCmapAxis = None
        """ index to the active HiC Map axis """

        self.marker = None
        self.markerColor = (1, 0, 0)  # Default Red
        self.markerLineWidth = 0.5    # Marker line width
        self.label_mij = None
        self.fontsize = 14

        self.figsize = PageSize[ePageType.A4]
        """Set page size of plot to A4 """

        self.PageType = ePageType.A4
        """Set page type enum """

        self.PageOrientation = 'Portrait'
        """Set page Orientation """

        # Horizontal and vertical space between plots
        self.horizontalSpace = None
        self.verticalSpace = None

        # Add options to toolbar. In GUI designer these options are outside of toolbars
        self.add_options_to_toolbar()

        # Connect menus and actions to methods
        self.add_menu_triggered_action()

        # Connect GUI options to methods
        self.connect_gui_to_options()
        self.connectGenomicDatasetOptions()

        # Embed matplotlib canvas inside scroll area
        self.embed_mpl_canvas_in_mainwindow()

        #TODO: temporary added to load map from command line
        if len(sys.argv) > 2:
            self.add_new_ccmap_axes()
            self.enable_gui_options()
            self.InitMapImage(sys.argv[2], 'ccmap')


    def closeEvent(self, event):
        """Close any additional opened dialog box

            Remove all opened files
        """
        if self.DialogAxisProps is not None:
            self.DialogAxisProps.close()
        if self.corrDialog is not None:
            self.corrDialog.close()

        # Close any opened files
        for key in self.filesOpened:
            self.filesOpened[key].close()

        if self.hiCmapAxes is not None:
            for hax in self.hiCmapAxes:
                del hax.ccmap
                if hax.genomicPlotAxes is not None:
                    for gax in hax.genomicPlotAxes:

                        # Close converter dialog
                        if gax.converterDialog is not None:
                            gax.converterDialog.hideMode = False
                            try:
                                gax.converterDialog.close()
                            except:
                                pass

                        if gax.txtFileHand is not None:
                            del gax.txtFileHand
                        del gax
                    del hax

    def add_menu_triggered_action(self):
        ''' This function connects the menus and actions to the methods.
            It is called only during initialization.

        .. note::

            * Menu action to add genomic subplots for each ccmap plot is added dynamically in addMenuActionsForGenomicDatasetSelector().


        '''
        # Connects File menu
        self.actionLoad_Hi_C_Map_file.triggered.connect(self.open_map_pyobj)
        self.actionAdd_Hi_C_Map.triggered.connect(self.open_map_pyobj)
        self.actionSave_plot.triggered.connect(self.save_plot)
        self.actionQuit.triggered.connect(self.close)
        self.actionPrint.triggered.connect(self.printPlot)

        # Create new action-group for page-size
        self.actionGroupPageSize = QActionGroup(self)
        # Add to action group
        self.actionPlotSizeLetter.setActionGroup(self.actionGroupPageSize)
        self.actionPlotSizeA5.setActionGroup(self.actionGroupPageSize)
        self.actionPlotSizeA4.setActionGroup(self.actionGroupPageSize)
        self.actionPlotSizeA3.setActionGroup(self.actionGroupPageSize)
        self.actionPlotSizeA2.setActionGroup(self.actionGroupPageSize)
        self.actionPlotSizeA1.setActionGroup(self.actionGroupPageSize)
        self.actionPlotSizeA0.setActionGroup(self.actionGroupPageSize)
        self.actionPlotSizeCustom.setActionGroup(self.actionGroupPageSize)

        # Connect action group to method
        self.actionGroupPageSize.triggered.connect(self.set_page_size)

        # Create new action-group for Page Orientation
        self.actionGroupPageOrientation = QActionGroup(self)
        self.actionPlotOrientationLandscape.setActionGroup(self.actionGroupPageOrientation)
        self.actionPlotOrientationPortrait.setActionGroup(self.actionGroupPageOrientation)
        self.actionGroupPageOrientation.triggered.connect(self.set_page_orientation)

        # user defined color dialog
        self.actionAddModifyColormap.triggered.connect( self.showUserDefinedColormapDialog )

        # Analysis menu
        self.actionCorrelationMaps.triggered.connect( self.analysisCorrelateMaps )

        # Help menu
        self.actionWhatIsThis.triggered.connect( QWhatsThis.enterWhatsThisMode )
        self.actionAbout.triggered.connect( self.on_about )

    def connect_gui_to_options(self):
        ''' This function connects all options to respective methods.
            Only used during initialization.
        '''
        # Connect axis tree widget and its right click action
        self.axisTreeWidget.itemClicked.connect(self.make_ccmap_active_on_click)
        self.axisTreeWidget.currentItemChanged.connect(self.make_ccmap_active_on_click)
        self.axisTreeWidget.customContextMenuRequested.connect(self.ShowRightClickMenuTreeWidget)

        # Reset maps
        self.reset_maps_button.clicked.connect(self.ResetMap)
        browserHelpers.add_icon_to_widget(self.reset_maps_button, 'resetArrow.png')

        self.contactMapNameCBox.currentTextChanged.connect( self.changeMapNames )

        # Main setting options
        # Hide tabbars
        tabbars = self.settingsTabWidget.findChildren(QTabBar)
        for tabbar in tabbars:
            tabbar.hide()
        self.settingsCBox.currentIndexChanged.connect( self.settingsTabWidget.setCurrentIndex )

        # Connect marker options
        self.markerCBox.currentIndexChanged.connect(self.draw_marker)
        self.markerColorButton.clicked.connect( self.set_marker_color )
        self.markerLineWidthSpinBox.valueChanged.connect( self.draw_marker )

        # Connect color scale options
        self.color_scale_spin_box.valueChanged.connect(self.scale_color_by_spinbox)
        self.color_scale_slider.valueChanged.connect(self.scale_color_by_slider)
        self.colorScaleTypeSelectorComboBox.currentIndexChanged.connect(self.update_color_scale_type_map)

        # Connect color range options
        self.lineEditMaxColorRange.setValidator(QDoubleValidator())
        self.lineEditMaxColorRange.returnPressed.connect( self.update_color_scale_type_map )
        self.lineEditMinColorRange.setValidator(QDoubleValidator())
        self.lineEditMinColorRange.returnPressed.connect( self.update_color_scale_type_map )
        browserHelpers.add_icon_to_widget(self.go_color_range_button, 'gotoPlay.png')
        self.go_color_range_button.clicked.connect(self.update_color_scale_type_map)
        browserHelpers.add_icon_to_widget(self.reset_color_range_button, 'resetArrow.png')
        self.reset_color_range_button.clicked.connect(self.reset_color_range)

        # Connect color-map and interpolation options
        self.colorMapsDictionary = browserHelpers.add_colormaps_to_combobox(self.cmapCBox)
        self.cmapCBox.currentIndexChanged.connect(self.change_color_map_types)
        self.interpolation = browserHelpers.get_interpolation_dict()
        self.interpolationCBox.currentIndexChanged.connect(self.change_interpolation_method)

        # Connect go to options
        self.gotoXbox.setValidator(QIntValidator())
        self.gotoXbox.returnPressed.connect(self.goto_xy)
        self.gotoYbox.setValidator(QIntValidator())
        self.gotoYbox.returnPressed.connect(self.goto_xy)
        self.goto_button.clicked.connect(self.goto_xy)

        # connect navigation buttons
        self.navUpButton.pressed.connect(self.navigate_to_up)
        self.navDownButton.pressed.connect(self.navigate_to_down)
        self.navLeftButton.pressed.connect(self.navigate_to_left)
        self.navRightButton.pressed.connect(self.navigate_to_right)
        self.navUpRightButton.pressed.connect(self.navigate_to_up_right)
        self.navUpLeftButton.pressed.connect(self.navigate_to_up_left)
        self.navDownRightButton.pressed.connect(self.navigate_to_down_right)
        self.navDownLeftButton.pressed.connect(self.navigate_to_down_left)

        # Connect zoom options
        self.zoomInButton.pressed.connect(self.do_zoom_in)
        self.zoomOutButton.pressed.connect(self.do_zoom_out)
        self.zoomBinsSpinBox.valueChanged.connect(self.do_zoom_in_out)

        # Connect horizontal and vertical toolbar options
        self.vspaceIncrementButton.clicked.connect( lambda: self.set_vertical_space_to_figure('plus') )
        self.vspaceDecrementButton.clicked.connect( lambda: self.set_vertical_space_to_figure('minus') )
        self.hspaceIncrementButton.clicked.connect( lambda: self.set_horizontal_space_to_figure('plus') )
        self.hspaceDecrementButton.clicked.connect( lambda: self.set_horizontal_space_to_figure('minus') )
        self.vspaceLineEdit.setValidator(QDoubleValidator())
        self.hspaceLineEdit.setValidator(QDoubleValidator())
        self.fitPlotsTightButton.clicked.connect( self.fits_plots_on_page )

    def add_options_to_toolbar(self):
        """Place GUI widgets with icons inside toolbars.

        """
        # Set icons to zoom toolbar
        browserHelpers.add_icon_to_widget(self.zoomInButton, 'zoomIn.png')
        browserHelpers.add_icon_to_widget(self.zoomOutButton, 'zoomOut.png')

        # Set icons to navigation toolbar
        browserHelpers.add_icon_to_widget(self.navUpButton, 'upArrow.png')
        browserHelpers.add_icon_to_widget(self.navDownButton, 'downArrow.png')
        browserHelpers.add_icon_to_widget(self.navLeftButton, 'leftArrow.png')
        browserHelpers.add_icon_to_widget(self.navRightButton, 'rightArrow.png')
        browserHelpers.add_icon_to_widget(self.navUpLeftButton, 'leftUpArrow.png')
        browserHelpers.add_icon_to_widget(self.navDownRightButton, 'rightDownArrow.png')
        browserHelpers.add_icon_to_widget(self.navDownLeftButton, 'leftDownArrow.png')
        browserHelpers.add_icon_to_widget(self.navUpRightButton, 'rightUpArrow.png')

        # Place navigation toolbar
        self.navToolBar.addSeparator()
        self.navToolBar.addWidget(self.navUpButton)
        self.navToolBar.addWidget(self.navDownButton)
        self.navToolBar.addWidget(self.navLeftButton)
        self.navToolBar.addWidget(self.navRightButton)
        self.navToolBar.addWidget(self.navUpLeftButton)
        self.navToolBar.addWidget(self.navDownRightButton)
        self.navToolBar.addWidget(self.navDownLeftButton)
        self.navToolBar.addWidget(self.navUpRightButton)
        self.navToolBar.addSeparator()
        self.navToolBar.addWidget(self.navStepsLabel)
        self.navToolBar.addWidget(self.navStepsSpinBox)

        # Place zoom toolbar
        self.zoomToolBar.addSeparator()
        self.zoomToolBar.addWidget(self.zoomLabel)
        self.zoomToolBar.addWidget(self.zoomInButton)
        self.zoomToolBar.addWidget(self.zoomOutButton)
        self.zoomToolBar.addSeparator()
        self.zoomToolBar.addWidget(self.zoomBinsLabel)
        self.zoomToolBar.addWidget(self.zoomBinsSpinBox)
        self.zoomLabel.setAlignment(Qt.AlignCenter)
        self.zoomBinsLabel.setAlignment(Qt.AlignCenter)

        # Place "Go To" toolbar
        self.gotoToolBar.addSeparator()
        self.gotoToolBar.addWidget(self.gotoLabel)
        self.gotoToolBar.addWidget(self.gotoXLabel)
        self.gotoToolBar.addWidget(self.gotoXbox)
        self.gotoToolBar.addWidget(self.gotoYLabel)
        self.gotoToolBar.addWidget(self.gotoYbox)
        self.gotoToolBar.addWidget(self.gotoSpacerLabel)
        self.gotoToolBar.addWidget(self.gotoSpaceCBox)
        self.gotoToolBar.addWidget(self.goto_button)
        browserHelpers.add_icon_to_widget(self.goto_button, 'gotoPlay.png')
        self.gotoToolBar.addSeparator()
        self.gotoLabel.setAlignment(Qt.AlignCenter)
        self.gotoXLabel.setAlignment(Qt.AlignCenter)
        self.gotoYLabel.setAlignment(Qt.AlignCenter)

        # Set icons to horizontal-vertical spacer tool bar
        browserHelpers.add_icon_to_widget(self.vspaceIncrementButton, 'vspaceIncrement.png')
        browserHelpers.add_icon_to_widget(self.vspaceDecrementButton, 'vspaceDecrement.png')
        browserHelpers.add_icon_to_widget(self.hspaceIncrementButton, 'hspaceIncrement.png')
        browserHelpers.add_icon_to_widget(self.hspaceDecrementButton, 'hspaceDecrement.png')
        browserHelpers.add_icon_to_widget(self.fitPlotsTightButton, 'fitPlotsTight.png')


        # Place horizontal-vertical spacer tool bar
        self.spacerToolBar.addSeparator()
        self.spacerToolBar.addWidget(self.vspaceIncrementButton)
        self.spacerToolBar.addWidget(self.vspaceLineEdit)
        self.spacerToolBar.addWidget(self.vspaceDecrementButton)
        self.spacerToolBar.addSeparator()
        self.spacerToolBar.addWidget(self.hspaceIncrementButton)
        self.spacerToolBar.addWidget(self.hspaceLineEdit)
        self.spacerToolBar.addWidget(self.hspaceDecrementButton)
        self.spacerToolBar.addSeparator()
        self.spacerToolBar.addWidget(self.fitPlotsTightButton)
        self.spacerToolBar.addSeparator()

    def create_mpl_canvas(self):
        """ Create a new mpl Figure and FigCanvas object.
        """

        # Set dpi
        self.dpi = 100

        # New figure object of given figure size
        self.figure = mplFigure.Figure(dpi=self.dpi, facecolor='white', figsize=self.figsize)

        # New canvas object connect with figure object
        self.canvas = FigureCanvas(self.figure)

        # Clear the figure
        self.figure.clear()

        # Set focus of mouse
        self.canvas.setFocusPolicy( Qt.ClickFocus )
        self.canvas.setFocus()

        # we define the widget as expandable
        #self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        # Attach right click menu to canvas
        self.canvas.setContextMenuPolicy(Qt.CustomContextMenu)
        self.canvas.customContextMenuRequested.connect(self.ShowRightClickMenuCanvas)

        #notify the system of updated policy
        #self.canvas.updateGeometry()

    def embed_mpl_canvas_in_mainwindow(self):
        """ This function creates new canvas and embed it into main window.
            This function is created because new canvas is needed when page size or orientation is changed.
            Therefore, this function is called at initialization and when self.resize_mpl_figure() is called.
        """

        # If already figure and canvas, particularly during resizing, delete old figure and canvas
        if self.figure is not None:
            del self.figure
        if self.canvas is not None:
            del self.canvas

        # Create new canvas
        self.create_mpl_canvas()

        # Remove old widget and add new canvas widget to scroll area
        self.scrollAreaPlot.takeWidget()
        self.scrollAreaPlot.setWidgetResizable(True)
        self.scrollAreaPlot.setWidget(self.canvas)

    def connect_mouse_events_to_canvas(self):
        """connect to all the events we need
        """
        self.cidpress = self.canvas.mpl_connect('button_press_event', self.press_on_canvas)
        self.cidrelease = self.canvas.mpl_connect('button_release_event', self.release_on_canvas)
        self.cidmotion = self.canvas.mpl_connect('motion_notify_event', self.motion_on_canvas)
        self.cidscroll = self.canvas.mpl_connect('scroll_event', self.scroll_on_canvas)
        # Bind the 'pick' event for clicking on one of the bars
        #
        #self.cidpick = self.canvas.mpl_connect('pick_event', self.on_pick)

    def disconnect_mouse_events_from_canvas(self):
        """Disconnect to event
        Useful when resizing the figure size because new canvas is generated, so events need to connect again.
        """
        if hasattr(self, 'cidpress'):
            self.canvas.mpl_disconnect(self.cidpress)
        if hasattr(self, 'cidrelease'):
            self.canvas.mpl_disconnect(self.cidrelease)
        if hasattr(self, 'cidmotion'):
            self.canvas.mpl_disconnect(self.cidmotion)
        if hasattr(self, 'cidscroll'):
            self.canvas.mpl_disconnect(self.cidscroll)

    def resize_mpl_figure(self):
        """Resize the figure using figsize variable
        """

        prev_ActiveHiCmapAxis = self.ActiveHiCmapAxis

        # Disconnect all events
        self.disconnect_mouse_events_from_canvas()

        # Destroy old canvas, generate and embed new canvas
        self.embed_mpl_canvas_in_mainwindow()

        if self.hiCmapAxes is None:  return

        mid_idx = len(self.hiCmapAxes)
        i = 0
        while i < len(self.hiCmapAxes):
            if i == 0:
                ax = self.figure.add_subplot(1, mid_idx, i+1)
                ax.set_aspect('equal', 'box-forced')
            else:
                ax = self.figure.add_subplot(1, mid_idx, i+1, sharey=self.hiCmapAxes[0].ax)
                ax.set_aspect('equal', 'box-forced')
            self.hiCmapAxes[i].update_axes(ax) # Replotted both ccmap and any subplots here
            self.ActiveHiCmapAxis = i
            i = i+1

        # Connect mouse interaction with canvas
        self.connect_mouse_events_to_canvas()
        self.ActiveHiCmapAxis = prev_ActiveHiCmapAxis

        # Tight layout at start and then turn off it
        self.figure.set_tight_layout(True)
        self.canvas.draw()
        self.figure.set_tight_layout(False)

        # Get horizontal and vertical space between subplots
        self.get_horizontal_vertical_space_from_figure()

    def fits_plots_on_page(self):
        # Tight layout at start and then turn off it
        self.figure.set_tight_layout(True)
        self.canvas.draw()
        self.figure.set_tight_layout(False)
        self.get_horizontal_vertical_space_from_figure()

    def set_page_size(self):
        """Set page size when user select new page size from menu

        At first figsize is determined and then figure is resized.

        """

        old_figsize = self.figsize

        NewPageType = None
        # Checking whether a new menu is clicked and checked
        if self.actionPlotSizeA5.isChecked() and self.PageType != ePageType.A5 :
            NewPageType = ePageType.A5
        if self.actionPlotSizeA4.isChecked() and self.PageType != ePageType.A4 :
            NewPageType =  ePageType.A4
        if self.actionPlotSizeA3.isChecked() and self.PageType != ePageType.A3 :
            NewPageType = ePageType.A3
        if self.actionPlotSizeA2.isChecked() and self.PageType != ePageType.A2 :
            NewPageType =  ePageType.A2
        if self.actionPlotSizeA1.isChecked() and self.PageType != ePageType.A1 :
            NewPageType = ePageType.A1
        if self.actionPlotSizeA0.isChecked() and self.PageType != ePageType.A0 :
            NewPageType = ePageType.A0
        if self.actionPlotSizeLetter.isChecked() and self.PageType != ePageType.Letter :
            NewPageType = ePageType.Letter
        if self.actionPlotSizeCustom.isChecked() :
            NewPageType = ePageType.Custom

        # Assign new page type
        self.PageType = NewPageType

        # Check portrait or landscape is checked in page orientation menu
        bPortrait = self.actionPlotOrientationPortrait.isChecked()

        # Assign new page size depends on the clicked page type menu and page orientation menu
        if self.PageType == ePageType.A5:
            if bPortrait:
                self.figsize = PageSize[ePageType.A5]
            else:
                self.figsize = PageSize[ePageType.A5][::-1]

        if self.PageType == ePageType.A4:
            if bPortrait:
                self.figsize = PageSize[ePageType.A4]
            else:
                self.figsize = PageSize[ePageType.A4][::-1]

        if self.PageType == ePageType.A3:
            if bPortrait:
                self.figsize = PageSize[ePageType.A3]
            else:
                self.figsize = PageSize[ePageType.A3][::-1]

        if self.PageType == ePageType.A2:
            if bPortrait:
                self.figsize = PageSize[ePageType.A2]
            else:
                self.figsize = PageSize[ePageType.A2][::-1]

        if self.PageType == ePageType.A1:
            if bPortrait:
                self.figsize = PageSize[ePageType.A1]
            else:
                self.figsize = PageSize[ePageType.A1][::-1]

        if self.PageType == ePageType.A0:
            if bPortrait:
                self.figsize = PageSize[ePageType.A0]
            else:
                self.figsize = PageSize[ePageType.A0][::-1]

        if self.PageType == ePageType.Letter:
            if bPortrait:
                self.figsize = PageSize[ePageType.Letter]
            else:
                self.figsize = PageSize[ePageType.Letter][::-1]

        # Show a dialog box for custom page size
        if self.PageType == ePageType.Custom:
            Dialog = browserHelpers.DialogCustomPlotsize(self.figsize)
            result = Dialog.exec_()
            if result == Dialog.Accepted:
                self.figsize = tuple(Dialog.figsize)

            if Dialog.PageOrientation == 'Portrait':
                self.actionPlotOrientationPortrait.setChecked(True)
            else:
                self.actionPlotOrientationLandscape.setChecked(True)

            Dialog.close()
            del Dialog

        if old_figsize != self.figsize:
            self.resize_mpl_figure()

    def set_page_orientation(self):
        """ Set page Orientation when user select an option from menu
        """
        old_figsize = self.figsize
        NewPageOrientation = None

        # Checking whether a new menu is clicked and checked
        if self.actionPlotOrientationLandscape.isChecked() and self.PageOrientation != 'Landscape' :
            NewPageOrientation = 'Landscape'
        if self.actionPlotOrientationPortrait.isChecked() and self.PageOrientation != 'Portrait' :
            NewPageOrientation = 'Portrait'

        self.figsize = self.figsize[::-1]

        if old_figsize != self.figsize:
            self.resize_mpl_figure()

        return

    def get_horizontal_vertical_space_from_figure(self):
        """To get horizontal and vertical space between plots from the figure instance
        """
        self.verticalSpace = self.figure.subplotpars.hspace
        self.horizontalSpace = self.figure.subplotpars.wspace

    def set_vertical_space_to_figure(self, operator):
        """Increase or Decrease vertical spacing between plots
        """
        if operator == 'plus':
            new_verticalSpace =  self.verticalSpace + float( self.vspaceLineEdit.text() )
        else:
            new_verticalSpace =  self.verticalSpace - float( self.vspaceLineEdit.text() )

        if new_verticalSpace < 0:
            new_verticalSpace = 0.0

        self.figure.subplots_adjust(hspace=new_verticalSpace)

        # Update separately for all genomic plots
        for axs in self.hiCmapAxes:
            axs.verticalSpace = new_verticalSpace

        self.canvas.draw()
        self.get_horizontal_vertical_space_from_figure()

    def set_horizontal_space_to_figure(self, operator):
        """Increase or Decrease horizontal spacing between plots
        """
        if operator == 'plus':
            new_horizontalSpace =  self.horizontalSpace + float( self.hspaceLineEdit.text() )
        else:
            new_horizontalSpace =  self.horizontalSpace - float( self.hspaceLineEdit.text() )

        if new_horizontalSpace < 0:
            new_horizontalSpace = 0.0

        self.figure.subplots_adjust(wspace=new_horizontalSpace)
        self.canvas.draw()
        self.get_horizontal_vertical_space_from_figure()

    def add_new_ccmap_axes(self, pathToMapFile):
        ''' Add new ccmap axes

            This function add new hic-map on the canvas. To add new plots, previous plots are at first removed.
            Then, all previous plots are added including the new ones.

        '''

        if self.hiCmapAxes is None:
            ''' This condition run When first time ccmap is loaded.
                Because, there is no plots in figure object, no need to clear figure.
                Just add new single axes object and plot hic-map in it.
            '''

            # New axis
            ax = self.figure.add_subplot(1, 1, 1)
            self.hiCmapAxes = []

            # Initialize and append new ccmap axis object into hic-map list
            self.hiCmapAxes.append(CCMAPAXIS(0, ax))
            self.hiCmapAxes[0].title = os.path.splitext(pathToMapFile)[0]

            # Add this ccmap axis to the tree widget
            self.hiCmapAxes[0].set_tree_widget_item(self.axisTreeWidget)
            self.axisTreeWidget.addTopLevelItem(self.hiCmapAxes[0].treeWidgetItem)
            self.axisTreeWidget.setCurrentItem(self.hiCmapAxes[0].treeWidgetItem)


            # Add menu QAction for genomic data selector
            self.addMenuActionsForGenomicDatasetSelector(0)

            # set index of active ccmap axis
            self.ActiveHiCmapAxis = 0

            # Enable options to add other ccmap plot, only once done
            self.actionAdd_Hi_C_Map.setEnabled(True)

        else:
            ''' When more than one ccmaps are already present.

                At first, disconnect mouse event from canvas. Clear the figure.
                Generate new axes object using index (location) on figure object. Plots previous maps to newly generated axes object.
                Re-connect mouse event to new axes object. Add new ccmap axes object to tree widget and set it as an active map.

                I found that upon addition of more than one maps, map became rectangular.
                Therefore, I have to use set_aspect to enforce square maps.

                I also found that upon addition of more than two maps, first map entirely shrink and not visible.
                Therefore, I enforce lim on x-axis (set_xlim) to maintain size of the maps.

            '''

            # disconnect mouse event from canvas. Clear the figure.
            self.disconnect_mouse_events_from_canvas()
            self.canvas.figure.clear()

            # Generating new axes object by re-indexing old plots ( below while loop works for old maps)
            mid_idx = len(self.hiCmapAxes) + 1
            i = 0
            while i < len(self.hiCmapAxes):
                if i == 0:
                    ''' First axes-plot does not share its axis with any other plots.
                        Therefore, for first one, axes is generated differently.
                    '''
                    ax = self.figure.add_subplot(1, mid_idx, i+1)
                    ax.set_aspect('equal', 'box-forced')
                    ax.set_xlim(0, self.binsDisplayed)
                else:
                    ''' Other axes-plots share its axis with first plot.
                    '''
                    ax = self.figure.add_subplot(1, mid_idx, i+1, sharey=self.hiCmapAxes[0].ax)
                    ax.set_aspect('equal', 'box-forced')
                    ax.set_xlim(0, self.binsDisplayed)

                # Update axes_props object with new axes object. Also, plot is updated here
                self.hiCmapAxes[i].update_axes(ax)
                i = i+1

            # Add new axes object for new ccmap
            ax = self.figure.add_subplot(1, mid_idx, i+1, sharey=self.hiCmapAxes[0].ax)
            ax.set_aspect('equal', 'box-forced')
            self.hiCmapAxes.append(CCMAPAXIS(i, ax))

            # generate title
            self.hiCmapAxes[i].title = os.path.splitext(pathToMapFile)[0]

            # set index of active ccmap
            self.ActiveHiCmapAxis = i

            # Add new map to tree widget
            self.hiCmapAxes[i].set_tree_widget_item(self.axisTreeWidget)
            self.axisTreeWidget.addTopLevelItem(self.hiCmapAxes[i].treeWidgetItem)
            self.axisTreeWidget.setCurrentItem(self.hiCmapAxes[i].treeWidgetItem)

            # Add menu QAction for genomic data selector
            self.addMenuActionsForGenomicDatasetSelector(i)

        #self.figure.set_tight_layout(True)
        self.connect_mouse_events_to_canvas()

    def make_ccmap_active_on_click(self, item, col):
        ''' Set active ccmap axes object upon either clicking on list in tree widget or on plots.

        A ccmapaxis object store the current status of all the options. Therefore, it can be used to modify
        all GUI options. When user clicked on inactive ccmap plot, this function reads all the options from ccmapaxis
        object and modify these options on GUI.

        Note that, we only have to change options display. However, several GUI objects are connected with methods, and any change in GUI
        triggers these methods. Therefore, I have came with a new variable in ccmapaxis, which I use to prevent replotting.
        '''
        # If no ccmap or only one ccmap.
        if self.hiCmapAxes is None: return

        # Getting index of axis selected in tree widget
        hidx = None
        gidx = None
        for axis in self.hiCmapAxes:
            if axis.treeWidgetItem is self.axisTreeWidget.currentItem():
                hidx = axis.index
                break
            if axis.genomicPlotAxes is not None:
                for gax in axis.genomicPlotAxes:
                    if gax.treeWidgetItem is self.axisTreeWidget.currentItem():
                        hidx = axis.index
                        gidx = gax.index
                        break


        # Other conditions
        if hidx is None and gidx is None:    return
        if self.hiCmapAxes[hidx].image is None: return

        if len(self.hiCmapAxes) > 1 and self.ActiveHiCmapAxis != hidx:
            # Change index of active ccmap
            self.ActiveHiCmapAxis = hidx
            aidx = self.ActiveHiCmapAxis

            # Do not plot because only status of GUI options are changed here.
            self.hiCmapAxes[aidx].doNotPlot = True

            # Set resolution and mapunit
            self.resolutionLineEdit.setText(self.hiCmapAxes[aidx].resolution)
            self.mapUnitLineEdit.setText(self.hiCmapAxes[aidx].mapUnit)

            # Color scale status in qcombobox
            self.colorScaleTypeSelectorComboBox.setCurrentText(self.hiCmapAxes[aidx].colorScaleStatus)

            minvalue, maxvalue = self.hiCmapAxes[aidx].colorRangeMinValue, self.hiCmapAxes[aidx].colorRangeMaxValue
            if self.hiCmapAxes[aidx].colorScaleStatus == 'Logarithm of map':
                minvalue = np.log(self.hiCmapAxes[aidx].colorRangeMinValue)
                maxvalue = np.log(self.hiCmapAxes[aidx].colorRangeMaxValue)

            # Update color scale slider and spinbox
            self.color_scale_spin_box.setRange(minvalue, maxvalue)
            self.color_scale_spin_box.setSingleStep(self.hiCmapAxes[aidx].color_scale_steps[1]-self.hiCmapAxes[aidx].color_scale_steps[0])
            self.color_scale_spin_box.setValue(self.hiCmapAxes[aidx].color_scale_spinbox_value)
            self.color_scale_slider.setValue(self.hiCmapAxes[aidx].color_scale_slider_value)

            # Update color range line edit boxes
            self.lineEditMinColorRange.setText( str(self.hiCmapAxes[aidx].colorRangeMinValue) )
            self.lineEditMaxColorRange.setText( str(self.hiCmapAxes[aidx].colorRangeMaxValue) )

            # Update color map types and interpolation method in combo box
            if isinstance(self.hiCmapAxes[aidx].colormap, mpl.colors.LinearSegmentedColormap):
                self.cmapCBox.setCurrentText(self.hiCmapAxes[aidx].colormap.name)
            else:
                self.cmapCBox.setCurrentText(self.hiCmapAxes[aidx].colormap)
            self.interpolationCBox.setCurrentText(self.hiCmapAxes[aidx].interpolation)

            # Toggle replotting True
            self.hiCmapAxes[aidx].doNotPlot = False

            # Try to change map name list in combo box
            self.changeMapNamesInComboBox(aidx)

        # make active genomic dataset plot
        if gidx is not None:
            self.makeGenomicSupPlotOptionsActive(hidx, gidx)
        else:
            self.makeGenomicSupPlotOptionsInactive()

    def enable_gui_options(self):
        """To enable all options in tools.

        Initially all options are disabled to prevent any exception and error because user may interact with any GUI option without any ccmap loaded.
        This may throw exception and errors.

        This method should be called after ccmap is loaded first time. All options will be enabled for interactions.
        """

        # Tool bars
        self.navToolBar.setEnabled(True)
        self.zoomToolBar.setEnabled(True)
        self.gotoToolBar.setEnabled(True)
        self.spacerToolBar.setEnabled(True)

        # Reset map option
        self.reset_maps_button.setEnabled(True)

        # Marker option
        self.markerCBox.setEnabled(True)
        self.markerColorButton.setEnabled(True)
        self.markerLineWidthSpinBox.setEnabled(True)

        # Color scale options
        self.reset_color_range_button.setEnabled(True)
        self.go_color_range_button.setEnabled(True)
        self.lineEditMaxColorRange.setEnabled(True)
        self.lineEditMinColorRange.setEnabled(True)
        self.color_scale_spin_box.setEnabled(True)
        self.color_scale_slider.setEnabled(True)
        self.colorScaleTypeSelectorComboBox.setEnabled(True)
        self.interpolationCBox.setEnabled(True)
        self.cmapCBox.setEnabled(True)

    def init_gui_options_for_ccmap(self):
        """Set color related GUI options during initialization
        """
        aidx = self.ActiveHiCmapAxis
        self.hiCmapAxes[aidx].doNotPlot = True

        # Set resolution and mapunit
        self.resolutionLineEdit.setText(self.hiCmapAxes[aidx].resolution)
        self.mapUnitLineEdit.setText(self.hiCmapAxes[aidx].mapUnit)

        # Set color in GUI widget
        self.lineEditMaxColorRange.setText(str(self.hiCmapAxes[aidx].colorRangeMaxValue))
        self.lineEditMinColorRange.setText(str(self.hiCmapAxes[aidx].colorRangeMinValue))

        # Options
        self.color_scale_spin_box.setRange(self.hiCmapAxes[aidx].colorRangeMinValue, self.hiCmapAxes[aidx].colorRangeMaxValue)
        self.color_scale_spin_box.setSingleStep(self.hiCmapAxes[aidx].color_scale_steps[1]-self.hiCmapAxes[aidx].color_scale_steps[0])
        self.color_scale_spin_box.setValue(self.hiCmapAxes[aidx].colorRangeMaxValue)
        self.color_scale_slider.setValue(100)
        self.hiCmapAxes[aidx].set_color_spinbox_slider_values(self.color_scale_spin_box.value(), self.color_scale_slider.value())

        if isinstance(self.hiCmapAxes[aidx].colormap, mpl.colors.LinearSegmentedColormap):
            self.cmapCBox.setCurrentText(self.hiCmapAxes[aidx].colormap.name)
        else:
            self.cmapCBox.setCurrentText(self.hiCmapAxes[aidx].colormap)

        self.interpolationCBox.setCurrentText(self.hiCmapAxes[aidx].interpolation)
        self.colorScaleTypeSelectorComboBox.setCurrentText(self.hiCmapAxes[aidx].colorScaleStatus)

        self.hiCmapAxes[aidx].doNotPlot = False

        # Change contact map list
        self.changeMapNamesInComboBox(aidx)

    def open_map_pyobj(self):
        ''' To load new python object when clicked on Load Contact Map object

        This also add new ccmap axes instance to self.hiCmapAxes list.

        '''
        Load = True

        # This portion is from previous version. Now, it does not work as ccmap object is par of CCMAPAXIS class.
        # Here kept for reference
        if hasattr(self, 'ccmap'):
            if self.ccmap is not None:
                msg = '''Hi-C map object is already loaded...
                         Clcik OK to remove previously loaded Python MAP object
                      '''
                msgBox = QMessageBox()
                msgBox.setWindowTitle('Warning')
                msgBox.setText('Hi-C map object is already loaded...')
                msgBox.setInformativeText('Clcik OK to remove previously loaded Python MAP object')
                msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
                msgBox.setDefaultButton(QMessageBox.Cancel);
                pressed_button = msgBox.exec_()
                if pressed_button == QMessageBox.Ok:
                    Load = True
                else:
                    Load = False

        # A dialog box will be displayed to select a ccmap file
        file_choices = "Contact Map File (*.ccmap *.gcmap);;Hi-C contact map File (*.hicmap);;Chromosomal contact map file (*.ccmap);;Genome Contact Map File (*.gcmap)"
        path = QFileDialog.getOpenFileName(self, 'Load Map File', '', file_choices)

        # if user select a file
        if path[0]:
            file_extension = os.path.splitext(path[0])[1]
            self.status_bar.clearMessage()
            self.status_bar.showMessage('Loading file : {0}' .format(path[0]))

            fileType = None
            if file_extension == '.ccmap' or file_extension == '.hicmap':
                fileType = 'ccmap'
            elif file_extension == '.gcmap':
                fileType = 'gcmap'
            else:
                raise IOError ('File not end with gcmap, ccmap or hicmap. File may not be compatible.')

            dialog = None
            if fileType == 'gcmap':
                dialog = browserHelpers.GCMapSelectorDialog(path[0])
                dialog.exec_()
                if dialog.result() != QDialog.Accepted:
                    dialog.close()
                    del dialog
                    return

            # Add new CCMAPAXIS instance and respective matplotlib axes instance
            self.add_new_ccmap_axes(path[0])

            # If first time, enable all options
            if self.hiCmapAxes is not None:
                if len(self.hiCmapAxes) == 1:
                    self.enable_gui_options()

            # Initialize plot
            if fileType == 'gcmap':
                self.InitMapImage(path[0], fileType, mapName=dialog.mapName, resolution=dialog.resolution)
                dialog.close()
                del dialog
            else:
                self.InitMapImage(path[0], fileType)

            # Add treewidget item to tree widget
            self.axisTreeWidget.setCurrentItem(self.hiCmapAxes[self.ActiveHiCmapAxis].treeWidgetItem)
            self.status_bar.showMessage('{0} loaded.' .format(path[0]))

    def InitMapImage(self, path, fileType, mapName=None, resolution=None):
        """Initialize map plot for first time
        """
        aidx = self.ActiveHiCmapAxis

        #Initialize the class with ccmap file
        self.hiCmapAxes[aidx].initialize(path, fileType, mapName=mapName, resolution=resolution, filesOpened=self.filesOpened)

        # Make matrix readable
        if isinstance(self.hiCmapAxes[aidx].ccmap, gmlib.ccmap.CCMAP):
            self.hiCmapAxes[aidx].ccmap.make_readable()

        # To get new limit
        if self.hiCmapAxes[aidx].ccmap.shape[0] > self.binsDisplayed and self.hiCmapAxes[aidx].ccmap.shape[1] > self.binsDisplayed:
            xrange = [1, self.binsDisplayed]
            yrange = [1, self.binsDisplayed]

        else:
            if self.hiCmapAxes[aidx].ccmap.shape[0] <= self.binsDisplayed:
                self.binsDisplayed = self.hiCmapAxes[aidx].ccmap.shape[0]

            if self.hiCmapAxes[aidx].ccmap.shape[1] <= self.binsDisplayed:
                self.binsDisplayed = self.hiCmapAxes[aidx].ccmap.shape[1]

            xrange = [1, self.binsDisplayed]
            yrange = [1, self.binsDisplayed]

        # Init GUI options
        self.init_gui_options_for_ccmap()

        # Try if dynamics resolutions interchange is possible
        self.tryEnableInterchangeResolutions()

        self.hiCmapAxes[aidx].rangeXY = (xrange, yrange)

        if len(self.hiCmapAxes) > 1:
            self.ResetMap()
        else:
            if self.interchangeableResolutions is not None:
                self.currentResolutionIndex = self.interchangeableResolutions.index( self.hiCmapAxes[aidx].ccmap.resolution )
            self.fits_plots_on_page()

    def tryEnableInterchangeResolutions(self):
        """ Try to enable dynamic interchange in resolution.

        If it is successful, self.interchangeableResolutions contains list of resolutions
        """
        allGCMap = True
        binsizes = []

        # First check if file format is ccmap or gcmap
        for hax in self.hiCmapAxes:
            if hax.fileType == 'ccmap':
                allGCMap = False
            else:
                binsizes.append( hax.ccmap.binsizes )

        if allGCMap:
            binSizeMatched = True
            common = None

            # If only one map is loaded
            if len(binsizes) == 1:
                common = binsizes[0]

            # If more than one, determine common available resolutions among all maps
            for i in range(len(binsizes) - 1 ):
                b1 = set( binsizes[i] )
                b2 = set( binsizes[i+1] )
                c = set.intersection( b1 , b2 )
                if not c:
                    binSizeMatched = False
                else:
                    if common is None:
                        common = c
                    else:
                        common = set.intersection( common, set(c) )
                        if common is None:
                            binSizeMatched = False

            if binSizeMatched:
                self.interchangeableResolutions = list( map(gmlib.util.binsizeToResolution, sorted( list(common) ) ) )
            else:
                self.interchangeableResolutions = None

        else:
            self.interchangeableResolutions = None

    def changeMapNamesInComboBox(self, aidx):
        """ Change map names in combo box depending on the active map
        """
        self.contactMapNameCBox.blockSignals(True)

        if self.hiCmapAxes[aidx].interchangeableCMapNames is not None:
            self.contactMapNameCBox.setEnabled(True)
            self.contactMapNameCBox.clear()
            self.contactMapNameCBox.addItems( self.hiCmapAxes[aidx].interchangeableCMapNames )
            self.contactMapNameCBox.setCurrentText( self.hiCmapAxes[aidx].ccmap.groupName )
        else:
            self.contactMapNameCBox.clear()
            self.contactMapNameCBox.setEnabled(False)

        self.contactMapNameCBox.blockSignals(False)

    def changeMapNames(self, newMapName):
        """ Changes chromosome in browser
        It at first check that change is possible than it changes the maps
        """
        aidx = self.ActiveHiCmapAxis

        success = True

        # Check if current resolution is present in the data
        if self.hiCmapAxes[aidx].fileType == 'gcmap':
            success = self.hiCmapAxes[aidx].ccmap.checkMapExist(mapName=newMapName, resolution=self.hiCmapAxes[aidx].resolution)
            if not success:
                msg = 'For {0}, {1} resolution map is not available.\nChange the resolution before changing map.'.format(newMapName, self.hiCmapAxes[aidx].resolution)
                browserHelpers.showWarningMessageBox(msg, self)

        # First check is all data is available in genomic subplot
        if self.hiCmapAxes[aidx].genomicPlotAxes is not None and success:
            for gax in self.hiCmapAxes[aidx].genomicPlotAxes:
                if not gax.changeDataByName(self, newMapName, change=False):
                    success = False

        # In case if same data name is not found in genomic plot. Do not change data name
        if not success:
            self.changeMapNamesInComboBox(aidx)
            return

        # Block signals temporarily
        self.contactMapNameCBox.blockSignals(True)

        # Change contact map data
        self.InitMapImage(self.hiCmapAxes[aidx].ccmap.hdf5.filename, self.hiCmapAxes[aidx].fileType,
                        mapName=newMapName, resolution=self.hiCmapAxes[aidx].resolution)
        self.reset_color_range() # Maintain color scaling type and reset scale

        # Change genomic subplot data
        if self.hiCmapAxes[aidx].genomicPlotAxes is not None:
            for gax in self.hiCmapAxes[aidx].genomicPlotAxes:
                gax.changeDataByName(self, newMapName, change=True)
                self.resetGenomicDataYScaleLimits(hidx=aidx, gidx=gax.index)
                gax.updatePlot()

        self.contactMapNameCBox.blockSignals(False)

    def update_color_scale_type_map(self):
        if self.hiCmapAxes is None: return
        aidx = self.ActiveHiCmapAxis

        if self.hiCmapAxes[aidx].ccmap.minvalue < 0.0 and self.colorScaleTypeSelectorComboBox.currentText() != 'Change color linearly':
            self.colorScaleTypeSelectorComboBox.blockSignals(True)
            self.colorScaleTypeSelectorComboBox.setCurrentText('Change color linearly')
            self.colorScaleTypeSelectorComboBox.blockSignals(False)
            return

        # Take options value from GUI
        if not self.hiCmapAxes[aidx].doNotPlot:
            minvalue = float(self.lineEditMinColorRange.text())
            maxvalue = float(self.lineEditMaxColorRange.text())

            # If user put wrong values
            if minvalue >= maxvalue:
                msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'Minimum value is equal or larger than maximum value.', QMessageBox.Ok, self)
                msgBox.exec_()
                msgBox.close()
                self.lineEditMinColorRange.setText(str(self.hiCmapAxes[aidx].colorRangeMinValue))
                self.lineEditMaxColorRange.setText(str(self.hiCmapAxes[aidx].colorRangeMaxValue))
                return

            # If user put minimum value as zero
            if minvalue == 0 and self.hiCmapAxes[aidx].colorScaleStatus == 'Logarithm of map':
                msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'Minimum value cannot be zero in case of logarithm', QMessageBox.Ok, self)
                msgBox.exec_()
                msgBox.close()
                self.lineEditMinColorRange.setText(str(self.hiCmapAxes[aidx].colorRangeMinValue))
                self.lineEditMaxColorRange.setText(str(self.hiCmapAxes[aidx].colorRangeMaxValue))
                return

            self.hiCmapAxes[aidx].colorRangeMinValue = minvalue
            self.hiCmapAxes[aidx].colorRangeMaxValue = maxvalue


        self.hiCmapAxes[aidx].colorScaleStatus = self.colorScaleTypeSelectorComboBox.currentText()
        self.canvas.draw()

        # Change options in GUI
        self.hiCmapAxes[aidx].doNotPlot = True # Plot already updated above, so do not re-plot when spinbox and slider value is changing
        self.color_scale_spin_box.setRange(self.hiCmapAxes[aidx].color_scale_steps[0], self.hiCmapAxes[aidx].color_scale_steps[-1])
        self.color_scale_spin_box.setSingleStep(self.hiCmapAxes[aidx].color_scale_steps[1]-self.hiCmapAxes[aidx].color_scale_steps[0])
        self.color_scale_spin_box.setValue(self.hiCmapAxes[aidx].color_scale_spinbox_value)
        self.color_scale_slider.setValue(100)
        self.hiCmapAxes[aidx].doNotPlot = False

    def reset_color_range(self):
        aidx = self.ActiveHiCmapAxis
        self.hiCmapAxes[aidx].colorRangeMinValue = self.hiCmapAxes[aidx].ccmap.minvalue
        self.hiCmapAxes[aidx].colorRangeMaxValue = self.hiCmapAxes[aidx].ccmap.maxvalue
        self.lineEditMaxColorRange.setText(str(self.hiCmapAxes[aidx].colorRangeMaxValue))
        self.lineEditMinColorRange.setText(str(self.hiCmapAxes[aidx].colorRangeMinValue))
        self.hiCmapAxes[aidx].scaleMinMaxForResolutionInterchange = False
        self.update_color_scale_type_map()

    def scale_color_by_spinbox (self):
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes is None: return
        if self.hiCmapAxes[aidx].image is None: return

        idx = ( np.abs(self.hiCmapAxes[aidx].color_scale_steps - self.color_scale_spin_box.value() )).argmin()
        self.color_scale_slider.setValue(idx)

    def scale_color_by_slider(self):
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes is None: return
        if self.hiCmapAxes[aidx].image is None: return


        if self.hiCmapAxes[aidx].ccmap is not None:
            color_value = self.hiCmapAxes[aidx].color_scale_steps[self.color_scale_slider.value()]
            self.color_scale_spin_box.setValue(color_value)
            self.hiCmapAxes[aidx].color_scale_spinbox_value = color_value
            self.hiCmapAxes[aidx].color_scale_slider_value = self.color_scale_slider.value()
            self.canvas.draw()

    def change_color_map_types(self, index):
        """ Change colormap from combo box
        """
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        self.hiCmapAxes[aidx].colormap =  self.colorMapsDictionary[int(index)]
        self.canvas.draw()

    def showUserDefinedColormapDialog(self):
        """ Open dialog to generate new user defined colormap
        """
        colormapDialog = browserHelpers.DialogUserColorMap()
        colormapDialog.setColorMapList(self.colorMapsDictionary)
        colormapDialog.accepted.connect( self.userDefinedColorMapDialogAccepted )
        colormapDialog.exec_()

    def userDefinedColorMapDialogAccepted(self):
        """ Get user defined colormap from dialog
        """
        dialog = self.sender()
        if dialog.resultColorInfo is None:
            return

        # Get color information
        colorInfo = dialog.resultColorInfo
        idx = self.cmapCBox.findText(colorInfo['name'], Qt.MatchExactly)
        if idx == -1:
            browserHelpers.add_external_colormap_to_combobox(self.cmapCBox, self.colorMapsDictionary, colorInfo)
        else:
            self.colorMapsDictionary[idx] = browserHelpers.colorInfoToSegmentDataColorMap(colorInfo)
            browserHelpers.change_colormap_icon_to_combobox(self.cmapCBox, idx, colorInfo)
            dialog.close()

        # update plot in browser
        idx = self.cmapCBox.findText(colorInfo['name'], Qt.MatchExactly)
        if self.cmapCBox.currentIndex() == idx:
            self.change_color_map_types(idx)
        else:
            self.cmapCBox.setCurrentIndex(idx)

    def change_interpolation_method(self, index):
        """ Change interpolation of plot
        """
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return
        self.hiCmapAxes[aidx].interpolation = self.interpolation[int(index)]
        self.canvas.draw()

    def draw_marker(self):
        """ Draw the marker on the plots
        """
        if self.hiCmapAxes is None: return
        axis_list = []
        for i in range(len(self.hiCmapAxes)):
            axis_list.append(self.hiCmapAxes[i].ax)
            if self.hiCmapAxes[i].genomicPlotAxes is not None:
                for j in range(len(self.hiCmapAxes[i].genomicPlotAxes)):
                    axis_list.append(self.hiCmapAxes[i].genomicPlotAxes[j].ax)

        if self.marker is not None:
            self.marker.disconnect()
            self.marker = None

        # Any change in Spinbox, change linewidth here, kept here as no need to write an additional function
        self.markerLineWidth = self.markerLineWidthSpinBox.value()

        if self.markerCBox.currentIndex() == 1:
            self.marker = malWidgets.MultiCursor(self.canvas, axis_list, horizOn=True, color=self.markerColor, linewidth=self.markerLineWidth)
        elif self.markerCBox.currentIndex() == 2:
            self.marker = malWidgets.MultiCursor(self.canvas, axis_list, horizOn=True, vertOn=False, color=self.markerColor, linewidth=self.markerLineWidth)
        elif self.markerCBox.currentIndex() == 3:
            self.marker = malWidgets.MultiCursor(self.canvas, axis_list, horizOn=False, color=self.markerColor, linewidth=self.markerLineWidth)
        else:
            pass

        self.canvas.draw()

    def set_marker_color(self):
        """ Set the color of marker
        """

        # Convert RGB tuple to QColor, QColor range from 0 to 255 while matplotlib color range from 0 to 1.
        qcolor = QColor.fromRgb( int(self.markerColor[0]*255), int(self.markerColor[1]*255), int(self.markerColor[2]*255), 255 )

        # QColorDialog open here
        colorDialog = QColorDialog(qcolor, self)
        colorDialog.setWindowTitle("Choose a Color")
        colorDialog.exec_()

        if colorDialog.result() == QDialog.Accepted:
            # Get new color from the user
            pickedColor = colorDialog.selectedColor()
            self.markerColor = (pickedColor.red()/255, pickedColor.green()/255, pickedColor.blue()/255)

            # Draw marker with new color
            self.draw_marker()

            # Set background color of button
            style = 'background-color: rgb({0}, {1}, {2});'.format(pickedColor.red(), pickedColor.green(), pickedColor.blue())
            self.markerColorButton.setStyleSheet(style)

        del colorDialog
        del qcolor

    def goto_xy(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        xword = str( self.gotoXbox.text() )
        yword = str( self.gotoYbox.text() )

        try:
            x = int(xword.replace(',', ''))
        except:
            return

        try:
            y = int(yword.replace(',', ''))
        except:
            return

        if not xword:
            msg = 'Please enter a value for X'
            QMessageBox.warning(self, "Warning!", msg)
            return
        if not yword:
            msg = 'Please enter a value for Y'
            QMessageBox.warning(self, "Warning!", msg)
            return

        mapUnitSize = gmlib.util.resolutionToBinsize(self.hiCmapAxes[aidx].mapUnit)

        if self.gotoSpaceCBox.currentIndex() == 1:
            x = int( x / self.hiCmapAxes[aidx].ccmap.binsize)
            y = int( y / self.hiCmapAxes[aidx].ccmap.binsize)

        else:
            x = int( x * mapUnitSize / self.hiCmapAxes[aidx].ccmap.binsize)
            y = int( y * mapUnitSize / self.hiCmapAxes[aidx].ccmap.binsize)


        '''
        xdiff = self.hiCmapAxes[aidx].xrange[1] - self.hiCmapAxes[aidx].xrange[0]
        ydiff = self.hiCmapAxes[aidx].yrange[1] - self.hiCmapAxes[aidx].yrange[0]

        div = self.hiCmapAxes[aidx].ccmap.binsize

        xrange = [ (x - self.hiCmapAxes[aidx].ccmap.xticks[0]/div) - (xdiff/2), (x - self.hiCmapAxes[aidx].ccmap.xticks[0]/div) + (xdiff/2) ]
        yrange = [ (y - self.hiCmapAxes[aidx].ccmap.yticks[0]/div) - (ydiff/2), (y - self.hiCmapAxes[aidx].ccmap.yticks[0]/div) + (ydiff/2) ]
        '''

        xrange = [x - int(self.binsDisplayed/2), x +  int(self.binsDisplayed/2) ]
        yrange = [y - int(self.binsDisplayed/2), y +  int(self.binsDisplayed/2) ]

        if xrange[0] < 0:
            xrange = [ 0, self.binsDisplayed ]

        if xrange[1] > self.hiCmapAxes[aidx].ccmap.shape[0]:
            xrange = [self.hiCmapAxes[aidx].ccmap.shape[0] - self.binsDisplayed, self.hiCmapAxes[aidx].ccmap.shape[0]]

        if yrange[0] < 0:
            yrange = [ 0, self.binsDisplayed ]

        if yrange[1] > self.hiCmapAxes[aidx].ccmap.shape[1]:
            yrange = [self.hiCmapAxes[aidx].ccmap.shape[1] - self.binsDisplayed, self.hiCmapAxes[aidx].ccmap.shape[1]]

        draw = True
        if xrange[0] < 0 or yrange[0] < 0:
            draw = False

        for i in range(len(self.hiCmapAxes)):
            if self.hiCmapAxes[i].ccmap.shape[1] < yrange[1] or self.hiCmapAxes[i].ccmap.shape[0] < xrange[1]:
                draw = False
        if draw:
            for i in range(len(self.hiCmapAxes)):
                self.hiCmapAxes[i].rangeXY = (xrange, yrange)
            self.canvas.draw()

    def do_zoom_in(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        input_steps = 10
        actual_steps = input_steps

        yrange = [self.hiCmapAxes[aidx].yrange[0]+actual_steps, self.hiCmapAxes[aidx].yrange[1]-actual_steps]
        xrange = [self.hiCmapAxes[aidx].xrange[0]+actual_steps, self.hiCmapAxes[aidx].xrange[1]-actual_steps]

        self.changeMapSize(xrange, yrange)

    def do_zoom_out(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        input_steps = 10
        actual_steps = input_steps

        yrange = [self.hiCmapAxes[self.ActiveHiCmapAxis].yrange[0]-actual_steps, self.hiCmapAxes[self.ActiveHiCmapAxis].yrange[1]+actual_steps]
        xrange = [self.hiCmapAxes[self.ActiveHiCmapAxis].xrange[0]-actual_steps, self.hiCmapAxes[self.ActiveHiCmapAxis].xrange[1]+actual_steps]

        self.changeMapSize(xrange, yrange)

    def do_zoom_in_out(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        yrange = None
        xrange = None

        yrange = self.hiCmapAxes[self.ActiveHiCmapAxis].yrange
        xrange = self.hiCmapAxes[self.ActiveHiCmapAxis].xrange

        # To prevent this function, when other functions try to set value in spinbox
        if self.binsDisplayed == int(self.zoomBinsSpinBox.value()): return

        curr_mid_x = int( (xrange[1] + xrange[0])/2 )
        curr_mid_y = int( (yrange[1] + yrange[0])/2 )

        target_mid_value = int(int(self.zoomBinsSpinBox.value())/2)

        #print(xrange, yrange)
        xrange = [curr_mid_x-target_mid_value, curr_mid_x+target_mid_value]
        yrange = [curr_mid_y-target_mid_value, curr_mid_y+target_mid_value]
        #print(xrange, yrange)

        self.changeMapSize(xrange, yrange)

    def navigate_to_up(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        input_steps = self.navStepsSpinBox.value()
        actual_steps = input_steps

        # To check if new range value is greater than the maximum in any of the ccmap
        # Here, input steps is not used because it will allow to calculate minimum steps in case of several ccmap
        for i in range(len(self.hiCmapAxes)):
            yrange = [self.hiCmapAxes[i].yrange[0]+actual_steps, self.hiCmapAxes[i].yrange[1]+actual_steps]

            # If yrange maximum is larger than ccmap shape, change the increment (steps) so that yrange maximum equal to ccmap shape.
            if yrange[1] >= self.hiCmapAxes[i].ccmap.shape[1]:
                actual_steps = self.hiCmapAxes[i].ccmap.shape[1] - self.hiCmapAxes[i].yrange[1]

        if actual_steps > 0:
            for i in range(len(self.hiCmapAxes)):
                yrange = [self.hiCmapAxes[i].yrange[0]+actual_steps, self.hiCmapAxes[i].yrange[1]+actual_steps]
                #print(self.hiCmapAxes[i].ccmap.shape, yrange, actual_steps)
                self.hiCmapAxes[i].rangeXY = (self.hiCmapAxes[i].xrange, yrange)
            self.canvas.draw()
        actual_steps = input_steps

    def navigate_to_down(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        input_steps = self.navStepsSpinBox.value()
        actual_steps = input_steps

        # To check if new range value is greater than the minimum in any of the ccmap
        for i in range(len(self.hiCmapAxes)):
            yrange = [self.hiCmapAxes[i].yrange[0]-actual_steps, self.hiCmapAxes[i].yrange[1]-actual_steps]
            if yrange[0] < 0:
                actual_steps = self.hiCmapAxes[i].yrange[0]

        if actual_steps > 0:
            for i in range(len(self.hiCmapAxes)):
                yrange = [self.hiCmapAxes[i].yrange[0]-actual_steps, self.hiCmapAxes[i].yrange[1]-actual_steps]
                #print(self.hiCmapAxes[i].ccmap.shape, yrange, actual_steps)
                self.hiCmapAxes[i].rangeXY = (self.hiCmapAxes[i].xrange, yrange)
            self.canvas.draw()

        actual_steps = input_steps

    def navigate_to_right(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        input_steps = self.navStepsSpinBox.value()
        actual_steps = input_steps

        # To check if new range value is greater than the maximum in any of the ccmap
        for i in range(len(self.hiCmapAxes)):
            xrange = [self.hiCmapAxes[i].xrange[0]+actual_steps, self.hiCmapAxes[i].xrange[1]+actual_steps]
            if xrange[1] >= self.hiCmapAxes[i].ccmap.shape[0]:
                actual_steps = self.hiCmapAxes[i].ccmap.shape[0] - self.hiCmapAxes[i].xrange[1]

        # If range is within the maximum limits and difference between current position and maximum limits is more than input steps, redraw all the plots
        if actual_steps > 0:
            for i in range(len(self.hiCmapAxes)):
                xrange = [self.hiCmapAxes[i].xrange[0]+actual_steps, self.hiCmapAxes[i].xrange[1]+actual_steps]
                #print(self.hiCmapAxes[i].ccmap.shape, xrange, actual_steps)
                self.hiCmapAxes[i].rangeXY = (xrange, self.hiCmapAxes[i].yrange)
            self.canvas.draw()

        actual_steps = input_steps

    def navigate_to_left(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        input_steps = self.navStepsSpinBox.value()
        actual_steps = input_steps

        # To check if new range value is greater than the minimum in any of the ccmap
        for i in range(len(self.hiCmapAxes)):
            xrange = [self.hiCmapAxes[i].xrange[0]-actual_steps, self.hiCmapAxes[i].xrange[1]-actual_steps]
            if xrange[0] < 0:
                actual_steps = self.hiCmapAxes[i].xrange[0]

        # If range is within the minimum limits and difference between current position and minimum limits is more than input steps, redraw all the plots
        if actual_steps > 0:
            for i in range(len(self.hiCmapAxes)):
                xrange = [self.hiCmapAxes[i].xrange[0]-actual_steps, self.hiCmapAxes[i].xrange[1]-actual_steps]
                #print(self.hiCmapAxes[i].ccmap.shape, xrange, actual_steps)
                self.hiCmapAxes[i].rangeXY = (xrange, self.hiCmapAxes[i].yrange)
            self.canvas.draw()

        actual_steps = input_steps

    def navigate_to_up_right(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        input_steps = self.navStepsSpinBox.value()
        actual_steps = input_steps

        # along Y-axis
        # To check if new range value is greater than the maximum in any of the ccmap
        # Here, input steps is not used because it will allow to calculate minimum steps in case of several ccmap
        for i in range(len(self.hiCmapAxes)):
            yrange = [self.hiCmapAxes[i].yrange[0]+actual_steps, self.hiCmapAxes[i].yrange[1]+actual_steps]

            # If yrange maximum is larger than ccmap shape, change the increment (steps) so that yrange maximum equal to ccmap shape.
            if yrange[1] >= self.hiCmapAxes[i].ccmap.shape[1]:
                actual_steps = self.hiCmapAxes[i].ccmap.shape[1] - self.hiCmapAxes[i].yrange[1]

        # Along X-axis
        # To check if new range value is greater than the maximum in any of the ccmap
        for i in range(len(self.hiCmapAxes)):
            xrange = [self.hiCmapAxes[i].xrange[0]+actual_steps, self.hiCmapAxes[i].xrange[1]+actual_steps]
            if xrange[1] >= self.hiCmapAxes[i].ccmap.shape[0]:
                actual_steps = self.hiCmapAxes[i].ccmap.shape[0] - self.hiCmapAxes[i].xrange[1]


        if actual_steps > 0:
            for i in range(len(self.hiCmapAxes)):
                yrange = [self.hiCmapAxes[i].yrange[0]+actual_steps, self.hiCmapAxes[i].yrange[1]+actual_steps]
                xrange = [self.hiCmapAxes[i].xrange[0]+actual_steps, self.hiCmapAxes[i].xrange[1]+actual_steps]
                #print(self.hiCmapAxes[i].ccmap.shape, xrange, yrange, actual_steps)
                self.hiCmapAxes[i].rangeXY = (xrange, yrange)
            self.canvas.draw()
        actual_steps = input_steps

    def navigate_to_up_left(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        input_steps = self.navStepsSpinBox.value()
        actual_steps = input_steps

        # along Y-axis
        # To check if new range value is greater than the maximum in any of the ccmap
        # Here, input steps is not used because it will allow to calculate minimum steps in case of several ccmap
        for i in range(len(self.hiCmapAxes)):
            yrange = [self.hiCmapAxes[i].yrange[0]+actual_steps, self.hiCmapAxes[i].yrange[1]+actual_steps]

            # If yrange maximum is larger than ccmap shape, change the increment (steps) so that yrange maximum equal to ccmap shape.
            if yrange[1] >= self.hiCmapAxes[i].ccmap.shape[1]:
                actual_steps = self.hiCmapAxes[i].ccmap.shape[1] - self.hiCmapAxes[i].yrange[1]

        # Along X-axis
        # To check if new range value is greater than the minimum in any of the ccmap
        for i in range(len(self.hiCmapAxes)):
            xrange = [self.hiCmapAxes[i].xrange[0]-actual_steps, self.hiCmapAxes[i].xrange[1]-actual_steps]
            if xrange[0] < 0:
                actual_steps = self.hiCmapAxes[i].xrange[0]

        if actual_steps > 0:
            for i in range(len(self.hiCmapAxes)):
                yrange = [self.hiCmapAxes[i].yrange[0]+actual_steps, self.hiCmapAxes[i].yrange[1]+actual_steps]
                xrange = [self.hiCmapAxes[i].xrange[0]-actual_steps, self.hiCmapAxes[i].xrange[1]-actual_steps]
                #print(self.hiCmapAxes[i].ccmap.shape, xrange, yrange, actual_steps)
                self.hiCmapAxes[i].rangeXY = (xrange, yrange)
            self.canvas.draw()
        actual_steps = input_steps

    def navigate_to_down_right(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        input_steps = self.navStepsSpinBox.value()
        actual_steps = input_steps

        # along Y-axis
        # To check if new range value is greater than the maximum in any of the ccmap
        # Here, input steps is not used because it will allow to calculate minimum steps in case of several ccmap
        # To check if new range value is greater than the minimum in any of the ccmap
        for i in range(len(self.hiCmapAxes)):
            yrange = [self.hiCmapAxes[i].yrange[0]-actual_steps, self.hiCmapAxes[i].yrange[1]-actual_steps]
            if yrange[0] < 0:
                actual_steps = self.hiCmapAxes[i].yrange[0]

        # Along X-axis
        # To check if new range value is greater than the maximum in any of the ccmap
        for i in range(len(self.hiCmapAxes)):
            xrange = [self.hiCmapAxes[i].xrange[0]+actual_steps, self.hiCmapAxes[i].xrange[1]+actual_steps]
            if xrange[1] >= self.hiCmapAxes[i].ccmap.shape[0]:
                actual_steps = self.hiCmapAxes[i].ccmap.shape[0] - self.hiCmapAxes[i].xrange[1]


        if actual_steps > 0:
            for i in range(len(self.hiCmapAxes)):
                yrange = [self.hiCmapAxes[i].yrange[0]-actual_steps, self.hiCmapAxes[i].yrange[1]-actual_steps]
                xrange = [self.hiCmapAxes[i].xrange[0]+actual_steps, self.hiCmapAxes[i].xrange[1]+actual_steps]
                #print(self.hiCmapAxes[i].ccmap.shape, xrange, yrange, actual_steps)
                self.hiCmapAxes[i].rangeXY = (xrange, yrange)
            self.canvas.draw()
        actual_steps = input_steps

    def navigate_to_down_left(self):
        if self.hiCmapAxes is None:  return
        aidx = self.ActiveHiCmapAxis
        if self.hiCmapAxes[aidx].image is None: return

        input_steps = self.navStepsSpinBox.value()
        actual_steps = input_steps

        # along Y-axis
        # To check if new range value is greater than the maximum in any of the ccmap
        # Here, input steps is not used because it will allow to calculate minimum steps in case of several ccmap
        # To check if new range value is greater than the minimum in any of the ccmap
        for i in range(len(self.hiCmapAxes)):
            yrange = [self.hiCmapAxes[i].yrange[0]-actual_steps, self.hiCmapAxes[i].yrange[1]-actual_steps]
            if yrange[0] < 0:
                actual_steps = self.hiCmapAxes[i].yrange[0]

        # Along X-axis
        # To check if new range value is greater than the maximum in any of the ccmap
        for i in range(len(self.hiCmapAxes)):
            xrange = [self.hiCmapAxes[i].xrange[0]-actual_steps, self.hiCmapAxes[i].xrange[1]-actual_steps]
            if xrange[0] < 0:
                actual_steps = self.hiCmapAxes[i].xrange[0]


        if actual_steps > 0:
            for i in range(len(self.hiCmapAxes)):
                yrange = [self.hiCmapAxes[i].yrange[0]-actual_steps, self.hiCmapAxes[i].yrange[1]-actual_steps]
                xrange = [self.hiCmapAxes[i].xrange[0]-actual_steps, self.hiCmapAxes[i].xrange[1]-actual_steps]
                #print(self.hiCmapAxes[i].ccmap.shape, xrange, yrange, actual_steps)
                self.hiCmapAxes[i].rangeXY = (xrange, yrange)
            self.canvas.draw()
        actual_steps = input_steps

    def ResetMap(self):
        self.binsDisplayed = 1000

        xrange = []
        yrange = []
        for i in range(len(self.hiCmapAxes)):
            if self.binsDisplayed > self.hiCmapAxes[i].ccmap.shape[0]:
                self.binsDisplayed = self.hiCmapAxes[i].ccmap.shape[0]
            if self.binsDisplayed > self.hiCmapAxes[i].ccmap.shape[1]:
                self.binsDisplayed = self.hiCmapAxes[i].ccmap.shape[1]

            if self.interchangeableResolutions is not None:
                idx = self.interchangeableResolutions.index( self.hiCmapAxes[i].ccmap.resolution )
                self.currentResolutionIndex =  max( idx , self.currentResolutionIndex )

        for aidx in range(len(self.hiCmapAxes)):
            if self.hiCmapAxes[aidx].ccmap.shape[0] > self.binsDisplayed and self.hiCmapAxes[aidx].ccmap.shape[1] > self.binsDisplayed:
                xrange = [1, self.binsDisplayed]
                yrange = [1, self.binsDisplayed]

            else:
                if self.hiCmapAxes[aidx].ccmap.shape[0] <= self.binsDisplayed:
                    self.binsDisplayed = self.hiCmapAxes[aidx].ccmap.shape[0]

                if self.hiCmapAxes[aidx].ccmap.shape[1] <= self.binsDisplayed:
                    self.binsDisplayed = self.hiCmapAxes[aidx].ccmap.shape[1]

                xrange = [1, self.binsDisplayed]
                yrange = [1, self.binsDisplayed]

            if self.interchangeableResolutions is not None:
                self.hiCmapAxes[aidx].ccmap.changeResolution(self.interchangeableResolutions[self.currentResolutionIndex])

            self.hiCmapAxes[aidx].rangeXY = (xrange, yrange)

        self.zoomBinsSpinBox.setValue(self.binsDisplayed)
        self.canvas.draw()

        # Tight layout at start and then turn off it
        self.figure.set_tight_layout(True)
        self.canvas.draw()
        self.figure.set_tight_layout(False)

        # Get horizontal and vertical space between subplots
        self.get_horizontal_vertical_space_from_figure()

    def show_value_at_cursor(self, event):
        ''' Show value at mouse pointer
        '''
        if self.hiCmapAxes is None:  return
        obj = self.get_ccmap_obj_axis_under_mouse_mpl(event)
        mi = obj.xrange[0] + int(event.xdata)
        mj = obj.yrange[0] + int(event.ydata)
        self.label_mij = obj.ax.text(event.xdata, event.ydata, '{0}, {1}, {2:.3f}' .format(obj.xticklabels[mi], obj.yticklabels[mj],obj.ccmap.matrix[mi][mj]) )
        self.canvas.draw()

    def press_on_canvas(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if self.hiCmapAxes is None:  return

        # Nothing do on right click
        if event.button == 3:   return

        # Get the ccmap axis which is under the mouse when it is pressed
        temp_ccmap_axis = self.get_ccmap_obj_axis_under_mouse_mpl(event)
        temp_subplot_axis = self.get_subplot_obj_axis_under_mouse_mpl(event)

        # if no ccmap axis under mouse when pressed, no need to do anything
        if temp_ccmap_axis is None and temp_subplot_axis is None: return

        # If mouse button is pressed then set the current ccmap axis as active
        if temp_ccmap_axis is not None:
            self.axisTreeWidget.setCurrentItem(temp_ccmap_axis.treeWidgetItem)

        if temp_subplot_axis is not None:
            self.axisTreeWidget.setCurrentItem(temp_subplot_axis.treeWidgetItem)

        self.press_on_plot = int(event.xdata), int(event.ydata)
        #self.show_value_at_cursor(event)
        self.canvas.draw()

    def motion_on_canvas(self, event):
        'on motion we will move the rect if the mouse is over us'

        if self.hiCmapAxes is None:  return
        temp_ccmap_axis = self.get_ccmap_obj_axis_under_mouse_mpl(event)
        temp_subplot_axes = self.get_subplot_obj_axis_under_mouse_mpl(event)

        if temp_ccmap_axis is None and temp_subplot_axes is None:
            self.status_bar.clearMessage()

            # Hide marker from plots
            if self.marker is not None:
                self.marker.visible = False
                self.marker.onmove(event)
                self.canvas.draw()
            return


        if temp_subplot_axes is not None:
            self.motion_on_canvas_on_subplots(event, temp_subplot_axes)
            return

        if event.inaxes is not temp_ccmap_axis.ax:  return


        # Display marker on plots
        if self.marker is not None:
            self.marker.visible = True
            self.marker.onmove(event)

        curr_pos_x = int(event.xdata)
        curr_pos_y = int(event.ydata)
        real_xcoor = temp_ccmap_axis.xticklabels[curr_pos_x + temp_ccmap_axis.xrange[0]]
        real_ycoor = temp_ccmap_axis.yticklabels[curr_pos_y + temp_ccmap_axis.yrange[0]]
        map_value = temp_ccmap_axis.ccmap.matrix[curr_pos_x + temp_ccmap_axis.xrange[0]][curr_pos_y + temp_ccmap_axis.yrange[0]]

        if temp_ccmap_axis.colorScaleStatus == 'Logarithm of map':
            if map_value != 0:
                map_value = np.log(map_value)

        mapUnitSize = gmlib.util.resolutionToBinsize(temp_ccmap_axis.mapUnit)
        self.status_bar.showMessage('|  Real X = {0:,}  |  Real Y = {1:,}  |  X = {2}  |  Y = {3}  |  Value = {4:.6f}  |' .format(
                                    int(mapUnitSize*real_xcoor), int(mapUnitSize*real_ycoor),
                                    real_xcoor, real_ycoor, map_value))

        # If mouse button is not pressed on the plots, return
        if self.press_on_plot is None: return

        # When mouse button is press and on hold (drag motion)
        xpress, ypress = self.press_on_plot
        dx = xpress - curr_pos_x
        dy = ypress - curr_pos_y
        #print('curr_pos_x=%f, xpress=%f, event.xdata=%f, dx=%f, curr_pos_x+dx=%f'%(curr_pos_x, xpress, event.xdata, dx, curr_pos_x+dx))
        #print(dx, dy)
        #self.rect.set_x(x0+dx)
        #self.rect.set_y(y0+dy)
        map_xrange = temp_ccmap_axis.xrange.copy()
        map_yrange = temp_ccmap_axis.yrange.copy()

        map_xrange = [ temp_ccmap_axis.xrange[0] + dx, temp_ccmap_axis.xrange[1] + dx ]
        map_yrange = [ temp_ccmap_axis.yrange[0] + dy, temp_ccmap_axis.yrange[1] + dy ]

        if map_xrange[0] < 0:
            return

        if map_xrange[1] > temp_ccmap_axis.ccmap.shape[0]:
            return

        if map_yrange[0] < 0:
            return

        if map_yrange[1] > temp_ccmap_axis.ccmap.shape[1]:
            return

        if (temp_ccmap_axis.xrange[1] - temp_ccmap_axis.xrange[0]) != (map_xrange[1] - map_xrange[0]):
            return

        if (temp_ccmap_axis.yrange[1] - temp_ccmap_axis.yrange[0]) != (map_yrange[1] - map_yrange[0]):
            return

        #print('on move', map_xrange, map_yrange, self.hiCmapAxes[0].ccmap.shape, self.hiCmapAxes[1].ccmap.shape, self.hiCmapAxes[0].xrange, self.hiCmapAxes[1].xrange)
        #print(self.press_on_plot)
        draw = True
        for i in range(len(self.hiCmapAxes)):
            if self.hiCmapAxes[i].ccmap.shape[1] < map_yrange[1] or self.hiCmapAxes[i].ccmap.shape[0] < map_xrange[1]:
                draw = False
        if draw:
            for i in range(len(self.hiCmapAxes)):
                self.hiCmapAxes[i].rangeXY = (map_xrange, map_yrange)
            self.canvas.draw()

        self.press_on_plot = int(event.xdata), int(event.ydata)

    def motion_on_canvas_on_subplots(self, event, temp_subplot_axes):
        'on motion we will move the rect if the mouse is over us'

        if event.inaxes is not temp_subplot_axes.ax:    return

        # Display marker on plots
        if self.marker is not None:
            self.marker.visible = True
            self.marker.onmove(event)

        curr_pos_x = int(event.xdata)
        curr_pos_y = int(event.ydata)
        real_xcoor = temp_subplot_axes.hiCmapAxis.xticklabels[curr_pos_x + temp_subplot_axes.hiCmapAxis.xrange[0]]
        value = temp_subplot_axes.dataArray[curr_pos_x + temp_subplot_axes.hiCmapAxis.xrange[0]]

        self.status_bar.showMessage('|  Real X = {0:,}  |  X = {1}  |  Y = {2:.6f} |' .format(
                                    int(temp_subplot_axes.hiCmapAxis.ccmap.binsize*real_xcoor), real_xcoor, value))

        # If mouse button is not pressed on the plots, return
        if self.press_on_plot is None: return

        # When mouse button is press and on hold (drag motion)
        xpress, ypress = self.press_on_plot
        dx = xpress - curr_pos_x

        map_xrange = [ temp_subplot_axes.hiCmapAxis.xrange[0] + dx, temp_subplot_axes.hiCmapAxis.xrange[1] + dx ]

        if map_xrange[0] < 0:
            return

        if map_xrange[1] > temp_subplot_axes.hiCmapAxis.ccmap.shape[0]:
            return

        if (temp_subplot_axes.hiCmapAxis.xrange[1] - temp_subplot_axes.hiCmapAxis.xrange[0]) != (map_xrange[1] - map_xrange[0]):
            return

        #print('on move', map_xrange, map_yrange, self.hiCmapAxes[0].ccmap.shape, self.hiCmapAxes[1].ccmap.shape, self.hiCmapAxes[0].xrange, self.hiCmapAxes[1].xrange)
        #print(self.press_on_plot)
        draw = True
        for i in range(len(self.hiCmapAxes)):
            if self.hiCmapAxes[i].ccmap.shape[0] < map_xrange[1]:
                draw = False
        if draw:
            for i in range(len(self.hiCmapAxes)):
                self.hiCmapAxes[i].rangeXY = (map_xrange, self.hiCmapAxes[i].yrange)
            self.canvas.draw()

        self.press_on_plot = int(event.xdata), int(event.ydata)

    def release_on_canvas(self, event):
        'on release we reset the press data'

        self.press_on_plot = None
        try:
            self.label_mij.remove()
        except:
            pass
        self.canvas.draw()

    def scroll_on_canvas(self, event):
        """ When user scroll on map
            * Zoom in
            * Zoom out
            * Do nothing
        """

        if self.hiCmapAxes is None:  return
        temp_ccmap_axis = self.get_ccmap_obj_axis_under_mouse_mpl(event)
        temp_subplot_axes = self.get_subplot_obj_axis_under_mouse_mpl(event)

        if temp_ccmap_axis is None and temp_subplot_axes is None: return

        if temp_subplot_axes is not None:
            temp_ccmap_axis = temp_subplot_axes.hiCmapAxis

        if event.inaxes != temp_ccmap_axis.ax and event.inaxes != temp_subplot_axes.ax: return

        if self.press_on_plot is not None: return

        map_xrange = temp_ccmap_axis.xrange.copy()
        map_yrange = temp_ccmap_axis.yrange.copy()

        if (event.step > 0):
            map_xrange[0] = temp_ccmap_axis.xrange[0] + 10
            map_xrange[1] = temp_ccmap_axis.xrange[1] - 10
            map_yrange[0] = temp_ccmap_axis.yrange[0] + 10
            map_yrange[1] = temp_ccmap_axis.yrange[1] - 10

        if (event.step < 0):
            map_xrange[0] = temp_ccmap_axis.xrange[0] - 10
            map_xrange[1] = temp_ccmap_axis.xrange[1] + 10
            map_yrange[0] = temp_ccmap_axis.yrange[0] - 10
            map_yrange[1] = temp_ccmap_axis.yrange[1] + 10

        # If new map size is fine, try to change the map size on the interface
        self.changeMapSize(map_xrange, map_yrange)

    def changeMapSize(self, map_xrange, map_yrange):
        """ It changes the size of the map on zoom in and zoom out.

        It also changes the resolution if possible
        """
        if map_xrange[0] < 0:
            map_xrange[0] = 0
        if map_yrange[0] < 0:
            map_yrange[0] = 0
        if map_xrange[1] > self.hiCmapAxes[self.ActiveHiCmapAxis].ccmap.shape[0]:
            map_xrange[1] = self.hiCmapAxes[self.ActiveHiCmapAxis].ccmap.shape[0]
        if map_yrange[1] > self.hiCmapAxes[self.ActiveHiCmapAxis].ccmap.shape[1]:
            map_yrange[1] = self.hiCmapAxes[self.ActiveHiCmapAxis].ccmap.shape[1]

        # Always make a square map. Below two conditions are implemented to restrain the map as a square.
        # Never be an rectangle
        # Idea here is that to decrease the longer axis and make it equal to shorter axis.
        # This way, zoom in and out at the end of any axis works fine.
        xlength = map_xrange[1] - map_xrange[0]
        ylength = map_yrange[1] - map_yrange[0]

        if xlength < ylength:
            if map_yrange[0] == 0:
                map_yrange[1] = map_yrange[1] - (ylength - xlength)
            elif map_yrange[1] ==  self.hiCmapAxes[self.ActiveHiCmapAxis].ccmap.shape[1]:
                map_yrange[0] = map_yrange[0] + (ylength - xlength)
            else:
                if (ylength - xlength) % 2 == 0:
                    map_yrange[0] = map_yrange[0] + int( (ylength - xlength)/2 )
                    map_yrange[1] = map_yrange[1] - int( (ylength - xlength)/2 )
                else:
                    map_yrange[0] = map_yrange[0] + math.ceil(  (ylength - xlength)/2 )
                    map_yrange[1] = map_yrange[1] - math.floor( (ylength - xlength)/2 )

        if ylength < xlength:
            if map_xrange[0] == 0:
                map_xrange[1] = map_xrange[1] - (xlength - ylength)
            elif map_xrange[1] ==  self.hiCmapAxes[self.ActiveHiCmapAxis].ccmap.shape[0]:
                map_xrange[0] = map_xrange[0] + (xlength - ylength)
            else:
                if (xlength - ylength) % 2 == 0:
                    map_xrange[0] = map_xrange[0] + int( (xlength - ylength)/2 )
                    map_xrange[1] = map_xrange[1] - int( (xlength - ylength)/2 )
                else:
                    map_xrange[0] = map_xrange[0] + math.ceil(  (xlength - ylength)/2 )
                    map_xrange[1] = map_xrange[1] - math.floor( (xlength - ylength)/2 )

        # Do not go below this zoom
        if (map_xrange[1] - map_xrange[0]) < self.zoomBinsSpinBox.minimum():
            return

        # Do not go above this zoom
        if (map_xrange[1] - map_xrange[0]) > self.zoomBinsSpinBox.maximum():
            return

        # To check if new xrange and yrange is available for all maps. If not do not update the plot.
        draw = True
        for i in range(len(self.hiCmapAxes)):
            if self.hiCmapAxes[i].ccmap.shape[1] < map_yrange[1] or self.hiCmapAxes[i].ccmap.shape[0] < map_xrange[1]:
                draw = False


        if draw:
            # If resolution interchange is allowed
            if self.interchangeableResolutions is not None:
                new_map_xrange, new_map_yrange = self.tryChangingResolutionsAll(map_xrange, map_yrange)

                # In case when change in resolution not possible, continue with original zoom
                if new_map_xrange is not None and new_map_xrange is not None:
                    map_xrange = new_map_xrange
                    map_yrange = new_map_xrange

            for i in range(len(self.hiCmapAxes)):
                self.hiCmapAxes[i].rangeXY = (map_xrange, map_yrange)
            self.binsDisplayed = map_xrange[1] - map_xrange[0]
            self.zoomBinsSpinBox.setValue(int(self.binsDisplayed))
            self.canvas.draw()

    def tryChangingResolutionsAll(self, map_xrange, map_yrange):
        """ Try to change resolutions of all maps given the threshold of map shape

        Also, change resolution of genomic dataset.

        """

        newNumberBins = map_xrange[1] - map_xrange[0]
        resolutionChanged = True
        # do fine -> coarse
        if newNumberBins > 1001:
            newResIndex = self.currentResolutionIndex + 1
            resolutionChanged = False

            # Determine if new resolution is available for all maps
            if newResIndex < len(self.interchangeableResolutions):
                resolutionChanged = True
                for i in range(len(self.hiCmapAxes)):
                    if gmlib.util.resolutionToBinsize(self.interchangeableResolutions[newResIndex]) > self.hiCmapAxes[i].ccmap.binsizes[-1]:
                        resolutionChanged = False
                        break

                    # Check if resolution changing is possible for genomic dataset
                    if self.hiCmapAxes[i].genomicPlotAxes is not None:
                        for gax in self.hiCmapAxes[i].genomicPlotAxes:
                            if not gax.changeResolution(self.interchangeableResolutions[newResIndex], change=False):
                                resolutionChanged = False
                                if gax.showWarningNoChangeResolution:
                                    msgBox = QMessageBox(QMessageBox.Warning, 'Warning',
                                                         'Not able to change resolution for {0} of {1}.\n\n Continuing without change in resolution...'
                                                         .format(self.axisTreeWidget.itemWidget(gax.treeWidgetItem, 0).text(),
                                                                 self.axisTreeWidget.itemWidget(self.hiCmapAxes[i].treeWidgetItem, 0).text()),
                                                                 QMessageBox.Ok, self)
                                    msgBox.exec_()
                                    msgBox.close()
                                    gax.showWarningNoChangeResolution = False

                                break

            # If new resolution is available for new maps, change the resolution and respective xrange and yrange
            if resolutionChanged:
                self.currentResolutionIndex = self.currentResolutionIndex + 1
                for i in range(len(self.hiCmapAxes)):
                    self.hiCmapAxes[i].ccmap.changeResolution(self.interchangeableResolutions[self.currentResolutionIndex])
                    self.hiCmapAxes[i].updatePropsForResolution()                   # Update some properties
                    self.hiCmapAxes[i].scaleMinMaxForResolutionInterchange = True
                    self.resolutionLineEdit.setText(self.hiCmapAxes[self.ActiveHiCmapAxis].resolution)  # Update resolution on GUI

                    # Change resolution of genomic dataset
                    if self.hiCmapAxes[i].genomicPlotAxes is not None:
                        for gax in self.hiCmapAxes[i].genomicPlotAxes:
                            gax.changeResolution(self.interchangeableResolutions[self.currentResolutionIndex])

                map_xrange[0] = int( np.floor(map_xrange[0]/2))
                map_xrange[1] = int( np.floor(map_xrange[1]/2))
                map_yrange[0] = int( np.floor(map_yrange[0]/2))
                map_yrange[1] = int( np.floor(map_yrange[1]/2))

        # do coarse -> fine
        if newNumberBins < 500 :
            newResIndex = self.currentResolutionIndex - 1
            resolutionChanged = False

            # Determine if new resolution is available for all maps
            if newResIndex >= 0:
                resolutionChanged = True
                for i in range(len(self.hiCmapAxes)):
                    if gmlib.util.resolutionToBinsize(self.interchangeableResolutions[newResIndex]) < self.hiCmapAxes[i].ccmap.binsizes[0]:
                        resolutionChanged = False
                        break

                    # Check if resolution changing is possible for genomic dataset
                    if self.hiCmapAxes[i].genomicPlotAxes is not None:
                        for gax in self.hiCmapAxes[i].genomicPlotAxes:
                            if not gax.changeResolution(self.interchangeableResolutions[newResIndex], change=False):
                                resolutionChanged = False
                                if gax.showWarningNoChangeResolution:
                                    msgBox = QMessageBox(QMessageBox.Warning, 'Warning',
                                                         'Not able to change resolution for {0} of {1}.\n\nContinuing without change in resolution...'
                                                         .format(self.axisTreeWidget.itemWidget(gax.treeWidgetItem, 0).text(),
                                                                 self.axisTreeWidget.itemWidget(self.hiCmapAxes[i].treeWidgetItem, 0).text()),
                                                                 QMessageBox.Ok, self)
                                    msgBox.exec_()
                                    msgBox.close()
                                    gax.showWarningNoChangeResolution = False
                                break

            # If new resolution is available for new maps, change the resolution and respective xrange and yrange
            if resolutionChanged:
                self.currentResolutionIndex = self.currentResolutionIndex - 1

                for i in range(len(self.hiCmapAxes)):
                    self.hiCmapAxes[i].ccmap.changeResolution(self.interchangeableResolutions[self.currentResolutionIndex])
                    self.hiCmapAxes[i].updatePropsForResolution()                   # Update some properties
                    self.hiCmapAxes[i].scaleMinMaxForResolutionInterchange = True
                    self.resolutionLineEdit.setText(self.hiCmapAxes[self.ActiveHiCmapAxis].resolution)  # Update resolution on GUI

                    # Change resolution of genomic dataset
                    if self.hiCmapAxes[i].genomicPlotAxes is not None:
                        for gax in self.hiCmapAxes[i].genomicPlotAxes:
                            gax.changeResolution(self.interchangeableResolutions[self.currentResolutionIndex])
                            gax.showWarningNoChangeResolution = True

                map_xrange[0] = map_xrange[0] * 2
                map_xrange[1] = map_xrange[1] * 2
                map_yrange[0] = map_yrange[0] * 2
                map_yrange[1] = map_yrange[1] * 2

        if resolutionChanged:
            return map_xrange, map_yrange
        else:
            return None, None

    def get_ccmap_obj_axis_under_mouse_mpl(self, event):
        if self.hiCmapAxes is not None:
            for i in range(len(self.hiCmapAxes)):
                contains, attrib = self.hiCmapAxes[i].ax.contains(event)
                if contains:
                    return self.hiCmapAxes[i]

    def get_ccmap_obj_axis_under_mouse_qt(self, qpoint):
        if self.hiCmapAxes is not None:
            cw, ch = self.canvas.get_width_height()
            for i in range(len(self.hiCmapAxes)):
                bbox = self.hiCmapAxes[i].ax.get_position()
                if bbox.contains(qpoint.x()/cw, 1-(qpoint.y()/ch)):
                    return self.hiCmapAxes[i]

    def save_plot(self):
        """ Save plot as a image
        """
        # Get list of all available formats
        formatGroups = self.canvas.get_supported_filetypes_grouped()

        # Make full list of formats
        file_choices = 'Image formats ('
        for desc in formatGroups:
            for imgFormat in formatGroups[desc]:
                file_choices += ' *.{0}'.format(imgFormat)
        file_choices += ')'

        # Make list format by each category
        for desc in formatGroups:
            file_choices += ';;{0} ('.format(desc)
            for imgFormat in formatGroups[desc]:
                file_choices += ' *.{0}'.format(imgFormat)
            file_choices += ')'

        # Just added for any formats
        file_choices += ";;All formats (*.*)"

        path = QFileDialog.getSaveFileName(self, 'Save image file', '', file_choices)

        if path[0]:
            self.figure.savefig(path[0], dpi=300)
            self.status_bar.showMessage('Saved to %s' % path[0], 2000)

    def on_about(self):
        """ Show about dialog
        """
        aboutDialog = browserHelpers.aboutBrowserDialog(self)
        aboutDialog.show()

    def printPlot(self):
        self.canvas.figure.set_dpi(300)
        self.canvas.draw()
        printer = QPrinter()
        printer.setResolution(300)
        printWidget= QPrintDialog(printer, self.canvas)

        if(printWidget.exec_() == QDialog.Accepted):
            p = self.canvas.grab(QRect(0, 0, int(300*self.figsize[0]), int(300*self.figsize[1])))
            painter = QPainter(printer)
            painter.setRenderHint(QPainter.HighQualityAntialiasing)
            painter.drawPixmap(0, 0, p)
            painter.end()

        self.canvas.figure.set_dpi(100)
        self.canvas.draw()

    def ShowRightClickMenuCanvas(self, qpoint):
        if self.hiCmapAxes is None:  return

        temp_ccmap_axis = self.get_ccmap_obj_axis_under_mouse_qt(qpoint)
        temp_subplot_axis = self.get_subplot_obj_axis_under_mouse_qt(qpoint)

        if temp_ccmap_axis is not None:
            self.menuRightClick = browserHelpers.menuRightClick()
            selectedItem = self.menuRightClick.exec_(QCursor.pos())
            if selectedItem == 0:   return
            self.processRightClickItem(selectedItem, temp_ccmap_axis)
            self.menuRightClick = None

        if temp_subplot_axis is not None:
            self.menuRightClick = browserHelpers.menuRightClick()
            selectedItem = self.menuRightClick.exec_(QCursor.pos())
            if selectedItem == 0:   return
            self.processRightClickItem(selectedItem, temp_subplot_axis)
            self.menuRightClick = None

    def ShowRightClickMenuTreeWidget(self, qpoint):
        if self.hiCmapAxes is None:  return
        self.menuRightClick = browserHelpers.menuRightClick()
        selectedItem = self.menuRightClick.exec_(QCursor.pos())
        if selectedItem == 0:   return

        temp_ccmap_axis = self.hiCmapAxes[self.ActiveHiCmapAxis]

        # Determine if genomic dataset subplot is clicked
        active_genomic_axis = None
        if temp_ccmap_axis.genomicPlotAxes is not None:
            for gax in temp_ccmap_axis.genomicPlotAxes:
                if gax.treeWidgetItem is self.axisTreeWidget.currentItem():
                    active_genomic_axis = gax
                    break

        if active_genomic_axis is None:
            self.processRightClickItem(selectedItem, temp_ccmap_axis)
        else:
            self.processRightClickItem(selectedItem, active_genomic_axis)

        self.menuRightClick = None

    def processRightClickItem(self, selectedItem, temp_axis):
        if self.menuRightClick is None: return

        for i in range(4):
            if self.menuRightClick.actionList[i] is selectedItem:

                if self.DialogAxisProps is not None:
                    if self.DialogAxisProps.result() == QDialog.Rejected or self.DialogAxisProps.result() == QDialog.Accepted:
                        self.DialogAxisProps.close()
                        self.DialogAxisProps = None

                if self.DialogAxisProps is None:
                    self.DialogAxisProps = browserHelpers.DialogAxisProps(self.canvas, i, temp_axis.ax, axesProps=temp_axis.axes_props)
                    self.DialogAxisProps.show()
                else:
                    if temp_axis.ax is not self.DialogAxisProps.axes:
                        self.DialogAxisProps.update_axes(temp_axis.ax, axesProps=temp_axis.axes_props)

                    self.DialogAxisProps.tabWidgetAxisProps.setCurrentIndex(i)
                    self.DialogAxisProps.activateWindow()

    def addMenuActionsForGenomicDatasetSelector(self, idx):
        """To add menu-actions to main menu to select genomic dataset for particular ccmap
        """
        #Add menu for respective ccmap and store it for future reference
        if self.menuActionsToAddGenomicDataset is None:

            # Remove a qaction
            self.menuAddGenomicDataset.removeAction(self.actionDummy)

            # Initialize list
            self.menuActionsToAddGenomicDataset = []
            self.menuAddGenomicDataset.setEnabled(True)

        # Add action menu
        self.hiCmapAxes[idx].qActionToAddGenomicDataset = self.menuAddGenomicDataset.addAction(self.hiCmapAxes[idx].title)
        self.menuActionsToAddGenomicDataset.append(self.hiCmapAxes[idx].qActionToAddGenomicDataset)

        # Connect action menu to open genomic file
        self.hiCmapAxes[idx].qActionToAddGenomicDataset.triggered.connect( lambda: self.openGenomicDatasetFile(self.hiCmapAxes[idx]) )

    def analysisCorrelateMaps(self):

        msg = 'This interface is in development.\nPlease use "gcMapExplorer corrBWcmaps" to calculate correlation.'
        guiHelpers.showWarningMessageBox(msg, self)

        '''
        self.corrDialog = browserHelpers.DialogCorrelationBetweenMaps()

        if self.hiCmapAxes is not None:
            for i in range(len(self.hiCmapAxes)):
                self.corrDialog.addHicmaps(self.hiCmapAxes[i].ccmap, self.hiCmapAxes[i].title)

        self.corrDialog.exec_()

        self.corrDialog.close()
        del self.corrDialog
        self.corrDialog = None
        '''

class GenomicDataPlotAxis:

    def __init__(self, index, hiCmapAxis):
        self.index = index             # Its index in self.hiCmapAxes[idx].genomicPlotAxes
        self.hiCmapAxis = hiCmapAxis   # reference to parent in list self.hiCmapAxes[idx]
        self.ax = None                 # matplotlib axes instance
        self.plot_type = None          # name of plot types
        self.plot = None               # matplotlib plot instance
        self.plotColor = (0, 0, 1)     # Color of plot, default Blue
        self.plotLineWidth = 0.5       # Vertical Line Width
        self.axes_props = None         # AxesProperties instance
        self.yticks = None             # yticks
        self.yticksFormat = 'plain'    # Ytick label formats - plain or scientific
        self.yticksDecimals = 2        # Decimals in y ticks formatting
        self.percentileValues = dict() # Store values at given percentile


        # Genomic dataset stuffs
        self.hdf5Hand = None           # HDF5Handler instance
        self.shownDataset = None       # Dictionary keyword for chromosome, resolution and data coarsening methods -> can be used to access array from hdf5Hand

        self.converterDialog = None    # dialog converter for bigwig, wig and bed file

        self.txtFileHand = None        # Genomic text file handler

        self.dataArray = None          # array presently plotted
        self.plotLocation = None       # Location of plot: top or bottom

        self.showWarningNoChangeResolution = True        # To show warning message when resolution cannot be shown. It should be only shown once

        self.treeWidgetItem = None      # tree widget item instance

        self.yscaleSlider = None        # value in slider
        self.yscaleSpinbox = None       # value in spinbox

        self._ylimit = None             # upper and lower limit in current plot
        self._yScaleSteps = None        # 100 steps between maximum and minimum allowed limit along y-axis.
        self._xlabel = None             # ylabel instance

    @property
    def ylimit(self):
        """Upper and lower limit of current plot.
        Permanent limits are stored in yScaleSteps, which can be changed by user using range box.
        """
        return self._ylimit

    @ylimit.setter
    def ylimit(self, value):
        """Upper and lower limit of current plot.
        Permanent limits are stored in yScaleSteps, which can be changed by user using range box.

        This is the only function through which plots can be updated.
        """
        if value[0] == value[1]:
            return

        if not self.hiCmapAxis.doNotPlot:
            self._ylimit = [ value[0], value[1] ]

            # Setting yticks based on ylimits, only change yticks when y-limits are changed
            ydiff = value[0] - value[1]
            yticks = np.linspace(self.ylimit[0], self.ylimit[1], 100)

            # Set tick label according to precision
            self.yticksDecimals = gmlib.util.locate_significant_digit_after_decimal(ydiff)
            if self.yticksDecimals > 3:
                self.yticksFormatStyle = 'sci'
                self.yticks = np.around(yticks, decimals=self.yticksDecimals+1)
            else:
                self.yticks = np.around(yticks, decimals=self.yticksDecimals+1)

            self.updatePlot()

    @property
    def yScaleSteps(self):
        return self._yScaleSteps

    @yScaleSteps.setter
    def yScaleSteps(self, value):
        """To change to upper and lower limit of y-scale slider and spinbox
           1) Change this variable -> change spinbox value => slider will update, ylimits will update and plot will be updated automatically.
           2) To reset scale: Change this variable to None -> change spinbox value => slider will update, ylimits will update and plot will be updated automatically.
        """
        if value is None:
            # Get minimum and maximum. Maximum is taken at 95 percentile in case minimum is not negative
            minimum = np.amin(self.dataArray)
            if minimum < 0:
                maximum = np.amax(self.dataArray)
            else:
                maximum = np.percentile(self.dataArray[self.dataArray[:] != 0 ], 95)
            self._yScaleSteps = np.linspace(minimum, maximum, 100)
        else:
            self._yScaleSteps = np.linspace(value[0], value[1], 100)

    @property
    def xlabel(self):
        return self._xlabel

    @xlabel.setter
    def xlabel(self, value):
        """ Also set xlabel in plot
        """
        self._xlabel = value
        if self.axes_props is not None:
            self.axes_props.xLabel['Text'] = value
            self.axes_props.set_to_axes()
        else:
            self.ax.set_xlabel(value, fontsize=12)


    def selectGenomicDataByH5Conversion(self, browser, filename, chromosome, filesOpened):
        """ This can be used to load a dataset by converting any other format
        like bigwig, wig or bed format to H5 format.

        NOTE: temporary H5 file should be deleted when browser is closed.

        """

        if self.converterDialog is None:
            self.converterDialog = h5Converter.DialogOther1DFormatLoader(browser)
            self.converterDialog.setModal(True)
            self.converterDialog.connect2Hide()
            self.converterDialog.connectToBrowserDataLoader(browser, self)
            self.converterDialog.setInputFile(filename)
            self.converterDialog.enableFileIndexing()
            self.converterDialog.setChromName(chromosome)
            self.converterDialog.setH5Name()
            self.converterDialog.show()

    def selectGenomicDataHdf5ByDialogBox(self, filename, filesOpened):
        """Dialog box to select genomic displayed_dataset
        A dialog box will be opened and content of hdf5 file will be read.
        User will have option to select chromosome, resolution and data.
        If user select dataset, define self.shownDataset and self.plotLocation.
        It also select data Aarray using setDataArray().

        Overall, calling this method, following variables will be available for use:

            * self.shownDataset
            * self.plotLocation
            * self.dataArray
            * self.ylimit
            * self.yScaleSteps

        """
        if filename not in filesOpened:
            self.hdf5Hand = gmlib.genomicsDataHandler.HDF5Handler(filename)
            filesOpened[filename] = self.hdf5Hand
        else:
            self.hdf5Hand = filesOpened[filename]

        dialog = browserHelpers.DialogGenomicsDataSelector(self.hdf5Hand, requestedBinsize=self.hiCmapAxis.ccmap.binsize)
        dialog.setChromosomeResolution(self.hiCmapAxis.ccmap.xlabel, gmlib.util.binsizeToResolution(self.hiCmapAxis.ccmap.binsize))
        dialog.exec_()

        if dialog.selected_data is not None:
            self.shownDataset = dialog.selected_data
            self.plotLocation = dialog.whereToPlot
            self.setDataArray()

        dialog.close()

    def readDataFromTextFile(self, filename):
        """Read data from a text file
        """
        dialog = browserHelpers.DialogTextFileSelector(inputFileName=filename)
        dialog.exec_()

        if dialog.inputFileName is not None:
            self.txtFileHand = gmlib.genomicsDataHandler.TextFileHandler(dialog.inputFileName, self.hiCmapAxis.ccmap.shape[0],
                                                     binsize=self.hiCmapAxis.ccmap.binsize, title=dialog.title)
            self.txtFileHand.readData()
            self.plotLocation = dialog.plotPosition
            self.setDataArray()

            dialog.close()

    def setDataArray(self):
        """set Data Array  and Y-scale steps using self.shownData

        It sets self.dataArray and self.yScaleSteps using self.shownData
        """
        if self.shownDataset is not None:
            chrom = self.shownDataset[0]
            res = self.shownDataset[1]
            rdata = self.shownDataset[2]
            self.dataArray = self.hdf5Hand.hdf5[chrom][res][rdata][:]

        if self.txtFileHand is not None:
            self.dataArray = self.txtFileHand.data

        # Get minimum and maximum. Maximum is taken at 95 percentile in case minimum is not negative
        minimum = np.amin(self.dataArray)
        if minimum < 0:
            maximum = np.amax(self.dataArray)
        else:
            maximum = np.percentile(self.dataArray[self.dataArray[:] != 0 ], 95)
        self._ylimit = [minimum, maximum]

        self.yScaleSteps = None

    def changeResolution(self, resolution, change=True):
        """ Change resolution of data.
        In case of dynamic change in resolution, this function try to change the resolution of data.

        Note: It does not change the original data. However, it changes self.dataArray.

        WARNING: This function is messy as there are several return statement.

        Parameters
        ----------
        change : bool
            If it is True, resolution will be changed. If it is False, function will just check that change to new resolution
            is whether possible. A way to know that change in resolution is possible before try to change the resolution and
            mess up everything afterward.

        Returns
        -------
        True (if successful) or False (not successful)

        """

        newBinsize = gmlib.util.resolutionToBinsize(resolution)
        minimumBinsize = None

        # When data is from a hdf5 file
        if self.hdf5Hand is not None:
            chrom = self.shownDataset[0]
            rdata = self.shownDataset[2]
            foundData = False

            # Search if new resolution is already inside the file, no need to downsample
            if self.hdf5Hand.hasDataName(chrom, resolution, rdata):
                # Change the data array in case change is required
                if change:
                    self.shownDataset = (chrom, resolution, rdata)                # Change shownDataset
                    self.dataArray = self.hdf5Hand.hdf5[chrom][resolution][rdata][:]
                    foundData = True
                return True

            # In case if not found, minimum resolution present for data is obtained for latter use
            if not foundData:
                minimumBinsize = np.amin(list(map(gmlib.util.resolutionToBinsize, self.hdf5Hand.getResolutionList(chrom, dataName=rdata))))

        # When data is from a text file
        if self.txtFileHand is not None:
            minimumBinsize = self.txtFileHand.binsize

        # In this case downsampling is not possible. therefore no success
        if newBinsize < minimumBinsize:
            return False

        # Try to determine the factor by which data will be downsampled
        # Note: factor should alyways be an integer, not a float.
        level = None
        if newBinsize > minimumBinsize:
            if newBinsize % minimumBinsize == 0:
                level = newBinsize / minimumBinsize

        # If downsampling factor is not a whole number, in that case abort the downsampling.
        if level is None:
            return False

        # Perform downsampling in case of known level
        success = False

        # When data is from a hdf5 file
        if self.hdf5Hand is not None:
            chrom = self.shownDataset[0]
            minResolution = gmlib.util.binsizeToResolution(minimumBinsize)
            rdata = self.shownDataset[2]

            success = True
            if change:
                self.dataArray = gmlib.ccmap.downSample1D(self.hdf5Hand.hdf5[chrom][minResolution][rdata][:], level=level)
                self.shownDataset = (chrom, resolution, rdata)                # Change shownDataset


        # When data is from a text file
        if self.txtFileHand is not None:
            if change:
                self.dataArray = gmlib.ccmap.downSample1D(self.txtFileHand.data, level=level)
            success = True

        if not success:
            return False
        else:
            return True

    def genChromNameList(self):
        """ Get the list of all possible datanames
        """
        chromNameList = None
        if self.hdf5Hand is not None:
            chromNameList = self.hdf5Hand.getChromList()
            chromNameList = gmlib.util.sorted_nicely(chromNameList)

        return chromNameList

    def changeDataByName(self, browser, name, change=False):
        success = False

        chrom = name
        res = self.shownDataset[1]
        rdata = self.shownDataset[2]

        if self.hdf5Hand is not None:
            if self.hdf5Hand.hasDataName(chrom, res, rdata):
                success = True
                if change:
                    self.shownDataset = (chrom, res, rdata)
                    self.setDataArray()
                    self.xlabel = self.hiCmapAxis.axes_props.xLabel['Text']

        # In case if other formats are loaded
        if self.converterDialog is not None and not success:
            self.hdf5Hand.close()
            self.converterDialog.connectToBrowseMapChanger(browser, chrom)
            self.converterDialog.setChromName(chrom)
            self.converterDialog.show()
            self.converterDialog.startProcess()

        return success

    def set_tree_widget_item(self, treeWidget):
        """Initialize tree widget item
        """
        title = None
        if self.hdf5Hand is not None:
            title = self.hdf5Hand.title

        if self.txtFileHand is not None:
            title = self.txtFileHand.title

        if title is not None:
            label = QLabel(title)
            label.setWordWrap(True)
            self.treeWidgetItem = QTreeWidgetItem()
            self.hiCmapAxis.treeWidgetItem.addChild(self.treeWidgetItem)
            treeWidget.setItemWidget(self.treeWidgetItem, 0, label)

    def updatePlot(self, ax=None):
        """Update plots

        It updates the plots. If a new axes is given, it uses new axes to plot
        Also updates axis properties.

        """
        if ax is not None:
            if self.ax is not None:
                # remove previous axis. Have to use because of KeyError. Does not impact anything
                try:
                    ax.figure.delaxes(self.ax)
                except KeyError:
                    pass
                del self.ax

            self.ax = ax
            if self.axes_props is not None:
                self.axes_props.axes = ax
                self.axes_props.set_to_axes()

        # to get start and end point of plot
        sidx = self.hiCmapAxis.xrange[0]
        eidx = self.hiCmapAxis.xrange[1]

        # If end point is larger than the data, do not plot
        if eidx > self.dataArray.shape[0]:  return

        # If plot is already there, clear it
        if self.plot is not None:
            self.ax.clear()
            del self.plot

        # If x-tick-locations are not set in parent ccmap axes, set here
        if self.hiCmapAxis.axes_props.xTickLocations is None:
            self.hiCmapAxis.update_axes_props()

        # plot
        #self.plot = self.ax.plot(self.hiCmapAxis.axes_props.xTickLocations, self.dataArray[sidx : eidx])

        # Vertical lines plot
        xticks, minvalue, maxvalue = self.get_vlines_min_max(self.hiCmapAxis.axes_props.xTickLocations, self.dataArray[sidx : eidx])
        self.plot = self.ax.vlines(xticks, minvalue, maxvalue, colors=self.plotColor, linewidths=self.plotLineWidth)

        # Set upper and lower limit along y-axis
        self.ax.set_ylim(self.ylimit[0], self.ylimit[1])

        # set y tick label formatting
        self.ax.ticklabel_format(style=self.yticksFormat, axis='y')

        # Set upper and lower limit along x-axis using parent ccmap axes
        self.ax.set_xlim(self.hiCmapAxis.axes_props.xTickLocations[0], self.hiCmapAxis.axes_props.xTickLocations[-1])

        # Update axis properties
        self.update_axes_props()

        # Update parent ccmap axes properties
        self.hiCmapAxis.update_axes_props()

    def get_vlines_min_max(self, xticks, data):
        """To get minimum and maximum values for each vertical line
        """
        #Index for all non-zeros
        idx = np.nonzero( data != 0 )

        # get signs of each entry and their index
        sign = np.sign(data)
        neg_idx = np.nonzero( sign == -1 )
        pos_idx = np.nonzero( sign == 1 )

        minvalue = np.zeros(data.shape)
        maxvalue = np.zeros(data.shape)

        minvalue[neg_idx] = data[neg_idx]
        maxvalue[pos_idx] = data[pos_idx]

        return xticks[idx], minvalue[idx], maxvalue[idx]

    def update_axes_props(self):
        """ To update axis properties.
        """
        # When not present, initialize new
        if self.axes_props is None:
            self.ax.set_xlabel(self.hiCmapAxis.axes_props.xLabel['Text'], fontsize=12)   # Set xlabel for first time
            self.axes_props = browserHelpers.AxesProperties(self.ax)
            self.axes_props.yTickLabel['Tick Intervals'] = 3

        # update y-tick-locations
        self.axes_props.yTickLocations = self.yticks
        self.axes_props.yTickLabelTexts = self.yticks

        # Hide x-label and x-ticklabel if not the lowest plots
        if self.hiCmapAxis.lowerMostGenomicPlotAxes != self.index:
            self.axes_props.xLabel['Show Label'] = 'none'
            self.axes_props.xTickLabel['Label Position'] = 'none'

        self.axes_props.set_to_axes()

class CCMAPAXIS:
    def __init__(self, index, ax):
        self.index = index                      # Current index to parent instance
        self.ax = ax                            # matplotlib axes instance
        self.axes_props = None                  # AxesProperties instance

        self.fileType = None                    # Type of file loaded ccmap or gcmap

        self.image = None                       # ax.imshow instance
        self.ccmap = None                       # CCMAP instance
        self.mapUnit = '10kb'                   # Unit of map in resolution
        self.resolution = None
        self.title = None                       # title of ccmap
        self.xrange = None                      # Current upper and lower limit of x-axis
        self.yrange = None                      # Current upper and lower limit of y-axis

        self.treeWidgetItem = None              # Tree widget item instance
        self.qActionToAddGenomicDataset = None  # action menu instance to add new genomic data plot as child

        # GUI Options
        self.colorRangeMaxValue = None          # Maximum value for color-map
        self.colorRangeMinValue = None          # Minimum value for color-map
        self.color_scale_steps = None           # Array of 100 steps between maximum and minimum
        self.color_scale_slider_value = None    # current value on slider
        self.marker = None

        #  Global color scaler
        self.scaleMinMaxForResolutionInterchange = False

        self.interchangeableCMapNames = None    # List of all contact maps available from file

        self.xticks = None                      # xtick locations -- all values between maximum and minimum
        self.yticks = None                      # yticks locations -- all values between maximum and minimum
        self.xticklabels = None                 # x-tick-labels for all xtick locations
        self.yticklabels = None                 # y-tick-labels for all ytick locations

        # Sub Genomic plots
        self.genomicPlotAxes = None             # List of child GenomicDataPlotAxis instances
        self.divider = None
        self.upperMostGenomicPlotAxes = None
        self.lowerMostGenomicPlotAxes = None

        # Maximum and minimum value of colormap
        self.vmin = None
        self.vmax = None

        # All GUI option are connected to plotting method. When a plot become active, it changes also GUI options, which in turn
        # plot again. To prevent this use this varaible to block plotting.
        self.doNotPlot = False

        # Properties that update plots. These are private variables. However, updating the public counterpart also updates plots
        self._rangeXY = None
        self._colorScaleStatus = 'Change color linearly'
        self._color_scale_spinbox_value = None
        self._colormap = 'gray_r'
        self._interpolation = 'none'
        self._verticalSpace = 0.2
        self._xlabel = None                      # label for x-axis
        self._ylabel = None                      # label for y-axis


    @property
    def rangeXY(self):
        return self._rangeXY

    @rangeXY.setter
    def rangeXY(self, value):
        """ Update xrange, yrange and plot
        """
        if self.doNotPlot:  return

        xrange, yrange = value

        replot = False
        if self.xrange is not None and self.yrange is not None:
            CurrXdiff = xrange[1] - xrange[0]
            PrevXdiff = self.xrange[1] - self.xrange[0]
            CurrYdiff = yrange[1] - yrange[0]
            PrevYdiff = self.yrange[1] - self.yrange[0]

            # This one for zoom in and out
            if CurrXdiff != PrevXdiff or CurrYdiff != PrevYdiff:
                replot = True
            else:
                replot = False

        self._rangeXY = (xrange, yrange)
        self.xrange = xrange
        self.yrange = yrange
        self.makeMapImage(replot=replot)

        if self.genomicPlotAxes is not None:
            for obj in self.genomicPlotAxes:
                obj.updatePlot()

    @property
    def color_scale_spinbox_value(self):
        return self._color_scale_spinbox_value

    @color_scale_spinbox_value.setter
    def color_scale_spinbox_value(self, value):
        """Update color scale spinbox value and plot
        """
        if self.doNotPlot:  return

        self._color_scale_spinbox_value = value
        self.makeMapImage(replot=True)

    @property
    def colorScaleStatus(self):
        return self._colorScaleStatus

    @colorScaleStatus.setter
    def colorScaleStatus(self, newColorScaleStatus):
        """Update color scale type and plot
        """
        if self.doNotPlot:  return

        if self.ccmap.minvalue < 0.0:
            self._colorScaleStatus == 'Change color linearly'
        else:
            self._colorScaleStatus = newColorScaleStatus

        minvalue = None
        maxvalue = None
        if self._colorScaleStatus == 'Logarithm of map':
            minvalue = np.log(self.colorRangeMinValue)
            maxvalue = np.log(self.colorRangeMaxValue)
            self.color_scale_steps = np.linspace(minvalue , maxvalue, 101, endpoint=True)
            self.color_scale_steps = np.linspace(self.color_scale_steps[1] , maxvalue, 101, endpoint=True)

        if self._colorScaleStatus == 'Change color logarithmically':
            minvalue = self.colorRangeMinValue
            maxvalue = self.colorRangeMaxValue
            self.color_scale_steps = np.linspace( np.log(minvalue), np.log(maxvalue), 101, endpoint=True)
            self.color_scale_steps = np.exp(self.color_scale_steps)

        if self._colorScaleStatus == 'Change color linearly':
            minvalue = self.colorRangeMinValue
            maxvalue = self.colorRangeMaxValue
            self.color_scale_steps = np.linspace( minvalue, maxvalue, 101, endpoint=True)

        self._color_scale_spinbox_value = maxvalue
        self.makeMapImage(replot=True)

    @property
    def colormap(self):
        return self._colormap

    @colormap.setter
    def colormap(self, value):
        """Update color map type and plot
        """
        if self.doNotPlot:  return

        if self._colormap != value:
            self._colormap = value
            self.makeMapImage(replot=True)

    @property
    def interpolation(self):
        return self._interpolation

    @interpolation.setter
    def interpolation(self, value):
        """Update interpolation method and plot
        """
        if self.doNotPlot:  return

        if self._interpolation != value:
            self._interpolation = value
            self.makeMapImage(replot=True)

    @property
    def verticalSpace(self):
        return self._verticalSpace

    @verticalSpace.setter
    def verticalSpace(self, value):
        """Update vertical space between plots
        """
        if value == self._verticalSpace:    return

        self._verticalSpace = value
        if self.genomicPlotAxes is not None:
            self.updateGenomicPlotAxes(resize=True)

    @property
    def xlabel(self):
        return self._xlabel

    @xlabel.setter
    def xlabel(self, value):
        """ Also set xlabel in plot
        """
        self._xlabel = value
        if self.axes_props is not None:
            self.axes_props.xLabel['Text'] = value
            self.axes_props.set_to_axes()
        else:
            self.ax.set_xlabel(value, fontsize=12)

    @property
    def ylabel(self):
        return self._ylabel

    @ylabel.setter
    def ylabel(self, value):
        """ ALso set ylabel to plots
        """
        self._ylabel = value
        if self.axes_props is not None:
            self.axes_props.yLabel['Text'] = value
            self.axes_props.set_to_axes()
        else:
            self.ax.set_ylabel(value, fontsize=12)

    def initialize(self, path, fileType, mapName=None, resolution=None, filesOpened=None):
        """Initialize the class after first call.
        Read maps and set several options/variables in accordance with maps
        """

        self.set_ccmap(path, fileType, mapName=mapName, resolution=resolution, filesOpened=filesOpened)

        self.colorRangeMinValue = self.ccmap.minvalue
        self.colorRangeMaxValue = self.ccmap.maxvalue

        self.xticklabels, self.yticklabels = self.ccmap.get_ticks()
        self.xticklabels = self.xticklabels/gmlib.util.resolutionToBinsize(self.mapUnit)
        self.yticklabels = self.yticklabels/gmlib.util.resolutionToBinsize(self.mapUnit)

        self.resolution = gmlib.util.binsizeToResolution( self.ccmap.binsize )

        if self.ccmap.xlabel is None:
            self.xlabel = 'Unknown' + ' [{0}]'.format(self.mapUnit)
        else:
            self.xlabel = self.ccmap.xlabel + ' [{0}]'.format(self.mapUnit)

        if self.ccmap.ylabel is None:
            self.ylabel = 'Unknown' + ' [{0}]'.format(self.mapUnit)
        else:
            self.ylabel = self.ccmap.ylabel + ' [{0}]'.format(self.mapUnit)

        self.color_scale_steps = np.linspace( self.colorRangeMinValue, self.colorRangeMaxValue, 101, endpoint=True)

    def updatePropsForResolution(self):
        """ Update ticklabels and resolution according to the new resolution
        """
        self.resolution = gmlib.util.binsizeToResolution( self.ccmap.binsize )

        self.xticklabels, self.yticklabels = self.ccmap.get_ticks()
        self.xticklabels = self.xticklabels/gmlib.util.resolutionToBinsize(self.mapUnit)
        self.yticklabels = self.yticklabels/gmlib.util.resolutionToBinsize(self.mapUnit)

    def set_color_spinbox_slider_values(self, spinbox, slider):
        """Used to directly assign the values for spinbox
        """
        self._color_scale_spinbox_value = spinbox
        self.color_scale_slider_value = slider

    def update_axes(self, ax):
        self.ax = None
        self.image = None
        self.ax = ax
        self.axes_props.axes = None
        self.axes_props.axes = ax
        self.axes_props.set_to_axes()
        self.makeMapImage(replot=True)
        self.updateGenomicPlotAxes(resize=True)

    def updateGenomicPlotAxes(self, resize=False):
        if self.genomicPlotAxes is None:    return

        if resize:
            self.divider = None

        if self.divider is None:
            self.divider = make_axes_locatable(self.ax)

        if resize:
            for obj in self.genomicPlotAxes:
                ax = self.divider.append_axes(obj.plotLocation, size=0.6, pad=self.verticalSpace, sharex=self.ax)
                obj.updatePlot(ax)
        else:
            ax = self.divider.append_axes(self.genomicPlotAxes[-1].plotLocation, size=0.6, pad=self.verticalSpace, sharex=self.ax)
            self.genomicPlotAxes[-1].updatePlot(ax)

        for obj in self.genomicPlotAxes:
            if self.lowerMostGenomicPlotAxes  != obj.index:
                obj.axes_props.xLabel['Show Label'] = 'none'
                obj.axes_props.xTickLabel['Label Position'] = 'none'
                obj.axes_props.set_to_axes()

    def remove_ccmap(self):
        ''' To initialize a new ccmap
        '''
        if self.ccmap is not None:
            del self.ccmap

        if self.image is not None:
            del self.image

        self.ax.clear()

        self.ccmap = None
        self.image = None
        self.xrange = None
        self.yrange = None
        self.xticklabels = None
        self.yticklabels = None

    def set_ccmap(self, path, fileType, mapName=None, resolution=None, filesOpened=None):
        if fileType == 'ccmap':
            self.ccmap = gmlib.ccmap.load_ccmap(path)
            self.fileType = 'ccmap'

            # In case if minimum value is zero, change it
            if self.ccmap.minvalue == 0:
                self.ccmap.make_readable()
                ma = np.ma.masked_equal(self.ccmap.matrix, 0.0, copy=False)
                self.ccmap.minvalue = ma.min()
                del ma
                self.ccmap.make_unreadable()
        else:
            if path not in filesOpened:
                filesOpened[path] = h5py.File(path, 'r')

            self.ccmap = gmlib.gcmap.GCMAP(filesOpened[path], mapName=mapName, resolution=resolution)
            self.fileType = 'gcmap'

            # In case if minimum value is zero, change it
            if self.ccmap.minvalue == 0:
                ma = np.ma.masked_equal(self.ccmap.matrix, 0.0, copy=False)
                self.ccmap.minvalue = ma.min()
                del ma

            # Make list of all contact maps
            self.tryEnableInterchangeCMapName()

    def set_tree_widget_item(self, treeWidget):
        label = QLabel(self.title)
        label.setWordWrap(True)
        self.treeWidgetItem = QTreeWidgetItem()
        treeWidget.addTopLevelItem(self.treeWidgetItem)
        treeWidget.setItemWidget(self.treeWidgetItem, 0, label)
        #self.treeWidgetItem =  QTreeWidgetItem([self.title], 0)

    def update_axes_props(self):

        if self.axes_props is None:
            xdiff = self.xrange[1] - self.xrange[0]
            ydiff = self.yrange[1] - self.yrange[0]
            xpos_tick = np.linspace(0, xdiff-1, 9, dtype=int)
            ypos_tick = np.linspace(0, ydiff-1, 9, dtype=int)
            xticks = np.linspace(self.xrange[0], self.xrange[1]-1, 9, dtype=int)
            yticks = np.linspace(self.yrange[0], self.yrange[1]-1, 9, dtype=int)
            self.ax.set_xticks(xpos_tick)
            self.ax.set_xticklabels(self.xticklabels[xticks], fontsize=12)
            self.ax.set_xlabel(self.xlabel, fontsize=12)
            self.ax.set_yticks(ypos_tick)
            self.ax.set_yticklabels(self.yticklabels[yticks], fontsize=12)
            self.ax.set_ylabel(self.ylabel, fontsize=12)

            self.axes_props = browserHelpers.AxesProperties(self.ax)

        else:
            xdiff = self.xrange[1] - self.xrange[0]
            ydiff = self.yrange[1] - self.yrange[0]
            xpos_tick = np.arange(0, xdiff)
            ypos_tick = np.arange(0, ydiff)

            xticks = np.arange(self.xrange[0], self.xrange[1], dtype=int)
            yticks = np.arange(self.yrange[0], self.yrange[1], dtype=int)
            self.axes_props.xTickLocations = xpos_tick
            self.axes_props.yTickLocations = ypos_tick
            self.axes_props.xTickLabelTexts = self.xticklabels[xticks]
            self.axes_props.yTickLabelTexts = self.yticklabels[yticks]

            # Hide x-label and x-ticklabel if other plot is present
            if self.lowerMostGenomicPlotAxes is not None:
                self.axes_props.xLabel['Show Label'] = 'none'
                self.axes_props.xTickLabel['Label Position'] = 'none'

            self.axes_props.set_to_axes()

    def set_vmin_vmax(self):
        minvalue = self.color_scale_steps[0]
        if self.color_scale_spinbox_value is None:
            maxvalue = self.color_scale_steps[-1]
        else:
            maxvalue = self.color_scale_spinbox_value

        if minvalue==0 and self.ccmap.minvalue != 0:
            minvalue =  self.colorRangeMinValue

        self.vmin = minvalue
        self.vmax = maxvalue

    def makeMapImage(self, replot=False):
        ''' Make image of map for the given range
        '''

        if self.doNotPlot:  return

        matrix = self.ccmap.matrix[self.yrange[0]:self.yrange[1], self.xrange[0]:self.xrange[1]]

        # Scale matrix
        if self.scaleMinMaxForResolutionInterchange:
            matrix = browserHelpers.scaleMatrix(matrix, self.colorRangeMinValue, self.colorRangeMaxValue)

        # Change log by masking zero and filled it with minimum value
        if self.colorScaleStatus == 'Logarithm of map':
            minvalue = self.color_scale_steps[0] - (self.color_scale_steps[1] - self.color_scale_steps[0])
            new_matrix = ma.log(matrix)
            matrix = new_matrix.filled(minvalue)

        if self.image is None:
            self.ax.clear()
            self.set_vmin_vmax()
            self.image = self.ax.imshow( matrix, vmin=self.vmin, vmax=self.vmax, origin='lower', aspect='equal', cmap=self.colormap, interpolation=self.interpolation, resample=True )

        elif replot:
            self.ax.clear()
            self.set_vmin_vmax()
            if self.image is not None:
                self.image = None

            self.image = self.ax.imshow( matrix, vmin=self.vmin, vmax=self.vmax, origin='lower', aspect='equal', cmap=self.colormap, interpolation=self.interpolation, resample=True )

        else:
            self.image.set_data( matrix )

        self.update_axes_props()

    def tryEnableInterchangeCMapName(self):
        """ Try to enable interchange in contact map name.

        If it is successful, self.interchangeableCMapNames contains list of chromosome
        """
        GCMap = True
        mapNameMatched = True
        chroms = []

        if self.fileType == 'ccmap':
            GCMap = False
        else:
            self.ccmap.genMapNameList()          # Generate only once, so calling several time does not do anything
            chroms.append( self.ccmap.mapNameList )

        # Append list from genomic dataset also
        if self.genomicPlotAxes is not None:
            for gax in self.genomicPlotAxes:

                # In case of converter dialog, do not generate common map
                if gax.converterDialog is not None:
                    continue

                dataNames = gax.genChromNameList()

                # In case only data names are possible
                if dataNames is None:
                    mapNameMatched = False
                else:
                    chroms.append( dataNames )

        # Now determine common map names between all datasets
        if GCMap and mapNameMatched:
            common = None

            # If only one map is loaded
            if len(chroms) == 1:
                common = set( chroms[0] )

            # If more than one, determine common available resolutions among all maps
            for i in range(len(chroms) - 1 ):
                b1 = set( chroms[i] )
                b2 = set( chroms[i+1] )
                c = set.intersection( b1 , b2 )
                if not c:
                    mapNameMatched = False
                else:
                    if common is None:
                        common = c
                    else:
                        common = set.intersection( common, c )
                        if common is None:
                            mapNameMatched = False

            if mapNameMatched:
                self.interchangeableCMapNames = gmlib.util.sorted_nicely( common )
            else:
                self.interchangeableCMapNames = None

        else:
            self.interchangeableCMapNames = None


if __name__=="__main__":
    main()
