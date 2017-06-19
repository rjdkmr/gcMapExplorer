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

import os
import numpy as np
import h5py

import tempfile

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.uic import loadUiType

import matplotlib as mpl
from matplotlib import font_manager as mplFontManager
from matplotlib import colors as mplColors
from matplotlib.ticker import AutoMinorLocator, MaxNLocator

import gcMapExplorer.lib as gmlib

from . import guiHelpers

# Determine absolute path to UIs directory. Relative path from this directory does not work.
DirToThisScript = os.path.dirname(os.path.abspath(__file__))
PathToUIs = os.path.join(DirToThisScript, 'UIs')

# Dialog box to change page size by custom length
pathToThisUI = os.path.join(PathToUIs, 'dialogCustomPageSize.ui')
Ui_DialogCustomPlotsize, DialogCustomPlotsizeBase = loadUiType(pathToThisUI)
class DialogCustomPlotsize(DialogCustomPlotsizeBase, Ui_DialogCustomPlotsize):
    def __init__(self, figsize):
        super(DialogCustomPlotsize, self).__init__()
        self.setupUi(self)

        # Validate input as double
        self.lineEditHeight.setValidator(QDoubleValidator())
        self.lineEditWidth.setValidator(QDoubleValidator())

        # Store figsize for later use
        self.figsize = list(figsize)

        # Display figsize in boxes
        self.lineEditWidth.setText(str(figsize[0]))
        self.lineEditHeight.setText(str(figsize[1]))

        # GEt orientation from figure size
        self.PageOrientation = self.check_orientation()

        # Connect to return button
        self.lineEditWidth.editingFinished.connect(self.change_figsize)
        self.lineEditHeight.editingFinished.connect(self.change_figsize)
        self.buttonArrow.clicked.connect(self.change_orientation)

    def change_figsize(self):
        self.figsize = (float(self.lineEditWidth.text()), float(self.lineEditHeight.text()))
        self.PageOrientation = self.check_orientation()

    def check_orientation(self):
        self.PageOrientation = None
        if self.figsize[1] >= self.figsize[0]:
            self.lineEditOrientation.setText('Portrait')
            self.PageOrientation = 'Portrait'
        else:
            self.lineEditOrientation.setText('Landscape')
            self.PageOrientation = 'Landscape'

    def change_orientation(self):
        self.figsize = self.figsize[::-1]
        self.PageOrientation = self.check_orientation()
        # Display figsize in boxes
        self.lineEditWidth.setText(str(self.figsize[0]))
        self.lineEditHeight.setText(str(self.figsize[1]))


# Dialog box to choose and load genomic data
pathToThisUI = os.path.join(PathToUIs, 'gcmapSelectorForBrowser.ui')
Ui_gcmapSelectorDialog, gcmapSelectorDialogBase = loadUiType(pathToThisUI)
class GCMapSelectorDialog(gcmapSelectorDialogBase, Ui_gcmapSelectorDialog):
    def __init__(self, filename, hdf5=None):
        super(GCMapSelectorDialog, self).__init__()
        self.setupUi(self)

        self.filename = filename    # Input gcmap file
        self.mapName = None         # Result
        self.resolution = None      # Result
        self.dictBinSizes = None    # Dictionary of binsizes available for maps

        # h5py File object
        self.hdf5 = hdf5
        if hdf5 is None:
            self.needToCloseHdf5 = True
        else:
            self.needToCloseHdf5 = False

        self.mapList = None                              # List of maps in hdf5 file
        self.dictMapListWidgetItems = None               # Dictionary for map listWidgetItem for respective map
        self.dictResolutionListWidgetItems = None        # Dictionary for resolution listWidgetItem for selected map listWidgetItem

        self.fileNameLineEdit.setText( self.filename )   # Set name of filename
        self.getMapList()                                # Generate list of maps for first time
        self.cancelButton.setFocus()                     # Set focus on cancel button as no map is selected yet

        # Set icon for browse file
        self.browseFileButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )

        # Connnecting widget to function
        self.fileNameLineEdit.editingFinished.connect( self.checkInputFile )

        self.mapListWidget.itemClicked.connect( self.displayResolutionList )
        self.resolutionListWidget.itemClicked.connect( self.displayAttributes )
        self.resolutionListWidget.itemClicked.connect( self.okButton.setFocus )

        self.browseFileButton.clicked.connect( self.browseFile )
        self.okButton.clicked.connect( lambda: self.onButtonClicked(self.okButton) )
        self.cancelButton.clicked.connect( lambda: self.onButtonClicked(self.cancelButton) )

        # Disable editing of table's cell
        self.mapAttrsTableWidget.setEditTriggers( QAbstractItemView.NoEditTriggers )

    def getMapList(self):
        """Read gcmap file, make a list of maps and display it in QListWidget
        """
        if self.hdf5 is None:
            self.hdf5 = h5py.File(self.filename)

        self.mapList = gmlib.util.sorted_nicely( list( self.hdf5.keys() ) )

        self.dictMapListWidgetItems = dict()
        for key in self.mapList:
            listItem = QListWidgetItem(key, self.mapListWidget)
            self.dictMapListWidgetItems[key] = listItem

    def displayResolutionList(self, item):
        """ Generate a resolution list for selected map
        """

        mapName = self.getMapForListWidget(item)

		# determining finest resolution map
        if self.dictBinSizes is None:
            self.dictBinSizes = dict()

        if mapName not in self.dictBinSizes:
            binsizes = []
            for key in self.hdf5[mapName]:
                if 'bNoData' not in key:
                    binsizes.append( self.hdf5[mapName][key].attrs['binsize'] )
            binsizes = sorted(binsizes)
            self.dictBinSizes[mapName] = binsizes

        if self.dictResolutionListWidgetItems is not None:
            for i in range(self.resolutionListWidget.count()):
                item = self.resolutionListWidget.takeItem( 0 )
                del item

        self.dictResolutionListWidgetItems = None
        self.dictResolutionListWidgetItems = dict()
        for binsize in self.dictBinSizes[mapName]:
            listItem = QListWidgetItem(gmlib.util.binsizeToResolution(binsize), self.resolutionListWidget)
            self.dictResolutionListWidgetItems[binsize] = listItem

        resolutionItem = self.dictResolutionListWidgetItems[ self.dictBinSizes[mapName][0] ]
        self.resolutionListWidget.setCurrentItem( resolutionItem )
        self.displayAttributes( resolutionItem )


    def displayAttributes(self, item):
        """When a Map is selected, display its attributes in table (QTableWidget)
        """
        mapName = self.getMapForListWidget( self.mapListWidget.currentItem() )
        resolution = self.getResolutionForListWidget( item )

        self.mapAttrsTableWidget.item(0, 0).setText( str(self.hdf5[mapName].attrs['xlabel']) )
        self.mapAttrsTableWidget.item(1, 0).setText( str(self.hdf5[mapName].attrs['ylabel']) )

        attrs = ['minvalue', 'maxvalue', 'binsize', 'xshape', 'yshape']
        for i in range(len(attrs)):
            self.mapAttrsTableWidget.item(i+2, 0).setText( str(self.hdf5[mapName][resolution].attrs[attrs[i]]) )

    def getMapForListWidget(self, item):
        """Get map name for given QListWidgetItem
        """
        mapName = None
        for key in self.mapList:
            if item is self.dictMapListWidgetItems[key]:
                mapName = key
                break
        return mapName

    def getResolutionForListWidget(self, item):
        """Get map name for given QListWidgetItem
        """

        mapName = self.getMapForListWidget( self.mapListWidget.currentItem() )
        resolution = None

        for binsize in self.dictBinSizes[mapName]:
            if item is self.dictResolutionListWidgetItems[binsize]:
                resolution = gmlib.util.binsizeToResolution(binsize)
                break
        return resolution

    def resetForNewFile(self):
        """Reset all widgets and class attributes.
        Call it when a new valid file is selected.
        """
        for i in range(self.mapListWidget.count()):
            item = self.mapListWidget.takeItem( 0 )
            del item

        self.mapList = None
        self.dictMapListWidgetItems = None

        for i in range(7):
            self.mapAttrsTableWidget.item(i, 0).setText('Select a Map')

        if self.needToCloseHdf5:
            self.hdf5.close()
        self.hdf5 = None
        self.cancelButton.setFocus()

    def checkInputFile(self):
        """Check input file exist. If file exist, read it, generate a list of map, and show it.
        """
        inputFile = str( self.fileNameLineEdit.text() )
        if inputFile and not os.path.isfile(inputFile):
            msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'File not found.', QMessageBox.Ok, self)
            msgBox.exec_()
            msgBox.close()

            self.fileNameLineEdit.clear()
            self.fileNameLineEdit.setText(self.filename)
        else:
            self.filename = inputFile
            self.resetForNewFile()
            self.getMapList()

    def browseFile(self):
        """Browse and select input text file. If file is selected, read it, generate a list of map, and show it.
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = "gcmap file (*.gcmap);;All files (*.*)"
        path = QFileDialog.getOpenFileName(self, 'Select a gcmap file', '', file_choices)
        if path[0]:
            self.fileNameLineEdit.setText( path[0] )
            self.resetForNewFile()
            self.getMapList()

    def onButtonClicked(self, button):
        """when ok/cancel button is clicked
        """
        if button is self.okButton:
            item = self.mapListWidget.currentItem()
            if item:
                self.mapName = self.getMapForListWidget(item)

                resItem = self.resolutionListWidget.currentItem()
                if resItem:
                    self.resolution = self.getResolutionForListWidget( resItem )
                else:
                    self.reject()

                self.accept()
            else:
                self.reject()
        else:
            self.reject()


# Dialog box to change axis properties
pathToThisUI = os.path.join(PathToUIs, 'dialogAxisProps.ui')
Ui_DialogAxisProps, QDialogAxisPropsBase = loadUiType(pathToThisUI)
class DialogAxisProps(QDialogAxisPropsBase, Ui_DialogAxisProps):
    def __init__(self, canvas, tabIdx, axes, axesProps=None, xTickLocations=None, xTickLabelTexts=None, yTickLocations=None, yTickLabelTexts=None):
        super(DialogAxisProps, self).__init__()
        self.setupUi(self)

        self.tabWidgetAxisProps.setCurrentIndex(tabIdx)
        self.axes = axes
        self.canvas = canvas

        self.init_font_matplotlib_qt()

        if axesProps is None:
            self.axesProps = AxesProperties(axes, xTickLocations, xTickLabelTexts, yTickLocations, yTickLabelTexts)
        else:
            self.axesProps = axesProps

        # Display X label properties
        self.displayXLabelPropsOnTabWidget()

        # Display Y label properties
        self.displayYLabelPropsOnTabWidget()

        # Display X tick label properties
        self.displayXTickLabelPropsOnTabWidget()

        # Display X tick label properties
        self.displayYTickLabelPropsOnTabWidget()

        # What happens when buttons are clicked
        self.xTickShowMinorTicksCBox.currentIndexChanged.connect(self.xMinorTicksController)
        self.yTickShowMinorTicksCBox.currentIndexChanged.connect(self.yMinorTicksController)
        self.okCancelApplyButtonBox.clicked.connect(self.uponButtonClicked)

    def update_axes(self, axes, axesProps, xTickLocations=None, xTickLabelTexts=None, yTickLocations=None, yTickLabelTexts=None):
        ''' If dialog box is already opened by user, use this function to update dialog box for given axes and axes properties
        '''
        self.axes = axes
        self.axesProps = axesProps

        # Display X label properties
        self.displayXLabelPropsOnTabWidget()

        # Display Y label properties
        self.displayYLabelPropsOnTabWidget()

        # Display X tick label properties
        self.displayXTickLabelPropsOnTabWidget()

        # Display X tick label properties
        self.displayYTickLabelPropsOnTabWidget()

    def init_font_matplotlib_qt(self):
        ''' Initialize font name list QFontComboBox.
        Also adds font from matplotlib default fonts
        '''
        # Font database of qt
        self.qtFontDatabase = QFontDatabase()

        # Getting path of matplotlib ttf directory
        mplDataPath = mpl.rcParams['datapath']
        mplFontPath = os.path.join(mplDataPath, 'fonts')
        mplTtfFontPath = os.path.join(mplFontPath, 'ttf')

        # Buinding ttf files list
        mplTtfList = []
        for f in os.listdir(mplTtfFontPath):
            if f.endswith('.ttf'):
                mplTtfList.append(os.path.join(mplTtfFontPath, f))

        # Getting name of each font from ttf files
        mplFontNameList = [mplFontManager.FontProperties(fname=fname).get_name() for fname in mplTtfList]

        # Getting name of each font from qt font database
        qtFontNameListQtype = self.qtFontDatabase.families()
        qtFontNameList =  [ str(qtFont) for qtFont in qtFontNameListQtype ]

        # Comapring and adding matplotlib font to qt font database
        for i in range(len(mplFontNameList)):
            matched = False
            for qtFont in qtFontNameList:
                if str(qtFont) == mplFontNameList[i]:
                    matched = True

            if not matched:
                self.qtFontDatabase.addApplicationFont(mplTtfList[i])

    def is_font_exist_in_mpl(self, fontname):
        result = None
        result = mplFontManager.findfont(mplFontManager.FontProperties(family=fontname))

        if mplFontManager.FontProperties(fname=result).get_name() == fontname:
            return True
        else:
            return False

    def displayXLabelPropsOnTabWidget(self):
        ''' Display X Label properties on Tab Widget
        '''
        # Setting X label in tab widget
        self.xLabelTitleLineEdit.setText(self.axesProps.xLabel['Text'])

        # Setting x label postion in Show Label combo box
        if self.axesProps.xLabel['Show Label'] == 'bottom':
            self.xLabelShowLabelCBox.setCurrentIndex(0)
        if self.axesProps.xLabel['Show Label'] == 'top':
            self.xLabelShowLabelCBox.setCurrentIndex(1)
        if self.axesProps.xLabel['Show Label'] == 'none':
            self.xLabelShowLabelCBox.setCurrentIndex(2)

        # Setting font in qfontcombobox
        qfont = self.qtFontDatabase.font(self.axesProps.xLabel['Font Name'], 'normal', self.axesProps.xLabel['Font Size'])
        self.xLabelFontNameCBox.setCurrentFont(qfont)

        # setting font size in tab widget
        self.xLabelFontSizeCBox.setValidator(QIntValidator())
        self.xLabelFontSizeCBox.setCurrentText(str(self.axesProps.xLabel['Font Size']))

        # padding
        self.xLabelPaddingSpinbox.setValue(self.axesProps.xLabel['padding'])

    def changeXLabelPropsOnAxes(self):
        ''' Changes X Label properties on axes plot
        '''
        # Check whether x-label is changed in tab-widget and then change it in plot
        if self.axesProps.xLabel['Text'] != self.xLabelTitleLineEdit.text():
            self.axes.get_xaxis().set_label_text(str(self.xLabelTitleLineEdit.text()))

        # Check whether font size is changed and then change it in plot
        if self.axesProps.xLabel['Font Size'] != int( self.xLabelFontSizeCBox.currentText() ):
            self.axes.get_xaxis().get_label().set_fontsize( int( self.xLabelFontSizeCBox.currentText() ) )

        # Check whether font name is changed in tab-widget and then chenge it in plot
        qfont = self.xLabelFontNameCBox.currentFont()
        if self.is_font_exist_in_mpl(str(qfont.family())):
            if self.axesProps.xLabel['Font Name'] != str(qfont.family()):
                self.axes.get_xaxis().get_label().set_fontname(str(qfont.family()))
        else:
            # resetting font in qfontcombobox if font is not present in matplotlib
            qfont = self.qtFontDatabase.font(self.axesProps.xLabel['Font Name'], 'normal', self.axesProps.xLabel['Font Size'])
            self.xLabelFontNameCBox.setCurrentFont(qfont)
            self.axes.get_xaxis().get_label().set_fontname(str(qfont.family()))

        # Check whether x-lable position is changed in widget and then change it
        if self.xLabelShowLabelCBox.currentIndex() == 0:
            self.axes.get_xaxis().get_label().set_visible(True)
            self.axes.get_xaxis().set_label_position('bottom')
        if self.xLabelShowLabelCBox.currentIndex() == 1:
            self.axes.get_xaxis().get_label().set_visible(True)
            self.axes.get_xaxis().set_label_position('top')
        if self.xLabelShowLabelCBox.currentIndex() == 2:
            self.axes.get_xaxis().get_label().set_visible(False)

        # Padding
        self.axes.xaxis.labelpad = self.xLabelPaddingSpinbox.value()

    def displayYLabelPropsOnTabWidget(self):
        ''' Display Y Label properties on Tab Widget
        '''
        # Setting Y label in tab widget
        self.yLabelTitleLineEdit.setText(self.axesProps.yLabel['Text'])

        # Setting Y label postion in Show Label combo box
        if self.axesProps.yLabel['Show Label'] == 'left':
            self.yLabelShowLabelCBox.setCurrentIndex(0)
        if self.axesProps.yLabel['Show Label'] == 'right':
            self.yLabelShowLabelCBox.setCurrentIndex(1)
        if self.axesProps.yLabel['Show Label'] == 'none':
            self.yLabelShowLabelCBox.setCurrentIndex(2)

        # Setting font in qfontcombobox
        qfont = self.qtFontDatabase.font(self.axesProps.yLabel['Font Name'], 'normal', self.axesProps.yLabel['Font Size'])
        self.yLabelFontNameCBox.setCurrentFont(qfont)

        # setting font size in tab widget
        self.yLabelFontSizeCBox.setValidator(QIntValidator())
        self.yLabelFontSizeCBox.setCurrentText(str(self.axesProps.yLabel['Font Size']))

        # padding
        self.yLabelPaddingSpinbox.setValue(self.axesProps.yLabel['padding'])

    def changeYLabelPropsOnAxes(self):
        ''' Changes Y Label properties on axes plot
        '''
        # Check whether x-label is changed in tab-widget and then change it in plot
        if self.axesProps.yLabel['Text'] != self.yLabelTitleLineEdit.text():
            self.axes.get_yaxis().set_label_text(str(self.yLabelTitleLineEdit.text()))

        # Check whether font size is changed and then change it in plot
        if self.axesProps.yLabel['Font Size'] != int( self.yLabelFontSizeCBox.currentText() ):
            self.axes.get_yaxis().get_label().set_fontsize( int( self.yLabelFontSizeCBox.currentText() ) )

        # Check whether font name is changed in tab-widget and then chenge it in plot
        qfont = self.yLabelFontNameCBox.currentFont()
        if self.is_font_exist_in_mpl(str(qfont.family())):
            if self.axesProps.yLabel['Font Name'] != str(qfont.family()):
                self.axes.get_yaxis().get_label().set_fontname(str(qfont.family()))
        else:
            # Setting font in qfontcombobox if font is not present in matplotlib
            qfont = self.qtFontDatabase.font(self.axesProps.yLabel['Font Name'], 'normal', self.axesProps.yLabel['Font Size'])
            self.yLabelFontNameCBox.setCurrentFont(qfont)
            self.axes.get_yaxis().get_label().set_fontname(str(qfont.family()))

        # Check whether y-lable position is changed in widget and then change it
        if self.yLabelShowLabelCBox.currentIndex() == 0:
            self.axes.get_yaxis().get_label().set_visible(True)
            self.axes.get_yaxis().set_label_position('left')
        if self.yLabelShowLabelCBox.currentIndex() == 1:
            self.axes.get_yaxis().get_label().set_visible(True)
            self.axes.get_yaxis().set_label_position('right')
        if self.yLabelShowLabelCBox.currentIndex() == 2:
            self.axes.get_yaxis().get_label().set_visible(False)

        # Padding
        self.axes.yaxis.labelpad = self.yLabelPaddingSpinbox.value()

    def displayXTickLabelPropsOnTabWidget(self):
        # Display tick label position
        if self.axesProps.xTickLabel['Label Position'] == 'bottom':
            self.xTickShowTickLabelsCBox.setCurrentIndex(0)
        if self.axesProps.xTickLabel['Label Position'] == 'top':
            self.xTickShowTickLabelsCBox.setCurrentIndex(1)
        if self.axesProps.xTickLabel['Label Position'] == 'both':
            self.xTickShowTickLabelsCBox.setCurrentIndex(2)
        if self.axesProps.xTickLabel['Label Position'] == 'none':
            self.xTickShowTickLabelsCBox.setCurrentIndex(3)

        # Display label rotation
        self.xTickLabelRotationCBox.setValidator(QIntValidator())
        self.xTickLabelRotationCBox.setCurrentText(str(self.axesProps.xTickLabel['Label Rotation']))

        # Setting font in qfontcombobox
        qfont = self.qtFontDatabase.font(self.axesProps.xTickLabel['Font Name'], 'normal', self.axesProps.xTickLabel['Font Size'])
        self.xTickFontNameCBox.setCurrentFont(qfont)

        # setting font size in tab widget
        self.xTickFontSizeCBox.setValidator(QIntValidator())
        self.xTickFontSizeCBox.setCurrentText(str(self.axesProps.xTickLabel['Font Size']))

        # padding
        self.xTickLabelPaddingSpinbox.setValue(self.axesProps.xTickLabel['padding'])

        # number of ticks
        self.xTickIntervalsLineEdit.setValidator(QIntValidator())
        self.xTickIntervalsLineEdit.setText(str(self.axesProps.xTickLabel['Tick Intervals']))

        # Ticks position: top, bottom, both, none
        if self.axesProps.xTickLabel['Tick Position'] == 'both':
            self.xTickShowTicksCBox.setCurrentIndex(0)
        if self.axesProps.xTickLabel['Tick Position'] == 'bottom':
            self.xTickShowTicksCBox.setCurrentIndex(1)
        if self.axesProps.xTickLabel['Tick Position'] == 'top':
            self.xTickShowTicksCBox.setCurrentIndex(2)
        if self.axesProps.xTickLabel['Tick Position'] == 'none':
            self.xTickShowTicksCBox.setCurrentIndex(3)

        # Major tick length
        self.xTickMajorTickLenLineEdit.setValidator(QDoubleValidator())
        self.xTickMajorTickLenLineEdit.setText(str(self.axesProps.xTickLabel['Major Ticks Length']))

        # Major tick width
        self.xTickMajorTickWidthLineEdit.setValidator(QDoubleValidator())
        self.xTickMajorTickWidthLineEdit.setText(str(self.axesProps.xTickLabel['Major Ticks Width']))

        # Show minor ticks
        if self.axesProps.xTickLabel['Minor Ticks Visible']:
            self.xTickShowMinorTicksCBox.setCurrentIndex(1)
        else:
            num_minorticks = len(self.axes.get_xaxis().get_minor_ticks()) / ( len(self.axes.get_xaxis().get_major_ticks()) - 1)
            self.xTickShowMinorTicksCBox.setCurrentIndex( num_minorticks )

        # Minor ticks height
        self.xTickMinorTickLenLineEdit.setValidator(QDoubleValidator())
        if self.axesProps.xTickLabel['Minor Ticks Length'] is not None:
            self.xTickMinorTickLenLineEdit.setText(str(self.axesProps.xTickLabel['Minor Ticks Length']))

        # Minor ticks width
        self.xTickMinorTickWidthLineEdit.setValidator(QDoubleValidator())
        if self.axesProps.xTickLabel['Minor Ticks Width'] is not None:
            self.xTickMinorTickWidthLineEdit.setText(str(self.axesProps.xTickLabel['Minor Ticks Width']))

    def changeXTickLabelPropsOnAxes(self):
        # Change tick label  and ticks position
        ticks = self.axes.get_xaxis().get_major_ticks()

        # Bottom tick label
        if self.xTickShowTickLabelsCBox.currentIndex() == 0:
            for tick in ticks:
                tick.label1On = True
                tick.label2On = False
        # top tick label
        if self.xTickShowTickLabelsCBox.currentIndex() == 1:
            for tick in ticks:
                tick.label1On = False
                tick.label2On = True
        # tick label at both side
        if self.xTickShowTickLabelsCBox.currentIndex() == 2:
            for tick in ticks:
                tick.label1On = True
                tick.label2On = True
        # hide tick label
        if self.xTickShowTickLabelsCBox.currentIndex() == 3:
            for tick in ticks:
                tick.label1On = False
                tick.label2On = False


        ticks = self.axes.get_xaxis().get_major_ticks()
        # Both ticks position
        if self.xTickShowTicksCBox.currentIndex() == 0:
            for tick in ticks:
                tick.tick1On = True
                tick.tick2On = True
        # bottom ticks position
        if self.xTickShowTicksCBox.currentIndex() == 1:
            for tick in ticks:
                tick.tick1On = True
                tick.tick2On = False
        # top ticks position
        if self.xTickShowTicksCBox.currentIndex() == 2:
            for tick in ticks:
                tick.tick1On = False
                tick.tick2On = True
        # hide ticks position
        if self.xTickShowTicksCBox.currentIndex() == 3:
            for tick in ticks:
                tick.tick1On = False
                tick.tick2On = False


        # Change Tick labels rotation
        labels = self.axes.get_xticklabels()
        for label in labels:
            label.set_rotation( int( self.xTickLabelRotationCBox.currentText() ) )


        # Check whether font name is changed in tab-widget and then chenge it in plot
        qfont = self.xTickFontNameCBox.currentFont()
        if self.is_font_exist_in_mpl(str(qfont.family())):
            for label in labels:
                label.set_fontname( str(qfont.family()) )
        else:
            # Setting font in qfontcombobox if font is not present in matplotlib
            qfont = self.qtFontDatabase.font(self.axesProps.xTickLabel['Font Name'], 'normal', self.axesProps.xTickLabel['Font Size'])
            self.xTickFontNameCBox.setCurrentFont(qfont)
            for label in labels:
                label.set_fontname( str(qfont.family()) )

        # font size of labels
        for label in labels:
            label.set_fontsize( int(self.xTickFontSizeCBox.currentText()) )

        #Setting tick-label padding
        self.axes.xaxis.set_tick_params(pad=self.xTickLabelPaddingSpinbox.value())

        #TODO number of ticks
        num_tick_intervals = int ( self.xTickIntervalsLineEdit.text() )
        if num_tick_intervals != self.axesProps.xTickLabel['Tick Intervals']:
            self.axesProps.modify_xtick_intervals(num_tick_intervals)
            self.xTickIntervalsLineEdit.setText( str(len(self.axes.get_xaxis().get_major_ticks())) )
            self.axesProps.xTickLabel['Tick Intervals'] = len(self.axes.get_xaxis().get_major_ticks())

        # Major tick lines: ticklines is a list of 2D Line objects
        ticklines = self.axes.get_xaxis().get_majorticklines()
        if ticklines:
            for tickline in ticklines:
                tickline.set_markersize( float( self.xTickMajorTickLenLineEdit.text()) )
                tickline.set_markeredgewidth( float( self.xTickMajorTickWidthLineEdit.text()) )


        # If minor ticks are enabled, show them and reset their postions according to major ticks
        if self.xTickShowMinorTicksCBox.currentIndex() != 0:

            # Setting minor locator
            minorLocator = AutoMinorLocator( self.xTickShowMinorTicksCBox.currentIndex() + 1 )
            self.axes.get_xaxis().set_minor_locator(minorLocator)

            if not str(self.xTickMinorTickLenLineEdit.text()) and not str(self.xTickMinorTickWidthLineEdit.text()):
                self.axes.tick_params(axis='x', which='minor', bottom=True, top=True)
            if str(self.xTickMinorTickLenLineEdit.text()) and not str(self.xTickMinorTickWidthLineEdit.text()):
                self.axes.tick_params(axis='x', which='minor', length = float( self.xTickMinorTickLenLineEdit.text() ) )
            if not str(self.xTickMinorTickLenLineEdit.text()) and str(self.xTickMinorTickWidthLineEdit.text()):
                self.axes.tick_params(axis='x', which='minor', width =  float( self.xTickMinorTickWidthLineEdit.text() ) )
            if str(self.xTickMinorTickLenLineEdit.text()) and str(self.xTickMinorTickWidthLineEdit.text()):
                self.axes.tick_params(axis='x', which='minor', length = float( self.xTickMinorTickLenLineEdit.text()), width =  float( self.xTickMinorTickWidthLineEdit.text()),  bottom=True, top=True)


            # Update length and width of minor tick lines in widget
            ticklines = self.axes.get_xaxis().get_minorticklines()
            self.axesProps.xTickLabel['Minor Ticks Length'] = ticklines[0].get_markersize()
            self.axesProps.xTickLabel['Minor Ticks Width'] = ticklines[0].get_markeredgewidth()
            self.xTickMinorTickLenLineEdit.setText(str(self.axesProps.xTickLabel['Minor Ticks Length']))
            self.xTickMinorTickWidthLineEdit.setText(str(self.axesProps.xTickLabel['Minor Ticks Width']))

            # Reset positions of minor ticks according to major ticks
            ticks = self.axes.get_xaxis().get_minor_ticks()
            # Bottom tick label
            if self.xTickShowTicksCBox.currentIndex() == 0:
                for tick in ticks:
                    tick.tick1On = True
                    tick.tick2On = True
            # top tick label
            if self.xTickShowTicksCBox.currentIndex() == 1:
                for tick in ticks:
                    tick.tick1On = True
                    tick.tick2On = False
            # tick label at both side
            if self.xTickShowTicksCBox.currentIndex() == 2:
                for tick in ticks:
                    tick.tick1On = False
                    tick.tick2On = True
            # hide tick label
            if self.xTickShowTicksCBox.currentIndex() == 3:
                for tick in ticks:
                    tick.tick1On = False
                    tick.tick2On = False

        # If minor ticks are disabled, remove them
        if self.xTickShowMinorTicksCBox.currentIndex() == 0:
            ticks = self.axes.get_xaxis().get_minor_ticks()
            if ticks:
                for tick in ticks:
                    tick.tick1On = False
                    tick.tick2On = False

    def displayYTickLabelPropsOnTabWidget(self):
        # Display tick label positionx
        if self.axesProps.yTickLabel['Label Position'] == 'left':
            self.yTickShowTickLabelsCBox.setCurrentIndex(0)
        if self.axesProps.yTickLabel['Label Position'] == 'right':
            self.yTickShowTickLabelsCBox.setCurrentIndex(1)
        if self.axesProps.yTickLabel['Label Position'] == 'both':
            self.yTickShowTickLabelsCBox.setCurrentIndex(2)
        if self.axesProps.yTickLabel['Label Position'] == 'none':
            self.yTickShowTickLabelsCBox.setCurrentIndex(3)

        # Display label rotation
        self.yTickLabelRotationCBox.setValidator(QIntValidator())
        self.yTickLabelRotationCBox.setCurrentText(str(self.axesProps.yTickLabel['Label Rotation']))

        # Setting font in qfontcombobox
        qfont = self.qtFontDatabase.font(self.axesProps.yTickLabel['Font Name'], 'normal', self.axesProps.yTickLabel['Font Size'])
        self.yTickFontNameCBox.setCurrentFont(qfont)

        # setting font size in tab widget
        self.yTickFontSizeCBox.setValidator(QIntValidator())
        self.yTickFontSizeCBox.setCurrentText(str(self.axesProps.yTickLabel['Font Size']))

        # padding
        self.yTickLabelPaddingSpinbox.setValue(self.axesProps.yTickLabel['padding'])

        # number of ticks
        self.yTickIntervalsLineEdit.setValidator(QIntValidator())
        self.yTickIntervalsLineEdit.setText(str(self.axesProps.yTickLabel['Tick Intervals']))

        # Ticks position: top, bottom, both, none
        if self.axesProps.yTickLabel['Tick Position'] == 'both':
            self.yTickShowTicksCBox.setCurrentIndex(0)
        if self.axesProps.yTickLabel['Tick Position'] == 'left':
            self.yTickShowTicksCBox.setCurrentIndex(1)
        if self.axesProps.yTickLabel['Tick Position'] == 'right':
            self.yTickShowTicksCBox.setCurrentIndex(2)
        if self.axesProps.yTickLabel['Tick Position'] == 'none':
            self.yTickShowTicksCBox.setCurrentIndex(3)

        # Major tick length
        self.yTickMajorTickLenLineEdit.setValidator(QDoubleValidator())
        self.yTickMajorTickLenLineEdit.setText(str(self.axesProps.yTickLabel['Major Ticks Length']))

        # Major tick width
        self.yTickMajorTickWidthLineEdit.setValidator(QDoubleValidator())
        self.yTickMajorTickWidthLineEdit.setText(str(self.axesProps.yTickLabel['Major Ticks Width']))

        # Show minor ticks
        if self.axesProps.yTickLabel['Minor Ticks Visible']:
            self.yTickShowMinorTicksCBox.setCurrentIndex(1)
        else:
            num_minorticks = len(self.axes.get_yaxis().get_minor_ticks()) / ( len(self.axes.get_yaxis().get_major_ticks()) - 1)
            self.yTickShowMinorTicksCBox.setCurrentIndex( num_minorticks )

        # Minor ticks height
        self.yTickMinorTickLenLineEdit.setValidator(QDoubleValidator())
        if self.axesProps.yTickLabel['Minor Ticks Length'] is not None:
            self.yTickMinorTickLenLineEdit.setText(str(self.axesProps.yTickLabel['Minor Ticks Length']))

        # Minor ticks width
        self.yTickMinorTickWidthLineEdit.setValidator(QDoubleValidator())
        if self.axesProps.yTickLabel['Minor Ticks Width'] is not None:
            self.yTickMinorTickWidthLineEdit.setText(str(self.axesProps.yTickLabel['Minor Ticks Width']))

    def changeYTickLabelPropsOnAxes(self):
        # Change tick label  and ticks position
        ticks = self.axes.get_yaxis().get_major_ticks()

        # Bottom tick label
        if self.yTickShowTickLabelsCBox.currentIndex() == 0:
            for tick in ticks:
                tick.label1On = True
                tick.label2On = False
        # top tick label
        if self.yTickShowTickLabelsCBox.currentIndex() == 1:
            for tick in ticks:
                tick.label1On = False
                tick.label2On = True
        # tick label at both side
        if self.yTickShowTickLabelsCBox.currentIndex() == 2:
            for tick in ticks:
                tick.label1On = True
                tick.label2On = True
        # hide tick label
        if self.yTickShowTickLabelsCBox.currentIndex() == 3:
            for tick in ticks:
                tick.label1On = False
                tick.label2On = False


        ticks = self.axes.get_yaxis().get_major_ticks()
        # Both ticks position
        if self.yTickShowTicksCBox.currentIndex() == 0:
            for tick in ticks:
                tick.tick1On = True
                tick.tick2On = True
        # bottom ticks position
        if self.yTickShowTicksCBox.currentIndex() == 1:
            for tick in ticks:
                tick.tick1On = True
                tick.tick2On = False
        # top ticks position
        if self.yTickShowTicksCBox.currentIndex() == 2:
            for tick in ticks:
                tick.tick1On = False
                tick.tick2On = True
        # hide ticks position
        if self.yTickShowTicksCBox.currentIndex() == 3:
            for tick in ticks:
                tick.tick1On = False
                tick.tick2On = False


        # Change Tick labels rotation
        labels = self.axes.get_yticklabels()
        for label in labels:
            label.set_rotation( int( self.yTickLabelRotationCBox.currentText() ) )


        # Check whether font name is changed in tab-widget and then chenge it in plot
        qfont = self.yTickFontNameCBox.currentFont()
        if self.is_font_exist_in_mpl(str(qfont.family())):
            for label in labels:
                label.set_fontname( str(qfont.family()) )
        else:
            # Setting font in qfontcombobox if font is not present in matplotlib
            qfont = self.qtFontDatabase.font(self.axesProps.yTickLabel['Font Name'], 'normal', self.axesProps.yTickLabel['Font Size'])
            self.yTickFontNameCBox.setCurrentFont(qfont)
            for label in labels:
                label.set_fontname( str(qfont.family()) )

        # font size of labels
        for label in labels:
            label.set_fontsize( int(self.yTickFontSizeCBox.currentText()) )

        #Setting tick-label padding
        self.axes.yaxis.set_tick_params(pad=self.yTickLabelPaddingSpinbox.value())

        #TODO number of ticks
        num_tick_intervals = int ( self.yTickIntervalsLineEdit.text() )
        if num_tick_intervals != self.axesProps.yTickLabel['Tick Intervals']:
            self.axesProps.modify_ytick_intervals(num_tick_intervals)
            self.yTickIntervalsLineEdit.setText( str(len(self.axes.get_yaxis().get_major_ticks())) )
            self.axesProps.yTickLabel['Tick Intervals'] = len(self.axes.get_yaxis().get_major_ticks())

        # Major tick lines: ticklines is a list of 2D Line objects
        ticklines = self.axes.get_yaxis().get_majorticklines()
        if ticklines:
            for tickline in ticklines:
                tickline.set_markersize( float( self.yTickMajorTickLenLineEdit.text()) )
                tickline.set_markeredgewidth( float( self.yTickMajorTickWidthLineEdit.text()) )


        # If minor ticks are enabled, show them and reset their postions according to major ticks
        if self.yTickShowMinorTicksCBox.currentIndex() != 0:

            # Setting minor locator
            minorLocator = AutoMinorLocator( self.yTickShowMinorTicksCBox.currentIndex() + 1 )
            self.axes.get_yaxis().set_minor_locator(minorLocator)

            if not str(self.yTickMinorTickLenLineEdit.text()) and not str(self.yTickMinorTickWidthLineEdit.text()):
                self.axes.tick_params(axis='y', which='minor', bottom=True, top=True)
            if str(self.yTickMinorTickLenLineEdit.text()) and not str(self.yTickMinorTickWidthLineEdit.text()):
                self.axes.tick_params(axis='y', which='minor', length = float( self.yTickMinorTickLenLineEdit.text() ) )
            if not str(self.yTickMinorTickLenLineEdit.text()) and str(self.yTickMinorTickWidthLineEdit.text()):
                self.axes.tick_params(axis='y', which='minor', width =  float( self.yTickMinorTickWidthLineEdit.text() ) )
            if str(self.yTickMinorTickLenLineEdit.text()) and str(self.yTickMinorTickWidthLineEdit.text()):
                self.axes.tick_params(axis='y', which='minor', length = float( self.yTickMinorTickLenLineEdit.text()), width =  float( self.yTickMinorTickWidthLineEdit.text()),  bottom=True, top=True)


            # Update length and width of minor tick lines in widget
            ticklines = self.axes.get_yaxis().get_minorticklines()
            self.axesProps.yTickLabel['Minor Ticks Length'] = ticklines[0].get_markersize()
            self.axesProps.yTickLabel['Minor Ticks Width'] = ticklines[0].get_markeredgewidth()
            self.yTickMinorTickLenLineEdit.setText(str(self.axesProps.yTickLabel['Minor Ticks Length']))
            self.yTickMinorTickWidthLineEdit.setText(str(self.axesProps.yTickLabel['Minor Ticks Width']))

            # Reset positions of minor ticks according to major ticks
            ticks = self.axes.get_yaxis().get_minor_ticks()
            # Bottom tick label
            if self.yTickShowTicksCBox.currentIndex() == 0:
                for tick in ticks:
                    tick.tick1On = True
                    tick.tick2On = True
            # top tick label
            if self.yTickShowTicksCBox.currentIndex() == 1:
                for tick in ticks:
                    tick.tick1On = True
                    tick.tick2On = False
            # tick label at both side
            if self.yTickShowTicksCBox.currentIndex() == 2:
                for tick in ticks:
                    tick.tick1On = False
                    tick.tick2On = True
            # hide tick label
            if self.yTickShowTicksCBox.currentIndex() == 3:
                for tick in ticks:
                    tick.tick1On = False
                    tick.tick2On = False

        # If minor ticks are disabled, remove them
        if self.yTickShowMinorTicksCBox.currentIndex() == 0:
            ticks = self.axes.get_yaxis().get_minor_ticks()
            if ticks:
                for tick in ticks:
                    tick.tick1On = False
                    tick.tick2On = False

    def uponButtonClicked(self, button):
        if button.text() == 'Apply':
            self.changeAxisPropsOnAxes()
        if button.text() == 'Cancel':
            self.done(QDialog.Rejected)
        if button.text() == 'OK':
            self.changeAxisPropsOnAxes()
            self.done(QDialog.Accepted)

    def changeAxisPropsOnAxes(self):
        self.changeXLabelPropsOnAxes()
        self.changeYLabelPropsOnAxes()
        self.changeXTickLabelPropsOnAxes()
        self.changeYTickLabelPropsOnAxes()
        self.canvas.draw()
        self.axesProps.get_from_axes()

    def xMinorTicksController(self):
        ''' To modify minor ticks properties when user enable or disable it
        '''
        if self.xTickShowMinorTicksCBox.currentIndex() != 0:
            self.xTickMinorTickLenLineEdit.setEnabled(True)
            self.xTickMinorTickWidthLineEdit.setEnabled(True)
        if self.xTickShowMinorTicksCBox.currentIndex() == 0:
            self.xTickMinorTickLenLineEdit.setEnabled(False)
            self.xTickMinorTickWidthLineEdit.setEnabled(False)

    def yMinorTicksController(self):
        ''' To modify minor ticks properties when user enable or disable it
        '''
        if self.yTickShowMinorTicksCBox.currentIndex() != 0:
            self.yTickMinorTickLenLineEdit.setEnabled(True)
            self.yTickMinorTickWidthLineEdit.setEnabled(True)
        if self.yTickShowMinorTicksCBox.currentIndex() == 0:
            self.yTickMinorTickLenLineEdit.setEnabled(False)
            self.yTickMinorTickWidthLineEdit.setEnabled(False)



# Dialog box to select and load genomic data
pathToThisUI = os.path.join(PathToUIs, 'genomicDataSelector.ui')
Ui_DialogGenomicsDataSelector, QDialogGenomicsDataSelectorBase = loadUiType(pathToThisUI)
class DialogGenomicsDataSelector(QDialogGenomicsDataSelectorBase, Ui_DialogGenomicsDataSelector):
    def __init__(self, hdf5Handle, requestedBinsize=None):
        super(DialogGenomicsDataSelector, self).__init__()
        self.setupUi(self)

        self.hdf5Handle = None
        self.closeFile = False

        self.data = None
        self.requestedBinsize = requestedBinsize
        self.resolutions = None
        self.dataNames = None
        self.coarse_methods_name = { 'min':'Minimum', 'max':'Maximum', 'median':'Median', 'amean':'Arithmatic Mean', 'gmean':'Geometric Mean', 'hmean':'Harmonic Mean' }
        self.coarse_methods_name_r = { 'Minimum':'min', 'Maximum':'max', 'Median':'median', 'Arithmatic Mean':'amean', 'Geometric Mean':'gmean', 'Harmonic Mean':'hmean' }
        self.selected_data = None
        self.whereToPlot = None

        if isinstance(hdf5Handle, gmlib.genomicsDataHandler.HDF5Handler):
            self.hdf5Handle = hdf5Handle
            if self.hdf5Handle.hdf5 is None:
                self.hdf5Handle.open()
        else:
            self.hdf5Handle = gmlib.genomicsDataHandler.HDF5Handler(hdf5Handle)
            self.hdf5Handle.open()
            self.closeFile = True

        # Add chromosomes in list
        self.chroms = []
        for key in self.hdf5Handle.hdf5:
            self.chroms.append(key)
        self.chromNameListWidget.addItems(gmlib.util.sorted_nicely(self.chroms))

        # Connect to interface
        self.chromNameListWidget.itemSelectionChanged.connect(self.chromosomeSelector)
        self.chromResolutionListWidget.itemSelectionChanged.connect(self.resolutionSelector)
        self.openAbortButton.clicked.connect(self.on_button_clicked)

    def setChromosome(self, chromosome):
        """ Set chromsome in list widget
        """
        if self.hdf5Handle.hasChromosome(chromosome):
            self.chromNameListWidget.setCurrentItem( self.chromNameListWidget.findItems(chromosome, Qt.MatchExactly)[0] )

    def setChromosomeResolution(self, chrom, res):
        """ Set both chromosome and resolution in list widget
        """
        if self.hdf5Handle.hasResolution(chrom, res):
            self.chromNameListWidget.setCurrentItem( self.chromNameListWidget.findItems(chrom, Qt.MatchExactly)[0] )
            self.chromResolutionListWidget.setCurrentItem( self.chromResolutionListWidget.findItems(res, Qt.MatchExactly)[0] )

    def chromosomeSelector(self):
        """ To display/change resolution list when a chromosome is selected by user
        """

        # Removing resolutions list
        if self.resolutions is not None:
            for i in range(len(self.resolutions)):
                dummy = self.chromResolutionListWidget.takeItem(0)
                del dummy

            self.resolutions = None

        # Removing dataNames list
        if self.dataNames is not None:
            for i in range(len(self.dataNames)):
                dummy = self.dataNameListWidget.takeItem(0)
                del dummy

            self.dataNames = None

        # Generating new resolution list
        self.resolutions = []
        currentText = self.chromNameListWidget.currentItem().text()
        for key in self.hdf5Handle.hdf5[currentText]:
            self.resolutions.append(key)

        # Sort from fine to coarse resolution
        binsizes = list( map(gmlib.util.resolutionToBinsize, self.resolutions)  )
        binsizes = sorted( binsizes )
        self.chromResolutionListWidget.addItems( list( map(gmlib.util.binsizeToResolution, binsizes)  ) )

    def resolutionSelector(self):
        """ To display/change resolution list when a chromosome is selected by user
        """

        if self.chromResolutionListWidget.currentItem() is None : return

        # Removing dataNames list
        if self.dataNames is not None:
            for i in range(len(self.dataNames)):
                dummy = self.dataNameListWidget.takeItem(0)
                del dummy

            self.dataNames = None

        # Generating dataNames list
        self.dataNames = []
        currentChrom = self.chromNameListWidget.currentItem().text()
        currentResolution = self.chromResolutionListWidget.currentItem().text()
        for key in self.hdf5Handle.hdf5[currentChrom][currentResolution]:
            if key in self.coarse_methods_name:
                self.dataNames.append(self.coarse_methods_name[key])
            else:
                self.dataNames.append(key)

        self.dataNameListWidget.addItems(sorted(self.dataNames))

    def on_button_clicked(self, button):
        """Do when Open/Abort button is clicked
        """
        if button.text() == "Open":
            if self.dataNameListWidget.currentItem() is not None:
                chroms = self.chromNameListWidget.currentItem().text()
                resolution = self.chromResolutionListWidget.currentItem().text()
                dataType = self.dataNameListWidget.currentItem().text()

                # In case if resolution of map not match with genomic dataset
                if self.requestedBinsize is not None:
                    requestedResolution = gmlib.util.binsizeToResolution(self.requestedBinsize)
                    if resolution != requestedResolution:
                        msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'Selected resolution "{0}" does not match with map resolution "{1}" .'
                                                     .format(resolution, requestedResolution),QMessageBox.Ok, self)
                        msgBox.exec_()
                        msgBox.close()
                        return

                if dataType in self.coarse_methods_name_r:
                    self.selected_data = (chroms, resolution, self.coarse_methods_name_r[dataType])
                else:
                    self.selected_data = (chroms, resolution, dataType)
                self.whereToPlot = self.plotAtComboBox.currentText().lower()
                self.accept()

        if button.text() == "Abort":
            self.reject()


# Dialog box to select and load genomic data from text file
pathToThisUI = os.path.join(PathToUIs, 'textFileSelector.ui')
Ui_DialogTextFileSelector, QDialogTextFileSelectorBase = loadUiType(pathToThisUI)
class DialogTextFileSelector(QDialogTextFileSelectorBase, Ui_DialogTextFileSelector):
    def __init__(self, inputFileName=None):
        super(DialogTextFileSelector, self).__init__()
        self.setupUi(self)

        self.inputFileName = inputFileName
        self.plotPosition = None
        self.title = str( self.titleLineEdit.text() )

        if self.inputFileName is not None:
            self.fileNameLineEdit.setText(self.inputFileName)

        self.openFileButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )

        self.fileNameLineEdit.editingFinished.connect(self.checkInputFile)
        self.titleLineEdit.editingFinished.connect( self.checkTitle )

        self.openFileButton.clicked.connect( self.browseFile )
        self.okCancelButtons.clicked.connect(self.onOkCancelButtonsClicked)

    def browseFile(self):
        """Browse and select input text file
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " text file (*.txt *.dat)"
        path = QFileDialog.getOpenFileName(self, 'Select a text file', '', file_choices)
        if path[0]:
            self.fileNameLineEdit.setText( path[0] )

    def checkInputFile(self):
        """Check input file exist
        """
        inputFile = str( self.fileNameLineEdit.text() )
        if inputFile and not os.path.isfile(inputFile):
            self.fileNameLineEdit.clear()

    def checkTitle(self):
        """Check title of the file"""
        title =  str( self.titleLineEdit.text() )
        if not title:
            msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'No title! Title should be present.', QMessageBox.Ok, self)
            msgBox.exec_()
            msgBox.close()
            self.titleLineEdit.setText( self.title )
        else:
            self.title = title

    def parseDisplayedOptions(self):
        """Parse all options displayed on dialog
        """
        # If input file name store it or return 0
        inputFileName = str( self.fileNameLineEdit.text() )
        if inputFileName:
            self.inputFileName = inputFileName
        else:
            msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'No input text file selected.', QMessageBox.Ok, self)
            msgBox.exec_()
            msgBox.close()
            self.inputFileName = None
            return 0

        # Store plot position
        if self.plotPositionCBox.currentIndex() == 0:
            self.plotPosition = 'bottom'
        else:
            self.plotPosition = 'top'

        return 1

    def onOkCancelButtonsClicked(self, button):
        """Do when final OK or Cancel buttons are clicked by user
        """
        if button.text() == 'OK':
            if self.parseDisplayedOptions():
                self.done(QDialog.Accepted)
            else:
                return
        if button.text() == 'Cancel':
            self.done(QDialog.Rejected)

# Dialog box to choose and load genomic data
pathToThisUI = os.path.join(PathToUIs, 'correlationMaps.ui')
Ui_DialogCorrelationMaps, QDialogCorrelationMapsBase = loadUiType(pathToThisUI)
class DialogCorrelationBetweenMaps(QDialogCorrelationMapsBase, Ui_DialogCorrelationMaps):
    def __init__(self):
        super(DialogCorrelationBetweenMaps, self).__init__()
        self.setupUi(self)

        self.hicmapList = None
        self.firstMap = None
        self.secondMap = None
        self.firstHicmap = None
        self.secondHicmap = None
        self.ignore_triangular = None
        self.diagonal_offset = None
        self.corrType = None
        self.blockSize = None
        self.workDir = None
        self.slideStepSize = None

        # Results
        self.corr = None
        self.aux = None

        # calculating Thread
        self.thread = None
        self.que = None

        # Set default icons for buttons
        self.set_icons()

        # Set default working directory
        self.setDefaultWorkingDirectory()

        # Set loogger handler to a new type
        self.loggerHandler = guiHelpers.pyQtLogHandler()
        self.loggerHandler.messageEmitted.connect( self.loggerPlainTextEdit.appendPlainText ) # Connect loggerHandler to log textbox


        # Set validator for lineedit boxes
        self.diagOffsetLineEdit.setValidator(QIntValidator())
        self.slideStepSizeLineEdit.setValidator(QIntValidator())

        # Connect browse buttons for input ccmap
        self.firstMapBrowsButton.clicked.connect(self.browseFirstMap)
        self.secondMapBrowsButton.clicked.connect(self.browseSecondMap)
        self.workDirBrowseButton.clicked.connect( self.browseWorkingDirectory )

        # check input ccmap line-edits for valid file
        self.firstMapInputLineEdit.editingFinished.connect(self.checkInputMapFiles)
        self.secondMapInputLineEdit.editingFinished.connect(self.checkInputMapFiles)
        self.workDirLineEdit.editingFinished.connect( self.checkWorkDirExist )

        # Connect other remaining buttons
        self.calculateButton.clicked.connect( self.calculate )
        self.stopButton.clicked.connect( self.terminateThread )
        self.closeButton.clicked.connect( self.closeDialog )
        self.saveOutputButton.clicked.connect( self.saveOutput )
        self.clearOutputButton.clicked.connect( self.outputTextEdit.clear )

    def addHicmaps(self, ccmap, title):
        """Add hicmaps
        """
        if self.hicmapList is None:
            self.hicmapList = []
            self.firstMapInputCBox.setItemText(0, title)
            self.secondMapInputCBox.setItemText(0, title)
        else:
            self.firstMapInputCBox.addItem(title)
            self.secondMapInputCBox.addItem(title)

        self.hicmapList.append(ccmap)

    def set_icons(self):
        """Set standard file icon to browse buttons
        """
        self.firstMapBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.secondMapBrowsButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )
        self.workDirBrowseButton.setIcon( self.style().standardIcon(QStyle.SP_DirOpenIcon) )

    def setDefaultWorkingDirectory(self):
        """Set default working directory and self.workDir
        """
        self.workDir = tempfile.gettempdir()
        self.workDirLineEdit.setText(self.workDir)

    def browseWorkingDirectory(self):
        """Browse and choose working directory
        """
        path = QFileDialog.getExistingDirectory(self, 'Select Working Directory', os.getcwd() )
        if path:
            self.workDirLineEdit.setText(path)

    def checkWorkDirExist(self):
        """To check if working directory exists
        """
        workDir = str( self.workDirLineEdit.text() )
        if workDir and os.path.isdir(workDir):
            self.workDir = workDir
        else:
            msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'Directoy: {0} not found!!'.format(workDir), QMessageBox.Ok, self)
            msgBox.exec_()
            msgBox.close()
            self.workDirLineEdit.setText( self.workDir )

    def browseFirstMap(self):
        """Browse and select first ccmap file
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " ccmap file (*.ccmap)"
        path = QFileDialog.getOpenFileName(self, 'Load File', '', file_choices)
        if path[0]:
            self.firstMapInputLineEdit.setText( path[0] )

    def browseSecondMap(self):
        """Browse and select second ccmap file
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " ccmap file (*.ccmap)"
        path = QFileDialog.getOpenFileName(self, 'Load File', '', file_choices)
        if path[0]:
            self.secondMapInputLineEdit.setText( path[0] )

    def checkInputMapFiles(self):
        """Check input ccmap file is exist
        """
        mapFile = str( self.firstMapInputLineEdit.text() )
        if mapFile and not os.path.isfile(mapFile):
            self.firstMapInputLineEdit.clear()

        mapFile = str( self.secondMapInputLineEdit.text() )
        if mapFile and not os.path.isfile(mapFile):
            self.secondMapInputLineEdit.clear()

    def parseDisplayedOptions(self):
        """Parse and store all displayed options in dialog box

        It should return 1 if everything is parsed successfully.

        """

        # Check for first input ccmap
        if self.firstMapInputCBox.currentText() == 'Not Available':
            firstMap = str( self.firstMapInputLineEdit.text() )
            if firstMap:
                self.firstMap = firstMap
            elif self.hicmapList is not None:
                pass
            else:
                msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'Please select first ccmap', QMessageBox.Ok, self)
                msgBox.exec_()
                msgBox.close()
                self.firstMap = None
                return 0

        # Check for second input ccmap
        secondMap = str( self.secondMapInputLineEdit.text() )
        if secondMap:
            self.secondMap = secondMap
        elif self.hicmapList is not None:
            pass
        else:
            msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'Please select second ccmap', QMessageBox.Ok, self)
            msgBox.exec_()
            msgBox.close()
            self.secondMap = None
            return 0

        # Check for ignore_triangular option
        if self.ignoreTriangularCBox.currentIndex() == 0:
            self.ignore_triangular = True
        else:
            self.ignore_triangular = False

        # Check for diagonal offset option
        diagonal_offset = int( self.diagOffsetLineEdit.text() )
        if diagonal_offset:
            self.diagonal_offset = diagonal_offset
        else:
            msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'No offset from Diagonal. Setting it to zero.', QMessageBox.Ok, self)
            msgBox.exec_()
            msgBox.close()
            self.diagonal_offset = 0

        # Check for correlation type option
        if self.corrTypeCBox.currentIndex() == 0:
            self.corrType = 'pearson'
        else:
            self.corrType = 'spearman'

        # Check for block-size option
        if self.blockGBox.isChecked():
            resolution = self.sizeBlockLineEdit.text()
            if resolution:
                try:
                    gmlib.util.resolutionToBinsize(resolution)
                except ValueError:
                    msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'Resolution should be in kb or mb.', QMessageBox.Ok, self)
                    msgBox.exec_()
                    msgBox.close()
                    return 0
                self.blockSize = resolution
        else:
            self.blockSize = None

        # Check for slide step size option
        if int(self.slideStepSizeLineEdit.text()) > 0:
            self.slideStepSize = int(self.slideStepSizeLineEdit.text())
        else:
            msgBox = QMessageBox(QMessageBox.Warning, 'Warning', 'Slide step size should be larger than zero', QMessageBox.Ok, self)
            msgBox.exec_()
            msgBox.close()
            self.slideStepSizeLineEdit.setText('1')
            return 0

        return 1

    def loadCCMAP(self):
        """Load hicmaps from input ccmap files
        """
        if self.firstMap is not None:
            self.firstHicmap = gmlib.ccmap.load_ccmap(self.firstMap)
            self.firstHicmap.make_readable()
        else:
            self.firstHicmap = self.hicmapList[ self.firstMapInputCBox.currentIndex() ]

        if self.secondMap is not None:
            self.secondHicmap = gmlib.ccmap.load_ccmap(self.secondMap)
            self.secondHicmap.make_readable()
        else:
            self.secondHicmap = self.hicmapList[ self.secondMapInputCBox.currentIndex() ]

    def calculate(self, checked):
        """Perform final calculation and display result in output box
        """
        # Check if user has given correct options otherwise return
        if not self.parseDisplayedOptions():    return

        if checked:
            # Clear all text edit boxes at the start of calculation
            self.outputTextEdit.clear()
            self.loggerPlainTextEdit.clear()

            # Load hicmaps
            self.loadCCMAP()

            # Keywords arguments for correlateHicMaps
            kwargs = { 'ignore_triangular':self.ignore_triangular, 'diagonal_offset':self.diagonal_offset, 'corrType':self.corrType,
                        'blockSize':self.blockSize, 'slideStepSize':self.slideStepSize, 'workDir':self.workDir, 'logHandler':self.loggerHandler }

            # Initialize a new thread and connect it
            self.thread = guiHelpers.qtThread(target=gmlib.cmstats.correlateCMaps, args=(self.firstHicmap, self.secondHicmap), kwargs=kwargs)
            self.thread.resultReady.connect( self.finishedCalculation )

            # Start the thread
            self.thread.start()

            # Change states of buttons on GUI
            self.stopButton.setEnabled(True)
            self.calculateButton.setEnabled(False)

    def finishedCalculation(self):
        """When calculation is finished
        """
        # Uncheck the calculate button
        self.calculateButton.setChecked(False)

        # If thread is still alive, particulary when calculation is properly finished but not terminated abruptly
        if self.thread is not None:

            # Write results to ouput text box
            corr, aux = self.thread.results
            if self.blockSize is None:
                self.outputTextEdit.appendPlainText('Correlation: {0}'.format(corr))
                self.outputTextEdit.appendPlainText('p value: {0}'.format(aux))
            else:
                for i in range(len(corr)):
                    if corr[i] != 0:
                        self.outputTextEdit.appendPlainText('{0}\t{1}'.format(int(aux[i]), corr[i]))

            # Close thread and remove from memory
            self.thread.wait()
            self.thread.quit()
            del self.thread
            self.thread = None

            # set states of button on GUI
            self.stopButton.setEnabled(False)
            self.calculateButton.setEnabled(True)

    def terminateThread(self):
        """When calculation is stopped abruptly
        """
        # disable stop button
        self.stopButton.setEnabled(False)
        self.loggerPlainTextEdit.appendPlainText('Wait, trying to stopping calculation...')      # Write message in log box
        QApplication.processEvents()                                                             # Update events on GUI
        if self.thread is None: return
        self.thread.requestInterruption()   # Request for thread interruption. Terminate does not work.
        self.thread.wait()                  # Wait to terminate the thread
        del self.thread                     # remove the thread completely from memory
        self.loggerPlainTextEdit.appendPlainText('Calculation stopped.')    # Update status in log box
        self.thread = None
        self.calculateButton.setEnabled(True)

    def closeDialog(self):
        """Close this dialog
        """
        self.reject()

    def saveOutput(self):
        """Save output as file from output box
        """
        # A dialog box will be displayed to select a text file and path will be stored in the cell
        file_choices = " text file (*.txt *.dat)"
        path = QFileDialog.getSaveFileName(self, 'Save Output File', '', file_choices)
        if path[0]:
            outfile = path[0]
            ext = os.path.splitext(outfile)[1]
            if not (ext == '.txt' or ext == '.dat'):
                outfile += '.txt'
            fout = open(outfile, 'w')
            fout.write(self.outputTextEdit.toPlainText())
            fout.close()


class AxesProperties:
    ''' Used to hold and change axis properties.
    '''
    def __init__(self, axes, xTickLocations=None, xTickLabelTexts=None, yTickLocations=None, yTickLabelTexts=None):
        self.axes = axes
        self.xTickLocations = xTickLocations
        self.xTickLabelTexts = xTickLabelTexts
        self.yTickLocations = yTickLocations
        self.yTickLabelTexts = yTickLabelTexts

        # Initialize axes proeprties dictionary
        self.xLabel = self.init_xlabel_dict()
        self.yLabel = self.init_yLabel_dict()
        self.xTickLabel = self.init_xTickLabel_dict()
        self.yTickLabel = self.init_yTickLabel_dict()

        self.get_from_axes()

    def init_xlabel_dict(self):
        xLabel = dict()
        xLabel['Text'] = None
        xLabel['Show Label'] = None
        xLabel['Font Name'] = None
        xLabel['Font Size'] = None
        xLabel['padding'] = None
        return xLabel

    def init_yLabel_dict(self):
        yLabel = dict()
        yLabel['Text'] = None
        yLabel['Show Label'] = None
        yLabel['Font Name'] = None
        yLabel['Font Size'] = None
        yLabel['padding'] = None
        return yLabel

    def init_xTickLabel_dict(self):
        xTickLabel = dict()
        xTickLabel['Label Position'] = None
        xTickLabel['Label Rotation'] = None
        xTickLabel['Font Name'] = None
        xTickLabel['Font Size'] = None
        xTickLabel['padding'] = None

        xTickLabel['Tick Intervals'] = None
        xTickLabel['Tick Position'] = None
        xTickLabel['Major Ticks Length'] = None
        xTickLabel['Major Ticks Width'] = None
        xTickLabel['Minor Ticks Visible'] = None
        xTickLabel['Minor Ticks Length'] = None
        xTickLabel['Minor Ticks Width'] = None
        return xTickLabel

    def init_yTickLabel_dict(self):
        yTickLabel = dict()
        yTickLabel['Label Position'] = None
        yTickLabel['Label Rotation'] = None
        yTickLabel['Font Name'] = None
        yTickLabel['Font Size'] = None
        yTickLabel['padding'] = None

        yTickLabel['Tick Intervals'] = None
        yTickLabel['Tick Position'] = None
        yTickLabel['Major Ticks Length'] = None
        yTickLabel['Major Ticks Width'] = None
        yTickLabel['Minor Ticks Visible'] = None
        yTickLabel['Minor Ticks Length'] = None
        yTickLabel['Minor Ticks Width'] = None
        return yTickLabel

    def modify_xtick_intervals(self, num_tick_intervals):

        # Do not apply setting when no input is given
        if self.xTickLocations is None: return

        # Original range shown on axis
        view_interval = self.axes.get_xaxis().get_view_interval()

        # Try to change interval with nice tick locations. Sometimes changes, some time not.
        major_locator = MaxNLocator(nbins=num_tick_intervals)
        self.axes.get_xaxis().set_major_locator(major_locator)
        self.axes.get_xaxis().set_view_interval(view_interval[0], view_interval[1])

        # Removing any xticks outside the original intervals
        xticks = list(self.axes.get_xticks())
        while(1):
            if xticks[0] <  view_interval[0]:
                del xticks[0]
            else:
                break

        while(1):
            if xticks[-1] >  view_interval[1]:
                del xticks[-1]
            else:
                break

        # Resetting ticks interval
        self.axes.set_xticks(xticks)

        # Change tick labels if required
        if self.xTickLabelTexts is None:
            self.axes.set_xticklabels(self.axes.get_xticks())
        else:
            # Check whther ticks locations and tick labels are compaitable
            if len(self.xTickLocations) != len(self.xTickLabelTexts):
                raise AssertionError ('Number of x-axis locations you provided do not match with number of xticklabels')

            # Retrive indices of tick labels corresponding to ticks and make list of labels
            xticks = self.axes.get_xticks()
            ticklabels = []
            for xtick in xticks:
                idx = np.argmin( np.abs(self.xTickLocations - xtick) )
                ticklabels.append(self.xTickLabelTexts[idx])

            self.axes.set_xticklabels(ticklabels)

    def modify_ytick_intervals(self, num_tick_intervals):


        # Original range shown on axis
        view_interval = self.axes.get_yaxis().get_view_interval()

        # Try to change interval with nice tick locations. Sometimes changes, some time not.
        major_locator = MaxNLocator(nbins=num_tick_intervals)
        self.axes.get_yaxis().set_major_locator(major_locator)
        self.axes.get_yaxis().set_view_interval(view_interval[0], view_interval[1])

        # Removing any xticks outside the original intervals
        yticks = list(self.axes.get_yticks())
        while(1):
            if yticks[0] <  view_interval[0]:
                del yticks[0]
            else:
                break

        while(1):
            if yticks[-1] >  view_interval[1]:
                del yticks[-1]
            else:
                break

        # Resetting ticks interval
        self.axes.set_yticks(yticks)

        # Do not apply setting when no input is given
        if self.yTickLocations is None: return

        # Change tick labels if required
        if self.yTickLabelTexts is None:
            self.axes.set_yticklabels(self.axes.get_yticks())
        else:
            # Check whther ticks locations and tick labels are compaitable
            if len(self.yTickLocations) != len(self.yTickLabelTexts):
                raise AssertionError ('Number of y-axis locations do not match with number of yticklabels')

            # Retrive indices of tick labels corresponding to ticks and make list of labels
            yticks = self.axes.get_yticks()
            ticklabels = []
            for ytick in yticks:
                idx = np.argmin( np.abs(self.yTickLocations - ytick) )
                ticklabels.append(self.yTickLabelTexts[idx])

            self.axes.set_yticklabels(ticklabels)

    def getXLabelPropsFromAxes(self):
        ''' Get X-label properties from the axes instant
        '''
        self.xLabel['Text'] = self.axes.get_xaxis().get_label_text()
        self.xLabel['Show Label'] = self.axes.get_xaxis().get_label_position()
        if not self.axes.get_xaxis().get_label().get_visible():
            self.xLabel['Show Label'] = 'none'
        self.xLabel['Font Name'] = self.axes.get_xaxis().get_label().get_fontname()
        self.xLabel['Font Size'] = int( self.axes.get_xaxis().get_label().get_fontsize() )
        self.xLabel['padding'] = self.axes.xaxis.labelpad

    def getYLabelPropsFromAxes(self):
        ''' Get Y-label properties from the axes instant
        '''
        self.yLabel['Text'] = self.axes.get_yaxis().get_label_text()
        self.yLabel['Show Label'] = self.axes.get_yaxis().get_label_position()
        if not self.axes.get_yaxis().get_label().get_visible():
            self.yLabel['Show Label'] = 'none'
        self.yLabel['Font Name'] = self.axes.get_yaxis().get_label().get_fontname()
        self.yLabel['Font Size'] = int( self.axes.get_yaxis().get_label().get_fontsize() )
        self.yLabel['padding'] = self.axes.yaxis.labelpad

    def getXTickLabelPropsFromAxes(self):
        # Tick position and Tick-label position
        ticks = self.axes.get_xaxis().get_major_ticks()

        # Getting ticks position
        if ticks[0].tick1On and ticks[0].tick2On:
            self.xTickLabel['Tick Position'] = 'both'
        if ticks[0].tick1On and not ticks[0].tick2On:
            self.xTickLabel['Tick Position'] = 'bottom'
        if not ticks[0].tick1On and ticks[0].tick2On:
            self.xTickLabel['Tick Position'] = 'top'
        if not ticks[0].tick1On and not ticks[0].tick2On:
            self.xTickLabel['Tick Position'] = 'none'

        # Getting tick-labels position
        if ticks[0].label1On and ticks[0].label2On:
            self.xTickLabel['Label Position'] = 'both'
        if ticks[0].label1On and not ticks[0].label2On:
            self.xTickLabel['Label Position'] = 'bottom'
        if not ticks[0].label1On and ticks[0].label2On:
            self.xTickLabel['Label Position'] = 'top'
        if not ticks[0].label1On and not ticks[0].label2On:
            self.xTickLabel['Label Position'] = 'none'

        #Getting tick-label padding
        self.xTickLabel['padding'] = ticks[0].get_pad()


        self.xTickLabel['Minor Ticks Visible'] = bool( self.axes.get_xaxis().get_minor_ticks() )
        # Get number of minor ticks if it is present
        if self.xTickLabel['Minor Ticks Visible']:
            self.xTickLabel['Minor Ticks Visible'] = len(self.axes.get_xaxis().get_minor_ticks()) / ( len(self.axes.get_xaxis().get_major_ticks()) -1 )

        # Major tick lines: ticklines is a list of 2D Line objects
        ticklines = self.axes.get_xaxis().get_majorticklines()
        if ticklines:
            self.xTickLabel['Major Ticks Length'] = ticklines[0].get_markersize()
            self.xTickLabel['Major Ticks Width'] = ticklines[0].get_markeredgewidth()

        # Minor tick lines: ticklines is a list of 2D Line objects
        ticklines = self.axes.get_xaxis().get_minorticklines()
        if ticklines:
            self.xTickLabel['Minor Ticks Length'] = ticklines[0].get_markersize()
            self.xTickLabel['Minor Ticks Width'] = ticklines[0].get_markeredgewidth()

        # Tick labels font-name, font-size and rotations
        if self.xTickLabel['Label Position'] != 'none':
            labels = self.axes.get_xticklabels()
            self.xTickLabel['Font Name'] = labels[0].get_fontname()
            self.xTickLabel['Font Size'] = int( labels[0].get_fontsize() )
            self.xTickLabel['Label Rotation'] = int( labels[0].get_rotation() )

        # Getting number of tick intervals
        if self.xTickLabel['Tick Intervals'] is None:
            self.xTickLabel['Tick Intervals'] = len(self.axes.get_xticks()) - 1

    def getYTickLabelPropsFromAxes(self):
        # Tick position and Tick-label position
        ticks = self.axes.get_yaxis().get_major_ticks()

        # Getting ticks position
        if ticks[0].tick1On and ticks[0].tick2On:
            self.yTickLabel['Tick Position'] = 'both'
        if ticks[0].tick1On and not ticks[0].tick2On:
            self.yTickLabel['Tick Position'] = 'left'
        if not ticks[0].tick1On and ticks[0].tick2On:
            self.yTickLabel['Tick Position'] = 'right'
        if not ticks[0].tick1On and not ticks[0].tick2On:
            self.yTickLabel['Tick Position'] = 'none'

        # Getting tick-labels position
        if ticks[0].label1On and ticks[0].label2On:
            self.yTickLabel['Label Position'] = 'both'
        if ticks[0].label1On and not ticks[0].label2On:
            self.yTickLabel['Label Position'] = 'left'
        if not ticks[0].label1On and ticks[0].label2On:
            self.yTickLabel['Label Position'] = 'right'
        if not ticks[0].label1On and not ticks[0].label2On:
            self.yTickLabel['Label Position'] = 'none'

        #Getting tick-label padding
        self.yTickLabel['padding'] = ticks[0].get_pad()

        self.yTickLabel['Minor Ticks Visible'] = bool( self.axes.get_yaxis().get_minor_ticks() )
        # Get number of minor ticks if it is present
        if self.yTickLabel['Minor Ticks Visible']:
            self.yTickLabel['Minor Ticks Visible'] = len(self.axes.get_yaxis().get_minor_ticks()) / ( len(self.axes.get_yaxis().get_major_ticks()) -1 )

        # Major tick lines: ticklines is a list of 2D Line objects
        ticklines = self.axes.get_yaxis().get_majorticklines()
        if ticklines:
            self.yTickLabel['Major Ticks Length'] = ticklines[0].get_markersize()
            self.yTickLabel['Major Ticks Width'] = ticklines[0].get_markeredgewidth()

        # Minor tick lines: ticklines is a list of 2D Line objects
        ticklines = self.axes.get_yaxis().get_minorticklines()
        if ticklines:
            self.yTickLabel['Minor Ticks Length'] = ticklines[0].get_markersize()
            self.yTickLabel['Minor Ticks Width'] = ticklines[0].get_markeredgewidth()

        # Tick labels font-name, font-size and rotations
        if self.yTickLabel['Label Position'] != 'none':
            labels = self.axes.get_yticklabels()
            self.yTickLabel['Font Name'] = labels[0].get_fontname()
            self.yTickLabel['Font Size'] = int( labels[0].get_fontsize() )
            self.yTickLabel['Label Rotation'] = int( labels[0].get_rotation() )

        # Getting number of tick intervals
        if self.yTickLabel['Tick Intervals'] is None:
            self.yTickLabel['Tick Intervals'] = len(self.axes.get_yticks()) - 1

    def setXLabelPropsToAxes(self):
        ''' Get X-label properties from the axes instant
        '''
        self.axes.get_xaxis().set_label_text(self.xLabel['Text'])
        if self.xLabel['Show Label'] == 'none':
            self.axes.get_xaxis().get_label().set_visible(False)
        else:
            self.axes.get_xaxis().set_label_position(self.xLabel['Show Label'])
        self.axes.get_xaxis().get_label().set_fontname(self.xLabel['Font Name'])
        self.axes.get_xaxis().get_label().set_fontsize(self.xLabel['Font Size'])
        self.axes.xaxis.labelpad = self.xLabel['padding']

    def setYLabelPropsToAxes(self):
        self.axes.get_yaxis().set_label_text(self.yLabel['Text'])
        if self.yLabel['Show Label'] == 'none':
            self.axes.get_yaxis().get_label().set_visible(False)
        else:
            self.axes.get_yaxis().set_label_position(self.yLabel['Show Label'])
        self.axes.get_yaxis().get_label().set_fontname(self.yLabel['Font Name'])
        self.axes.get_yaxis().get_label().set_fontsize(self.yLabel['Font Size'])
        self.axes.yaxis.labelpad = self.yLabel['padding']

    def setXTickLabelPropsToAxes(self):
        # Setting minor ticks
        if self.xTickLabel['Minor Ticks Visible']:
            # Setting minor locator
            minorLocator = AutoMinorLocator( self.xTickLabel['Minor Ticks Visible'] + 1 )
            self.axes.get_xaxis().set_minor_locator(minorLocator)

        # Tick position and Tick-label position
        ticks = self.axes.get_xaxis().get_major_ticks()
        if self.xTickLabel['Minor Ticks Visible']:
            ticks = ticks + self.axes.get_xaxis().get_minor_ticks()

        if self.xTickLabel['Tick Position'] == 'both':
            for tick in ticks:
                tick.tick1On = True
                tick.tick2On = True
        if self.xTickLabel['Tick Position'] == 'bottom':
            for tick in ticks:
                tick.tick1On = True
                tick.tick2On = False
        if self.xTickLabel['Tick Position'] == 'top':
            for tick in ticks:
                tick.tick1On = False
                tick.tick2On = True
        if self.xTickLabel['Tick Position'] == 'none':
            for tick in ticks:
                tick.tick1On = False
                tick.tick2On = False

        ticks = self.axes.get_xaxis().get_major_ticks()
        if self.xTickLabel['Label Position'] == 'both':
            for tick in ticks:
                tick.label1On = True
                tick.label2On = True
        if self.xTickLabel['Label Position'] == 'bottom':
            for tick in ticks:
                tick.label1On = True
                tick.label2On = False
        if self.xTickLabel['Label Position'] == 'top':
            for tick in ticks:
                tick.label1On = False
                tick.label2On = True
        if self.xTickLabel['Label Position'] == 'none':
            for tick in ticks:
                tick.label1On = False
                tick.labelOn = False

        #Setting tick-label padding
        self.axes.xaxis.set_tick_params(pad=self.xTickLabel['padding'])

        # Major tick lines: ticklines is a list of 2D Line objects
        ticklines = self.axes.get_xaxis().get_majorticklines()
        if ticklines:
            for tickline in ticklines:
                tickline.set_markersize(self.xTickLabel['Major Ticks Length'])
                tickline.set_markeredgewidth(self.xTickLabel['Major Ticks Width'])


        # Tick labels font-name, font-size and rotations
        labels = self.axes.get_xticklabels()
        for label in labels:
            label.set_fontname(self.xTickLabel['Font Name'])
            label.set_fontsize(self.xTickLabel['Font Size'])
            label.set_rotation(self.xTickLabel['Label Rotation'])

        # Getting number of tick intervals
        self.modify_xtick_intervals(self.xTickLabel['Tick Intervals'])

        # Minor tick lines: ticklines is a list of 2D Line objects
        ticklines = self.axes.get_xaxis().get_minorticklines()
        if ticklines:
            for tickline in ticklines:
                tickline.set_markersize(self.xTickLabel['Minor Ticks Length'])
                tickline.set_markeredgewidth(self.xTickLabel['Minor Ticks Width'])

    def setYTickLabelPropsToAxes(self):
        # Setting minor ticks
        if self.yTickLabel['Minor Ticks Visible']:
            # Setting minor locator
            minorLocator = AutoMinorLocator( self.yTickLabel['Minor Ticks Visible'] + 1 )
            self.axes.get_yaxis().set_minor_locator(minorLocator)

        # Tick position and Tick-label position
        ticks = self.axes.get_yaxis().get_major_ticks()
        if self.yTickLabel['Minor Ticks Visible']:
            ticks = ticks + self.axes.get_yaxis().get_minor_ticks()

        if self.yTickLabel['Tick Position'] == 'both':
            for tick in ticks:
                tick.tick1On = True
                tick.tick2On = True
        if self.yTickLabel['Tick Position'] == 'left':
            for tick in ticks:
                tick.tick1On = True
                tick.tick2On = False
        if self.yTickLabel['Tick Position'] == 'right':
            for tick in ticks:
                tick.tick1On = False
                tick.tick2On = True
        if self.yTickLabel['Tick Position'] == 'none':
            for tick in ticks:
                tick.tick1On = False
                tick.tick2On = False

        ticks = self.axes.get_yaxis().get_major_ticks()
        if self.yTickLabel['Label Position'] == 'both':
            for tick in ticks:
                tick.label1On = True
                tick.label2On = True
        if self.yTickLabel['Label Position'] == 'left':
            for tick in ticks:
                tick.label1On = True
                tick.label2On = False
        if self.yTickLabel['Label Position'] == 'right':
            for tick in ticks:
                tick.label1On = False
                tick.label2On = True
        if self.yTickLabel['Label Position'] == 'none':
            for tick in ticks:
                tick.label1On = False
                tick.labelOn = False

        #Setting tick-label padding
        self.axes.yaxis.set_tick_params(pad=self.yTickLabel['padding'])

        # Major tick lines: ticklines is a list of 2D Line objects
        ticklines = self.axes.get_yaxis().get_majorticklines()
        if ticklines:
            for tickline in ticklines:
                tickline.set_markersize(self.yTickLabel['Major Ticks Length'])
                tickline.set_markeredgewidth(self.yTickLabel['Major Ticks Width'])


        # Tick labels font-name, font-size and rotations
        labels = self.axes.get_yticklabels()
        for label in labels:
            label.set_fontname(self.yTickLabel['Font Name'])
            label.set_fontsize(self.yTickLabel['Font Size'])
            label.set_rotation(self.yTickLabel['Label Rotation'])

        # Getting number of tick intervals
        self.modify_ytick_intervals(self.yTickLabel['Tick Intervals'])

        # Minor tick lines: ticklines is a list of 2D Line objects
        ticklines = self.axes.get_yaxis().get_minorticklines()
        if ticklines:
            for tickline in ticklines:
                tickline.set_markersize(self.yTickLabel['Minor Ticks Length'])
                tickline.set_markeredgewidth(self.yTickLabel['Minor Ticks Width'])

    def get_from_axes(self):
        # getting axes properties
        self.getXLabelPropsFromAxes()
        self.getYLabelPropsFromAxes()
        self.getXTickLabelPropsFromAxes()
        self.getYTickLabelPropsFromAxes()

    def set_to_axes(self):
        self.setXLabelPropsToAxes()
        self.setYLabelPropsToAxes()
        self.setXTickLabelPropsToAxes()
        self.setYTickLabelPropsToAxes()


class menuRightClick(QMenu):
    def __init__(self, ):
        super(menuRightClick, self).__init__()
        self.actionList = dict()
        self.addSubMenus()

    def addSubMenus(self):
        self.actionXlabel = self.addAction("Edit X-axis Label...")
        self.actionList[0] = self.actionXlabel
        self.actionYlabel = self.addAction("Edit Y-axis Label...")
        self.actionList[1] = self.actionYlabel
        self.actionXtickLabel = self.addAction("Edit X-ticks Label...")
        self.actionList[2] = self.actionXtickLabel
        self.actionYtickLabel = self.addAction("Edit Y-ticks Label...")
        self.actionList[3] = self.actionYtickLabel

# Dialog box to change page size by custom length
# Dialog box to choose and load genomic data
pathToThisUI = os.path.join(PathToUIs, 'aboutBrowser.ui')
Ui_aboutBrowserDialog, aboutBrowserDialogBase = loadUiType(pathToThisUI)
class aboutBrowserDialog(aboutBrowserDialogBase, Ui_aboutBrowserDialog):
    def __init__(self, parent=None):
        super(aboutBrowserDialog, self).__init__(parent=parent)
        self.setupUi(self)
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.okButton.clicked.connect( self.close )
        self.textBrowser.setOpenExternalLinks(True)


# Dialog box to change page size by custom length
# Dialog box to choose and load genomic data
pathToThisUI = os.path.join(PathToUIs, 'userColorMapDialog.ui')
Ui_DialogUserColorMap, DialogUserColorMapBase = loadUiType(pathToThisUI)
class DialogUserColorMap(DialogUserColorMapBase, Ui_DialogUserColorMap):
    def __init__(self, parent=None):
        super(DialogUserColorMap, self).__init__(parent=parent)
        self.setupUi(self)
        self.setAttribute(Qt.WA_DeleteOnClose)

        self.qdSpinBoxToButton = dict()
        """ dictionary from QDoubleSpinBox to QPushButton in QTableWidget"""

        self.buttonToqdSpinBox = dict()
        """ dictionary from QPushButton to QDoubleSpinBox in QTableWidget"""

        self.resultColorInfo = None
        """ Contain resulted color information dictionary

            Structure of colorInfo dictionary

            colorInfo ----- 'name'
                |
                └------ 'colors'
                            |
                            ├------ values -> color
                            |
                            ├------ values -> color
                            |
                            :
                            └------ values -> color
        """

        self.standardColorMaps = None
        """ Contain lsit of standard color maps """

        self.userColorMaps = dict()
        """ Dictionary of user defined color maps. It is used to reload table
        when user select different color maps in QComboBox."""

        # Hide tabbars
        tabbars = self.addModifyWidget.findChildren(QTabBar)
        for tabbar in tabbars:
            tabbar.hide()

        # Initialize table
        self.tableWidget.verticalHeader().setMinimumWidth(30)
        self.initTable()
        self.previewColorMap()

        self.addModifyCBox.currentIndexChanged.connect( self.addModifyWidget.setCurrentIndex )
        self.addRowTableButton.clicked.connect( self.addRowToTable )
        self.removeRowTableButton.clicked.connect( self.removeRowFromTable )
        self.saveTableButton.clicked.connect( self.saveColorMap )
        self.loadTableButton.clicked.connect( self.loadColorMapFromJson )
        self.cmapListCBox.currentIndexChanged.connect( self.loadColorInfoFromSpinbox )
        self.cmapListCBox.activated.connect( self.loadColorInfoFromSpinbox )

        self.helpButton.clicked.connect( QWhatsThis.enterWhatsThisMode )

        self.okCancelButtonBox.rejected.connect( self.reject )
        self.okCancelButtonBox.accepted.connect( self.dialogAccepted )

    def initTable(self):
        """ Initialize table with three colors
        """
        colors = ['#0000FF', '#00FF00', '#FF0000']
        value = [0, 0.5, 1]
        for r in [0, 1, 2]:
            self.tableWidget.insertRow(r)
            if self.tableWidget.item(r, 0) is None:
                self.addContentsToRow(r, color=colors[r], value=value[r])

    def addRowToTable(self):
        """Add a new row with QDoubleSpinBox and QPushButton
        at the end of table widget
        """
        row = self.tableWidget.rowCount()
        self.tableWidget.insertRow(row)
        self.addContentsToRow(row)

    def addContentsToRow(self, row, color=None, value=None):
        """ Add QDoubleSpinBox and QPushButton to cells of input row in table
        """
        # Add Spin box
        qsbox = QDoubleSpinBox(self)
        qsbox.setRange(0, 1)
        qsbox.setDecimals(3)
        qsbox.setSingleStep(0.05)
        if value is not None:
            qsbox.setValue(value)
        qsbox.valueChanged.connect( self.previewColorMap )
        self.tableWidget.setCellWidget( row, 0, qsbox )

        # Add Button
        qButton = QPushButton(self)
        qButton.setFlat(True)
        self.tableWidget.setCellWidget( row, 1, qButton )
        if color is not None:
            style = 'background-color: {0}; border: none;'.format(color)
            qButton.setStyleSheet(style)
        qButton.clicked.connect( self.chooseColorByButton )

        # Modify dictionaries
        self.qdSpinBoxToButton[qsbox] = qButton
        self.buttonToqdSpinBox[qButton] = qsbox

    def removeRowFromTable(self):
        """Remove a selected row from color table
        """
        # Get selected cell
        row, col = getSelectedRowColumnFromTable(self.tableWidget)

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

        # determine QDoubleSpinBox and QPushButton
        qsbox = self.tableWidget.cellWidget(row, 0)
        qButton = self.tableWidget.cellWidget(row, 1)

        # Remove row, QDoubleSpinBox, and QPushButton
        self.qdSpinBoxToButton.pop(qsbox)
        self.buttonToqdSpinBox.pop(qButton)
        self.tableWidget.removeRow(row)

    def clearTable(self):
        """ Clear entire table
        """
        while self.tableWidget.rowCount() != 0:
            # determine QDoubleSpinBox and QPushButton
            qsbox = self.tableWidget.cellWidget(0, 0)
            qButton = self.tableWidget.cellWidget(0, 1)

            # Remove row, QDoubleSpinBox, and QPushButton
            self.qdSpinBoxToButton.pop(qsbox)
            self.buttonToqdSpinBox.pop(qButton)
            self.tableWidget.removeRow(0)

    def setColorMapList(self, colormapList):
        """ Set color map list from browser
        """
        self.standardColorMaps = []
        self.cmapListCBox.blockSignals(True)   # Block signal fo QComboBox

        for colormap in colormapList:
            # Check if it is a LinearSegmentedColormap.
            # If a user added the colormap it is  LinearSegmentedColormap
            # otherwise a string
            if isinstance(colormapList[colormap], mplColors.LinearSegmentedColormap):
                # Find index in QComboBox
                idx = self.cmapListCBox.findText(colormapList[colormap].name, Qt.MatchExactly)

                # Only add when this colormap when it is not present in QComboBox
                if idx == -1:
                    colorInfo = segmentDataColorMapToColorInfo(colormapList[colormap])
                    self.userColorMaps[colormapList[colormap].name] = colorInfo
                    self.cmapListCBox.addItem(colormapList[colormap].name)
            else:
                self.standardColorMaps.append(colormapList[colormap])

        self.cmapListCBox.blockSignals(False)

    def previewColorMap(self):
        """ Update preview of color map
        """
        colorDict = self.getColorDictFromTable()
        style = 'background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0'
        for value in colorDict:
            style += ', stop:{0} {1}'.format(value, colorDict[value])
        style += ');'
        self.previewColorMapWidget.setStyleSheet(style)

    def getColorDictFromTable(self):
        """ Get color dictionary from table
        """
        valuesToSBox = dict()
        colorDict = dict()
        for qsbox in self.qdSpinBoxToButton:
            valuesToSBox[qsbox.value()] = qsbox

        values = list(valuesToSBox.keys())
        values = np.sort(values)

        for i in range(len(values)):
            button = self.qdSpinBoxToButton[valuesToSBox[values[i]]]
            color = button.palette().color(button.backgroundRole())
            colorDict[values[i]] = color.name()

        if values[0] != 0.0:
            colorDict[0.0] = colorDict[values[0]]

        if values[-1] != 1.0:
            colorDict[1.0] = colorDict[values[-1]]

        return colorDict

    def chooseColorByButton(self):
        """ Choose color of a button and set it
        """
        if self.sender() == 0:
            return
        # Get the button from where signal was originated
        button = self.sender()

        # Get the background-color of button
        qcolor = button.palette().color(button.backgroundRole())

        # QColorDialog open here
        colorDialog = QColorDialog(qcolor, self)
        colorDialog.setWindowTitle("Choose a Color")
        colorDialog.exec_()

        if colorDialog.result() == QDialog.Accepted:
            # Get new color from the user
            pickedColor = colorDialog.selectedColor()

            # Set new color to button
            style = 'background-color: {0}; border: none;'.format(pickedColor.name())
            button.setStyleSheet(style)

            # Update preview color map
            self.previewColorMap()

    def getColorInformation(self):
        """ Get color-information dictionary by parsing displayed widgets.

        Structure of colorInfo dictionary

        colorInfo ----- 'name'
            |
            └------ 'colors'
                        |
                        ├------ values -> color
                        |
                        ├------ values -> color
                        |
                        :
                        └------ values -> color


        """
        if self.addModifyCBox.currentText() == 'Add':
            name = self.cmapNameLineEdit.text()
            if not name:
                msgBox = QMessageBox(self)
                msgBox.setWindowModality(Qt.WindowModal)
                msgBox.setWindowTitle('Warning')
                msgBox.setIcon(QMessageBox.Information)
                msgBox.setText(' No name !!')
                msgBox.setStandardButtons(QMessageBox.Cancel)
                msgBox.exec_()
                msgBox.close()
                self.cmapNameLineEdit.setFocus(True)
                return

        if self.addModifyCBox.currentText() == 'Modify':
            name = self.cmapListCBox.currentText()

        if name in self.standardColorMaps:
            msgBox = QMessageBox(self)
            msgBox.setWindowModality(Qt.WindowModal)
            msgBox.setWindowTitle('Warning')
            msgBox.setIcon(QMessageBox.Information)
            msgBox.setText(' colormap "{0}" already present in standard maps.\n Please choose another name.'.format(name))
            msgBox.setStandardButtons(QMessageBox.Cancel)
            msgBox.exec_()
            msgBox.close()
            self.cmapNameLineEdit.setFocus(True)
            return

        colorInfo = dict()
        colorInfo['name'] = name
        colorInfo['colors'] = self.getColorDictFromTable()


        return colorInfo

    def saveColorMap(self):
        """ Save the color map as a json file

        Format of file is similar to colorInfo dictionary.

        colorInfo ----- 'name'
            |
            └------ 'colors'
                        |
                        ├------ values -> color
                        |
                        ├------ values -> color
                        |
                        :
                        └------ values -> color

        """
        # Get color information from table
        colorInfo = self.getColorInformation()
        if colorInfo is None:
            return

        # A dialog box will be displayed to select a file and path will be stored in the cell
        file_choices = " json file (*.json);;All Files(*.*)"
        path = QFileDialog.getSaveFileName(self, 'Select or Create File', '', file_choices, options=QFileDialog.DontConfirmOverwrite)
        if path[0]:
            outDir = os.path.dirname( path[0] )
            baseName = os.path.basename( path[0] )
            extension = os.path.splitext( baseName )[1]

            if not (extension == '.json'):
                outName = os.path.join(outDir, baseName+'.json')
            else:
                outName = path[0]

            fout =  open( outName, "w" )
            json.dump(colorInfo, fout, indent=4, separators=(',', ':'))
            fout.close()

    def loadColorMapFromJson(self):
        """ Load color map from a json file
        """
        # A dialog box will be displayed to select a json file
        file_choices = " json file (*.json);;All Files(*.*)"
        path = QFileDialog.getOpenFileName(self, 'Load File', '', file_choices)
        if not path[0]:
            return

        inputFile = path[0]

        try:
            fin = open( inputFile, "r" )
        except:
            return

        try:
            with open( inputFile, "r" ) as fin:
                jsonColorInfo = json.load( fin )
        except:
            msgBox = QMessageBox(self)
            msgBox.setWindowModality(Qt.WindowModal)
            msgBox.setWindowTitle('Warning')
            msgBox.setIcon(QMessageBox.Information)
            msgBox.setText(' Not able to read from file.')
            msgBox.setStandardButtons(QMessageBox.Cancel)
            msgBox.exec_()
            msgBox.close()
            self.cmapNameLineEdit.setFocus(True)
            return

        # Changes value to float
        colorInfo = dict()
        try:
            colorInfo['name'] = jsonColorInfo['name']
            colorInfo['colors'] = dict()
            for key in jsonColorInfo['colors']:

                # Check if color is readable
                color = jsonColorInfo['colors'][key]
                if not mplColors.is_color_like(color):
                    raise(ValueError)

                colorInfo['colors'][float(key)] = color

        except:
            msgBox = QMessageBox(self)
            msgBox.setWindowModality(Qt.WindowModal)
            msgBox.setWindowTitle('Error')
            msgBox.setIcon(QMessageBox.Information)
            fmt = """File format is not compatible. See Below an example:

                        {
                            "colors":{
                                "0.0":"#0000ff",
                                "0.5":"#00ff00",
                                "1.0":"#ff0000"
                            },
                            "name":"blue-green-red"
                        """
            msgBox.setText(fmt)
            msgBox.setStandardButtons(QMessageBox.Cancel)
            msgBox.exec_()
            msgBox.close()
            return

        self.loadColorInfoToTable(colorInfo)

    def loadColorInfoFromSpinbox(self):
        """ Change colormap preview and load table freshly
        """
        name = self.cmapListCBox.currentText()
        if name:
            colorInfo = self.userColorMaps[name]
            self.clearTable()  # Clear table
            self.loadColorInfoToTable(colorInfo)  # Load new table

    def loadColorInfoToTable(self, colorInfo):
        """ Load color information from dictonary to table widget
        """
        # Set name
        self.cmapNameLineEdit.setText(colorInfo['name'])

        # Get and sort stop values
        values = list( colorInfo['colors'].keys() )
        values = np.sort(values)

        # Add row in table and set QDoubleSpinBox and QPushButton
        for i in range(len(values)):
            r = self.tableWidget.rowCount()
            if i > r-1 :
                self.addRowToTable()

            qsbox = self.tableWidget.cellWidget(i, 0)
            qsbox.setValue( values[i] )

            color = colorInfo['colors'][values[i]]
            qButton = self.tableWidget.cellWidget(i, 1)
            style = 'background-color: {0}; border: none;'.format(color)
            qButton.setStyleSheet(style)

            # Update colormap preview
            self.previewColorMap()

    def dialogAccepted(self):
        """ When OK button is clicked
        """
        colorInfo = self.getColorInformation()
        if colorInfo is not None:
            self.resultColorInfo = colorInfo
            self.accept()


def add_colormaps_to_combobox(cbox):
    cmap_lists = [ 'afmhot_r', 'autumn_r', 'bone_r', 'cool_r', 'copper_r',
                  'gist_heat_r', 'gray_r', 'hot_r', 'pink_r', 'spring_r', 'summer_r', 'winter_r', 'Blues',
                  'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu', 'PuBuGn', 'PuRd',
                  'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd',
                  'afmhot', 'autumn', 'bone', 'cool', 'copper', 'gist_heat', 'gray', 'hot', 'pink',
                   'spring', 'summer', 'winter', ]

    cmap_icons = ['afmhot_r.ico', 'autumn_r.ico', 'bone_r.ico', 'cool_r.ico', 'copper_r.ico',
                  'gist_heat_r.ico', 'gray_r.ico', 'hot_r.ico', 'pink_r.ico', 'spring_r.ico', 'summer_r.ico', 'winter_r.ico', 'blues.ico',
                  'blue_green.ico', 'blue_purple.ico', 'green_blue.ico', 'greens.ico', 'greys.ico', 'oranges.ico', 'oranges_red.ico', 'purple_blue.ico', 'purple_blue_green.ico',
                  'purple_red.ico', 'purples.ico', 'red_purples.ico', 'reds.ico', 'yellow_green.ico', 'yellow_green_blue.ico', 'yellow_orange_brown.ico', 'yellow_orange_red.ico',
                  'afmhot.ico', 'autumn.ico', 'bone.ico', 'cool.ico', 'copper.ico', 'gist_heat.ico', 'gray.ico', 'hot.ico', 'pink.ico',
                  'spring.ico', 'summer.ico', 'winter.ico',  ]

    icon_path = os.path.join(DirToThisScript, 'icons')
    icon_path = os.path.join(icon_path, 'cmaps')

    cmap = dict()
    for i in range(len(cmap_lists)):
        cmap[i] = cmap_lists[i]

    for i in range(len(cmap_lists)):
        qicon = QIcon()
        icon_file = os.path.join(icon_path, cmap_icons[i])
        qicon.addFile(icon_file)
        cbox.addItem(qicon, cmap_lists[i])

    cbox.setCurrentIndex(6)

    return cmap

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

def segmentDataColorMapToColorInfo(colormap):
    """ make a colormap information dictionary from
    segmented Data colormap object.

    Structure of ouput colorInfo dictionary

    colorInfo ----- 'name'
        |
        └------ 'colors'
                    |
                    ├------ values -> color
                    |
                    ├------ values -> color
                    |
                    :
                    └------ values -> color

    """
    colorInfo = dict()
    colorInfo['name'] = colormap.name
    colorInfo['colors'] = dict()
    for i in range(len(colormap._segmentdata['blue'])):
        value = colormap._segmentdata['blue'][i][0]
        r = colormap._segmentdata['red'][i][1]
        g = colormap._segmentdata['green'][i][1]
        b = colormap._segmentdata['blue'][i][1]
        a = colormap._segmentdata['alpha'][i][1]
        colorInfo['colors'][value] = mplColors.rgb2hex((r, g, b, a))

    return colorInfo

def colorInfoToSegmentDataColorMap(colorInfo):
    """ make a segmented Data colormap object from a colormap dictionary
    """
    cdict = colorInfo['colors']
    keys = list( cdict.keys())
    keys = np.sort(keys)
    color_list = []
    for i in range(len(keys)):
        color_list.append( (keys[i], cdict[keys[i]]) )

    colormap = mplColors.LinearSegmentedColormap.from_list(colorInfo['name'], color_list)
    return colormap

def add_external_colormap_to_combobox(cbox, colorMapList, colorInfo):
    """ Add user defined colormap
    """
    # Generate new pixmap for icons
    width = 146
    height = 12
    pixmap = QPixmap(width, height);
    pixmap.fill(Qt.transparent);

    # Generate gradient and fill it in pixmap
    gradient = QLinearGradient(0, 0, width, height);
    for v in colorInfo['colors']:
        qtcolor = QColor(colorInfo['colors'][v])
        gradient.setColorAt(v, qtcolor)
    painter = QPainter(pixmap)
    painter.fillRect(0, 0, width, height, gradient)
    painter.end()

    # Add gradient into combobox and cmap dictionary
    qicon = QIcon(pixmap)
    colorMapList[cbox.count()] = colorInfoToSegmentDataColorMap(colorInfo)
    cbox.addItem(qicon, colorInfo['name'])


def change_colormap_icon_to_combobox(cbox, index, colorInfo):
    """ Add user defined colormap
    """
    # Generate new pixmap for icons
    width = 146
    height = 12
    pixmap = QPixmap(width, height);
    pixmap.fill(Qt.transparent);

    # Generate gradient and fill it in pixmap
    gradient = QLinearGradient(0, 0, width, height);
    for v in colorInfo['colors']:
        qtcolor = QColor(colorInfo['colors'][v])
        gradient.setColorAt(v, qtcolor)
    painter = QPainter(pixmap)
    painter.fillRect(0, 0, width, height, gradient)
    painter.end()

    # Add gradient into combobox and cmap dictionary
    qicon = QIcon(pixmap)
    cbox.setItemIcon(index, qicon)

def scaleMatrix(matrix, vmin, vmax):
    # Masked the matrix for zero values
    maskedMatrix = np.ma.masked_equal(matrix, 0.0, copy=False)

    scaledMatrix = (maskedMatrix - maskedMatrix.min()) / (maskedMatrix.max() - maskedMatrix.min())

    scaledMatrix = scaledMatrix * (vmax - vmin) + vmin

    return scaledMatrix.filled(0.0)


def add_icon_to_widget(widget, icon_file):
    qicon = QIcon()
    icon_path = os.path.join(DirToThisScript, 'icons')
    icon_file = os.path.join(icon_path, icon_file)
    qicon.addFile(icon_file)
    widget.setIcon(qicon)

def HLine():
    frame = QFrame()
    frame.setFrameShape(QFrame.HLine)
    frame.setFrameShadow(QFrame.Sunken)
    return frame

def VLine():
    frame = QFrame()
    frame.setFrameShape(QFrame.VLine)
    frame.setFrameShadow(QFrame.Sunken)
    return frame

def get_interpolation_dict():
    method = ['none', 'nearest', 'bilinear', 'sinc', 'lanczos', 'bicubic', 'gaussian', 'spline16', 'spline36', 'hanning', 'hamming',
              'hermite', 'kaiser', 'quadric', 'catrom', 'bessel', 'mitchell']

    interpolation = dict()
    for i in range(len(method)):
        interpolation[i] = method[i]

    return interpolation

def showWarningMessageBox(msg, qwidget):
    msgBox = QMessageBox(QMessageBox.Warning, 'Warning', msg, QMessageBox.Ok, qwidget)
    msgBox.exec_()
    msgBox.close()
