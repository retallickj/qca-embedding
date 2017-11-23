#!/usr/bin/env python

# -----------------------------------
# Name: gui_settings.py
# Desc: GUI settings for embedder application
# Author: Jake Retallick
# Created: 2015.11.25
# Modified: 2015.11.25
# Licence: Copyright 2015
# -----------------------------------

import sys

# import Qt based on installed version
from gui.pyqt_import import importPyQt
QtGui, QtCore = importPyQt('QtGui', 'QtCore')

# MAIN WINDOW SETTINGS
WIN_DX = 1400   # width of the main window
WIN_DY = 800    # height of the main window
WIN_X0 = 100    # x-offset of the main window
WIN_Y0 = 100    # y-offset of the main window

ICO_SIZE = 50           # icon size
ICO_DIR = './gui/ico/'   # icon directory

BUTTON_SIZE = 25    # size of buttons

# QCA CELL PARAMETERS
CELL_SEP = 70               # Space between QCA cells
CELL_SIZE = .8*CELL_SEP     # Size of QCA cell
CELL_ALPHA = int(1.*255)   #

# --colors
RED = QtGui.QColor(255, 0, 0)
BLUE = QtGui.QColor(0, 0, 255)
QCA_COL = {'default': QtGui.QColor(255, 255, 255),
           'inactive': QtGui.QColor(100, 100, 100),
           'output': QtGui.QColor(255, 255, 0),
           'input': QtGui.QColor(0, 200, 0),
           'fixed': QtGui.QColor(255, 165, 0),
           'clicked': QtGui.QColor(0, 150, 150)}

DOT_RAD = 0.25*CELL_SIZE

# --qca pen
CELL_PEN_WIDTH = max(1, int(0.05*CELL_SIZE))    #
CELL_PEN_COLOR = QtGui.QColor(180, 180, 180)    #
TEXT_PEN_WIDTH = max(1, int(0.05*CELL_SIZE))    #
TEXT_PEN_COLOR = QtGui.QColor(0, 0, 0)          #
INT_PEN_STYLE = {'strong': QtCore.Qt.SolidLine,
                 'weak': QtCore.Qt.DashLine}
INT_PEN_WIDTH = 3

# --qca magnification
MAX_MAG = 5             # maximum magnification
MIN_MAG = 0.1           # minimum magnification
MAG_STEP = 0.1          # magnification step
MAG_WHEEL_FACT = 0.2    # factor for wheel zoom

QCA_CANVAS_OFFSET = 0.3

# CHIMERA PARAMETERS

CHIMERA_TILE_SIZE = 80
CHIMERA_NODE_RAD = 0.07
CHIMERA_PEN_WIDTH = 0.01*CHIMERA_TILE_SIZE
CHIMERA_NODE_OFFSET = 0.12
CHIMERA_NODE_DELTA = 0.20
CHIMERA_EDGE_WIDTH = 0.02*CHIMERA_TILE_SIZE
CHIMERA_FONT_SIZE = 0.12*CHIMERA_TILE_SIZE
CHIMERA_LABEL_OFFSET = 0.05*CHIMERA_TILE_SIZE
CHIMERA_LABEL_COLOR = QtGui.QColor(255, 0, 0)
CHIMERA_DEFAULT_FILE = '../bin/ch16.txt'

CHIMERA_COL = {'tile': QtGui.QColor(220, 220, 220),
               'tile-selected': QtGui.QColor(170, 170, 170),
               'active': QtGui.QColor(255, 255, 255),
               'inactive': QtGui.QColor(255, 255, 255, 0),
               'used': QtGui.QColor(0, 200, 0),
               'clicked': QtGui.QColor(230, 0, 0),
               'local': QtGui.QColor(0, 0, 200)}
