#!/usr/bin/env python

# -----------------------------------
# Name: embedder.py
# Desc: Main loop for QCA embedder application
# Author: Jake Retallick
# Created: 2015.11.25
# Modified: 2016.03.02
# Licence: Copyright 2015
# -----------------------------------

from PyQt4 import QtGui
from gui.mw_embedder import MainWindow
import sys


def main():
    '''Main loop which initialises embedder application'''

    app = QtGui.QApplication(sys.argv)

    w = MainWindow()
    w.show()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
