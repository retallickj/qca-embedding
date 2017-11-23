#!/usr/bin/env python

# -----------------------------------
# Name: embedder.py
# Desc: Main loop for QCA embedder application
# Author: Jake Retallick
# Created: 2015.11.25
# Modified: 2016.03.02
# Licence: Copyright 2015
# -----------------------------------

import sys

# import Qt based on installed version
from gui.pyqt_import import importPyQt
QtWidgets = importPyQt('QtWidgets')

from gui.mw_embedder import MainWindow



def main():
    '''Main loop which initialises embedder application'''

    app = QtWidgets.QApplication(sys.argv)

    w = MainWindow()
    w.show()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
