#!/usr/bin/env python

from PyQt4 import QtGui
from core.parse_qca import parse_qca_file
from core.auxil import qca_to_coef, hash_problem

import sys
from pprint import pprint
from random import shuffle
import numpy as np

#
#class QCACell(QtGui.QWidget):
#    ''' '''
#    
#    def __init__(self, parent=None):
#        super(QCACell, self).__init__(parent)
#
#class Canvas(QtGui.QWidget):
#    ''' '''
#
#    def __init__(self, parent=None):
#        super(Canvas, self).__init__(parent)
#
#
#    def paintEvent(self, e):
#
#class MainWindow(QtGui.QMainWindow):
#    ''' '''
#
#    def __init__(self):
#        ''' '''
#        super(MainWindow, self).__init__()
#        self.initUI()
#
#    def initUI(self):
#
#        self.scrollarea = QtGui.QScrollArea()
#        self.setCentralWidget(self.scrollarea)
#
#        self.canvas = Canvas(self)



def mix_seriation(h, J):
    '''Apply random invariants to h and J'''

    # shuffle
    inds = range(len(h))
    shuffle(inds)

    K = np.random.rand()
    hs = 1-2*(np.random.rand()<.5)
    h_ = h[inds]*K*hs
    J_ = J[inds, :][:, inds]*K

    return h_, J_, K, hs, inds

def main(fname):
    ''' '''

    try:
        cells, spacing, zones, J, fb = parse_qca_file(fname, one_zone=True)
    except:
        print('Failed to load qca file')
        return

    h, J = qca_to_coef(cells, spacing, J, adj='full')

    h /= np.max(np.abs(J))
    J /= np.max(np.abs(J))

    for _ in range(100):
        h_, J_, K, hp, inds = mix_seriation(h, J)
        hval, K_, hp_, inds_ = hash_problem(h_, J_)

        if False:
            print('K: {0:.4f}\t hp: {1}\ninds: {2}'.format(K, hp, inds))
            print('K: {0:.4f}\t hp: {1}\ninds: {2}'.format(K_, hp_, inds_))

        print('hash val: {0}'.format(hval))

#    app = QtGui.QApplication(sys.argv)
#
#    w = MainWindow()
#    w.show()
#
#    sys.exit(app.exec_())

if __name__ == '__main__':

    try:
        fname = sys.argv[1]
    except:
        print('No QCAD file given...')
        sys.exit()

    main(fname)
