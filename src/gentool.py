#!/usr/bin/env python
#encoding:  utf-8

from PyQt5.QtWidgets import qApp, QMainWindow, QAction, QApplication,\
        QGraphicsView, QGraphicsScene, QGraphicsEllipseItem, QGraphicsLineItem
from PyQt5.QtGui import QIcon, QBrush, QPen, QColor
from PyQt5.QtCore import Qt, QSize, QPointF

import os

SF = 100

class Node(QGraphicsEllipseItem):

    col = QColor(255,255,255)
    D = 0.3*SF
    z = 2

    def __init__(self, n=0, m=0, parent=None):
        super(Node, self).__init__(parent=parent)

        self.setZValue(self.z)
        self.setRect(0, 0, self.D, self.D)
        self.setBrush(QBrush(self.col))
        self.set_center(n,m)
        
    def set_center(self, n, m):
        self.n, self.m = n, m
        self.setPos(n*SF-.5*self.D, m*SF-.5*self.D)


class SnapDot(Node):

    col = QColor(150,150,0)
    D = 0.2*SF
    z = 3


class Edge(QGraphicsLineItem):

    pen = QPen()
   

class DesignView(QGraphicsView):

    bgcol = QColor(50,50,50)

    def __init__(self, parent=None):
        super(DesignView, self).__init__(parent)
        self.scene().setBackgroundBrush(QBrush(self.bgcol))
        self.setMouseTracking(True)

        self.snapdot = SnapDot()
        self.scene().addItem(self.snapdot)

        self.nodes = {}

        self.cpos = QPointF(0,0)    # position cache
        
    def mouseMoveEvent(self, e):
        pos = self.mapToScene(e.pos())
        delta = (pos - self.cpos).manhattanLength()
        if delta > SF:
            self._nearest_snap(pos)
            self.cpo = pos

    def mouseReleaseEvent(self, e):
        if e.button() == Qt.LeftButton:
            self._add_node()
        if e.button() == Qt.RightButton:
            self._rem_node()

    # internal methods

    def _nearest_snap(self, pos):
        x,y = pos.x(), pos.y()
        n0, m0 = [round(v/SF) for v in [x,y]]
        self.snapdot.set_center(n0,m0)

    def _add_node(self):
        n, m = self.snapdot.n, self.snapdot.m
        if (n,m) not in self.nodes:
            print('Adding Node: ({0},{1})'.format(n,m))
            node = Node(n, m)
            self.scene().addItem(node)
            self.nodes[(n,m)] = node

    def _rem_node(self):
        n, m = self.snapdot.n, self.snapdot.m
        if (n,m) in self.nodes:
            print('Removing Node: ({0},{1})'.format(n,m))
            node = self.nodes.pop((n,m))
            self.scene().removeItem(node)
            
        
        



class GenTool(QMainWindow):

    res_dir = os.path.join('.', 'res')
    
    def __init__(self):
        ''' '''
        super(GenTool, self).__init__()

        self.res = lambda s: os.path.join(self.res_dir,s)

        self.init_gui()
        self.init_actions()
        
    def init_gui(self):

        self.resize(1200,800)
        self.setWindowTitle('Graph Generation Tool')

        self.toolbar = self.addToolBar('Tools')
        self.toolbar.setIconSize(QSize(40,40))

        self.scene = QGraphicsScene()
        self.view = DesignView(self.scene)
        self.setCentralWidget(self.view)

        

    def init_actions(self):

        save_act = QAction(QIcon(self.res('save.png')), 'Save', self)
        save_act.triggered.connect(self.save)

        load_act = QAction(QIcon(self.res('load.png')), 'Load', self)
        load_act.triggered.connect(self.load)

        self.toolbar.addAction(save_act)
        self.toolbar.addAction(load_act)

    # actions

    def save(self):
        pass

    def load(self):
        pass

    def keyPressEvent(self, e):
        if e.key() == Qt.Key_Q:
            qApp.quit()



if __name__ == '__main__':

    import sys
    app = QApplication(sys.argv)

    w = GenTool()
    w.show()

    sys.exit(app.exec_())

    


