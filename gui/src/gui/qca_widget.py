#!/usr/bin/env python

# -----------------------------------
# Name: qca_widget.py
# Desc: QCAWidget class definition
# Author: Jake Retallick
# Created: 2015.11.25
# Modified: 2015.11.25
# Licence: Copyright 2015
# -----------------------------------

from PyQt4 import QtGui, QtCore, QtSvg
import gui_settings as settings
from core.parse_qca import parse_qca_file
from core.auxil import convert_to_full_adjacency, convert_to_lim_adjacency,\
    prepare_convert_adj, CELL_FUNCTIONS, getOutputs

import numpy as np
import matplotlib.pyplot as plt


class QCACellWidget(QtGui.QWidget):
    '''QCA Cell Widget'''

    cm = plt.get_cmap('bwr')

    def __init__(self, parent, cell, spacing=1., offset=[0, 0]):
        ''' '''
        super(QCACellWidget, self).__init__(parent)

        self.parent = parent
        self.qca_widget = parent.parent

        # cell parameters
        self.x = settings.CELL_SEP*(cell['x']-offset[0])*1./spacing
        self.y = settings.CELL_SEP*(cell['y']-offset[1])*1./spacing

        self.cell_pen = QtGui.QPen(settings.CELL_PEN_COLOR)
        self.cell_pen.setWidth(settings.CELL_PEN_WIDTH)

        self.text_pen = QtGui.QPen(settings.TEXT_PEN_COLOR)
        self.text_pen.setWidth(settings.TEXT_PEN_WIDTH)

        self.color = None

        self.type = cell['cf']
        self.qdots = []
        self.num = cell['number']
        for qd in cell['qdots']:
            x = settings.CELL_SEP*(qd['x']-offset[0])*1./spacing
            y = settings.CELL_SEP*(qd['y']-offset[1])*1./spacing
            self.qdots.append([x, y])

        # flags for cell type (simplifies access later)
        self.fixed = self.type == CELL_FUNCTIONS['QCAD_CELL_FIXED']
        self.driver = self.type == CELL_FUNCTIONS['QCAD_CELL_INPUT']
        self.output = self.type == CELL_FUNCTIONS['QCAD_CELL_OUTPUT']
        self.normal = not (self.fixed or self.driver)

        if self.fixed:
            self.pol = cell['pol']
        else:
            self.pol = None

        self.setGeometry(0, 0, settings.CELL_SIZE, settings.CELL_SIZE)
        self.clicked = False
        self.mouse_pos = None

    def set_color(self, col=None):
        self.color = col

    def get_color(self):
        '''Determine the background color of the QCA Cell'''

        if self.color is not None:
            color = self.color
        else:
            # check if selected
            if self.clicked:
                return settings.QCA_COL['clicked']
            if self.type == 0:
                color = settings.QCA_COL['default']
            elif self.type == 1:    # input
                color = settings.QCA_COL['input']
            elif self.type == 2:    # output
                color = settings.QCA_COL['output']
            elif self.type == 3:    # fixed
    #            color = settings.QCA_COL['fixed']
                color = settings.RED if self.pol > 0 else settings.BLUE
            else:
                print('Invalid cell type')
                color = settings.QCA_COL['default']

        if type(color) is tuple:
            color = [int(255*c) for c in color]
            color[3] = settings.CELL_ALPHA
            color = QtGui.QColor(*color)
        return color

    def paintEvent(self, e):
        ''' '''
        painter = QtGui.QPainter()
        painter.begin(self)
        self.drawCell(painter)
        painter.end()

    def drawCell(self, painter):
        '''Redraw cell widget'''

        painter.setPen(self.cell_pen)
        painter.setBrush(self.get_color())

        painter.drawRect(self.geometry())

        # write cell label
        painter.setPen(self.text_pen)
        painter.setFont(QtGui.QFont('Decorative', 10))
        painter.drawText(self.geometry(), QtCore.Qt.AlignCenter, str(self.num))

    def mousePressEvent(self, e):
        ''' '''
        self.mouse_pos = e.pos()
        self.parent.mousePressEvent(e)

    def mouseReleaseEvent(self, e):
        ''' '''
        if self.mouse_pos is not None:
            diff = e.pos()-self.mouse_pos
            if max(abs(diff.x()), abs(diff.y())) < settings.CELL_SIZE:
                self.qca_widget.onClick(self.num)
        self.mouse_pos = None
        self.parent.mouseReleaseEvent(e)

class Canvas(QtGui.QWidget):
    '''Canvas to draw QCA cells'''

    def __init__(self, parent):
        ''' '''
        super(Canvas, self).__init__()
        self.parent = parent
        self.scaling = 1.

        self.w = 1.
        self.h = 1.

    def rescale(self, zoom=True, f=1.):
        ''' '''
        step = f*settings.MAG_STEP
        if zoom:
            self.scaling = min(settings.MAX_MAG, self.scaling + step)
        else:
            self.scaling = max(settings.MIN_MAG, self.scaling - step)
        geo = self.geometry()
        geo.setWidth(self.w*self.scaling)
        geo.setHeight(self.h*self.scaling)
        self.setGeometry(geo)
        self.update()

    def moveCells(self):
        ''' '''
        for cell in self.parent.cells:
            _x = cell.x*self.scaling
            _y = cell.y*self.scaling
            _size = settings.CELL_SIZE*self.scaling
            cell.setGeometry(_x, _y, _size, _size)

    def drawCells(self, painter):
        ''' '''
        for cell in self.parent.cells:
            cell.drawCell(painter)

    def drawConnections(self, painter):
        ''' '''
        J0 = self.parent.J0

        # store cell centers
        X = []
        Y = []
        for cell in self.parent.cells:
            geo = cell.geometry()
            X.append(geo.x()+.5*geo.width())
            Y.append(geo.y()+.5*geo.height())
        # draw all non-zero interactions
        pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 255))
        pen.setWidth(max(1, settings.INT_PEN_WIDTH*self.scaling))
        for i in xrange(J0.shape[0]-1):
            for j in xrange(i+1, J0.shape[0]):
                if abs(J0[i, j]) > 0.:
                    pen.setStyle(settings.INT_PEN_STYLE[
                        'strong' if abs(J0[i, j]) > 0.5 else 'weak'])
                    painter.setPen(pen)
                    painter.drawLine(X[i], Y[i], X[j], Y[j])

    def paint(self, painter):
        ''' '''
        self.moveCells()
        self.drawConnections(painter)
        self.drawCells(painter)

    # Interrupts

    def paintEvent(self, e):
        ''' '''
        painter = QtGui.QPainter()
        painter.begin(self)
        self.paint(painter)
        painter.end()

    def mouseDoubleClickEvent(self, e):
        ''' '''
        pass
        # determine which cell was clicked


class QCAWidget(QtGui.QScrollArea):
    '''Widget for viewing QCA circuits'''

    def __init__(self, parent=None, filename=None):
        ''' '''
        super(QCAWidget, self).__init__(parent)

        self.parent = parent

        # parameters
        self.cells = []         # list of qca cells
        self.spacing = 1.       # cell-cell spacing value
        self.J = np.zeros(0)    # cell-cell interaction matrix
        self.J0 = self.J        # J with only included interactions

        # mouse tracking
        self.mouse_pos = None
        self.zoom_flag = False

        # initialise UI
        self.initUI()

        if filename is not None:
            self.updateCircuit(filename)

    def initUI(self):
        ''' '''

        # create main widget
        self.canvas = Canvas(self)
        self.canvas.setGeometry(0, 0, 0, 0)
        self.setWidget(self.canvas)

    def updateCircuit(self, filename, full_adj):
        ''' '''

        try:
            cells, spacing, J = parse_qca_file(filename)
        except Exception as e:
            print(e.message)
            print('Failed to load QCA File...')
            return

        # save filename for reference
        self.filename = filename

        # forget old circuit
        for cell in self.cells:
            cell.setParent(None)
        self.cells = []

        # update J coefficients
        self.J = J

        # set up adjacency conversion variables
        Js, T, A, DX, DY = prepare_convert_adj(cells, spacing, J)
        self.convert_vars = {'Js': Js,
                             'T': T,
                             'A': A,
                             'DX': DX,
                             'DY': DY}

        self.outputs = getOutputs(cells)

        # set adjacency type
        self.setAdjacency(full_adj, update=False)  # sets full_adj and J0

        # find span and offset of circuit: currently inefficient
        x_min = min([cell['x'] for cell in cells])
        x_max = max([cell['x'] for cell in cells])
        y_min = min([cell['y'] for cell in cells])
        y_max = max([cell['y'] for cell in cells])

        o = settings.QCA_CANVAS_OFFSET
        span = [x_max-x_min, y_max-y_min]
        offset = [x_min-o*span[0], y_min-o*span[1]]

        # update circuit constants
        self.spacing = spacing
        self.offset = offset

        # update size and scaling of canvas
        factor = (1+2*o)*settings.CELL_SEP*1.1/spacing
        self.canvas.scaling = 1.
        self.canvas.w = (1+span[0])*factor
        self.canvas.h = (1+span[1])*factor
        self.canvas.setGeometry(0, 0,
                                self.canvas.w, self.canvas.h)

        # add new cells
        for cell in cells:
            self.addCell(cell)

        self.canvas.update()

    def setAdjacency(self, full_adj, update=True):
        ''' '''
        self.full_adj = full_adj
        if full_adj:
            self.J0 = convert_to_full_adjacency(self.J, **self.convert_vars)
        else:
            self.J0 = convert_to_lim_adjacency(self.J, **self.convert_vars)

        self.J0 /= np.max(np.abs(self.J0))

        if update:
            self.canvas.update()

    def setPols(self, pols):
        '''set the cell colours based on the given cell polarizations'''

        scale = 0.7

        if pols is None:
            for cell in self.cells:
                cell.setColor(None)
        else:
            for i, pol in pols.items():
                cell = self.cells[i]
                if pol is None:
                    cell.set_color(None)
                else:
                    p = .5*(1+scale*pol)
                    cell.set_color(cell.cm(p))

        self.canvas.update()


    def addCell(self, cell):
        ''' '''
        cell = QCACellWidget(self.canvas, cell,
                             spacing=self.spacing, offset=self.offset)
        self.cells.append(cell)
        cell.show()

    def onClick(self, cell):
        '''Response to clicking on one the QCA cells'''
        if self.cells[cell].driver or self.cells[cell].fixed:
            return

        if self.parent.active_embedding != -1:
            self.selectCell(cell)
            # make changes in chimera widget
            embedding = self.parent.embeddings[self.parent.active_embedding]
            self.parent.chimera_widget.selectNodes(embedding, cell)

    def selectCell(self, num):
        ''' '''
        cell = self.cells[num]
        if cell.driver or cell.fixed or cell.clicked:
            return

        for c2 in self.cells:
            if c2.clicked:
                c2.clicked = False
                self.canvas.update(c2.geometry())

        cell.clicked = True
        self.canvas.update(cell.geometry())

    def prepareCircuit(self):
        '''Return needed parameters for embedding'''

        return self.J0, self.cells

    def save_svg(self, fname):
        '''Write the QCA circuit to an svg file'''

        generator = QtSvg.QSvgGenerator()
        generator.setFileName(fname)
        generator.setSize(self.canvas.size())
        generator.setViewBox(self.canvas.rect())

        painter = QtGui.QPainter()
        painter.begin(generator)
        self.canvas.paint(painter)
        painter.end()

    # interrupts

    def mousePressEvent(self, e):
        '''On left click drag circuit, on right click highlight cell'''
        self.mouse_pos = e.pos()

    def mouseMoveEvent(self, e):
        ''' '''
        if self.mouse_pos is not None:
            diff = e.pos()-self.mouse_pos
            self.mouse_pos = e.pos()
            self.verticalScrollBar().setValue(
                self.verticalScrollBar().value()-diff.y())
            self.horizontalScrollBar().setValue(
                self.horizontalScrollBar().value()-diff.x())

    def mouseReleaseEvent(self, e):
        '''On mouse release, forget old mouse position to avoid
        jumping'''
        self.mouse_pos = None

    def keyPressEvent(self, e):
        ''' '''
        if e.key() == QtCore.Qt.Key_Minus:
            self.canvas.rescale(zoom=False)
        elif e.key() == QtCore.Qt.Key_Plus:
            self.canvas.rescale(zoom=True)
        elif e.key() == QtCore.Qt.Key_Control:
            self.zoom_flag = True
        self.parent.keyPressEvent(e)

    def keyReleaseEvent(self, e):
        ''' '''
        if e.key() == QtCore.Qt.Key_Control:
            self.zoom_flag = False

    def wheelEvent(self, e):
        '''Scrolling options'''

        if self.zoom_flag:
            self.canvas.rescale(zoom=e.delta() > 0, f=settings.MAG_WHEEL_FACT)
        else:
            super(QCAWidget, self).wheelEvent(e)
