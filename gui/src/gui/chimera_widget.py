#!/usr/bin/env python

# -----------------------------------
# Name: chimera_widget.py
# Desc: ChimeraWidget class definition
# Author: Jake Retallick
# Created: 2015.11.25
# Modified: 2015.11.25
# Licence: Copyright 2015
# -----------------------------------

import sys

# import Qt based on installed version
from gui.pyqt_import import importPyQt
QtGui, QtWidgets, QtCore, QtSvg = importPyQt('QtGui', 'QtWidgets', 'QtCore', 'QtSvg')

import gui.gui_settings as settings
from core.chimera import load_chimera_file, linear_to_tuple


class ChimeraNode(QtWidgets.QWidget):
    '''Graph node representing a qubit on D-Wave's processor'''

    def __init__(self, tile, h, l, active):
        '''inputs:  tile    - ChimeraTile object containing the node
                    h       - flag for horizontal qubit
                    l       - index of qubit in tile (0...3)
        '''

        super(ChimeraNode, self).__init__(tile)

        self.tile = tile        # Tile containing the node

        self.cell = None            # cell with model containing the node
        self.embedding_ind = None   # embedding index
        self.used = False           # flag if node is taken by an embedding
        self.clicked = False
        self.local = False

        self.active = active    # flag determines whether node is active

        self.h = h
        self.l = l

        self.mouse_pos = None

        # determine center position within tile
        if l < 2:
            self.z = settings.CHIMERA_NODE_OFFSET + \
                settings.CHIMERA_NODE_DELTA*l
        else:
            self.z = 1-settings.CHIMERA_NODE_OFFSET - \
                settings.CHIMERA_NODE_DELTA*(3-l)

    def getColor(self):
        '''Return node color'''

        if not self.active:
            color = settings.CHIMERA_COL['inactive']
        if self.local:
            if self.clicked:
                color = settings.CHIMERA_COL['clicked']
            else:
                color = settings.CHIMERA_COL['local']
        elif self.used:
            color = settings.CHIMERA_COL['used']
        else:
            color = settings.CHIMERA_COL['active']
        return color

    def drawNode(self, painter):
        '''Draw the node within the tile'''

        # determine center point of node
        if self.h:
            x = .5*self.tile.width()
            y = self.z*self.tile.height()
        else:
            y = .5*self.tile.height()
            x = self.z*self.tile.width()

        # shift to top left corner
        x = x-settings.CHIMERA_NODE_RAD*self.tile.width()
        y = y-settings.CHIMERA_NODE_RAD*self.tile.width()

        size = 2*settings.CHIMERA_NODE_RAD*self.tile.width()

        self.setGeometry(x, y, size, size)
        rect = QtCore.QRect(x+self.tile.x(), y+self.tile.y(),
                            size, size)
        pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 255))
        pen.setWidth(settings.CHIMERA_PEN_WIDTH)
        brush = QtGui.QBrush(self.getColor())

        painter.setPen(pen)
        painter.setBrush(brush)

        painter.drawEllipse(rect)

    def setUsed(self, ind, cell):
        ''' '''
        self.used = True
        self.embedding_ind = ind
        self.cell = cell

    def mousePressEvent(self, e):
        ''' '''
        self.mouse_pos = e.pos()
        self.tile.parent.mousePressEvent(e)

    def mouseReleaseEvent(self, e):
        if self.mouse_pos is not None:
            diff = e.pos() - self.mouse_pos
            if max(abs(diff.x()), abs(diff.y())) < self.width():
                self.tile.parent.onNodeClick(self)
        self.mouse_pos = None


class ChimeraTile(QtWidgets.QWidget):
    ''' '''

    def __init__(self, parent, m, n, adj=None):
        '''Tile in the chimera graph'''

        super(ChimeraTile, self).__init__(parent)

        self.parent = parent    # ChimeraWidget
        self.m = m              # tile row
        self.n = n              # tile column
        self.adj = adj          # full circuit adjacency dict

        self.mouse_pos = None
        self.nodes = {}
        self.selected = False

        # initialise nodes
        for h in [True, False]:
            for l in range(4):
                key = (self.m, self.n, h, l)
                active = len(adj[key]) > 0
                if active:
                    node = ChimeraNode(self, h, l, active)
                    self.nodes[(h, l)] = node
                    node.show()

    def mousePressEvent(self, e):
        '''On mouse click, store mouse position to check release'''
        self.mouse_pos = e.pos()
        self.parent.mousePressEvent(e)

    def mouseReleaseEvent(self, e):
        '''On mouse release, if same tile, select tile for embedding'''
        if e.button() == QtCore.Qt.RightButton:
            self.parent.releaseSelection()
        elif self.mouse_pos is not None:
            diff = e.pos()-self.mouse_pos
            if max(abs(diff.x()), abs(diff.y())) < self.width():
                # same tile release, select tile pass to
                self.parent.onTileClick(self.m, self.n)
        self.parent.mouseReleaseEvent(e)

    def getColor(self):
        '''Get the color of the tile'''
        if self.selected:
            return settings.CHIMERA_COL['tile-selected']
        else:
            return settings.CHIMERA_COL['tile']

    def drawTile(self, painter):
        '''Draw tile background'''

        # outline
        pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 255))
        brush = QtGui.QBrush(self.getColor())

        painter.setPen(pen)
        painter.setBrush(brush)

        painter.drawRect(self.geometry())

    def drawLabel(self, painter):

        pen = QtGui.QPen(settings.CHIMERA_LABEL_COLOR)
        painter.setPen(pen)
        painter.setFont(QtGui.QFont('Decorative',
                                    settings.CHIMERA_FONT_SIZE))
        painter.drawText(
            self.x()+settings.CHIMERA_LABEL_OFFSET,
            self.y()+self.height()-settings.CHIMERA_LABEL_OFFSET,
            '{0}:{1}'.format(self.m, self.n))

    def drawNodes(self, painter):
        '''Draw nodes within the tile'''

        for key in self.nodes:
            node = self.nodes[key]
            node.drawNode(painter)


class Canvas(QtWidgets.QWidget):
    ''' '''

    def __init__(self, parent):
        '''Initialise Canvas'''
        super(Canvas, self).__init__(parent)
        self.parent = parent
        self.scaling = 1.

    def drawTiles(self, painter):
        '''Draw nodes and connectors within tiles'''

        for key in self.parent.tiles:
            m, n = key
            tile = self.parent.tiles[key]
            tile.drawTile(painter)

    def drawConnectors(self, painter):  # UGLY
        '''Draw connectors between tiles'''

        pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 255))
        pen.setWidth(settings.CHIMERA_EDGE_WIDTH)
        painter.setPen(pen)

        # for each tile, draw internal edges and all edges to higher
        # indexed tiles
        for m1, n1 in self.parent.tiles:
            tile1 = self.parent.tiles[(m1, n1)]
            dxy = settings.CHIMERA_NODE_RAD*tile1.width()
            for h1, l1 in tile1.nodes:
                node1 = tile1.nodes[(h1, l1)]
                x1 = tile1.x() + node1.x() + dxy
                y1 = tile1.y() + node1.y() + dxy
                for m2, n2, h2, l2 in self.parent.adj[(m1, n1, h1, l1)]:
                    # only draw internal connections once
                    if m1 == m2 and n1 == n2 and h2:
                        node2 = tile1.nodes[(h2, l2)]
                        x2 = tile1.x() + node2.x() + dxy
                        y2 = tile1.y() + node2.y() + dxy
                        painter.drawLine(x1, y1, x2, y2)
                    # only draw external connections once
                    elif m2 > m1 or n2 > n1:
                        tile2 = self.parent.tiles[(m2, n2)]
                        node2 = tile2.nodes[(h2, l2)]
                        x2 = tile2.x() + node2.x() + dxy
                        y2 = tile2.y() + node2.y() + dxy
                        painter.drawLine(x1, y1, x2, y2)

    def drawNodes(self, painter):
        ''' '''

        for key in self.parent.tiles:
            m, n = key
            tile = self.parent.tiles[key]
            tile.drawNodes(painter)

    def drawLabels(self, painter):
        ''' '''

        for key in self.parent.tiles:
            m, n = key
            tile = self.parent.tiles[key]
            tile.drawLabel(painter)

    def paint(self, painter):
        ''' '''

        # draw tile background
        self.drawTiles(painter)

        # draw connectors
        self.drawConnectors(painter)

        # draw nodes
        self.drawNodes(painter)

        # draw tile labels
        self.drawLabels(painter)

    def paintEvent(self, e):
        ''' '''

        painter = QtGui.QPainter()
        painter.begin(self)
        self.paint(painter)
        painter.end()


class ChimeraWidget(QtWidgets.QScrollArea):
    '''Widget for viewing QCA circuits'''

    def __init__(self, parent=None):
        ''' '''
        super(ChimeraWidget, self).__init__(parent)

        self.parent = parent

        # parameters
        self.shift = False          # flag for shift pressed
        self.clicked_tile = None    # pair of indices for selected tile
        self.active_range = None    # range of tiles in active graph
        self.tiles = {}
        self.adj = {}

        self.mouse_pos = None

        self.initUI()

#        self.updateChimera(settings.CHIMERA_DEFAULT_FILE)

    def initUI(self):
        '''Initialise UI'''

        self.layout = QtWidgets.QGridLayout(self)
        self.canvas = Canvas(self)
        self.setWidget(self.canvas)
        self.canvas.setLayout(self.layout)

    def updateChimera(self, filename):
        '''Process a chimera specification file and update the widget'''

        try:
            M, N, adj = load_chimera_file(filename)
        except IOError:
            print('Failed to load given file...')
            return

        # forget old grid layout
        for tile in self.tiles:
            self.tiles[tile].setParent(None)
        self.tiles = {}

        while self.layout.count():
            item = self.layout.takeAt(0)
            item.widget().deleteLater()

        # resize canvas
        width = N*settings.CHIMERA_TILE_SIZE
        height = M*settings.CHIMERA_TILE_SIZE

        self.canvas.setGeometry(0, 0, width, height)

        # convert adjacency dict to tuple format
        adj = {linear_to_tuple(k, M, N):\
            [linear_to_tuple(a, M, N) for a in adj[k]] for k in adj}

        self.M = M
        self.N = N
        self.adj = adj

        for m in range(M):
            for n in range(N):
                tile = ChimeraTile(self, m, n, adj=adj)
                self.tiles[(m, n)] = tile
                self.layout.addWidget(tile, m, n)
                tile.show()

        self.canvas.update()

    def onTileClick(self, m, n):
        '''If a tile is clicked and shift-flag False, start subgraph
        select. If a tile is clicked and shift-flag True, end subgraph
        select.'''

        if self.clicked_tile is not None and self.shift:
            # unselected all other tiles
            for key in self.tiles:
                self.tiles[key].selected = False
            # select subgraph
            if self.clicked_tile != (m, n):
                self.setActiveRange(self.clicked_tile, (m, n))
                self.canvas.update()
                self.click_tile = None
        else:
            # first corner of subgraph
            self.clicked_tile = (m, n)
            # unselected all other tiles
            for key in self.tiles:
                self.tiles[key].selected = False
            # update selected subgraph
            self.setActiveRange((m, n), None)
            self.canvas.update()

    def unclickNodes(self):
        '''Reset the clicked and local flags of all nodes'''
        for m, n in self.tiles:
            tile = self.tiles[(m, n)]
            redraw = False
            for h, l in tile.nodes:
                node = tile.nodes[(h, l)]
                if node.local:
                    node.local = False
                    node.clicked = False
                    redraw = True
            if redraw:
                self.canvas.update(tile.geometry())

    def resetNodes(self, embedding):
        '''Reset the nodes of an embedding to non-used'''

        models = embedding.models
        active_range = embedding.active_range

        tiles = set()
        for c1 in models:
            for m, n, h, l in models[c1]:
                # correct tile for active range
                m += active_range['M'][0]
                n += active_range['N'][0]
                # set node as clicked
                tiles.add((m, n))
                self.tiles[(m, n)].nodes[(h, l)].used=False
                self.tiles[(m, n)].nodes[(h, l)].clicked=False
                self.tiles[(m, n)].nodes[(h, l)].local=False
        for m, n in tiles:
            self.canvas.update(self.tiles[(m, n)].geometry())

    def onNodeClick(self, node):
        '''On node click'''

        if not node.used:
            return

        cell = node.cell
        embedding = self.parent.embeddings[node.embedding_ind]

        # make changes to QCA widget
        if self.parent.active_embedding != node.embedding_ind:
            self.parent.switchEmbedding(node.embedding_ind, color=False)

        # highlight corresponding cell in QCAWidget
        self.parent.qca_widget.selectCell(cell)

        # highlight nodes of embedding
        self.selectNodes(embedding, cell)

    def selectNodes(self, embedding, cell):
        '''Select the nodes corresponding to a cell in an embedding'''

        active_range = embedding.active_range

        # unset all other nodes
        self.unclickNodes()

        # set all nodes of the circuit as local
        tiles = set()
        for c1 in embedding.models:
            for m, n, h, l in embedding.models[c1]:
                # correct tile for active range
                m += active_range['M'][0]
                n += active_range['N'][0]
                # set node as clicked
                tiles.add((m, n))
                self.tiles[(m, n)].nodes[(h, l)].local=True
                if c1 == cell:
                    self.tiles[(m, n)].nodes[(h, l)].clicked=True

        for m, n in tiles:
            self.canvas.update(self.tiles[(m, n)].geometry())

    def setActiveRange(self, tile1, tile2):
        '''Set the active range of the chimera graph'''

        # find range of selected tiles
        if tile2 is None:
            ry = [tile1[0], tile1[0]+1]
            rx = [tile1[1], tile1[1]+1]
        else:
            ry = [min(tile1[0], tile2[0]), max(tile1[0], tile2[0])+1]
            rx = [min(tile1[1], tile2[1]), max(tile1[1], tile2[1])+1]

        # for now just select tiles
        for m in range(ry[0], ry[1]):
            for n in range(rx[0], rx[1]):
                self.tiles[(m, n)].selected = True

        self.active_range = {'M': ry,
                             'N': rx}

    def releaseSelection(self):
        ''' '''

        for key in self.tiles:
            self.tiles[key].selected = False

        self.clicked_tile = None
        self.active_range = None

        self.canvas.update()

    def getActiveGraph(self):
        '''Return the adjacency matrix of active range of the chimera graph'''

        if self.active_range is None:
            return self.M, self.N, self.adj, {'M': [0, self.M], 'N': [0, self.N]}

        M = self.active_range['M'][1] - self.active_range['M'][0]
        N = self.active_range['N'][1] - self.active_range['N'][0]

        # define check functions
        tile_check = lambda m, n, h, l: \
            m >= self.active_range['M'][0] and\
            m < self.active_range['M'][1] and\
            n >= self.active_range['N'][0] and\
            n < self.active_range['N'][1]

        node_check = lambda m, n, h, l: tile_check(m, n, h, l) and\
            len(self.adj[(m, n, h, l)]) > 0 and\
            not self.tiles[(m, n)].nodes[(h, l)].used

        # include only nodes within active range
        adj = {k1: [] for k1 in self.adj if tile_check(*k1)}

        for k1 in adj:
            if node_check(*k1):
                adj[k1] = [k2 for k2 in self.adj[k1] if node_check(*k2)]

        # offset adjacency list to zero
        dm = self.active_range['M'][0]
        dn = self.active_range['N'][0]

        offset = lambda m, n, h, l: (m-dm, n-dn, h, l)

        adj = {offset(*k1): [offset(*k2) for k2 in adj[k1]] for k1 in adj}

        return M, N, adj, self.active_range

    # ADD EMBEDDING

    def addEmbedding(self, embedding, ind):
        '''Add an embedding with given index'''

        models = embedding.models
        active_range = embedding.active_range

        tiles = set()
        for cell in models:
            for m, n, h, l in models[cell]:
                m += active_range['M'][0]
                n += active_range['N'][0]
                self.tiles[(m, n)].nodes[(h, l)].setUsed(ind, cell)

        for tile in tiles:
            self.canvas.update(self.tiles[tile].geometry())

    def save_svg(self, fname):
        '''Write the cchimera graph to an svg file'''

        generator = QtSvg.QSvgGenerator()
        generator.setFileName(fname)
        generator.setSize(self.canvas.size())
        generator.setViewBox(self.canvas.rect())

        painter = QtGui.QPainter()
        painter.begin(generator)
        self.canvas.paint(painter)
        painter.end()

    # EVENT HANDLERS

    def keyPressEvent(self, e):
        ''' '''
        scroll_delta = settings.CHIMERA_TILE_SIZE
        if e.key() == QtCore.Qt.Key_Shift:
            self.shift = True
        elif e.key() == QtCore.Qt.Key_Left:
            self.horizontalScrollBar().setValue(
                self.horizontalScrollBar().value() - scroll_delta)
        elif e.key() == QtCore.Qt.Key_Right:
            self.horizontalScrollBar().setValue(
                self.horizontalScrollBar().value() + scroll_delta)
        elif e.key() == QtCore.Qt.Key_Up:
            self.verticalScrollBar().setValue(
                self.verticalScrollBar().value() - scroll_delta)
        elif e.key() == QtCore.Qt.Key_Down:
            self.verticalScrollBar().setValue(
                self.verticalScrollBar().value() + scroll_delta)

        self.parent.keyPressEvent(e)

    def keyReleaseEvent(self, e):
        '''Reset key flags'''

        self.shift = False

    def mousePressEvent(self, e):
        '''On left click drag circuit, on right click highlight cell'''
        self.mouse_pos = e.pos()

    def mouseMoveEvent(self, e):
        ''' '''
        return  # disable for now.... cause problems
        if self.mouse_pos is not None:
            diff = e.pos()-self.mouse_pos
            self.mouse_pos = e.pos()
            self.verticalScrollBar().setValue(
                self.verticalScrollBar().value()-diff.y())
            self.horizontalScrollBar().setValue(
                self.horizontalScrollBar().value()-diff.x())

    def mouseReleaseEvent(self, e):
        '''On mouse release, forget old mouse position to avoid
        jumping. On right click release unselect everything'''
        self.mouse_pos = None
