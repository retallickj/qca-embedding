#!/usr/bin/env python

# -----------------------------------
# Name: Main Window widget for embedder application
# Author: Jake Retallick
# Created: 2015.11.25
# Modified: 2015.11.25
# Licence: Copyright 2015
# -----------------------------------

import sys

# import Qt based on installed version
from gui.pyqt_import import importPyQt
QtGui, QtWidgets, QtCore = importPyQt('QtGui', 'QtWidgets', 'QtCore')

import os

import gui.gui_settings as settings
from gui.qca_widget import QCAWidget
from gui.chimera_widget import ChimeraWidget
from core.classes import Embedding, get_embedder_flags
from core.chimera import tuple_to_linear
import traceback

class MainWindow(QtWidgets.QMainWindow):
    '''Main Window widget for embedder application'''

    def __init__(self):
        '''Create the Main Window widget'''

        super(MainWindow, self).__init__()
        self.initUI()

    def initUI(self):
        '''Initialise the UI'''

        # default parameters
        self.qca_dir = os.getcwd()
        self.embed_dir = os.getcwd()
        self.chimera_dir = os.getcwd()
        self.coef_dir = os.getcwd()
        self.svg_dir = os.getcwd()

        # functionality parameters

        self.chimera_file = ''      # relative path to chimera file
        self.qca_active = False     # True when QCAWidget set
        self.full_adj = True        # True when using full adjacency
        self.embed_method = 'dense' # Value of embedding method (Default: Dense)
        self.tile_style = 0         # tile style

        self.embeddings = {}        # list of embeddings
        self.active_embedding = -1  # index of active embedding
        self.embedding_count = 0    # next embedding index
        self.embedding_actions = {}
        self.embedding_menus = {}

        self.coupling_strength = 1. # relative strength of internal couplers

        # main window parameters
        geo = [settings.WIN_X0, settings.WIN_Y0,
               settings.WIN_DX, settings.WIN_DY]
        self.setGeometry(*geo)
        self.setWindowTitle('QCA Embedder')

        self.statusBar()

        # build the menu
        self.init_menubar()

        # build the toolbar
        self.init_toolbar()

        # set up the main layout
        hbox = QtWidgets.QHBoxLayout()

        # QCA widget placeholder
        self.qca_widget = QCAWidget(self)

        # Chimera widget
        self.chimera_widget = ChimeraWidget(self)
        self.chimera_file = os.path.relpath(settings.CHIMERA_DEFAULT_FILE)
        self.chimera_widget.updateChimera(self.chimera_file)
        self.action_save_chimera_svg.setEnabled(True)

        hbox.addWidget(self.qca_widget, stretch=4)
        hbox.addWidget(self.chimera_widget, stretch=4)

        main_widget = QtWidgets.QWidget(self)
        main_widget.setLayout(hbox)
        self.setCentralWidget(main_widget)

    def init_menubar(self):
        ''' '''

        menubar = self.menuBar()

        file_menu = menubar.addMenu('&File')
        tool_menu = menubar.addMenu('&Tools')
        menubar.addSeparator()
        self.embeddings_menu = menubar.addMenu('&Embeddings')
        self.embeddings_menu.setEnabled(False)

#        view_menu = menubar.addMenu('&View')
#        help_menu = menubar.addMenu('&Help')

        # construct actions

        # Loading methods

        qcaFileAction = QtWidgets.QAction(
            QtGui.QIcon(settings.ICO_DIR+'qca_file.png'),
            'Open QCA file...', self)
        qcaFileAction.triggered.connect(self.load_qca_file)

        embedFileAction = QtWidgets.QAction(
            QtGui.QIcon(settings.ICO_DIR+'open_embed.png'),
            'Open EMBED file...', self)
        embedFileAction.triggered.connect(self.load_embed_file)

        chimeraFileAction = QtWidgets.QAction(
            QtGui.QIcon(settings.ICO_DIR+'chimera_file.png'),
            'Open chimera file...', self)
        chimeraFileAction.triggered.connect(self.load_chimera_file)

        # Saving methods

        self.action_save_embedding = QtWidgets.QAction('Save active embedding...', self)
        self.action_save_embedding.triggered.connect(self.save_active_embedding)
        self.action_save_embedding.setEnabled(False)

        self.action_save_all = QtWidgets.QAction('Save EMBED file...', self)
        self.action_save_all.triggered.connect(self.save_all_embeddings)
        self.action_save_all.setEnabled(False)

        self.action_export_coefs = QtWidgets.QAction('Export coef file...', self)
        self.action_export_coefs.setIcon(
            QtGui.QIcon(settings.ICO_DIR+'upload.png'))
        self.action_export_coefs.setStatusTip('Create coefficient files...')
        self.action_export_coefs.triggered.connect(self.export_coefs)
        self.action_export_coefs.setEnabled(False)

        # SVG exporting

        self.action_save_qca_svg = QtWidgets.QAction('Save qca widget as SVG...', self)
        self.action_save_qca_svg.triggered.connect(self.save_qca_svg)
        self.action_save_qca_svg.setEnabled(False)

        self.action_save_chimera_svg = QtWidgets.QAction('Save chimera widget as SVG...', self)
        self.action_save_chimera_svg.triggered.connect(self.save_chimera_svg)
        self.action_save_chimera_svg.setEnabled(False)

        # exit

        exitAction = QtWidgets.QAction('Exit', self)
        exitAction.setShortcut('Ctrl+W')
        exitAction.triggered.connect(self.close)

        # Tool menu

        self.action_dense_embed_flag = QtWidgets.QAction('Dense', self)
        self.action_dense_embed_flag.triggered.connect(self.switch_dense_embed)

        self.action_layout_embed_flag = QtWidgets.QAction('Layout-Aware', self)
        self.action_layout_embed_flag.triggered.connect(self.switch_layout_embed)

        self.action_heur_embed_flag = QtWidgets.QAction('Heuristic', self)
        self.action_heur_embed_flag.triggered.connect(self.switch_heur_embed)

        tile_func_ab = lambda: self.set_tile_style(0)
        tile_func_a = lambda: self.set_tile_style(-1)
        tile_func_b = lambda: self.set_tile_style(1)

        self.action_tile_AB_flag = QtWidgets.QAction('AB', self)
        self.action_tile_AB_flag.triggered.connect(tile_func_ab)
        self.action_tile_AB_flag.setEnabled(False)

        self.action_tile_A_flag = QtWidgets.QAction('A', self)
        self.action_tile_A_flag.triggered.connect(tile_func_a)
        self.action_tile_A_flag.setEnabled(True)

        self.action_tile_B_flag = QtWidgets.QAction('B', self)
        self.action_tile_B_flag.triggered.connect(tile_func_b)
        self.action_tile_B_flag.setEnabled(True)

        self.action_set_coupling = QtWidgets.QAction('Coupling Strength...', self)
        self.action_set_coupling.triggered.connect(self.set_coupling)
        self.action_set_coupling.setEnabled(True)

        file_menu.addAction(qcaFileAction)
        file_menu.addAction(embedFileAction)
        file_menu.addAction(chimeraFileAction)
        file_menu.addSeparator()
        file_menu.addAction(self.action_save_embedding)
        file_menu.addAction(self.action_save_all)
        file_menu.addAction(self.action_export_coefs)
        file_menu.addSeparator()
        file_menu.addAction(self.action_save_qca_svg)
        file_menu.addAction(self.action_save_chimera_svg)
        file_menu.addSeparator()
        file_menu.addAction(exitAction)

        embedder_menu = tool_menu.addMenu('Embedding method')

        embedders = get_embedder_flags()

        self.embed_method = 'dense'
        self.action_dense_embed_flag.setEnabled(False)

        if embedders['dense']:
            embedder_menu.addAction(self.action_dense_embed_flag)

        if embedders['layout']:
            embedder_menu.addAction(self.action_layout_embed_flag)

        if embedders['heur']:
            embedder_menu.addAction(self.action_heur_embed_flag)

        print('Using embedder: {0}'.format(str(self.embed_method).upper()))

        tile_style_menu = tool_menu.addMenu('Tile style')
        tile_style_menu.addAction(self.action_tile_AB_flag)
        tile_style_menu.addAction(self.action_tile_A_flag)
        tile_style_menu.addAction(self.action_tile_B_flag)

        tool_menu.addAction(self.action_set_coupling)

    def init_toolbar(self):
        ''' '''

        toolbar = QtWidgets.QToolBar()
        toolbar.setIconSize(QtCore.QSize(settings.ICO_SIZE, settings.ICO_SIZE))
        self.addToolBar(QtCore.Qt.LeftToolBarArea, toolbar)

        # construct actions
        action_qca_file = QtWidgets.QAction(self)
        action_qca_file.setIcon(
            QtGui.QIcon(settings.ICO_DIR+'qca_file.png'))
        action_qca_file.setStatusTip('Open QCA file...')
        action_qca_file.triggered.connect(self.load_qca_file)

        action_embed_file = QtWidgets.QAction(self)
        action_embed_file.setIcon(
            QtGui.QIcon(settings.ICO_DIR+'embed_file.png'))
        action_embed_file.setStatusTip('Open embedding file...')
        action_embed_file.triggered.connect(self.load_embed_file)

        action_chimera_file = QtWidgets.QAction(self)
        action_chimera_file.setIcon(
            QtGui.QIcon(settings.ICO_DIR+'chimera_file.png'))
        action_chimera_file.setStatusTip('Open chimera file...')
        action_chimera_file.triggered.connect(self.load_chimera_file)

        self.action_switch_adj = QtWidgets.QAction(self)
        self.action_switch_adj.setIcon(
            QtGui.QIcon(settings.ICO_DIR+'lim_adj.png'))
        self.action_switch_adj.setStatusTip('Switch to Limited Adjacency...')
        self.action_switch_adj.triggered.connect(self.switch_adjacency)
        self.action_switch_adj.setEnabled(False)

        self.action_embed = QtWidgets.QAction(self)
        self.action_embed.setIcon(
            QtGui.QIcon(settings.ICO_DIR+'embed.png'))
        self.action_embed.setStatusTip('Embed diplayed circuit...')
        self.action_embed.triggered.connect(self.embed_circuit)
        self.action_embed.setEnabled(False)

        self.action_del_embed = QtWidgets.QAction(self)
        self.action_del_embed.setIcon(
            QtGui.QIcon(settings.ICO_DIR+'del-embed.png'))
        self.action_del_embed.setStatusTip('Delete active embedding...')
        self.action_del_embed.triggered.connect(self.removeEmbedding)
        self.action_del_embed.setEnabled(False)

        toolbar.addAction(action_qca_file)
#        toolbar.addAction(action_embed_file)
        toolbar.addAction(action_chimera_file)
        toolbar.addAction(self.action_switch_adj)
        toolbar.addAction(self.action_embed)
        toolbar.addAction(self.action_del_embed)
        toolbar.addAction(self.action_export_coefs)

    def reset(self):
        '''Delete all embeddings and reset counters'''

        for ind in self.embeddings:
            self.active_embedding = ind
            self.removeEmbedding()
        self.active_embedding = -1
        self.embedding_count = 0

    def switch_adjacency(self):
        ''' '''
        if self.qca_active:
            self.full_adj = not self.full_adj
            ico_file = 'lim_adj.png' if self.full_adj else 'full_adj.png'
            sub_message = 'Limited' if self.full_adj else 'Full'
            self.action_switch_adj.setIcon(
                QtGui.QIcon(settings.ICO_DIR+ico_file))
            self.action_switch_adj.setStatusTip(
                'Switch to {0} Adjacency...'.format(sub_message))
            self.qca_widget.setAdjacency(self.full_adj)

    def switch_embedder(self, embed_method):
        if embed_method=='dense':
            self.switch_dense_embed()
        if embed_method=='layout':
            self.switch_layout_embed()
        if embed_method=='heur':
            self.switch_heur_embed()

    def switch_dense_embed(self):
        self.action_dense_embed_flag.setEnabled(False)
        self.action_layout_embed_flag.setEnabled(True)
        self.action_heur_embed_flag.setEnabled(True)
        self.embed_method = 'dense'

    def switch_layout_embed(self):
        self.action_dense_embed_flag.setEnabled(True)
        self.action_layout_embed_flag.setEnabled(False)
        self.action_heur_embed_flag.setEnabled(True)
        self.embed_method = 'layout'

    def switch_heur_embed(self):
        self.action_dense_embed_flag.setEnabled(True)
        self.action_layout_embed_flag.setEnabled(True)
        self.action_heur_embed_flag.setEnabled(False)
        self.embed_method = 'heuristic'

    def set_tile_style(self, style):
        ''' '''

        self.tile_style = style
        self.action_tile_AB_flag.setEnabled(True)
        self.action_tile_A_flag.setEnabled(True)
        self.action_tile_B_flag.setEnabled(True)
        if style==0:
            self.action_tile_AB_flag.setEnabled(False)
        elif style<0:
            self.action_tile_A_flag.setEnabled(False)
        else:
            self.action_tile_B_flag.setEnabled(False)

    def apply_tile_style(self, chimera_adj):
        ''' '''

        if self.tile_style == 0:
            return chimera_adj, 4

        check = lambda qb: (qb[3]%2) == (0 if self.tile_style < 0 else 1)

        chim_adj = {qb1: [qb2 for qb2 in chimera_adj[qb1] if check(qb2)]
            for qb1 in chimera_adj if check(qb1)}

        return chim_adj, 2

    def set_coupling(self):
        '''Change the coupling strength of couplers within vertex models'''

        # create popup dialog
        val, ok = QtWidgets.QInputDialog.getDouble(self, 'Dailog',
            'Coupling Strength:', value=self.coupling_strength)

        if ok and val > 0:
            self.coupling_strength = val


    def embed_circuit(self):
        '''Run embedding on displayed circuit into selected chimera
        sub-graph'''

        print('Running embedding...')

        try:
            # get chimera sub-graph
            M, N, chimera_adj, active_range = self.chimera_widget.getActiveGraph()

            # apply tile style
            chimera_adj, L = self.apply_tile_style(chimera_adj)

            # get qca parameters
            J, cells = self.qca_widget.prepareCircuit()

            # embedding object
            embedding = Embedding(self.qca_widget.filename)
            embedding.set_embedder(self.embed_method)
            embedding.set_chimera(chimera_adj, active_range, M, N)
            embedding.set_qca(J, cells, self.full_adj, self.qca_widget.spacing)

            # run embedding
            try:
                embedding.run_embedding()
            except Exception as e:
                if type(e).__name__ == 'KeyboardInterrupt':
                    print('Embedding interrupted...')
                else:
                    print('Something went wrong...')
                    print (traceback.print_exc())
                return
        except:
            print('\nUnexpected crash in embedding... possible disjoint graph')
            return

        if embedding.good:
            self.addEmbedding(embedding)
        else:
            print('Embedding failed...')

    def addEmbedding(self, embedding):
        '''Add an embedding object'''

        # enable relevant actions
        if len(self.embeddings) == 0:
            self.embeddings_menu.setEnabled(True)
            self.action_save_all.setEnabled(True)
            self.action_export_coefs.setEnabled(True)

        # get label for embedding in embedding menu
        label = os.path.basename(embedding.qca_file)

        # create new sub-menu if needed
        if label not in self.embedding_menus:
            self.embedding_menus[label] = self.embeddings_menu.addMenu(label)

        # create action for menu
        ind = int(self.embedding_count)
        func = lambda: self.switchEmbedding(ind)
        action = QtWidgets.QAction(str(self.embedding_count), self)
        action.triggered.connect(func)

        # add action to sub-menu
        self.embedding_menus[label].addAction(action)

        # store action for access/deletion
        self.embedding_actions[self.embedding_count] = action

        # add embedding to list of embeddings
        self.embeddings[self.embedding_count] = embedding

        # add embedding to chimera
        self.chimera_widget.addEmbedding(embedding, self.embedding_count)

        # set as active embedding
        self.switchEmbedding(ind)

        # update embedding_count
        self.embedding_count += 1

    def removeEmbedding(self):
        ''' '''

        if self.active_embedding == -1:
            return

        ind = self.active_embedding

        if ind not in self.embeddings:
            print('Attempted to delete a non-existing embedding...')
            return

        # special case if active embedding
        if ind == self.active_embedding:
            self.active_embedding = -1
            self.action_del_embed.setEnabled(False)

        embedding = self.embeddings.pop(ind)
        label = os.path.basename(embedding.qca_file)

        # clean nodes of embedding
        self.chimera_widget.resetNodes(embedding)

        # delete embedding object
        del(embedding)

        # delete action from sub-menu
        self.embedding_menus[label].removeAction(self.embedding_actions[ind])

        # delete action
        self.embedding_actions.pop(ind)

        # delete sub-menu if no more elements
        if self.embedding_menus[label].isEmpty():
            menu_action = self.embedding_menus[label].menuAction()
            self.embeddings_menu.removeAction(menu_action)
            self.embedding_menus.pop(label)

        # disable embeddings_menu if no embeddings
        if len(self.embeddings) == 0:
            self.action_save_all.setEnabled(False)
            self.action_export_coefs.setEnabled(False)
            self.embeddings_menu.setEnabled(False)

    def switchEmbedding(self, ind, color=True):
        '''Switch active embedding'''

        if ind in self.embeddings:
            # reanable embedding action
            if self.active_embedding != -1:
                self.embedding_actions[self.active_embedding].setEnabled(True)
            else:
                self.action_save_embedding.setEnabled(True)

            # allow deletion of active embedding
            self.action_del_embed.setEnabled(True)

            # disable new active embedding action
            self.embedding_actions[ind].setEnabled(False)

            # update active embedding
            self.active_embedding = ind
            self.qca_widget.updateCircuit(self.embeddings[ind].qca_file,
                                          self.embeddings[ind].full_adj)
            if self.embeddings[ind].full_adj != self.full_adj:
                self.switch_adjacency()

            if self.embeddings[ind].embed_method != self.embed_method:
                self.switch_embedder(self.embeddings[ind].embed_method)

            # default coloring
            if color:
                # color nodes, no cell selected (assume no -1 cell)
                self.chimera_widget.selectNodes(self.embeddings[ind], -1)

    # FILE IO

    def create_coef_pol_file(self, fname, pol_sets):
        ''' '''

        try:
            fp = open(fname, 'w')
        except:
            print('Failed to open file: {0}'.format(fname))
            raise IOError

        for pol_set in pol_sets:
            for ind in pol_set:
                pols = pol_set[ind]
                fp.write('{0}: {1}\n'.format(ind, str(pols)))
            fp.write('\n')

        fp.close()

    def save_coef_file(self, hq, Jq, fname):
        ''' '''

        try:
            fp = open(fname, 'w')
        except:
            print('Failed to open file: {0}'.format(fname))
            raise IOError

        # chimera size
        M, N, L = self.chimera_widget.M, self.chimera_widget.N, 4
        Nqbits = 2*M*N*L

        fp.write('{0}\n'.format(Nqbits))

        # h parameters
        for qb in sorted(hq):
            qbl = tuple_to_linear(qb, M, N, L=L, index0=False)
            if hq[qb] != 0:
                fp.write('{0} {1} {2}\n'.format(qbl, qbl, round(hq[qb],3)))

        # J parameters
        for qb1, qb2 in sorted(Jq):
            qbl1 = tuple_to_linear(qb1, M, N, L=L, index0=False)
            qbl2 = tuple_to_linear(qb2, M, N, L=L, index0=False)
            if Jq[(qb1, qb2)] != 0:
                fp.write('{0} {1} {2}\n'.format(qbl1, qbl2,
                         round(Jq[(qb1, qb2)], 3)))

        fp.close()

    def create_embed_file(self, fname):
        '''Create an info file for all current embeddings'''

        try:
            fp = open(fname, 'w')
        except:
            print('Failed to open file: {0}'.format(fp))
            raise IOError

        # header
        chim_file = os.path.relpath(self.chimera_file, os.path.dirname(fname))
        fp.write('chimera_file: {0}\n\n'.format(chim_file))

        # embedding files
        for ind in self.embeddings:
            fp.write('{0}: {0}.txt\n'.format(ind))

        fp.close()

    def save_active_embedding(self):
        '''save the active embedding'''

        if self.active_embedding == -1:
            print('Trying to save nothing....should not have happened')
            return

        fname = self.getSaveFileName('Save active embedding', self.embed_dir)

        if not fname:
            return

        embedding = self.embeddings[self.active_embedding]
        self.save_embedding(embedding, fname)

    def save_embedding(self, embedding, fname):
        '''Save a single embedding to file'''

        try:
            fp = open(fname, 'w')
        except:
            print('Failed to open file: {0}'.format(fp))
            raise IOError

        # chimera file
        chim_file = os.path.relpath(self.chimera_file, os.path.dirname(fname))
        fp.write('chimera_file: {0}\n'.format(chim_file))

        # qca file
        qca_file = os.path.relpath(embedding.qca_file, os.path.dirname(fname))
        fp.write('qca_file: {0}\n\n'.format(qca_file))

        # adjacency and embedding type
        fp.write('full_adj: {0}\n'.format(embedding.full_adj))
        fp.write('embed_method: {0}\n\n'.format(embedding.embed_method))

        # chimera parameters
        fp.write('M: {0}\n'.format(embedding.M))
        fp.write('N: {0}\n'.format(embedding.N))
        fp.write('L: {0}\n'.format(embedding.L))
        fp.write('M0: {0}\n'.format(embedding.active_range['M'][0]))
        fp.write('N0: {0}\n\n'.format(embedding.active_range['N'][0]))

        # cell models
        for cell in embedding.models:
            fp.write('{0}: {1}\n'.format(cell,
                     ';'.join(str(qb) for qb in embedding.models[cell])))

        fp.close()

    def save_all_embeddings(self):
        '''Save all embeddings to a directory with an embed (summary) file'''

        # prompt for directory name
        dir_name = self.getExistingDirectory('Create/Select empty directory...', self.embed_dir)

        if not dir_name:
            return

        # convert to standard form
        dir_name = os.path.join(os.path.normpath(dir_name), '')

        # update embed home directory
        self.embed_dir = dir_name

        # if directory does not exist, create it
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        # if files exist in directory, prompt user
        files = [f for f in os.listdir(dir_name)\
            if os.path.isfile(os.path.join(dir_name, f))]

        if len(files) > 0:
            reply = QtWidgets.QMessageBox.question(self, 'Message',
            'This directory already contains content that will be deleted. Do\
            you want to continue?',
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.Cancel,
            QtWidgets.QMessageBox.Cancel)

            if reply == QtWidgets.QMessageBox.Yes:
                for f in files:
                    os.remove(os.path.join(dir_name, f))
            else:
                return

        try:
            # create embed file
            self.create_embed_file(os.path.join(dir_name, 'summary.embed'))

            # save each embedding
            for ind in self.embeddings:
                fname = os.path.join(dir_name, '{0}.txt'.format(ind))
                self.save_embedding(self.embeddings[ind], fname)
        except IOError:
            print('Failed to save embeddings...')
            return

    def load_embedding(self, fname, chimera_adj):
        ''' '''

        # create embedding
        embedding = Embedding()
        try:
            embedding.from_file(fname, self.chimera_file, chimera_adj)
            # set qca parameters
            self.qca_widget.updateCircuit(embedding.qca_file, embedding.full_adj)
            J, cells = self.qca_widget.prepareCircuit()
            embedding.set_qca(J, cells, embedding.full_adj)
        except:
            print('Failed to load embedding')
            return

        self.addEmbedding(embedding)

    def load_embed_file(self):
        '''Prompt filename for embed file'''

        fname = self.getOpenFileName('Select Embedding File', self.embed_dir,
                'EMBED (*.embed);; All files (*)')

        if not fname:
            return

        # update embed home directory
        fdir = os.path.dirname(fname)
        self.embed_dir = fdir

        try:
            fp = open(fname, 'r')
        except:
            print('Failed to open file: {0}'.format(fname))
            return None

        # parse file
        info = {}
        inds = []
        for line in fp:
            if '#' in line or len(line)<3:
                continue
            key, data = [x.strip() for x in line.split(':')]
            info[key] = data
            if key.isdigit():
                inds.append(int(key))
        fp.close()

        # delete all embeddings
        self.reset()
        ndir = os.path.dirname(fname)
        chim_file = os.path.normpath(os.path.join(ndir, info['chimera_file']))
        self.chimera_widget.updateChimera(chim_file)
        M, N, chimera_adj, active_range = self.chimera_widget.getActiveGraph()
        for ind in inds:
            ndir = os.path.dirname(fname)
            fn = os.path.normpath(os.path.join(ndir, info[str(ind)]))
            self.load_embedding(fn, chimera_adj)

        if not self.qca_active:
            self.qca_active = True
            self.action_embed.setEnabled(True)
            self.action_switch_adj.setEnabled(True)

    def load_qca_file(self):
        '''Prompt filename for qca file'''

        fname = self.getOpenFileName('Select QCA File', self.qca_dir)

        if not fname:
            return

        # update qca home directory
        fdir = os.path.dirname(fname)
        self.qca_dir = fdir

        self.qca_widget.updateCircuit(fname, self.full_adj)

        # disable old embedding
        self.chimera_widget.unclickNodes()
        if self.active_embedding != -1:
            self.embedding_actions[self.active_embedding].setEnabled(True)
        self.active_embedding = -1
        self.action_save_embedding.setEnabled(False)

        if not self.qca_active:
            self.qca_active = True
            self.action_embed.setEnabled(True)
            self.action_switch_adj.setEnabled(True)
            self.action_save_qca_svg.setEnabled(True)

    def load_chimera_file(self):
        '''Prompt filename for chimera structure'''

        fname = self.getOpenFileName('Select Chimera File', self.chimera_dir)

        if not fname:
            return

        # update chimera home directory
        fdir = os.path.dirname(fname)
        self.chimera_dir = fdir

        self.chimera_widget.updateChimera(fname)
        self.chimera_file = os.path.relpath(fname)

    def export_coefs(self):
        '''Determine the smallest set of files which need to be produced to
        allow for all unique input combinations to all independent embeddings.
        Save each to a file'''

        # prompt for directory name
        dir_name = self.getExistingDirectory('Create/Select empty directory...', self.coef_dir)

        if not dir_name:
            return

        # convert to standard form
        dir_name = os.path.join(os.path.normpath(dir_name), '')

        # update embed home directory
        self.coef_dir = dir_name

        # if directory does not exist, create it
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        # if files exist in directory, prompt user
        files = [f for f in os.listdir(dir_name)\
            if os.path.isfile(os.path.join(dir_name, f))]

        if len(files) > 0:
            reply = QtWidgets.QMessageBox.question(self, 'Message',
            'This directory already contain content that will be deleted. Do\
            you want to continue?',
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.Cancel,
            QtWidgets.QMessageBox.Cancel)

            if reply == QtWidgets.QMessageBox.Yes:
                for f in files:
                    os.remove(os.path.join(dir_name, f))
            else:
                return

        # generate list of all polarizations
        driver_pols = {ind: self.embeddings[ind].generate_driver_pols()
            for ind in self.embeddings}
        print(driver_pols)
        Npol = max([len(pols) for pols in driver_pols.values()])
        for ind in driver_pols:
            N = len(driver_pols[ind])
            if N == 0:
                driver_pols[ind] = [[]]*Npol
            else:
                driver_pols[ind] *= Npol/N
        pol_sets = []
        for n in range(Npol):
            pol_sets.append({ind: driver_pols[ind][n] for ind in driver_pols})

        # find the coefficients for each set of polarizations
        HQs = []
        JQs = []
        J_inner = -self.coupling_strength
        for pol_set in pol_sets:
            hq = {}
            Jq = {}
            for ind in self.embeddings:
                embedding = self.embeddings[ind]
                pols = pol_set[ind]
                if len(pols)==0:
                    pols = [[]]
                h, J = embedding.generate_coefs(pols, J_inner=J_inner)
                # add new parameters
                for q in h:
                    hq[q] = h[q]
                for q1, q2 in J:
                    Jq[(q1, q2)] = J[(q1, q2)]

            HQs.append(hq)
            JQs.append(Jq)

        try:
            # create pol file
            self.create_coef_pol_file(os.path.join(dir_name,
                                                   'pols.info'), pol_sets)

            # create each coef file
            for i in range(Npol):
                hq, Jq = HQs[i], JQs[i]
                fname = os.path.join(dir_name, 'coefs{0}.txt'.format(i))
                self.save_coef_file(hq, Jq, fname)

        except IOError:
            print('Failed to save coefficient files...')
            return

    # VECTOR GRAPHICS

    def save_qca_svg(self):
        ''' '''

        fname = self.getSaveFileName('Save SVG file...', self.svg_dir)

        if fname:
            self.svg_dir = os.path.dirname(fname)

        self.qca_widget.save_svg(fname)

    def save_chimera_svg(self):
        ''' '''

        fname = self.getSaveFileName('Save SVG file...', self.svg_dir)

        if fname:
            self.svg_dir = os.path.dirname(fname)

        self.chimera_widget.save_svg(fname)

    # EVENT HANDLING

    def closeEvent(self, e):
        '''Handle main window close event'''

        reply = QtWidgets.QMessageBox.question(
            self, 'Message', 'Are you sure you want to quit?',
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.Cancel,
            QtWidgets.QMessageBox.Cancel)

        if reply == QtWidgets.QMessageBox.Yes:
            e.accept()
        else:
            e.ignore()

    def keyPressEvent(self, e):
        ''' '''
        if e.key() == QtCore.Qt.Key_E:
            self.embed_circuit()

    # HACKS to fix filename getters in PyQt5

    def getSaveFileName(self, msg, rdir, flt=''):

        x = QtWidgets.QFileDialog.getSaveFileName(self, msg, rdir, flt)
        try:
            x = next(iter(x))
        except:
            pass
        return str(x)

    def getOpenFileName(self, msg, rdir, flt=''):
        x = QtWidgets.QFileDialog.getOpenFileName(self, msg, rdir, flt)
        try:
            x = next(iter(x))
        except:
            pass
        return str(x)

    def getExistingDirectory(self, msg, rdir):
        return str(QtWidgets.QFileDialog.getExistingDirectory(self, msg, rdir))
