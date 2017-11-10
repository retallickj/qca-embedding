#!/usr/bin/env python

#---------------------------------------------------------
# Name: classes.py
# Purpose: File constaining functional classes
# Author:	Jacob Retallick
# Created: 12.01.2015
# Last Modified: 12.01.2015
#---------------------------------------------------------

import re
import os
import traceback
import numpy as np

from core.chimera import tuple_to_linear, linear_to_tuple
import core.core_settings as settings
from core.dense_embed.assign import assign_parameters

# try to import different embedding methods
embedders = {'dense': True,
             'layout': True,
             'heur': True}
try:
    from core.dense_embed.embed import denseEmbed, setChimera
    from core.dense_embed.convert import convertToModels
except Exception as e:
    print('Could not load dense embedding method...')
    embedders['dense'] = False

try:
    from layout_embed.embed import layoutEmbed, layoutConfiguration
    from layout_embed.embed import setProblem, setTarget, layoutToModels
except Exception as e:
    print('Could not load layout embedding method...')
    embedders['layout'] = False
    print (traceback.print_exc())

try:
    from dwave_sapi2.embedding import find_embedding
except Exception as e:
    print('Could not load heuristic embedding method...')
    embedders['heur'] = False

# echo available embedder methods
print('Emedder options:')
for key, flag in embedders.iteritems():
    print('\t{0}:\t{1}'.format(key.upper(), 'Enabled' if flag else 'Disabled'))

SABOTAGE = False

def get_embedder_flags():
    return embedders

class Embedding:
    '''Container class for an embedding'''

    def __init__(self, qca_file=None):
        '''Initiate an embedding object'''

        # flags and general parameters
        self.qca_file = qca_file    # name of qca file embedded
        self.embed_file = None      # file to save embedding into
        self.coef_dir = None        # directory to save coefficient files into
        self.good = False           # flag for good embedding

        self.full_adj = True        # flag for full adjacency
        self.embed_method = 'dense' # value for for embedding method

        self.dense_trials = 1   # number of allowed dense placement trials

        # QCA and Chimera structure
        self.qca_adj = {}       # adjacency dict for qca circuit (all cells)
        self.chimera_adj = {}   # adjacency dict for chimera graph
        self.M = None           # number of rows in Chimera graph
        self.N = None           # number of columns in Chimera graph
        self.L = None
        self.active_range = {}  # for mapping between suba dn full graphs

        # QCA cell identifications
        self.cells = {}         # dict of cells for reference
        self.drivers = set()    # list of driver cells
        self.fixed = set()      # list of fixed cells
        self.normal = set()     # list of normal and output cells

        # embedding parameters
        self.fixed_pols = {}    # polarizations of fixed cells
        self.models = {}        # models for each qca cell

        # export parameters
        self.J = np.zeros(0)    # array of cell interactions (all cells)
        self.spacing = 1.       # spacing in layout input

    # EMBEDDING METHODS

    def run_embedding(self):
        ''' '''

        if self.embed_method == 'dense':
            self.run_dense_embedding()
        elif self.embed_method == 'layout':
            self.run_layout_embedding()
        elif self.embed_method == 'heuristic':
            self.run_heur_embedding()

        nqbits = sum(len(x) for x in self.models.values())
        ncells = len(self.models)

        print('Used {0} qbits for {1} cells'.format(nqbits, ncells))

    def run_dense_embedding(self, full_adj=True):
        '''Setup and run the Dense Placement algorithm'''

        # update embedding type in case direct call
        self.embed_method = 'dense'

        # format embedding parameters
        setChimera(self.chimera_adj, self.M, self.N, self.L)
        active_cells, qca_adj = self.get_reduced_qca_adj()

        # run a number of embedding and choose the best
        embeds = []
        for trial in xrange(settings.DENSE_TRIALS):
            print('Trial {0}...'.format(trial)),
            try:
                cell_map, paths = denseEmbed(qca_adj, write=False)
                print('success')
            except Exception as e:
                if type(e).__name__ == 'KeyboardInterrupt':
                    raise KeyboardInterrupt
                print('failed')
                continue
            embeds.append((cell_map, paths))

        if len(embeds) == 0:
            self.good = False
            return

        # sort embedding by number of qubits used (total path length)
        ind = -1 if SABOTAGE else 0
        cell_map, paths = sorted(embeds,
                                 key=lambda x: sum([len(p) for p in x[1]]))[ind]
        self.good = True

        # get cell models
        print('Converting to models...')
        models, max_model = convertToModels(paths, cell_map)
        print('done')

        self.models = {k: models[k]['qbits'] for k in models}

    def run_layout_embedding(self):
        '''Setup and run the Layout-Aware Placement algorithm'''

        # update embedding type in case direct call
        self.embed_method = 'layout'

        stats = {}
        configuration = {}
        configuration['M'] = self.M
        configuration['N'] = self.N
        configuration['L'] = self.L
        configuration['CIRCUIT'] = self.qca_file.split('/')[-1]
        layoutConfiguration(configuration)

        active_cells, qca_adj = self.get_reduced_qca_adj()

        try:
            setProblem(qca_adj, self.cells, self.spacing)
            setTarget(self.chimera_adj)
            self.good, cell_map = layoutEmbed(configuration, stats)
            if self.good:
                print('Layout-Aware Embedding Successful')
        except Exception as e:
            self.good = False
            if type(e).__name__ == 'KeyboardInterrupt':
                raise KeyboardInterrupt
            print('Layout-Aware Embedding Failed')
            print (traceback.print_exc())
            return

        self.models = layoutToModels(cell_map)

    def run_heur_embedding(self, full_adj=True):
        '''Setup and run the Heuristic algorithm'''

        # update embedding type in case direct call
        self.embed_method = 'heur'

        active_cells, qca_adj = self.get_reduced_qca_adj()
        S_size = len(qca_adj)
        A_size = len(self.chimera_adj)

        # construct S, the problem adjacency edge list
        S = set()
        smap = {c:i for i,c in enumerate(active_cells)}
        for c1, adj in qca_adj.items():
            for c2 in adj:
                if c1<c2:
                    S.add((smap[c1],smap[c2]))

        # construct A, the chimera adjacency edge list
        A = set()
        for qb1 in self.chimera_adj:
            for qb2 in self.chimera_adj[qb1]:
                l1 = tuple_to_linear(qb1, self.M, self.N, L=self.L, index0=True)
                l2 = tuple_to_linear(qb2, self.M, self.N, L=self.L, index0=True)
                A.add((l1, l2))

        try:
            print 'Running heuristic embedding'
            #models = find_embedding(S, S_size, A, A_size)
            models = find_embedding(S, A)
        except Exception as e:
            print(e.message())

        print 'Embedding finished'
        self.good = len(models) == S_size

        # map models to standard format
        mapper = lambda ind: linear_to_tuple(ind, self.M, self.N,
                                             L=self.L, index0=True)
        self.models = {active_cells[i]: [mapper(c) for c in model]
            for i,model in enumerate(models)}


    # PARAMETER ACCESS

    def set_embedder(self, method='dense'):
        '''Set embedder type'''
        self.embed_method = method

    def set_chimera(self, adj, active_range, M, N, L=4):
        '''Set the Chimera graph to embed into'''

        self.M = M
        self.N = N
        self.L = L
        self.active_range = active_range
        self.chimera_adj = adj

    def set_qca(self, J, cells, full_adj=True, spacing=1.):
        '''Set up the qca structure'''

        self.full_adj = full_adj

        self.qca_adj = {i: J[i].nonzero()[0].tolist()
            for i in xrange(J.shape[0])}

        # driver cells
        self.drivers = set([i for i in self.qca_adj if cells[i].driver])

        # fixed cells
        self.fixed = set([i for i in self.qca_adj if cells[i].fixed])

        # normal cells
        self.normal = set([i for i in self.qca_adj if cells[i].normal])

        # output cells
        self.outputs = set([i for i in self.qca_adj if cells[i].output])

        # store cell information
        self.cells  = {i: cells[i] for i in self.qca_adj}

        # fixed cell polarizations
        self.fixed_pols = {cell: cells[cell].pol for cell in self.fixed}

        # coupler weights
        self.J = J

        # QCA layout spacing
        self.spacing = spacing


    def generate_driver_pols(self, full=True):
        '''generate a list of all unique driver polarization sets'''

        N = len(self.drivers)   # number of driver cells
        if N == 0:
            return [[]]
        count = 2**N if full else 2**(N-1)
        bins = [format(i, '#0{0}b'.format(N+2))[2::] for i in range(2**N)]
        pols = [[(int(x)*2-1) for x in b] for b in bins]
        return pols

    def generate_coefs(self, pols, J_inner=-1):
        '''Generate h and J parameters for a given set of driver polarizations'''

        # set up input polarization vector
        P = np.asmatrix(np.zeros([len(self.qca_adj),1], dtype=float))

        # fixed cell contributions
        for cell in self.fixed:
            P[cell, 0] = self.fixed_pols[cell]

        # driver cell contributions
        drivers = sorted(self.drivers)  # list of drivers
        for i in range(len(drivers)):
            P[drivers[i], 0] = pols[i]

        # get h coefficients
        h = np.round(np.asmatrix(self.J)*P, 3)

        # get reduce J matrix
        inds = sorted(self.normal)
        h = {ind: float(h[ind]) for ind in inds}
        J = {}
        for i1 in inds:
            J[i1] = {}
            for i2 in inds:
                if self.J[i1, i2] != 0:
                    J[i1][i2] = np.round(self.J[i1, i2], 3)

        # compute qubit parameters
        hq, Jq = assign_parameters(h, J, self.models, self.chimera_adj,
                                   flip_J=True, J_inner=J_inner)

        # correct for tile offset
        M0, N0 = self.active_range['M'][0], self.active_range['N'][0]
        mapping = lambda q: (q[0]+M0, q[1]+N0, q[2], q[3])
        hq = {mapping(q): hq[q] for q in hq}
        Jq = {(mapping(q1), mapping(q2)): Jq[(q1, q2)] for q1, q2 in Jq}

        return hq, Jq

    def get_reduced_qca_adj(self):
        '''Get a reduced form of qca_adj only for non-driver/fixed cells'''

        # check function for membership in drivers or fixed
        check = lambda cell: cell not in self.drivers.union(self.fixed)

        reduced_adj = {c1: [c2 for c2 in self.qca_adj[c1] if check(c2)]
            for c1 in self.qca_adj if check(c1)}

        return sorted(reduced_adj), reduced_adj

    # FILE IO

    def from_file(self, fname, chimera_file, chimera_adj):
        '''Load an embedder object from file relative to main directory'''

        try:
            fp = open(fname, 'r')
        except:
            print('Failed to open file: {0}'.format(fname))
            raise IOError

        # parse file
        info = {}
        cells = []
        for line in fp:
            if '#' in line or len(line) < 3:
                continue
            key, data = [x.strip() for x in line.split(':')]
            info[key] = data
            if key.isdigit():
                cells.append(key)
        fp.close()

        # process info
        ndir = os.path.dirname(fname)
        chim_file = os.path.normpath(os.path.join(ndir, info['chimera_file']))

        if not os.path.samefile(chim_file, chimera_file):
            print('Chosen embedding is not native to this chimera graph')
            return

        self.qca_file = os.path.normpath(os.path.join(ndir, info['qca_file']))

        # flags
        self.full_adj = info['full_adj'] == 'True'
        #TODO: Not tested with other methods
        self.embed_method = info['embed_method']

        # chimera parameters
        self.M = int(info['M'])
        self.N = int(info['N'])
        self.L = int(info['L'])

        M0 = int(info['M0'])
        N0 = int(info['N0'])
        self.active_range = {'M': [M0, M0+self.M],
                             'N': [N0, N0+self.N]}

        # refine chimera adjacency
        tile_check = lambda m, n, h, l: \
            m >= self.active_range['M'][0] and\
            m < self.active_range['M'][1] and\
            n >= self.active_range['N'][0] and\
            n < self.active_range['N'][1]

        chimera_adj = {k1: [k2 for k2 in chimera_adj[k1] if tile_check(*k2)]
            for k1 in chimera_adj if tile_check(*k1)}

        offset = lambda m, n, h, l: (m-M0, n-N0, h, l)
        chimera_adj = {offset(*k1): [offset(*k2) for k2 in chimera_adj[k1]]
            for k1 in chimera_adj}
        self.chimera_adj = chimera_adj

        # get models
        models = {}
        regex = re.compile('[0-9]+')
        str_to_tuple = lambda s: tuple([int(x) for x in regex.findall(s)])
        for cell in cells:
            models[int(cell)] = [str_to_tuple(s) for s in info[cell].split(';')]
        self.models = models
