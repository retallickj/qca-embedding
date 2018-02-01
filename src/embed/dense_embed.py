#!/usr/bin/env python
# encoding: utf-8

'''
Classes for embedding graphs onto the Chimera architecture using the Dense
Placement algorithm.
'''

from __future__ import print_function   # for verbose functionality

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2018-01-30'      # last update

from core.utility import dget, range_product
from collections import defaultdict
import networkx as nx
import numpy as np


# exceptions

class PlacementError(Exception):
    '''Custom Exception for all failure modes of the Embedding algorithm'''
    def __init__(self, message='', mode=0):
        super(PlacementError, self).__init__(message)
        self.mode = mode

class SourceError(Exception):
    '''Custom Exception for problems with the source graph'''
    pass


# structs

def echo_struct(T, sep=' :: '):
    '''Read out all the class and instance variables of a class as name=val pairs'''
    return sep.join('{0}={1}'.format(attr, getattr(T, attr)) for attr in dir(T)
                if not callable(getattr(T, attr)) and not attr.startswith('__'))

class CellPars:
    '''All relevant parameters for a given cell during embedding'''

    qbit = None     # qbit assigned to the cell, None if no qbit assigned
    placed = False  # True if cell has been placed
    nadj = -1       # current number of unplaced adjacent cells

    def __init__(self, nadj=None):
        if isinstance(nadj, int) and nadj >= 0:
            self.nadj = nadj

    def __str__(self):
        return echo_struct(self)


class QubitPars:
    '''All relevant parameters for a given qubit during embedding'''

    cell = None         # cell the qubit is assigned to, else None
    reserved = False    # qbit is reserved for pending adjacency
    taken = False       # qbit taken for routing
    prox = set()        # set of assigned adjacent qubits
    c_in = set()        # set of free internal qubits
    c_out = -1          # number of free external qubits
    paths = set()       # paths containing the qubit

    def __str__(self):
        return echo_struct(self)

# logging

class Logger:
    '''Handler for writing the embedder status to either stdout or a log file'''

    def __init__(self, fn=None, append=False):
        '''Attempt to initialise a Logger instance for the given log file path.
        If no filepath is given, outputs to stdout.'''

        self.fn = fn
        self.fp = None
        self.append = bool(append)

    def start(self):
        '''Start the log'''

        if self.fn is None:
            self.fp = None
        else:
            try:
                self.fp = open(self.fn, 'a' if self.append else 'w')
            except:
                print('Failed to open logger for file: {0}'.format(self.fn))
                self.fp = None

    def end(self):
        '''End the log'''

        if self.fp is not None:
            self.fp.close()

    def log(self, s):
        '''Write the given string to the log. No endline is added.'''

        if self.fp is not None:
            self.fp.write(s)
        else:
            print(s, end='')

# main embedder class

class DenseEmbedder:
    '''Handler class for running the dense placement embedding algorithm.
    Attempts to find an embedding of the source graph into a target Chimera
    graph.'''

    # seeding parameters
    INIT_TRIALS = 5         # number of attempts to find initial seed per trial
    FIRST_PROB_POW = 4.     # power for firstCell method
    FIRST_QBIT_SIG = .3     # gaussian dist. range for tile select
    FIRST_QBIT_ATTEMPTS = 5 # number of attempt to find starting qubit

    # search costs
    IN_TILE_COST = 1        # new qubit in the same tile
    OUT_TILE_COST = 1.9     # new qubit in different tile
    EDGE_REP_COST = 0.5     # linear cost for tiles closer to processor edge

    # seam ranking costs
    SEAM_QBIT_COST = 30     # assigned qubit moved to broken qubit
    SEAM_PATH_COST = 40     # path broken
    SEAM_EXT_COST = 0       # path extension
    SEAM_DIST_COST = 30     # distance between seam and average qubit position

    MAX_SEARCH_COUNT = 3    # maximum number of multisource search attempts

    def __init__(self, chimera=None, logfile=None, append=False):
        '''Initialise a DenseEmbedder instance.

        parameters:
            chimera : optional Chimera graph
            logfile : optional filename if verbose status written to a log
            append  : True if should append log to logfile
        '''

        self.set_logger(logfile, append)

        if chimera is not None:
            self.set_chimera(chimera)
        else:
            self.chimera = None

    def set_logger(self, logfile=None, append=False):
        '''Override the logger'''

        self._logger = Logger(logfile, append)

    def set_chimera(self, chimera):
        '''Update the Chimera graph for the Embedder

        input:
            chimera : Chimera graph. All nodes in the graph must have a 'tup'
                      field which is equivalent to a 4-tuple of form (m,n,h,l)
                      describing the qubit in the Chimera graph (see core.chimera)
        '''

        # check that graph is valid
        mins, maxs, L = self._check_chimera(chimera)
        if L is None:
            print('Invalid chimera format')
            return

        # store chimera and properties
        self.chimera = self._remap_chimera(chimera)
        self.origin = mins          # top-left tile
        self.M, self.N = [1+b-a for a,b in zip(mins, maxs)]
        self.L = L

    def embed(self, source, **kwargs):
        '''Attempt to find an embedding of the given source graph onto the
        current Chimera graph.

        inputs:
            source  : source graph, can be any networkx.Graph or derived

        kwargs:
            verbose : if True, verbose logging for the embedding
            ntrials : number of embedding attempts
            best    : if True, return only the most economical embedding
        '''

        # check that chimera is defined
        if self.chimera is None:
            print('Chimera Graph has not been specified')
            return []

        # handle default arguments
        verbose = dget(kwargs, 'verbose', False, mp=bool)
        ntrials = dget(kwargs, 'ntrials', 1, mp=int, chk=lambda n: n>0)
        best = dget(kwargs, 'best', False, mp=bool)

        # validate optional arguments
        self.vlog = self.log if verbose else lambda *a, **k: None

        # initialize solver for source graph
        try:
            self._initialise(source)
        except SourceError:
            return []

        # main loop
        solutions = []
        for trial in range(ntrials):
            try:
                self._embed_trial()
                # embedding should now be stored in the class, query it
                solutions.append(self._get_embedding())
            except PlacementError as e:
                continue
            except KeyboardInterrupt:
                print('Keyboard Interrupt: returning completed solutions')
                break
            except:
                print('Unknown error... aborting')
                return None

        if best:
            return min(solutions, key=self._embed_eval())
        else:
            return solutions


    def _initialise(self, source):
        '''Prepare the Embedder for the given source adjacency graph'''

        # start the Logger
        self._logger.start()

        self.vlog('Initialising embedder for source graph... ')

        # check source
        if not self._check_source(source):
            raise SourceError('Invalid source graph')
        self.source = source

        self.vlog('done\n')


    def _reset(self):
        '''Reset the Embedder for a new trial'''

        self.vlog('Resetting embedder... ')

        # set up cell parameters
        self._cps = {cell: CellPars(deg) for cell, deg in self.source.degree()}

        # set up qubit parameters
        self._qps = {qbit: QubitPars() for qbit in self.chimera}

        self._paths = {}    # paths between cell qubits
        self._tile_occ = {k: defaultdict(int) for k in 'mn'}    # qubits per row/col
        self._vacancy = [-1,]*4     # number of free cols/rows of tiles [L,R,D,U]

        self.vlog('done\n')


    def _first_cell(self):
        '''Pick the first cell to place'''


        self.log('Selecting first cell... ')
        keys = self._cps.keys()
        worth = lambda key: pow(self._cps[key].nadj, self.FIRST_PROB_POW)
        probs = np.array([worth(key) for key in keys], dtype=float)
        cell = np.random.choice(keys, p=probs/np.sum(probs))
        self.log('done\n')

        return cell


    def _first_qbit(self, cell):
        '''Find the qubit seed for the first cell'''

        self.log('Attempting to select first qubit... ')

        # 1D Gaussian PDFs for each tile axis
        X = [x0+np.arange(k) for x0, k in zip(self.origin, [self.M, self.N])]
        Z = [np.linspace(-.5*(k-1), .5*(k-1), k) for k in [self.M, self.N]]
        PDFs = [np.exp(-z*z/(2*self.FIRST_QBIT_SIG)) for z in Z]
        PDFs = [pdf/np.sum(pdf) for pdf in PDFs]

        # make attempts to pick first qubit
        Chim = self.chimera     # shorthand
        for attempt in range(self.FIRST_QBIT_ATTEMPTS):
            # pick tile according to the PDFs
            m,n = [np.random.choice(x, p=pdf) for x,pdf in zip(X,PDFs)]
            # randomly loop through qubits in tile to check if any have enough
            # neighbours to satisfy cell connectivity
            order = [(h,l) for h,l in range_product(2, self.L)]
            np.random.shuffle(order)
            for h,l in order:
                qbit = (m,n,h,l)
                if qbit in Chim and Chim.degree(qbit) >= self._cps[cell].nadj:
                    self.log('done\n')
                    return qbit

        self.log('failed all {0} attempts\n'.format(self.FIRST_QBIT_ATTEMPTS))
        return None

    def _assign_qbit(self, cell, qbit):
        '''Associate the given qubit and cell'''

        # decrement nadj for all adjacent cells
        for neigh in self.source.neighbors(cell):
            self._cps[neigh].nadj -= 1

        # associate cell and qubit
        self._cps[cell].qbit = qbit
        self._qps[qbit].cell = cell

        # update tile_occ and vacancy
        if not self._qps[qbit].reserved:    # not already counted
            self._update_tile_occ([qbit])
            self._set_vacancy()

        # update parameters
        self._cps[cell].placed = True
        self._qps[qbit].taken = True
        self._qps[qbit].assigned = True

        # make adjacent qubits aware of placed cell (for reserve check)
        for qb in self.chimera.neighbors(qbit):
            self._qps[qbit].prox.add(qbit)

        # disable qubit from Routing framework
        Routing.disable_qubits([qbit])



    def _assign_paths(self, paths):
        '''Set all relevant flags for the given routed paths'''
        pass

    def _reserve_qbits(self, qbits):
        '''For each of a set of qubits, check if the adjacent qubits need to
        be reserved in order to satisfy connectivity constraints. Reserve as
        appropriate.'''
        pass

    def _place_cell(self, cell):
        '''Attempt to find a placement for the given cell.

        returns:
            qbit    : label of qubit to assign to cell
            paths   : list of paths from previously placed neighbours to the
                      assigned qubit
        '''
        pass

    def _update_tile_occ(self, qbits, dec=False):
        '''Either add or remove the set of qubits from the tile row/col
        occupancy tally'''

        x = -1 if dec else -
        for qb in qbits:
            self._tile_occ['m'][qb[0]] += x
            self._tile_occ['n'][qb[1]] += x

    def _set_vacancy(self):
        '''Update the number of free rows/cols from the tile occupancy'''

        # compute left/right vacancy
        left = min(k for k,c in self._tile_occ['n'].items() if c>0))
        right = max(k for k,c in self._tile_occ['n'].items() if c>0))
        left,right = [x - self.origin[1] for x in [left,right]]

        # compute down/up vacancy
        down = min(k for k,c in self._tile_occ['m'].items() if c>0))
        up = max(k for k,c in self._tile_occ['m'].items() if c>0))
        down,up = [x - self.origin[0] for x in [down,up]]

        self._vacancy = [left, right, down ,up]

    def _embed_trial(self):
        '''Run a single embedding trial'''

        # Reset the embedder for the new trial
        self._reset()

        # attempt to find initial seed
        for init_trial in range(self.INIT_TRIALS):
            cell = self._first_cell()       # seed cell
            qbit = self._first_qbit(cell)   # seed qbit
            if qbit is None:
                break
        else:
            raise PlacementError('Failed to place initial seed')

        # assign seed cell to seed qubit
        self._assign_qbit(cell, qbit)
        self._reserve_qbits([qbit])

        # current list of cells to place
        do_now = sorted(self.source.neighbors(cell),
                        key=lambda c: self._cps[c].nadj, reverse=True)
        do_next = set()

        # if do_now is non-empty, there are still cells to place
        while do_now:

            self.log('toDo: {0}'.format(str(do_now)))

            # place each cell in do_now
            for cell in do_now:

                # find qbit and paths from placed cells
                qbit, paths = self._place_cell(cell)

                # abort on failed placement
                if qbit is None:
                    raise PlacementError('Failed to place cell {0}'.format(cell))

                # assign qubit and paths
                self._assign_qbit(cell, qbit)
                self._assign_paths(paths)

                # reservations
                self._reserve_qbits([qbit])

                # add unplaced neighbours to do_next
                for c2 in self.source.neighbors(cell):
                    if not (self._cps[c2].placed or c2 in do_now):
                        do_next.add(c2)

            # update do_now and clear do_next
            do_now = sorted(do_next, key=lambda c: self._cps[c].nadj, reversed=True)
            do_next.clear()

        # post processing
        pass

    def log(self, s):
        '''Write the given string to the logger'''
        try:
            self._logger.log(s)
        except:
            pass

    def _check_source(self, source):
        '''Check that a source graph is formatted correctly and contains all
        necessary information. A valid source graph must be connected and have
        maximum degree L+2'''

        try:
            assert isinstance(source, nx.Graph), \
                'source must be networkx.Graph or derived'
            assert nx.is_connected(source), \
                'source must be connected'
            assert max(d for n,d in source.degree()) <= self.L+2, \
                'source must have maximum degree of L+2'
            return True
        except AssertionError as e:
            print('Invalid source: {0}'.format(e.message))
            return False

    @staticmethod
    def _check_chimera(chimera):
        '''Attempt to extract the chimera parameters from a given Chimera graph.
        Checks that the Chimera graph is valid. If valid, returns a tuple:
        (mins, maxs, lmax). Otherise returns a tuple of Nones.

        returns:
            mins    : smallest (m,n) values in the Chimera graph
            maxs    : largest (m,n) values in the Chimera graph
            L       : detected number of qubits per half-tile.

        note: L estimation is equal to the largest detected l+1 value so for
        accurate linear indexing at least one qubit must have the maximum
        l value.
        '''

        try:
            # check tup formats
            mins, maxs = [0,0], [0,0]   # range of m,n values
            lmax = -1                   # largest l value
            for node, tup in chimera.nodes('tup'):
                assert tup is not None, 'Missing tup parameter'
                assert isinstance(tup[0], int) and tup[0]>=0, 'Invalid tile m'
                assert isinstance(tup[1], int) and tup[1]>=0, 'Invalid tile n'
                assert isinstance(tup[3], int) and tup[3]>=0, 'Invalid tile l'
                mins = [min(a,b) for a,b in zip(mins, tup[:2])]
                maxs = [max(a,b) for a,b in zip(maxs, tup[:2])]
                lmax = max(lmax, tup[3])
            assert min(mins) >= 0, 'Invalid tile range'
            assert lmax >= 0, 'Invalid l range'
            return mins, maxs, lmax+1
        except AssertionError as e:
            print(e.message)
            return (None,)*3

    @staticmethod
    def _remap_chimera(chimera):
        '''Make a shallow copy of the Chimera graph and change all the names of
        the nodes to the 4-tuple'''

        mapping = lambda x: tuple(chimera.node[x]['tup'])
        return nx.relabel_nodes(chimera, mapping)
