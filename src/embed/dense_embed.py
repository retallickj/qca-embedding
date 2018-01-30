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

from ..core.utility import dget
import .logger as lg
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

class CellPars:
    '''All relevant parameters for a given cell during embedding'''
    qbit = None     # qbit assigned to the cell, None if no qbit assigned
    placed = False  # True if cell has been placed
    nadj = -1       # current number of unplaced adjacent cells
    def __init__(self, nadj=None):
        if isinstance(nadj, int) and nadj >= 0:
            self.nadj = nadj

class QubitPars:
    '''All relevant parameters for a given qubit during embedding'''
    cell = None         # cell the qubit is assigned to, else None
    reserved = False    # qbit is reserved for pending adjacency
    taken = False       # qbit taken for routing
    prox = set()        # set of assigned adjacent qubits
    c_in = set()        # set of free internal qubits
    c_out = -1          # number of free external qubits
    paths = set()       # paths containing the qubit

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
    FIRST_PROB_POW = 4      # power for firstCell method
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
        if lmax:
            print('Invalid chimera format')
            return

        # store chimera and properties
        self.chimera = chimera
        self.origin = mins                              # top-left tile
        self.span = [1+b-a for a,b in zip(mins, maxs)]  # number of columns/rows
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

        # handle default arguments
        verbose = dget(kwargs, 'verbose', False, mp=bool)
        ntrials = dget(kwargs, 'ntrials', 1, mp=int, chk=lambda n: n>0)
        best = dget(kwargs, 'best', False, mp=bool)

        # validate optional arguments
        self.vlog = self.log if verbose else lambda *a, **k: None

        # initialize solver for source graph
        try:
            self._initialise(source):
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

        self.vlog('Initialising embedder for source graph...')

        # check source
        if not self._check_source(source):
            raise SourceError('Invalid source graph')
        self.source = source

        self.vlog('done')


    def _reset(self):
        '''Reset the Embedder for a new trial'''

        # set up cell parameters
        self._cps = {cell: CellPars(deg) for cell, deg in self.source.degree()}

        # set up qubit parameters
        # NOTE: may have to use tuple keys rather than the qbit labels
        self._qps = {qbit: QubitPars() for qbit in self.chimera}

        # paths between cell qubits
        self._paths = {}


    def _first_cell(self):
        '''Pick the first cell to place'''

        self.log('Selecting first cell...')
        keys = self._cps.keys()
        worth = lambda key: pow(self._cps[key].nadj, self.FIRST_PROB_POW)
        probs = np.array([worth(key) for key in keys])
        cell = np.random.choice(keys, p=probs)
        self.log('done')


    def _first_qbit(self, cell):
        '''Find the qubit seed for the first cell'''

        self.log('Attempting to select first cell')


        self.log('done')

    def _assign_qbit(self, cell, qbit):
        '''Associate the given qubit and cell'''
        pass

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


    def _embed_trial(self):
        '''Run a single embedding trial'''

        # Reset the embedder for the new trial
        self._reset()

        # attempt to find initial seed
        for init_trial in range(self.INIT_TRIALS):
            cell = self._first_cell()       # seed cell
            qbit = self._first_qbit(cell)   # seed qbit
            if qbit:
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

    @staticmethod
    def _check_source(source):
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
                assert isinstance(tup[0], int) and tup[0]>=0, 'Invalid tile m'
                assert isinstance(tup[1], int) and tup[1]>=0, 'Invalid tile n'
                assert isinstance(tup[3], int) and tup[3]>=0, 'Invalid tile l'
                mins = [min(a,b) for a,b in zip(mins, tup[:2])]
                maxs = [max(a,b) for a,b in zip(maxs, tup[:2])]
                lmax = max(lmax, tup[3])
            assert min(mins) >= 0, 'Invalid tile range'
            assert lmax >= 0, 'Invalid l range'
            return mins, maxs, lmax+1
        except AssertionError:
            return (None,)*3
