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
__date__        = '2018-02-01'      # last update

from core.utility import dget, range_product
from collections import defaultdict
import networkx as nx
import numpy as np
from heapq import heappush, heappop
from copy import copy as cp
from itertools import chain

from .router import Router

# exceptions

class SeamError(Exception):
    '''Failure due to seam opening'''
    pass

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
    reserve = set()     # set of qbits reserved by this qbit
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

        self._router = Router()

        # setup for all possible extender algorithms
        self.extenders = {'dijkstra': self._extend_dijkstra}

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


    # internal methods

    def _initialise(self, source):
        '''Prepare the Embedder for the given source adjacency graph'''

        # start the Logger
        self._logger.start()

        self.vlog('Initialising embedder for source graph... ')

        # check source
        if not self._check_source(source):
            raise SourceError('Invalid source graph')
        self.source = source

        # set up router
        self._router.initialise(source)

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


    # seed

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


    # handling

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
        self._router.disable_nodes([qbit])

    def _assign_paths(self, paths):
        '''Set all relevant flags for the given routed paths'''

        reserve_check = set()   # set of qbits to check later for reserves
        for path in paths:
            key = (self._qps[qb].cell for qb in [path[0], path[-1]])
            for qbit in path:
                self._qps[qbit].taken = True
                self._qps[qbit].paths.add(key)
                # if qbit has any prox, flag for later reserve check
                if self._qps[qbit].prox:
                    reserve_check |= self._qps[qbit].prox
            self._update_tile_occ(path[1:-1])
            self._router.disable_nodes(path[1:-1])
            self._paths[key] = path

        self._set_vacancy()
        self._reserve_qubits(reserve_check)

    def _pop_reserve(self, qbit):
        '''Clears the reserve of the given qubit and unsets each reserve qbits
        reserved flag. Returns a copy of the reserve'''

        reserve = self._qps[qbit].reserve.copy()
        if reserve:
            for qb in reserve:
                self._qps[qb].reserved = False
            self._update_tile_occ(reserve, dec=True)
            self._qps[qbit].reserve.clear()
        return reserve

    def _poll_neighbors(self, qbit):
        '''Get a list of unreserved neighbors for the given qbit and update the
        qbits internal and external counts'''

        unreserved = []
        self._qps[qbit].c_in = set()
        self._qps[qbit].c_out = 0
        for qb in self.chimera.neighbors(qbit):
            if not (self._qps[qb].taken or self._qps[qb].reserved):
                unreserved.append(qb)
                if qb[:2] == qbit[:2]:
                    self._qps[qbit].c_in.add(qb)
                else:
                    self._qps[qbit].c_out += 1
        return unreserved

    def _push_reserve(self, qbit, reserve):
        '''Adds a reserve to the given qubit. Returns a set of qbits which
        may now need an updated reserve'''
        res_check = set()
        self._qps[qbit].reserve = reserve.copy()
        for qb in reserve:
            self._qps[qb].reserved = True
            res_check |= self._qps[qb].prox
        self._update_tile_occ(reserve)
        return res_check

    def _reserve_qbits(self, qbits):
        '''For each of a set of qubits, check if the adjacent qubits need to
        be reserved in order to satisfy connectivity constraints. Reserve as
        appropriate.'''

        if not qbits:
            return

        for qbit in qbits:
            cell = self._qps[qbit].cell
            # get number of adjacent unplaced cells for associated cell
            if cell is None:
                raise KeyError('Qubit {0} not assigned to a cell'.format(qbit))
            nadj = self._cps[cell].nadj

            old_reserve = self._pop_reserve(qbit)   # previous reserved neighbors
            unreserved = self._poll_neighbors(qbit) # new unreserved neighbors

            # reserve only if there are only just enough unreserved qbits
            if nadj == len(unreserved):
                res_check = self._push_reserve(qbit, unreserved)
                #NOTE this might actually be != rather than ==
                if old_reserve == unreserved:
                    self._reserve_qbits(res_check - set(qbits))
            elif nadj > len(reserved):
                raise PlacementError('No free qubits for {0}'.format(cell))

        self._set_vacancy()

    def _forget_qbit(self, qbit, check=True):
        '''Release a given qubit. Returns a list of all connected paths which
        should be forgotten before continuing with the embedding'''

        cell = self._qps[qbit].cell
        if cell is None:
            raise KeyError('Qubit {0} not assigned to a cell'.format(qbit))

        self._update_tile_occ([qbit], dec=True)

        self._qps[qbit].cell = None
        self._qps[qbit].assigned = False
        self._qps[qbit].taken = False

        self._cps[cell].qbit = None
        self._cps[cell].placed = False

        for c2 in self._source.neighbor(cell):
            self._cps[c2] += 1

        self._pop_reserve(qbit)

        if check:
            self._set_vacancy()

        for qb in self.chimera.neighbors(qbit):
            self._qps[qb].prox.remove(qbit)

        self._router.enable_nodes([qbit])

        paths = self._qps[qbit].paths.copy()
        self._qps[qbit].clear()
        return paths

    def _forget_path(self, key, check=True):
        '''Free up qubits of the path with the given key and update appropriate
        flags. If check, also check qubit reservations for nearby qubits.'''

        path = self._paths.pop(key)

        # remove path from the end nodes
        for qb in [path[0], path[-1]]:
            if key in self._qps[qb].paths:
                self._qps[qb].paths.remove(key)

        res_check = set()
        for qbit in path[1:-1]:
            self._qps[qbit].taken = False
            self._qps[qbit].paths.clear()
            if self._qps[qbit].prox:
                res_check |= self._qps[qbit].prox
        self._router.enable_nodes(path[1:-1])
        self._update_tile_occ(path[1:-1], dec=True)
        if check:
            self._reserve_qbits(res_check)


    # multi-source search

    def _extend_dijkstra(self, src):
        '''Generator for Dijkstra search extension'''

        visited = {}
        for qb in self._chimera:
            visited[qb] = self._qps[qb].taken or self._qps[qb].reserved
        costs = {qb: 2*len(self._chimera) for qb in self._chimera}
        queue = [src]
        costs[src]=0

        # search loop
        while queue:
            # pop and yield the lowest cost node in the heap
            qbit = heappop(queue)
            if visited[qbit]:
                continue
            yield qbit
            visited[qbit] = True

            # update neighbor costs
            for qb in self._chimera[qbit]:
                if not visited[qb]:
                    dcost = self.IN_TILE_COST if qb[:2]==qbit[:2] else self.OUT_TILE_COST
                    dcost += self.EDGE_REP_COST*max(map(abs,[qb[0]-.5*(self.M-1),
                                                             qb[1]-.5*(self.N-1)]))
                    costs[qb] = min(costs[qb], costs[qbit]+dcost)
                    heappush(queue, qb)

    def _setup_extend_gen(self, src, extender):
        '''Setup the extension generator from the given source qubit'''

        for qb in self._qps[src].reserve:
            self._qps[qb].reserved = False

        gen = extender(src)
        next(gen)   # burn src qbit and initialise

        for qb in self._qps[src].reserve:
            self._qps[qb].reserved = True

        return gen

    def _next_mss_candidates(self, extend, visits, forb):
        '''Get the next set of MSS candidates: those reached by all extenders'''

        cands = []
        for src, ext in extend:
            try:
                node = next(ext)
            except StopIteration:
                raise KeyError('No MSS solution found...')
            visits[node] += 1
            if visits[node] == len(extend) and node not in forb:
                cands.append(node)
        mcands = map(lambda x: [self._suitability(x, extend.keys()), x], cands)
        cands = sorted(mcands, reverse=True)
        cands = list(filter(lambda x: x[0] >= adj, cands))
        return cands

    def _multi_source_search(self, srcs, adj, forb=set(), typ='dijkstra'):
        '''Attempt to find the lowest cost suitable qubits with free paths to
        the given sources

        inputs:
            srcs    : iterable of sources to search from
            adj     : required number of free neighbors at target
            forb    : set of forbidden nodes to exclude from search
            typ     : type of search expansion, Dijkstra search by default

        return:
            cands   : iterable of all candidate qbits
        '''

        try:
            extender = self.extenders[typ]
        except:
            raise PlacementError('Invalid MSS extender type')

        # setup extender from each src
        extend = {src: self._setup_extend_gen(src, extender) for src in srcs}
        visits = {qb: 0 for qb in self._chimera}

        # main search loop
        while True:
            try:
                cands = self._next_mss_candidates(extend, visits, forb)
                if cands:
                    return cands
            except KeyError:
                return None     # search finished and failed

    # cell placement

    def _suitability(self, qbit, srcs=[]):
        '''Determine the effective number of free adjacent qubits. Affected
        by the number of mutual free qubits for in-tile assigned qubits'''

        # for each src reserve that intersects the qbit's neighbourhood: 1 pt
        nbh = set(self._chimera[qbit])
        suit = sum(bool(nbh & self._qps[src].reserve) for src in srcs)

        for qb in nbh:
            if qb in srcs:
                suit += 1
            elif not (self._qps[qb].taken or self._qps[qb].reserved):
                suit += 1
                if qb[:2] == qbit[:2]:
                    c_in.add(qb)

        # account for negotitating free qubits with other in-tile qubits
        m,n,h,l = qbit
        qbs = [(m,n,h,_l) for _l in range(self.L) if _l != l]
        qbs = filter(lambda x: x in self._qps and self._qps[x.assigned], qbs)

        for qb in qbs:
            cell = self._qps[qb].cell
            suit -= max(0, self._qps[qb].nadj - self._qps[qb].c_out -
                            len(self._qps[qb].c_in - c_in))

        return suit

    def _suitable_qbit(self, adj_qbs, avl, forb, reserves):
        '''Attempt to find a suitable qbit '''
        pass

    def _place_cell(self, cell):
        '''Attempt to find a placement for the given cell.

        returns:
            qbit    : label of qubit to assign to cell
            paths   : list of paths from previously placed neighbours to the
                      assigned qubit
        '''

        self.log('Placing cell: {0}\n'.format(cell))

        seam_flag = False   # open seam next loop iteration
        qbit = None         # assigned cell, loop escape

        # qubits used for placed adjacent cells
        adj_qbs = [self._cps[c].qbit for c in self._source[cell]
                            if self._cps[c].placed]
        avl = len(self._source[cell])
        reserves = {qb: pars.reserve for qb,pars in self._qps.items()}

        # multisource search parameters
        forb = set()        # list of forbidden qbits for MSS
        search_count = 0    #  counter for number of failed searches

        # placement loop
        while qbit is None:

            if seam_flag:
                seam_flag = False
                adj_qbs = self._open_seam() # adj_qbs can change on seam opening

            # attempt to find suitable qbit and handle result
            qbit = self._suitable_qbit(adj_qbs, avl, forb, reserves)
            if qbit is None:
                search_count += 1
                if search_count >= self.MAX_SEARCH_COUNT:
                    seam_flag, search_count = True, 0
                    forb.clear()
                else:
                    forb.update(map(lambda x: x[1], qbits))

        paths = cp(self._router.paths.values())
        self._router.disable_nodes(set(chain(*paths)))

        self.log('Placed on qubit: {0}\n\n'.format(qbit))

        return qbit, paths



            # get list of candidate qbits
            qbits = self._multi_source_search(adj_qbits, avl, forb=forb)
            if not qbits:
                seam_flag = True
                continue
            self.log('Found {0} candidate qubits: {1}\n'.format(
                        len(qbits), ','.join(qb[1] for qb in qbits)))

            for suit, qb in qbits:
                if self._route(qb, adj_qbs, reserves):
                    qbit = qb
                    break
            else:
                self.log('No suitable qbit found\n')
                search_count += 1
                if search_count >= self.MAX_SEARCH_COUNT:
                    seam_flag = True
                    search_count = 0
                    forb.clear()
                else:
                    forb.update(map(lambda x: x[1], qbits))


    # seam opening

    def _available_seams(self, qbits):
        '''Identify possible seams to split given a list of qbits to be routed.
        Seams are candidates if they arenext to one of the qubits and can be
        opened'''

        seams = set()

        # NOTE: could clean this up a bit... 4 similar blocks
        for qbit in qbits:

            # open left
            if self._vacancy[0] > 0 and qbit[1] > self._vacancy[0]:
                if qbit[1] > 1:
                    seams.add(((1, qbit[1]), False))
                seams.add(((1, qbit[1]+1), False))

            # open right
            if self._vacancy[1] > 0 and qbit[1] < ((self.N-1)-self._vacancy[1]):
                seams.add(((1, qbit[1]), True))
                if qbit[1] < (self.N-2):
                    seams.add(((1, qbit[1]+1), True))

            # open down
            if self._vacancy[2] > 0 and qbit[0] > self._vacancy[2]:
                if qbit[0] > 1:
                    seams.add(((0, qbit[0]), False))
                seams.add(((0, qbit[0]+1), False))

            # open up
            if self._vacancy[3] > 0 and qbit[0] < ((self.M-1)-self._vacancy[3]):
                seams.add(((0, qbit[0]), True))
                if qbit[0] < (self.M-2):
                    seams.add(((0, qbit[0]+1), True))

        return seams



    def _open_seam(self, adj_qbits):
        '''Identify and attempt to open a seam to allow for new routes'''

        self.log('Running Seam Opening Routine\n')

        if not any(self._vacancy):
            self.log('No vacant columns/rows to open\n\n')
            raise SeamError('Out of room')

        seams = self._available_seams(adj_qbits)

        # analyze seams
        seam_dicts = map(lambda s: self._gen_seam_dict(s, adj_qbits), seams)
        seam_dicts = filter(None, seam_dicts)

        if len(seam_dicts) == 0:
            self.log('No suitable seams detected\n')
            raise SeamError('No suitable seams')

        # select seam to open
        seam_dict = self._select_seam(seam_dicts)
        self.log('current vacancy: {0}\n'.format(str(self._vacancy)))
        self.log('selected seam {0} :: {1}\n'.format(
                    str(seam_dict['sm']), str(seam_dict['dr'])))

        if not self._split_seam(**seam_dict):
            self.log('Failed to open seam\n\n')
            raise SeamError('Failed to split seam')
        self.log('Seam successfully opened...\n')
        # update adjacent qubits
        adj_qbits = [self._cps[c].qbit for c in self._cps if self_cps[c].placed]
        self.log('New adjacent qbits: {0}\n'.format(str(adj_qbits)))

        return adj_qbits



    # other, for now

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
