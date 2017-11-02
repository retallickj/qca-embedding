#---------------------------------------------------------
# Name: embed.py
# Purpose: Main loop for running Dense Placement embedding.
# Author:	Jacob Retallick
# Created: 06.08.2014
#---------------------------------------------------------

#######################################################################
### IMPORTS ###

from __future__ import division

from bisect import bisect
from random import random, shuffle
from copy import copy as cp
import numpy as np
import sys
import os   # for avoiding file overwriting
import re
import itertools

import routing as Routing


#######################################################################
### GLOBALS ###

# default Chimera parameters
M = 8   # number of tile rows
N = 8	  # number of tile columns
L = 4	  # number of qubits per half tile


### working variables

# source variables
_numAdj = {}        # number of unplaced adjacent cells, source keyed
_numAdj2 = {}       # number of cells at a 2-distance
_source = {}        # adjacency list for cell connectivity

# target variables
_qubits = {}        # source keyed dict of each cell's qubit
_cells = {}         # 4-tup keyed dict of each qbits' cell
_qbit_paths = {}    # 4-tup keyed dict of paths containing each qbit
_qbitAdj = {}       # 4-tup keyed dict of adjacent qbit lists
_vacancy = []       # number of free columns/rows of tiles [L,R,D,U]
_tile_occ = {}      # number of used qbits in each tile row/column, 4-tup keyed
_reserved = {}      # 4-tup keyed dict of reserved qbit lists

# flags
_cell_flags = {}    # source keyed flag dict for each cell
_qbit_flags = {}    # 4-tup keyed flag dict for each qubit
_reserved = {}      # source keyed dict of sets of reserved adjacent qubits

_paths = {}	     # source keyed dict of all paths in the embedding


### handles

FIRST_PROB_POW = 4          # power for firstCell probabilistic method
FIRST_QBIT_SIG = .3         # gaussian dist. tile select standard dev.
FIRST_QBIT_ATTEMPTS = 5     # number of attempts to find starting qbit

# search costs
IN_TILE_COST = 1        # cost for a new qubit in the same tile
OUT_TILE_COST = 1.9     # cost for a new qubit outside the same tile
EDGE_REP_COST = 0.5     # linear cost for tiles closer to the processor edge

# seam ranking weights
SEAM_QBIT_COST = 30      # cost for an assigned qbit to a broken qbit
SEAM_PATH_COST = 40      # cost for breaking a path
SEAM_EXT_COST = 0       # cost for each extended connection
SEAM_DIST_COST = 30      # cost for the distance between the seam and av. qbit

MAX_SEARCH_COUNT = 3    # maximum number of additional times to run the
                        # multisource search algorithm before failure asserted

VERBOSE = False
WRITE = True

WRITE_DIR = '../sols/temp/'
WRITE_PATH = None
#ROUTE_PATH = 'routing' if WRITE else ''
ROUTE_PATH = None
PORT_PATH = '../bin/temp/port'


#######################################################################
#######################################################################
## LOGGING

LOGGING = False
LOG_PATH = '../bin/logs/log'
_fp_log = None


# checked
def writeSol(fp):
    '''write the solution to the given file pointer'''

    global _qubits, _paths, M, N, L, _qbitAdj, _qbit_flags

#    nq = len(_qubits)
#    npq = len(filter(None, _qubits.values()))

    # header
    fp.write('<header>\n')
    fp.write('M = %d\n' % M)
    fp.write('N = %d\n' % N)
    fp.write('L = %d\n' % L)
    fp.write('</header>\n\n')

    # disabled qubit list
    qbs = itertools.product(range(M), range(N), [0, 1], range(L))
    dis_qbits = [qb for qb in qbs if not _qbitAdj[qb]]

    fp.write('<dis_qbits>\n')
    for qb in dis_qbits:
        fp.write('%s\n' % str(qb))
    fp.write('</dis_qbits>\n')

    # write qbit assignments
    fp.write('\n\n<qubits>\n')
    for cell in sorted(_qubits):
        fp.write('%s : %s\n' % (str(cell), str(_qubits[cell])))
    fp.write('</qubits>\n')

    # write paths
    fp.write('\n\n<paths>\n')
    for ends in sorted(_paths):
        path = _paths[ends]
        fp.write('%s : %s\n' % ('; '.join(map(str, ends)),
                                '; '.join(map(str, path))))
    fp.write('</paths>')

    num_qbits = 0
    for qbit in _qbit_flags:
        if _qbit_flags[qbit]['taken']:
            num_qbits += 1

    print 'Used %d qubits...' % num_qbits


# checked
def portSol(ext=None):
    '''write the current solution to file: open and close file so update is
    immediate'''

    log('\n\n'+'*'*40 + '\n')
    log('Porting solution for analysis...\n')

    try:
        fname = PORT_PATH
        if not ext is None:
            fname += str(ext)
        fp = open(fname, 'w')
    except IOError:
        print('Failed to open port file... likely invalid directory path')
        return None

    writeSol(fp)
    fp.close()

    log('\nComplete...\n')
    log('*'*40 + '\n')


# checked
def initLog():
    '''Initialise log file'''
    global _fp_log

    if LOGGING:
        _fp_log = open(LOG_PATH, 'w')


# checked
def killLog():
    '''Close log file'''
    global _fp_log

    if LOGGING:
        _fp_log.close()


# checked
def log(txt):
    '''Log if LOGGING flag'''
    global _fp_log

    if LOGGING:
        if VERBOSE:
            print txt,
        _fp_log.write(txt)


# checked
def logSol(cell_map, paths):
    '''Log and print formatted solution'''
    global _qbit_flags

    if LOGGING:
        log('\n\n'+'*'*30+'\n\n')

        log('Cell Mapping:\n\n')
        for cell in cell_map:
            log('%s: %s\n' % (str(cell), str(cell_map[cell])))

        log('\n\nPaths:\n\n')
        for key in paths:
            c1, c2 = key
            path = paths[key]
            log('%s <> %s : \t %s\n' % (str(c1), str(c2), str(path)))

        # number of qubits used
        nq = 0
        for key in _qbit_flags:
            if _qbit_flags[key]['taken']:
                nq += 1
        log('Number of qubit: %d' % nq)


# checked
def formatSol():
    '''Format the solution for ease of interpretation

    output: paths (dict)    : (cell,cell) keyed dictionary of qubit paths
            cell_map (dict) : cell keyed dictionary of assigned qubits
    '''

    global _paths, _qubits, _cells

    # generate inverse cell map
    cell_map = cp(_qubits)

    # generate formatted paths dictionary
    paths = {}
    for key in _paths:
        path = _paths[key]
        if _qubits[key[0]] != path[0]:
            path = path[::-1]
        paths[key] = path

    return cell_map, paths


#######################################################################
#######################################################################
### FORMATTING and CONVERSION ###


# checked, complete
def indexToLinear(tup, index0=False):
    ''' convert a 4-tuple index to a linear qubit index. Tuple format
    tup=(row,col,horiz?,index) with row,col,index starting from the
    bottom left tile/qubit.
    '''

    global M, N, L

    qpr = 2*N*L     # qbits per row
    qpt = 2*L       # qbits per tile

    return (0 if index0 else 1) + qpr*tup[0]+qpt*tup[1]+L*tup[2]+tup[3]


# checked, complete
def indexToTuple(index, index0=False):
    ''' converts a linear qubit index to a 4-tuple index. '''

    global M, N, L

    qpr = 2*N*L     # qbits per row
    qpt = 2*L       # qbits per tile

    if not index0:
        index -= 1

    row, rem = divmod(index, qpr)
    col, rem = divmod(rem, qpt)
    horiz, ind = divmod(rem, L)

    return (row, col, horiz, ind)


#######################################################################
#######################################################################
### INITIALIZATION METHODS ###

# checked, possibly include pro-processing for wire attraction
def initialize(source):
    '''Initialise embedding solver'''

    global _source, M, N, L

    initLog()

    log('Starting embedding\n\n')

    log('Target Conditions:\n')
    log('M:%d\t N:%d\t L:%d\n\n' % (M, N, L))

    _source = source

    # set trial dependent parameters
    reset()

    log('\n\n')


# checked
def reset():
    '''Reset all trial specific parameters'''

    global _numAdj, _numAdj2, _qbitAdj, _tile_occ, _source
    global _reserved, _cells, _qubits, _qbit_paths, _paths, _vacancy

    log('Generating numAdj...')
    # generate numAdj from source
    _numAdj = {key: len(_source[key]) for key in _source}
    log('complete\n')

    log('Generating numAdj2...')
    # generate numAdj2
    _numAdj2 = {}   # sum of numAdj over each cell adjacent to key
    for key in _source:
        _numAdj2[key] = sum([_numAdj[adj] for adj in _source[key]])
    log('complete\n')

    # initialise flags and reserved dictionaries
    initFlags()

    log('Initializing routing algorithm')
    # configure routing algorithm
    Routing.initialize(_qbitAdj)

    # initialise _tile_occ
    _tile_occ = {}
    _tile_occ['r'] = [0 for _ in xrange(M)]
    _tile_occ['c'] = [0 for _ in xrange(N)]

    # set _vacancy... need to have placed a cell first so set to -1
    _vacancy = [-1, -1, -1, -1]

    _reserved, _cells, _qbit_paths, _paths, _qubits = {}, {}, {}, {}, {}

    # _reserved sets
    for key in _qbitAdj:
        _reserved[key] = set()
        _cells[key] = None
        _qbit_paths[key] = set()

    for key in _source:
        _qubits[key] = None

    _paths = {}

# checked, complete
def setChimera(chimera_adj, m, n, l):
    '''
    updates the Chimera graph size.

    inputs: m (int) : number of tile rows
            n (int) : number of tile columns
            l (int) : number of horizontal or vertical qubits per tile

    outputs: none
    '''

    global _qbitAdj, M, N, L

    M, N, L = m, n, l

    _qbitAdj = chimera_adj

    # sort each keyed list
    for key in _qbitAdj:
        _qbitAdj[key].sort()


# checked
def initFlags():
    '''Initialise *_flags and reserved dicts'''
    global _numAdj, _qbitAdj, _cell_flags, _qbit_flags

    log('Setting flags...')

    # initialise cell_flags

    for key in _numAdj:
        d = {}
        d['placed'] = False     # cell has assigned qubit
        _cell_flags[key] = d

    # initialise qbit_flags
    for key in _qbitAdj:
        d = {}
        d['taken'] = False      # qubit is used for routing
        d['reserved'] = False   # qubit is reserved for adjacency
        d['assigned'] = False   # qubit is assigned to cell
        d['prox'] = set()       # qubit adjacent to assigned qbit
        d['c_in'] = set()       # set of free adjacent internal qubits
        d['c_out'] = 0          # number of free adjacent external qubits
        _qbit_flags[key] = d

    log('complete\n')


# checked, complete
def setTileOcc(qbits, dec=False):
    '''
    '''
    global _tile_occ

    for qbit in qbits:
        if dec:
            _tile_occ['r'][qbit[0]] -= 1
            _tile_occ['c'][qbit[1]] -= 1
        else:
            _tile_occ['r'][qbit[0]] += 1
            _tile_occ['c'][qbit[1]] += 1


# checked, complete
def setVacancy():
    '''update and set _vacancy'''

    global _tile_occ, _vacancy

    # compute left/right vacancy
    occupied = map(lambda x: x > 0, _tile_occ['c'])
    left = occupied.index(True)
    right = occupied[::-1].index(True)

    # compute bottom/top vacancy
    occupied = map(lambda x: x > 0, _tile_occ['r'])
    bot = occupied.index(True)
    top = occupied[::-1].index(True)

    _vacancy = [left, right, bot, top]


# checked, complete
def firstCell(M1=False):
    '''returns the first cell to be placed'''
    global _numAdj, _numAdj2

    log('Selecting first cell\n')

    if len(_numAdj.keys()) == 1:
        log('\tselected cell 0\n\n')
        return 0

    # create adjacency worths for each cell
    worth = {key: (_numAdj[key], _numAdj2[key]) for key in _numAdj}
    # sort by decreasing worth
    order = sorted(_numAdj.keys(), key=lambda x: worth[x])[::-1]

    ### method 1: max adj

    if M1:
        log('\trunning method 1: max adjacency...')
        # determine how many cells have the maximum worth
        num_max = worth.values().count(worth[order[0]])
        # randomly select one of these cells
        i = int(random()*num_max)
        cell = order[i]
        log('done\n')

    ### method 2: fully probabilistic
    # probability is ~ numAdj**POW for some power

    else:
        log('\trunning method 2: probabilistic...')
        # give a probability score for each cell
        probs = {key: pow(worth[key][0], FIRST_PROB_POW) for key in worth}
        # normalise and compute comparison values
        total_prob = sum(probs.values())
        comps = [0]
        for key in order:
            probs[key] /= total_prob
            comps.append(comps[-1]+probs[key])
        # randomly select starting key
        i = max(bisect(comps, random()), 1)

        cell = order[i-1]
        log('done\n')

    log('\tselected cell %d\n\n' % cell)
    return cell


# checked, complete
def firstQubit(cell, M1=False):
    '''Selects the qubit corresponding to the first cell'''
    global _qbitAdj, _numAdj

    log('Selecting first Qubit\n')
    qb = None
    adj = _numAdj[cell]

    ### method 1: middle cell

    if M1:
        log('\trunning method 1: middle tile... ')
        # select candidate tile(s)
        n, m = [N//2], [M//2]
        if N % 2 == 0:
            n.append(N//2-1)
        if M % 2 == 0:
            m.append(M//2-1)
        tiles = [(_n, _m) for _n in n for _m in m]

        # shuffle tiles
        shuffle(tiles)

        for tile in tiles:
            r, c = tile
            # try to find suitable qubit
            order = [(h, i) for h in xrange(2) for i in xrange(L)]
            shuffle(order)
            for h, i in order:
                qbit = (r, c, h, i)
                if len(_qbitAdj[qbit]) >= adj:
                    qb = qbit
                    break
            else:
                continue
            break

        log('done\n')

    ### method 2: Gaussian dist

    else:
        log('\trunning method 2: gaussian dist... ')

        if N % 2:   # if odd rows
            Y = np.arange(-(N//2), N//2+1)
        else:
            Y = .5+np.arange(-(N//2), N//2)

        if M % 2:   # if odd rows
            X = np.arange(-(M//2), M//2+1)
        else:
            X = .5+np.arange(-(M//2), M//2)

        # generate probabilities
        CDF = []
        for ax in [X, Y]:
            Z = np.exp(-ax*ax/(2*FIRST_QBIT_SIG))
            Z /= np.sum(Z)
            cdf = [0.]
            for z in Z:
                cdf.append(cdf[-1]+z)
            CDF.append(cdf)

        # attempt to find qubit
        attempt = 0
        while attempt < FIRST_QBIT_ATTEMPTS:
            attempt += 1
            # pick tile
            r = max(bisect(CDF[0], random()), 1)-1
            c = max(bisect(CDF[1], random()), 1)-1
            # pick qubit
            order = [(h, i) for h in xrange(2) for i in xrange(L)]
            shuffle(order)
            for h, i in order:
                qbit = (r, c, h, i)
                if len(_qbitAdj[qbit]) >= adj:
                    qb = qbit
                    break
            else:
                continue
            break

        log('done\n')

    if qb is None:
        log('\n***Failed to identify a suitable qubit')
        return None

    log('\tselected qbit: %s\n\n' % (str(qb)))
    return qb


#######################################################################
#######################################################################
### UPDATE METHODS ###


# checked
def assignQubit(cell, qbit):
    '''Assign a qubit to QCA cell'''

    global _numAdj, _source, _qubits, _cells
    global _cell_flags, _qbit_flags, _qbitAdj

    # decrement numAdj for each edjacent cell
    for adj in _source[cell]:
        _numAdj[adj] -= 1

    # set _qubits and _cells
    _qubits[cell] = qbit
    _cells[qbit] = cell

    # update tile_occ and vacancy
    if not _qbit_flags[qbit]['reserved']:
        setTileOcc([qbit])
        setVacancy()

    # update flags
    _cell_flags[cell]['placed'] = True

    _qbit_flags[qbit]['taken'] = True
    _qbit_flags[qbit]['assigned'] = True

    # make adjacent qubits aware of place cell (for reserved check)
    for qb in _qbitAdj[qbit]:
        _qbit_flags[qb]['prox'].add(qbit)

    Routing.disableQubits([qbit])


# unchecked
def assignPaths(paths):
    '''Flag and assign routed paths'''

    global _qbit_flags, _qbit_paths, _paths, _cells

    reserve_check = set()

    for path in paths:
        # get end points as key
        key = tuple(map(lambda x: _cells[x], [path[0], path[-1]]))
        for qbit in path:
            # take qbit
            _qbit_flags[qbit]['taken'] = True
            _qbit_paths[qbit].add(key)
            # if qbit is prox, flag for later reserved check
            if _qbit_flags[qbit]['prox']:
                reserve_check.update(_qbit_flags[qbit]['prox'])
        # update tile_occ
        setTileOcc(path[1:-1])
        Routing.disableQubits(path[1:-1])
        _paths[key] = path

    # update vacancy
    setVacancy()

    # check for reservations
    try:
        reserveQubits(list(reserve_check))
    except KeyError:
        raise KeyError('Qubit reservation failed during path assignment')


# unchecked
def reserveQubits(qbits):
    '''for each qbit in qbits, check if adjacent qubits should be
    reserved. Reserve if appropriate.
    '''
    global _cells, _numAdj, _qbit_flags, _reserved

    if not qbits:
        return

    #log('\n\nReserving locals qbits for: %s\n' % map(str, qbits))

    for qbit in qbits:
        #log('\nchecking qbit %s\n' % str(qbit))
        # cell properties
        try:
            cell = _cells[qbit]
            if cell is None:
                raise KeyError
        except KeyError:
            raise KeyError('Qbit %s is not assigned to a cell...'
                           % str(qbit))
        num_adj = _numAdj[cell]

        #log('Required adjacency: %d\n' % num_adj)

        # wipe old reservations
        old_res = set()
        if _reserved[qbit]:
            #log('releasing qbits: %s: \n' % map(str, _reserved[qbit]))
            for qb in _reserved[qbit]:
                _qbit_flags[qb]['reserved'] = False
            setTileOcc(_reserved[qbit], dec=True)
            old_res = cp(_reserved[qbit])
            _reserved[qbit].clear()

        # get list of all adjacent unreserved qubits and count qbit type
        qbs = []
        _qbit_flags[qbit]['c_in'] = set()
        _qbit_flags[qbit]['c_out'] = 0
        for q in _qbitAdj[qbit]:
            if not (_qbit_flags[q]['taken'] or _qbit_flags[q]['reserved']):
                qbs.append(q)
                if q[0:2] == qbit[0:2]:
                    _qbit_flags[qbit]['c_in'].add(q)
                else:
                    _qbit_flags[qbit]['c_out'] += 1

        # if exact amount of adjacent qubits available, reserve
        if num_adj == len(qbs):
            # log('Reserving all free qbits for qbit %s\n' % str(qbit))
            # reserve all adjacent qubits
            res_check = set()
            for qb in qbs:
                _qbit_flags[qb]['reserved'] = True
                _reserved[qbit].add(qb)
                res_check.update(_qbit_flags[qb]['prox'])
            setTileOcc(qbs)
            # if reserved qubits changed, check local qubits for reservations
            if old_res == _reserved[qbit]:
                reserveQubits(res_check-set(qbits))

        # check for insufficient qubits
        elif num_adj > len(qbs):
            raise KeyError('Insufficent free qubits for cell %s' % str(cell))

    setVacancy()


# unchecked
def forgetQubit(qbit, check=True):
    ''' release a given qubit. Returns a list of all connected paths which
    should be forgotten before continuing with the embedding
    '''

    global _cells, _cell_flags, _numAdj, _source, _qubits
    global _qbit_flags, _reserved, _qbit_paths, _qbitAdj

    try:
        cell = _cells[qbit]
    except KeyError:
        log('Qbit has not been assigned to any cell')
        raise KeyError('Qbit has not been assigned to any cell')

    qbs = set([qbit])    # set of qbits to decrement from _tile_occ

    _cells[qbit] = None
    _qubits[cell] = None

    # update flags
    _cell_flags[cell]['placed'] = False

    _qbit_flags[qbit]['assigned'] = False
    _qbit_flags[qbit]['taken'] = False

    # refresh source parameters
    for adj in _source[cell]:
        _numAdj[adj] += 1

    # clear reserved list
    for qb in _reserved[qbit]:
        _qbit_flags[qb]['reserved'] = False
    qbs.update(_reserved[qbit])
    _reserved[qbit].clear()

    # get list of paths connected to qbit
    paths = cp(_qbit_paths[qbit])
    _qbit_paths[qbit].clear()

    # update _tile_occ
    setTileOcc(qbs, dec=True)

    if check:
        setVacancy()

    for qb in _qbitAdj[qbit]:
        _qbit_flags[qb]['prox'].remove(qbit)

    Routing.enableQubits([qbit])

    return paths


# unchecked
def forgetPath(key, check=True):
    '''Free up qubits of the path with the given key and update appropriate
    flags. If check==True, also check qubit reservations for nearby qubits.
    '''
    global _paths, _qbit_flags, _qbit_paths

    reserve_check = set()
    path = cp(_paths[key])
    _paths.pop(key)

    if key in _qbit_paths[path[0]]:
        _qbit_paths[path[0]].remove(key)
    if key in _qbit_paths[path[-1]]:
        _qbit_paths[path[-1]].remove(key)

    for qbit in path[1:-1]:
        _qbit_flags[qbit]['taken'] = False
        _qbit_paths[qbit].clear()
        if _qbit_flags[qbit]['prox']:
            reserve_check.update(_qbit_flags[qbit]['prox'])
    Routing.enableQubits(path[1:-1])
    setTileOcc(path[1:-1], dec=True)
    if check:
        reserveQubits(reserve_check)
        setVacancy()


#######################################################################
#######################################################################
### MULTI-SOURCE SEARCH ###


# checked, modify cost scheme if necessary
def extend_Dijkstra(src):
    '''Generator for Dijkstra search extension'''
    global _qbitAdj, _qbit_flags, M, N

    BIG_VAL = 2*len(_qbitAdj)   # large value for initial node cost

    # initialise
    visited = {}
    for qbit in _qbitAdj:
        if _qbit_flags[qbit]['taken'] or _qbit_flags[qbit]['reserved']:
            visited[qbit] = True
        else:
            visited[qbit] = False

    costs = {qbit: BIG_VAL for qbit in _qbitAdj}
    next_qb = set()
    next_qb.add(src)
    costs[src] = 0

    # tree growth loop
    while next_qb:

        # pick lowest cost qbit and yield
        qbit = sorted(next_qb, key=lambda x: costs[x])[0]
        next_qb.remove(qbit)
        yield qbit

        # mark as visited
        visited[qbit] = True

        # update costs of all unvisited adjacent nodes
        for qb in _qbitAdj[qbit]:
            if not visited[qb]:
                # add cost increment
                dcost = IN_TILE_COST if qb[0:2] == qbit[0:2] else OUT_TILE_COST
                # include edge repulsion
                dcost += EDGE_REP_COST*max(map(abs, [qb[0]-.5*(M-1),
                                                     qb[1]-.5*(N-1)]))
                costs[qb] = min(costs[qb], costs[qbit]+dcost)
                next_qb.add(qb)


# unchecked...implement later
def extend_Astar():
    '''
    '''

    pass


# checked, complete
def multiSourceSearch(srcs, adj, forb=set(), typ='Dijkstra'):
    '''
    Attempt to find the lowest cost suitable point with free paths to
    the given sources
    '''

    global _qbitAdj, _reserved

    # create path extension generator

    if typ.upper() == 'DIJKSTRA':
        extend_func = extend_Dijkstra
    else:
        extend_func = extend_Astar

    extend = {}
    # initialise generator for each source
    for src in srcs:
        # release local reserved qbits (only for current src)
        for qb in _reserved[src]:
            _qbit_flags[qb]['reserved'] = False

        # create generator
        extend[src] = extend_func(src)
        next(extend[src])   # burn src qbit and initialise

        # reset local reserved abits (forbid for other srcs)
        for qb in _reserved[src]:
                _qbit_flags[qb]['reserved'] = True

    # set visit counts for each qbit
    visits = {qbit: 0 for qbit in _qbitAdj}

    # search loop

    while True:

        cands = []  # candidate nodes

        ## extend each source tree

        for src in srcs:

            # extend
            try:
                node = next(extend[src])
            except StopIteration:     # break if no more nodes to visit
                return None

            # increment visited node count
            visits[node] += 1

            # if node visited from all sources add as candidate
            if visits[node] == len(srcs):
                if node in forb:
                    continue
                cands.append(node)

        ## check for suitable candidate

        # sort by suitability
        cands = sorted(map(lambda x: [suitability(x, srcs), x], cands))[::-1]

        # filter out unsuitable qbits
        cands = filter(lambda x: x[0] >= adj, cands)

        # select qbit
        if cands:
            return cands    # return all candidate qbit


#######################################################################
#######################################################################
### EMBEDDING SUB-ALGORITHMS ###

# checked
def checkSol():
    '''Check that embedding solution is valid'''
    global _qubits, _paths, _source, _qbitAdj

    log('\n\nChecking Solution...\n\n')

    check = False

    if not all(map(lambda x: not x is None, _qubits.values())):
        print _paths
        print _qubits
        raise KeyError('Not all cells were assigned a qubit')

    for c1 in _source:
        for c2 in _source[c1]:
            if not ((c1, c2) in _paths or (c2, c1) in _paths):
                print('No path between %s and %s' % (str(c1), str(c2)))
                check = True
    if check:
        raise KeyError('Not all paths were placed')

    # check all path connections are available and count uses of non cell qbits
    uses = {q: 0 for q in _qbitAdj}
    for path in _paths.values():
        for i in xrange(len(path)-1):
            q1, q2 = path[i: i+2]
            if not q2 in _qbitAdj[q1]:
                print('No coupler available for %s to %s' % (str(q1), str(q2)))
                check = True
            if i > 0:
                uses[q1] += 1

    if check:
        raise KeyError('Not all couplers available')

    for q in uses:
        if uses[q] > 1:
            print('Qubit %s used in %d paths' % (str(q), uses[q]))

    if any(map(lambda x: x > 1, uses.values())):
        raise KeyError('Qubit shared among multiple paths')


# always changing... tentatively done
def suitability(qbit, srcs=[]):
    '''Determine the effective number of free adjacent qubits. Effected by the
    number of mutual free qubits for in-tile assigned qubits'''
    global _reserved, _qbit_flags, _numAdj, _cells, L

    s = 0
    res_flags = {src: 0 for src in srcs}

    c_in = set()
    for qb in _qbitAdj[qbit]:
        # check for special cases for end points
        for src in srcs:
            if qb in _reserved[src]:
                res_flags[src] = 1
        if qb in srcs:
            s += 1
        # general bulk condition
        elif not (_qbit_flags[qb]['taken'] or _qbit_flags[qb]['reserved']):
            s += 1
            # note as free internal qbit if so
            if qb[0:2] == qbit[0:2]:
                c_in.add(qb)

    # account for negotiating free qubits with other in-tile qubits
    r, c, h, l0 = qbit
    qbs = [(r, c, h, l) for l in xrange(L) if l != l0]
    qbs = filter(lambda x: x in _cells and _qbit_flags[x]['assigned'], qbs)

    for qb in qbs:
        cell = _cells[qb]
        s -= max(0, _numAdj[cell] - _qbit_flags[qb]['c_out'] -
                 len(_qbit_flags[qb]['c_in']-c_in))

    return s+sum(res_flags.values())


def placeCell(cell):
    '''Attempt to find a suitable qbit to place input cell on.

    inputs: cell(int)		: source index of cell to place

    output: qbit(tuple)		: 4-tup for qbit to assign to cell
            paths(list)		: list of paths from placed cell qubits to
                              qbit
    '''

    global _qbitAdj, _source, _qubits
    global _cell_flags, _qbit_flags, _reserved, _vacancy

    log('\n'+'#'*30+'\n')
    log('Placing cell: %s\n' % str(cell))

    ### Initialise

    seam_flag = False
    qbit = None

    # find qubits for placed adjacent cells
    adj_qbits = [_qubits[c] for c in _source[cell] if _cell_flags[c]['placed']]
    log('Adjacent qbits: %s\n' % str(adj_qbits))

    # find required availability of target qbit
    avb = len(_source[cell])
    log('Required availability: %d\n' % avb)

    # multisourcesearch parameters
    forb = set()        # list of forbidden qbit for multisourcesearch
    search_count = 0    # counter for number of failed searches

    # every time a seam is opened, we should check to see if there is a
    # better qubit to consider
    while qbit is None:

        ### Open Seam

        if seam_flag:
            log('Running Seam Opening Routine\n')
            seam_flag = False

            # check for vacancies
            if not any(_vacancy):
                log('No vacant columns/rows to open\n\n')
                raise KeyError('Out of room')

            # find available seams
            seams = availableSeams(adj_qbits)

            # analyse seams
            seam_dicts = map(lambda s: genSeamDict(s, adj_qbits), seams)
            seam_dicts = filter(None, seam_dicts)

            if len(seam_dicts) == 0:
                log('No suitable seams detected\n')
                return None, []

            # select seam to open
            seam_dict = selectSeam(seam_dicts)
            log('current vacancy: %s\n' % str(_vacancy))
            log('selected seam %s :: %s\n' % (str(seam_dict['sm']),
                                              str(seam_dict['dr'])))

            # open seam
            success = openSeam(**seam_dict)

            if not success:
                log('Failed to open seam\n\n')
                return None, []

            # update adjacent qubits
            log('Seam successfully opened... \n')
            adj_qbits = [_qubits[c] for c in _source[cell]
                         if _cell_flags[c]['placed']]
            log('New adjacent qbits: %s\n' % str(adj_qbits))

        ### Pick qubit to assign

        # run multisource search method, get list of candidate qbits
        qbits = multiSourceSearch(adj_qbits, avb, forb=forb)

        # check if found
        if not qbits:
            log('multiSourceSearch failed\n')
            seam_flag = True
            continue

        log('Found %d candidate qubits: %s \n'
            % (len(qbits), map(lambda x: str(x[1]), qbits)))

        # check each candidate qbit from multisourcesearch in order
        for qbit in qbits:

            suit, qbit = qbit

            log('Trying qbit: %s with suitability %s ...'
                % (str(qbit), str(suit)))

            ### Find paths

            routes = [[qb, qbit] for qb in adj_qbits]
            end_points = list(set([it for rt in routes for it in rt]))
            # find best consistent paths
            cost = Routing.Routing(routes, _reserved, writePath=ROUTE_PATH)

            # check successful routing

            if cost >= Routing.COST_BREAK:
                log('routing failed...\n')
                # disable end points
                Routing.disableQubits(end_points)
                continue

            log('\t success\n')
            break
        else:
            log('No suitable qbit found\n')
            qbit = None
            search_count += 1
            if search_count >= MAX_SEARCH_COUNT:
                seam_flag = True
                search_count = 0
                forb.clear()
            else:
                forb.update(map(lambda x: x[1], qbits))

    # get paths
    paths = cp(Routing.getPaths().values())

    # disable path qubits
    qbs = list(set([it for path in paths for it in path]))
    Routing.disableQubits(qbs)

    log('Placed on qubit: %s\n\n' % str(qbit))
    
    return qbit, paths


### SEAM OPENING ###

## seam format: seam= [(vert,ind),dir]
#	horz	-> flag for horizontal or vertical seam
# 	ind		-> index of seam from bottom-left
# 	dir		-> flag for opening direction; True opens away from b-l


# unchecked
def availableSeams(qbits):
    ''' returns a list of possible seams to split given a list of qbits. Seams
    are candidates if they are next to one of the qubits and can be'''

    global _vacancy, N, M, _tile_occ

    seams = set()

    for qbit in qbits:
        # left/right seams

        if _vacancy[0] > 0 and qbit[1] > _vacancy[0]:     # open left
            if qbit[1] > 1:
                seams.add(((1, qbit[1]), False))
            seams.add(((1, qbit[1]+1), False))

        if _vacancy[1] > 0 and qbit[1] < ((N-1)-_vacancy[1]):  # open right
            seams.add(((1, qbit[1]), True))
            if qbit[1] < (N-2):
                seams.add(((1, qbit[1]+1), True))

        # down/up seams

        if _vacancy[2] > 0 and qbit[0] > _vacancy[2]: 	# open down
            if qbit[0] > 1:
                seams.add(((0, qbit[0]), False))
            seams.add(((0, qbit[0]+1), False))

        if _vacancy[3] > 0 and qbit[0] < ((M-1)-_vacancy[3]):  # open up
            seams.add(((0, qbit[0]), True))
            if qbit[0] < (M-2):
                seams.add(((0, qbit[0]+1), True))

    return seams


def prepSeam(seam):
    '''prepare parameters for seam opening'''

    global _qbit_flags,  _paths

    sm, dr = seam

    # local functions

    def check_qb(qb):
        '''check if qb lies on the mobile side of seam'''

        if dr:  # increasing index
            return qb[sm[0]] >= sm[1]
        else:   # decreasing index
            return qb[sm[0]] < sm[1]

    # list of qubits to be moved
    qbits = []
    for qb in filter(check_qb, _qbit_flags):
        if _qbit_flags[qb]['assigned']:
                qbits.append(qb)

    # dict of connectors to be 'moved' for each path
    paths = {}
    num_halfs = {}

    for key in _paths:
        path = _paths[key]
        d = []
        nf, nh, nn = 0, 0, 0
        for i in xrange(1, len(path)):
            connect = [path[i-1], path[i]]
            n1, n2 = map(check_qb, connect)
            if n1 and n2:       # completely on mobile side of seam
                d.append([connect, 'full'])
                nf += 1
            elif n1 or n2:      # cut by seam
                d.append([connect, 'half', n1])
                nh += 1
            else:               # not affected by seam
                d.append([connect, 'none'])
                nn += 1
        if nf+nh > 0:
            paths[key] = d
            num_halfs[key] = nh

    return qbits, paths, num_halfs


# unchecked
def genSeamDict(seam, adj_qbits):
    '''generate a dict of seam opening parameters'''

    global _qbit_flags, _paths, _qbitAdj, _source, _cells, _qbit_paths

    ## local functions

    sm, dr = seam

    def target_qb(qb):
        '''returns the target qubit for a given qubit under the specified
        seam and direction'''

        q = list(qb)
        if dr:
            q[sm[0]] += 1
        else:
            q[sm[0]] -= 1
        return tuple(q)

    qbits, paths, num_halfs = prepSeam(seam)
    if len(qbits) + len(paths) == 0:
        return None

    ## check for qbit conflicts
    qbit_confs = []
    for qb in qbits:
        tg = target_qb(qb)
        # throw exception unless target qbit exists and is suitable
        try:
            assert len(_qbitAdj[tg]) >= len(_source[_cells[qb]])
        except KeyError:
            qbit_confs.append(qb)

    ## check for path conflicts
    map_conn = lambda conn: map(target_qb, conn)
    path_confs = set()
    for key in paths:
        full_conns = filter(lambda x: x[1] == 'full', paths[key])
        new_conns = map(lambda x: map_conn(x[0]), full_conns)
        for conn in new_conns:
            if not conn[1] in _qbitAdj[conn[0]]:
                path_confs.add(key)
                break

    # clean up qbits and connects
    for qb in qbit_confs:
        qbits.remove(qb)
    for key in path_confs:
        paths.pop(key)

    # give score based on nearest seam to mean qbit
    mean_qbit = reduce(lambda x, y: (x[0]+y[0], x[1]+y[1]), adj_qbits)
    mean_qbit = map(lambda x: x/len(adj_qbits), mean_qbit)

    seam_dist = abs(mean_qbit[sm[0]]-sm[1]-.5)

    ### compute cost

    # number of qbit and paths to repair
    nq, np = map(len, [qbit_confs, path_confs])
    # number of paths to extend
    ne = 0
    for key in paths:
        ne += num_halfs[key]

    cost = nq*SEAM_QBIT_COST+np*SEAM_PATH_COST+ne*SEAM_EXT_COST
    cost += seam_dist*SEAM_DIST_COST

    seam_dict = {'sm': sm,
                 'dr': dr,
                 'cost': cost,
                 'qbits': qbits,
                 'paths': paths,
                 'qb_conf': qbit_confs,
                 'pt_conf': path_confs,
                 'tg_fn': target_qb,
                 'par': [nq, np, ne, seam_dist]}

    return seam_dict


def moveQbit(qbit, tg_fn):
    ''' Move a qubit under seam opening. Handles forgetting and
    replacing the qubit'''

    global _cells, _qubits

    old_qb, new_qb = qbit, tg_fn(qbit)
    cell = _cells[old_qb]

    log('qbit: %s \t -> \t %s\n' % (str(old_qb), str(new_qb)))
    # forget old qubit and flags,  old paths sohuld be handled elsewhere
    forgetQubit(old_qb)

    # map to new qubit
    assignQubit(cell, new_qb)


def newPath(key, path, tg_fn):
    '''create new path'''

    new_path = []
    # print path
    for i in xrange(len(path)):
        conn, typ = path[i][0:2]
        new_conn = map(tg_fn, conn)
        if typ == 'full':
            if i == 0:
                new_path.append(new_conn[0])
            new_path.append(new_conn[1])
        elif typ == 'half':
            if path[i][2]:
                if i == 0:
                    new_path.append(new_conn[0])
                new_path.append(conn[0])
                new_path.append(conn[1])
            else:
                if i == 0:
                    new_path.append(conn[0])
                new_path.append(conn[1])
                new_path.append(new_conn[1])
        elif typ == 'none':
            if i == 0:
                new_path.append(conn[0])
            new_path.append(conn[1])
        else:
            log('Invalid coupler type in movePath\n')
            raise ValueError('Invalid coupler')

    return new_path


def movePath(key, path, tg_fn):
    ''' Move a path under seam opening. Handles forgetting and replacing'''

    global _paths

    forgetPath(key)

    new_path = []
    log('Moving path: %s \n' % str(key))
    # print path
    for i in xrange(len(path)):
        conn, typ = path[i][0:2]
        new_conn = map(tg_fn, conn)
        if typ == 'full':
            if i == 0:
                new_path.append(new_conn[0])
            new_path.append(new_conn[1])
        elif typ == 'half':
            if path[i][2]:
                if i == 0:
                    new_path.append(new_conn[0])
                new_path.append(conn[0])
                new_path.append(conn[1])
            else:
                if i == 0:
                    new_path.append(conn[0])
                new_path.append(conn[1])
                new_path.append(new_conn[1])
        elif typ == 'none':
            if i == 0:
                new_path.append(conn[0])
            new_path.append(conn[1])
        else:
            log('Invalid coupler type in movePath\n')
            raise ValueError('Invalid coupler')
    # print new_path
    assignPaths([new_path])


# unchecked
def openSeam(sm, dr, cost, qbits, paths, qb_conf, pt_conf, tg_fn, par):
    ''' recursively open seam'''

    global _cells, _numAdj, M, N, L, _qbitAdj, _cell_flags, _qubits, _reserved
    ## erase conflicts

    log('Broken qubits: \n %s \n\n' % str(qb_conf))

    # erase qbit conflicts and update broken paths
    log('Erasing %d broken qubits...\n' % len(qb_conf))
    cell_conf = map(lambda x: _cells[x], qb_conf)
    for qb in qb_conf:
        new_paths = forgetQubit(qb)
        log('%s: %s\n' % (str(qb), str(new_paths)))
        pt_conf.update(new_paths)
        for pt in new_paths:
            if pt in paths:
                paths.pop(pt)
    log(' done\n')

    # erase path conflicts
    log('Erasing %d broken paths...\n' % len(pt_conf))
    for path in pt_conf:
        log('path: %s\n' % str(path))
        forgetPath(path)
    log(' done\n')

    # retain path keys between good qbits
    pt_rp = []
    for key in pt_conf:
        if _cell_flags[key[0]]['placed'] and _cell_flags[key[1]]['placed']:
            pt_rp.append(key)

    ## STORE OLD VALUES AND WIPE SEAM TARGET REGION

    # map of new qbits for each cell
    qbit_dict = {_cells[qb]: tg_fn(qb) for qb in qbits}

    # map of new paths
    path_dict = {key: newPath(key, paths[key], tg_fn) for key in paths}

    # wipe old qubits and paths
    log('Wiping old values...\n')
    log('\t qbits...')
    for qb in qbits:
        forgetQubit(qb, False)
    log('done\n')
    log('\t paths...')
    for key in paths:
        forgetPath(key, False)
    log('done\n')

    # assign new qubits and paths
    log('Assigning new values...\n')
    for cell in qbit_dict:
        assignQubit(cell, qbit_dict[cell])

    assignPaths(path_dict.values())

    # update all reserved qubits
    log('updating all reservations...\n')
    reserveQubits(filter(None, _qubits.values()))

    ## repair broken paths

    # only place paths between moved qubits
    log('Attempting to replace broken routes between moved qubits\n')
    routes = []
    for pt in pt_rp:
        rt = map(lambda x: _qubits[x], pt)
        routes.append(rt)
    cost = Routing.Routing(routes, _reserved,  writePath=ROUTE_PATH)

    # check successful routing
    if cost >= Routing.COST_BREAK:
        log('routing failed...\n')
        raise KeyError('Routing failed for paths between moved qbits in \
        seam opening... fix code later')

    # get paths
    log('Routing successfull, preparing to assign new paths\n')
    fixed_paths = cp(Routing.getPaths().values())

    # disable path qubits
    log('Disabling new paths\n')
    qbs = list(set([it for path in fixed_paths for it in path]))
    Routing.disableQubits(qbs)

    # assign paths and update reservations
    log('Assigning new paths \n')
    assignPaths(fixed_paths)

    ## repair qbit placements, should automatically deal with paths

    # order qbits to be placed by decreasing adjacency
    log('Preparing new qbit placements \n')
    cell_ord = sorted(cell_conf, key=lambda x: -_numAdj[x])

    for cell in cell_ord:
        log('cell: %s ...' % str(cell))
        new_qb, new_paths = placeCell(cell)     # recursive call

        # abort on failed placement
        if new_qb is None:
            return False
            raise KeyError('Failed cell %s placement in seam \
            opening' % str(cell))

        # assign qubit and paths
        assignQubit(cell, new_qb)
        assignPaths(new_paths)

        # handle reservations
        reserveQubits([new_qb])

    return True


def selectSeam(seam_dicts):
    '''select which seam to open based on minimum cost. If multiple
    seams have the same cost randomly select one'''

    seams = sorted(seam_dicts, key=lambda x: x['cost'])

    log('Seam candidates:\n')
    for seam in seams:
        log('%s :: %5s   cost = %.1f \t %s\n'
            % (str(seam['sm']), str(seam['dr']),
               seam['cost'], str(seam['par'])))

    cands = filter(lambda x: x['cost'] == seams[0]['cost'], seams)

    return cands[int(random()*len(cands))]




# NEW METHODS FOR WIRE SHORTENING
def shorten_wire_paths():
    '''find the wires. Find a qbit path for each wire that most resemble's
    the length of the wire path. To do this, find the shortest possible path
    and lengthen the path until it is an appropriate length'''

    global _paths
    global _qubits
    global _qbit_flags

    wires, nodes = get_wires()
    wires = sorted(wires, key=lambda x: len(x))

    all_paths = dict(_paths)
    all_qbs = dict(_qubits)

    pre_ex = sum(map(lambda x: len(x)-2, _paths.values()))

    print "There are %d cells, there are %d in qbit flags" %(len(_cells), len(_qbit_flags))
    # clear qbit flags except the nodes (which do not move)
    qb_nodes = []
    for c in nodes:
        qb_nodes.append(_qubits[c])
    for qb in _cells:
        _qbit_flags[qb]['taken'] = qb in qb_nodes

    update = True

    qb_wires = []
    elongate = []

    for wire in wires:

        # get the shortest path for the wire
        qb_wire = get_wire_path(_qubits[wire[0]], _qubits[wire[-1]], qb_nodes)

        # set flags, add qb_wire to correct list
        # qb_wires are for wires that are already the best length
        # elongate are for wires that need to be lengthened
        if qb_wire:
            for qb in qb_wire:
                _qbit_flags[qb]['taken'] = True

            if len(qb_wire) < len(wire):
                elongate.append((wire, qb_wire))
            else:
                qb_wires.append((wire, qb_wire))

        else:
            for qb in wire:
                if qb in nodes:
                    continue
                _qbit_flags[all_qbs[qb]]['taken'] = False

    # lengthen wires taht need lenthening
    for w, qb_w in elongate:
        qb_w = lengthen(qb_w, w, qb_nodes, 1)
        if qb_w is None:
            update = False
        qb_wires.append((w, qb_w))

    # updatek paths and qbs
    if update:
        for w, qb_w in qb_wires:
            update_mapping(qb_w, w, nodes, all_paths, all_qbs)

    # ensure correct flags
    used_qb = []
    for path in all_paths.values():
        used_qb.extend(path)

    for qb in _qbit_flags:
        if _qbit_flags[qb]['taken']:
            if qb not in used_qb:
                _qbit_flags[qb]['taken'] = False

    post_ex = sum(map(lambda x: len(x)-2, all_paths.values()))

    # if improvements are found, save them
    if post_ex < pre_ex and update:
        _paths = all_paths
        _qubits = all_qbs
        return (pre_ex - post_ex)
    # otherwise restore older correct mappings
    else:
        print "Removing long paths failed, using earlier embedding"
        used_qb = []
        for path in _paths.values():
            used_qb.extend(path)
        for qb in _qbit_flags:
            _qbit_flags[qb]['taken'] = qb in used_qb
        return 0


def get_wire_path(head, end_qb, nodes):
    '''get the best qbit path for the qca circuit wire'''
    visited = []
    complete, paths = extend_path([head], end_qb, nodes, visited)
    if complete:
        return paths

    # get new paths until one reaches the end (simple djikstra's)
    while paths:
        curr_path = paths.pop(0)
        complete, new_paths = extend_path(curr_path, end_qb, nodes, visited)
        if complete:
            return new_paths

        added = False
        for i in xrange(len(paths)):
            if new_paths and len(paths[i]) > len(new_paths[0]):
                added = True
                for p in new_paths:
                    paths.insert(i, p)
                break

        if not added:
            for p in new_paths:
                paths.append(p)


def extend_path(curr_path, end_qb, nodes, visited):
    '''given a current path, returns the a list of the possible next paths'''

    head = curr_path[-1]
    paths = []
    next_qbs = _qbitAdj[head]
    arrived = False
    # append the curr path with all possibilities
    for qb in next_qbs:
        if arrived:
            if not _qbit_flags[qb]['taken']:
                return (True, curr_path + [qb])
        elif qb == end_qb:
            if qb in nodes or not _qbit_flags[qb]['taken']:
                return (True, curr_path + [qb])
            else:
                arrived = True
        elif not _qbit_flags[qb]['taken'] and \
            qb not in curr_path and qb not in visited:
            visited.append(qb)
            paths.append(curr_path + [qb])

    return (False, paths)


def lengthen(qb_wire, wire, nodes, count):
    '''increase the length of a path'''

    global _qbit_flags

    n = len(wire) - len(qb_wire)
    i = 0
    # method to increase path length by 4
    while n >= 3 and i+1 < len(qb_wire):
        qb = qb_wire[i]

        # if the next qb is perpendicular move a tile over,
        # switch orientation twice, then move back to the original tile
        if qb_wire[i+1][2] != qb[2]:

            if qb[0] + (qb[2]^1) >= 8 or qb[1] + qb[2] >= 8:
                i += 1
                continue

            qb1 = (qb[0] + (qb[2]^1), qb[1] + qb[2], qb[2], qb[3])

            if _qbit_flags[qb1]['taken']:

                if qb[0] + (-1)*(qb[2]^1) < 0 or qb[1] + (-1)*qb[2] < 0:
                    i += 1
                    continue

                qb1 = (qb[0] + (-1)*(qb[2]^1), qb[1] + (-1)*qb[2], qb[2], qb[3])

                if _qbit_flags[qb1]['taken']:
                    i += 1
                    continue

            found = False
            for j in xrange(L):
                qb2 = (qb1[0], qb1[1], qb1[2]^1, j)
                if not _qbit_flags[qb2]['taken']:
                    found = True
                    break

            if not found:
                i += 1
                continue

            found = False
            for j in xrange(L):
                qb3 = (qb1[0], qb1[1], qb1[2], j)
                qb4 = (qb[0], qb[1], qb[2], j)
                if not _qbit_flags[qb3]['taken'] and qb3 != qb1 \
                    and not _qbit_flags[qb4]['taken'] and qb4 != qb:

                    found = True
                    break

            if not found:
                i += 1
                continue

            qb_wire.insert(i + 1, qb4)
            qb_wire.insert(i + 1, qb3)
            qb_wire.insert(i + 1, qb2)
            qb_wire.insert(i + 1, qb1)

            _qbit_flags[qb1]['taken'] = True
            _qbit_flags[qb2]['taken'] = True
            _qbit_flags[qb3]['taken'] = True
            _qbit_flags[qb4]['taken'] = True

            n -= 4

        # if the next qb is parallel switch orientation twice,
        # move to the correct tile, switch orientation again
        elif qb_wire[i+1][2] == qb[2]:
            found = False
            for j in xrange(L):
                qb1 = (qb[0], qb[1], qb[2]^1, j)
                if not _qbit_flags[qb1]['taken']:
                    found = True
                    break

            if not found:
                i += 1
                continue

            found = False
            for j in xrange(L):
                qb2 = (qb[0], qb[1], qb[2], j)

                if qb_wire[i+1][0] == qb[0]:
                    m = qb_wire[i+1][1] - qb[1]
                else:
                    m = qb_wire[i+1][0] - qb[0]

                if qb[0] + (m)*(qb[2]^1) < 0 or qb[1] + (m)*qb[2] < 0 \
                    or qb[0] + (m)*(qb[2]^1) >= 8 or qb[1] + (m)*qb[2] >= 8:
                    i += 1
                    continue
                qb3 = (qb[0] + (m)*(qb[2]^1), qb[1] + (m)*qb[2], qb[2], j)

                if not _qbit_flags[qb2]['taken'] and \
                    not _qbit_flags[qb3]['taken'] and qb2 != qb:

                    found = True
                    break

            if not found:
                i += 1
                continue

            found = False
            for j in xrange(L):
                qb4 = (qb3[0], qb3[1], qb3[2]^1, j)
                if not _qbit_flags[qb4]['taken']:
                    found = True
                    break

            if not found:
                i += 1
                continue

            qb_wire.insert(i + 1, qb4)
            qb_wire.insert(i + 1, qb3)
            qb_wire.insert(i + 1, qb2)
            qb_wire.insert(i + 1, qb1)

            _qbit_flags[qb1]['taken'] = True
            _qbit_flags[qb2]['taken'] = True
            _qbit_flags[qb3]['taken'] = True
            _qbit_flags[qb4]['taken'] = True

            n -= 4

    i = 0
    # method to increase path length by 2
    while n > 0 and i+1 < len(qb_wire):
        qb = qb_wire[i]

        # must be in the same tile
        # switch orientation twice
        if qb_wire[i+1][2] != qb[2]:
            found = False
            for j in xrange(L):
                qb1 = (qb[0], qb[1], qb[2]^1, j)
                if not _qbit_flags[qb1]['taken']:
                    found = True
                    break

            if not found:
                i += 1
                continue

            found = False
            for j in xrange(L):
                qb2 = (qb[0], qb[1], qb[2], j)
                if not _qbit_flags[qb2]['taken'] and qb2 != qb:
                    found = True
                    break

            if not found:
                i += 1
                continue

            qb_wire.insert(i + 1, qb2)
            qb_wire.insert(i + 1, qb1)

            _qbit_flags[qb1]['taken'] = True
            _qbit_flags[qb2]['taken'] = True

            n -= 2

        elif qb_wire[i+1][2] == qb[2]:
            i += 1

    # if unable to lengthen enough, create a new different shortest path,
    # and try lengthening again
    if n > 0:
        holder_qb = []
        if count > len(qb_wire):
            return None
        for j in xrange(0,count):
            if qb_wire[j] in nodes:
                continue
            holder_qb.append(qb_wire[j])
        for qb in qb_wire[count:]:
            if qb_wire[j] in nodes:
                continue
            _qbit_flags[qb]['taken'] = False

        new_qb_wire = get_wire_path(qb_wire[0], qb_wire[-1], nodes)
        if not new_qb_wire:
            return None

        for qb in holder_qb:
            _qbit_flags[qb]['taken'] = False
        for qb in new_qb_wire:
            _qbit_flags[qb]['taken'] = True

        qb_wire = lengthen(new_qb_wire, wire, nodes, count + 1)

    return qb_wire


def update_mapping(qb_wire, wire, nodes, all_paths, all_qbs):
    '''update paths and qbs'''

    # remove old unused paths
    for i in xrange(1, len(wire)):
        if (wire[i-1], wire[i]) in all_paths:
            del all_paths[(wire[i-1], wire[i])]
        elif (wire[i], wire[i-1]) in all_paths:
            del all_paths[(wire[i], wire[i-1])]
        else:
            print "ERROR: Can't find path in paths"

    # add new paths and qbs
    if len(qb_wire) > len(wire):
        for i in xrange(len(wire)):
            all_qbs[wire[i]] = qb_wire[i]
            if i > 0 and i < len(wire) - 1:
                all_paths[(wire[i-1],wire[i])] = [qb_wire[i-1], qb_wire[i]]

        all_paths[(wire[len(wire) - 2], wire[len(wire) - 1])] = qb_wire[len(wire) - 2:]
    else:
        for i in xrange(len(qb_wire)):
            all_qbs[wire[i]] = qb_wire[i]
            if i > 0 and i < len(wire) - 1:
                all_paths[(wire[i-1],wire[i])] = [qb_wire[i-1], qb_wire[i]]

        all_paths[(wire[len(wire) - 2], wire[len(wire) - 1])] = qb_wire[len(wire) - 2:]


def get_wires():
    '''get the wires (defined as segments where every cell only has 2
    connections) of a qca circuit'''

    wires = []
    cells = _qubits.keys()
    nodes = []

    # find the first node (cell with more than 2 connections) and start there
    curr_c = None
    for cell in cells:
        if len(_source[cell]) > 2:
            curr_c = cell
            cells.remove(cell)
            break

    # for all cells left
    next_c = None
    curr_wire = []
    while cells:
        # add to current wire
        curr_wire.append(curr_c)

        # if we reach a node, save current wire, start a new one
        adj_cells = _source[curr_c]
        if len(adj_cells) > 2:
            nodes.append(curr_c)
            if len(curr_wire) > 1:
                wires.append(curr_wire)
                curr_wire = [curr_c]

        # find the next cell connected to this cell
        next_c = None
        for cell in adj_cells:
            if cell in cells:
                next_c = cell
                break

        # if next cell is already in a wire, add the end node to the wire
        if next_c is None:
            for cell in adj_cells:
                if cell in nodes and cell not in curr_wire:
                    curr_wire.append(cell)

        # find a new node to start a new wire from
        if next_c is None:
            for node in nodes:
                for pcell in _source[node]:
                    if pcell in cells:
                        wires.append(curr_wire)
                        curr_wire = [node]
                        next_c = pcell
                        break
                if next_c:
                    break

        curr_c = next_c
        cells.remove(next_c)

    curr_wire.append(curr_c)
    wires.append(curr_wire)
    return wires, nodes


#######################################################################
#######################################################################
### MAIN ###


# unchecked
def denseEmbed(source, write=False):
    '''
    Attempts to find an embedding of the source graph into a global
    target Chimera graph.

    inputs:	source(dict)	: adjacency dict: source graph

    outputs: qubits (dict)	: source node indexed mapping of assigned
                             qubits.
             routes (dict)	: (node1,node2) indexed dictionary of qubit
                             routes; node1 < node2
             info (dict)	: dsescribe later ...
    '''

    global _source, _cell_flags, _numAdj, _qbitAdj

    ### INITIALIZE ###

    initialize(source)

    ### INITIAL SEED ###

    # select first cell
    cell = firstCell()

    # select first qubit
    qbit = firstQubit(cell)

    if qbit is None:
        #print 'No suitable first qubit found'
        sys.exit()

    # take first qubit
    assignQubit(cell, qbit)

    # handle reservations
    reserveQubits([qbit])

    # update do* lists
    doNow = sorted(source[cell], key=lambda x: -_numAdj[x])
    doNext = set()

    ### GENERAL PLACEMENT LOOP ###

    # if doNow is non-empty, there are still cells to place
    while doNow:

        log('\n\n')
        log('*'*50 + '\n')
        log('*'*50 + '\n')
        log('*'*50 + '\n')
        log('toDo: %s\n' % str(doNow))
        # place each cell in doNow

        for cell in doNow:

            # find qbit and paths from placed cells
            qbit, paths = placeCell(cell)

            # abort on failed placement
            if qbit is None:
                print 'No placement of cell %s found' % str(cell)
                raise

            # assign qubit and paths
            assignQubit(cell, qbit)
            assignPaths(paths)

            # handle reservations
            reserveQubits([qbit])

            # add unplaced adjacent cells to doNext
            for c2 in source[cell]:
                if not (_cell_flags[c2]['placed'] or c2 in doNow):
                    doNext.add(c2)

        # update doNow and clear doNext

        doNow = sorted(doNext, key=lambda x: -_numAdj[x])
        doNext.clear()

    # post processing path shortening
#    shorten_wire_paths()

    checkSol()
    log('\n\n***Embedding complete\n\n')
    cell_map, paths = formatSol()

    logSol(cell_map, paths)
    killLog()

    if WRITE and write:
        print 'writing solution',
        try:
            if not (WRITE_PATH is None):
                fname = WRITE_PATH
                fp = open(fname, 'w')
            elif WRITE_DIR:
                regex = re.compile('^sol[0-9]+$')   # default name format
                old_sols = filter(regex.match, os.listdir(WRITE_DIR))
                old_ext = map(lambda x: int(x[3::]), old_sols)  # 3 for 'sol'
                old_max = max(old_ext) if len(old_ext) > 0 else -1
                fname = WRITE_DIR + 'sol%d' % (old_max+1)
            else:
                fname = None
                raise IOError('No valid file destination')

            if not fname is None:
                fp = open(fname, 'w')
                writeSol(fp)
                fp.close()
        except IOError as e:
            print e.message
            print 'Invalid filename: %s' % fname

    return cell_map, paths
