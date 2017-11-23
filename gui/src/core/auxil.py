#!/usr/bin/env python

#---------------------------------------------------------
# Name: auxil.py
# Purpose: Auxiliary function for use in embedder application
# Author: Jacob Retallick
# Created: 2015.11.26
# Last Modified: 2015.11.26
#---------------------------------------------------------

import numpy as np
from copy import deepcopy as dc
from heapq import heappop, heappush
from collections import defaultdict

from itertools import combinations

# physical parameters
eps0 = 8.85412e-12  # permittivity of free space
epsr = 12.          # relative permittivity
q0 = 1.602e-19      # elementary charge

CELL_FUNCTIONS = {'QCAD_CELL_NORMAL': 0,
                  'QCAD_CELL_INPUT': 1,
                  'QCAD_CELL_OUTPUT': 2,
                  'QCAD_CELL_FIXED': 3}

CELL_MODES = {'QCAD_CELL_MODE_NORMAL': 0,
              'QCAD_CELL_MODE_CROSSOVER': 1,
              'QCAD_CELL_MODE_VERTICAL': 2,
              'QCAD_CELL_MODE_CLUSTER': 3}

R_MAX = 1.8             # maximum interaction range
STRONG_THRESH = 0.3     # threshold for qualifying as a strong interaction
D_ERR = 0.2             # allowed error in DX, DY equality


### GENERAL FUNCTIONS


def pinch(string, pre, post):
    '''selects the string between the first instance of substring (pre) and
    the last instance of substring (post)'''
    return string.partition(pre)[2].rpartition(post)[0]


def zeq(x, y, err):
    '''returns True if x and y differ by at most err'''
    return abs(x-y) < err


def gen_pols(n):
    '''Generate all possible polarizations for n cells'''
    if n <= 0:
        return []
    return [tuple(2*int(x)-1 for x in format(i, '#0{0}b'.format(n+2))[2:])
            for i in range(pow(2, n))]


### QCA CELL PROCESSING

def getOutputs(cells):
    return [i for i,cell in enumerate(cells) if cell['cf'] == CELL_FUNCTIONS['QCAD_CELL_OUTPUT']]

def getEk(c1, c2, DR=2):
    '''Compute the kink energy using the qdot positions (in nm) of two cells.
    If the cell displacement is greater than DR return False'''

    # check cell-cell range
    dx = c1['x']-c2['x']
    dy = c1['y']-c2['y']
    if dx*dx+dy*dy > DR*DR:
        return False

    qdots_1 = c1['qdots']
    qdots_2 = c2['qdots']

    # compute displacements

    x1 = [qd['x'] for qd in qdots_1]
    y1 = [qd['y'] for qd in qdots_1]

    x2 = [qd['x'] for qd in qdots_2]
    y2 = [qd['y'] for qd in qdots_2]

    X1 = np.array([x1, y1]).T.reshape([4, 1, 2])
    X2 = np.array([x2, y2]).T.reshape([1, 4, 2])

    R = np.sqrt(np.sum(pow(X1-X2, 2), axis=2))

    if np.min(R) == 0:
        print('qdot overlap detected')
        return 0.

    # QCADesigner orders qdots either CW or CCW so same and diff configurations
    # are always alternating indices.

    Q = np.array([1, -1, 1, -1])    # template for charge arrangement

    Q = np.outer(Q, Q)

    Ek = -1e9*q0*np.sum((Q/R))/(8*np.pi*eps0*epsr)

    return Ek

def prepare_convert_adj(cells, spacing, J):
    '''Prepares useful variables for converting from the parse_qca J matrix to
    a reduced adjacency form.

    outputs:    Js  : J scaled by the nearest neighbour interaction of two
                     non-rotated cells.
                T   : Array of cell-cell types for each element of J
                        1  -> non-rotated - non-rotated
                        0  -> non-rotated - rotated
                        -1 -> rotated - rotated
                A   : Enumeration of interaction types
                        3  -> adapter
                        0  -> T=0
                        1  -> nearest neighbour
                        -1 -> diagonal
                        2  -> dx,dy = (2,0) or (0,2)
                        -2 -> dx,dy = (2,1) or (1,2)
                DX  : X displacements in grid-spacings
                DY  : Y displacements in grid-spacings
    '''

    # scale J by the kink energy of two non-rotated adjacent cells
    E_nn = np.max(np.abs(J))
    Js = np.round(J/E_nn, 4)

    # determine interaction type of each element of J:
    #   1  -> non-rotated - non-rotated
    #   0  -> non-rotated - rotated
    #   -1 -> rotated - rotated

    rot = [cell['rot'] for cell in cells]   # array of rotated flags
    rot = 1*np.array(rot).reshape([-1, 1])

    T = 1-(rot+rot.T)
    #T = T.astype(int)

    # get displacements between each cell

    X = np.array([cell['x'] for cell in cells]).reshape([-1, 1])
    Y = np.array([cell['y'] for cell in cells]).reshape([-1, 1])

    DX = (X.T - X)/spacing
    DY = (Y - Y.T)/spacing

    #

    N = len(cells)
    A = np.zeros([len(cells), len(cells)], dtype=int)

    for i in range(N-1):
        for j in range(i+1, N):
            dx, dy = abs(DX[i, j]), abs(DY[i, j])
            # first check for adapter condition
            if zeq(min(dx, dy), 0.5, D_ERR) and zeq(max(dx, dy), 1, D_ERR):
                A[i, j] = 3
            elif T[i, j] == 0:
                continue
            # (1, 0) or (0, 1) interaction
            elif zeq(dx+dy, 1, D_ERR):
                A[i, j] = 1
            # (1, 1) interaction
            elif zeq(dx, 1, D_ERR) and zeq(dy, 1, D_ERR):
                A[i, j] = -1
            elif zeq(min(dx, dy), 0, D_ERR) and zeq(max(dx, dy), 2, D_ERR):
                A[i, j] = 2
            elif zeq(min(dx, dy), 1, D_ERR) and zeq(max(dx, dy), 2, D_ERR):
                A[i, j] = -2

    A += A.T

    return Js, T, A, DX, DY


def new_identify_inverts(A):

    # an inverting interaction is one between two cells A = -1 that share
    # no A=1 neighbours

    # list of A=1 neighbours for each cell
    Ns = [set(a.nonzero()[0].tolist()) for a in A==1]

    # set overlap
    S = np.zeros(A.shape, dtype=bool)
    for i,j in combinations(range(A.shape[0]),2):
        if Ns[i].intersection(Ns[j]):
            S[i,j] = S[j,i] = True

    # loop through all inverter candidates (will handle (i,j) and (j,i))
    is_inv = (S==False)*(A==-1)

    return is_inv


def identify_inverters(A):
    '''Identify inverter cells'''

    # an inverter is a cell with at most one strong interaction and two weak
    # interactions, each having one strong interaction
    invs = {}
    N = A.shape[0]

    # count number of strong interactions for each cell
    num_strong = [np.count_nonzero(A[i,:] == 1) for i in range(N)]

    for i in range(N):
        if num_strong[i] <= 1:   # check if an inverter
            adj = [j for j in range(N) if A[i, j] == -1 and num_strong[j] == 1]
            if len(adj) == 2:
                for k in range(N):
                    if A[adj[0], k] == 1 and A[adj[1], k] == 1:
                        break
                else:
                    invs[i] = adj

    return invs


def identify_xovers(A):
    '''Identify all rotated crossover cells in a circuit'''

    # xover condition: two cells have A=2 and no path of A:1,1
    N = A.shape[0]
    cands = [(i,j) for i in range(N-1) for j in range(i+1, N) if A[i,j] == 2]

    xovers = []
    for cand in cands:
        for k in range(N):
            if abs(A[cand[0], k]) == 1 and abs(A[cand[1], k]) == 1:
                break
        else:
            xovers.append(cand)

    return xovers


def convert_to_full_adjacency(J, Js, T, A, DX, DY):
    '''Convert the J matrix from parse_qca to include only full adjacency
    interactions'''

    xovers = identify_xovers(A)
    N = A.shape[0]

    F = np.ones([N, N])

    for i in range(N-1):
        for j in range(i+1, N):
            if (i, j) in xovers or A[i, j] in [1, -1, 3]:
                continue
            else:
                F[i, j] = 0
                F[j, i] = 0

    return J*F


def convert_to_lim_adjacency(J, Js, T, A, DX, DY):
    '''Convert the J matrix from parse_qca to include only limited adjacency
    interactions'''

    # start with full adjacency representation
    Js = np.array(convert_to_full_adjacency(J, Js, T, A, DX, DY))

    # get inverters and their included diagonal interactions
    #invs = identify_inverters(A)
    is_inv = new_identify_inverts(A)

    # clear all diagonal interactions for non-inverters
    N = A.shape[0]
    for i in range(N-1):
        for j in range(i+1, N):
            if A[i, j] != -1 or is_inv[i,j]:
                continue
            else:
                Js[i, j] = Js[j, i] = 0.

    return J*(Js != 0)

def qca_to_coef(cells, spacing, J, adj=None):
    '''construct h and J matrix from parse_qca parameters. J matrix includes
    all cell interactions, including driver and fixed'''

    # convert J matrix for appropriate adjacency

    if adj is not None:
        Js, T, A, DX, DY = prepare_convert_adj(cells, spacing, J)
        if adj=='full':
            J_ = convert_to_full_adjacency(J, Js, T, A, DX, DY)
        elif adj=='lim':
            J_ = convert_to_lim_adjacency(J, Js, T, A, DX, DY)
        else:
            J_ = dc(J)
    else:
        J_ = dc(J)

    # separate indices of driver/fixed cells from normal/output cells
    dinds, inds = [], []
    CFs = [CELL_FUNCTIONS['QCAD_CELL_NORMAL'], CELL_FUNCTIONS['QCAD_CELL_OUTPUT']]
    for i in range(len(cells)):
        if cells[i]['cf'] in CFs:
            inds.append(i)
        else:
            dinds.append(i)

    # isolate active J matrix
    Js = J_[inds, :][:, inds]

    # interaction matrix between driver/fixed cells and normal/output cells
    Jx = J_[inds, :][:, dinds]

    # compute h parameters
    try:
        pols = np.array([cells[i]['pol'] for i in dinds])
    except:
        print('One of the driver/fixed cells has no specified polarization')
        pols = np.zeros([len(dinds),], dtype=float)

    h = np.dot(Jx, pols).reshape([-1,])

    return h, Js

def bfs_order(x, A):
    '''Determine modified bfs order for degrees x and adjacency matrix A'''
    N = len(x)
    adj = [np.nonzero(A[i])[0] for i in range(N)]
    dx = 1./(N*N)
    C = defaultdict(int)
    pq = [min([(y, i) for i, y in enumerate(x)]),]
    visited = [False]*N
    order = []
    while len(order)<N:
        y, n = heappop(pq)
        if visited[n]:
            continue
        order.append(n)
        visited[n] = True
        for n2 in [a for a in adj[n] if not visited[a]]:
            x[n2] -= C[y]*dx
            heappush(pq, (x[n2], n2))
        C[y] += 1
    return order

def hash_mat(M):
    '''Return a hash value for a numpy array'''
    M.flags.writeable = False
    val = hash(M.data)
    M.flags.writeable = True
    return val

def hash_problem(h, J, res=3):
    '''Generate a hash for a given (h,J) pair. Decimal resolution of h and J
    should be specified by res argument.'''

    # regularize formating
    h = np.array(h).reshape([-1,])
    J = np.array(J)

    # assert sizes
    N = h.size
    assert J.shape == (N, N), 'J and h are of different sizes'

    # normalise
    K = np.max(np.abs(J))
    h /= K
    J /= K

    # reduce resolution to avoid numerical dependence of alpha-centrality
    h = np.round(h, res)
    J = np.round(J, res)

    # compute alpha centrality coefficients for positive parity
    evals = np.linalg.eigvalsh(J)
    lmax = max(np.abs(evals))
    alpha = .99/lmax

    cent = np.linalg.solve(np.eye(N)-alpha*J, h)

    # enumerate centrality coefficients up to given percent resolution
    val, enum = min(cent), 0
    for c, i in sorted([(c,i) for i,c in enumerate(cent)]):
        f = (c-val)/abs(val)
        if f > 1e-10:
            val, enum = c, enum+1
        cent[i] = enum

    # candidate parities
    s = np.sign(np.sum(h))
    hps = [+1, -1] if s == 0 else [s]

    hash_vals = {}
    inds = {}

    for hp in hps:
        inds[hp] = bfs_order(cent*hp, J != 0)
        h_ = ((10**res)*h[inds]*hp).astype(int)
        J_ = ((10**res)*J[inds,:][:, inds]).astype(int)
        hash_vals[hp] = hash((hash_mat(h_), hash_mat(J_)))

    return hash_vals, K, hps, inds
