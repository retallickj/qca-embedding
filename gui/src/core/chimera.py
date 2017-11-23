#!/usr/bin/env python

#---------------------------------------------------------
# Name: chimera.py
# Purpose: Functions for handling the chimera graph
# Author:	Jacob Retallick
# Created: 11.30.2015
# Last Modified: 11.30.2015
#---------------------------------------------------------

import numpy as np
from itertools import product

# constants
L = 4   # number of qubits per half tile


def linear_to_tuple(ind, M, N, L=4, index0=False):
    '''Convert the linear index of a qubit in an (N, M, L) processor to
    tuple format'''

    qpr = 2*N*L     # qbits per row
    qpt = 2*L       # qbits per tile

    if not index0:
        ind -= 1

    row, rem = divmod(ind, qpr)
    col, rem = divmod(rem, qpt)
    horiz, ind = divmod(rem, L)

    return (row, col, horiz, ind)


def tuple_to_linear(tup, M, N, L=4, index0=False):
    '''Convert a tuple format index of a qubit in an (N, M, L) processor
    to linear format'''

    qpr = 2*N*L     # qbits per row
    qpt = 2*L       # qbits per tile

    return (0 if index0 else 1) + qpr*tup[0]+qpt*tup[1]+L*tup[2]+tup[3]


def load_chimera_file(filename):
    '''Load a chimera graph from an edge specification file'''

    try:
        fp = open(filename, 'r')
    except:
        print('Failed to open file: {0}'.format(filename))
        raise IOError

    # get number of qubits and number of connectors
    num_qbits, num_conns = [int(x) for x in fp.readline().split()]

    adj = {i: [] for i in range(1, num_qbits+1)}

    for line in fp:
        a, b = [int(x) for x in line.strip().split()]
        adj[a].append(b)
        adj[b].append(a)

    # processor size
    M = int(np.sqrt(num_qbits/(2*L)))
    N = M

    fp.close()

    return M, N, adj


def generate_chimera_adj(M, N, L=4):
    '''Generate a full chimera adjacency dict of size MxN'''

    adj = {}

    # intra-tile connections
    for m, n, h, l in product(range(M), range(N), range(2), range(L)):
        adj[(m, n, h, l)] = [(m,n,1-h,x) for x in range(L)]

    # horizontal inter-tile connections
    for m, n, l in product(range(M), range(N-1), range(L)):
        adj[(m, n, 1, l)].append((m, n+1, 1, l))
        adj[(m, n+1, 1, l)].append((m, n, 1, l))

    # vertical inter-tile connections
    for m, n, l in product(range(M-1), range(N), range(L)):
        adj[(m, n, 0, l)].append((m+1, n, 0, l))
        adj[(m+1, n, 0, l)].append((m, n, 0, l))

    return adj
