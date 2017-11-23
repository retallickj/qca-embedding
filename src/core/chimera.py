#!/usr/bin/env python
# encoding: utf-8

'''
Chimera hardware connectivity graph structure and related methods
'''

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2017-10-31'      # last update


from collections import namedtuple
from graph import Graph
from numpy import sqrt,ceil,random
from itertools import product

# index conversion methods

class Tup(namedtuple('Tup', 'm n h l')):
    '''4-tuple format: (m,n,h,l)

    paramters:
       m   : tile row
       n   : tile column
       h   : qubit is horizontal
       l   : qubit subindex within tile: 0 <= l < L
    '''
    __slots__ = ()


def linearToTuple(n, M, N, L=4, index0=False):
    '''Convert linear index of a qubit in an (N, M, L) processor to
    (m,n,h,l) 4-tuple format.

    inputs:
        n       : linear index
        M       : number of tile rows
        N       : number of tile columns
        L       : number of qubits per half-tile
        index0 : smallest linear index is 0
    '''

    qpr = 2*N*L     # qbits per row
    qpt = 2*L       # qbits per tile

    if not index0:
        ind -= 1

    row, rem = divmod(ind, qpr)
    col, rem = divmod(rem, qpt)
    horiz, ind = divmod(rem, L)

    return Tup(row, col, horiz, ind)

def tupleToLinear(tup, M, N, L=4, index0=False):
    '''Convert 4-tuple (m,n,h,l) index of a qubit in an (N, M, L) processor to
    linear format

    inputs:
        tup     : 4-tuple index
        M       : number of tile rows
        N       : number of tile columns
        L       : number of qubits per half-tile
        index0 : smallest linear index is 0
    '''

    qpr = 2*N*L     # qbits per row
    qpt = 2*L       # qbits per tile

    return (0 if index0 else 1) + qpr*tup[0]+qpt*tup[1]+L*tup[2]+tup[3]


# Chimera structure

class Tile:
    '''Container for a the qubits in a physical tile'''

    def __init__(self, n, m, qbs=[]):
        self.n, self.m, self.qbs = n, m, qbs

class Chimera(Graph):
    '''Specialized Graph subclass for the Chimera hardware connectivity graph.'''

    L = 4   # number of qubits per half-tile
    M = 12  # number of tile rows
    N = M   # number of tile columns

    NB = 1  # number of header lines in Chimera graph file

    def __init__(self, fn=None, qb_yield=1.):
        '''Construct a new Chimera graph instance

        inputs:
            fn          : Chimera graph file
            qb_yield    : Qubit yield if not from file, #qbs/(2*M*N*L)
        '''

        super(Chimera, self).__init__()

        if fn is not None:
            self.fromFile(fn)
        else:
            self.__generateChimera(qb_yield)

        # fill tiles with qubit indices for convenience
        self.__fillTiles()

    # File I/O

    def fromFile(self, fn):
        '''Load the Chimera graph from file'''
        # construct graph from file
        super(Chimera,self).fromFile(fn, nb=self.NB, mp=int)

        # Chimera parameters
        self.M = self.N = int(round( sqrt(len(self)/(2*self.L)) , 0))

    def toFile(self, fn):
        '''Save the Chimera graph to file'''

        M, N, L = self.M, self.N, self.L
        with open(fn, 'w') as fp:
            # header content
            fp.write('{0} {1}\n'.format(2*M*N*L, self.number_of_edges()))
            # collect list of included couplers
            edges = list(self.edges)
            for i,j in sorted(edges, key=lambda x: x[::-1]):
                fp.write('{0} {1}\n'.format(i,j))


    def __fillTiles(self):
        '''Store the qubits included in each tile for convenience'''

        self.tiles = [[Tile(m,n) for n in range(self.N)] for m in range(self.N)]
        for m,n in product(range(self.M), range(self.N)):
            tile = self.tiles[m][n]
            qb0 = tupleToLinear((m,n,0,0), self.M, self.N, self.L)
            tile.qbs = [qb for qb in range(qb0, qb0+2*self.L) if qb in self.nodes]

    def __generateChimera(self, qb_yield=1.):
        '''Randomly generate the Chimera graph with the given qubit yield. The
        number of qubits in the resulting graph is qb_yield*(2*M*N*L)'''

        self.clear()
        self.tiles = None

        # generate a complete Chimera graph
        qpr = 2*self.N*self.L   # qbits per row
        qpt = 2*self.L          # qbits per tile

        # add all qubits
        for n in range(2*self.M*self.N*self.L):
            self.add_node(n+1)

        # add all couplers
        for m,n in product(range(self.M),range(self.N)):
            qbmap = lambda h,l: tupleToLinear((m,n,h,l), self.M, self.N, self.L)
            # coupler within the tile
            for l1, l2 in product(range(self.L), range(self.L)):
                self.add_edge(qbmap(0,l1), qbmap(1,l2))
            # couplers between rows
            if m<self.M-1:
                for l in range(self.L):
                    qb = qbmap(0,l)
                    self.add_edge(qb, qb+qpr)
            if n<self.N-1:
                for l in range(self.L):
                    qb = qbmap(1,l)
                    self.add_edge(qb, qb+qpt)

        # select broken qubits
        n_bqb = int(round((1-qb_yield)*self.size(),0)) # number of broken qubits
        br_qbs = random.choice(self.nodes.keys(),n_bqb, replace=False).tolist()

        # remove broken qubits
        for qb in br_qbs:
            self.remove_node(qb)


if __name__ == '__main__':

    import sys

    try:
        fn = sys.argv[1]
    except:
        print('No chimera file given, generating chimera')
        fn = None

    chimera = Chimera(fn)
    chimera.toFile('../temp/chimera2')
