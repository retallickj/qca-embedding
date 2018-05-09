#!/usr/bin/env python
# encoding: utf-8

'''
Chimera hardware connectivity graph structure and related methods
'''

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2017-11-28'  # last significant update

from graph import Graph
from utility import range_product, dget

from numpy import sqrt, ceil
from numpy.random import choice
from collections import namedtuple, defaultdict
from itertools import chain, product


class Tup(namedtuple('Tup', 'm n h l')):
    '''4-tupleformat: (m,n,h,l)

    parameters:
        m   : tile row
        n   : tile column
        h   : qubit is horizontal
        l   : qubit subindex within tile 0 <= l < L
    '''
    __slots__ = ()

def linear_to_tuple(ind, M, N, L=4, index0=False):
    '''Convert linear index of a ubit in an (M, N, L) processor to a
    (m, n, h, l) 4-tuple

    inputs:
        ind     : linear index
        M       : number of tile rows
        N       : number of tile columns
        L       : number of qubits per half-tile
        index0 : smallest linear index is 0
    '''

    m, rem = divmod(ind if index0 else ind-1, 2*N*L)
    n, rem = divmod(rem, 2*L)
    h, l   = divmod(rem, L)

    return Tup(m, n, h, l)

def tuple_to_linear(tup, M, N, L=4, index0=False):
    '''Convert 4-tuple (m,n,h,l) of a qubit in an (M, N, L) processor to a
    linear index

    inputs:
        tup     : 4-tuple index
        M       : number of tile rows
        N       : number of tile columns
        L       : number of qubits per half-tile
        index0  : smallest linear index is 0
    '''

    return (not index0) + 2*N*L*tup[0] + 2*L*tup[1] + L*tup[2] + tup[3]


class Tile(namedtuple('Tile', 'm n qbs')):
    ''' Container for tile information and contained qubits

    parameters:
        m       : row index of tile
        n       : column index of tile
        qbs     : unsorted set of qubits in the tile, as node labels
    '''
    __slots__ = ()


class Chimera(Graph):
    '''Specialized Graph subclass for the Chimera hardware connectivity graph'''

    M = 16  # number of tile rows on chip
    N = M   # number of tile columns on chip
    L = 4   # number of qubits per half-tile

    # TODO: change this to a pattern for the Chimera file header
    NB = 1  # number of header lines in a Chimera graph file

    def __init__(self, G=None, fn=None, qb_yield=1., **kwargs):
        '''Construct a new Chimera graph instance

        inputs:
            G           : optional nx.Graph to construct from
            fn          : Chimera graph file
            qb_yield    : Qubit yield if not from file, #qbs/(2*M*N*L)

        kwargs:
            M   : number of tile rows
            N   : number of tile columns
            L   : number of qubits per half tile
        '''

        super(Chimera, self).__init__(G)

        self.tiles = defaultdict(dict)

        if fn is not None:
            self.from_file(fn)
        else:
            self._generate_chimera(qb_yield, **kwargs)


    def clear(self):
        '''Clear all Graph content'''
        super(Chimera, self).clear()
        self.M, self.N, self.L = Chimera.M, Chimera.N, Chimera.L

        self.tiles.clear()

    def subgraph(self, qbs):
        '''Get a subgraph of the Chimera graph'''
        G = super(Chimera, self).subgraph(qbs)
        print(len(G))

    def subchimera(self, m_set, n_set):
        '''Get a subgraph of the Chimera graph containing only qubits in the tiles
        specified by the m and n sets. Qubits maintian their numbering and
        4-tuple information.

        inputs:
            m_range : iterable of tile rows, integer gets mapped to [n]
            n_range : iterable of tile columns, integer gets mapped to [n]

        usage:
            # subgraph of qubits in tiles with m in [1,2,3] and n in [2,3,4]
            sub = chimera.subchimera(range(1,4), range(2,5))

            # subgraph of qubits in tiles with m in [2] and n in [5,6]
            sub = chimera.subchimera(2, [5,6])
            sub = chimera.subchimera([2], [5,6])
        '''

        # format integer inputs
        def inform(x):
            try:
                iter(x)
            except TypeError:
                x = [x]
            return x
        m_set, n_set = inform(m_set), inform(n_set)

        # get set of all qubits to include, could one-line this with chain
        qbs = set()
        for m,n in product(m_set,n_set):
            qbs |= self.tiles[m][n].qbs     # union shorthand

        return self.subgraph(qbs)

    # File I/O

    def from_file(self, fn):
        '''Load the Chimera graph from file'''

        # construct graph from file
        super(Chimera, self).from_file(fn, nb=self.NB, mp=int)

        # Chimera parameters
        self.L = 4
        self.M = self.N = int(round(sqrt(len(self)/(2*self.L)), 0))

        # assign node tuple fields
        for node in self:
            self.add_node(node, tup=linear_to_tuple(node, self.M, self.N, self.L))

        self._fill_tiles()

    def to_file(self, fn):
        '''Save the Chimera graph to file'''

        M, N, L = self.M, self.N, self.L    # help verbosity a bit
        with open(fn, 'w') as fp:
            # header content
            fp.write('{0} {1}\n'.format(2*M*N*L, self.number_of_edges()))
            # included couplers
            for i, j in sorted(self.edges, key=lambda x: x[::-1]):
                fp.write('{0} {1}\n'.format(i,j))

    # Helper Functions

    def _fill_tiles(self):
        '''Populate the tiles with their containined qubits'''

        M, N, L = self.M, self.N, self.L
        self.tiles.clear()

        # default tiles
        for m,n in range_product(M,N):
            self.tiles[m][n] = Tile(m, n, set())

        # add qubits to tiles
        for qb, data in self.nodes.data():
            m,n,h,l = data['tup']
            self.tiles[m][n].qbs.add(qb)


    def _generate_chimera(self, qb_yield=1., **kwargs):
        '''Randomly generate the Chimera graph with the given qubit yield. The
        number of qubits in the resulting graph is qb_yield*(2*M*N*L)

        kwargs:
            M   : number of tile rows
            N   : number of tile columns
            L   : number of qubits per half tile
        '''

        self.clear()

        # Chimera parameters
        self.M = dget(kwargs, 'M', self.M, int)
        self.N = dget(kwargs, 'N', self.M, int) # N == M by default
        self.L = dget(kwargs, 'L', self.L, int)
        M, N, L = self.M, self.N, self.L        # help verbosity a bit

        self.name = 'Chimera({0}, {1}, {2})'.format(M, N, L)

        # add all qubits
        for n in range(2*M*N*L):
            self.add_node(n+1, tup=linear_to_tuple(n+1, M, N, L))

        # add all couplers
        for m,n in range_product(M, N):
            mp = lambda h,l: tuple_to_linear((m,n,h,l), M, N, L)
            # inter-tile couplers
            self.add_edges_from((mp(0,u), mp(1,v)) for u,v in range_product(4,4))
            # horizontal couplers
            if n < N-1:
                for l in range(L):
                    qb = mp(1,l)
                    self.add_edge(qb, qb+2*L)
            # vertical couplers
            if m < M-1:
                for l in range(L):
                    qb = mp(0,l)
                    self.add_edge(qb, qb+2*N*L)

        # select broken qubits
        n_bqb = int(round((1-qb_yield)*len(self), 0))
        br_qbs = sorted(choice(self, n_bqb, replace=False))

        # remove broken qubits
        self.remove_nodes_from(br_qbs)

        self._fill_tiles()


    # automatic parameter detection

    def _detect_parameters(self):
        '''Detect the chimera parameters from the Graph structure.

        method:
            The set, Q, of absolute differences between the node indices of each
            edge is a subset of {1, 2, ..., 2L-1, 2L, 2LN}. We assume that
            min(M,N) > 1'''

        # get estimate of lower bound on L
        trials, L = 0, 1
        while trials < 10:
            qb = choice(self)   # randomly select a qubit in the graph
            G = self_expand_to_tile(qb) # subset of nodes in qubit tile
            L = max(L, ceil(.5*(max(G)-min(G)+1)))
            trials += 1

        if max(self) <= 2*L:
            # believe the chimera graph is a single tile
            return 1, 1, int(max(L, ceil(.5*max(self))))

    def _expand_to_tile(self, qb):
        '''Attempt to find a subset of the nodes contained in the tile of the
        given qubit'''

        # find the largest biconnected subgraph of the nodes 2 edges away from qb
        S = set(chain.from_iterable(self[k] for k in self[qb]))
        Gs = list(nx.biconnected_components(self.subgraph(S)))
        if Gs:
            return max(Gs, key=lambda x:len(x))
        else:
            return self.subgraph([qb])


if __name__ == '__main__':

    import sys

    try:
        fn = sys.argv[1]
    except:
        print('No chimera file given, generating chimera')
        fn = None

    chimera = Chimera(fn=fn, qb_yield=1.)
    # chimera.to_file('../temp/ch16b.txt')

    G = chimera.subchimera([1,2], 2)

    from pprint import pprint
    pprint(G)
