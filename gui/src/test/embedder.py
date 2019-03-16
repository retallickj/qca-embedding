#!/usr/bin/env python
# encoding: utf-8

import sys
assert sys.version_info < (3,), 'Embedding methods only support Python 2'

from traceback import print_exc

from chimera import Chimera

from itertools import product, combinations

_embedders = {'dense': True,
             'heur': True}

try:
    from dense_embed.embed import denseEmbed, setChimera
    from dense_embed.convert import convertToModels
except Exception as e:
    print('Could not load dense embedding method')
    _embedders['dense'] = False
    print(print_exc())

try:
    from dwave_sapi2.embedding import find_embedding
except Exception as e:
    print('Could not load heuristic embedding method')
    _embedders['heur'] = False
    print(print_exc())



class Embedder:

    zero = 1e-3  # zero threshold value


    def __init__(self, method):
        '''Initialise an embedding instance

        input:
            method  : embedder type, must be in embedders
        '''

        assert method in _embedders, 'Invalid embedding method'
        assert _embedders[method], 'Given embedding method not available'

        self.method = method


    def setChimera(self, M, N=None, L=4, qb_yield=1):
        '''Set the adjacency structure for an (M,N=M,L) chimera graph with the
        given yield'''

        if N is None:
            N = M

        self.chimera = Chimera(M=M, N=N, L=L, qb_yield=qb_yield)


    def run(h, J, ntrials=1):
        '''Run the embedder for the given (h,J) problem. Performs multiple
        trials if n_trials set'''

        self.h = np.array(h).reshape([-1,])
        self.J = np.array(J)

        self.N = len(self.h)
        assert self.J.shape == (self.N,)*2, "J matrix doesn't match size of h"

        self.J /= np.max(np.abs(self.J))

    def _run_dense(self, h, J, ntrials):
        '''Setup and run the Dense Placement algorithm'''
        pass

    def _run_heur(self, h, J, ntrials):
        '''Setup and run the Heuristic algorithm'''

        S = set()
        for ij in combinations(range(N),2):
            if abs(J[ij]) > self.zero:
                S.add(ij)


if __name__ == '__main__':

    embedder = Embedder('dense')
