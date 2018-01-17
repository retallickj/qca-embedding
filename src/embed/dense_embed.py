#!/usr/bin/env python
# encoding: utf-8

'''
Classes for embedding graphs onto the Chimera architecture using the Dense
Placement algorithm.
'''

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2018-01-17'      # last update


class DenseEmbedder:
    '''Handler class for running the dense placement embedding algorithm.
    Attempts to find an embedding of the source graph into a target Chimera
    graph.'''

    def __init__(self, chimera=None):
        '''Initialise a DenseEmbedder instance.

        parameters:
            chimera : optional Chimera graph
        '''
        pass

    def set_chimera(self, chimera):
        '''Update the Chimera graph for the Embedder

        input:
            chimera : Chimera graph. Must be of type Chimera (see core.chimera)
        '''
        pass

    def embed(self, source, **kwargs):
        '''Attempt to find an embedding of the given source graph onto the
        current Chimera graph.

        inputs:
            source  : source graph, can be any networkx.Graph or derived

        kwargs:
            verbose : if True, verbose output for the embedding
            ntrials : number of embedding attempts
            best    : if True, return only the most economical embedding
        '''

        pass
