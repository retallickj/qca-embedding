#!/usr/bin/env python
# encoding: utf-8

'''
Classes for describing general graph structures with optional location
information. Base for the source graph structure in the embedding GUI.
'''

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2017-11-06'      # last update

from copy import copy
import networkx as nx

class Graph(nx.Graph):
    '''General base class for qubit graphs.'''

    file_delim = ' '    # graph file delimiter

    def __init__(self, G=None, fn=None, nb=0, mp=int):
        '''Initialise a Graph instance. If fn is give, initializes the Graph
        from a graph file (see Graph.fromFile for details).

        parameters:
            G   : optional nx.Graph to construct from
            fn  : filename
            nb  : number of header lines to burn before graph information
            mp  : map to apply to node labels
        '''

        super(Graph, self).__init__(G)

        if fn is not None:
            try:
                self.fromFile(fn, nb, mp)
            except AssertionError as e:
                print('Graph initialization failed with error:\n{0}'.format(e.message))
                self.clear()

    def clear(self):
        '''Reset all data in the Graph'''
        super(Graph, self).clear()  # clear nx.Graph data


    def from_file(self, fn, nb=0, mp=int):
        '''Populate graph from file. Each line of the file must be a tuple of the
        form '<i> <j> <val> <x> <y>' with only the first three arguments needed.
        If i==j, then val describes the node value and x and y are optional
        position information. Otherwise, val is the edge weight and x and y are
        ignored.

        inputs:
            fn  : filename
            nb  : number of header lines to burn before graph information
            mp  : map to apply to node labels
        '''

        assert isinstance(nb, int) and nb >= 0, 'Invalid burn count'

        # reset any existing graph information
        self.clear()

        # read graph file
        try:
            with open(fn, 'r') as fp:
                # burn header
                for n in range(nb):
                    fp.readline()
                self.__read_file(fp, mp)
        except IOError:
            print('Failed to read graph source file: {0}'.format(fn))


    def subgraph(self, nodes, copy=False):
        '''Return the sub-graph composed of the given list of nodes. The nodes
        in the subgraph are deep-copies of those of the full graph. The returned
        graph must be cast to the appropriate derived class.'''

        G0 = nx.Graph(super(Graph, self).subgraph(nodes))
        G = self.__new__(type(self))
        G.__init__(G0)

        return G


    # private methods

    def __read_file(self, fp, mp):
        '''Read the node information from a file pointer

        inputs:
            fp  : Python-style file pointer in a 'with' block
            mp  : node label map
        '''

        for line in fp:
            data = line.split(self.file_delim)
            if len(data)==2:
                data.append(0)
            if len(data)==3:
                data += [None,]*2
            elif len(data) != 5:
                print('Invalid line format: {0} ... skipping'.format(line))
                continue
            i,j = [mp(d) for d in data[:2]]
            w,x,y = [float(d) if d is not None else None for d in data[2:]]

            # add nodes
            for k in [i,j]:
                self.add_node(k)    # does nothing if node already exists

            # add node or edge information
            if i==j:
                self.add_node(i, weight=w, x=x, y=y)
            else:
                self.add_edge(i, j, weight=w)


if __name__ == '__main__':

    import sys
    from pprint import pprint

    try:
        fn = sys.argv[1]
    except:
        print('No file given')
        sys.exit()

    G = Graph()
    G.from_file(fn, 1)
