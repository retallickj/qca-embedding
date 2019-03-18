#!/usr/bin/env python
#encoding: utf-8

'''
'''

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '1.0'
__date__        = '2019-03-16'

import networkx as nx
import re

class GraphException(Exception):
    '''Defined errors during Graph operations'''
    pass


class Graph(nx.Graph):
    '''General base class for qubit graphs'''


    def __init__(self, G=None, fn=None, nb=0, mp=int):
        '''Initialise a Graphj instance. If fn is given, initialize the Graph
        from a graph file (see Graph.from_file for details)
        
        inputs:
            G   : optional nx.Graph to construct from
            fn  : filename of graph file
            nb  : number of lines to burn in the graph file
            mp  : map to apply to node labels
        '''

        super(Graph, self).__init__(G)

        if fn is not None:
            try:
                self.from_file(fn, nb=nb, mp=mp)
            except GraphException as e:
                print('Graph initialization failed:\n{0}'.format(e.message))
                self.clear()


    def clear(self):
        '''Reset all data in the Graph'''
        super(Graph, self).clear()


    def from_file(self, fn, nb=0, mp=int):
        '''Populate graph from file. Each line of the file must be a tuple of
        the form "<i> <j> <val> <x> <y>" with <x> and <y> optional. Any line
        which does not match

        valid format:

            i i val x y :   node with label mp(i), value float(val), and 
                            position (float(x), float(y))

            i i val : node without position information

            i j val : edge between nodes mp(i) and mp(j) with value float(val)

            i j : node or edge with val=0

        inputs:
            fn  : filename
            nb  : number of lines to burn
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
                self._read_file(fp, mp)
        except IOError:
            raise GraphException(
                'Failed to read graph source file: {0}'.format(fn))

    def to_file(self, fn):
        '''Write teh graph to file'''

        try:
            with open(fn, 'w' ) as fp:
                self._write_file(fp)
        except IOError:
            print('Failed to write graph  to file: {0}'.format(fn))

    def subgraph(self, nodes, copy=False):
        '''Return the sub-graph composed of the gioven list of nodes. The nodes
        in the subgraph are deep-copies of those of the full graph. The returned
        graph must be cast to the appropriate derived class.'''

        # nx.Graph version of subgraph
        G0 = nx.Graph(super(Graph, self).subgraph(nodes))

        # forward cast to relevant nx.Graph derived class
        G = self.__new__(type(self))
        G.__init__(G0)

        return G

    # private methods

    def _read_file(self, fp, mp):
        '''Read the node information from a file pointer
        
        inputs:
            fp  : file pointer
            mp  : node label map 
        '''

        for line in fp:
            data =  re.findall('[\-\+\w\.]+', line)
            if len(data) == 2:
                data.append(0)
            if len(data) == 3:
                data += [None,]*2
            elif len(data) != 5:
                print('Invalid lione format: {0} .. skipping'.format(line))
                continue
            i,j = [mp(d) for d in data[:2]]
            w,x,y = [float(d) if d is not None else None for d in data[2:]]

            # add nodes
            for k in [i,j]:
                self.add_node(k)    # does nothing if node already exists

            # add node or edge information
            if i==j:
                self.add_node(i, weigh=w, x=x, y=y)
            else:
                self.add_edge(i, j, weight=w)

    def _write_file(self, fp):
        '''Write node information to file
        
        inputs:
            fp  : file pointer
        ''' 

        def extract(obj, key, dft):
            return obj[key] if key in obj else dft
                
        for n in self.nodes:
            node = self.node[n]
            val = extract(node,'val',0)
            x,y = [extract(node, k, None) for k in 'xy']
            s = '{0} {0} {1}'.format(n,val)
            if x is not None and y is not None:
                s += ' {0} {1}'.format(x, y)
            fp.write(s+'\n')

        for (n,m) in self.edges:
            edge = self.edges[(n,m)]
            val = extract(edge, 'val',0)
            s = '{0} {1} {2}'.format(n,m,val)
            fp.write(s+'\n')




if __name__ == '__main__':

    if False:
        G = Graph()

        G.add_node(1, val=-1, x=0, y=1)
        G.add_node(2, val=1, x=-1, y=0)
        G.add_node(3, val=1, x=0, y=-1)
        G.add_node(4, val=0, x=0, y=0)
        G.add_node(5, val=0, x=1, y=0)

        for n in [1,2,3,5]:
            G.add_edge(n,4,val=-1)
        for n,m in [(1,2),(2,3),(3,5),(1,5)]:
            G.add_edge(n,m,val=0.2)

        import os
        fn = os.path.join('cache', 'temp_graph.txt')
        G.to_file(fn)
    else:
        import sys
        fn = sys.argv[1]
        G = Graph(fn=fn, mp=int)
        import pdb; pdb.set_trace()

        

