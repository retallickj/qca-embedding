#!/usr/bin/env python
# encoding: utf-8

'''
Classes for describing general graph structures with optional location
information. Base for the source graph structure in the embedding GUI.
'''

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2017-10-31'      # last update

from copy import copy

class Node(object):
    '''Node class for Graph construction

    Parameters:
        n   : name or label for the node
        adj : adjacency dict with adj[m] the edge weight between n and m
        v   : value contained in the node
        x   : optional x location information
        y   : optional y location information
    '''

    def __init__(self, n, adj=dict(), v=0., x=None, y=None):
        '''Construct a Node instance'''
        self.n, self.adj, self.v, self.x, self.y = n, adj, v, x, y

    def copy(self):
        '''Return a deep copy of the node, assuming 'adj' is a dict of floats and
        and other parameters are copyable '''
        n, v, x, y = [copy(z) for z in [self.n, self.v, self.x, self.y]]
        return Node(n, dict(self.adj), v, x, y)

    def echo(self):
        '''Human readable node information for debug'''

        from pprint import pprint

        print('Node: {0}::{1}'.format(self.n,self.v))
        if self.x is not None:
            print('Location: {0}::{1}'.format(self.x,self.y))
        print('Adjacency:')
        pprint(self.adj)
        print('\n')
        n.echo()



class Graph(object):
    '''General base class for qubit graphs'''

    file_delim = ' '    # graph file delimiter

    def __init__(self):
        '''Initialise a graph instance'''

        self.nodes = {}     # dictionary of nodes

    def fromFile(self, fn, nb=0, mp=int):
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
        self.nodes = {}

        try:
            with open(fn, 'r') as fp:
                # burn header
                for n in range(nb):
                    fp.readline()
                self.__readFile(fp, mp)
        except IOError:
            print('Failed to read graph source file: {0}'.format(fn))


    def addNode(self, n, v=0., x=None, y=None):
        '''Add an unconnected node to the Graph'''

        assert n not in self.nodes, 'Node already exists'
        self.nodes[n] = Node(n, v=v, x=x, y=y)

    def addEdge(self, n, m, v=0.):
        '''Add an edge between two existing nodes'''

        assert n in self.nodes and m in self.nodes, 'Missing node'
        self.nodes[n].adj[m] = v
        self.nodes[m].adj[n] = v

    # graph methods

    def size(self):
        '''Return the number of nodes in the graph'''
        return len(self.nodes)

    def degree(self):
        '''Return the degree of the graph'''
        return max(len(n.adj) for n in self.nodes)

    def number_of_nodes(self):
        '''Return the number of nodes in the graph'''
        return self.size()

    def number_of_edges(self):
        '''Return the number of edges in the graph'''
        return int(sum(len(self.nodes[n].adj) for n in self.nodes)//2)

    def subgraph(self, nodes):
        '''Return the sub-graph composed of the given list of nodes.
        The nodes in the subgraph are deep-copies of those of the full graph.'''

        G = Graph()

        # get copy of the included nodes
        G.nodes = {k: self.nodes[k].copy() for k in nodes if k in self.nodes}

        # remove extra edges
        for k, node in G.nodes.items():
            node.adj = {j: v for j,v in node.adj.items() if j in G.nodes}

        return G

    # private methods

    def __readFile(self, fp, mp):
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
            v,x,y = [float(d) if d is not None else None for d in data[2:]]

            # add nodes if they do not exist
            for k in [i,j]:
                if k not in self.nodes:
                    self.nodes[k] = Node(k, adj={})

            # add node or edge information
            if i==j:
                node = self.nodes[i]
                node.v, node.x, node.y = v, x, y
            else:
                self.nodes[i].adj[j] = v
                self.nodes[j].adj[i] = v

if __name__ == '__main__':

    import sys
    from pprint import pprint

    try:
        fn = sys.argv[1]
    except:
        print('No file given')
        sys.exit()

    G = Graph()
    G.fromFile(fn, 1)

    for k,n in G.nodes.items():
        print('\n')
        n.echo()
