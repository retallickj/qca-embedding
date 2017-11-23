#!/usr/bin/python

# -----------------------------------------------------------------------------
# Name: assign.py
# Purpose:  qubit and coupler parameter assignment
# Author:   Jacob Retallick
# Created:  24.06.2015
# Last modified: 24.06.2015
# -----------------------------------------------------------------------------

import numpy as np
import networkx as nx
import itertools
import re

# ASSIGNMENT PARAMETERS

strategies = ['maximal',    # Use as many couplers as possible
              'minimal'     # Use as few couplers as possible
              ]

strategy = 'maximal'
assert strategy.lower() in strategies, 'Invalid edge selection strategy'


def sort_dict(d):
    '''Sort the keys of a dict by d[key]'''
    return zip(*sorted([(d[k], k) for k in d]))[1]


def pinch(s, s1, s2):
    return s.partition(s1)[2].rpartition(s2)[0]


# Edge selection strategies


def maximal_couplers(subgraphs, edges):
    '''Use as many edges as possible between and within subgraphs'''
    return subgraphs, edges


def minimal_couplers(subgraphs, edges):
    '''Use the fewest possible number of couplers between and within
    subgraphs'''

    N = len(subgraphs)

    # map each subgraph to its minimum spanning tree
    subgraphs = [nx.minimum_spanning_tree(subgraph) for subgraph in subgraphs]

    # for each tree, find a root node and store the shortest path to each
    # node as a cost metric.
    costs = {}
    for tree in subgraphs:
        # identify the root
        path_lengths = nx.shortest_path_length(tree)
        root_weights = {k: sum(path_lengths[k].values()) for k in path_lengths}
        root = sort_dict(root_weights)[0]
        # assign path lengths as node costs
        for node in path_lengths[root]:
            costs[node] = path_lengths[root][node]

    # for each pair of subgraphs, keep the inter-subgraph edge with the
    # minimum total cost of its end nodes
    nodes = sorted(subgraphs.keys())
    for i in range(N-1):
        q1 = nodes[i]
        for j in range(i+1, N):
            q2 = nodes[j]
            edge_costs = {e: costs[e[0]]+costs[e[1]] for e in edges[(q1, q2)]}
            edges[(q1, q2)] = sort_dict(edge_costs)[0]

    return subgraphs, edges


#


def write_to_file(hq, Jq, fn):
    '''Write the parameters to file.

    inputs: hq  : dictionary of h parameters for used qubits
            Jq  : dictionary of J parameters for used couplers
            fn  : file path
    '''

    try:
        fp = open(fn, 'w')
    except IOError:
        print('Invalid file path... {0}'.format(fn))
        return None

    fp.write('<parameters>\n')

    # write the h parameters
    fp.write('<h>\n')
    for qbit in hq:
        fp.write('\t{0}:\t{1}\n'.format(qbit, hq[qbit]))
    fp.write('</h>\n\n')

    # write the J parameters
    fp.write('<J>\n')
    for coupler in Jq:
        fp.write('\t{0}:\t{1}\n'.format(coupler, Jq[coupler]))
    fp.write('</J>\n\n')

    fp.write('</parameters>')


def read_from_file(fn):
    '''Read the h and J parameters of an embedding from file'''

    try:
        fp = open(fn, 'w')
    except IOError:
        print('Invalid file path... {0}'.format(fn))
        return None

    # reading flags
    reading_h = False
    reading_J = False

    # regex
    re_start = re.compile('^<\w*>$')
    re_end = re.compile('^</\w*>$')

    hq = {}
    Jq = {}

    for line in fp:
        if '#' in line or len(line) < 3:
            continue
        if re_start.match(line.strip()):
            key = pinch(line, '<', '>').strip()
            if key == 'h':
                reading_h = True
            elif key == 'J':
                reading_J = True
        elif re_end.match(line.strip()):
            key = pinch(line, '</', '>').strip()
            if key == 'h':
                reading_h = False
            elif key == 'J':
                reading_J = False
        elif reading_h:
            qbit, h = line.strip().split(':')
            qbit = int(qbit)
            hq[qbit] = float(h)
        elif reading_J:
            coupler, J = line.strip().split(':')
            coupler = tuple(map(int, pinch(coupler, '(', ')').split(',')))
            Jq[coupler] = float(J)

    return hq, Jq


def partition_graph(G, parts):
    '''Partition graph G into a set of disjoint subgraphs given by a list of
    lists of node labels for each subgraph

    inputs: G       : Graph object to be partitioned
            parts   : dict of lists of subgraph nodes. Must have same labels as
                     in G
    '''

    N = len(parts)  # number of partitions

    # get partition subgraphs
    subgraphs = {}
    for node in parts:
        part = parts[node]
        subgraph = G.subgraph(part)
        if len(subgraph) < len(part):
            conflicts = [n for n in part if not n in subgraph.nodes()]
            raise KeyError('Invalid indices given: {0}'.format(conflicts))
        subgraphs[node] = subgraph

    # get list of edges between each subgraph
    edges = {}
    nodes = sorted(parts.keys())
    for i1 in range(N-1):
        n1 = nodes[i1]
        for i2 in range(i1+1, N):
            n2 = nodes[i2]
            edges[(n1, n2)] = []
            for u, v in list(itertools.product(parts[n1], parts[n2])):
                edge = sorted([u, v])
                if G.has_edge(*edge):
                    edges[(n1, n2)].append(edge)

    return subgraphs, edges


def convert_to_parameters(h, J, subgraphs, edges, J_inner):
    '''Construct parameter dictionaries from the problem h and J coefficients
    and the determined subgraphs and edges'''

    assert len(h) == len(subgraphs),\
        'Mismatch between problem nodes and subgraphs'

    N = len(h)

    hq = {}     # h parameter for each qubit used: keys are integers
    Jq = {}     # J parameter for each coupler: key format (u, v) with u < v

    # scale factor for handling general J_inner
    h_max = max(abs(x) for x in h.values())
    j_max = max(max(abs(x) for x in j.values()) for j in J.values())
    hj_max = max(h_max, j_max)
    scale = float(max(hj_max, abs(J_inner)))
    print('J_inner set to: {0}'.format(J_inner))

    # handle internal subgraph parameters
    for node in h:
        subgraph = subgraphs[node]
        for qbit in subgraph.nodes():
            hq[qbit] = h[node]*1./subgraph.number_of_nodes()/scale
        for q1, q2 in subgraph.edges():
            Jq[tuple(sorted((q1, q2)))] = J_inner/scale

    # handle inter-subgraph parameters
    nodes = sorted(subgraphs.keys())
    for i in range(N-1):
        n1 = nodes[i]
        for j in range(i+1, N):
            n2 = nodes[j]
            if n2 in J[n1]:
                for q1, q2 in edges[(n1, n2)]:
                    Jq[tuple(sorted((q1, q2)))] = J[n1][n2]*1./len(edges[(n1, n2)])/scale

    return hq, Jq


def assign_parameters(h, J, qbits, chimera, flip_J=False, J_inner=-1):
    '''Determine the h and J coefficients for a given embedding problem given
    a list of qubits for each problem node and an adjacency list for the
    target chimera structure (with qbit labels as verticed)

    inputs: h       : dict of on-site terms for each problem node
            J       : dict of coupling terms between each problem node
            qbits   : list of qbits for each problem node
            chimera : adjacency list structure for the target chimera graph
            flip_J   : flag for flipping the sign of J
    '''

    # check that h and J are normalised
    max_hj = max(np.max(np.abs(h.values())),
                 max([max(J[node].values()) for node in J]))
    if max_hj == 0:
        print('Invalid problem statement. All zero parameters')
        return None
    if max_hj != 1:
        print('Normalizing h and J by factor {0}'.format(max_hj))
        h = {node: h[node]/max_hj for node in h}
        J = {n1: {n2: J[n1][n2]/max_hj for n2 in J[n1]} for n1 in J}

    # flip J signs if flagged
    if flip_J:
#        print('Flipping signs of J coefficients')
        J = {n1: {n2: -J[n1][n2] for n2 in J[n1]} for n1 in J}

    # build chimera graph
    G_chimera = nx.Graph(chimera)

    # get subgraphs and edge lists for problem node qbit lists
    try:
        subgraphs, edges = partition_graph(G_chimera, qbits)
    except KeyError as e:
        print('Qbit label error during Chimera graph partition...')
        print(e.message)
        return None, None

    # remove unwanted edges
    if strategy.lower() == 'maximal':
        subgraphs, edges = maximal_couplers(subgraphs, edges)
    elif strategy.lower() == 'minimal':
        subgraphs, edges = minimal_couplers(subgraphs, edges)
    else:
        print('No valid edge selection strategy given...')
        return None, None

    # convert subgraphs and edges to parameter dictionaries
    hq, Jq = convert_to_parameters(h, J, subgraphs, edges, J_inner=J_inner)

    # return parameter
    return hq, Jq
