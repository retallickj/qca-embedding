#!/usr/bin/env python

#---------------------------------------------------------
# Name: parse_qca.py
# Purpose: Parsing functions for QCADesigner files
# Author: Jacob Retallick
# Created: 2015.10.22
# Last Modified: 2015.10.22
#---------------------------------------------------------

# NOTE
# the original parse script no longer seems to work (change in networkx?)
# for the purposes of the embedder, we don't need to consider clock zones so
# I have simplified the parseing script to remove that functionality.

import re

import networkx as nx
import numpy as np
from core.auxil import getEk, CELL_FUNCTIONS, CELL_MODES
from itertools import combinations, chain

from pprint import pprint

## general global parameters

R_MAX = 2.1         # max cell-cell interaction range (rel to grid spacing)
EK_THRESH = 1e-3    # threshold for included Ek, relative to max(abs(Ek))
X_ROUND = 4         # places to round to when deciding if cell is rotated


### FILE PROCESSING

def build_hierarchy(fn):
    '''Build a dict hierarchy containing all objects, their parameters, and
    childen.'''

    fp = open(fn, 'r')

    linemap = lambda s: s.replace(',', '.')

    # general re expression. may need to change if future format changes
    re_start = re.compile('^\[.+\]$')
    re_term = re.compile('^\[#.+\]$')

    hier = {'label': 'Hierarchy', 'children': [], 'vars': {}}

    key_stack = ['Hierarchy']  # stack of active keys, pop of top of stack
    dict_stack = [hier]   # stack of corresponding dict objects.

    line_cnt = 0
    for line in fp:
        line = linemap(line)
        line_cnt += 1
        line = line.strip()  # remove endline and possible whitespace

        # must check object termination first
        if re_term.match(line):
            key = line[2:-1]
            if key_stack[-1] == key:
                d = dict_stack.pop()
                key_stack.pop()
                try:
                    dict_stack[-1]['children'].append(d)
                except:
                    print('Somehow over-popped dict_stack...')
                    return None
            else:
                print('Start-end mismatch in line {0}'.format(line_cnt))
                return None

        # for a new object, create a new dict template
        elif re_start.match(line):
            key = line[1:-1]
            key_stack.append(key)
            d = {'label': key, 'children': [], 'vars': {}}
            dict_stack.append(d)

        # otherwise check for new variable to add to most recent dict
        else:
            if '=' in line:
                var, val = line.split('=')
                dict_stack[-1]['vars'][var] = val
    fp.close()

    return hier


def proc_hierarchy(hier):
    '''Process the extracted data hierarchy to extract useful information. In
    the current information, we are interested in the overall cell grid spacing
    (for deciding on the range of included cell) and the properties of each
    cell in the circuit'''

    cells = []
    spacing = None

    # hierarchy should only have two children: VERSION and TYPE:DESIGN. The
    # former might be useful in later implentations for selecting formatting
    # options but for now all we care about are the DESIGN objects

    hier = [child for child in hier['children']
            if child['label'] == 'TYPE:DESIGN'][0]

    # for now assert that there can be only one cell layer, no vertical x-over
    layers = [child for child in hier['children']
              if child['label'] == 'TYPE:QCADLayer']

    # isolate cell layers
    cell_layers = [layer for layer in layers if layer['vars']['type'] == '1']

    # merge cell layers, will lead to qdot conflict if vertical x-over
    cell_dicts = [layer['children'] for layer in cell_layers]
    cell_dicts = list(chain.from_iterable(cell_dicts))

    # get grid-spacing (average cell bounding box)
    cx = float(cell_dicts[0]['vars']['cell_options.cxCell'])
    cy = float(cell_dicts[0]['vars']['cell_options.cyCell'])

    spacing = np.sqrt(cx*cy)
    # create cell objects
    cells = []

    for cd in cell_dicts:
        cell = {}

        # cell type
        cell['cf'] = CELL_FUNCTIONS[cd['vars']['cell_function']]
        cell['cm'] = CELL_MODES[cd['vars']['cell_options.mode']]
        cell['clk'] = int(cd['vars']['cell_options.clock'])

        # just for show sol
        cell['cx'] = float(cd['vars']['cell_options.cxCell'])
        cell['cy'] = float(cd['vars']['cell_options.cyCell'])

        # position, first child will be the QCADesignObject
        design_object = cd['children'][0]
        cell['x'] = float(design_object['vars']['x'])
        cell['y'] = float(design_object['vars']['y'])

        # quantum dots
        qdot_dicts = [child for child in cd['children']
                      if child['label'] == 'TYPE:CELL_DOT']

        qdots = []
        for d in qdot_dicts:
            dot = {}
            dot['x'] = float(d['vars']['x'])
            dot['y'] = float(d['vars']['y'])
            dot['c'] = float(d['vars']['charge'])
            qdots.append(dot)

        cell['qdots'] = qdots

        # determine if cell is rotated, will have three x values
        x = set([round(dt['x'], X_ROUND) for dt in qdots])
        if len(x) == 3:
            cell['rot'] = True
        elif len(x) == 2:
            cell['rot'] = False
        else:
            print('Could not decide cell rotation')
            cell['rot'] = False

        # keep track of polarization if cell is fixed: don't rely on labels
        if cell['cf'] == CELL_FUNCTIONS['QCAD_CELL_FIXED']:
            pol = qdots[0]['c']+qdots[2]['c']-qdots[1]['c']-qdots[3]['c']
            pol /= qdots[0]['c']+qdots[2]['c']+qdots[1]['c']+qdots[3]['c']
            cell['pol'] = pol

        cells.append(cell)

    return cells, spacing


## CIRCUIT PROCESSING

def build_J(cells, spacing, r_max=R_MAX):
    '''Build the J matrix for the given circuit. Restricts the interaction
    distance to r_max but does not apply any adjacency contraints'''

    N = len(cells)

    # contruct connectivvity matrix
    J = np.zeros([N, N], dtype=float)
    DR = r_max*spacing
    for i,j in combinations(range(N), 2):
        Ek = getEk(cells[i], cells[j], DR=DR)
        if Ek:
            J[i,j] = J[j,i] = Ek

    # remove very weak interactions
    J = J*(np.abs(J) >= np.max(np.abs(J)*EK_THRESH))

    return J

def zone_cells(cells, spacing, show=False):
    '''Split cells into clock zones. Distinguishes disjoint zones with the
    same zone index'''

    N = len(cells)  # number of cells

    # construct connectivity matrix
    J = np.zeros([N, N], dtype=float)
    DR = R_MAX*spacing
    for i in range(N-1):
        for j in range(i+1, N):
            Ek = getEk(cells[i], cells[j], DR=DR)
            if Ek:
                J[i, j] = Ek
                J[j, i] = Ek

    # remove very weak interactions
    J = J * (np.abs(J) >= np.max(np.abs(J)*EK_THRESH))

    # make full cell connectivity Graph
    G = nx.Graph(J)

#    if show:
#        plt.figure(0)
#        plt.clf()
#        nx.draw_graphviz(G)
#        plt.show()

    # get indices for each clock index
    clk = [cell['clk'] for cell in cells]
    clk_ind = list(set(clk))    # will sort by default
    inds = [[i for i, x in enumerate(clk) if x == ind] for ind in clk_ind]

    # split graph into sub-graphs with the same clock indices
    sub_G = {ind: G.subgraph(inds[ind]) for ind in clk_ind}

    # split disconnected components for each label graph
    sub_ind = {ind: list(nx.connected_components(sub_G[ind]))
               for ind in clk_ind}

    ## find zone order

    # create abstract zone connectivity graph
    G = nx.DiGraph()
    # nodes
    for clk in clk_ind:
        for i in range(len(sub_ind[clk])):
            key = (clk, i)
            G.add_node(key, inds=sub_ind[clk][i])
    # edges
    for clk in clk_ind:
        adj_clk = 3 if clk == 0 else clk-1
        if not adj_clk in sub_ind:
            continue
        for i in range(len(sub_ind[clk])):
            k1 = (clk, i)
            for j in range(len(sub_ind[adj_clk])):
                k2 = (adj_clk, j)
                if np.any(J[G.node[k1]['inds'], :][:, G.node[k2]['inds']]):
                    G.add_edge(k2, k1)

#    if show:
#        plt.figure(1)
#        plt.clf()
#        nx.draw_graphviz(G)
#        plt.show()

    # find input nodes, have no predecessors
    predecs = {n: len(G.predecessors(n)) for n in G.nodes_iter()}
    inputs = [ky for ky, val in predecs.items() if val == 0]

    # expand from inputs
    visited = {key: False for key in G.nodes()}
    nodes = inputs
    order = [nodes]
    while nodes:
        new_nodes = set()
        for node in nodes:
            new_nodes.update(G.successors(node))
            visited[node] = True
        # remove already visited nodes from new nodes
        new_nodes = [node for node in new_nodes if not visited[node]]
        nodes = new_nodes
        if nodes:
            order.append(nodes)

    # find feedback interactions
    feedback = {}
    for n in G.nodes_iter():
        for p in G.predecessors(n):
            pshell = 0
            nshell = 0
            pzone = 0
            nzone = 0
            for shell in order:
                if p in shell:
                    pshell = order.index(shell)
                    pzone = shell.index(p)
                if n in shell:
                    nshell = order.index(shell)
                    nzone = shell.index(n)
            if pshell > nshell:
                if (pshell,pzone) in feedback:
                    feedback[(pshell,pzone)].append((nshell,nzone))
                else:
                    feedback[(pshell,pzone)] =  [(nshell,nzone)]

    # reformat order list to contain zone indices
    form_func = lambda n: sub_ind[n[0]][n[1]]
    order = [[form_func(zone) for zone in shell] for shell in order]

    return order, J, feedback


def reorder_cells(cells, J, flipy=False):
    '''Renumber cells by position rather than the default QCADesigner placement
    order. Cells ordered by the tuple (zone, y, x)'''

    keys = {}
    ysgn = -1 if flipy else 1

    # assign sortable tuples for each cell
    for ind, cell in enumerate(cells):
        keys[ind] = (ysgn*cell['y'], cell['x'])

    order = list(zip(*sorted([(keys[i], i) for i in keys])))[1]

    # relabel cells and reorder the J matrix
    cells = [cells[i] for i in order]
    J = J[order, :][:, order]

    for i in range(len(cells)):
        cells[i]['num'] = i
        cells[i]['number'] = i

    return cells, J


## MAIN FUNCTION

def parse_qca_file(fn, verbose=False):
    '''Parse a QCADesigner file to extract cell properties. Returns an ordered
    list of cells, the QCADesigner grid spacing in nm, a list structure of the
    indices of each clock zone (propogating from inputs), and a coupling matrix
    J which contains the Ek values for cells within a radius of R_MAX times the
    grid spacing'''

    # build data hierarchy
    hier = build_hierarchy(fn)

    # extract useful information from data hierarchy
    cells, spacing = proc_hierarchy(hier)

    if verbose:
        print('Parsed QCA file...')

    for cell in cells:
        cell['clk'] = 0

    # construct J matrix
    J = build_J(cells, spacing)

    # reorder cells by zone and position
    cells, J = reorder_cells(cells, J)

    return cells, spacing, J
