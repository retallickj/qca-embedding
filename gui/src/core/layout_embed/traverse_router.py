'''
Created on Jun 16, 2017

@author: JosePinilla
'''
import sys
import random
import networkx as nx
import matplotlib.pyplot as plt
from virtualConflict import solveChains

VERBOSE = True
PLOT = True

M = 8
N = 8
L = 4
SEED = None
COST_THRESHOLD = 100000.0

# ROUTER FLAGS
LENGTH_COST = True
NEIGHBORHOOD = 'cross'
RANDOMIZE_CELLS = True
LENGTH_PRIORITY = True
RANDOMIZE_CANDIDATES = True


# History cost parameter
DELTA_H = 0.0
# Base cost constants
#DEGREE
BASE_A = 0.0
MAX_DEG = 6
#SCOPE
BASE_B = 0.0
#CONCENTRATION
BASE_C = 0.0
# Sharing cost parameter
DELTA_P = 0.0


## ROUTING RESOURCES
_alpha_p = 0            # Present-sharing cost scaler
_max_edge = sys.maxsize

_QCA = None
_Chimera =  None
_RGraph = None

def plotGraph(G):
    print("PLOT GRAPH")
    plt.clf()
    nx.draw(G,with_labels=True)
    plt.show()

def ripUpAll():
    
    ##### QCA Graph attributes
    ##########################

    # Initialize edges to be routed
    for u,v in _QCA.edges():
        _QCA.edges[u,v]['routed'] = False
        _QCA.edges[u,v]['path'] = []
    
    # Initialize node to be searched
    for node in _QCA:
        _QCA.nodes[node]['main'] = None
        _QCA.nodes[node]['qubits'].clear()
        _QCA.nodes[node]['routed'] = False
        
    
    ##### Chimera Graph attributes
    ##############################

    for node in _Chimera:
        _RGraph.nodes[node]['path'][:] = []
        _RGraph.nodes[node]['cells'].clear()
        _RGraph.nodes[node]['mapped'].clear()

def confChimeraGraph():
    '''
    
    :param chimera_adj:
    '''
    
    global M,N,L
    
    global _Chimera
    
    for node in _Chimera:
        
        # initialize empty set of cells assigned to qubit
        _Chimera.nodes[node]['cells'] = set()
        # initialize empty list of paths assigned to qubit
        _Chimera.nodes[node]['path'] = []
        # initialize empty dictionary of paths:set([cells])
        _Chimera.nodes[node]['mapped'] = {}
                
        # assign qubit tile from tuple(row,col,horiz,ind)
        row = node[0]
        col = node[1]
        _Chimera.nodes[node]['tile'] = (row, col)
        
        # NEGOTIATED CONGESTION    
        # historical cost
        _Chimera.nodes[node]['history'] = 1.0
        # Degree/2 (In directed graph)
        _Chimera.nodes[node]['degree'] = _Chimera.degree(node)/2
        
def confQCAGraph():
    '''
    
    :param QCA:
    '''
    
    global _QCA
    
    labels = {}
    for node in _QCA:
        # Initialize empty main qubits
        _QCA.nodes[node]['main'] = None
        # Initialize empty embedding on all cells
        _QCA.nodes[node]['qubits'] = set()
        # Initialize cell as not routed
        _QCA.nodes[node]['routed'] = False
        # Initialize cell model size
        _QCA.nodes[node]['size'] = 1
        # Initialize cell degree to not recompute
        degree = _QCA.degree[node]
        _QCA.nodes[node]['degree'] = degree 
        # Initialize cell for search
        _QCA.nodes[node]['routed'] = False
        # Initialize cell priority
        _QCA.nodes[node]['priority'] = degree
                       
        # Relabel QCA nodes to C<i>
        label = 'C' + str(node)
        labels[node] = label
    
    # Initialize size of patha
    nx.set_edge_attributes(_QCA, 0.0, 'size')
    
    # relabel nodes C<i> to avoid conflicts in routing graph F
    nx.relabel_nodes(_QCA, labels, copy=False)


def setupRGraph():
    '''
    
    :param chimera_adj:
    :param QCA:
    '''
    global _RGraph
    # Configure Chimera nodes and edges
    confChimeraGraph()
    # Configure QCA nodes and edges    
    confQCAGraph()
    # Make a disconnected copy of G
    V_QCA = nx.DiGraph()
    V_QCA.add_nodes_from(_QCA.nodes(data=False))
    # QCA U Chimera
    _RGraph = nx.compose(V_QCA, _Chimera)
    
    if PLOT:
        plotGraph(_QCA)
        plotGraph(_Chimera)
        plotGraph(_RGraph)

def neighbourTiles(tile):
    tile_x = tile % N
    tile_y = (tile-tile_x) / N

    n_tile = (tile - N)     if (tile_y > 0)      else   None
    s_tile = (tile + N)     if (tile_y < M-1)    else   None
    w_tile = (tile - 1)     if (tile_x > 0)      else   None
    e_tile = (tile + 1)     if (tile_x < N-1)    else   None
    
    nw_tile = (tile - N - 1)  if (tile_y > 0    and tile_x > 0)    else None
    ne_tile = (tile - N + 1)  if (tile_y > 0    and tile_x < N-1)  else None
    se_tile = (tile + N + 1)  if (tile_y < M-1  and tile_x < N-1)  else None
    sw_tile = (tile + N - 1)  if (tile_y < M-1  and tile_x > 0)    else None
    

    return n_tile, s_tile, w_tile, e_tile, nw_tile, ne_tile, se_tile, sw_tile

def getQubits(tile):
    '''
    List of qubits in the given tile index
    :param tile: Index of tile [ 0 : (M*N-1) ]
    '''
    
    tile_x = tile % N
    tile_y = (tile-tile_x) / N
    
    qubits = []
    for index in range(0,L):
        qubit_tuple_1 = (tile_y, tile_x, 0, index)
        qubit_tuple_2 = (tile_y, tile_x, 1, index)
        
        if _RGraph.nodes[qubit_tuple_1]['alive']:
            qubits.append(qubit_tuple_1)
        if _RGraph.nodes[qubit_tuple_2]['alive']:
            qubits.append(qubit_tuple_2)
            
    return qubits

def getCandidates(cell):
    '''
    
    :param F:
    :param cell:
    '''
    
    tile = _QCA.nodes[cell]['tile']
    n_tile, s_tile, w_tile, e_tile,  nw_tile, ne_tile, se_tile, sw_tile = neighbourTiles(tile)


    # Get tiles

    candidate_tiles = []
    
    if NEIGHBORHOOD=='tile':        candidate_tiles.append( tile ) 
        
    elif NEIGHBORHOOD=='cross':     candidate_tiles.extend( [tile, n_tile, s_tile, w_tile, e_tile] )

    elif NEIGHBORHOOD=='2-block':   candidate_tiles.extend( [tile, n_tile, s_tile, w_tile, e_tile, ne_tile, se_tile, nw_tile, sw_tile] )
    
    elif NEIGHBORHOOD=='block':
    
        cell_x = _QCA.nodes[cell]['cell']['x']
        cell_y = _QCA.nodes[cell]['cell']['y']
        tile_x = tile % N
        tile_y = (tile-tile_x) / N
        
        offset_x = cell_x - tile_x
        offset_y = cell_y - tile_y
        
        candidate_tiles.append( tile )
        
        if (offset_x > 0.5) or (w_tile == None):
            if (offset_y > 0.5) or (s_tile == None):
                candidate_tiles.extend( [e_tile, n_tile, ne_tile] )
            else:
                candidate_tiles.extend([e_tile, s_tile, se_tile])
        if (offset_x <= 0.5) or (e_tile == None):
            if (offset_y > 0.5) or (s_tile == None):
                candidate_tiles.extend([w_tile, n_tile, nw_tile])
            else:
                candidate_tiles.extend([w_tile, s_tile, sw_tile])

    # Get qubits from tiles
    candidates = []
    # Add qubits from candidate tiles        
    for t in candidate_tiles:
        if (t!=None):   candidates.extend( getQubits(t) )
        
    return candidates

def initSearch(source_cell):
    '''
     
    :param F:
    :param queue:
    :param source_cell:
    '''
    
    # initialize graph search
    nx.set_node_attributes(_RGraph, COST_THRESHOLD + 1, 'cost')
    nx.set_node_attributes(_RGraph, None, 'parent')
    nx.set_node_attributes(_RGraph, False, 'visited')

    _RGraph.nodes[source_cell]['level'] = 0
    _RGraph.nodes[source_cell]['cost'] = 0.0

def getCost(node, neighbor, source_cell, target_cell):
    '''
    
    :param node:
    :param neighbor:
    :param source_cell:
    :param target_cell:
    '''

    # Found target
    if (neighbor==target_cell):
        return 0.0

    cells = _RGraph.nodes[neighbor]['cells']
    # Qubit is assigned to target or source
    if (target_cell in cells):
        return 0.0
    if (source_cell in cells):
        return 0.0
    
    edge_cost = 0.0
    scope_cost = 0.0
    if target_cell:
        # Edge length cost
        path_len = _QCA.edges[source_cell,target_cell]['size']
        edge_cost = path_len / _max_edge
        # Scope cost
        node_tile = _RGraph.nodes[node]['tile']
        neighbor_tile = _RGraph.nodes[neighbor]['tile']   
        scope_cost = 0.0 if (node_tile==neighbor_tile) else BASE_B
    
    # Degree Cost
    degree_q = _RGraph.nodes[neighbor]['degree']
    degree_cost = (1 - (degree_q/MAX_DEG))
    
    # Base Cost (b_n)
    base_cost =  1 + BASE_A*(degree_cost) + scope_cost

    # Present-sharing Cost (p_n)    
    embedding = set()
    for path in _RGraph.nodes[neighbor]['path']:
        embedding.add( tuple(_RGraph.nodes[neighbor]['mapped'][path]) )
    k = len(embedding)
    sharing_cost = 1.0 + k * _alpha_p

    # History cost (h_n)    
    history_cost = _RGraph.nodes[neighbor]['history']
    
    # Node Cost (c_n = b_n * h_n * p_n)
    node_cost =  base_cost * sharing_cost * history_cost
    
    
    # Total cost (C_n) 
    if LENGTH_COST:
        cost = edge_cost*base_cost + (1.0-edge_cost)*node_cost
    else:
        cost = node_cost
    
    return cost

def BFS(source_cell, target_cell, queue):
    '''
    
    :param source_cell:
    :param target_cell:
    '''
     
    node = queue.pop()
    while (node != target_cell):
        
        for neighbor in _RGraph.neighbors(node):

            if not _RGraph.nodes[neighbor]['visited']:
                
                # Calculate cost of using node in path
                cost = getCost(node, neighbor, source_cell, target_cell)
                
                path_cost = _RGraph.nodes[node]['cost'] + cost
                
                if ( path_cost < _RGraph.nodes[neighbor]['cost'] ):
                    # Assign cost to node 
                    _RGraph.nodes[neighbor]['cost'] = path_cost
                    # Set path parent of node
                    _RGraph.nodes[neighbor]['parent'] = node
                    # Add to queue
                    queue.add(neighbor)
        
        # Set node as visited
        _RGraph.nodes[node]['visited'] = True
        # Get min cost node
        node = min(queue, key=lambda node: _RGraph.nodes[node]['cost'])
        # Remove node from queue 
        queue.remove(node)
        
    assert(node==target_cell)
    
def traceback(source_cell, target_cell):
    '''
    
    :param source_cell:
    :param target_cell:
    '''
    
    path = []
    node = _RGraph.nodes[target_cell]['parent']
    while (node!=source_cell):
        path.insert(0, node)
        node = _RGraph.nodes[node]['parent']
    
    return path
          
def populateQubits(path, source_cell, target_cell):
    '''
    
    :param path: [source_main, <virtual_cells>, target_main]
    :param source_cell:
    :param target_cell:
    '''
    
    global _RGraph, _QCA
    
    # Use sorted edge to have consistency in dictionaries
    edge = tuple(sorted([source_cell,target_cell]))

    # Write path to QCA graph
    _QCA.edges[source_cell,target_cell]['path'] = list(path)
    _QCA.edges[target_cell,source_cell]['path'] = list(path)[::-1]
    
    # Assign cells to main qubits
    source_main = path.pop(0)
    target_main = path.pop()

    # Add new edge to main qubit    
    _RGraph.nodes[source_main]['path'].append(edge)
    _RGraph.nodes[source_main]['mapped'][edge] = [source_cell]
    
    # Embed target cell in target main     
    _QCA.nodes[target_cell]['main']  =  target_main
    _QCA.nodes[target_cell]['qubits'].add(target_main)
    _RGraph.nodes[target_main]['path'].append(edge)
    _RGraph.nodes[target_main]['mapped'][edge] = [target_cell]
    _RGraph.nodes[target_main]['cells'].add(target_cell)
    
    ###################################################################################
    # Embed cells in chain qubits
    for qubit in path:
        # Look through path mappings to find if source or target are in qubit
        path_cell = None
        for pair in _RGraph.nodes[qubit]['path']:
            # Determine which cell (if any) is in path 
            mapping = set(_RGraph.nodes[qubit]['mapped'][pair])
            in_path = mapping.intersection(edge)
            
            # Modify qubit mapping according to cell in path
            if(in_path):
                path_cell = in_path.pop()
                remove_cell = [cell for cell in pair if cell!=path_cell][0]
                _RGraph.nodes[qubit]['mapped'][pair] = [path_cell]
                _RGraph.nodes[qubit]['mapped'][edge] = [path_cell]
                _QCA.nodes[remove_cell]['qubits'].discard(qubit)
                _RGraph.nodes[qubit]['cells'].discard(remove_cell)
                break
        
        # Neither source or target cell were in qubit    
        if (not path_cell):
            _RGraph.nodes[qubit]['mapped'][edge] = [source_cell, target_cell]
            _RGraph.nodes[qubit]['cells'].add(source_cell)
            _RGraph.nodes[qubit]['cells'].add(target_cell)
            _QCA.nodes[source_cell]['qubits'].add(qubit)
            _QCA.nodes[target_cell]['qubits'].add(qubit)

        # Add path to qubit            
        _RGraph.nodes[qubit]['path'].append(edge)
    ###################################################################################
    
def joinTargetCell(source_cell, target_cell, target_main):
    '''
    
    :param target_cell:
    :param target_main:
    '''

    global _RGraph
    
    if not target_main:
        candidates = getCandidates(target_cell)
        for qubit in candidates:
            if source_cell not in _RGraph.nodes[qubit]['cells']: 
                _RGraph.add_edge(qubit, target_cell)
    else:
        _RGraph.add_edge(target_main, target_cell)


def splitTargetCell(target_cell):
    global _RGraph
    
    to_remove = [n for n in _RGraph.in_edges(target_cell)]
    _RGraph.remove_edges_from(to_remove)      

def populateConflictedQubit(qubit, source_cell, target_cell):
    '''
    
    :param qubit:
    :param source_cell:
    :param target_cell:
    '''
    global _RGraph, _QCA
    
    # Use sorted edge to have consistency in dictionaries
    edge = tuple(sorted([source_cell,target_cell]))
    
    # Write path to corresponding QCA edge
    _QCA.edges[source_cell,target_cell]['path'] = [qubit]
    _QCA.edges[target_cell,source_cell]['path'] = [qubit]
    
    # Assign cells to main qubits
    _RGraph.nodes[qubit]['path'].append(edge)
    _RGraph.nodes[qubit]['mapped'][edge] = [source_cell, target_cell]
    

def expandPaths(source_cell, source_main):
    '''
    All the nodes on paths to previously found sinks are expanded automatically and added to the queue
    :param F:
    :param source_cell:
    :param target_cell:
    '''
    
    # set with intersection of qubits in paths
    queue = set()

    # Expand on main qubit
    _RGraph.nodes[source_main]['cost'] = 0.0
    _RGraph.nodes[source_main]['parent'] = source_cell
    queue.add(source_main)
    
    # Expand on qubits used in all paths    
    tree = [edge for edge in _RGraph.nodes[source_main]['path'] if source_cell in edge]
    for edge in tree:

        neighbor = [x for x in edge if x!=source_cell].pop()
        path = _QCA.edges[source_cell,neighbor]['path']
        parent = source_main

        for qubit in path[1:]:
            if (source_cell in _RGraph.nodes[qubit]['cells']):
                # Expand qubit
                _RGraph.nodes[qubit]['cost'] = 0.0
                _RGraph.nodes[qubit]['parent'] = parent
                # add to BFS queue
                queue.add(qubit)
                # Parent of next qubit
                parent = qubit
            else:
                break
    

    return queue

def routingTree(source_cell, targets, iteration):
    '''
    
    :param source_cell:
    :param targets:
    '''
    
    global _QCA, _RGraph
    
    for target_cell in targets:
        
        if VERBOSE: print('==============================Routing ' + str(source_cell) + ' ' + str(target_cell))

        source_main = _QCA.nodes[source_cell]['main']
        target_main = _QCA.nodes[target_cell]['main']

        if (source_main == target_main):
            if VERBOSE: print([source_main])
            # Populate conflicted qubit
            populateConflictedQubit(source_main, source_cell, target_cell)
            
        else:
            # Initialize qubit costs
            initSearch(source_cell)
            # Connect target cell to candidate qubits    
            joinTargetCell(source_cell, target_cell, target_main)
            # Initialize graph nodes and queue
            queue = expandPaths(source_cell, source_main)
            # Breadth First Search of target cell over routing graph
            BFS(source_cell, target_cell, queue)
            # Trace back the lowest cost path
            path = traceback(source_cell, target_cell)
            # Disconnect target cell from qubits
            splitTargetCell(target_cell)
    
            if VERBOSE: print(path)
            
            # Popoluate qubits
            populateQubits(path, source_cell, target_cell)
    
        
        # Mark edge as routed
        _QCA.edges[source_cell,target_cell]['routed'] =  True
    
def updateHistoryCost():
    '''
    
    :param F: Routing graph
    :param conflicts:
    '''
    
    global _RGraph
    
    conflict_qubits = set()
    conflict_cells = set()
    conflict_paths = set()
    
    conflicts = 0
    occ_qubits = 0
    # Count qubits with conflicts
    for qubit in _Chimera:
        embedded = set()
        paths = _RGraph.nodes[qubit]['path']
        for path in paths:
            cells = tuple(_RGraph.nodes[qubit]['mapped'][path])
            embedded.add(cells)
        
        num_cells = len(embedded)    
        if (num_cells>1):
            conflicts = conflicts+1 
            conflict_qubits.add(qubit)
            for mapping in embedded:
                for cell in [mapping]:
                    conflict_cells.add(cell)
            conflict_paths.update(paths)
        
        if (num_cells>0):    
            occ_qubits = occ_qubits + 1
            _RGraph.nodes[qubit]['history'] = _RGraph.nodes[qubit]['history'] + (DELTA_H*num_cells)    
        

    
    if VERBOSE:
        print('######################## OCCUPIED QUBITS')
        print(occ_qubits)
        if conflict_qubits:
            print('######################## CONFLICTS')
            for qubit in conflict_qubits:
                print(str(qubit) + ' ' + str([x for x in _RGraph.nodes[qubit]['cells']]))

            
    legal = not bool(conflicts)
    return legal


def getTargets(source_cell):
    
    targets = []
    
    for neighbor in _QCA.neighbors(source_cell):
        if not _QCA.edges[neighbor,source_cell]['routed']:
            targets.append(neighbor)
    
    return targets


def measureChains():
    '''
    
    '''
    global _max_edge
    
    for cell in _QCA:
        size = len(_QCA.nodes[cell]['qubits'])
        _QCA.nodes[cell]['size'] = size
        _QCA.nodes[cell]['priority'] = size
    
    _max_edge = 0.0    
    for edge in _QCA.edges():
        u,v = edge
        path_len = float(len(_QCA.edges[u,v]['path']))
        _QCA.edges[u,v]['size'] = path_len 
        _QCA.edges[v,u]['size'] = path_len 
        if path_len > _max_edge:
            _max_edge = path_len

def embedFirst(cell):
    
    candidates = getCandidates(cell)
    if RANDOMIZE_CANDIDATES:
        random.Random(SEED).shuffle(candidates)
    
    qubit = min(candidates , key=lambda candidate: getCost(cell, candidate, cell, None))
    
    _QCA.nodes[cell]['main']  =  qubit
    _QCA.nodes[cell]['qubits'].add(qubit)
    _RGraph.nodes[qubit]['cells'].add(cell)


def sortCells(cells, randomize, priority):
        
    cells_list = list(cells)
    
    if randomize:
        random.Random(SEED).shuffle(cells_list)
    
    if priority:
        return sorted(cells_list, key=lambda key: _QCA.nodes[key]['priority'], reverse=True)
    else:
        return cells_list


def negotiatedCongestion():
    '''
    
    '''
    global _alpha_p, _QCA
    
    # Router statistics
    stats = {}
    
    # Initialize present-sharing congestion cost
    _alpha_p = 0.0

    # Initialize congestion
    legal = False; itry = 0
    while ((not legal) and (itry < 100)):
        
        if VERBOSE: print('########### ROUTER ITERATION: ' + str(itry))
        
        ripUpAll()
        
        # Initialize cell queue
        sorted_cells = sortCells(_QCA.nodes(), False, True)
        firstCell = sorted_cells[0]
        embedFirst(firstCell)
        
        # Traverse cells        
        cell_queue = set([firstCell])
        while cell_queue:
        
            # Get next source cell to route
            source_cell = sortCells(cell_queue, RANDOMIZE_CELLS, LENGTH_PRIORITY)[0]
            #source_cell = list(cell_queue)[0]
            # Sinks are QCA neighbours        
            targets = getTargets(source_cell)
            
            # Get Routing Tree dictionary
            routingTree(source_cell, targets, itry)
            
            # Mark routed cell
            _QCA.nodes[source_cell]['routed'] =  True
            cell_queue.remove(source_cell)
            # Next cell
            cell_queue.update([cell for cell in targets if not _QCA.nodes[cell]['routed']])

        # Determine conflict costs and legality        
        legal = updateHistoryCost()
        # Measure cell chains
        measureChains()
        # Increase present-sharing cost
        _alpha_p = _alpha_p + DELTA_P
        # Increase iteration counter
        itry = itry + 1
    
    # Router result statistics
    stats['ROUTER_ITERATIONS'] = itry
    
    return legal, stats

def parseConfiguration(configuration):
    
    global M, N, L, SEED
     
    M, N, L = configuration['M'], configuration['N'], configuration['L']
    SEED = configuration['SEED']
    
    global VERBOSE, PLOT
    
    PLOT    = configuration['PLOT']
    VERBOSE = configuration['VERBOSE']
    
    global BASE_A, BASE_B, BASE_C, DELTA_P, DELTA_H
    
    BASE_A =       configuration['routing']['BASE_A']
    BASE_B =       configuration['routing']['BASE_B']
    BASE_C =       configuration['routing']['BASE_C']
    DELTA_P =      configuration['routing']['DELTA_P'] 
    DELTA_H =      configuration['routing']['DELTA_H']
    
    global SIMPLE, LENGTH_COST, NEIGHBORHOOD, LENGTH_PRIORITY, RANDOMIZE_CELLS, RANDOMIZE_CANDIDATES  
    
    SIMPLE =                configuration['routing']['SIMPLE']
    LENGTH_COST =           configuration['routing']['LENGTH_COST']
    NEIGHBORHOOD =          configuration['routing']['NEIGHBORHOOD']
    LENGTH_PRIORITY =       configuration['routing']['LENGTH_PRIORITY']
    RANDOMIZE_CELLS =       configuration['routing']['RANDOMIZE_CELLS']
    RANDOMIZE_CANDIDATES =  configuration['routing']['RANDOMIZE_CANDIDATES']

def legalize(Chimera, QCA, bins, tiles, configuration):
    
    global _Chimera, _QCA, _RGraph
    
    # Parse routing parameters
    parseConfiguration(configuration)
    
    # Parse routing structures 
    _Chimera = Chimera.to_directed() 
    _QCA = QCA
    
    # Create routing graph
    setupRGraph()
    
    # Negotiated congestion router
    good, stats = negotiatedCongestion()
    
    # Split chains to minimize max chain length 
    if good:
        solveChains(_RGraph, _QCA, configuration)

    return good, stats