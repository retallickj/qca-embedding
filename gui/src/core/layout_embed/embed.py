'''
Created on Jan 23, 2017

@author: JosePinilla
'''

import os
import sys
import json
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

##### DIFFUSION
from diffusion import diffusion
##### ROUTING
from traverse_router import legalize


VERBOSE =  False
PLOT = False

# Chimera and Dimensions
Chimera = None
N = 12
M = 12
L = 4

# Problem and Dimensions
QCA = None
X = 0
Y = 0
W = 1000
H = 1000

DIFFUSION = True
CIRCUIT = None

def plotLocations(QCA):
    pos = {}
    for c in QCA:
        cell = QCA.node[c]['cell']
        pos[c] = (cell['x'] , cell['y'])
    plt.figure(0)
    plt.clf()
    nx.draw(QCA, pos=pos, with_labels=True)
    
    plt.grid('on')
    plt.axis('on')
    plt.axis([0,N,0,M])
    x_ticks = np.arange(0, N) # steps are width/width = 1 without scaling
    y_ticks = np.arange(0, M)
    plt.xticks(x_ticks)
    plt.yticks(y_ticks)
    plt.gca().invert_yaxis()
    if WRITE:
        if not os.path.exists('./plot/'):
            os.makedirs('./plot/')
        plt.savefig('./plot/plt_' + CIRCUIT + '0.png')
    if PLOT:
        plt.show()


def partitionScale(QCA):
    '''
    Scale the QCA layout to the dimensions of the Chimera
    and assign cells to partitions/tiles/regions. Store new
    scaled layout locations
    :param QCA: QCA graph
    :param num_regions_x: columns of tiles/regions
    :param num_regions_y: rows of tiles/regions
    '''
    
    num_regions_x, num_regions_y = N, M
    
    # Cell container for each tile/region
    num_regions = num_regions_x*num_regions_y
    bins = {key : [] for key in xrange( 0,  num_regions)}
    
    # Size of a region
    qca_w_scale = W / float(num_regions_x)
    qca_h_scale = H / float(num_regions_y)  

    # Assign partitions and scale to chimera size
    for node in QCA:
        # Cell data from graph
        cell = QCA.nodes[node]['cell']
        
        # Scale
        scaled_x = (cell['x'] - X) / qca_w_scale
        scaled_y = (cell['y'] - Y) / qca_h_scale
        cell['x'] = scaled_x
        cell['y'] = scaled_y
        
        # Assign tile/region
        cell_region_x = int(np.floor( scaled_x ))
        cell_region_y = int(np.floor( scaled_y ))
        region = cell_region_x + cell_region_y*num_regions_x
        QCA.node[node]['tile'] = region
        bins[region].append(node)
    
    # Plot QCA Graph   
    if PLOT or WRITE: plotLocations(QCA)
        
    return bins

def layoutToModels(cell_map):
    
    models = {}
    
    for node in cell_map:
         
        if VERBOSE:
            print ("CELL " + str(node))
            print (cell_map[node]['qubits'])
         
        models[int(node[1:])] = cell_map[node]['qubits']
    
    return models

def getSupply(Chimera):
    
    supply = {}
    
    for m in xrange(M):
        for n in xrange(N):
            tile = n + N*m
            supply[tile] = []
            for h in xrange(2):
                for i in xrange(L):
                    qubit_tuple = (m,n,h,i)
                    if qubit_tuple in Chimera:
                        if Chimera.neighbors(qubit_tuple):
                            Chimera.node[qubit_tuple]['alive'] = True
                            supply[tile].append(qubit_tuple) 
                        else: 
                            Chimera.node[qubit_tuple]['alive'] = False
    return supply
        

def defaultConf ():
    DEFAULT_CONF = {}
    
    DEFAULT_CONF['PLOT'] = False
    DEFAULT_CONF['STATS'] = True
    DEFAULT_CONF['WRITE'] = False
    DEFAULT_CONF['VERIFY'] = False
    DEFAULT_CONF['VERBOSE'] = False
    DEFAULT_CONF['WRITELP'] = False
    
    DEFAULT_CONF['SEED'] =                              None
    
    DEFAULT_CONF['D_MAX'] =                             8.0
    DEFAULT_CONF['MAX_DEG'] =                           6.0
    
    DEFAULT_CONF['diffusion'] = {}
    DEFAULT_CONF['diffusion']['ENABLE'] =               True
    DEFAULT_CONF['diffusion']['C_LIM'] =                8.0
    DEFAULT_CONF['diffusion']['DELTA_T'] =              0.3
    DEFAULT_CONF['diffusion']['D_SCALER'] =             3.0
    DEFAULT_CONF['diffusion']['DIAGONAL'] =             0.5
    DEFAULT_CONF['diffusion']['VARIANCE_GROUP'] =       3.0
    DEFAULT_CONF['diffusion']['VARIANCE_THR'] =         0.01
     
    DEFAULT_CONF['routing'] = {}

    DEFAULT_CONF['routing']['SIMPLE'] =                 False
    DEFAULT_CONF['routing']['LENGTH_COST'] =            True
    DEFAULT_CONF['routing']['NEIGHBORHOOD'] =           '2-block'
    DEFAULT_CONF['routing']['LENGTH_PRIORITY'] =        True
    DEFAULT_CONF['routing']['RANDOMIZE_CELLS'] =        True
    DEFAULT_CONF['routing']['RANDOMIZE_CANDIDATES'] =   False

    DEFAULT_CONF['routing']['BASE_A'] =                 1.0
    DEFAULT_CONF['routing']['BASE_B'] =                 0.2
    DEFAULT_CONF['routing']['BASE_C'] =                 0.0
    DEFAULT_CONF['routing']['DELTA_H'] =                0.1
    DEFAULT_CONF['routing']['DELTA_P'] =                0.45
    
    return DEFAULT_CONF

def layoutConfiguration(conf, test_conf={}):

    TEMP_CONF = defaultConf()
    
    # Import class I and II dictionaries into configuration
    for k_I, v_I in test_conf.items():
        if not isinstance(v_I, dict):
            TEMP_CONF.update({ k_I : v_I })
        else:
            if k_I in TEMP_CONF:
                for k_II, v_II in v_I.items():
                    TEMP_CONF[k_I].update({ k_II : v_II })
                
    conf.update(TEMP_CONF)

def setProblem(problem_adj, nodes_loc, spacing):
    '''
    Create graph from adjacency and attributes of input problem.
    :param problem_adj: Adjacency dictionary
    :param node_loc: Attribute dictionary including 'x' and 'y' node location attributes.
    :param spacing: Diameter or size of input problem node
    '''
    
    global QCA
    global X, Y, W, H

    QCA = nx.DiGraph()

    qca_w1 = sys.maxint
    qca_w2 = 0 
    qca_h1 = sys.maxint
    qca_h2 = 0

    # Iterate only through active cells
    for cell in problem_adj:
        cell_x = nodes_loc[cell].x
        cell_y = nodes_loc[cell].y
        degree =  len(problem_adj[cell])
        
        cell_dict = {'x' : cell_x, 'y' : cell_y, 'degree' : degree}
        QCA.add_node(cell,cell=cell_dict)
        
        edges_in    = [(cell,neigh) for neigh in problem_adj[cell]]
        edges_out   = [(neigh,cell) for neigh in problem_adj[cell]]
        
        QCA.add_edges_from(edges_in)
        QCA.add_edges_from(edges_out)
        
        if ( cell_x < qca_w1):
            qca_w1 = cell_x
        elif (cell_x > qca_w2):
            qca_w2 = cell_x
        if ( cell_y < qca_h1):
            qca_h1 = cell_y
        elif ( cell_y > qca_h2):
            qca_h2 = cell_y
          

    X = qca_w1-(spacing/2)
    Y = qca_h1-(spacing/2)
    W = qca_w2-qca_w1+(spacing)
    H = qca_h2-qca_h1+(spacing)
    
def setTarget(chimera_adj):
    '''
    Create the target graph for minor-embedding.
    :param chimera_adj: Adjacenecy matrix of Chimera graph. TODO: Use DWave_networkx
    '''
    
    global Chimera
    
    # Create Chimera graph structure 
    Chimera = nx.Graph(chimera_adj)
    

def parseConfiguration(configuration):
    '''
    
    :param configuration: Configuration dictionary for embedding algorithm
    '''
    
    global M, N, L
     
    M, N, L = configuration['M'], configuration['N'], configuration['L']
    
    global VERBOSE, PLOT, WRITE
    
    VERBOSE =   configuration['VERBOSE']
    PLOT =      configuration['PLOT']
    WRITE =     configuration['WRITE']
    
    global CIRCUIT, DIFFUSION
    
    CIRCUIT = configuration['CIRCUIT']
    DIFFUSION = configuration['diffusion']['ENABLE']
    
def save_problem_file(filename, problem_adj, node_loc):
    '''
    Store problem adjacency and attributes in JSON format
    :param filename: JSON file to store input problem
    :param problem_adj: Adjacency dictionary
    :param node_loc: Attribute dictionary including 'x' and 'y' node location attributes.
    '''
    
    with open(filename, 'a') as data_file:
        
        node_loc_dict = {}
        for node in node_loc:
            node_x = node_loc[node].x
            node_y = node_loc[node].y
            node_dict = {'x':node_x, 'y':node_y}
            node_loc_dict[node] = node_dict 
            

        json.dump(problem_adj, data_file)
        data_file.write('\n')
        json.dump(node_loc_dict, data_file)    

def layoutEmbed(configuration, stats):
    '''
    Layout-aware embedding 
    :param configuration: Input layout embedding configuration dictionary 
    :param stats: Output statistics dictionary of embedding attempt
    '''
        
    parseConfiguration(configuration)
    
    # Parse Chimera to collect alive qubits per tiles
    tiles = getSupply(Chimera)
    
    # QCA graph layout partition into bins
    if VERBOSE: print("Partition")
    bins = partitionScale(QCA)
    
    # Diffusion to make a denser clustering
    if VERBOSE: print("Diffusion")
    if DIFFUSION:
        stats_diffusion = diffusion(Chimera, QCA, bins, tiles, configuration)
        stats.update(stats_diffusion)

    # Return true if embedding was successful
    if VERBOSE: print("Routing")    
    good, stats_legal = legalize(Chimera, QCA, bins, tiles, configuration)
    stats.update(stats_legal) 
    
    cell_map =  QCA.node
    
    return good, cell_map
