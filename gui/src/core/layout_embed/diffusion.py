'''
Created on Mar 27, 2017

@author: JosePinilla
'''
import numpy as np
import networkx as nx
from math import floor, sqrt
import matplotlib.pyplot as plt

M = None
N = None
L = None
C_LIM = None            # Concentration limit
D_QCA = None            # Diffusivity
S_MAX =  None           # Maximum supply
DIAGONAL = None


DELTA_T = None          # time step size of diffusion
VARIANCE_THR = None      # target variance
VARIANCE_GROUP = None    # number of iterations used to calculate variance
D_SCALER = None          # value to scale the density metric (# of cells * scaler) / (# of qubits)
CIRCUIT = None


PLOT = True
VERBOSE = False
WRITE = True


_QCA = None
_bins = None
_concentration =  None
_supply = None

def plotLocations(iteration):
    pos = {}
    for c in _QCA:
        cell = _QCA.nodes[c]['cell']
        pos[c] = (cell['x'] , cell['y'])
    plt.figure(0)
    plt.clf()
    nx.draw(_QCA, pos=pos, with_labels=True)

    plt.grid('on')
    plt.axis('on')
    plt.axis([0,N,0,M])
    x_ticks = np.arange(0, N) # steps are width/width = 1 without scaling
    y_ticks = np.arange(0, M)
    plt.xticks(x_ticks)
    plt.yticks(y_ticks)
    plt.gca().invert_yaxis()
    if WRITE:
        plt.savefig('./plot/plt_' + CIRCUIT + str(iteration) + '.png')
    if PLOT:
        plt.show()


def partition(iteration):

    global _concentration


    regions = range( 0, N*M )

    # Plot QCA Graph
    if PLOT or WRITE: plotLocations(iteration)


    concentration = dict.fromkeys( regions, 0.0 )

    overpopulated = False
    for node in _QCA:
        cell = _QCA.nodes[node]['cell']

        # assign cell to region
        cell_tile_x = int(floor(cell['x'] ))
        cell_tile_y = int(floor(cell['y'] ))
        tile = cell_tile_x + cell_tile_y*N
        _QCA.nodes[node]['tile'] = tile
        _bins[tile].append(node)
        C_tile = (concentration[tile] + 1.0)
        concentration[tile] =  C_tile
        S_tile = _supply[tile]
        if C_tile/S_tile > C_LIM/S_MAX:
            overpopulated = True


    _concentration =  {k: v for k, v in concentration.items()}
    # for boundaries
    _concentration[-1] = C_LIM

    return overpopulated

def printConcentration():

    print('##################CONCENTRATION')
    for j in range(M):
        print([_concentration[i+j*N]  for i in  range(N)])

def measureDispersion():
    '''
    Measure dispersion
    :param QCA:
    '''

    dist_accum = 0.0

    center_x = N/2.0
    center_y = M/2.0

    for node in _QCA:
        x1 = _QCA.nodes[node]['cell']['x']

        y1 = _QCA.nodes[node]['cell']['y']

        dist_accum = dist_accum + pow(x1-center_x,2) + pow(y1-center_y,2)

    disperse = dist_accum / len(_QCA)

    if VERBOSE:
        print("DISPERSE:" + str(disperse))

    return disperse


def measureSparsity():
    '''
    Measure Sparsity
    :param concentration:
    '''
    dist_accum = 0

    for edge in _QCA.edges():
        x1 = _QCA.nodes[edge[0]]['cell']['x']
        x2 = _QCA.nodes[edge[1]]['cell']['x']
        y1 = _QCA.nodes[edge[0]]['cell']['y']
        y2 = _QCA.nodes[edge[1]]['cell']['y']

        dist_accum = dist_accum + pow(x1-x2,2) + pow(y1-y2,2)

    sparsity = dist_accum #* max_n

    if VERBOSE:
        print("SPARSITY:" + str(sparsity))

    return sparsity

def isSpread(history):

    ###### Calculate variance

    # Get mean value
    mean = sum(history) / VARIANCE_GROUP

    # Identify trend
    increasing = True
    diff_accum = 0.0
    prev_val = 0.0
    for value in history:
        sq_diff = pow(value - mean, 2)
        diff_accum = diff_accum + sq_diff
        if (value < prev_val):
            increasing = False
        prev_val = value

    variance =  (diff_accum/VARIANCE_GROUP)
    std_dev = sqrt(variance)


    spread = (std_dev/mean > VARIANCE_THR)

    return spread, increasing

def neighbourTiles(tile):
    tile_x = tile % N
    tile_y = (tile-tile_x) / N

    n_tile = (tile - N)     if (tile_y > 0)      else   -1
    s_tile = (tile + N)     if (tile_y < M-1)    else   -1
    w_tile = (tile - 1)     if (tile_x > 0)      else   -1
    e_tile = (tile + 1)     if (tile_x < N-1)    else   -1

    nw_tile = (tile - N - 1)  if (tile_y > 0    and tile_x > 0)    else -1
    ne_tile = (tile - N + 1)  if (tile_y > 0    and tile_x < N-1)  else -1
    se_tile = (tile + N + 1)  if (tile_y < M-1  and tile_x < N-1)  else -1
    sw_tile = (tile + N - 1)  if (tile_y < M-1  and tile_x > 0)    else -1


    return n_tile, s_tile, w_tile, e_tile, nw_tile, ne_tile, se_tile, sw_tile


def layoutCost():

    #return measureSparsity()
    return measureDispersion()


def velocityVectors(D):
    V = {}
    num_tiles = N*M
    for tile in range(num_tiles):

        # Tile concentration
        C_tile   = _concentration[tile]
        if (C_tile==0): continue

        # Neighboring tiles
        n_tile, s_tile, w_tile, e_tile,  nw_tile, ne_tile, se_tile, sw_tile = neighbourTiles(tile)

        D_tile =  D[tile]

        Dn_tile = D[n_tile]
        Ds_tile = D[s_tile]
        Dw_tile = D[w_tile]
        De_tile = D[e_tile]

        Dnw_tile = D[nw_tile]
        Dne_tile = D[ne_tile]
        Dse_tile = D[se_tile]
        Dsw_tile = D[sw_tile]

        D_E = De_tile + (Dne_tile + Dse_tile)*DIAGONAL
        D_W = Dw_tile + (Dnw_tile + Dsw_tile)*DIAGONAL
        D_N = Dn_tile + (Dnw_tile + Dne_tile)*DIAGONAL
        D_S = Ds_tile + (Dsw_tile + Dse_tile)*DIAGONAL

        V_H = -(D_E - D_W) / (2.0*D_tile)

        V_V = -(D_S - D_N) / (2.0*D_tile)

        V[tile] = (V_H, V_V)

    V[-1] = 0.0
    return V


def getAttractors(tile, n,m):
    '''
    Get attractor tiles
    :param tile:
    '''

    n_tile, s_tile, w_tile, e_tile, nw_tile, ne_tile, se_tile, sw_tile = neighbourTiles(tile)

    dh = n - (N/2)

    dv = m - (M/2)

    # if cell is to the right of center
    if (dh > 0.0):
        attr_tile_x = w_tile
        # if cell is below center
        if (dv < 0.0):
            attr_tile_y = s_tile
            attr_tile_xy = sw_tile
        else:
            attr_tile_y = n_tile
            attr_tile_xy = nw_tile
    # if cell is to the right of center
    else:
        attr_tile_x = e_tile
        # if cell is below center
        if (dv < 0.0):
            attr_tile_y = s_tile
            attr_tile_xy = se_tile
        else:
            attr_tile_y = n_tile
            attr_tile_xy = ne_tile

    return attr_tile_x, attr_tile_y, attr_tile_xy

def moveCells(D):

    for m in range(M):
        for n in range(N):

            tile = n + m*N
            cells =  _bins[tile]
            attr_tile_x, attr_tile_y, attr_tile_xy = getAttractors(tile, n, m)

            for cell in cells:

                D_tile = D[tile]

                D_attr_x  = D[attr_tile_x]
                D_attr_y  = D[attr_tile_y]
                D_attr_xy = D[attr_tile_xy]

                N_x = - ( ( C_LIM / S_MAX ) - ( D_attr_x + (D_attr_xy/2.0) ) ) / (2.0*D_tile)

                N_y = - ( ( C_LIM / S_MAX ) - ( D_attr_y + (D_attr_xy/2.0) ) ) / (2.0*D_tile)

                c_x = _QCA.nodes[cell]['cell']['x']
                c_y = _QCA.nodes[cell]['cell']['y']

                lx_cell = ((2.0 * c_x) / N ) - 1
                ly_cell = ((2.0 * c_y) / M ) - 1

                ################################################################################
                ################ NEW LOCATION
                ################################################################################

                ###### X velocity
                v_h = N_x * lx_cell * (1-D_QCA)

                delta_x =  v_h*DELTA_T
                c_x =  c_x + delta_x
                # cap X value inside tiles
                c_x = min(c_x,float(N)-0.0001)
                c_x = max(0,c_x)

                ###### Y velocity
                v_v = N_y * ly_cell * (1-D_QCA)

                delta_y =  v_v*DELTA_T
                c_y =  c_y + delta_y
                # cap Y value inside tiles
                c_y = min(c_y,float(M)-0.0001)
                c_y = max(0,c_y)

                _QCA.nodes[cell]['cell']['x'] = c_x
                _QCA.nodes[cell]['cell']['y'] = c_y

            # Reset cells in bin until next partition
            _bins[tile][:] = []

def getDensities():
    '''

    '''

    D = {}
    C_accum = 0.0
    D_accum = 0.0
    occupancy = 0
    num_tiles = N*M

    for tile in range(num_tiles):

        C_tile = _concentration[tile]
        S_tile = _supply[tile]
        if C_tile:
            occupancy =  occupancy + 1
            if S_tile:
                D_tile = C_tile/S_tile
                D[tile] =  D_tile
            else:
                D_tile = -1
                D[tile] = -1

            # Global metrics
            C_accum = C_accum + C_tile
            D_accum = D_accum + D_tile
        else:
            D[tile] = 0


    # Boundary condition
    D[-1] = 1
    #Global metrics
    C_avg = C_accum / occupancy
    D_avg = D_accum / occupancy

    return D, C_avg, D_avg

def parseConfiguraion(configuration):

    global M, N, L

    M, N, L = configuration['M'], configuration['N'], configuration['L']

    global CIRCUIT

    CIRCUIT = configuration['CIRCUIT']

    global PLOT, VERBOSE, WRITE

    PLOT =          configuration['PLOT']
    VERBOSE =       configuration['VERBOSE']
    WRITE =         configuration['WRITE']

    global C_LIM, DELTA_T, D_SCALER, DIAGONAL, VARIANCE_GROUP, VARIANCE_THR

    C_LIM =             configuration['diffusion']['C_LIM']
    DELTA_T =           configuration['diffusion']['DELTA_T']
    D_SCALER =          configuration['diffusion']['D_SCALER']
    DIAGONAL =          configuration['diffusion']['DIAGONAL']
    VARIANCE_GROUP =    configuration['diffusion']['VARIANCE_GROUP']
    VARIANCE_THR =      configuration['diffusion']['VARIANCE_THR']


def diffusion(QCA, bins, tiles, configuration):

    global _QCA, _bins, _concentration, _supply
    global D_QCA, S_MAX

    # Parse diffusion parameters
    parseConfiguraion(configuration)

    # Parse diffusion structures
    _QCA =  QCA
    _bins =  bins
    _concentration = {tile:len(bins[tile]) for tile in bins}
    _concentration[-1] = 1

    # Set target density (D_QCA)
    D_QCA = min( (len(_QCA)*D_SCALER) / (N*M*L*2.0), 1.0)
    S_MAX = L*2.0
    
    _supply = {tile:float(len(tiles[tile])) for tile in tiles}
    _supply[-1] = 1
    # Measure initial cost
    cost = layoutCost()

    cost_history = [cost] * int(VARIANCE_GROUP)

    i = 0
    diffuse = True
    while (diffuse and (i<100)):

        # Density Matrix
        if VERBOSE: printConcentration()

        ########## Refresh tile densities
        D, C_avg, D_avg = getDensities()

        ########## Calculate cell locations
        moveCells(D)

        # Keep count of iterations
        i = i + 1

        ########## Measure costs
        overpopulated = partition(i)
        cost = layoutCost()

        # FIFO
        cost_history.pop(0)
        cost_history.append(cost)
        spread, increasing = isSpread(cost_history)

        ########## Diffusion condition
        diffuse = (overpopulated or spread) and not increasing

    stats = {}
    stats['DIFFUSION_ITERATIONS'] = i
    stats['AVG_CONCENTRATION'] =  C_avg
    stats['AVG_DENSITY'] =  D_avg

    return stats
