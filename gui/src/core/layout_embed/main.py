'''
Created on Oct 24, 2017

@author: JosePinilla
'''
######################################################################################################
######################################################################################################
########################################## STANDALONE ################################################ 
import sys
import json
import traceback
import numpy as np
import dwave_networkx as dnx
from embed import layoutEmbed, layoutToModels
from embed import layoutConfiguration, setProblem, setTarget

M = 8
N = 8
L = 4
SPACING = 18.0 

class JSONDecoder(json.JSONDecoder):
    '''
    Decode JSON strings into components
    '''
    def decode(self, s):
        result = super(JSONDecoder, self).decode(s)
        return self._decode(result)

    def _decode(self, o):
        try:
            isUnicode = isinstance(o, unicode) 
        except:
            isUnicode = isinstance(o, str) # Python3 renamed unicode to str
        if isUnicode: 
            try:
                return int(o)
            except ValueError:
                return o
        elif isinstance(o, dict):
            return {self._decode(k): self._decode(v) for k, v in o.items()}
        elif isinstance(o, list):
            return [self._decode(v) for v in o]
        else:
            return o

class Node():
    def __init__(self,x=0,y=0):
        self.x = x
        self.y = y


def load_problem_file(filename, spacing):
    '''
    Reads in input problem description
    
    :param filename: JSON input problem file with adjacency dictionary and attributes dictionary
        {<node>:<neighbors>}
        {<node>:{'x':X, 'y':Y, ...}} 
    :param spacing: Diameter or size of input problem node
    :returns problem_adj: Adjacency dictionary
    :returns node_loc: Attribute dictionary including 'x' and 'y' node location attributes.
    '''

    try:
        with open(filename, 'r') as data_file:    
            problem_adj, node_dict = (json.loads(line,cls=JSONDecoder) for line in data_file)
    except:
        print('Failed to open file: {0}'.format(filename))
        raise IOError
    
    node_loc = {}
    for node in node_dict:
        node_x = int(node_dict[node]['x'])
        node_y = int(node_dict[node]['y'])
        node_obj = Node(node_x,node_y)
        node_loc[node] = node_obj
        
    return problem_adj, node_loc
        
def linear_to_tuple(qubit,N,M):
    
    qpr = 2*N*L     # qbits per row
    qpt = 2*L       # qbits per tile

    qubit -= 1
    row, rem = divmod(qubit, qpr)
    col, rem = divmod(rem, qpt)
    horiz, ind = divmod(rem, L)

    qubit_tuple = (row, col, horiz, ind) 

    return qubit_tuple


def load_chimera_file(filename):
    '''Load a chimera graph from an edge specification file'''
    
    try:
        fp = open(filename, 'r')
    except:
        print('Failed to open file: {0}'.format(filename))
        raise IOError
        
    # Get number of qubits and number of connectors
    num_qbits, num_conns = [int(x) for x in fp.readline().split()]
    
    # Get dimensions of tiles
    # Allows nxm, n!=m for testing
    #M, N = [int(x) for x in fp.readline().split()]
    
    adj = {i: [] for i in range(1, num_qbits+1)}
                    
    for line in fp:
        a, b = [int(x) for x in line.strip().split()]
        adj[a].append(b)
        adj[b].append(a)
    
    # Processor size
    # LEGACY: only allows m*m processors.
    M = int(np.sqrt(num_qbits/(2*L)))
    N = M
    
    adj_tup = {linear_to_tuple(k,N,M):\
            [linear_to_tuple(v,N,M) for v in adj[k]] for k in adj}
    
    fp.close()

    return M, N, adj_tup


def main():
    '''Setup and run the Layout-Aware Placement algorithm'''
    global M,N,L
    global X,Y,W,H 

    if len( sys.argv) == 3:
        problem_file = sys.argv[1]
        chimera_file = sys.argv[2]
    else:
        raise ValueError("main.py <problem_file> <chimera_file>")

    problem_adj, node_loc =  load_problem_file(problem_file, SPACING)
    M, N, chimera_adj = load_chimera_file(chimera_file)
    
    stats = {}
    configuration = {}
    configuration['M'] = M
    configuration['N'] = N
    configuration['L'] = L
    configuration['CIRCUIT'] = 'problem'
    layoutConfiguration(configuration)
    
    try:
        setProblem(problem_adj, node_loc, SPACING)
        setTarget(chimera_adj)
        good, cell_map = layoutEmbed(configuration, stats)
        if good:
            print('Layout-Aware Embedding Successful')
    except Exception as e:
        good = False
        if type(e).__name__ == 'KeyboardInterrupt':
            raise KeyboardInterrupt
        print('Layout-Aware Embedding Failed')
        print (traceback.print_exc())
        return
    
    print( layoutToModels(cell_map) ) 

if __name__ == '__main__':
    main()