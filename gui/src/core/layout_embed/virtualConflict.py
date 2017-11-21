'''
Created on Mar 11, 2017

@author: JosePinilla
'''

import pulp

WRITE = False
VERBOSE = False
PLOT = False


_RGraph = None
_QCA = None

def parseSystem():
    '''
    
    :param RGraph: Routing graph
    :param QCA: QCA graph
    :param paths: dictionary of resulting paths from routing. key  = QCA edge : value = qubits in path  
    '''
    
    shared = {}
    fixed = {}
    sharing = False
    for cell in _QCA:
        fixed[cell] = 0
        shared[cell] = set()
        cell_qubits = _QCA.nodes[cell]['qubits'] 
        for qubit in cell_qubits:
            qubit_shared = len(_RGraph.nodes[qubit]['cells']) > 1
            if(not qubit_shared):
                fixed[cell] = fixed[cell] + 1 

    
    chains = {}
    chain_lengths = {}
    
    for edge in _QCA.edges():
        chains[edge] = []
        chain_lengths[edge] = 0
        
        for qubit in _QCA.edges[edge]['path']:
            qubit_shared = len(_RGraph.nodes[qubit]['cells']) > 1 
            if(qubit_shared):
                chains[edge].append(qubit)
                chain_lengths[edge] = chain_lengths[edge] + 1
                 
                cell_S = edge[0]
                cell_T = edge[1]
                var_S = "Q"+str(cell_S)+str(cell_S)+str(cell_T)
                var_T = "Q"+str(cell_T)+str(cell_S)+str(cell_T)
                
                shared[cell_S].add(var_S)
                shared[cell_T].add(var_T)
                sharing =  True        

    return fixed, shared, chains, chain_lengths, sharing
                
                
def setupProblem(fixed, shared, chains, chain_lengths):
    '''
    
    :param RGraph: Routing graph
    :param fixed: key = cell, value = number of fixed qubits assigned to cell 
    :param shared: key = cell, value = set of shared qubits assigned to cell
    :param chains:
    :param chain_lengths:
    '''
    
    prob = pulp.LpProblem("Solve Chains",pulp.LpMinimize)



    # OBJECTIVE FUNCTION
    Z = pulp.LpVariable('Z',lowBound=0,cat='Integer')
    prob += Z, "OBJ"
    
    # VARIABLES
    Lpvars = {}
    var_map = {}
    # Each variable represents the number of virtual
    # qubits corresponding to the cells assigned to the chain
    for chain in chains:
        if chains[chain]:
            # Name of variables are Q<source><edge> and Q<target><edge>
            cell_S = chain[0]
            cell_T = chain[1]
            var_S = "Q"+str(cell_S)+str(cell_S)+str(cell_T)
            var_T = "Q"+str(cell_T)+str(cell_S)+str(cell_T)

            # create variables
            Lpvars[var_S] = pulp.LpVariable(var_S, lowBound=0, cat='Integer')
            Lpvars[var_T] = pulp.LpVariable(var_T, lowBound=0, cat='Integer')
    
            var_map[chain] = {}
    
            # mapping from path and cell to variable name for post-solve parsing
            var_map[chain][cell_S] = var_S
            var_map[chain][cell_T] = var_T
                
            # Equations to constraint number of virtual qubits on chain
            prob += Lpvars[var_S] + Lpvars[var_T] == chain_lengths[chain], "ZeroSum" + str(var_S) + str(var_T) 
        
    # Balancing equations for cells in multiple chains
    for cell in shared:
        if shared[cell]:
            prob += fixed[cell] + sum(  Lpvars[x] for x in list(shared[cell]) ) <= Z, "Models" + str(cell)
            
    
    
    return prob, var_map


def assignQubits(LpSolution, var_map):
    '''
    Format the solution given by PuLP into
    dictionaries key : path (source = s_sol, target = t_sol)
    :param LpSolution:
    '''

    if VERBOSE:
        print(LpSolution)

    for path in var_map:
        
        cell_S = path[0]
        cell_T = path[1]
        
        var_S = var_map[path][cell_S]
        limit_S = LpSolution[ var_S ]
        
        # join nodes to source until limit is reached
        # other nodes belong to target
        
        joined = 0        
        for qubit in _QCA.edges[cell_S,cell_T]['path']:
            qubit_shared = len(_RGraph.nodes[qubit]['cells']) > 1
            if (qubit_shared):
                if (joined < limit_S):
                    _QCA.nodes[cell_T]['qubits'].discard(qubit)
                    _RGraph.nodes[qubit]['cells'].discard(cell_T)
                    joined = joined + 1
                else:
                    _QCA.nodes[cell_S]['qubits'].discard(qubit)
                    _RGraph.nodes[qubit]['cells'].discard(cell_S)


def parseConfiguration(configuration):
    
    global M, N, L
     
    M, N, L = configuration['M'], configuration['N'], configuration['L']
    
    global WRITE, VERBOSE, PLOT
    
    PLOT    = configuration['PLOT']
    VERBOSE = configuration['VERBOSE']
    WRITE    = configuration['WRITELP']
    
    
    

def solveChains(RGraph, QCA, configuration):
    '''
    
    :param RGraph:
    :param QCA:
    :param paths:
    '''

    global _QCA, _RGraph

    # Parse solver parameters
    parseConfiguration(configuration)
    
    # Parse routing structures
    _RGraph = RGraph
    _QCA = QCA

    fixed, shared, chains, chain_lengths, sharing = parseSystem()
    
    if sharing:
        prob, var_map = setupProblem(fixed, shared, chains, chain_lengths)
        
        if WRITE:
            prob.writeLP("SHARING.lp")       
        
        prob.solve(solver=pulp.GLPK_CMD(msg=VERBOSE))
        
        # read solution
        LpSolution = {}
        for v in  prob.variables():
            LpSolution[v.name] = v.varValue
            
        
        assignQubits(LpSolution, var_map)
    
    