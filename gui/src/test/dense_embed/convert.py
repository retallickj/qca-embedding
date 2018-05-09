#------------------------------------------------------------------------------
# Name: convert.py
# Purpose:  Conversion between qbits and paths to vertex models: minimizes
#           maximum model size
#
# Author:   Jacob Retallick
# Created:  20.02.2015
# Last modified: 06.05.2015
#------------------------------------------------------------------------------

#from scipy.optimize import linprog     # need scipy.__version__ >= 0.15.1
import pulp     # python LP and MIP classes and binding
from random import shuffle
from copy import copy as cp
import os    # only for deleting LP solution file


### HANDLES ###

USE_DEFAULT = False     # flag for using default PuLP solver.
USE_LPR = False         # flag for using LP-relaxation. Current not implemented

### GLOBALS ###

_adj = {}       # adjacency dict for the chain structure
_chains = {}    # end-point keyed dict of lists of chain nodes
_end_occ = {}   # end_point keyed dict of lists of local chains

LP_FILE = './lp-sol.lp'     # temp file for lp problem file
DEL_LP_FILE = True          # delete temp file after use
WRITE_MAX = True                # write larger model size in header


####################################################################
### WRITE FUNCTION


def writeSol(models, max_model, fname):
    '''write models to given filename'''

    try:
        fp = open(fname, 'w')
    except:
        print 'Failed to open file for writing: %s' % fname
        return False

    # header
    if WRITE_MAX:
        fp.write('max model: %.2f\n' % max_model)

    for cell in models:
        qbits = models[cell]['qbits']
        fp.write('%d\t:\t%s\n' % (cell, str(qbits)))

    fp.close()
    return True


#####################################################################
### PROBLEM FORMULATION/FORMATTING


def preProc(paths, qbits):
    '''Generates the adjacency structure of the chains'''

    global _adj

    # inverse map
    qbit_map = {qbits[key]: key for key in qbits}

    # initialize adjacency dictionary
    keys = set(reduce(lambda x, y: list(x)+list(y), paths))
    _adj = {key: {} for key in keys}

    # path keys may or may not be sorted, use ends of paths for keyes
    for path in paths.values():
        key = tuple(map(lambda x: qbit_map[x], [path[0], path[-1]]))
        # store directional path
        _adj[key[0]][key[1]] = path[1:-1]
        _adj[key[1]][key[0]] = path[-2:0:-1]

    return paths.keys()


def mergeChains(ch1, ch2, start=True):
    '''Combine two chains and return the junction qbit'''

    global _chains, _end_occ, _adj

    # pop old chains and generate new chain
    if start:
        _end_occ.pop(ch1[0])
        _end_occ[ch1[1]].remove(ch1)
        if ch2[0] == ch1[0]:    # chains meet at starts
            _end_occ[ch2[1]].remove(ch2)
            new_nodes = _chains.pop(ch2)[::-1] + [ch1[0]] + _chains.pop(ch1)
            new_key = (ch2[1], ch1[1])
        else:   # end of ch2 meets start of ch 1
            _end_occ[ch2[0]].remove(ch2)
            new_nodes = _chains.pop(ch2) + [ch1[0]] + _chains.pop(ch1)
            new_key = (ch2[0], ch1[1])
    else:
        _end_occ.pop(ch1[1])
        _end_occ[ch1[0]].remove(ch1)
        if ch2[0] == ch1[1]:  # end of ch1 meets start of ch2
            _end_occ[ch2[1]].remove(ch2)
            new_nodes = _chains.pop(ch1) + [ch1[1]] + _chains.pop(ch2)
            new_key = (ch1[0], ch2[1])
        else:   # chains meet at ends
            _end_occ[ch2[0]].remove(ch2)
            new_nodes = _chains.pop(ch1) + [ch1[1]] + _chains.pop(ch2)[::-1]
            new_key = (ch1[0], ch2[0])

    _chains[new_key] = new_nodes

    # add new chain to _end_occ[ends]
    for key in new_key:
        _end_occ[key].append(new_key)

    return new_key


def addNode(ch, start=True):
    ''' Add the adjacent chain to the current chain '''

    global _chains, _end_occ, _adj

    if start:
        # generate new chain and key
        old_end = _chains[ch][0] if _chains[ch] else ch[1]
        new_end = filter(lambda x: x != old_end, _adj[ch[0]])[0]
        new_nodes = [ch[0]] + _chains.pop(ch)
        new_key = (new_end, ch[1])
        # update _end_occ
        _end_occ[ch[1]].remove(ch)
        _end_occ.pop(ch[0])
        _end_occ[ch[1]].append(new_key)
    else:
        # generate new chain and key
        old_end = _chains[ch][-1] if _chains[ch] else ch[0]
        new_end = filter(lambda x: x != old_end, _adj[ch[1]])[0]
        new_nodes = _chains.pop(ch) + [ch[1]]
        new_key = (ch[0], new_end)
        # update _end_occ
        _end_occ[ch[0]].remove(ch)
        _end_occ.pop(ch[1])
        _end_occ[ch[0]].append(new_key)

    if not new_end in _end_occ:
        _end_occ[new_end] = []
    _end_occ[new_end].append(new_key)
    _chains[new_key] = new_nodes

    return new_key


def extendChain(ch, start=True):
    '''Extend either the start or end of a chain by one node. Merge two
    connecting chains if connected by a node with 2 connections (wire)'''

    global _chains, _end_occ, _adj

    i = 0 if start else 1

    in_wire = lambda x: len(_adj[x]) == 2 

    # Ears are K_3 subgraphs with only one vertex connected to another vertex not in K_3
    in_ear = lambda x: len(set(_adj[x[0]]) & set(_adj[x[1]])) > 0 


    while in_wire(ch[i]) and not in_ear(ch):

        # assert _end_occ consistency
        n = len(_end_occ[ch[i]])
        try:
            assert n in [1, 2], 'Lost track of _end_occ: %d' % n
            assert ch in _end_occ[ch[i]], 'End point lost'
        except AssertionError:
            return ch

        # check for merge condition
        if n == 2:
            ch2 = filter(lambda x: x != ch, _end_occ[ch[i]])[0]
            ch = mergeChains(ch, ch2, start=start)
        else:   # add next node from start
            ch = addNode(ch, start=start)

    return ch


def extendChains(chains):
    '''Extend all chains until they terminate at high adjacency nodes '''

    global _adj, _chains, _end_occ

    _chains = {chain: [] for chain in chains}

    term_chains = {}
    end_points = set(reduce(lambda x, y: list(x)+list(y), _chains))

    # create list of chains ending at each point
    _end_occ = {key: [] for key in end_points}

    for key in chains:
        _end_occ[key[0]].append(key)
        _end_occ[key[1]].append(key)

    while _chains:

        # select some key of _chains
        ch = _chains.keys()[0]

        # extend the start of the chain to termination
        ch = extendChain(ch, start=True)

        # extend the end of the chain to termination
        ch = extendChain(ch, start=False)

        # store and pop finished chain
        if not ch in term_chains:
            term_chains[ch] = []
        term_chains[ch].append(_chains.pop(ch))
        for end in ch:
            _end_occ[end].remove(ch)
            if len(_end_occ[end]) == 0:
                _end_occ.pop(end)

    return term_chains


#######################################################################
### SOLVERS


def formatProblem(extended_chains, qbits):
    '''Construct a problem dict object with useful values'''

    global _adj

    # switch to lists for consistent indexing
    keys = []
    chains = []
    for key in extended_chains:
        for ch in extended_chains[key]:
            keys.append(key)
            chains.append(ch)

    # make list of chain node lengths
    node_lens = map(len, chains)

    # record indices for chains containing nodes
    l = len(node_lens)
    ind_nodes = [i for i in xrange(l) if node_lens[i] > 0]
    ind_none = [i for i in xrange(l) if node_lens[i] == 0]

    # generate lists of terminating nodes
    term_nodes = list(set(reduce(lambda x, y: list(x) + list(y), keys)))
    end_lists = {node: {'n': [], 'm': []} for node in term_nodes}

    for i in xrange(l):
        key = keys[i]
        end_lists[key[0]]['n'].append(i)
        end_lists[key[1]]['m'].append(i)

    # generate list of ordered internal qubits
    chain_qbits = []
    for i in xrange(l):
        key = keys[i]
        chain = []
        # add end points to nodes list
        nodes = cp(chains[i])
        nodes.insert(0, key[0])
        nodes.append(key[1])
        for i in xrange(len(nodes)-1):
            c1, c2 = nodes[i:i+2]
            # add connecting chain and end qbit
            chain += _adj[c1][c2] + [qbits[c2]]
        # forget end qbit
        chain.pop()
        chain_qbits.append(chain)

    # store chain lengths
    chain_lens = map(len, chain_qbits)

    # set up constraint parameters
    prob_dict = {'node_lens': node_lens,        # number of nodes in each chain
                 'chain_lens': chain_lens,      # number of qbits in each chain
                 'chain_qbits': chain_qbits,    # list of qbits in each chain
                 'keys': keys,                  # list of chain keys
                 'ind_nodes': ind_nodes,        # indices for chains with nodes
                 'ind_none': ind_none,          # indices for chains w/o nodes
                 'end_lists': end_lists,        # node keyed dict of n,m inds
                 'chains': chains,              # list of chains
                 'qbits': qbits}                # cell keyed qbit map

    return prob_dict


def allocateChain(seq_n, seq_m, seq_left, chain, chain_qb, models):
    ''' '''

    mn = len(seq_left)*1./len(chain)
    a = int(mn)
    alpha = int(round((a+1-mn)*len(chain)))
    sizes = [a]*alpha + [a+1]*(len(chain)-alpha)
    shuffle(sizes)  # randomly order group sizes

    k = 0
    for i in xrange(len(chain)):
        cell = chain[i]
        qbs = seq_left[k:k+sizes[i]]
        k += sizes[i]
        models[cell]['qbits'] += qbs


def solToModels(sol, prob_dict):
    '''Convert optimized solution parameters to vertex-models'''

    global _adj

    # assign parameters to shorter variable names
    qbits = prob_dict['qbits']  # cell -> qbit map
    keys = prob_dict['keys']    # chain keys
    chains = prob_dict['chains']    # list of nodes in each chain
    chain_qbits = prob_dict['chain_qbits']  # list of qbits in each chain
    ind_nodes = prob_dict['ind_nodes']  # indices of chains with nodes
    l = len(keys)   # number of chains

    if False:   # for debugging
        print '*'*30
        print 'Chains:\n'
        for i in xrange(len(chains)):
            ch = chains[i]
            n, m = int(round(sol['n'][i])), int(round(sol['m'][i]))
            print 'n: %d \tm: %d' % (n, m)
            print 'nodes: %s' % str(ch)
            print 'ends: %s' % str(keys[i])
            print 'qubits: %s\n' % str(chain_qbits[i])

    # initialize model dicts
    models = {}
    for cell in qbits:
        models[cell] = {'qbits': [], 'int_coup': [], 'ext_coup': []}

    # get set of end points
    end_points = set()
    for k in keys:
        end_points.add(k[0])
        end_points.add(k[1])

    # first handle long chains
    for i in xrange(l):
        key = keys[i]   # end points
        chain = chains[i]
        chain_qb = chain_qbits[i]

        # solution parameters
        n = int(round(sol['n'][i]))
        m = int(round(sol['m'][i]))

        # ordered chain sequences
        seq_n = chain_qb[:n]    # qbits in start model
        seq_m = chain_qb[-m:] if m > 0 else []  # qbits in end model
        seq_left = chain_qb[n:-m] if m > 0 else chain_qb[n:]    # left over

        # add seq_n couplers to key[0]
        models[key[0]]['qbits'] += seq_n
        seq_n.insert(0, qbits[key[0]])
        if len(seq_n) > 1:
            for j in xrange(len(seq_n)-1):
                models[key[0]]['int_coup'].append((seq_n[j], seq_n[j+1]))

        # add seq_m couplers to key[1]
        models[key[1]]['qbits'] += seq_m
        seq_m.append(qbits[key[1]])
        if len(seq_m) > 1:
            for j in xrange(len(seq_m)-1):
                models[key[1]]['int_coup'].append((seq_m[j], seq_m[j+1]))

        # key[0] external coupler
        q1 = seq_n[0]
        q2 = seq_left[0] if seq_left else seq_m[0]
        models[key[0]]['ext_coup'].append((q1, q2))

        # key[1] external coupler
        q2 = seq_m[-1]
        q1 = seq_left[-1] if seq_left else seq_n[-1]
        models[key[0]]['ext_coup'].append((q1, q2))

        # allocate remaining qubits to chain nodes
        if i in ind_nodes:
            allocateChain(seq_n, seq_m, seq_left, chain, chain_qb, models)

    # assign default models to the remaining cells
    for cell in models:
        if not models[cell]['qbits'] or cell in end_points:
            models[cell]['qbits'].append(qbits[cell])

    return models


def solveLP(prob_dict, verbose):
    '''Solve the optimation problem using either Mixed Integer Programming
    or LP-Relaxation'''

    K = len(prob_dict['keys'])
    N = prob_dict['node_lens']
    M = prob_dict['chain_lens']
    end_lists = prob_dict['end_lists']

    # initialize LP solver

    prob = pulp.LpProblem('LP-prob', sense=pulp.LpMinimize)

    ### set up problem

    ## parameter dicts for LpVariable generation
    if USE_LPR:
        par_dict = {'lowBound': 0, 'upBound': None, 'cat': pulp.LpContinuous}
    else:
        par_dict = {'lowBound': 0, 'upBound': None, 'cat': pulp.LpInteger}

    mu_dict = {'lowBound': 0, 'upBound': None, 'cat': pulp.LpContinuous}

    # variable to be optimized
    n = [pulp.LpVariable('n%d' % i, **par_dict) for i in xrange(K)]  # starts
    m = [pulp.LpVariable('m%d' % i, **par_dict) for i in xrange(K)]  # end
    mu = pulp.LpVariable('mu', **mu_dict)   # max model size

    ## objective function

    prob += mu

    ## constraints

    # end group sizes
    for node in end_lists:
        ns = map(lambda x: n[x], end_lists[node]['n'])
        ms = map(lambda x: m[x], end_lists[node]['m'])
        prob += pulp.lpSum(ns) + pulp.lpSum(ms) + 1 <= mu

    # average chain sizes
    for i in prob_dict['ind_nodes']:
        prob += M[i] - n[i] - m[i] <= mu*N[i]
        prob += n[i] + m[i] <= M[i] - N[i]

    # chain without internal nodes must have n+m==M
    for i in prob_dict['ind_none']:
        prob += n[i] + m[i] == M[i]

    # write lp file and solve
    try:
        prob.writeLP(LP_FILE)
        fn = LP_FILE
    except:
        print 'Invalid LP filename'
        prob.writeLP('./lp-sol.lp')
        fn = './lp-sol.lp'
    if verbose:
        print 'Saved to file: %s' % fn

    print '\n\n'
    if USE_DEFAULT:
        status = prob.solve()
    else:
        status = prob.solve(solver=pulp.GLPK_CMD())
    print '\n\n'

    # check solution status
    if verbose:
        print(pulp.LpStatus[status])

    # store solution values
    sol = {}
    sol['n'] = map(lambda v: pulp.value(v), n)
    sol['m'] = map(lambda v: pulp.value(v), m)
    sol['mu'] = pulp.value(mu)  # for consistency check

    # delete solution file
    if DEL_LP_FILE:
        if verbose:
            print 'Deleting file: %s' % fn
        os.remove(fn)

    if verbose:
        print sol['m']
        print sol['n']

    models = solToModels(sol, prob_dict)

    if verbose:
        for cell in models:
            print 'c: %s \t :: %s' % (str(cell), str(models[cell]['qbits']))

    return models, sol['mu']


#######################################################################
### MAIN

def convertToModels(paths, qbits, verbose=False):
    '''After running dense placement, takes the dict of paths and the dict
    defining the associated cell for each qubit and finds a pseudo-optimal
    vertex-model representation of each cell

    inputs: paths (dict)    : (qb1,qb2) keyed dict of qbit chains between
                             assigned qbits.
            qbits (dict)    : qbit keyed dict of qbits for each cell.

    outputs: models (dict)  : cell indexed dict of model dictionaries.

    Model dictionary format:
        key := cell label
        elements:   'qbits'     := list of qbits in model
                    'int_coup'  := list of internal couplers in model
                    'ext_coup'  := cell keyed dict of external couplers to
                                   adjacent cell models.

    CURRENT IMPLEMENTATION ONLY CONSIDERS model['qbits'] AS THAT IS ALL THAT
    IS NEEDED TO MATCH OUTPUT FROM D-WAVE'S HEURISTIC ALGORITHM. SHOULD
    IMPLEMENT COUPLERS FOR COMPLETENESS.
    '''

    global _adj

    if verbose:
        print '\n\nStarting conversion...'
        print 'Largest path size: %d' % (max(map(len, paths.values()))-2)

    ## generate adjacency dictionary from paths
    all_paths = preProc(paths, qbits)

    ## get list of long chains
    long_chains = filter(lambda x: len(_adj[x[0]][x[1]]) > 0, all_paths)

    models = {}

    if long_chains:
        # extend long chains
        extended_chains = extendChains(long_chains)
        # generate problem dictionary
        prob_dict = formatProblem(extended_chains, qbits)
        # solve for optimal model parameters
        models, max_model = solveLP(prob_dict, verbose=verbose)
        if models is None:
            print('Error occurred in model optimization...')
            return None, -1
    else:
        if verbose:
            print 'No extra qubits used... assigning default models'
        # generate model parameters
        for cell in qbits:
            qbit = qbits[cell]
            model = {}
            model['qbits'] = [qbit]
            model['int_coup'] = []
            model['ext_coup'] = {}
            for c2 in _adj[cell]:
                model['ext_coup'][c2] = tuple(sorted([qbit, qbits[c2]]))
            models[cell] = model
        max_model = 1

    return models, max_model
