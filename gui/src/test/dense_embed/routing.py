#---------------------------------------------------------
# Name: routing.py
# Purpose: Code for routing qubit chains between assigned qubits.
# Author:	Jacob Retallick
# Created: 06.08.2014
# Last Modified: 06.05.2015
#---------------------------------------------------------

from math import exp
from copy import deepcopy as cp

### GLOBALS ###

# VARIABLES

_paths = []         # potential paths during each iteration
_allPaths = {}      # path for each route
_curr_used = {}     # flag for path inclusion of a qbit (by any path)

_is_shared = {}     # flag for qbit sharing between paths
_is_used = {}       # count of number of paths through qbit
_active = {}        # flag for active qbits

_qbitAdj = {}       # adjacency dict for qubits

_sharing_cost = 1.0    # cost scaler assigned for shared qbits
_hist_cost = {}        # persistent cost

# CONSTANTS

COST_BASE = 1.0             # base cost for adding a new qubit
COST_HISTORY = 1.0          # increment in hist_cost for a qubit being shared
COST_BREAK = 9999           # cost threshold for detecting a failed routing
COST_DISABLED = int(1e8)    # arbitrary high cost for disabled qubit

INC_SHARING = 1.0       # increment in sharing cost per iteration
BREAK_SHARING = 20      # threshold for sharing cost to assert failed routing
RATE_FORGET = 0.001     # forget rate for hist_cost: as exp(-RATE_FORGET)


def initialize(qbitAdj):
    '''Initialise routing solver. Only call once per embedding trial'''
    global _paths, _allPaths, _curr_used, _is_shared
    global _is_used, _active, _qbitAdj, _hist_cost, _sharing_cost

    _qbitAdj = cp(qbitAdj)

    _paths, _allPaths = [], {}
    _curr_used, _is_shared, _is_used, _active, _hist_cost = {}, {}, {}, {}, {}

    _sharing_cost = 1.0

    for key in qbitAdj:
        _curr_used[key] = False
        _is_shared[key] = False
        _is_used[key] = 0
        _active[key] = True
        _hist_cost[key] = 0


def resetFlags():
    '''Reset is_shared and is_used flags'''

    global _is_shared, _is_used, _allPaths

    for key in _is_shared:
        _is_shared[key] = False
        _is_used[key] = 0


def updateFlags():
    '''Update is_shared and is_used flags'''
    global _allPaths, _is_shared, _is_used

    resetFlags()
    end_lists = {key: 0 for key in _is_shared}

    # update used count
    for path in _allPaths.values():
        end_lists[path[0]] += 1
        end_lists[path[-1]] += 1
        for key in path:
            _is_used[key] += 1

    # account for path ends
    for key in end_lists:
        if end_lists[key] > 1:
            _is_used[key] -= (end_lists[key]-1)

    # flag shared qubits
    for key in _is_shared:
        if _is_used[key] > 1:
            _is_shared[key] = True


def resetData():
    '''Reset trial specific data'''
    global _paths, _is_used, _qbitAdj, _sharing_cost

    _is_used = {k: 0 for k in _qbitAdj}
    _paths = []
    _sharing_cost = 1.0


def genHist():
    '''Update hist_cost'''
    global _hist_cost, _is_shared

    for key in _hist_cost:
        if _hist_cost[key] > 0:
            _hist_cost[key] *= exp(-RATE_FORGET)
        if _is_shared[key]:
            _hist_cost[key] += COST_HISTORY


def sortPaths():
    '''Sort paths by cost '''
    global _paths

    if len(_paths) == 0:
        #print 'InternalRoutingError: no detected paths for indexOfPath'
        return -1

    # sort paths by cost
    _paths.sort(key=lambda x: x[0])


def nodeCost(qbit):
    '''Calculate cost of given node'''
    global _is_used, _is_shared, _hist_cost, _sharing_cost

    if _is_shared[qbit]:
        return _is_used[qbit]*(COST_BASE+_hist_cost[qbit])*_sharing_cost
    else:
        return _is_used[qbit]*(COST_BASE+_hist_cost[qbit])


def expandPath(unavailable = []):
    '''Expand lowest cost path'''
    global _paths, _is_used, _curr_used, _qbitAdj

    # select expanding path, remove path from list of paths
    path = _paths.pop(0)
    # possible extensions
    extensions = [qb for qb in _qbitAdj[path[-1]] if not _curr_used[qb]]

    new_paths = []
    for qbit in extensions:
        _curr_used[qbit] = True
        # new path
        temp_new = cp(path)
        temp_new.append(qbit)
        # new cost
        _is_used[qbit] += 1            # cost including new qbit
        temp_new[0] += nodeCost(qbit)
        _is_used[qbit] -= 1            # qbit only used if goal reached
        # add extended paths to list of paths
        new_paths.append(temp_new)

    _paths += new_paths
    return new_paths


def bestPath(route, reserved):
    '''Determine the best path for the given route'''
    global _paths, _curr_used

    if route[0] == route[1]:
        print 'bestPath ERROR: start is same as goal!'
        return 0

    # initialise path
    new_path, goal = [0, route[0]], route[1]    # first value of path is cost
    _paths.append(new_path)     # add to paths
    check = False

    # free reserved qbits
    for qb in reserved[route[0]].union(reserved[route[1]]):
        _curr_used[qb] = False

    _curr_used[route[1]] = False    # free last qbit

    while not check and _paths:
        sortPaths()        # sort paths by increasing cost
        new_paths = expandPath()    # expand cheapest path
        # check for goal
        for path in new_paths:
            if path[-1] == goal:
                check = True
                output = path[1::]
                break

    if not check:
        print 'bestPath ERROR: ran out of paths... unconnected graph?'
        return 0

    # reset values
    _paths = []
    for key in _curr_used:
        _curr_used[key] = False

    return output


# not implemented
def writeToFile(writePath):
    '''Write routing status to file'''
    global _allPaths

    print('writeToFile for routing is not implemented')
    return

    fp = open(writePath, 'a')    # append to file
    fp.write('new:\n')

    for path in _allPaths.values():
        for qbit in path:
            pass


def enableQubits(qbits):
    '''enable given qubits if inactive'''
    global _active, _hist_cost

    #print('\n***Enabling qubits: \n%s\n' % str(qbits))
    for qbit in qbits:
        # only enable inactive qubits in the current target graph
        if qbit in _active and _active[qbit] is False:
            _active[qbit] = True
            _hist_cost[qbit] = 0


def disableQubits(qbits):
    '''disable given qubits if active'''
    global _active, _hist_cost

    #print('\n***Disabling qubits: \n%s\n' % str(qbits))
    for qbit in qbits:
        # only disable active qubits in the current target graph
        if qbit in _active and _active[qbit] is True:
            _active[qbit] = False
            _hist_cost[qbit] = COST_DISABLED


def resetQubits():
    '''Reset active status and hist_cost of all qubits'''

    for qbit in _active:
        _active[qbit] = True
        _hist_cost[qbit] = 0


def getPaths():
    '''Return paths'''
    return _allPaths

def Routing(routes, reserved, writePath=''):
    '''Run routing algorithm. Find the lowest cost mutual paths for the list
    of routes to facilitate. Special consideration is given to reserved qubits
    so they must be passed as an input.'''
    global _flags, _paths, _allPaths, _sharing_cost, _curr_used

    ## Reset Data
    resetData()

    ## Negotiated Congestion and Routing

    _is_shared[_is_shared.keys()[0]] = True    # set loop condition

    rt_set = set([it for rt in routes for it in rt])    # list of route qbits
    res_qbits = set()   # set of reserved qubits.
    for s in reserved.values():
        res_qbits.update(s)

    # enable end qubits for routes
    enableQubits(rt_set)

    # iteration loop
    while any(_is_shared.values()):

        # release flags
        resetFlags()
        _allPaths = {}

        for i in xrange(len(routes)):
            # mark off all end-points as used
            for rt in routes:
                _curr_used[rt[0]] = True
                _curr_used[rt[1]] = True
            for qb in res_qbits:
                _curr_used[qb] = True
            rt = routes[i]
            _allPaths[i] = bestPath(rt, reserved)  # find best path for route
            updateFlags()                # update flags

        genHist()   # update hist_cost

        # write to file
        if writePath:
            writeToFile(writePath)

        # update sharing cost
        _sharing_cost += INC_SHARING
        if _sharing_cost > BREAK_SHARING:
            break

    ## Handle end conditions

    # No route found
    if _sharing_cost > BREAK_SHARING:
            #print 'Routing BREAK ERROR: the routing timed out'
            return COST_BREAK

    # Compute total paths cost
    cost = sum([nodeCost(x) for path in _allPaths.values() for x in path])

    # disable end points
    disableQubits(rt_set)

    return cost
