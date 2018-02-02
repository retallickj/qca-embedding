#!/usr/bin/env python
# encoding: utf-8

'''
Negotiated congestion rougin algorithm for the dense placement algorithm
'''


__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2018-01-31'      # last update

from collections import defaultdict
from heapq import heappush, heappop
import numpy as np


class NodePars:
    '''All relevant parameters for a node in the routing'''
    users = 0       # number of paths using this node
    marked = False  # node marked for use in a path
    active = True   # node is available for use in routing
    hist_cost = 0   # history cost of the node

class Router:
    '''Router for the Dense Placement embedder'''

    COST_BASE = 1.0             # base cost for adding a new qubit
    COST_HISTORY = 1.0          # increment in hist_cost for a qubit being shared
    COST_BREAK = 9999           # cost threshold for detecting a failed routing
    COST_DISABLED = int(1e8)    # arbitrary high cost for disabled qubit

    INC_SHARING = 1.0       # increment in sharing cost per iteration
    BREAK_SHARING = 20      # threshold for sharing cost to assert failed routing
    RATE_FORGET = 0.001     # forget rate for hist_cost: as exp(-RATE_FORGET)

    def __init__(self):
        ''' '''

        self._forget_fact = np.exp(-self.RATE_FORGET)

    def initialise(self, source):
        '''Set up the Router for a given Graph'''

        # assume at this point source has already been checked
        self.source = source

        # parameters
        self._sharing_cost = 1.0
        self._pars = {node: NodePars() for node in self.source}

        # path storage
        self._cpaths = []   # candidate paths for each iteration
        self.paths = {}     # paths for each route, (start,end) indexed

    def enable_nodes(self, nodes):
        '''Enable a given set of nodes if inactive'''

        for node in nodes:
            if node in self._pars and not self._pars[node].active:
                self._pars[node].active = True
                self._pars[node].hist_cost = 0

    def disable_nodes(self, nodes):
        '''Disable a given set of nodes if active'''

        for node in nodes:
            if node in self._pars and self._pars[node].active:
                self._pars[node].active = False
                self._pars[node].hist_cost = self.COST_DISABLED

    def route(self, routes, reserved):
        '''Run negotiated congestion for the set of given routes. Special
        consideration is given to reserved qubits so they must be pass as an
        input.'''

        # reset data
        self.reset_data()

        rt_set = set([node for rt in routes for node in rt])    # list of route nodes
        res_nodes = set()   # set of all reserved nodes
        for s in reserved.values():
            res_nodes.update(s)

        # enable route ends
        self.enable_nodes(rt_set)

        # main iteration looped
        iter_count = 0
        while iter_count==0 or max(self._pars[n].users for n in self._pars)>1:
            iter_count += 1

            # release flags and forget any working paths
            self._reset_users()
            self.paths = {}

            for rt in routes:
                # mark off all end-points and reserved nodes as used
                for node in rt_set|res_nodes:
                    self._pars[end].marked = True
                self.paths[rt] = self._best_path(rt, reserved)
                self._update_pars()

            # update history costs
            self._gen_hist()

            # update sharing cost, end if too high
            self._sharing_cost += self.INC_SHARING
            if self._sharing_cost > self.BREAK_SHARING:
                break

        # handle end conditions
        if self._sharing_cost > self.BREAK_SHARING:
            return self.COST_BREAK

        # compute total cost of all routes
        cost = sum(self._nodes_cost(path) for rt, path in self.paths.items())

        # disable end points
        self.disable_nodes(rt_set)

        return cost



    # internal methods

    def _reset_users(self):
        '''Reset the number of users for each node'''

        for key in self._pars:
            self._pars[key].users = 0

    def _update_pars(self):
        '''Update the user counts'''

        self._reset_users()
        ends = defaultdict(int)

        # user count, paths include the ends
        for rt, path in self.paths.items():
            ends[rt[0]] += 1
            ends[rt[1]] += 1
            for node in path:
                self._pars[node].users += 1

        # only count at most one end per each node
        for node, count in ends.items():
            if count > 1:
                self._pars[node] -= (count-1)

    def _reset_data(self):
        '''Reset trial specific data'''

        self._cpaths = []
        self._sharing_cost = 1.0
        for node in self._pars:
            self._pars[node].users = 0

    def _gen_hist(self):
        '''Update history cost'''

        for node, par in self._pars.items():
            if par.hist_cost > 0:
                par.hist_cost *= self._forget_fact
            if par.users > 1:
                par.hist_cost += self.COST_HISTORY

    def _sort_cpaths(self):
        '''Sort the candidate paths by cost'''
        pass

    def _node_cost(self, node, inc=False):
        '''Calculate the cost of the given node. If inc is True, calculates the
        cost if one extra path included the node.'''

        par = self._pars[node]
        cost = (par.users+inc)*(self.COST_BASE+par.hist_cost)
        if par.users>1:
            cost *= self._sharing_cost
        return cost

    def _nodes_cost(self, nodes):
        '''Calculate the cost of a set of nodes'''
        return sum(self._node_cost(node) for node in nodes)

    def _expand_path(self):
        '''Get the list of paths extended from the current cheapest.'''

        path = heappop(self._cpaths)            # pop lowest cost path
        nghb = self.source.neighbors(path[-1])  # neighbors of path end node

        # add possible extensions with updated path costs
        new_paths = []
        for node in filter(lambda nd: not self._pars[nd].marked, nghb):
            self._pars[node].marked = True
            tpath = path + [node]   # new path
            tpath[0] += self._node_cost(node, inc=True)
            new_paths.append(tpath)
        return new_paths

    def _best_path(self, route, reserved):
        '''Determine the best path for the given route. Reserved[node] should
        give the set of nodes which have been reserved for routes to node.'''

        # add start of path
        heappush(self._cpaths, [0, route[0])  # first path value is the cost

        # free reserved nodes for routing
        for node in reserved[route[0]].union(reserved[route[1]]):
            self._pars[node].marked = False
        self._pars[route[1]].marked = False     # free end node

        # search for cheapest path for route
        check = True
        while check and self._cpaths:
            for path in self._expand_path():
                if path[-1] == route[1]:
                    best = path[1::]    # best path without cost
                    check = False
                    break
                else:
                    heappush(self._cpaths, path)

        self._cpaths = []
        for node, par in self._pars:
            par.marked = False

        return best

    def _reset_nodes(self):
        '''Reset active status and hist_cost of all nodes'''

        for node, par in self._pars:
            par.active = True
            par.hist_cost = 0
