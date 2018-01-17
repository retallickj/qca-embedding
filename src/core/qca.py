#!/usr/bin/env python
# encoding: utf-8

'''
Specialised QCA Graph subclass, methods for parsing QCADesigner files, and
related methods
'''

from __future__ import print_function   # needed for verbose print

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2017-11-06'      # last update

from collections import namedtuple, defaultdict
from itertools import combinations
import bisect
import numpy as np
import networkx as nx

from pprint import pformat
import re

from graph import Graph


def dget(d, key, default=None, mp=lambda x:x):
    '''Useful defaulted dict-like accessor method'''
    try:
        return mp(d[key])
    except:
        return default


# QCADesigner file parsing

# the following constants need to be used by local classes but should not
# be directly visible externally

# cell-type flags in QCACell
_CFs = {'normal':0, 'input':1, 'output':2, 'fixed':3}

# cell-type tags tags in QCADesigner files
_QCAD_CFs = {'QCAD_CELL_NORMAL': _CFs['normal'],
             'QCAD_CELL_INPUT':  _CFs['input'],
             'QCAD_CELL_OUTPUT': _CFs['output'],
             'QCAD_CELL_FIXED':  _CFs['fixed']}

# tags for different object types
_TAGS = {'design':  'TYPE:DESIGN',
         'cell':     'TYPE:QCADCell',
         'cell_obj': 'TYPE:QCADDesignObject',
         'dot':      'TYPE:CELL_DOT'}

# keys for cell values
_VALS = {'cx': 'cell_options.cxCell',
         'cy':  'cell_options.cyCell',
         'cf':  'cell_function'}

class ParserNode:
    '''Nested Node structure for QCADesigner's specific html-like format'''

    _tag_charset = '[a-zA-Z0-9_:]'
    _val_charset = '[a-zA-Z0-9_\.\-\+]'

    _rgx_in  = re.compile('\[({0}+)\]'.format(_tag_charset))      # start tag
    _rgx_out = re.compile('\[#({0}+)\]'.format(_tag_charset))     # end tag
    _rgx_val = re.compile('({0}+)=({0}+)'.format(_val_charset))   # key=val pair

    def __init__(self, fp, tag=None):
        '''Build the nested node hierarchy from the current file pointer

        inputs:
            fp  : current file pointer
            tag : tag of current Node
        '''

        self.tag = tag
        self.children = []
        self.d = {}
        self.rgxm = None

        for line in fp:
            if self._rgx_match(self._rgx_in, line, 1):
                self.children.append(ParserNode(fp, tag=self.rgxm))
            elif self._rgx_match(self._rgx_out, line, 1):
                if self.rgxm != self.tag:
                    print('tag mismatch: {0} :: {1}'.format(self.tag, self.rgxm))
                break
            else:
                m = self._rgx_val.match(line)
                if m.lastindex==2:
                    self.d[m.group(1)] = m.group(2)


    def __getitem__(self, key):
        '''Get the keyed Node value'''
        return self.d[key] if key in self.d else None

    def _rgx_match(self, rgx, s, n=0):
        '''March s with an re pattern and store the n^th group'''
        m = rgx.match(s)
        self.rgxm = m.group(n) if (m and m.lastindex >= n) else None
        return self.rgxm is not None

    def echo(self, npad=0):
        '''Recursive echo'''

        prefix = '*'*npad+' '
        def vprint(s):
            print('{0} {1}'.format(prefix, s))

        vprint('::NODE::')
        vprint('tag: {0}'.format(self.tag))
        vprint('nchildren: {0}'.format(len(self.children)))
        vprint('::fields')
        for key, val in self.d.items():
            vprint('{0} : {1}'.format(key, val))

        for child in self.children:
            child.echo(npad+1)

    def get_child(self, tag):
        '''Attempt to get the first child of the node with the given tag'''

        for child in self.children:
            if child.tag==tag:
                return child
        else:
            return None

    def extract_all_nested(self, tag):
        '''Get a list all Nodes at any depth below the node with the given
        tag. Will only find the firt instance of the tag in each branch'''

        if self.tag == tag:
            return [self,]

        nodes = []
        for child in self.children:
            nodes += child.extract_all_nested(tag)
        return nodes

    @staticmethod
    def parse(fname):
        '''Parse a QCADesigner file and return the ParserNode head. If file is
        invalid, return None

        input:
            fname   : filename of QCADesigner file
        '''

        head = None
        try:
            with open(fname, 'r') as fp:
                head = ParserNode(fp, 'root')
        except Exception as e:
            print('Failed to parse QCADesigner file with error:\n\t{0}'.format(e.message))
            head=None
        return head

# QCACircuit structure

class QCACircuit(Graph):
    '''General container for all relevant information about a QCA circuit.'''

    QCADot = namedtuple('QCADot', ['x', 'y', 'q'])

    CFs = _CFs          # make cell functions available through QCACircuit prefix
    R0 = 2              # maximum range of interaction in normalized units
    Q0 = 1.60218e-19    # elementary charge

    def __init__(self, head=None, fname=None, verbose=False):
        '''Initialise a QCACircuit from an optional QCADesigner file parser'''

        super(QCACircuit, self).__init__()

        self.vprint = print if verbose else lambda *a, **k : None
        self.adj_masks =  {None:    self.__nullMask,
                          'full':   self.__fullMask,
                          'lim':    self.__limMask}

        self.cells = []         # ordered list through insertion
        self.__cellkeys = []    # sorting keys for self.cells

        # default cell value for sorting
        self.metric = lambda c: (self.nodes[c]['y'], self.nodes[c]['x'])

        self.clists = {'normal': [], 'input': [], 'output': [], 'fixed': []}

        if fname is not None:
            head = ParserNode.parse(fname)

        # extract cell data from node structure
        if isinstance(head, ParserNode):
            self._from_head(head)
            self.vprint('Circuit information extracted from circuit parser')
        else:
            self.spacing = 1.
            self.vprint('No valid parser data found')


    def add_cell(self, x, y, scale=True, cf='normal', pol=0, rot=False):
        '''Add a new cell to the QCACircuit at the given location. No check is
        made for overlap with an existing cell.

        inputs:
            x       : cell x position
            y       : cell y position

            optional arguments

            scale   : x and y should be scaled by the scaling factor
            cf      : cell function
            pol     : cell polarization, if non-zero, overwrites cf to 'fixed'
            rot     : cell is rotated
        '''

        try:
            x, y = float(x), float(y)
        except:
            print('Invalid cell coordinates, must be castable to float')
            return

        if scale:
            x, y = x*self.spacing, y*self.spacing

        assert cf in _CFs, 'Invalid cell type'
        cf = _CFs[cf]
        pol = float(pol)
        if pol != 0:
            cf = _CFs['fixed']
        rot = bool(rot)

        dots = self._to_dots(x, y, rot, pol)

        self._add_cell(x, y, cf, pol, rot, dots)

    def subgraph_pn(self, adj=None, pn=None):
        '''Construct a subgraph excluding the input and fixed polarization cells.
        Node weights will be the cell biases for the given input index pn.

        input:
            adj : adjacency type, in self.adj_masks
            pn  : input polarization index
        '''

        inds, h, J = self.compute_coefs(adj=adj, pn=pn)

        G = nx.from_numpy_matrix(np.round(J,3))
        for k,(i,_h) in enumerate(zip(inds, h)):
            c = self.nodes[self.cells[i]]
            G.add_node(k, ind=i, weight=round(_h,3), x=c['x'], y=c['y'])

        return G


    def compute_coefs(self, adj=None, pn=None, gam=None):
        '''compute the normalized coefficients for the current QCACircuit under
        the given adjacency type and input cell polarization index

        inputs:
            adj : adjacency type
            pn  : polarization index := binary polarization flag
            gam : transverse field
        '''

        assert adj in self.adj_masks, 'Invalid adjacency type'

        if len(self)==0:
            if gam is None:
                return np.zeros(0), np.zeros([0,0])
            else:
                return np.zeros(0), np.zeros([0,0]), np.zeros(0)

        # convert all edges to ordered adjacency array
        J = np.asarray(nx.to_numpy_matrix(self, nodelist=self.cells, dtype=float))
        J /= np.max(np.abs(J))

        # enforce adjacency
        J *= self.adj_masks[adj]()

        # only want normal and output cell coefficients
        inds = [i for i,c in enumerate(self.cells)
                    if self.nodes[c]['cf'] in [_CFs['normal'], _CFs['output']]]
        J = J[inds,:]

        # biases
        self.process()  # don't feel like being clever

        # input cell polarizations
        pols = np.zeros([len(self),], dtype=float)
        for i in self.clists['fixed']:
            pols[i] = self.nodes[self.cells[i]]['pol']
        if isinstance(pn, int):
            input_pols = self.int_to_pol(pn, len(self.clists['input']))
            for i, p in zip(self.clists['input'], input_pols):
                pols[i] = p

        h = np.dot(J, pols)
        J = J[:,inds]

        if gam is None:
            return inds, h, J
        else:
            return inds, h, J, np.ones(h.shape, dtype=float)*gam

    def process(self):
        '''Process useful parameters from the QCACircuit'''

        self.clists = {'normal': [], 'input': [], 'output': [], 'fixed': []}
        cfmap = {n:k for k,n in _CFs.items()}
        for i, cell in enumerate(self.cells):
            self.clists[cfmap[self.nodes[cell]['cf']]].append(i)

        self.vprint(pformat(self.clists))

    def reorder_cells(self, metric=None):
        '''reorder the cells by the given metric. The metric needs to take a
        QCACell and return some sortable value'''

        # default metric
        if metric is None:
            metric = self.metric
        elif isinstance(metric, 'function'):
            try:
                metric(self.cells[0])
                self.metric = metric
            except:
                print('Metric failed test case')
                return

        self.cells.sort(key=self.metric)
        self.__cellkeys = [self.metric(c) for c in self.cells]

    @staticmethod
    def int_to_pol(n, l):
        '''convert an integer into a polarization list of length L'''

        assert isinstance(l, int), 'Invalid binary size type'
        if l==0: return []
        b = format(int(n)%pow(2,l), '0{0}b'.format(l))
        return [2*int(p)-1 for p in b]

    # internal methods

    def _from_head(self, head):
        '''Load the QCACircuit from the head Node of a parsed QCADesigner file'''

        design = head.get_child(_TAGS['design'])
        assert design is not None, 'Invalid QCADesigner file format'

        cell_nodes = head.extract_all_nested(_TAGS['cell'])
        assert len(cell_nodes)>0, 'Empty QCADesigner file'

        # get cell-cell spacing
        cx, cy = [float(cell_nodes[0][_VALS[k]]) for k in ['cx', 'cy']]
        self.spacing = np.sqrt(cx*cy)

        # extract cell content
        for node in cell_nodes:
            self._add_from_node(node)

    def _add_from_node(self, node):
        '''Extract cell parameters from a ParserNode describing a QCACell and
        add the cell to the QCACircuit'''

        # cell position from the cell's design object
        cell_obj = node.get_child(_TAGS['cell_obj'])
        if cell_obj is None:
            return None
        x, y = [float(cell_obj[k]) for k in ['x', 'y']]

        # cell parameters
        cf = _QCAD_CFs[node[_VALS['cf']]]

        # quantum dot locations
        dots = node.extract_all_nested(_TAGS['dot'])
        assert len(dots)==4, 'Invalid dot layout'

        dots = [self.QCADot(float(d['x']), float(d['y']), float(d['charge'])/self.Q0) for d in dots]

        pol = round((dots[0].q+dots[2].q-dots[1].q-dots[3].q)/2, 5)
        rot = len(set(round(d.x, 2) for d in dots))==3

        self._add_cell(x, y, cf, pol, rot, dots)

    def _add_cell(self, x, y, cf, pol, rot, dots):
        '''handler for adding a new cell to the circuit'''

        n = len(self)
        self.add_node(n, weight=0, x=x, y=y, cf=cf, pol=pol, rot=rot, dots=dots)

        # add the interactions (edges)
        for m in self.cells:
            Ek = self._compute_Ek(m, n)
            if Ek:
                self.add_edge(n, m, weight=-Ek)

        self._insort(n)

    def _to_dots(self, x, y, rot, pol):
        '''Compute suitable dots for the given QCACell parameters'''

        # dot locations
        if rot:
            dd = .5*self.spacing/np.sqrt(2)
            X = [x+dx for dx in [-dd, 0, dd, 0]]
            Y = [y+dy for dy in [0, -dd, 0, dd]]
        else:
            dd = .25*self.spacing
            X = [x+dx for dx in [dd, dd, -dd, -dd]]
            Y = [y+dy for dy in [-dd, dd, dd, -dd]]
        Q = [.5+dq for dq in [pol, -pol, pol, -pol]]

        return [self.QCADot(x, y, q) for x,y,q in zip(X,Y,Q)]

    def _compute_Ek(self, k1, k2):
        '''Compute the kink energy between two QCACells. Computation is only
        up to a normalization constant.'''

        c1, c2 = [self.nodes[k] for k in [k1,k2]]

        # distance between cells
        if pow(c2['x']-c1['x'],2)+pow(c2['y']-c1['y'],2) > pow(self.R0*self.spacing,2):
            return 0

        # 4x4 matrix of distances between dots in either cell
        x1, y1 = zip(*[(d.x, d.y) for d in c1['dots']])
        x2, y2 = zip(*[(d.x, d.y) for d in c2['dots']])
        X1 = np.array([x1,y1]).T.reshape([4,1,2])/self.spacing
        X2 = np.array([x2,y2]).T.reshape([1,4,2])/self.spacing
        R = np.sqrt(np.sum(pow(X1-X2,2), axis=2))

        if np.min(R) == 0:
            print('qdot overlap detected')
            return 0

        Q = np.array([1, -1, 1, -1])
        Q = np.outer(Q,Q)
        Ek = -np.sum(Q/R)

        return Ek

    def _insort(self, k):
        '''Insert node k into the sorted list of cells'''

        assert k in self, 'Invalid cell name'
        cell = self.nodes[k]

        key = self.metric(k)                        # sorting key
        ind = bisect.bisect(self.__cellkeys, key)   # insert location

        self.cells.insert(ind, k)
        self.__cellkeys.insert(ind, key)


    # adjacency masks

    @staticmethod
    def identify_inters(A, inv=True):
        '''Identify either all inverter or xover interactions from the
        interaction array'''

        N = A.shape[0]
        a = -1 if inv else 2
        cands = [(i,j) for i,j in combinations(range(N),2) if A[i,j]==a]

        # inverter condition, A=-1 and no A:1,1 path
        out = defaultdict(list)
        for i,j in cands:
            for k in range(N):
                if A[i,k]==A[k,j]==1:
                    break
            else:
                out[i].append(j)
        return out

    def _cell_interactions(self):
        '''Characterise all the cell-cell interactions. Used by the masks:

        0   : no interaction, rot-nrot
        1   : (0,1) or (1,0)    nearest
        -1  : (1,1)             diagonal
        2   : (0,2) or (2,0)    next nearest
        -2  : (2,1) or (1,2)    next diagonal
        3   : adapter interaction, half offset

        '''

        X = np.array([self.nodes[c]['x'] for c in self.cells]).reshape([-1,1])
        Y = np.array([self.nodes[c]['y'] for c in self.cells]).reshape([-1,1])

        DX = np.abs(X.T-X)/self.spacing
        DY = np.abs(Y.T-Y)/self.spacing

        # cell interaction matrix
        A = np.zeros(DX.shape, dtype=int)
        zeq = lambda x,y : abs(x-y)<0.2     # approximate equality: |x-y| < 0.2

        for i,j in combinations(range(A.shape[0]),2):
            c1, c2 = [self.nodes[self.cells[k]] for k in [i,j]]
            dx, dy = DX[i,j], DY[i,j]
            if zeq(min(dx,dy), .5) and zeq(max(dx,dy), 1):  A[i,j] = 3
            elif c1['rot'] != c2['rot']:                    A[i,j] = 0
            elif zeq(dx+dy, 1):                             A[i,j] = 1
            elif zeq(dx, 1) and zeq(dy, 1):                 A[i,j] = -1
            elif zeq(min(dx,dy),0) and zeq(max(dx,dy),2):   A[i,j] = 2
            elif zeq(min(dx,dy),1) and zeq(max(dx,dy),2):   A[i,j] = -2
        return A + A.T

    def __maskFill(self, crit=lambda i,j:True):
        '''Generate the mask using a defined criterion. The mask is a 2D array
        the same size as J with True wherever crit evaluates to nonzero'''

        N = len(self)
        M = np.ones([N,N], dtype=bool)

        for i,j in combinations(range(N), 2):
            if crit(i,j):
                continue
            else:
                M[i,j] = M[j,i] = False
        return M

    def __nullMask(self):
        '''Keep all interactions'''
        return self.__maskFill()

    def __fullMask(self):
        '''Keep up to diagonal and xover interactions'''

        A = self._cell_interactions()
        xovers = self.identify_inters(A, inv=False)
        crit = lambda i,j: j in xovers[i] or A[i,j] in [1, -1, 3]
        return self.__maskFill(crit=crit)

    def __limMask(self):
        '''Keep up to necssary inverting and xover interactions'''

        A = self._cell_interactions()
        xovers = self.identify_inters(A, inv=False)
        invs = self.identify_inters(A, inv=True)
        crit = lambda i,j: j in xovers[i]+invs[i] or A[i,j] in [1, 3]
        return self.__maskFill(crit=crit)


    # echoing

    def echo_cells(self, dots=False):
        for cell in self.cells:
            self.echo_cell(self.nodes[cell], dots=dots)

    @classmethod
    def echo_cell(cls, cell, dots=False):
        cfmap = {k:i for i,k in cls.CFs.items()}
        print('Cell:  x:{0:>6.2f}  ::  y:{1:>6.2f} \t {2:<10} pol:{3:>6.2f} \t rot:{4}'.format(
                    cell['x'], cell['y'], cfmap[cell['cf']], cell['pol'], cell['rot']))
        if dots:
            print('Dots: ')
            for d in cell['dots']:
                print('       x:{0:>6.2f}  ::  y:{1:>6.2f} \t charge:{2:>4.1f}'.format(d.x, d.y, d.q))
            print('\n')


if __name__ == '__main__':

    import sys
    from pprint import pprint

    try:
        fn = sys.argv[1]
    except:
        print('missing QCADesigner file')
        sys.exit()

    circuit = QCACircuit(fname=fn, verbose=False)

    G = circuit.subgraph_pn(adj='full', pn=2)

    for n, d in G.nodes.items():
        print(n, d)

    for e, d in G.edges.items():
        print(e, d)
