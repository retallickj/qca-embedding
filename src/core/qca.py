#!/usr/bin/env python
# encoding: utf-8

'''
Specialised QCA Graph subclass, methods for parsing QCADesigner files, and
related methods
'''

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2017-10-31'      # last update


from __future__ import print_function   # needed for verbose print method
import numpy as np
from collections import namedtuple, defaultdict
from itertools import combinations

from pprint import pprint, pformat
import re

from graph import Graph, Node


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
_QCAD_CFs = {'QCAD_CELL_NORMAL': CFs['normal'],
             'QCAD_CELL_INPUT':  CFs['input'],
             'QCAD_CELL_OUTPUT': CFs['output'],
             'QCAD_CELL_FIXED':  CFs['fixed']}

# tags for different object types
_tags = {'design':  'TYPE:DESIGN',
         'cell':     'TYPE:QCADCell',
         'cell_obj': 'TYPE:QCADDesignObject',
         'dot':      'TYPE:CELL_DOT'}

# keys for cell values
_vals = {'cx': 'cell_options.cxCell',
         'cy':  'cell_options.cyCell',
         'cf':  'cell_function'}

class ParserNode:
    '''Nested Node structure for QCADesigner's specific html-like format'''

    tag_charset = '[a-zA-Z0-9_:]'
    val_charset = '[a-zA-Z0-9_\.\-\+]'

    rgx_in  = re.compile('\[({0}+)\]'.format(tag_charset))      # start tag
    rgx_out = re.compile('\[#({0}+)\]'.format(tag_charset))     # end tag
    rgx_val = re.compile('({0}+)=({0}+)'.format(val_charset))   # key=val pair

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
            if self.rgxMatch(self.rgx_in, line, 1):
                self.children.append(ParserNode(fp, tag=self.rgxm))
            elif self.rgxMatch(self.rgx_out, line, 1):
                if self.rgxm != self.tag:
                    print('tag mismatch: {0} :: {1}'.format(self.tag, self.rgxm))
                break
            else:
                m = self.rgx_val.match(line)
                if m.lastindex==2:
                    self.d[m.group(1)] = m.group(2)


    def __getitem__(self, key):
        '''Get the keyed Node value'''
        return self.d[key] if key in self.d else None

    def rgxMatch(self, rgx, s, n=0):
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

    def getChild(self, tag):
        '''Attempt to get the first child of the node with the given tag'''

        for child in self.children:
            if child.tag==tag:
                return child
        else:
            return None

def parse(fname):
    '''Parse a QCADesigner file and return the ParserNode head. If file is
    invalid, return None.

    input:
        fname   : filename of QCADesigner file.
    '''

    head = None
    with open(fname, 'r') as fp:
        head = ParserNode(fp, 'root')
    return head


# QCACircuit structure

class QCACell(Node):
    '''Specialized Node structure for describing QCA cells

    parameters:

        Node parameters

        n       : cell label
        adj     : adjacency dict of neighbouring cells
        v       : value contained in node
        x       : cell x location
        y       : cell y location

        Cell parameters

        cf      : cell function
        pol     : cell polarization
        rot     : cell is rotated
        dots    : list of cell quantum dot objects
    '''

    CFs = _CFs  # make cell-types externally visible through class prefix
    QCADot = namedtuple('QCADot', ['x y q']) # quantum dot position and charge
    Q0 = 1.60218e-19    # elementary charge for qdot charge normalization

    def __init__(self, n, adj=dict(), v=0., x=None, y=None,
                    cf=_CFs['normal'], pol=0, rot=False, dots=[]):
        '''Initialise a QCACell instance'''
        super(QCACell, self).__init__(n, adj, v, x, y)
        self.cf, self.pol, self.rot, self.dots = cf, pol, rot, dots

    @staticmethod
    def toDots(x, y, rot, pol, scale=1.):
        '''Compute suitable quantum dots for a QCA cell with a given scale.

        input:
            x       : x position of cell center
            y       : y position of cell center
            rot     : cell is rotated
            pol     : cell polarization -> dot charges
            scale   : side length of QCA cell
        '''

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

    @classmethod
    def fromParserNode(cls, n, node):
        '''Construct a new QCACell with the given label from a ParserNode'''

        # cell position from the cell's design object
        cell_obj = node.getChild(_tags['cell_obj'])
        if cell_obj is None:
            return None
        x,y = [float(cell_obj[k]) for k in 'xy']

        # cell parameters
        cf = _QCAD_CFs[node[_vals['cf']]]

        # quantum dot locations
        dots = node.extractAllNested(_tags['dot'])
        assert len(dots)==4, 'Invalid dot layout'
        dots = [cls.QCADot(float(d['x']), float(d['y']), float(d['charge'])/cls.Q0) for d in dots]

        pol = round((dots[0].q+dots[2].q-dots[1].q-dots[3].q)/2, 5)
        rot = len(set(round(d.x, 2) for d in dots))==3

        return QCACell(n, {}, 0, x, y, cf, pol, rot, dots)




class QCACircuit(Graph):
    '''General container for all relevant information about a QCA circuit.'''

    CFs = QCACell.CFs
    QCAD_CFs = {'QCAD_CELL_NORMAL': CFs['normal'],
                'QCAD_CELL_INPUT':  CFs['input'],
                'QCAD_CELL_OUTPUT': CFs['output'],
                'QCAD_CELL_FIXED':  CFs['fixed']}

    tags = {'design': 'TYPE:DESIGN',
            'cell': 'TYPE:QCADCell',
            'cell_obj': 'TYPE:QCADDesignObject',
            'dot': 'TYPE:CELL_DOT'}

    vals = {'cx': 'cell_options.cxCell',
            'cy': 'cell_options.cyCell',
            'cf': 'cell_function'}

    R0 = 2              # maximum range of interaction in normalized units


    def __init__(self, head=None, fname=None, verbose=False):
        '''Initialise a QCACircuit from an optional QCADesigner file parser'''

        super(QCACircuit, self).__init__()

        self.vprint = print if verbose else lambda *a, **k : None
        self.adj_masks =  {None:    self.__nullMask,
                          'full':   self.__fullMask,
                          'lim':    self.__limMask}

        self.cells = self.nodes     # dictionary of cells
        self.order = []

        self.clists = {'normal': [], 'input': [], 'output': [], 'fixed': []}

        if fname is not None:
            try:
                head = parse(fname)
            except:
                print('Failed to parse QCA data from file: {0}'.format(fname))
                head = None

        # extract cell data from node structure
        if isinstance(head, ParserNode):
            self.__fromHead(head)
            self.vprint('Circuit information extracted from circuit parser')
        else:
            self.vprint('No valid parser data found')


    def addCell(self, x, y, scale=True, cf='normal', pol=0, rot=False):
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

        dots = self.__toDots(x, y, rot, pol)

        n = len(self.cells)
        self.cells[n] = QCACell(n, {}, 0., x, y, cf, pol, rot, dots)

    def computeCoefs(self, adj=None, pn=None, gam=None):
        '''compute the normalized coefficients for the current QCACircuit under
        the given adjacency type and input cell polarization index'''

        assert adj in self.adj_masks, 'Invalid adjacency type'

        if len(self.cells)==0:
            if gam is None:
                return np.zeros(0), np.zeros([0,0])
            else:
                return np.zeros(0), np.zeros([0,0]), np.zeros(0)

        # get the kink energies for all cell interactions
        J = np.zeros([len(self.cells),]*2, dtype=float)
        for i,j in combinations(range(J.shape[0]),2):
            J[i,j]=J[j,i] = -self.computeEk(self.cells[i], self.cells[j])
        J /= np.max(np.abs(J))

        # enforce adjacency
        J *= self.adj_masks[adj]()

        # only want normal and output cell coefficients
        inds = [i for i,c in enumerate(self.cells)
                    if c.cf in [self.CFs['normal'], self.CFs['output']]]
        J = J[inds,:]

        # biases
        self.process()  # don't feel like being clever

        # input cell polarizations
        pols = np.zeros([len(self.cells),], dtype=float)
        for i in self.clists['fixed']:
            pols[i] = self.cells[i].pol
        if isinstance(pn, int):
            input_pols = self.intToPol(pn, len(self.clists['input']))
            for i, p in zip(self.clists['input'], input_pols):
                pols[i] = p

        h = np.dot(J, pols)
        J = J[:,inds]

        if gam is None:
            return h, J
        else:
            return h, J, np.ones(h.shape, dtype=float)*gam

    def process(self):
        '''Process useful parameters from the QCACircuit'''

        self.clists = {'normal': [], 'input': [], 'output': [], 'fixed': []}
        cfmap = {n:k for k,n in self.CFs.items()}
        for i, cell in self.cells.items():
            self.clists[cfmap[cell.cf]].append(i)

        self.vprint(pformat(self.clists))

    def reorderCells(self, metric=None):
        '''reorder the cells by the given metric. The metric needs to take a
        QCACell and return some sortable value'''

        # default metric
        if metric is None:
            metric = lambda cell: (cell.y, cell.x)

        self.cells.sort(key=metric)

    def computeEk(self, c1, c2):
        '''Compute the kink energy between two QCACells. Computation is only
        up to a normalization constant.'''

        # distance between cells
        if pow(c2.x-c1.x,2)+pow(c2.y-c1.y,2) > pow(self.R0*self.spacing,2):
            return 0

        # 4x4 matrix of distances between dots in either cell
        x1, y1 = zip(*[(d.x, d.y) for d in c1.dots])
        x2, y2 = zip(*[(d.x, d.y) for d in c2.dots])
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

    @staticmethod
    def intToPol(n, l):
        '''convert an integer into a polarization list of length L'''

        assert isinstance(l, int), 'Invalid binary size type'
        if l==0: return []
        b = format(int(n)%pow(2,l), '0{0}b'.format(l))
        return [2*int(p)-1 for p in b]

    # internal methods

    def __fromHead(self, head):
        '''Load the QCACircuit from the head Node of a parsed QCADesigner file'''

        design = head.getChild(self.tags['design'])
        assert design is not None, 'Invalid QCADesigner file format'

        cell_nodes = head.extractAllNested(self.tags['cell'])
        assert len(cell_nodes)>0, 'Empty QCADesigner file'

        # get cell-cell spacing
        cx, cy = [float(cell_nodes[0][self.vals[k]]) for k in ['cx', 'cy']]
        self.spacing = np.sqrt(cx*cy)

        # extract cell content
        for node in cell_nodes:
            cell = self.__nodeToCell(node)
            if cell:
                self.cells.append(cell)

    def __nodeToCell(self, node):
        '''Extract the relevant parameters from a ParserNode describing a
        QCACell'''

        # cell position from the cell's design object
        cell_obj = node.getChild(self.tags['cell_obj'])
        if cell_obj is None:
            return None
        x, y = [float(cell_obj[k]) for k in ['x', 'y']]

        # cell parameters
        cf = self.QCAD_CFs[node[self.vals['cf']]]

        # quantum dot locations
        dots = node.extractAllNested(self.tags['dot'])
        assert len(dots)==4, 'Invalid dot layout'

        dots = [self.QCADot(float(d['x']), float(d['y']), float(d['charge'])/self.Q0) for d in dots]

        pol = round((dots[0].q+dots[2].q-dots[1].q-dots[3].q)/2, 5)
        rot = len(set(round(d.x, 2) for d in dots))==3

        return QCACell(x, y, cf, pol, rot, dots)

    def __toDots(self, x, y, rot, pol):
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


    # adjacency masks

    @staticmethod
    def identifyInters(A, inv=True):
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

    def __cellInteractions(self):
        '''Characterise all the cell-cell interactions. Used by the masks:

        0   : no interaction, rot-nrot
        1   : (0,1) or (1,0)    nearest
        -1  : (1,1)             diagonal
        2   : (0,2) or (2,0)    next nearest
        -2  : (2,1) or (1,2)    next diagonal
        3   : adapter interaction, half offset

        '''

        X = np.array([c.x for c in self.cells]).reshape([-1,1])
        Y = np.array([c.y for c in self.cells]).reshape([-1,1])

        DX = np.abs(X.T-X)/self.spacing
        DY = np.abs(Y.T-Y)/self.spacing

        # cell interaction matrix
        A = np.zeros(DX.shape, dtype=int)
        zeq = lambda x,y : abs(x-y)<0.2     # approximate equality

        for i,j in combinations(range(A.shape[0]),2):
            dx, dy = DX[i,j], DY[i,j]
            if zeq(min(dx,dy), .5) and zeq(max(dx,dy), 1): A[i,j] = 3
            elif self.cells[i].rot != self.cells[j].rot: continue
            elif zeq(dx+dy, 1): A[i,j] = 1
            elif zeq(dx, 1) and zeq(dy, 1): A[i,j] = -1
            elif zeq(min(dx,dy),0) and zeq(max(dx,dy),2): A[i,j] = 2
            elif zeq(min(dx,dy),1) and zeq(max(dx,dy),2): A[i,j] = -2
        return A + A.T

    def __maskFill(self, crit=lambda i,j:True):
        '''Generate the mask using a defined criterion, crit return True/False'''

        N = len(self.cells)
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

        A = self.__cellInteractions()
        xovers = self.identifyInters(A, inv=False)
        crit = lambda i,j: j in xovers[i] or A[i,j] in [1, -1, 3]
        return self.__maskFill(crit=crit)

    def __limMask(self):
        '''Keep up to necssary inverting and xover interactions'''

        A = self.__cellInteractions()
        xovers = self.identifyInters(A, inv=False)
        invs = self.identifyInters(A, inv=True)
        crit = lambda i,j: j in xovers[i]+invs[i] or A[i,j] in [1, 3]
        return self.__maskFill(crit=crit)


    # echoing

    def echoCells(self, dots=False):
        for cell in self.cells:
            self.echoCell(cell, dots=dots)

    @classmethod
    def echoCell(cls, cell, dots=False):
        assert isinstance(cell, cls.QCACell)
        cfmap = {k:i for i,k in cls.CFs.items()}
        print('Cell:  x:{0:>6.2f}  ::  y:{1:>6.2f} \t {2:<10} pol:{3:>6.2f} \t rot:{4}'.format(cell.x, cell.y, cfmap[cell.cf], cell.pol, cell.rot))
        if dots:
            print('Dots: ')
            for d in cell.dots:
                print('       x:{0:>6.2f}  ::  y:{1:>6.2f} \t charge:{2:>4.1f}'.format(d.x, d.y, d.q))
            print('\n')
