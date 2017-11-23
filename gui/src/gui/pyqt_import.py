#!/usr/bin/env python

'''
Handler for importing PyQt submodules independent of PyQt version. GUI code is
assumed to be written for Qt5 and certain requested submodules will be remapped
accordingly.
'''

__author__  = 'Jake Retallick'
__date__    = '2017-11-22'


from importlib import import_module

__pyqt_mods = ['PyQt5', 'PyQt4']

maps = {}
maps['PyQt4'] = {'QtWidgets' : 'QtGui'}
maps['PyQt5'] = {}

def _remap(pyqt, submod):
    mp = maps[pyqt.__name__]
    return mp[submod] if submod in mp else submod

def importPyQt(*submods):
    '''Returns a list of the requested submodule for PyQtX for your latest
    local version of PyQt. Input should either be a list of submodule names or
    a position based set of *args. Missing/Invalid submodules will return None
    at that site in the list.

    usage:
        QtGui, QtCore = importPyQt('QtGui', 'QtCore')
        QtGui, QtCore, QtSvg = importPyQt(['QtGui', 'QtCore', 'QtSvg'])
    '''

    if len(submods)==1 and isinstance(submods[0], list):
        submods = submods[0]

    # attempt to load a PyQtX module into pyqt
    for m in __pyqt_mods:
        try:
            pyqt = import_module(m)
        except ImportError:
            continue
        break
    else:
        print('One of [{0}] must be installed'.format(' '.join(__pyqt_mods)))
        return [None for _ in submods]

    # load all the submodules
    out = []
    for submod in submods:
        try:
            sm = import_module('{0}.{1}'.format(pyqt.__name__, _remap(pyqt, submod)))
        except ImportError:
            print('Failed to import submodule: {0}'.format(submod))
            sm = None
        out.append(sm)

    return out[0] if len(out)==1 else out
