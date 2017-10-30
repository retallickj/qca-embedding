import json
import os

import numpy as np


class Solution:
    '''General class describing an Ising solution'''

    def __init__(self, fname=None):
        '''Initialise a Solution object'''

        # parameters
        self.nodes = []         # ordered labels for problem nodes
        self.J = None           # interaction parameters
        self.h = None           # biases
        self.hash = None        # hash key
        self.fname = None       # file path
        self.qca_fname = None   # file path to qca file

        self.states = []    # observed outcomes
        self.energies = []  # energies for the observed outcomes
        self.probs = []      # probability of each outcome

        if fname:
            self.from_json(fname)

    def from_json(self, fname):
        '''Load solution from json file'''

        try:
            fp = open(fname, 'r')
        except IOError:
            print('Failed to load file: {0}'.format(fname))
            raise IOError

        data = json.load(fp)
        fp.close()

        # parameters
        self.nodes = data['nodes']
        self.J = np.array(data['J'])
        self.h = np.array(data['h']).reshape([-1,])
        self.hash = data['hash']
        self.fname = fname
        if 'qca_name' in data:
            self.qca_fname = data['qca_fname']

        self.states = data['states']
        self.energies = data['energies']
        self.probs = data['probs']

    def to_file(self, fname=None):
        '''Save solution to json file'''

        data = {}
        data['nodes'] = self.nodes
        data['J'] = self.J
        data['h'] = self.h
        data['hash'] = self.hash
        data['qca_fname'] = self.qca_fname

        data['states'] = self.states
        data['energies'] = self.energies
        data['probs'] = self.probs

        if fname is None:
            fname = self.fname

        # build save location
        dir_name = os.path.dirname(fname)
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        # output to json file
        fp = open(fname, 'w')
        json.dump(data, fp)
        fp.close()

    def model_maj(self, models, ind_map=None):
        '''Reduce given models to logical qubits'''
        pass
