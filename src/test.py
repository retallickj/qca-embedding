#!/usr/bin/env python
# encoding: utf-8

'''
Test script for the Dense Embedder
'''

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2018-01-30'      # last update


from core.qca import QCACircuit
from core.chimera import Chimera
from embed.dense_embed import DenseEmbedder

import sys

def main(qca_file):

    circuit = QCACircuit(fname=qca_file, verbose=False)
    source = circuit.subgraph_pn(adj='full', pn=0)
    chimera = Chimera(M=4, N=4, L=4)

    embedder = DenseEmbedder(chimera=chimera, logfile=None)
    embedder.embed(source, verbose=True)


if __name__ == '__main__':

    try:
        qca_file = sys.argv[1]
    except:
        print('No QCA file given')
        sys.exit()

    main(qca_file)
