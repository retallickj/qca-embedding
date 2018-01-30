#!/usr/bin/env python
# encoding: utf-8

'''
Useful utility functions
'''

__author__      = 'Jake Retallick'
__copyright__   = 'MIT License'
__version__     = '2.0'
__date__        = '2017-11-28'  # last significant update


from itertools import product

def range_product(*counts):
    '''Compact shorthand for iterating over multiple integer ranges.

    note:
        range_product(i,j,k,...) is equivalent to
        itertools.product(range(i), range(j), range(k), ...)

    usage:
        for i,j,k in range_product(5,3,6):
            print(i,j,k)    # i looped in range(5), j in range(3), etc.

        for i,j in range_product(2, 7):
            print(i,j)      # i looped in range(2), j in range(7)

    '''

    for x in product(*[range(k) for k in counts]):
        yield x


def dget(dict_, key, default=None, mp=lambda x:x, chk=lambda x: True):
    '''Defaulted accessor for dictionary elements with an optional mapping.

    inputs:
        dict_   : dictionary to access
        key     : dictionary key
        default : default value if key is not found in the dictionary
        mp      : map to apply to dictionary value if found
        chk     : check statement that must be passed else default

    usage:

    d = {'a': '5', 'test': 'val'}
    dget(d, 'test', 2)      # returns 'val'
    dget(d, 'a', 0, int)    # returns int('5')
    dget(d, 'b', 'foo')     # returns 'foo'
    '''
    if key in dict_:
        x = mp(dict_[key])
        return x if chk(x) else default
    return default

def bounded(n, lo, hi):
    '''Check if lo <= n <= hi'''
    return n>=lo and n<=hi
