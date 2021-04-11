import numpy as np
import sympy as sp
import sympy.physics.quantum as spq

#symbolic ket
#symket = 

#actual kets
zero_ket = np.array([1, 0], dtype=np.complex_)
one_ket = np.array([0, 1], dtype=np.complex_)
plus_ket = np.array([1 / np.sqrt(2), 1 / np.sqrt(2)], dtype=np.complex_)
minus_ket = np.array([1 / np.sqrt(2), -1 / np.sqrt(2)], dtype=np.complex_)

X = {
    'name':'X',
    'rules': [{
        'match':[spq.Ket(0)],
        'replace':[spq.Ket(1)]
    },
    {
        'match':[spq.Ket(1)],
        'replace':[spq.Ket(0)]
    }]
}

H = {
    'name':'H',
    'rules': [{
        'match':[spq.Ket(0)],
        'replace':[spq.Ket(0)/sp.sqrt(2) + spq.Ket(1)/sp.sqrt(2)]
    },
    {
        'match':[spq.Ket(1)],
        'replace':[spq.Ket(0)/sp.sqrt(2) - spq.Ket(1)/sp.sqrt(2)]
    }]
}

SWAP = {
    'name':'SWAP', 
    'rules': [{
        'match':[spq.Ket(0), spq.Ket(1)],
        'replace':[spq.Ket(1), spq.Ket(0)]
    },
    {
        'match':[spq.Ket(1), spq.Ket(0)],
        'replace':[spq.Ket(0), spq.Ket(1)]
    }]
} 

CX = {
    'name':'CX',
    'rules': [{
        'match':[spq.Ket(1), spq.Ket(1)],
        'replace':[spq.Ket(1), spq.Ket(0)]
    },
    {
        'match':[spq.Ket(1), spq.Ket(0)],
        'replace':[spq.Ket(1), spq.Ket(1)]
    }]
}

CCX = {
    'name':'CCX',
    'rules': [{
        'match':[spq.Ket(1), spq.Ket(1), spq.Ket(1)],
        'replace':[spq.Ket(1), spq.Ket(0)]
    },
    {
        'match':[spq.Ket(1), spq.Ket(1), spq.Ket(0)],
        'replace':[spq.Ket(1), spq.Ket(1)]
    }]
}