import numpy as np
import sympy as sp
import sympy.physics.quantum as spq

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

'''MS = {
    'name':'Mølmer–Sørensen gate',
    'sym_map': [spq.Ket('e'), spq.Ket('g')],
    'rules': [{
        'match':[0, 0],
        'replace':[{
                'weight': 1/spq.sqrt(2),
                'kets': [0, 0]
            },{
                'weight': sp.I/spq.sqrt(2),
                'kets': [1, 1]
            }
        ]},{
        'match':[0, 1],
        'replace':[{
                'weight': 1/spq.sqrt(2),
                'kets': [0, 1]
            },{
                'weight': -1*sp.I/spq.sqrt(2),
                'kets': [1, 0]
            }
        ]},{
        'match':[1, 0],
        'replace':[{
                'weight': 1/spq.sqrt(2),
                'kets': [1, 0]
            },{
                'weight': -1*sp.I/spq.sqrt(2),
                'kets': [0, 1]
            }
        ]},{
        'match':[1, 1],
        'replace':[{
                'weight': 1/spq.sqrt(2),
                'kets': [1, 1]
            },{
                'weight': sp.I/spq.sqrt(2),
                'kets': [0, 0]
            }
        ]},
    ]
}'''