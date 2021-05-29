import numpy as np
import sympy as sp
import sympy.physics.quantum as spq

#TODO get rid of sym_map
X = {
    'name':'X gate',
    'sym_map': [spq.Ket('0'), spq.Ket('1')],
    'params_map':[],
    'rules': [{
        'match':[spq.Ket('0')],
        'replace':[{
                'weight': 1,
                'kets': [spq.Ket('1')]
            }
        ]},{
        'match':[spq.Ket('1')],
        'replace':[{
                'weight': 1,
                'kets': [spq.Ket('1')]
            }
        ]}
    ]
}

MS = {
    'name':'Mølmer–Sørensen gate',
    'rules': [{
        'match':[spq.Ket('e'), spq.Ket('e')],
        'replace':[{
                'weight': 1/sp.sqrt(2),
                'kets': [spq.Ket('e'), spq.Ket('e')]
            },{
                'weight': sp.I/sp.sqrt(2),
                'kets': [spq.Ket('g'), spq.Ket('g')]
            }
        ]},{
        'match':[spq.Ket('e'), spq.Ket('g')],
        'replace':[{
                'weight': 1/sp.sqrt(2),
                'kets': [spq.Ket('e'), spq.Ket('g')]
            },{
                'weight': -1*sp.I/sp.sqrt(2),
                'kets': [spq.Ket('g'), spq.Ket('e')]
            }
        ]},{
        'match':[spq.Ket('g'), spq.Ket('e')],
        'replace':[{
                'weight': 1/sp.sqrt(2),
                'kets': [spq.Ket('g'), spq.Ket('e')]
            },{
                'weight': -1*sp.I/sp.sqrt(2),
                'kets': [spq.Ket('e'), spq.Ket('g')]
            }
        ]},{
        'match':[1, 1],
        'replace':[{
                'weight': 1/sp.sqrt(2),
                'kets': [spq.Ket('g'), spq.Ket('g')]
            },{
                'weight': sp.I/sp.sqrt(2),
                'kets': [spq.Ket('e'), spq.Ket('e')]
            }
        ]}
    ]
}