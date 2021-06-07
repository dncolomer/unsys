import numpy as np
import sympy as sp
import sympy.physics.quantum as spq

#TODO get rid of sym_map
X = {
    'name':'X gate',
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

CX = {
    'name':'CX gate',
    'rules': [{
        'match':[spq.Ket('1'),spq.Ket('0')],
        'replace':[{
                'weight': 1,
                'kets': [spq.Ket('1'),spq.Ket('1')]
            }
        ]},{
        'match':[spq.Ket('1'),spq.Ket('1')],
        'replace':[{
                'weight': 1,
                'kets': [spq.Ket('1'),spq.Ket('0')]
            }
        ]}
    ]
}

H = {
    'name':'H gate',
    'rules': [{
        'match':[spq.Ket('0')],
        'replace':[{
                'weight': 1,
                'kets': [spq.Ket('0')/sp.sqrt(2) + spq.Ket('1')/sp.sqrt(2)]
            }
        ]},{
        'match':[spq.Ket('1')],
        'replace':[{
                'weight': 1,
                'kets': [spq.Ket('0')/sp.sqrt(2) - spq.Ket('1')/sp.sqrt(2)]
            }
        ]}
    ]
}

MS = {
    'name':'Mølmer–Sørensen gate',
    'rules': [{
        'match':[spq.Ket('0'), spq.Ket('0')],
        'replace':[{
                'weight': 1/sp.sqrt(2),
                'kets': [spq.Ket('0'), spq.Ket('0')]
            },{
                'weight': sp.I/sp.sqrt(2),
                'kets': [spq.Ket('1'), spq.Ket('1')]
            }
        ]},{
        'match':[spq.Ket('0'), spq.Ket('1')],
        'replace':[{
                'weight': 1/sp.sqrt(2),
                'kets': [spq.Ket('0'), spq.Ket('1')]
            },{
                'weight': -1*sp.I/sp.sqrt(2),
                'kets': [spq.Ket('1'), spq.Ket('0')]
            }
        ]},{
        'match':[spq.Ket('1'), spq.Ket('0')],
        'replace':[{
                'weight': 1/sp.sqrt(2),
                'kets': [spq.Ket('1'), spq.Ket('0')]
            },{
                'weight': -1*sp.I/sp.sqrt(2),
                'kets': [spq.Ket('0'), spq.Ket('1')]
            }
        ]},{
        'match':[spq.Ket('1'), spq.Ket('1')],
        'replace':[{
                'weight': 1/sp.sqrt(2),
                'kets': [spq.Ket('1'), spq.Ket('1')]
            },{
                'weight': sp.I/sp.sqrt(2),
                'kets': [spq.Ket('0'), spq.Ket('0')]
            }
        ]}
    ]
}