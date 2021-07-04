import numpy as np
import sympy as sp
import sympy.physics.quantum as spq

X = [{
    'match' : [spq.Ket('0')],
    'replace' : [spq.Ket('1')]
},{
    'match' : [spq.Ket('1')],
    'replace' : [spq.Ket('0')]
}]

H = [{
    'match' : [spq.Ket('0')],
    'replace' : [(spq.Ket('0')/sp.sqrt(2)) + (spq.Ket('1')/sp.sqrt(2))]
},{
    'match' : [spq.Ket('1')],
    'replace' : [(spq.Ket('0')/sp.sqrt(2)) - (spq.Ket('1')/sp.sqrt(2))]
}]

CX = [{
    'match' : [spq.Ket('1'),spq.Ket('0')],
    'replace' : [spq.Ket('1'),spq.Ket('1')]
},{
    'match' : [spq.Ket('1'),spq.Ket('1')],
    'replace' : [spq.Ket('1'),spq.Ket('0')]
}]