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


 #1 If I know I am in a particular top-level hyperedge then I know the rest of the elements in the hyper edge are true
 #2 


#If gate is not entangling it can be broken down to qubit level operations
#Example Hadamard Layer
0+1
0
1

HLayer

#local match-replace
0 -> 0+1
1 -> 0-1

#non-local match-replace
1 -> {0 -> 1, 1 -> 0}

0+1 / 0-1
#lcomb
(0 / 1) / 0-1

#fmatch on ctrl
#but target qubit lcomb
(0 / 1) / (0 / -1)

(0 / 1) / (1 / -0)

#lcomnb match splits nodes
#fmatch on a nonsplit node on a non-local repl on a fmatch hit creates a hyperedge?




0+1 -> 0 / 1 - 0+1 / 0-1 -> 0






# assume 1 qubit node per hyperedge at max
# assume a mapping from index to qubit i.e. [q0,q1]
CX = [{
    'match' : [spq.Ket('1'),spq.Ket('0')], #equal match or linearly combined match
    'replace' : [spq.Ket('1'),spq.Ket('1')]
},{
    'match' : [spq.Ket('1'),spq.Ket('1')], #equal match or linearly combined match
    'replace' : [spq.Ket('1'),spq.Ket('0')]
}]

CCX = [{
    'match' : [spq.Ket('1'),spq.Ket('1'),spq.Ket('0')], #equal match or linearly combined match
    'replace' : [spq.Ket('1'),spq.Ket('1'),spq.Ket('1')]
},{
    'match' : [spq.Ket('1'),spq.Ket('1'),spq.Ket('1')], #equal match or linearly combined match
    'replace' : [spq.Ket('1'),spq.Ket('1'),spq.Ket('0')]
}]

q0  q1
0   0 (Init)
+   0 (H q0)


