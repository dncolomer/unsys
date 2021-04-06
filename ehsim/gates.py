import numpy as np

# Gates:
# fmt: off
X = np.array([
    [0, 1],
    [1, 0]
], dtype=np.complex_)

Z = np.array([
    [1, 0],
    [0, -1]
], dtype=np.complex_)

I = np.array([
    [1, 0],
    [0, 1]
], dtype=np.complex_)

H = np.array(
    [
        [1 / np.sqrt(2), 1 / np.sqrt(2)],
        [1 / np.sqrt(2), -1 / np.sqrt(2)]
    ],
    dtype=np.complex_,
)

CX = np.array(
    [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]
    ], dtype=np.complex_
)

CH = np.array(
    [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1 / np.sqrt(2), 1 / np.sqrt(2)],
        [0, 0, 1 / np.sqrt(2), -1 / np.sqrt(2)],
    ],
    dtype=np.complex_,
)

CZ = np.array(
    [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, -1]
    ], dtype=np.complex_
)

II = np.array(
    [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ], dtype=np.complex_
)

CCZ = np.array(
    [
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, -1],
    ],
    dtype=np.complex_,
)

CCX = np.array(
    [
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 1, 0],
    ],
    dtype=np.complex_,
)

SWAP = np.array([
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1]
], dtype=np.complex_)
# fmt: on

# Pre: assumes all nodes are comp. basis elements of Z {0 or 1}
# Post: we leave the system in the Z basis
'''def X_rule(qubit):
    #Check all nodes of the qubit
    # 0 -> 1
    # 1 -> 0

# Pre: assumes all nodes are comp. basis elements of Z {0 or 1}
# Post: we leave the system in the Z basis
def H_rule(qubit):
    #Check all nodes of the qubit
    #Create branches
    # 0 -> (0),(1)
    # 1 -> (0),(-1)

def rewrite(rules,input_map=["q0","q1","q2"]):
    #process the rules

CCCX {[1,1,1,0] -> [1,1,1,1],[1,1,1,0] -> [1,1,1,1]}

CX [1,0] -> [1,1], [1,1] -> [1,0]

+- -> [0,0],[0,1],-[1,0],-[1,1]

[0,0],[0,1],-[1,1],-[1,0]

X = [{
	'match':[np.array([1, 0], dtype=np.complex_)],
    'replace':[np.array([0, 1], dtype=np.complex_)]
},
{
	'match':[np.array([0, 1], dtype=np.complex_)],
    'replace':[np.array([1, 0], dtype=np.complex_)]
}]

system = Hypergraph(1)
system.rewrite(X,["q0"])

Note to self, the actual rules should also do 
pattern matching so we can implement Multicontrl gates more generically

+ -> 0 / 1 -> 0 / 1 / 0 / -1 -> (after factoring?) -> 0
'''
