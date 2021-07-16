from itertools import combinations

import numpy as np
import sympy as sp
import sympy.physics.quantum as spq
import string

# global variables
guid = 0

def getUID(prefix="id#"):
    global guid

    guid = guid + 1
    return prefix + str(guid - 1)

# qudit is a unique label identifying the qudit system the state belongs to
class BasisState:
    def __init__(self, qudit, ket=spq.Ket('0'), amplitude=1):
        self.qudit = qudit
        self.ket = ket
        self.amplitude = amplitude

#pure hanges the interpretation of qubits with multiple states
class StateCorrelation:
    def __init__(self, weight):
        self.uid = getUID("c#")
        self.state_map = {}
    
    #s is a Basis State
    def addState(self, s):
        if (s.qudit not in self.state_map.keys()):
            self.state_map[s.qudit] = [] 
  
        self.state_map[s.qudit].append(s)

class QuditSystem:
    def __init__(self, nb_qudits, dim, symbolic= False):
        global guid

        self.dim = dim
        self.qudits = {}
        self.correlations = {}
        self.nb_qudits = nb_qudits

        #init qudits
        for i in range(0, nb_qudits):
            qudit_label = "q" + str(i)
            self.qudits[qudit_label] = []

        for qudit in self.qudits:
            if (symbolic):
                j = 0
                while j < dim:
                    amplitude = sp.symbols(qudit_label+'_'+str(j))
                    state = BasisState(qudit, spq.Ket(j), amplitude)

                    self.qudits[qudit].append(state)
                    
                    j += 1
            else:
                q = BasisState(qudit)
                self.qudits[qudit].append(q)
        
    def solveInterference():
        pass

    def solveEntanglement():
        pass