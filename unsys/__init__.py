from itertools import combinations

import numpy as np
import sympy as sp
import sympy.physics.quantum as spq
import string

# global variables
qudit_nb = 0

def getUID(prefix="id"):
    global qudit_nb

    qudit_nb = qudit_nb + 1
    return prefix + str(qudit_nb - 1)

#state_map is a map with the qubit labels as keys and their spq state expression as value
class State:
    def __init__(self, state_map):
        self.state_map = state_map
        self.uid = getUID("state#")

#pure hanges the interpretation of qubits with multiple states
class Qudit:
    def __init__(self, label, pure= True):
        self.uid = getUID("qudit#")
        self.label = label
        self.pure = pure

class QuditSystem:
    def __init__(self, nb_qudits, d, initial_states= None, symbolic= False):
        global qudit_nb

        self.dimension = d
        self.qudits = {}
        self.states = {}
        self.nb_qudits = nb_qudits

        if (initial_states is None):
            qudit_nb = 0
            for i in range(0, nb_qudits):
                qudit_label = "q" + str(i)
                q = Qudit(qudit_label)                        
                self.qudits[q.uid] = q

                if (symbolic):
                    sv_size = d
                    i = 0
                    while i < sv_size:
                        sym = sp.symbols(qudit_label+'_'+str(i))
                        s = sym * spq.Ket(i)
                        
                        state = State({q.uid: s})
                        self.states[state.uid] = state
                        
                        i += 1
                else:
                    s = spq.Ket('0')
                    state = State({q.uid: s})
                    self.states[state.uid] = state
        else:
            for qudit_label in initial_states:
                q = Qudit(qudit_label)
                state = State({q.uid: initial_states[qudit_label]})
                self.qudits[q.uid] = q
                self.states[state.uid] = state
           
'''
###############################################################
# UTILS
###############################################################


    def isSuperposed(self):
        kets = self.getKets()

        if (len(kets) > 1):
            return True

        return False
    
    def getAmplitude(self, ket):
        return self.value.coeff(ket)

    def getKets(self):
        v = self.value
        syms = list(v.free_symbols)
        kets = []
        for sym in syms:
            if (hasattr(sym, 'label')):
                kets.append(sym)
        
        return kets
    
    def getExpressionKets(self, e):
        v = e
        syms = list(v.free_symbols)
        kets = []
        for sym in syms:
            if (hasattr(sym, 'label')):
                kets.append(sym)
        
        return kets

    #Assumption syngle qudit qudits
    #TODO should we have qudits and not expressions as input?
    def valueEq(self, s):
        s1 = self.value
        s2 = s
        
        kets1 = self.getExpressionKets(s1)
        kets2 = self.getExpressionKets(s2)
        
        if (set(kets1) != set(kets2)):
            return False
        
        ref1_arg = sp.arg(s1.coeff(kets1[0]))
        ref2_arg = sp.arg(s2.coeff(kets1[0]))
        for ket in kets1:
            s1_abs = sp.Abs(s1.coeff(ket))
            s1_arg = sp.arg(s1.coeff(ket))

            s2_abs = sp.Abs(s2.coeff(ket))
            s2_arg = sp.arg(s2.coeff(ket))

            s1_relphase = s2_relphase = 0
            if (ket != kets1[0]):
                s1_relphase =  sp.Abs(ref1_arg - s1_arg)
                s2_relphase =  sp.Abs(ref2_arg - s2_arg)

            if (s1_relphase != s2_relphase or s1_abs != s2_abs):
                return False
        
        return True

    def getQuditqudits(self, qudit):
        uids = []
        for n in self.qudits:
            if self.qudits[n].qudit == qudit:
                uids.append(self.qudits[n].uid)

        return uids

    def getQuditquditInCorrelation(self, qudit, correlation_uid):
        for n in self.qudits:
            if self.qudits[n].qudit == qudit:
                if (self.qudits[n].correlation_uid is None and correlation_uid is None) or (
                    self.qudits[n].correlation_uid == correlation_uid
                ):
                    return self.qudits[n].uid

        return None

    def getQuditCorrelations(self, qudit):
        uids = []
        for e in self.correlations:
            qudit_uids = self.correlations[e].qudit_uids
            for qudit_uid in qudit_uids:
                if self.qudits[qudit_uid].qudit == qudit:
                    uids.append(e)

        return uids

    def correlatequdit(self, qudit_uid, correlation_uid):
        # assume correlation and qudit exist
        if self.qudits[qudit_uid].correlation_uid == None:
            self.qudits[qudit_uid].correlation_uid = correlation_uid
            self.correlations[correlation_uid].qudit_uids.append(qudit_uid)

    def deletequdit(self, qudit_uid):
        # assuming qudits only belong to one element
        if qudit_uid in self.qudits.keys():
            e_uid = self.qudits[qudit_uid].correlation_uid

            self.qudits.pop(qudit_uid)

            if e_uid in self.correlations.keys():
                i = self.correlations[e_uid].qudit_uids.index(qudit_uid)
                self.correlations[e_uid].qudit_uids.pop(i)

    def deleteCorrelation(self, correlation_uid):
        try:
            for nid in self.correlations[correlation_uid].qudit_uids:
                self.qudits.pop(nid)

                self.correlations.pop(correlation_uid)
        except KeyError:
            pass

    def copyCorrelation(self, correlation_uid):
        copy_e = quditCorrelation(self.correlations[correlation_uid].weight)
        self.correlations[copy_e.uid] = copy_e

        for n in self.correlations[correlation_uid].qudit_uids:
            qudit = self.qudits[n]
            new_n = qudit(qudit.qudit,self.dimension)
            new_n.value = qudit.value
            self.qudits[new_n.uid] = new_n

            self.correlatequdit(new_n.uid, copy_e.uid)
        
        return copy_e.uid

###############################################################
# REWRITE
###############################################################

    #qudit_map=["q0","q1","q2"]
    def rewrite(self,rules,qudit_map):
        for rule in rules:
            for i in enumerate(qudit_map):
                # Get Qudit Sets
        


###############################################################
# MEASURE
###############################################################

    def measure(self, qudits):
        #TODO
        pass

    #TODO Calculate loss as well
    def postSelect(self, qudits, s):
        # TODO
        pass
'''