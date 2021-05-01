from itertools import combinations

import numpy as np
import sympy as sp
import sympy.physics.quantum as spq
import string

# global variables
state_nb = 0

digs = string.digits + string.ascii_letters

def getUID(prefix="id"):
    global state_nb

    state_nb = state_nb + 1
    return prefix + str(state_nb - 1)

class State:
    def __init__(self, qudit, d, initial_value= None, symbolic= False):
        self.uid = getUID("s")
        self.qudit = qudit
        self.dimension = d

        #TODO consistency checks with regards to initial value and dimensionality, etc.
        self.value = initial_value
        if (self.value is None):
            self.value = spq.Ket('0')
            if (symbolic):
                sv_size = d
                i = 0
                s = None
                while i < sv_size:
                    sym = sp.symbols(self.uid+'_'+str(i))
                    if (i == 0):
                        s = sym * spq.Ket(i)
                    else:
                        s += sym * spq.Ket(i)
                    
                    i += 1
                    
                self.value = s

        self.measured = False
        self.correlation_uid = None
        self.replaced = False

    def isSuperposed():
        nbKets = self.getKets()

        if (len(kets) > 1):
            return True

        return False
    
    def getAmplitude(ket):
        return self.value.coeff(ket)

    def getKets():
        v = self.value
        syms = list(v.free_symbols)
        kets = []
        for sym in syms:
            if (hasattr(sym, 'label')):
                kets.append(sym)
        
        return kets

    #Assumption syngle qudit states
    #TODO should we have states and not expressions as input?
    def valueEq(self, s):
        s1 = self.value
        s2 = s
        
        kets1 = s1.getKets()
        kets2 = s2.getKets()
        
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

class StateCorrelation:
    def __init__(self, weight):
        self.state_uids = []
        self.uid = getUID("sc")

        self.weight = weight

# sv_expr is in format Ket('0,1,2') |0,1,2>
# The comma separated format is needed for higher dimensions to be still
# readable and interpreted as integers
# TODO support things like |0,1,-2>, etc.
class StateSystem:
    def __init__(self, nb_qudits, d, sv_expr= None, symbolic= False):
        global state_nb

        self.dimension = d
        self.states = {}
        self.correlations = {}
        self.quditLabels = []

        self._gate_log = []
        self.nb_qudits = nb_qudits

        state_nb = 0

        for i in range(0, nb_qudits):
            self.quditLabels.append("q" + str(i))

        if sv_expr is not None:
            #TODO check for consistency with d
            syms = list(sv_expr.free_symbols)
            for sym in syms:
                if (hasattr(sym, 'label')):
                    lbl = str(sym.label[0])
                    qubitlbls = lbl.split(',')
                    if (len(qubitlbls) == nb_qudits):
                        coeff = sv_expr.coeff(sym)

                        sc = StateCorrelation(coeff)
                        self.correlations[sc.uid] = sc

                        j = 0
                        while j < nb_qudits:
                            state = State("q"+str(j),self.dimension,symbolic=symbolic)
                            state.value = spq.Ket(qubitlbls[j])

                            self.states[state.uid] = state

                            self.correlateState(state.uid, sc.uid)
                            j += 1
                    else:
                        print("The provided wavefunction is inconsistent with the system's number of qudits.")
                        print("The affected Kets in the wavefunction have been ignored.")
        else:
            sc = StateCorrelation(1)
            self.correlations[sc.uid] = sc

            for i in range(0, nb_qudits):
                state = State("q" + str(i), self.dimension, symbolic=symbolic)
                self.states[state.uid] = state

                self.correlateState(state.uid,sc.uid)

###############################################################
# UTILS
###############################################################

    def getQubitStates(self, qubit):
        uids = []
        for n in self.states:
            if self.states[n].qudit == qubit:
                uids.append(self.states[n].uid)

        return uids

    def getQubitStateInCorrelation(self, qubit, correlation_uid):
        for n in self.states:
            if self.states[n].qudit == qubit:
                if (self.states[n].correlation_uid is None and correlation_uid is None) or (
                    self.states[n].correlation_uid == correlation_uid
                ):
                    return self.states[n].uid

        return None

    def getQubitCorrelations(self, qubit):
        uids = []
        for e in self.correlations:
            state_uids = self.correlations[e].state_uids
            for state_uid in state_uids:
                if self.states[state_uid].qudit == qubit:
                    uids.append(e)

        return uids

    def correlateState(self, state_uid, correlation_uid):
        # assume correlation and state exist
        if self.states[state_uid].correlation_uid == None:
            self.states[state_uid].correlation_uid = correlation_uid
            self.correlations[correlation_uid].state_uids.append(state_uid)

    def moveCorrelation(self, state_uid, src_correlation_uid, target_correlation_uid):
        # assign to target
        self.states[state_uid].correlation_uid = target_correlation_uid
        #pop from src
        self.correlations[src_correlation_uid].state_uids.pop(self.correlations[src_correlation_uid].state_uids.index(state_uid))
        #append to target
        self.correlations[target_correlation_uid].state_uids.append(state_uid)

    def deleteState(self, state_uid):
        # assuming states only belong to one element
        if state_uid in self.states.keys():
            e_uid = self.states[state_uid].correlation_uid

            self.states.pop(state_uid)

            if e_uid in self.correlations.keys():
                i = self.correlations[e_uid].state_uids.index(state_uid)
                self.correlations[e_uid].state_uids.pop(i)

    def deleteCorrelation(self, correlation_uid):
        for nid in self.correlations[correlation_uid].state_uids:
            self.states.pop(nid)

        self.correlations.pop(correlation_uid)

    def copyCorrelation(self, correlation_uid):
        copy_e = StateCorrelation(self.correlations[correlation_id].weight)
        self.correlations[copy_e.uid] = copy_e

        for n in self.correlations[correlation_uid]:
            state = self.states[n]
            new_n = State(state.qudit,self.dimension)
            new_n.value = state.value

            self.correlateState(new_n.uid, copy_e.uid)

        pass

    def composedCorrelations(self, qubits):
        eids_base = []

        for i,q in enumerate(qubits):
            if (i == 0):
                eids_base = self.getQubitCorrelations(q)
            else:
                if (set(eids_base) != set(self.getQubitCorrelations(q))):
                    return None

        return eids_base

    def areComposed(self, qubits):
        eids_base = []

        for i,q in enumerate(qubits):
            if (i == 0):
                eids_base = self.getQubitCorrelations(q)
            else:
                if (set(eids_base) != set(self.getQubitCorrelations(q))):
                    return False

        return True

###############################################################
# MERGE STATES
###############################################################

    def mergeQubitState(self, qubit):
        new_e = StateCorrelation(1)
        self.correlations[new_e.uid] = new_e

        new_n = State(qubit,self.dimension)
        e_delete = []
        for i,eid in enumerate(self.getQubitCorrelations(qubit)):
            if (len(self.correlations[eid].state_uids) != 1):
                print("Can't merge qubit "+qubit)
                return None
            else:
                w = self.correlations[eid].weight
                n = self.states[self.getQubitStateInCorrelation(eid)]

                if (i == 0):
                    new_n.value = n.value * w
                else:
                    new_n.value = new_n.value + (n.value * w)
                
                e_delete.append(eid)
        
        for eid in e_delete:
            self.deleteCorrelation(eid)

        self.correlateState(new_n.uid,new_e.uid)

        return new_e.uid  

    #This is where we merge single system states into one correlation superposing the comp. basis 
    def mergeQubitStates(self, qubits):
        for q in qubits:
            self.mergeQubitState(q)

###############################################################
# SPLIT STATES
###############################################################

    #This is where we split single system states into one correlation per computational base 
    def splitQubitState(self, qubit):
        correlation_ids = self.getQubitCorrelations(qubit)
        for eid in correlation_ids:
            state = self.getQubitStateInCorrelation(qubit,eid)
            if (state.isSuperposed()):
                kets = state.getKets()
                for ket in kets:
                    eid_copy = self.copyCorrelation(eid)
                    ampl = state.getAplitude(ket)

                    #Update Correlation Copy
                    state_copy_uid = self.getQubitStateInCorrelation(qubit,eid_copy)
                    self.states[state_copy_uid].value = ket
                    self.correlations[eid_copy].weight *= ampl
        
        #Clean-up old correlations
        for corr_id in correlation_ids:
            self.deleteCorrelation(corr_id)

    #This is where we split single system states into one correlation per computational base 
    def splitQubitsStates(self, qubits):
        for q in qubits:
            self.splitQubitState(q)

###############################################################
# SIMPLIFY
###############################################################    

    def canSimplifyCorrelations(self, base_e, cand_e, qubits):
        for q in qubits:
            base_n = self.getQubitStateInCorrelation(q,base_e)
            cand_n = self.getQubitStateInCorrelation(q,cand_e)

            if (not self.states[base_n].valueEq(self.states[cand_n].value)):
                return False

        return True

    def simplifyCorrelations(self, base_e, cand_e, qubits):
        self.correlations[base_e].weight = self.correlations[base_e].weight + self.correlations[cand_e].weight  
        self.deleteCorrelation(cand_e)

    def simplifyRec(self, eids, qubits):
        if (len(eids)) <= 1:
            return eids

        for i, base_e in enumerate(eids):
            for j, cand_e in enumerate(eids):
                if (i != j and self.canSimplifyCorrelations(base_e, cand_e, qubits)):
                    self.simplifyCorrelations(base_e, cand_e, qubits)
                    
                    #Update the list of correlation ids
                    eids.remove(cand_e)

                    return self.simplifyRec(eids, qubits)
        
        return eids

    #returns list of new correlations
    def simplifyQubitCorrelations(self, qubits):
        # Check if we can decompose
        if (self.areComposed(qubits)):
            eids_base = self.composedCorrelations(qubits)
            return self.simplifyRec(eids_base, qubits)
        else:
            print("Error: Qubits are not composed.")
        
        return []


###############################################################
# COMPOSE
###############################################################
    
    def composeCorrelations(self, src_g, target_g, qubits):
        for src_e in src_g:
            for target_e in target_g:
                for state_id in self.correlations[src_e].state_uids:
                    self.moveCorrelation(state_id,src_e,target_e)

            self.correlations[target_e].weight = self.correlations[target_e].weight * self.correlations[src_e].weight  
            self.deleteCorrelation(src_e)

    def canComposeCorrelations(self, base_g, cand_g, qubits):
        return not bool(set(base_g) & set(cand_g))

    def composeRec(self, correlation_groups, qubits): 
        print(correlation_groups)    
        base_group = []
        for i,el in enumerate(correlation_groups):
            if (i == 0):
                base_group = el
            else:
                if (self.canComposeCorrelations(base_group[i], base_group[i-1], qubits)):
                    correlation_id = self.compose(base_group[i], base_group[i-1], qubits)
                    correlation_groups.pop(i)
                    correlation_groups.pop(i-1)
                    correlation_groups.append([correlation_id])

                    return self.composeRec(correlation_groups, qubits)
        
        return correlation_groups
    
    def composeQubitCorrelations(self, qubits):
        correlation_groups = []
        for q in qubits:
            correlation_groups.append(self.getQubitCorrelations(q))
        
        return self.composeRec(correlation_groups,qubits)

###############################################################
# DECOMPOSE
###############################################################

    def decomposeCorrelations(self, base_e, cand_e, qubits):
        new_e = StateCorrelation(1)
        self.correlations[new_e.uid] = new_e

        for q in qubits:
            base_n = self.getQubitStateInCorrelation(q,base_e)
            cand_n = self.getQubitStateInCorrelation(q,cand_e)

            self.moveCorrelation(base_n,base_e,new_e.uid)
            self.deleteState(cand_n)
            

    def canDecomposeCorrelations(self, base_e, cand_e, qubits):
        for q in qubits:
            base_n = self.getQubitStateInCorrelation(q,base_e)
            cand_n = self.getQubitStateInCorrelation(q,cand_e)

            if (not self.states[base_n].valueEq(self.states[cand_n].value)):
                return False

        return True

    def decomposeRec(self, eids, qubits):
        if (len(eids)) <= 1:
            return eids

        for i, base_e in enumerate(eids):
            for j, cand_e in enumerate(eids):
                if (i != j and self.canDecomposeCorrelations(base_e, cand_e, qubits)):
                    self.decomposeCorrelations(base_e, cand_e, qubits)
                    
                    #Update the list of correlation ids
                    eids.remove(base_e)
                    eids.remove(cand_e)

                    return self.decomposeRec(eids, qubits)
        
        return eids

    #returns list of new correlations
    def decomposeQubits(self, qubits):
        # Check if we can decompose
        if (self.areComposed(qubits)):
            eids_base = self.composedCorrelations(qubits)
            return self.decomposeRec(eids_base, qubits)
        else:
            print("Error: Qubits are not composed.")
        
        return []

###############################################################
# REWRITE
###############################################################

    # match [1,1,0]
    # correlation {uid = ...}
    # map ["q0","q1","q2"]
    def isMatch(self,match,correlation,qubit_map):
        euid = None
        if (correlation is not None):
            euid = self.correlations[correlation].uid

        for i,m in enumerate(match):
            q = qubit_map[i]

            if (self.getQubitStateInCorrelation(q,euid) is None):
                return False
            else:
                state = self.states[self.getQubitStateInCorrelation(q,euid)]
                if (not state.valueEq(m)):
                    return False

        return True
    
    def replaceMatch(self,correlation,replace,qubit_map):
        euid = None
        if (correlation is not None):
            euid = self.correlations[correlation].uid

        for i,m in enumerate(replace):
            q = qubit_map[i]

            if (not self.states[self.getQubitStateInCorrelation(q,euid)].replaced):
                self.states[self.getQubitStateInCorrelation(q,euid)].value = m

                m = sp.simplify(m)
                m = sp.expand(m)

                self.states[self.getQubitStateInCorrelation(q,euid)].replaced = True

    #qubit_map=["q0","q1","q2"]
    def rewrite(self,rules,qubit_map,params_map=[]):
        if (self.areComposed(qubit_map)):
            #set all replacement tracking t False
            for n in self.states:
                self.states[n].replaced = False

            for rule in rules['rules']:
                #find matches in correlations
                for e in self.correlations:
                    if (self.isMatch(rule['match'],e,qubit_map)):
                        self.replaceMatch(e,rule['replace'],qubit_map)
                
                #find matches in system
                if (self.isMatch(rule['match'],None,qubit_map)):
                        self.replaceMatch(None,rule['replace'],qubit_map)
        else:
            print("Can't rewrite. Qubits are not composed")

###############################################################
# MEASURE
###############################################################

    def measure(self, qubits):
        # iterate over each hyper correlation the qubit is in
        self.splitQubitsStates(qubits)

        for qubit in qubits:
            states = self.getQubitStates(qubit)
            for state_uid in states:
                self.states[state_uid].measured = True

    # TODO calculate how much we are omitting (like Quirk does)
    # Calculate loss
    def postSelectZ(self, qubits, s):
        self.measure(qubits)
        for qubit in qubits:
            stateIds = self.getQubitStates(qubit)
            for stateId in stateIds:
                state = self.states[stateId]
                if (not state.valueEq(s)):
                    self.deleteCorrelation(state.correlation_uid)