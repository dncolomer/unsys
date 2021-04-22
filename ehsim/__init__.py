from itertools import combinations

import numpy as np
import sympy as sp
import sympy.physics.quantum as spq
import string

# global variables
node_nb = 0

digs = string.digits + string.ascii_letters

def int2base(x, base):
    if x < 0:
        sign = -1
    elif x == 0:
        return digs[0]
    else:
        sign = 1

    x *= sign
    digits = []

    while x:
        digits.append(digs[int(x % base)])
        x = int(x / base)

    if sign < 0:
        digits.append('-')

    digits.reverse()

    return ''.join(digits)

def getUID(prefix="id"):
    global node_nb

    node_nb = node_nb + 1
    return prefix + str(node_nb - 1)

class State:
    def __init__(self, qudit, d, initial_value= None, symbolic= True):
        self.uid = getUID("n")
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
                    i_base_d = int2base(i,d)
                    sym = sp.symbols(self.uid+'_'+str(i_base_d))
                    if (i == 0):
                        s = sym * spq.Ket(i_base_d)
                    else:
                        s += sym * spq.Ket(i_base_d)
                    
                    i += 1
                    
                self.value = s

        self.measured = False
        self.edge_uid = None
        self.replaced = False

class StateCorrelation:
    def __init__(self, weight, uid=None):
        self.node_uids = []
        if uid is None:
            self.uid = getUID("edge")
        else:
            self.uid = uid
        self.weight = weight

#sv_expr is in format Ket('012')
class StateSystem:
    def __init__(self, nb_qudits, d, sv_expr= None, symbolic= True, record_gates= True):
        global node_nb

        self.dimension = d
        self.nodes = {}
        self.edges = {}
        self.quditLabels = []

        self._record_gates = record_gates
        self._gate_log = []
        self.nb_qudits = nb_qudits

        #TODO better system for unique IDs
        node_nb = 0

        if sv_expr is not None:
            sv_size = d**nb_qudits
            i = 0
            while i < sv_size:
                i_base_d = int2base(i,d)
                coeff = sv_expr.coeff(i_base_d)
                e = StateCorrelation(coeff)

                j = 0
                while j < nb_qudits:
                    n = State("q"+str(j),self.dimension)
                    self.quditLabels.append("q"+str(j))

                    n.value = coeff*Ket(i_base_d[j])
                    self.addNodeToEdge(n.uid, e.uid)
                    j += 1
                
                i += 1
        else:
            for i in range(0, nb_qudits):
                self.quditLabels.append("q" + str(i))

                node = State("q" + str(i), self.dimension, symbolic=symbolic)
                self.nodes[node.uid] = node

###############################################################
# UTILS
###############################################################

    def stateEq(self, s1, s2):
        s10_abs = sp.Abs(s1.coeff(spq.Ket('0')))
        s10_arg = sp.arg(s1.coeff(spq.Ket('0')))

        s20_abs = sp.Abs(s2.coeff(spq.Ket('0')))
        s20_arg = sp.arg(s2.coeff(spq.Ket('0')))

        s11_abs = sp.Abs(s1.coeff(spq.Ket('1')))
        s11_arg = sp.arg(s1.coeff(spq.Ket('1')))

        s21_abs = sp.Abs(s2.coeff(spq.Ket('1')))
        s21_arg = sp.arg(s2.coeff(spq.Ket('1')))

        s1_relphase =  sp.Abs(s11_arg - s10_arg)
        s2_relphase =  sp.Abs(s21_arg - s20_arg)

        return s1_relphase == s2_relphase and s10_abs == s20_abs and s11_abs == s21_abs

    def _record(self, qubits, gate: np.ndarray, controls: list = None):
        """
        Make a note of which gates have been applied to which qubits.

        Arguments:
            qubits (str | list[str]): A list of qubits or a single qubit name
            gate: The gate that was applied
            controls (optional): A list of control qubits

        Returns:
            None

        """
        if not self._record_gates:
            return
        # If this gate acts on a single qubit, save it as an array of length=1.
        # This will make the gate log types consistent.
        if isinstance(qubits, str):
            qubits = [qubits]
        self._gate_log.append(
            {"gate": gate, "qubits": qubits, "controls": controls or []}
        )

    def __len__(self):
        return self._num_qubits

    def getQubitNodeIds(self, qubit):
        uids = []
        for n in self.nodes:
            if self.nodes[n].qudit == qubit:
                uids.append(self.nodes[n].uid)

        return uids

    def getQubitNodeIdInEdge(self, qubit, edge_uid):
        for n in self.nodes:
            if self.nodes[n].qudit == qubit:
                if (self.nodes[n].edge_uid is None and edge_uid is None) or (
                    self.nodes[n].edge_uid == edge_uid
                ):
                    return self.nodes[n].uid

        return None

    def getQubitNodeId(self, qubit, edge):
        uids = []
        for n in self.edges[edge].node_uids:
            if self.nodes[n].qudit == qubit:
                return self.nodes[n].uid

        return None

    def getQubitEdgeIds(self, qubit):
        uids = []
        for e in self.edges:
            node_uids = self.edges[e].node_uids
            for node_uid in node_uids:
                if self.nodes[node_uid].qudit == qubit:
                    uids.append(e)

        return uids

    def addNodeToEdge(self, node_uid, edge_uid):
        # assume edge and node exist
        if self.nodes[node_uid].edge_uid == None:
            self.nodes[node_uid].edge_uid = edge_uid
            self.edges[edge_uid].node_uids.append(node_uid)

    def moveNodeToEdge(self, node_uid, src_edge_uid, target_edge_uid):
        # assign to target
        self.nodes[node_uid].edge_uid = target_edge_uid
        #pop from src
        self.edges[src_edge_uid].node_uids.pop(self.edges[src_edge_uid].node_uids.index(node_uid))
        #append to target
        self.edges[target_edge_uid].node_uids.append(node_uid)

    def deleteNode(self, node_uid):
        # assuming nodes only belong to one element
        if node_uid in self.nodes.keys():
            e_uid = self.nodes[node_uid].edge_uid

            self.nodes.pop(node_uid)

            if e_uid in self.edges.keys():
                i = self.edges[e_uid].node_uids.index(node_uid)
                self.edges[e_uid].node_uids.pop(i)

    def deleteEdge(self, edge_uid):
        for nid in self.edges[edge_uid].node_uids:
            self.nodes.pop(nid)

        self.edges.pop(edge_uid)

    def copyEdge(self, edge_uid):
        copy_e = StateCorrelation(self.edges[edge_id].weight)
        self.edges[copy_e.uid] = copy_e

        for n in self.edges[edge_uid]:
            node = self.nodes[n]
            new_n = State(node.qudit,self.dimension)
            new_n.value = node.value

            self.addNodeToEdge(new_n.uid, copy_e.uid)

        pass

    def composedEdges(self, qubits):
        eids_base = []

        for i,q in enumerate(qubits):
            if (i == 0):
                eids_base = self.getQubitEdgeIds(q)
            else:
                if (set(eids_base) != set(self.getQubitEdgeIds(q))):
                    return None

        return eids_base

    def areComposed(self, qubits):
        eids_base = []

        for i,q in enumerate(qubits):
            if (i == 0):
                eids_base = self.getQubitEdgeIds(q)
            else:
                if (set(eids_base) != set(self.getQubitEdgeIds(q))):
                    return False

        return True

###############################################################
# MERGE STATES
###############################################################

    def mergeQubitState(self, qubit):
        new_e = StateCorrelation(1)
        new_n = State(qubit,self.dimension)
        e_delete = []
        for i,eid in enumerate(self.getQubitEdgeIds(qubit)):
            if (len(self.edges[eid].node_uids) != 1):
                print("Can't simplify qubit "+qubit)
                return None
            else:
                w = self.edges[eid].weight
                n = self.nodes[self.getQubitNodeIdInEdge(eid)]

                if (i == 0):
                    new_n.value = n.value * w
                else:
                    new_n.value = new_n.value + (n.value * w)
                
                e_delete.append(eid)
        
        for eid in e_delete:
            self.deleteEdge(eid)

        self.addNodeToEdge(new_n.uid,new_e.uid)

        return new_e.uid  

    #This is where we merge single system states into one edge superposing the comp. basis 
    def mergeQubitStates(self, qubits):
        for q in qubits:
            self.mergeQubitState(q)

###############################################################
# SPLIT STATES
###############################################################

    #This is where we split single system states into one edge per computational base 
    def splitQubitState(self, qubit):
        edge_ids = self.getQubitEdgeIds(qubit)
        for eid in edge_ids:
            node = self.getQubitNodeIdInEdge(qubit,eid)
            if (not self.stateEq(sqp.Ket('0'),node.value) and not self.stateEq(sqp.Ket('1'),node.value)):
                eid_copy = self.copyEdge(eid)
                coeff_0 = self.nodes[self.edges[eid]].value.coeff(spq.Ket('0'))
                coeff_1 = self.nodes[self.edges[eid]].value.coeff(spq.Ket('1'))
                self.nodes[self.edges[eid]].value = spq.Ket('0')
                self.nodes[self.edges[eid_copy]].value = spq.Ket('1')

                self.edges[eid].weight *= coeff_0
                self.edges[eid_copy].weight *= coeff_1 

    #This is where we split single system states into one edge per computational base 
    def splitQubitsStates(self, qubits):
        for q in qubits:
            self.splitQubitState(q)

###############################################################
# SIMPLIFY
###############################################################    

    def canSimplifyEdges(base_e, cand_e, qubits):
        for q in qubits:
            base_n = self.getQubitNodeIdInEdge(q,base_e)
            cand_n = self.getQubitNodeIdInEdge(q,cand_e)

            if (not self.stateEq(self.nodes[base_n].value, self.nodes[cand_e].value)):
                return False

        return True

    def simplifyEdges(base_e, cand_e, qubits):
        self.edges[base_e].weight = self.edges[base_e].weight + self.edges[cand_e].weight  
        self.deleteEdge(cand_e)

    def simplifyRec(eids, qubits):
        if (len(eids)) <= 1:
            return eids

        for i, base_e in enumerate(eids):
            for j, cand_e in enumerate(eids):
                if (i != j and self.canSimplifyEdges(base_e, cand_e, qubits)):
                    self.simplifyEdges(base_e, cand_e, qubits)
                    
                    #Update the list of edge ids
                    eids.remove(cand_e)

                    return self.simplifyRec(eids, qubits)
        
        return eids

    #returns list of new edges
    def simplifyQubits(self, qubits):
        # Check if we can decompose
        if (self.areComposed(qubits)):
            eids_base = self.composedEdges(qubits)
            return simplifyRec(eids_base, qubits)
        else:
            print("Error: Qubits are not composed.")
        
        return []


###############################################################
# COMPOSE
###############################################################
    
    def composeEdges(src_g, target_g, qubits):
        for src_e in src_g:
            for target_e in target_g:
                for node_id in self.edges[src_e].node_uids:
                    self.moveNodeToEdge(node_id,src_e,target_e)

            self.edges[target_e].weight = self.edges[target_e].weight * self.edges[src_e].weight  
            self.deleteEdge(src_e)

    def canComposeEdges(base_g, cand_g, qubits):
        return not bool(set(base_g) & set(cand_g))

    def composeRec(self, edge_groups, qubits):     
        base_group = []
        for i,el in enumerate(edge_groups):
            if (i == 0):
                base_group = base_group[el]
            else:
                if (self.canComposeEdges(base_group[i],base_group[i-1])):
                    edge_id = self.compose(base_group[i],base_group[i-1])
                    edge_groups.pop(i)
                    edge_groups.pop(i-1)
                    edge_groups.append([edge_id])

                    return self.composeRec(edge_groups, qubits)
        
        return edge_groups
    
    def composeQubits(qubits):
        edge_groups = []
        for q in qubits:
            edge_groups.append(getQubitEdgeIds(q))
        
        return self.composeRec(edge_groups,qubits)

###############################################################
# DECOMPOSE
###############################################################

    def decomposeEdges(base_e, cand_e, qubits):
        new_e = StateCorrelation(1)
        for q in qubits:
            base_n = self.getQubitNodeIdInEdge(q,base_e)
            cand_n = self.getQubitNodeIdInEdge(q,cand_e)

            self.moveNodeToEdge(base_n,base_e,new_e)
            self.deleteNode(cand_n)
            

    def canDecomposeEdges(base_e, cand_e, qubits):
        for q in qubits:
            base_n = self.getQubitNodeIdInEdge(q,base_e)
            cand_n = self.getQubitNodeIdInEdge(q,cand_e)

            if (not self.stateEq(self.nodes[base_n].value, self.nodes[cand_e].value)):
                return False

        return True

    def decomposeRec(eids, qubits):
        if (len(eids)) <= 1:
            return eids

        for i, base_e in enumerate(eids):
            for j, cand_e in enumerate(eids):
                if (i != j and self.canDecomposeEdges(base_e, cand_e, qubits)):
                    self.decomposeEdges(base_e, cand_e, qubits)
                    
                    #Update the list of edge ids
                    eids.remove(base_e)
                    eids.remove(cand_e)

                    return self.decomposeRec(eids, qubits)
        
        return eids

    #returns list of new edges
    def decomposeQubits(self, qubits):
        # Check if we can decompose
        if (self.areComposed(qubits)):
            eids_base = self.composedEdges(qubits)
            return decomposeRec(eids_base, qubits)
        else:
            print("Error: Qubits are not composed.")
        
        return []

###############################################################
# REWRITE
###############################################################

    # match [1,1,0]
    # edge {uid = ...}
    # map ["q0","q1","q2"]
    def isMatch(self,match,edge,qubit_map):
        euid = None
        if (edge is not None):
            euid = self.edges[edge].uid

        for i,m in enumerate(match):
            q = qubit_map[i]

            if (self.getQubitNodeIdInEdge(q,euid) is None):
                return False
            else:
                node = self.nodes[self.getQubitNodeIdInEdge(q,euid)]
                if (not self.stateEq(node.value,m)):
                    return False

        return True
    
    def replaceMatch(self,edge,replace,qubit_map):
        euid = None
        if (edge is not None):
            euid = self.edges[edge].uid

        for i,m in enumerate(replace):
            q = qubit_map[i]

            if (not self.nodes[self.getQubitNodeIdInEdge(q,euid)].replaced):
                self.nodes[self.getQubitNodeIdInEdge(q,euid)].value = m

                m = sp.simplify(m)
                m = sp.expand(m)

                self.nodes[self.getQubitNodeIdInEdge(q,euid)].replaced = True

    #qubit_map=["q0","q1","q2"]
    def rewrite(self,rules,qubit_map,params_map=[]):
        if (self.areComposed(qubit_map)):
            #set all replacement tracking t False
            for n in self.nodes:
                self.nodes[n].replaced = False

            for rule in rules['rules']:
                #find matches in edges
                for e in self.edges:
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
        # iterate over each hyper edge the qubit is in
        for q in qubits:
            edge_ids = self.getQubitEdgeIds(q)
            if len(edge_ids) == 0:
                # if the qubit is in the comp. basis then we just flag it as measured = True
                # remve any global phase
                # Note: Quirk keeps track of the phase so that a statevector (assuming measurement deferred can still be shown)
                # self.nodes[node.uid].value = self.correctPhase(node.value)
                node = self.nodes[self.getQubitNodeIdInEdge(q, None)]
                self.nodes[node.uid].measured = True
                if not self.stateEq(node.value, spq.Ket('1')) and not self.stateEq(
                    node.value, spq.Ket('0')
                ):
                    # if the qubit is not then we need to split the edge into 2.
                    # One where the node will be in the 0 state + measured = True
                    # One where the node will be in the 1 state + measured = True
                    self.splitEdgeZ(None, node.qudit)
            else:
                for e_id in edge_ids:
                    # if the qubit is in the comp. basis then we just flag it as measured = False
                    # remve any global phase
                    # Note: Quirk keeps track of the phase so that a statevector (assuming measurement deferred can still be shown)
                    # self.nodes[node.uid].value = self.correctPhase(node.value)
                    node = self.nodes[self.getQubitNodeIdInEdge(q, e_id)]
                    self.nodes[node.uid].measured = True
                    if not self.stateEq(node.value, spq.Ket('1')) and not self.stateEq(
                        node.value, spq.Ket('0')
                    ):
                        # if the qubit is not then we need to split the edge into 2.
                        # One where the node will be in the 0 state + measured = True
                        # One where the node will be in the 1 state + measured = True
                        self.splitEdgeZ(e_id, node.qudit)

        return

    # TODO calculate how much we are omitting (like Quirk does)
    def postSelectZ(self, qubits, state):
        self.measure(qubits)
        loss = 0
        for qubit in qubits:
            nodeIds = self.getQubitNodeIds(qubit)
            for nodeId in nodeIds:
                node = self.nodes[nodeId]
                if (not self.stateEq(node.value, state)):
                    if (node.edge_uid is not None):
                        edge = self.edges[node.edge_uid]
                        '''loss = loss + (
                            edge.weight.real ** 2 + edge.weight.imag ** 2
                        )'''

                    self.deleteEdge(node.edge_uid)

        return loss