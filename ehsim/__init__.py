from itertools import combinations

import numpy as np
import sympy as sp
import sympy.physics.quantum as spq

# global variables
node_nb = 0

def getUID(prefix="id"):
    global node_nb

    node_nb = node_nb + 1
    return prefix + str(node_nb - 1)

class Node:
    def __init__(self, qubit, state=None, symbolic=True):
        self.uid = getUID("n")
        self.qubit = qubit

        self.state = state
        if (state is None):
            self.state = spq.Ket(0)
            if (symbolic):
                self.state = sp.symbols(self.uid+'_0')*spq.Ket(0) + sp.symbols(self.uid+'_1')*spq.Ket(1)

        self.measured = False
        self.edge_uid = None
        self.replaced = False

class Hyperedge:
    def __init__(self, weight, uid=None):
        self.node_uids = []
        if uid is None:
            self.uid = getUID("edge")
        else:
            self.uid = uid
        self.weight = weight

class Hypergraph:
    # TODO add support for register level operation in a later version
    def __init__(self, nb_qubits: int, sv: list = None, record_gates: bool = True, symbolic: bool = True):
        """
        Create a new Hypergraph simulator.

        Arguments:
            nb_qubits (int): The number of qubits to simulate
            sv (list): A list of state vectors to initialize with
            record_gates (bool: True): Whether to keep a log of gates as they
                are applied. This incurs a slight performance penalty but
                enables some useful tools like visualization.

        """
        global node_nb

        self.nodes = {}
        self.edges = {}
        self.qubitLabels = []

        sv = sv or []

        self._record_gates = record_gates
        self._gate_log = []
        self._num_qubits = nb_qubits

        node_nb = 0

        if len(sv) > 0:
            #TODO
            pass
        else:
            for i in range(0, nb_qubits):
                self.qubitLabels.append("q" + str(i))

                node = Node("q" + str(i), symbolic=symbolic)
                self.nodes[node.uid] = node

    def __len__(self):
        return self._num_qubits

    def getQubitNodeIds(self, qubit):
        uids = []
        for n in self.nodes:
            if self.nodes[n].qubit == qubit:
                uids.append(self.nodes[n].uid)

        return uids

    def getQubitNodeIdInEdge(self, qubit, edge_uid):
        for n in self.nodes:
            if self.nodes[n].qubit == qubit:
                if (self.nodes[n].edge_uid is None and edge_uid is None) or (
                    self.nodes[n].edge_uid == edge_uid
                ):
                    return self.nodes[n].uid

        return None

    def getQubitNodeId(self, qubit, edge):
        uids = []
        for n in self.edges[edge].node_uids:
            if self.nodes[n].qubit == qubit:
                return self.nodes[n].uid

        return None

    def getQubitEdgeIds(self, qubit):
        uids = []
        for e in self.edges:
            node_uids = self.edges[e].node_uids
            for node_uid in node_uids:
                if self.nodes[node_uid].qubit == qubit:
                    uids.append(e)

        return uids

    def addNodeToEdge(self, node_uid, edge_uid):
        # assume edge and node exist
        if self.nodes[node_uid].edge_uid == None:
            self.nodes[node_uid].edge_uid = edge_uid
            self.edges[edge_uid].node_uids.append(node_uid)

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

    # We must conbine them in distributabliy
    def combineEdges(self, a_edge_ids, b_edge_ids):
        for a_id in a_edge_ids:
            for b_id in b_edge_ids:
                a_edge = self.edges[a_id]
                b_edge = self.edges[b_id]

                e = Hyperedge(a_edge.weight * b_edge.weight)
                self.edges[e.uid] = e

                # recrerate the nodes of a inside the new edge
                for n_id in a_edge.node_uids:
                    p = Node(self.nodes[n_id].qubit, self.nodes[n_id].state)
                    p.measured = self.nodes[n_id].measured
                    self.nodes[p.uid] = p
                    self.addNodeToEdge(p.uid, e.uid)

                # recrerate the nodes of b inside the new edge
                for n_id in b_edge.node_uids:
                    p = Node(self.nodes[n_id].qubit, self.nodes[n_id].state)
                    p.measured = self.nodes[n_id].measured
                    self.nodes[p.uid] = p
                    self.addNodeToEdge(p.uid, e.uid)

        # Cleanup and delete the old nodes and edges
        for e_id in a_edge_ids:
            self.deleteEdge(e_id)

        for e_id in b_edge_ids:
            self.deleteEdge(e_id)

    def stateEq(self, s1, s2):
        s10_abs = sp.Abs(s1.coeff(spq.Ket(0)))
        s10_arg = sp.arg(s1.coeff(spq.Ket(0)))

        s20_abs = sp.Abs(s2.coeff(spq.Ket(0)))
        s20_arg = sp.arg(s2.coeff(spq.Ket(0)))

        s11_abs = sp.Abs(s1.coeff(spq.Ket(1)))
        s11_arg = sp.arg(s1.coeff(spq.Ket(1)))

        s21_abs = sp.Abs(s2.coeff(spq.Ket(1)))
        s21_arg = sp.arg(s2.coeff(spq.Ket(1)))

        s1_relphase =  sp.Abs(s11_arg - s10_arg)
        s2_relphase =  sp.Abs(s21_arg - s20_arg)

        return s1_relphase == s2_relphase and s10_abs == s20_abs and s11_abs == s21_abs

    # Transforms the full system
    def toStateVector(self):
        pass

    def canMergeEdges(self, euid1, euid2):
        nuids1 = self.edges[euid1].node_uids
        nuids2 = self.edges[euid2].node_uids
        diff = 0

        for nuid1 in nuids1:
            for nuid2 in nuids2:
                if (self.nodes[nuid1].qubit == self.nodes[nuid2].qubit and not self.stateEq(self.nodes[nuid1].state,self.nodes[nuid2].state)):
                    diff = diff + 1
                    if (diff > 1):
                        return False

        return diff <= 1

    #preconditin: we assume canMerge
    def mergeEdges(self, euid1, euid2):
        nuids1 = self.edges[euid1].node_uids
        nuids2 = self.edges[euid2].node_uids

        for nuid1 in nuids1:
            for nuid2 in nuids2:
                if (self.nodes[nuid1].qubit == self.nodes[nuid2].qubit and not self.stateEq(self.nodes[nuid1].state,self.nodes[nuid2].state)):
                    self.nodes[nuid1].state = (self.nodes[nuid1].state * self.edges[euid1].weight) + (self.nodes[nuid2].state * self.edges[euid2].weight)
                    #TODO update euid1 weight
        
        self.deleteEdge(euid2)

    def factorRec(self, edge_ids = [], measured_qubits=[], steps=None):
        if steps is not None and steps == 0:
            return edge_ids

        if len(m.edge_ids()) <= 1:
            return edge_ids

        if steps is not None and steps > 0:
            steps = steps - 1

        base_edge = edge_ids[0]
        for i, cand_e in enumerate(edge_ids):
            if (i > 0 and self.canMergeEdges(base_edge, cand_e)):
                self.mergeEdges(base_edge, cand_e)
                edge_ids.pop(0)
                edge_ids.pop(i)

                return factorRec(edge_ids,measured_qubits,steps)
        
        return edge_ids

    def splitAllEdgesZ(self, qubit):
        a_edge_ids = self.getQubitEdgeIds(qubit)
        if (len(a_edge_ids) != 0):
            for eid in a_edge_ids:
                self.splitEdgeZ(eid,a)
        else:
            self.splitEdgeZ(None,qubit)

    def splitEdgeZ(self, edge_uid, qubit):
        node = self.nodes[self.getQubitNodeIdInEdge(qubit, edge_uid)]

        if (not self.stateEq(node.state,spq.Ket(0)) and not self.stateEq(node.state,spq.Ket(1))):
            if edge_uid is None:
                # TODO This is a bit dirty ;)
                # create an edge for the 1 component
                edge = Hyperedge(1)
                self.edges[edge.uid] = edge

                # populate initial edge
                for n_id in self.nodes:
                    self.addNodeToEdge(n_id, edge.uid)
            else:
                edge = self.edges[edge_uid]

            # Create Edge for the 0 component
            e = Hyperedge(edge.weight * node.state.coeff(spq.Ket(0)))
            self.edges[e.uid] = e

            # recrerate the nodes of a inside the new edge
            for n_id in edge.node_uids:
                state = self.nodes[n_id].state
                if n_id == node.uid:
                    state = spq.Ket(0)

                p = Node(self.nodes[n_id].qubit, state)
                p.measured = self.nodes[n_id].measured
                self.nodes[p.uid] = p
                self.addNodeToEdge(p.uid, e.uid)

            # Update current edge to reflect the 1 component
            edge.weight = edge.weight * node.state.coeff(spq.Ket(1))
            for n_id in edge.node_uids:
                if n_id == node.uid:
                    self.nodes[n_id].state = spq.Ket(1)

    # TODO Measure a set of qubits
    # TODO Revise preconditions for the op methods
    # TODO Update th rest of the codebase to account for non measured nodes
    # TODO probably don't allow for any gates on non measured nodes (maybe bitflips?)
    # TODO we lose the measurement labels when factoring
    def measure(self, qubits):
        # iterate over each hyper edge the qubit is in
        for q in qubits:
            edge_ids = self.getQubitEdgeIds(q)
            if len(edge_ids) == 0:
                # if the qubit is in the comp. basis then we just flag it as measured = True
                # remve any global phase
                # Note: Quirk keeps track of the phase so that a statevector (assuming measurement deferred can still be shown)
                # self.nodes[node.uid].state = self.correctPhase(node.state)
                node = self.nodes[self.getQubitNodeIdInEdge(q, None)]
                self.nodes[node.uid].measured = True
                if not self.stateEq(node.state, one_ket) and not self.stateEq(
                    node.state, zero_ket
                ):
                    # if the qubit is not then we need to split the edge into 2.
                    # One where the node will be in the 0 state + measured = True
                    # One where the node will be in the 1 state + measured = True
                    self.splitEdgeZ(None, node.qubit)
            else:
                for e_id in edge_ids:
                    # if the qubit is in the comp. basis then we just flag it as measured = False
                    # remve any global phase
                    # Note: Quirk keeps track of the phase so that a statevector (assuming measurement deferred can still be shown)
                    # self.nodes[node.uid].state = self.correctPhase(node.state)
                    node = self.nodes[self.getQubitNodeIdInEdge(q, e_id)]
                    self.nodes[node.uid].measured = True
                    if not self.stateEq(node.state, one_ket) and not self.stateEq(
                        node.state, zero_ket
                    ):
                        # if the qubit is not then we need to split the edge into 2.
                        # One where the node will be in the 0 state + measured = True
                        # One where the node will be in the 1 state + measured = True
                        self.splitEdgeZ(e_id, node.qubit)

        return

    def normalizeHypergraph(self):
        #Edges
        w = []
        for euid in self.edges:
            w.append(self.edges[euid].weight)
        
        if (len(w) > 0):
            s_w = sp.matrices.Matrix([w])
            #the norm is the square root of the dot product of the vector with itself
            norm = sp.sqrt(s_w.dot(s_w))

            for euid in self.edges:
                self.edges[euid].weight = sp.simplify(self.edges[euid].weight / norm)
        
        #Nodes
        for n in self.nodes:
            state = sp.simplify(self.nodes[n].state)
            k0 = state.coeff(spq.Ket(0))
            k1 = state.coeff(spq.Ket(1))

            if (k0 != 0 and k1 != 0):
                s_w = sp.matrices.Matrix([k0, k1])
                norm = sp.sqrt(s_w.dot(s_w))

                self.nodes[n].state = k0*spq.Ket(0)/norm + k1*spq.Ket(1)/norm


    # TODO calculate how much we are omitting (like Quirk does)
    def postSelectZ(self, qubits, state=zero_ket):
        self.measure(qubits)
        loss = 0
        for qubit in qubits:
            nodeIds = self.getQubitNodeIds(qubit)
            for nodeId in nodeIds:
                node = self.nodes[nodeId]
                if not self.stateEq(node.state, state):
                    if node.edge_uid is not None:
                        edge = self.edges[node.edge_uid]
                        loss = loss + (
                            edge.weight.real ** 2 + edge.weight.imag ** 2
                        )

                    self.deleteEdge(node.edge_uid)

        return loss

    def expandQubits(self, a, b):
        #if one of the qubits is not entangled but the other is we need to add the qubit to all corresponding edges before we start with the operation (sort of decompress)
        #preprocessing NOT tested
        a_edge_ids = self.getQubitEdgeIds(a)
        b_edge_ids = self.getQubitEdgeIds(b)
        
        if (len(a_edge_ids) == 0 and len(b_edge_ids) != 0):
            a_node_ids = self.getQubitNodeIds(a)
            a_node = self.nodes[a_node_ids[0]]
            for edge_id in b_edge_ids:
                #create node
                p = Node(a_node.qubit,a_node.state)
                #add nodes to hgraph
                self.nodes[p.uid] = p
                #add nodes to edge
                self.addNodeToEdge(p.uid,edge_id)
            
            #delete original node
            #did nt belong t any edge anyway
            self.deleteNode(a_node.uid)
        
        elif (len(b_edge_ids) == 0 and len(a_edge_ids) != 0):
            b_node_ids = self.getQubitNodeIds(b)
            b_node = self.nodes[b_node_ids[0]]
            for edge_id in a_edge_ids:
                #create node
                p = Node(b_node.qubit,b_node.state)
                #add nodes to hgraph
                self.nodes[p.uid] = p
                #add nodes to edge
                self.addNodeToEdge(p.uid,edge_id)
            
            #delete original node
            #did nt belong t any edge anyway
            self.deleteNode(b_node.uid)
        
        #if both entangled but in different edges we need to do some preprocessing as well
        shared_edges = list(set(a) & set(b))

        #tbh, if edges are not shared they, the intersection must be empty. Any other set-up it's just not a valid state
        if (len(shared_edges) == 0):
            self.combineEdges(a_edge_ids,b_edge_ids)

    # Factor a specific set of entangled qubits
    # TODO add a verbose mode that explains a bit more what's going on (what's being merged)
    def factorQubits(self, qubits, steps=None, verbose=False):
        # preprocessing: we expand all the affected qubits to avoid global phase comparison issues
        for index, qubit in enumerate(qubits):
            if index > 0:
                self.expandQubits(qubits[index - 1], qubits[index])

        #Health check on input
        measured_qubits = []
        edge_ids = []
        for q in qubits:
            q_edge_ids = self.getQubitEdgeIds(q)
            if (len(edge_ids) == 0):
                edge_ids = self.getQubitEdgeIds(q)
            else:
                if (set(edge_ids) != set(q_edge_ids)):
                    print("The given qubits can't be factored")
                    return None

        # call recursive part
        return self.factorRec(edge_ids= edge_ids, measured_qubits= measured_qubits, steps= steps)
        
    #
    # match [1,1,0]
    # edge {uid = ...}
    # map ["q0","q1","q2"]
    def isMatch(self,match,edge,qubit_map):
        euid = None
        if (edge is not None):
            euid = self.edges[edge].uid

        for i,m in enumerate(match):
            q = qubit_map[i] #TODO check for error

            if (self.getQubitNodeIdInEdge(q,euid) is None):
                return False
            else:
                node = self.nodes[self.getQubitNodeIdInEdge(q,euid)]
                if (not self.stateEq(node.state,m)):
                    return False

        return True
    
    def replaceMatch(self,edge,replace,qubit_map):
        euid = None
        if (edge is not None):
            euid = self.edges[edge].uid

        for i,m in enumerate(replace):
            q = qubit_map[i] #TODO check for error

            if (not self.nodes[self.getQubitNodeIdInEdge(q,euid)].replaced):
                self.nodes[self.getQubitNodeIdInEdge(q,euid)].state = m

                m = sp.simplify(m)
                m = sp.expand(m)

                self.nodes[self.getQubitNodeIdInEdge(q,euid)].replaced = True

    #qubit_map=["q0","q1","q2"]
    def rewrite(self,rules,qubit_map,params_map=[]):
        #TODO need to refactor this too account for rules and qubit mappings
        #self._record(qubit, rules['name'])

        #set all replacement tracking t False
        for n in self.nodes:
            self.nodes[n].replaced = False

        #Expand & split mapped qubits
        for index, qubit in enumerate(qubit_map):
            if index > 0:
                self.expandQubits(qubit_map[index - 1], qubit_map[index])
        
        for qubit in qubit_map:
            self.splitAllEdgesZ(qubit)

        for rule in rules['rules']:
            #find matches in edges
            # TODO if th edge has been touched previously shouldnt match anymore
            for e in self.edges:
                if (self.isMatch(rule['match'],e,qubit_map)):
                    self.replaceMatch(e,rule['replace'],qubit_map)
            
            #find matches in system
            if (self.isMatch(rule['match'],None,qubit_map)):
                    self.replaceMatch(None,rule['replace'],qubit_map)
