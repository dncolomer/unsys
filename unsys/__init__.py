from itertools import combinations
from networkx.classes.function import neighbors

import numpy as np
import sympy as sp
import sympy.physics.quantum as spq
import hypernetx as hnx
from hnxwidget import HypernetxWidget

# global variables
guid = 0

def getUID(prefix="id#"):
    global guid

    guid = guid + 1
    return prefix + str(guid - 1)

class QuditSystem:
    def __init__(self, nb_qudits, dim, symbolic= False, hypergraph= None):
        global guid

        self.dim = dim
        self.nb_qudits = nb_qudits

        if (hypergraph is not None):
            self.hypergraph = hypergraph
        else:
            #init qudits
            qudits = []

            for qudit in range(nb_qudits):
                if (symbolic):
                    j = 0
                    state = 0
                    while j < dim:
                        coeff = sp.symbols('coeff_'+str(j))
                        state += coeff*spq.Ket(j)
                        
                        j += 1
                    qudits.append(hnx.Entity(getUID(),props={'qudit':qudit, 'state':state}))
                else:
                    qudits.append(hnx.Entity(getUID(),props={'qudit':qudit, 'state':spq.Ket(0)}))
            
            #Initialize Hypergraph
            syslabel = getUID('system#')
            self.hypergraph = hnx.Hypergraph({syslabel: hnx.Entity(syslabel,elements=qudits,props={'coeff':1})})

    def stateIntersection(self,state1,state2):
        s1_symbols = state1.free_symbols
        s2_symbols = state2.free_symbols

        return s1_symbols.intersection(s2_symbols)
        
    def draw(self):
        node_labels = {}
        for node in self.hypergraph.nodes():
            node_labels[node.uid] = 'qudit ' + str(node.props['qudit']) + ' | ' + str(node.props['state'])

        edge_labels = {}
        for edge in self.hypergraph.edges():
            edge_labels[edge.uid] = str(edge.uid) + ' | coeff: ' + str(edge.props['coeff'])

        hnx.draw(self.hypergraph,node_labels=node_labels,edge_labels=edge_labels)

    def drawBloch(self,qudit):
        pass

    def getStatevector(self,qudit):
        pass
    
    # subsystems is the list of qudits to consider
    def drawStatevector(self,subsystems):
        pass

    def getQuditProbabilities(self,qudit):
        mProb = []
        j = 0
        while j < self.dim:
            mProb[j] = 0
        
        return mProb

    # subsystems is the list of qudits to consider
    def drawProbabilities(self,subsystems):
        pass

    def getQuditNodes(self,qudit):
        node_list = []
        for e in self.hypergraph.nodes():
            if (e.props['qudit'] == qudit):
                node_list.append(e.uid)
        
        return node_list

    def cascadeNodeRemoval(self,hedges):
        if (len(hedges) == 0):
            return self.hypergraph

        next_hedges = []

        for h in hedges:
            hedge = hedges[h]

            #clean-up as many nodes as possible
            for n in hedge.children:
                node = hedge.children[n]
                if (len(node.memberships) <= 2): #all nodes are members of at least the system hedge
                    self.hypergraph.remove_node(node)
            
            #neighbour hedges
            edge_neighbors = self.hypergraph.edge_neighbors(hedge)
            for e in edge_neighbors:
                edge = edge_neighbors[e]
                if (len(edge.intersection(hedge)) == len(edge.children)):
                    #is nested: update coeff
                    edge.props['coeff'] = 0
                else:
                    #not nested
                    #remove intersecting nodes
                    self.hypergraph.remove_nodes(edge.intersection(hedge))
                    #queue edge
                    next_hedges.append(edge)

        return self.cascadeNodeRemoval(next_hedges) 

    # subsystems = list of qudits
    def postSelect(self,qudit,state):
        hg = self.hypergraph
        #delete all nodes where states don't overlap (share no computational basis elements)
        for e in self.hypergraph.nodes():
            state_intersection = self.stateIntersection(state,e.props['state'])
            if (e.props['qudit'] == qudit):
                if (len(state_intersection) == 0):
                    #fetch all the hyperedges this node belongs to
                    hedges = e.memberships

                    #remove node from hypergraph
                    self.hypergraph = self.hypergraph.remove_node(e.uid)

                    #update the rest of the hypergraph
                    self.cascadeNodeRemoval(hedges)
                    
                else:
                    #we colapse the noded state to <state>
                    #REVIEW
                    e.props['state'] = state


