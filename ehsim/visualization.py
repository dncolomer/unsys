import json
import numpy as np
import sympy as sp
import sympy.physics.quantum as spq

from ehsim.gates import H, X, SWAP, CX, CCX

_gate_names = [[H, "H"], [X, "X"]]
zero_ket = np.array([1, 0], dtype=np.complex_)
one_ket = np.array([0, 1], dtype=np.complex_)
plus_ket = np.array([1 / np.sqrt(2), 1 / np.sqrt(2)], dtype=np.complex_)
minus_ket = np.array([1 / np.sqrt(2), -1 / np.sqrt(2)], dtype=np.complex_)

#TODO This needs to be refactred in order to support rules
def quirkExport(hypergraph):
    """
    Expert a hypergraph simulation to Quirk.

    """
    if not hypergraph._record_gates:
        raise ValueError(
            "You cannot export a hypergraph to quirk unless `record_gates` was set to True in the Hypergraph constructor."
        )
    base_url = "https://algassert.com/quirk#circuit={}"

    gates = [
        [1 for _ in range(len(hypergraph))] for _ in range(len(hypergraph._gate_log))
    ]

    for i, gate in enumerate(hypergraph._gate_log):
        for (a, n) in _gate_names:
            if gate["gate"] is a:
                for q in gate["qubits"]:
                    gates[i][int(q[1:])] = n
                for q in gate["controls"]:
                    gates[i][int(q[1:])] = "•"
        else:
            raise ValueError(
                "Unsupported gate. (This is just a prototype export for now!)"
            )

    return base_url.format(json.dumps({"cols": gates}, ensure_ascii=False))

def cytoscapeExport(hypergraph):
    elements = {}
    elements["nodes"] = []
    elements["edges"] = []
    for i in hypergraph.nodes:
        # set precision to avoid awkward -zeros
        s = []
        for state in hypergraph.nodes[i].state:
            s.append(np.around(state, decimals=3))

        measured = ""
        if hypergraph.nodes[i].measured:
            measured = "[M]"

        lbl = (
            str(hypergraph.nodes[i].uid)
            + "  ["
            + str(hypergraph.nodes[i].qubit)
            + "]  "
            + str(np.around(s, 3))
            + " "
            + measured
        )
        nn = {}
        nn["data"] = {}
        nn["data"]["id"] = lbl

        if hypergraph.nodes[i].edge_uid is not None:
            euid = hypergraph.nodes[i].edge_uid
            l_euid = euid + "  " + str(np.around(hypergraph.edges[euid].weight, 3))
            nn["data"]["parent"] = l_euid

        elements["nodes"].append(nn)

    print(json.dumps(elements))

def simplifiedState(hypergraph, state):
    #This can be expensive!
    state = sp.simplify(state)
    state = sp.expand(state)

    if (hypergraph.stateEq(state,spq.Ket(0))):
        return "|0>"

    if (hypergraph.stateEq(state,spq.Ket(1))):
        return "|1>"

    if (hypergraph.stateEq(state,spq.Ket(0)/sp.sqrt(2) + spq.Ket(1)/sp.sqrt(2))):
        return "|+>"

    if (hypergraph.stateEq(state,spq.Ket(0)/sp.sqrt(2) - spq.Ket(1)/sp.sqrt(2))):
        return "|->"
    
    return str(state)

def normalizeState(state):
    k0 = state.coeff(spq.Ket(0))
    k1 = state.coeff(spq.Ket(1))

    s_w = sp.matrices.Matrix([k0, k1])
    norm = sp.sqrt(s_w.dot(s_w))

    if (k0 != 0 or k1 != 0):
        state = k0*spq.Ket(0)/norm + k1*spq.Ket(1)/norm
    
    return state

def print_raw(hypergraph, simplify=True, normalize=False, subs=[]):
    print("      ".join(str(ql) for ql in hypergraph.qubitLabels))
    print("------".join("--" for ql in hypergraph.qubitLabels))  
    phelper = []

    if (normalize):
        hypergraph.normalizeHypergraph()
    
    for e in hypergraph.edges:
        phelper = []
        for ql in hypergraph.qubitLabels:
            nid = hypergraph.getQubitNodeIdInEdge(ql,e)

            if (nid is not None):
                replaced = ' '
                if (hypergraph.nodes[nid].replaced):
                    replaced = '*'

                state = hypergraph.nodes[nid].state
                for sub in subs:
                    if (sp.symbols(sub) in state.free_symbols):
                        state = state.subs(sp.symbols(sub),subs[sub])
                
                if (normalize):
                    state = normalizeState(state)

                if (simplify):
                    phelper.append(simplifiedState(hypergraph,state))
                else:
                    phelper.append(state)
            else:
                phelper.append("N/A")
        
        if (len(phelper) != 0):
            amp = hypergraph.edges[e].weight
            for sub in subs:
                if (sp.symbols(sub) in amp.free_symbols):
                    amp = amp.subs(sp.symbols(sub),subs[sub])

            phelper.append("weight: "+str(amp))

        print("     ".join(str(x) for x in phelper))
    
    phelper = []
    for ql in hypergraph.qubitLabels:
        nid = hypergraph.getQubitNodeIdInEdge(ql,None)
        systemEmpty = True

        if (nid is not None):
            systemEmpty = False
            replaced = ' '
            if (hypergraph.nodes[nid].replaced):
                replaced = '*'

            state = hypergraph.nodes[nid].state
            for sub in subs:
                if (sp.symbols(sub) in state.free_symbols):
                    state = state.subs(sp.symbols(sub),subs[sub])
            
            if (normalize):
                state = normalizeState(state)

            if (simplify):
                phelper.append(simplifiedState(hypergraph,state))
            else:
                phelper.append(state)
        else:
            phelper.append("N/A")

    if (not systemEmpty):
        if (len(phelper) != 0):
                phelper.append("Not Entangled")

        print("     ".join(str(x) for x in phelper))
    
    print("      ".join("  " for ql in hypergraph.qubitLabels))
    #print(hypergraph.toStateVector())
    print("------".join("--" for ql in hypergraph.qubitLabels))
    print("      ".join("  " for ql in hypergraph.qubitLabels))
    print("      ".join("  " for ql in hypergraph.qubitLabels))