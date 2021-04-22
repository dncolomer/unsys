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
def quirkExport(state_system):
    """
    Expert a state_system simulation to Quirk.

    """
    if not state_system._record_gates:
        raise ValueError(
            "You cannot export a state_system to quirk unless `record_gates` was set to True in the Hypergraph constructor."
        )
    base_url = "https://algassert.com/quirk#circuit={}"

    gates = [
        [1 for _ in range(len(state_system))] for _ in range(len(state_system._gate_log))
    ]

    for i, gate in enumerate(state_system._gate_log):
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

def cytoscapeExport(state_system):
    elements = {}
    elements["nodes"] = []
    elements["edges"] = []
    for i in state_system.nodes:
        # set precision to avoid awkward -zeros
        s = []
        for state in state_system.nodes[i].value:
            s.append(np.around(state, decimals=3))

        measured = ""
        if state_system.nodes[i].measured:
            measured = "[M]"

        lbl = (
            str(state_system.nodes[i].uid)
            + "  ["
            + str(state_system.nodes[i].qubit)
            + "]  "
            + str(np.around(s, 3))
            + " "
            + measured
        )
        nn = {}
        nn["data"] = {}
        nn["data"]["id"] = lbl

        if state_system.nodes[i].edge_uid is not None:
            euid = state_system.nodes[i].edge_uid
            l_euid = euid + "  " + str(np.around(state_system.edges[euid].weight, 3))
            nn["data"]["parent"] = l_euid

        elements["nodes"].append(nn)

    print(json.dumps(elements))

# Transforms the full system
def toStateVector(self):
    pass

def simplifiedState(state_system, state):
    #This can be expensive!
    state = sp.simplify(state)
    state = sp.expand(state)

    if (state_system.stateEq(state,spq.Ket('0'))):
        return "|0>"

    if (state_system.stateEq(state,spq.Ket('1'))):
        return "|1>"

    if (state_system.stateEq(state,spq.Ket('0')/sp.sqrt(2) + spq.Ket('1')/sp.sqrt(2))):
        return "|+>"

    if (state_system.stateEq(state,spq.Ket('0')/sp.sqrt(2) - spq.Ket('1')/sp.sqrt(2))):
        return "|->"
    
    state = sp.simplify(state)
    
    return str(state)

def printState(state, simplify=True, subs=[]):
    #TODO
    print(state.value)

def printStateSystem(state_system, simplify=True, subs=[]):
    print("      ".join(str(ql) for ql in state_system.qubitLabels))
    print("------".join("--" for ql in state_system.qubitLabels))  
    phelper = []
    
    for e in state_system.edges:
        phelper = []
        for ql in state_system.qubitLabels:
            nid = state_system.getQubitNodeIdInEdge(ql,e)

            if (nid is not None):
                replaced = ' '
                if (state_system.nodes[nid].replaced):
                    replaced = '*'

                state = state_system.nodes[nid].value
                for sub in subs:
                    if (sp.symbols(sub) in state.free_symbols):
                        state = state.subs(sp.symbols(sub),subs[sub])

                if (simplify):
                    phelper.append(simplifiedState(state_system,state))
                else:
                    phelper.append(state)
            else:
                phelper.append("N/A")
        
        if (len(phelper) != 0):
            amp = state_system.edges[e].weight
            for sub in subs:
                if (sp.symbols(sub) in amp.free_symbols):
                    amp = amp.subs(sp.symbols(sub),subs[sub])

            phelper.append("weight: "+str(amp))

        print("     ".join(str(x) for x in phelper))
    
    phelper = []
    for ql in state_system.qubitLabels:
        nid = state_system.getQubitNodeIdInEdge(ql,None)
        systemEmpty = True

        if (nid is not None):
            systemEmpty = False
            replaced = ' '
            if (state_system.nodes[nid].replaced):
                replaced = '*'

            state = state_system.nodes[nid].value
            for sub in subs:
                if (sp.symbols(sub) in state.free_symbols):
                    state = state.subs(sp.symbols(sub),subs[sub])

            if (simplify):
                phelper.append(simplifiedState(state_system,state))
            else:
                phelper.append(state)
        else:
            phelper.append("N/A")

    if (not systemEmpty):
        if (len(phelper) != 0):
                phelper.append("Not Entangled")

        print("     ".join(str(x) for x in phelper))
    
    print("      ".join("  " for ql in state_system.qubitLabels))
    #print(state_system.toStateVector())
    print("------".join("--" for ql in state_system.qubitLabels))
    print("      ".join("  " for ql in state_system.qubitLabels))
    print("      ".join("  " for ql in state_system.qubitLabels))