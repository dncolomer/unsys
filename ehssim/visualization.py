import json
import numpy as np
import hypernetx as hnx

from .gates import *

_gate_names = [[H, "H"], [X, "X"]]


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
            l_euid = euid + "  " + str(np.around(hypergraph.edges[euid].amplitude, 3))
            nn["data"]["parent"] = l_euid

        elements["nodes"].append(nn)

    print(json.dumps(elements))


def draw(hypergraph):
    s = hnx.Entity("system", elements=[], amplitude=1)
    hg = hnx.Hypergraph()
    hg.add_edge(s)

    empty_system = True

    if len(hypergraph.nodes) == 0:
        print("Empty Hypergraph")

    for i in hypergraph.edges:
        edge = hnx.Entity(i + "  " + str(np.around(hypergraph.edges[i].amplitude, 3)))

        hg.add_edge(edge)

    for i in hypergraph.nodes:
        # set precision to avoid awkward -zeros
        s = []
        for state in hypergraph.nodes[i].state:
            s.append(np.around(state, decimals=3))

        measured = ""
        if hypergraph.nodes[i].measured:
            measured = "[M]"

        node = hnx.Entity(
            str(hypergraph.nodes[i].uid)
            + "  ["
            + str(hypergraph.nodes[i].qubit)
            + "]  "
            + str(np.around(s, 3))
            + " "
            + measured
        )

        if hypergraph.nodes[i].edge_uid is None:
            hg.add_node_to_edge(node, "system")
            empty_system = False
        else:
            euid = hypergraph.nodes[i].edge_uid
            l_euid = euid + "  " + str(np.around(hypergraph.edges[euid].amplitude, 3))
            hg.add_node_to_edge(node, l_euid)

    if empty_system:
        hg.remove_edge("system")

    """hnx.drawing.rubber_band.draw(hg,
            node_labels=getNodeLabels(hg),
            edge_labels=getEdgeLabels(hg),
            node_radius=2.0,
            label_alpha=1.0
    )"""

    hnx.drawing.two_column.draw(hg)