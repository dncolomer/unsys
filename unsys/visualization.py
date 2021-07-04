import json
import numpy as np
import sympy as sp
import sympy.physics.quantum as spq

def printQuditSystem(quditSystem):
    for s in quditSystem.states:
        print("------------" + s + "------------")
        for sm in quditSystem.states[s].state_map:
            print(sm + ": " + str(quditSystem.states[s].state_map[sm]))

'''
def printStateSystem(state_system, subs=[]):
    print("      ".join(str(ql) for ql in state_system.quditLabels))
    print("------".join("--" for ql in state_system.quditLabels))  
    phelper = []
    
    for e in state_system.correlations:
        phelper = []
        for ql in state_system.quditLabels:
            nid = state_system.getQuditStateInCorrelation(ql,e)

            if (nid is not None):
                replaced = ' '
                if (state_system.states[nid].replaced):
                    replaced = '*'

                state = state_system.states[nid].value
                for sub in subs:
                    if (sp.symbols(sub) in state.free_symbols):
                        state = state.subs(sp.symbols(sub),subs[sub])

                if (state_system.states[nid].measured):
                    phelper.append(str(state) + ' M')
                else:
                    phelper.append(state)
            else:
                phelper.append("N/A")
        
        if (len(phelper) != 0):
            amp = state_system.correlations[e].weight
            for sub in subs:
                if (sp.symbols(sub) in amp.free_symbols):
                    amp = amp.subs(sp.symbols(sub),subs[sub])

            phelper.append("weight: "+str(amp))

        print("     ".join(str(x) for x in phelper))
    
    print("      ".join("  " for ql in state_system.quditLabels))
    #print(state_system.toStateVector())
    print("------".join("--" for ql in state_system.quditLabels))
    print("      ".join("  " for ql in state_system.quditLabels))
    print("      ".join("  " for ql in state_system.quditLabels))
'''