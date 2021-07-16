import json
import numpy as np
import sympy as sp
import sympy.physics.quantum as spq

def printQuditSystem(quditSystem):
    for q in quditSystem.qudits:
        print("---- Qudit " + q + " ----")
        for s in quditSystem.qudits[q]:
            print("Amplitude: " + str(s.amplitude))
            print("Ket: " + str(s.ket))
            print("***")

        print("-------------------------")
    
    for cuid in quditSystem.correlations:
        print("---- Correlation " + cuid + " ----")
        for c in quditSystem.correlations[cuid]:
            for q in c.state_map:
                for s in c.state_map[q]:
                    print("Amplitude: " + str(s.amplitude))
                    print("Ket: " + str(s.ket))
                    print("Qudit: " + str(s.qudit))
                    print("***")

        print("-------------------------")