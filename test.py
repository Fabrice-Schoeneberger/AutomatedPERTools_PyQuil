from main import executor, get_backend
from pyquil import Program
from pyquil.gates import H, CNOT, Z, MEASURE, S, X, Y, I, RX, RZ, FENCE
from pyquil.quilbase import Declare
import numpy as np


# Here the backend for the simulation is prepared
from pyquil.quantum_processor import NxQuantumProcessor
from pyquil.noise import NoiseModel
from pyquil import get_qc
import networkx as nx
from pyquil.quil import DefGate

s_dagger = np.array([[ 1+0j,  0+0j],
                   [ 0+0j,  0-1j]])
# Get the Quil definition for the new gate
s_dagger_definition = DefGate("S_DAGGER", s_dagger)
# Get the gate constructor
S_DAGGER = s_dagger_definition.get_constructor()

backend = get_qc("5q-qvm") #str(n)+'q-qvm' #args.backend
isa = backend.to_compiler_isa()
backend_qubits = sorted(int(k) for k in isa.qubits.keys())
# By default every qubit on our fake backend is connected with every other qubit. This breaks the algorhythm. So we need to change the topology
# So this line here changes the topology of the fake backend. The condition in the end determines it
edges = [(q1, q2) for q1 in backend_qubits for q2 in backend_qubits if abs(q1-q2) == 1]
# Build the NX graph
topo = nx.from_edgelist(edges)
# You would uncomment the next line if you have disconnected qubits
# topo.add_nodes_from(qubits)
nx_quantum_processor = NxQuantumProcessor(topo)
backend.compiler.quantum_processor = nx_quantum_processor

#         ┌───┐                ░ ┌───┐      ░ ┌───┐ ┌───┐            ░ ┌─┐   
#q_0 -> 0 ┤ X ├────────────■───░─┤ Y ├──■───░─┤ Z ├─┤ X ├────────────░─┤M├───
#         ├───┤┌───┐┌───┐┌─┴─┐ ░ ├───┤┌─┴─┐ ░ ├───┤┌┴───┴┐┌───┐┌───┐ ░ └╥┘┌─┐
#q_1 -> 1 ┤ H ├┤ S ├┤ I ├┤ X ├─░─┤ X ├┤ X ├─░─┤ I ├┤ Sdg ├┤ H ├┤ I ├─░──╫─┤M├
#         └───┘└───┘└───┘└───┘ ░ └───┘└───┘ ░ └───┘└─────┘└───┘└───┘ ░  ║ └╥┘
# meas: 2/══════════════════════════════════════════════════════════════╩══╩═
#                                                                       0  1 
#{'11': 2, '00': 44, '10': 44, '01': 934}

prog = Program()
prog += Declare("ro", "BIT", 2)
prog += X(0)
prog += H(1)
prog += S(1)
prog += I(1)
prog += CNOT(0,1)
prog += Y(0)
prog += X(1)
prog += CNOT(0,1)
prog += Z(0)
prog += I(1)
prog += X(0)
prog += s_dagger_definition
prog += S_DAGGER(1)
prog += H(1)
prog += I(1)
prog += MEASURE(0, ("ro", 0))
prog += MEASURE(1, ("ro", 1))

prog2 = Program()
prog2 += Declare("ro", "BIT", 2)
prog2 += X(0)
prog2 += H(1)
prog2 += S(1)
#prog2 += I(1)
prog2 += CNOT(0,1)
prog2 += s_dagger_definition
prog2 += S_DAGGER(1)
prog2 += H(1)
prog2 += MEASURE(0, ("ro", 0))
prog2 += MEASURE(1, ("ro", 1))
print(prog2)

prog4 = Program()
prog4 += Declare("ro", "BIT", 2)
prog4 += CNOT(0,1)
prog4 += S(0)
prog4 += H(0)
prog4 += MEASURE(0, ("ro", 0))
prog4 += MEASURE(1, ("ro", 1))
#print(prog4)

twoqubit_errorops = ['IX']
twoqubit_errorprobs = [0.05]
twoqubit_error_template = [(op, p) for op,p in zip(twoqubit_errorops, twoqubit_errorprobs)]+[("II", 1-sum(twoqubit_errorprobs))]
#(noise_model, twoqubit_error_template, singlequbit_error_template) = get_noise_model()
#prog.define_noisy_gate("X", [0], [identity])
#prog.define_noisy_gate("CNOT", [1,0], [pauli_IX])
#apply_noise_model(prog, backend, twoqubit_error_template, "CNOT")

counts = executor([prog], backend, 10024, apply_noise=True)
print(counts)
counts = executor([prog4], backend, 1024, apply_noise=True)
print(counts)
