from pyquil import get_qc, Program
from pyquil.gates import H, CNOT, Z, MEASURE, S, X, Y, I, FENCE
from pyquil.api import local_forest_runtime
from pyquil.quilbase import Declare
import pyquil.paulis as Pauli
from pyquil.paulis import PauliTerm
import networkx as nx
from pyquil.quantum_processor import NxQuantumProcessor
from pyquil.noise import decoherence_noise_with_asymmetric_ro


def main():
    n = 4
    prog = Program(
        Declare("ro", "BIT", n),
        H(0),
        CNOT(0, 1),
        H(3),
        MEASURE(0, ("ro", 0)),
        MEASURE(1, ("ro", 1)),
    )
    qvm = get_qc("8q-qvm")
    isa = qvm.to_compiler_isa()
    qubits = sorted(int(k) for k in isa.qubits.keys())
    edges = [(q1, q2) for q1 in qubits for q2 in qubits if q2-q1 == 1 or q2-q1 == -1]
    # Build the NX graph
    topo = nx.from_edgelist(edges)
    # You would uncomment the next line if you have disconnected qubits
    # topo.add_nodes_from(qubits)
    quantum_processor = NxQuantumProcessor(topo)
    #quantum_processor.noise_model = decoherence_noise_with_asymmetric_ro(quantum_processor.to_compiler_isa())  # Optional
    qvm.compiler.quantum_processor = quantum_processor
    g = qvm.quantum_processor.qubit_topology().to_directed().edges()
    print(g)
    #print(g)
    #print([pair for pair in qvm.quantum_processor.qubit_topology().to_directed().edges() if any(p in [0,1,2,3] for p in pair)])
    #print(qvm.quantum_processor.qubit_topology().nodes)
    #print(cm)
    #qc = get_qc('1q-qvm')  # You can make any 'nq-qvm' this way for any reasonable 'n'
    executable = qvm.compile(prog)
    #print(prog.get_qubit_indices())
    #print()
    #print(executable.get_qubit_indices())
    result = qvm.run(executable)
    bitstrings = result.get_register_map()['ro']
    # %%
    # for p in prog.instructions:
        #print(p._InstructionMeta__name)
        #print(p.name)
        #if p._InstructionMeta__name == 'Declare':
            #continue
        #if hasattr(p, "qubits"):
            #print([q.index for q in p.qubits])
        #elif hasattr(p, "qubit"):
            #print([p.qubit.index])
        #else:
            #print(p)
            #raise Exception("No qubits")
    #print(prog.get_qubit_indices())
    results = qvm.run(qvm.compile(prog.wrap_in_numshots_loop(10))).get_register_map()['ro']
    #print(results)
    from collections import defaultdict
    def dicthis(array):
        dic = dict()
        for a in array:
            if tuple(k for k in a) in dic:
                dic[tuple(k for k in a)] += 1
            else:
                dic[tuple(k for k in a)] = 1
        return dic
    #print(dicthis(results))
    #print(prog)

    #%%
    prog += H(0)
    #print(prog)
    pauli = PauliTerm.from_list([("X", 0),("Y", 2), ("Z", 3), ("X", 4)])
    pauli1 = PauliTerm.from_list([(p,i) for i, p in enumerate("IX")])
    pauli2 = PauliTerm.from_list([(p,i) for i, p in enumerate("XX")])
    for i in pauli1:
        print(i)
    pauliliste= []
    pauliliste.append(PauliTerm.from_list([(p,i) for i, p in enumerate("XX")]))
    pauliliste.append(PauliTerm.from_list([(p,i) for i, p in enumerate("YX")]))
    pauliliste.append(PauliTerm.from_list([(p,i) for i, p in enumerate("IY")]))
    #print(pauli.pauli_string(range(1+max(pauli.get_qubits()))))
    #print(str(pauli))
    #print(Pauli.check_commutation(pauli_list=[pauli2], pauli_two=pauli1))
    # %%
    def conjugate_pauli_with_cliffords(pauli_term: PauliTerm, program: Program) -> PauliTerm:
        # Iterate through the program's instructions in reverse order
        for instruction in reversed(program.instructions):
            if instruction.name == "H":
                for qubit in instruction.qubits:
                    if pauli_term[qubit.index] == "X" or pauli_term[qubit.index] == "Z":
                        pauli_term *= PauliTerm("Y", qubit.index)
            elif instruction.name == "S" or instruction.name == "S^-1":
                for qubit in instruction.qubits:
                    if pauli_term[qubit.index] == "X" or pauli_term[qubit.index] == "Y":
                        pauli_term *= PauliTerm("Z", qubit.index)
            elif instruction.name == "CNOT":
                copy_term = pauli_term.copy()
                control_qubit = instruction.qubits[0].index
                target_qubit = instruction.qubits[1].index
                if copy_term[control_qubit] == "X" or copy_term[control_qubit] == "Y":
                    pauli_term *= PauliTerm("X", target_qubit)
                if copy_term[target_qubit] == "Z" or copy_term[target_qubit] == "Y":
                    pauli_term *= PauliTerm("Z", control_qubit)
        return pauli_term
    
    prog = Program()
    #prog += CNOT(0,1)
    prog += CNOT(0,1)
    pauli = PauliTerm.from_list([(p,i) for i, p in enumerate("XY")])
    #print(pauli.pauli_string(range(1+max(pauli.get_qubits()))))
    #print(conjugate_pauli_with_cliffords(pauli, prog).pauli_string(range(1+max(pauli.get_qubits()))))
    # %%









from sys import platform
if platform == "linux" or platform == "linux2":
    main()
    # linux
elif platform == "darwin":
    raise Exception("IOS Not supported")
    # OS X
elif platform == "win32":
    with local_forest_runtime():
        main()
    # Windows...
