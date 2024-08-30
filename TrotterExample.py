def main():
    import argparse
    parser = argparse.ArgumentParser()
        
    # Definiere ein Argument
    parser.add_argument('--plusone', '-p', help='Takes Neighboring qubits into account', default=False, action='store_true')
    parser.add_argument('--sum', '-s', help='Same as -p and turns sumation on over neighboring qubits', default=False, action='store_true')
    parser.add_argument('--pntsamples', type=int, help='How many samples in PNT? Default: 16', default=16)
    parser.add_argument('--pntsinglesamples', type=int, help='How many single samples in PNT? Default: 100', default=100)
    parser.add_argument('--persamples', type=int, help='How many samples in PER? Default: 100', default=100)
    parser.add_argument('--shots', type=int, help='How many shots? Default: 1024', default=1024)
    parser.add_argument('--backend', type=str, help='Which backend to use? Default: FakeVigoV2', default="5q-qvm")
    parser.add_argument('--setqubits', type=int, nargs='+', help='Which qubits to use?: Default: 0123 and transpile')

    #  Parse die Argumente
    args = parser.parse_args()
    # %%
    from pyquil import get_qc, Program
    from pyquil.gates import H, CNOT, Z, MEASURE, S, X, Y, I, RX, RZ, FENCE
    from pyquil.quilbase import Declare
    import pyquil.paulis as Pauli
    from pyquil.paulis import PauliTerm
    import networkx as nx
    from pyquil.quantum_processor import NxQuantumProcessor
    from pyquil.noise import decoherence_noise_with_asymmetric_ro
    from matplotlib import pyplot as plt
    import os
    import sys
    import numpy as np
    import json
    import time
    tim = time.localtime()
    print("%s.%s. %s:%s" % (tim.tm_mday, tim.tm_mon, tim.tm_hour, tim.tm_min))

    folder = os.getcwd()
    while not folder.endswith("pyquil_program"):
        folder = os.path.dirname(folder)
    sys.path.append(os.path.join(folder, "pauli_lindblad_per"))

    from tomography.experiment import SparsePauliTomographyExperiment as tomography
    from primitives.pauli import PyQuilPauli

    plt.style.use("ggplot")

    # %%
    # Zugriff auf die Variablen
    # Here the backend for the simulation is prepared
    backend = get_qc(args.backend) #str(n)+'q-qvm'
    isa = backend.to_compiler_isa()
    backend_qubits = sorted(int(k) for k in isa.qubits.keys())
    # By default every qubit on our fake backend is connected with every other qubit. This breaks the algorhythm. So we need to change the topology
    # So this line here changes the topology of the fake backend. The condition in the end determines it
    edges = [(q1, q2) for q1 in backend_qubits for q2 in backend_qubits if abs(q1-q2) == 1]
    # Build the NX graph
    topo = nx.from_edgelist(edges)
    # You would uncomment the next line if you have disconnected qubits
    # topo.add_nodes_from(qubits)
    quantum_processor = NxQuantumProcessor(topo)
    quantum_processor.noise_model = decoherence_noise_with_asymmetric_ro(quantum_processor.to_compiler_isa())  # Optional
    backend.compiler.quantum_processor = quantum_processor
    # Here the backend without noise is prepared for later comparison
    perfect_quantum_processor = NxQuantumProcessor(topo)
    perfect_backend = get_qc(args.backend)
    perfect_backend.compiler.quantum_processor = perfect_quantum_processor

    qubits = [0,1,2,3] #[9,10,11,12] for MelbourneV2
    num_qubits = 4
    if args.setqubits != None:
        if len(args.setqubits) != 4:
            raise Exception("Must be 4 qubits when given")
        qubits = args.setqubits
        num_qubits = len(backend.quantum_processor.qubit_topology().nodes)
        
    tomography_connections = args.plusone
    sum_over_lambda = args.sum
    if sum_over_lambda:
        tomography_connections = True

    pntsamples = args.pntsamples
    pntsinglesamples = args.pntsinglesamples
    persamples = args.persamples
    shots = args.shots

    namebase = "" 
    print("Arguments where set as:")
    for arg_name, arg_value in vars(args).items():
        if arg_name == "setqubits" and arg_value == None:
            arg_value = "[0,1,2,3]_and_transpile"
        print("\t%s: %s" % (arg_name, str(arg_value)))
        namebase += str(arg_value) + "_"
    print("Namebase will be: " + namebase)
    #%%
    def trotterLayer(h,J,dt,n):
        trotterLayer = Program(Declare("ro", "BIT", num_qubits))
        for q in qubits:
            trotterLayer += RX(dt*4*h, q)
        for i, j in [(qubits[2*i], qubits[2*i+1]) for i in range(n)]:
            trotterLayer += CNOT(i,j)
        for q in [qubits[2*i+1] for i in range(n)]:
            trotterLayer += RZ(-4*J*dt, q)
        for i, j in [(qubits[2*i], qubits[2*i+1]) for i in range(n)]:
            trotterLayer += CNOT(i,j)
        for i, j in [(qubits[2*i+1], qubits[2*i+2]) for i in range(n-1)]:
            trotterLayer += CNOT(i,j)
        for q in [qubits[2*i+1] for i in range(n-1)]:
            trotterLayer += RZ(-4*J*dt, q)
        for i, j in [(qubits[2*i+1], qubits[2*i+2]) for i in range(n-1)]:
            trotterLayer += CNOT(i,j)
        return trotterLayer

    def maketrotterCircuit(s):
        tL = trotterLayer(h=1, J=-0.15, dt=0.2, n=2)
        trotterCircuit = Program(Declare("ro", "BIT", num_qubits))
        for i in range(s):
            trotterCircuit += tL
            #trotterCircuit += FENCE() #the fence is not currently supported by quilc
        return trotterCircuit

    circuits = [maketrotterCircuit(i) for i in range(1,15)]
    #print(circuits[1])

    qubits = set()
    for circuit in circuits:
        for q in backend.compile(circuit).get_qubit_indices():
            qubits.add(q)
        break #temporary
    print("Qubits set to ", qubits)
    

    # %%
    def executor(circuits):
        def take_counts(array):
            dic = dict()
            for a in array:
                if tuple(k for k in a) in dic:
                    dic[tuple(k for k in a)] += 1
                else:
                    dic[tuple(k for k in a)] = 1
            return dic
        counts = []
        for circuit in circuits:
            result = backend.run(backend.compile(circuit.wrap_in_numshots_loop(shots=shots))).get_register_map()['ro']
            counts.append(take_counts(result))
        return counts

    # %%
    print("initialize experiment")
    experiment = tomography(circuits = circuits, inst_map = backend_qubits, backend = backend, tomography_connections=tomography_connections, sum_over_lambda=sum_over_lambda)

    print("generate circuits")
    experiment.generate(samples = 1, single_samples = 1, depths = [2,4,8,16])
    tim = time.localtime()
    print("%s.%s. %s:%s" % (tim.tm_mday, tim.tm_mon, tim.tm_hour, tim.tm_min))
    # %%
    print("run experiment")
    experiment.run(executor)
    # %%
    print("analyse experiment")
    noisedataframe = experiment.analyze()

    # %%
    perexp = experiment.create_per_experiment(circuits)

    # %%
    noise_strengths = [0,0.5,1,2]
    expectations = []
    for q in qubits:
        expect = "I"*(len(backend_qubits)) #15
        expect = expect[:q] + 'Z' + expect[q+1:]
        expectations.append("".join(expect))
    print("do PER runs")
    tim = time.localtime()
    print("%s.%s. %s:%s" % (tim.tm_mday, tim.tm_mon, tim.tm_hour, tim.tm_min))
    perexp.generate(expectations = expectations, samples = 4, noise_strengths = noise_strengths)

    # %%
    print("Run PER")
    perexp.run(executor)

    # %%
    circuit_results = perexp.analyze()













from pyquil.api import local_forest_runtime
if __name__ == "__main__":
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
