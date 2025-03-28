import logging
logger = logging.getLogger("experiment")

import time
tim = time.time()
laststring = ""
last_time_string = ""
def print_time(printstring=""):
    global tim, laststring, last_time_string
    import time
    local_tim = time.localtime()
    new_tim = time.time()
    time_difference = new_tim-tim
    tim = new_tim
    days = int(time_difference // (24 * 3600))
    time_difference %= (24 * 3600)
    hours = int(time_difference // 3600)
    time_difference %= 3600
    minutes = int(time_difference // 60)
    seconds = int(time_difference % 60)
    st = f" Time taken: {days:02}:{hours:02}:{minutes:02}:{seconds:02}"
    if not (laststring == "" and last_time_string == ""):
        print(laststring+last_time_string+st)    
    while len(printstring) < 24:
        printstring += " "
    printstring+="\t"
    laststring = printstring
    local_time_string = "%02d.%02d. %02d:%02d:%02d" % (local_tim.tm_mday, local_tim.tm_mon, local_tim.tm_hour, local_tim.tm_min, local_tim.tm_sec)
    last_time_string = local_time_string
    print(printstring+local_time_string+" Time taken: --:--:--:--", end="\r")

def make_initial_Circuit(qubits, num_qubits, backend, n):
    from pyquil import Program
    from pyquil.gates import H, CNOT, Z, MEASURE, S, X, Y, I, RX, RZ, FENCE
    from pyquil.quilbase import Declare
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
        tL = trotterLayer(h=1, J=-0.15, dt=0.2, n=n)
        trotterCircuit = Program(Declare("ro", "BIT", num_qubits))
        for i in range(s):
            trotterCircuit += tL
            #trotterCircuit += FENCE() #the fence is not currently supported by quilc
        return trotterCircuit

    return [maketrotterCircuit(i) for i in range(1,15)]

def make_initial_Circuit2(backend):
    from pyquil import Program
    from pyquil.gates import H, CNOT, Z, MEASURE, S, X, Y, I, RX, RZ, FENCE
    from pyquil.quilbase import Declare
    prog = Program()
    prog += Declare("ro", "BIT", 2)
    prog += CNOT(0,1)
    #prog += CNOT(0,1)
    return [prog]

def get_backend(args, return_perfect=False, return_backend_qubits=False):
    # Here the backend for the simulation is prepared
    from pyquil.quantum_processor import NxQuantumProcessor
    from pyquil.noise import decoherence_noise_with_asymmetric_ro
    from pyquil import get_qc
    import networkx as nx
    backend = get_qc(str(args.num_qubits) + "q-qvm") #str(n)+'q-qvm' #args.backend
    isa = backend.to_compiler_isa()
    backend_qubits = sorted(int(k) for k in isa.qubits.keys())
    if return_backend_qubits:
        return len(backend_qubits)
    # By default every qubit on our fake backend is connected with every other qubit. This breaks the algorhythm. So we need to change the topology
    # So this line here changes the topology of the fake backend. The condition in the end determines it
    edges = [(q1, q2) for q1 in backend_qubits for q2 in backend_qubits if abs(q1-q2) == 1]
    # Build the NX graph
    topo = nx.from_edgelist(edges)
    # You would uncomment the next line if you have disconnected qubits
    # topo.add_nodes_from(qubits)
    quantum_processor = NxQuantumProcessor(topo)
    backend.compiler.quantum_processor = quantum_processor
    return backend
    
def get_noise_model():
    #twoqubit_errorops = ['YZ', 'IY', 'YY', 'XY']
    #twoqubit_errorprobs = [0.008802700270751796, 0.0032989083407153896, 0.01917444731546973, 0.019520575974201874]
    twoqubit_errorops = ['IX']
    twoqubit_errorprobs = [0.05]
    twoqubit_error_template = [(op, p) for op,p in zip(twoqubit_errorops, twoqubit_errorprobs)]+[("II", 1-sum(twoqubit_errorprobs))]
    return (twoqubit_error_template, [])

def apply_noise_model(prog, backend, error_template, gate):
    import numpy as np
    # Defining the Pauli matrices
    Paulis = {}
    Paulis['I'] = np.array([[1.+0.j, 0.+0.j],
                [0.+0.j, 1.+0.j]])
    Paulis['X'] = np.array([[0.+0.j, 1.+0.j],
                [1.+0.j, 0.+0.j]])
    Paulis['Y'] = np.array([[0.+1.j, 0.+0.j],
                [0.+0.j, 0.-1.j]])
    Paulis['Z'] = np.array([[1.+0.j, 0.+0.j],
                [0.+0.j, -1.+0.j]])
    
    # Define the gates, that can be called
    Gates = {}
    Gates['CNOT'] = np.array([[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 1],
                    [0, 0, 1, 0]])
    
    # Draft the errornes gates
    noise_model = []
    for op, p in error_template:
        errorgate = np.array([[1]])
        for char in op:
            errorgate = np.kron(errorgate, Paulis[char])
        noise_model.append(np.sqrt(p)*np.matmul(errorgate, Gates[gate]))

    # Apply the gates to the circuit
    if np.log2(np.linalg.matrix_rank(errorgate)) == 2:
        for i, j in backend.qubit_topology().edges:
            prog.define_noisy_gate(gate, [i,j], noise_model)
            prog.define_noisy_gate(gate, [j,i], noise_model)
    elif np.log2(np.linalg.matrix_rank(errorgate)) == 1:
        for i in backend.qubits():
            prog.define_noisy_gate(gate, [i], noise_model)

def take_counts(array, backend):
    dic = dict()
    for a in array:
        a = "".join(list(str(k) for k in a))
        while len(a) < len(backend.qubits()):
            a += "0"
        if a in dic:
            dic[a] += 1
        else:
            dic[a] = 1
    return dic

def executor(circuits, backend, shots, apply_noise=True, noise_model=None):
    #import multiprocessing, os, pickle, uuid
    #manager = multiprocessing.Manager()
    #counts = manager.list()
    #lock = multiprocessing.Lock()
    counts = []
    for circuit in circuits:
        if apply_noise:
            apply_noise_model(circuit, backend, get_noise_model()[0], "CNOT")
        #process = multiprocessing.Process(target=executer_final, args=(circuit, backend, shots, counts, lock))
        #process.start()
        #while len(multiprocessing.active_children()) > multiprocessing.cpu_count():
            #print("test?",end='\r')
        #    pass
        result = backend.run((circuit.wrap_in_numshots_loop(shots=shots))).get_register_map()['ro']
        counts.append(take_counts(result, backend))
    #while len(multiprocessing.active_children()) > 1:
        #print("Apply Cross Talk Noise %s/%s" % (i+1, len_circuits), "; Number of active Threads: %03s" % len(multiprocessing.active_children()), "; List length: %s" % len(new_circuits), end='\r')
    #    pass
    #manager = None
    #counts = list(counts)
    for i, r in enumerate(counts):
        if not 1024 in r.values():
            circuit = circuits[i]
            liste = [str(inst) for inst in circuit if not "PRAGMA" in str(inst) and not "DEFGATE" in str(inst) and not "DECLARE" in str(inst)]
            #print(liste)
            #print(r)
    return counts

def executer_final(circuit, backend, shots, counts, lock):
    result = backend.run((circuit.wrap_in_numshots_loop(shots=shots))).get_register_map()['ro']
    with lock:
        counts.append(take_counts(result, backend))

def calculate_with_simple_backend(circuits, shots, persamples, backend, qubits, n, apply_cross_talk=False):
    res = []
    if apply_cross_talk:
        circuits = apply_cross_talk_proxy(circuits, backend)
    for circ in circuits:
        qc = circ.copy()
        qc.measure_all()
        count = executor(qc, backend, shots*persamples)
        count = {tuple(int(k) for k in key):count[key] for key in count.keys()}
        tot = 0
        for key in count.keys():
            num = sum([(-1)**bit for i, bit in enumerate(key) if len(key)-1-i in qubits])
            tot += num*count[key]
        res.append(tot/(shots*persamples*n*2))
    return res

# %% Cross Talk Noise Functions
def apply_cross_talk_proxy(circuits, backend):
    return circuits #TODO: Implement this

def circuit_to_layers(qc):
    layers = []
    inst_list = [inst for inst in qc if not inst.ismeas()] 

    #pop off instructions until inst_list is empty
    while inst_list:

        circ = qc.copy_empty() #blank circuit to add instructions
        layer_qubits = set() #qubits in the support of two-qubit clifford gates

        for inst in inst_list.copy(): #iterate through remaining instructions

            #check if current instruction overlaps with support of two-qubit gates
            #already on layer
            if not layer_qubits.intersection(inst.support()):
                circ.add_instruction(inst) #add instruction to circuit and pop from list
                inst_list.remove(inst)

            if inst.weight() >= 2:
                layer_qubits = layer_qubits.union(inst.support()) #add support to layer

        if circ: #append only if not empty
            layers.append(circ)

    return layers

# %% Start main(), Define the optional arguments
def main():
    import argparse
    parser = argparse.ArgumentParser()
        
    # Define an Argument
    parser.add_argument('--pntsamples', type=int, help='How many samples in PNT? Default: 16', default=16)
    parser.add_argument('--pntsinglesamples', type=int, help='How many single samples in PNT? Default: 100', default=100)
    parser.add_argument('--persamples', type=int, help='How many samples in PER? Default: 100', default=100)
    parser.add_argument('--shots', type=int, help='How many shots? Default: 1024', default=1024)
    parser.add_argument('--num_qubits', type=int, help='Define how many qubits the backend should have. Layout: Line? Default: 5', default=5)
    parser.add_argument('--cross', '-c', help='Simulates Cross Talk Noise', default=False, action='store_true')
    parser.add_argument('--onlyTomography', help='Only does the tomography and then ends the program', default=False, action='store_true')
    parser.add_argument('--setqubits', type=int, nargs='+', help='Which qubits to use?: Default: 0123 and transpile')
    parser.add_argument('--depths', type=int, nargs='+', help='Decide the depths of the pnt-samples. Default: [2,4,8,16]')
    parser.add_argument('--foldername_extra', type=str, help='Attach something to the end of the foldernamebase', default="")

    #  Parse die Argumente
    args = parser.parse_args()

    # %% imports and system path appending
    import pickle
    import os
    import json

    import os
    import sys
    i = 0
    #parentfolder = "pyquil_program"
    folder = os.getcwd()
    import sys
    print(sys.version)
    while not "pauli_lindblad_per" in os.listdir(folder):
        folder = os.path.dirname(folder)
        i+=1
        if i == 50:
            raise Exception("pauli_lindblad_per not found. Please make sure it is in this or a parent folder")
    sys.path.append(os.path.join(folder, "pauli_lindblad_per"))


    from tomography.experiment import SparsePauliTomographyExperiment as tomography

    # %% Decipher the Arguments

    backend = get_backend(args)

    qubits = [0,1,2,3] #[9,10,11,12] for MelbourneV2
    if args.setqubits != None:
        if len(args.setqubits) != 4:
            raise Exception("Must be 4 qubits when given")
        qubits = args.setqubits
    
    depths = [2,4,8,16]
    if args.depths != None:
        depths = args.depths
    do_cross_talk_noise = args.cross
    onlyTomography = args.onlyTomography

    pntsamples = args.pntsamples
    pntsinglesamples = args.pntsinglesamples
    persamples = args.persamples
    shots = args.shots

    # %% Make the initial Circuits
    print("")
    print("Make Circuits")
    n = 2
    
    circuits = make_initial_Circuit2(backend)
    qubits = set()
    for circuit in circuits:
        for q in backend.compile(circuit).get_qubit_indices():
            qubits.add(q)
    print("Qubits set to ", qubits)
    # %% Make namebase
    namebase = "" 
    print("Arguments where set as:")
    for arg_name, arg_value in vars(args).items():
        if arg_name == "setqubits" and arg_value == None:
            arg_value = str(qubits)
        if arg_name == "depths" and arg_value == None:
            arg_value = str(depths)
        print("\t%s: %s" % (arg_name, str(arg_value)))
        namebase += str(arg_value) + "_"
    namebase = namebase[:-1]
    os.makedirs(namebase, exist_ok=True)
    print("Namebase will be: " + namebase)
    namebase += "/"
    #circuits[0].draw()

    # %% initialize experiment
    print_time("initialize experiment")
    experiment = tomography(circuits = circuits, inst_map = [i for i in range(get_backend(args, return_backend_qubits=True))], backend = backend)

    import multiprocessing
    # %% generate PNT circuits
    print_time("generate circuits")
    experiment.generate(samples = pntsamples, single_samples = pntsinglesamples, depths = depths)

    # %% run PNT experiment
    print_time("run experiment")
    experiment.run(executor, shots, do_cross_talk=do_cross_talk_noise, apply_cross_talk=apply_cross_talk_proxy)

    # %% analyse PNT experiment. End Tomography Only
    print_time("analyse experiment")
    noisedataframe = experiment.analyze()
    # %% Save all the data. End Tomography Only
    print_time("Saving data")
    #with open(namebase + "experiment.pickle", "wb") as f:
    #    processor = experiment._procspec._processor
    #    experiment._procspec._processor = None
    #    pickle.dump(experiment, f)
    #    experiment._procspec._processor = processor

    coeffs_dict_list = []
    infidelities_list = []
    for layer in experiment.analysis.get_all_layer_data():
        coeffs_dict = dict(layer.noisemodel.coeffs)
        infidelities = {term: 1-layer._term_data[term].fidelity for term in layer._term_data}
        coeffs_dict_list.append(coeffs_dict)
        infidelities_list.append(infidelities)
    os.makedirs("server_run_collection/"+namebase, exist_ok=True)
    with open("server_run_collection/" + namebase + "coeffs.pickle", "wb") as f:
        pickle.dump(coeffs_dict_list, f)
    with open("server_run_collection/" + namebase + "infidelities.pickle", "wb") as f:
        pickle.dump(infidelities_list, f)
    (twoqubit_error_template, singlequbit_error_template) = get_noise_model()
    with open("server_run_collection/" + namebase + "noise_model.pickle", "wb") as f:
       pickle.dump((twoqubit_error_template, singlequbit_error_template), f)
    #with open("server_run_collection/" + namebase + "circuits.pickle", "wb") as f:
    #    pickle.dump(circuits, f)
    #process.join()
    if onlyTomography:
        print_time("Tomography Ended")
        print("")
        return
    # %% create per experiment
    print_time("Create PER Experiment")
    perexp = experiment.create_per_experiment(circuits)

    # %% generate per experiment and noise strengths
    noise_strengths = [0,0.5,1,2]
    expectations = []
    for q in qubits:
        expect = "I"*(get_backend(args, return_backend_qubits=True)) #15
        expect = expect[:q] + 'Z' + expect[q+1:]
        expectations.append("".join(reversed(expect)))
    print_time("do PER runs")
    perexp.generate(expectations = expectations, samples = persamples, noise_strengths = noise_strengths)
    # %% Run PER
    print_time("Run PER")
    if len(multiprocessing.active_children()) > 1:
        raise Exception("Too many children")
    perexp.run(executor, shots, do_cross_talk=do_cross_talk_noise, apply_cross_talk=apply_cross_talk_proxy)

    # %% Analyze PER and Delete Pickled PERrun Data
    print_time("Analyze PER")
    circuit_results = perexp.analyze()

    print_time("Delete Pickled PERrun Data")
    perexp.delete_pickles()

    # %% Extract data
    results_errors = []
    results_at_noise = []
    results_at_noise_errors = []
    results = []

    for run in circuit_results:
        tot = 0
        tot_error = 0
        tot_at_noise = [0 for _ in range(len(noise_strengths))]
        tot_at_noise_errors = [0 for _ in range(len(noise_strengths))]
        for op in expectations:
            #get the full per results
            expec = run.get_result(op).expectation
            tot += expec/len(expectations)

            #get the corresponding fit-errors
            expec_error = run.get_result(op).expectation_error
            tot_error += expec_error/len(expectations)

            #get the value at the different noise levels
            expec_at_noise = run.get_result(op).get_expectations()
            for i in range(0,len(tot_at_noise)):
                tot_at_noise[i] += expec_at_noise[i]/len(expectations)

            expec_at_noise_error = [run.get_result(op).get_std_of_strengths(strength) for strength in noise_strengths]
            for i in range(0,len(tot_at_noise)):
                tot_at_noise_errors[i] += expec_at_noise_error[i]/len(expectations)

            

        results.append(tot)
        results_errors.append(tot_error)
        results_at_noise.append(tot_at_noise)
        results_at_noise_errors.append(tot_at_noise_errors)

    savi = {}
    savi["results"] = results
    savi["results_errors"] = results_errors
    savi["results_at_noise"] = results_at_noise
    savi["results_at_noise_errors"] = results_at_noise_errors

    # %% Calculate unmitigated error and without error
    noisyresult = calculate_with_simple_backend(circuits, shots, persamples, backend, qubits, n, apply_cross_talk=do_cross_talk_noise)
    res = calculate_with_simple_backend(circuits, shots, persamples, get_backend(args, return_perfect=True), qubits, n)

    savi["noisyresult"] = noisyresult
    savi["res"] = res
    with open(namebase + '_arrays.json', 'w') as file:
        json.dump(savi, file)
    # %% Plot Expectation vs Noise Strength and save results

    with open(namebase + "circuit_results.pickle", "wb") as f:
        pickle.dump(circuit_results, f)

    print_time()
    print("")
    """ for i in range (len(expectations)):
        ax = circuit_results[0].get_result(expectations[i]).plot()
        plt.title("Expectation vs Noise Strength " + expectations[i])
        plt.xlabel("Noise Strength")
        plt.ylabel("Expectation")
        plt.savefig(namebase+"_Expectation_vs_Noise_Strength_" + expectations[i] + ".png") """

# %% Start Program
if __name__ == "__main__":
    print_time("Starting")
    from sys import platform
    if platform == "linux" or platform == "linux2":
        main()
        # linux
    elif platform == "darwin":
        raise Exception("IOS Not supported")
        # OS X
    elif platform == "win32":
        from pyquil.api import local_forest_runtime
        with local_forest_runtime():
            main()
        # Windows...
