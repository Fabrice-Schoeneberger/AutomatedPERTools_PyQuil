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
    return (twoqubit_error_template, None, None)

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

def executor(circuits, backend, shots, noise_model=None):
    #import multiprocessing, os, pickle, uuid
    #manager = multiprocessing.Manager()
    #counts = manager.list()
    #lock = multiprocessing.Lock()
    counts = []
    for circuit in circuits:
        if noise_model != None:
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

def find_used_qubits(circuits, backend):
    qubits = set()
    for circuit in circuits:
        for q in backend.compile(circuit).get_qubit_indices():
            qubits.add(q)
    return qubits

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

    # PNT Arguments
    parser.add_argument('--pntsamples', type=int, help='How many samples in PNT? Default: 16', default=16)
    parser.add_argument('--pntsinglesamples', type=int, help='How many single samples in PNT? Default: 100', default=100)
    parser.add_argument('--depths', type=int, nargs='+', help='Decide the depths of the pnt-samples. Default: [2,4,8,16]')
    parser.add_argument('--onlyTomography', help='Only does the tomography and then ends the program', default=False, action='store_true')
    # PER Arguments
    parser.add_argument('--persamples', type=int, help='How many samples in PER? Default: 1000', default=1000)
    parser.add_argument('--noise_strengths', type=float, nargs='+', help='Decide the noise strengths at which PER will run. Default: [0.5,1,2]')
    parser.add_argument('--expectations', type=str, nargs='+', help='Decide the expectation values which whill be measured at PER. Default: Z on all used qubits')
    # General Arguments
    parser.add_argument('--shots', type=int, help='How many shots per circuit? Default: 1024', default=1024)
    parser.add_argument('--num_qubits', type=int, help='Define how many qubits the simulated backend should have. Layout: Line? Default: How ever many the biggest circuit has', default=5)
    parser.add_argument('--circuitfilename', type=str, help='Set the name of the file, that contains the function making you circuit. Default: circuit.py', default="circuits")
    parser.add_argument('--circuitfunction', type=str, help='Set the name of the function, that return your circuits. The function can only take a backend as input and has to return an array of circuits. Default: make_initial_Circuit', default="make_initial_Circuit")
    parser.add_argument('--cross', '-c', help='Simulates Cross Talk Noise', default=False, action='store_true')
    parser.add_argument('--make_plots', help='If true, generates standard PER and vZNE plots', default=False, action='store_true')
    parser.add_argument('--setqubits', type=int, nargs='+', help='Which qubits on the backend should be used? When set tries to fullfill request, but does not garantie. Default: 0123 and transpile')
    parser.add_argument('--foldername_extra', type=str, help='Attach something to the end of the foldernamebase', default="")
    parser.add_argument('--real_backend', help='Toggels the use of a real IBM Backend. You will be prompted an extra input to confirm this. Use QiskitRuntimeService.save_account(channel="ibm_quantum", token="<your_token>") first', default=False, action='store_true')

    #  Parse the Arguments
    args = parser.parse_args()

    # %% imports and system path appending
    import pickle
    from matplotlib import pyplot as plt
    import os
    import json
    import uuid
    import numpy as np
    

    import os
    import sys
    i = 0
    folder = os.path.dirname(os.path.abspath(__file__))
    while not "pauli_lindblad_per" in os.listdir(folder):
        folder = os.path.dirname(folder)
        i+=1
        if i == 50:
            raise Exception("pauli_lindblad_per not found. Please make sure it is in this or a parent folder")
    sys.path.append(os.path.join(folder, "pauli_lindblad_per"))


    from tomography.experiment import SparsePauliTomographyExperiment as tomography

    plt.style.use("ggplot")
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
    noise_strengths = [0.5,1,2]
    if args.noise_strengths != None:
        noise_strengths = args.noise_strengths
    do_cross_talk_noise = args.cross
    onlyTomography = args.onlyTomography

    pntsamples = args.pntsamples
    pntsinglesamples = args.pntsinglesamples
    persamples = args.persamples
    make_plots = args.make_plots
    shots = args.shots
    if args.expectations != None:
        expectations = args.expectations
    else: 
        qubits = [0,1,2,3]
        expectations = []
        for q in qubits:
            expect = "I"*4
            expect = expect[:q] + 'Z' + expect[q+1:]
            expectations.append("".join(reversed(expect)))
    (noise_model, twoqubit_error_template, singlequbit_error_template) = get_noise_model()

    # %% Make the initial Circuits
    print("")
    print("Make Circuits")
    n = 2

    if args.circuitfilename.endswith(".py"):
        circuitfilename = args.circuitfilename[:-3]
    else:
        circuitfilename = args.circuitfilename
    import importlib.util
    spec = importlib.util.spec_from_file_location(circuitfilename, f"{circuitfilename}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    get_circuit = getattr(module, args.circuitfunction)
    circuits = get_circuit(backend)
    qubits = find_used_qubits(circuits, backend)
    max_qubits = 0
    for circuit in circuits: 
        max_qubits = max(max_qubits, max(qubits))
    if args.num_qubits == 0:
        args.num_qubits = max_qubits
        backend = get_backend(args)
        circuits = get_circuit(backend)
    elif max_qubits > get_backend(args, return_backend_qubits=True):
        raise Exception("Backend has to few qubits for a circuit. " +str(get_backend(args, return_backend_qubits=True))+"/"+str(max_qubits)+ " given. Give more qubits with --num_qubits x")
    
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
        if str(arg_value) != "":
            namebase += str(arg_value) + "_"
    print("Backend name is: "+ str(backend.name))
    namebase += "_" + str(backend.name)
    namebase = namebase[:-1]
    #os.makedirs(namebase, exist_ok=True)
    print("Namebase will be: " + namebase)
    namebase += "/"

    # %% initialize experiment
    print_time("initialize experiment")
    experiment = tomography(circuits = circuits, inst_map = [i for i in range(get_backend(args, return_backend_qubits=True))], backend = backend)
    # %% generate PNT circuits
    print_time("generate circuits")
    experiment.generate(samples = pntsamples, single_samples = pntsinglesamples, depths = depths)

    # %% run PNT experiment
    print_time("run experiment")
    experiment.run(executor, shots, do_cross_talk=do_cross_talk_noise, apply_cross_talk=apply_cross_talk_proxy, noise_model=noise_model)

    # %% analyse PNT experiment.
    print_time("analyse experiment")
    noisedataframe = experiment.analyze()
    # %% Save all the data. End Tomography Only
    print_time("Saving data")
    id = str(uuid.uuid4())
    coeffs_dict_list = []
    infidelities_list = []
    cliff_layer_list = []
    for layer in experiment.analysis.get_all_layer_data():
        coeffs_dict = dict(layer.noisemodel.coeffs)
        infidelities = {term: 1-layer._term_data[term].fidelity for term in layer._term_data}
        coeffs_dict_list.append(coeffs_dict)
        infidelities_list.append(infidelities)
        cliff_layer_list.append(layer.layer._cliff_layer)
    os.makedirs("automatedPERrun_collection/"+ namebase, exist_ok=True)
    os.makedirs("automatedPERrun_collection/"+ namebase + "coeffs/", exist_ok=True)
    os.makedirs("automatedPERrun_collection/"+ namebase + "infidelities/", exist_ok=True)
    os.makedirs("automatedPERrun_collection/"+ namebase + "cliff_layer/", exist_ok=True)
    with open("automatedPERrun_collection/" + namebase + "coeffs/" + str(id) + ".pickle", "wb") as f:
        pickle.dump(coeffs_dict_list, f)
    with open("automatedPERrun_collection/" + namebase + "infidelities/" + str(id) + ".pickle", "wb") as f:
        pickle.dump(infidelities_list, f)
    with open("automatedPERrun_collection/" + namebase + "cliff_layer/" + str(id) + ".pickle", "wb") as f:
        pickle.dump(cliff_layer_list, f)
    with open("automatedPERrun_collection/" + namebase + "noise_model.pickle", "wb") as f:
        pickle.dump((noise_model, twoqubit_error_template, singlequbit_error_template), f)
    with open("automatedPERrun_collection/" + namebase + "circuits.pickle", "wb") as f:
        pickle.dump(circuits, f)
        
    if onlyTomography:
        print_time("Tomography Ended")
        print("")
        return
    # %% create per experiment
    print_time("Create PER Experiment")
    perexp = experiment.create_per_experiment(circuits)

    # %% generate per experiment and noise strengths
    print_time("do PER runs")
    perexp.generate(expectations = expectations, samples_per_circuit = persamples, noise_strengths = noise_strengths)
    # %% Run PER
    print_time("Run PER")
    perexp.run(executor, shots, do_cross_talk=do_cross_talk_noise, apply_cross_talk=apply_cross_talk_proxy, noise_model=noise_model)

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
    noisyresult = calculate_with_simple_backend(circuits, shots, persamples, backend, qubits, n, noise_model, apply_cross_talk=do_cross_talk_noise)
    res = calculate_with_simple_backend(circuits, shots, persamples, get_backend(args, return_perfect=True), qubits, n, noise_model)

    savi["noisyresult"] = noisyresult
    savi["res"] = res

    # %% Plot PER
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#E6194B', '#3CB44B', '#0082C8', '#F58231', '#911EB4', '#FFD700', '#46F0F0', '#F032E6', '#A9A9A9']
    fontsize =20
    markers = ['D', 'X', 'p', '*', 'D', 'X', 'p', '*', 'D', 'X', 'p', '*']

    with open(namebase + 'PER_data/'+id+'.json', 'w') as file:
        json.dump(savi, file)
    os.makedirs(os.path.join("automatedPERrun_collection", str(persamples)+"_"+str(noise_strengths)), exist_ok=True)
    os.makedirs(os.path.join("automatedPERrun_collection", str(persamples)+"_"+str(noise_strengths), "rel_value"), exist_ok=True)
    os.makedirs(os.path.join("automatedPERrun_collection", str(persamples)+"_"+str(noise_strengths), "rel_error"), exist_ok=True)
    def make_PER_plots():
        for j in range(len(noise_strengths)+2): # 0 vZNE results, 1-x results at noise level
            plt.figure(figsize=(20,12))
            for i, noise in enumerate(noise_strengths):
                if j != 0 and j!= 1 and noise == noise_strengths[j-2]:
                    plt.errorbar(range(1,15), [res[i] for res in results_at_noise], fmt=markers[i]+'--', yerr=[res[i] for res in results_at_noise_errors], capsize=12*fontsize/20, label='N='+str(noise), color= colors[i], zorder=1, markersize=10*fontsize/20, linewidth=2*fontsize/20)
                else:
                    plt.plot(range(1,15), [res[i] for res in results_at_noise], markers[i], label='N='+str(noise), color= colors[i], zorder=2, markersize=10*fontsize/20)
                
            plt.plot(range(1,15), res, 'o:', label="Trotter Simulation", color= colors[len(noise_strengths)], zorder=2, markersize=10*fontsize/20, linewidth=2*fontsize/20)
            plt.plot(range(1,15), noisyresult, 'o', label="Unmitigated", color= colors[len(noise_strengths)+1], zorder=2, markersize=10*fontsize/20)
            if j == 1:
                plt.errorbar(range(1,15), results, yerr=[[np.abs(res) for res in results_errors]],fmt='s--', capsize=12*fontsize/20, label="PER", color= colors[len(noise_strengths)+2], zorder=1, markersize=10*fontsize/20, linewidth=2*fontsize/20)
            else:
                plt.plot(range(1,15), results, 's', label="PER", color= colors[len(noise_strengths)+2], zorder=2, markersize=10*fontsize/20)

            plt.ylim([-2.1,2.1])
            plt.legend(fontsize=fontsize, loc='upper right')
            plt.title("Trotter Simulation with PER", fontsize=fontsize*1.2)
            plt.xlabel("Trotter Steps", fontsize=fontsize)
            plt.ylabel("Z Magnetization", fontsize=fontsize)
            plt.tick_params(axis='both', which='major', labelsize=fontsize)
            plt.savefig("automatedPERrun_collection/"+str(persamples)+"_"+str(noise_strengths)+"/Trotter_Sim_PER"+str(j)+".png")
    def make_vZNE_fit_plots(i):
        fontsize = 30
        def plot(perdat):
            """Plots the expectation values against an exponential fit.
            """
            fig, ax = plt.subplots(figsize=(14, 8))
            ax.errorbar(list(sorted(perdat["data"].keys())), [perdat["data"][s]/perdat["counts"][s] for s in list(sorted(perdat["data"].keys()))], yerr=[np.std(perdat["dataStatistic"][strength]) for strength in list(sorted(perdat["data"].keys()))],  linestyle = "None", marker = "o", color = "tab:blue", capsize=12*fontsize/20, markersize=12*fontsize/20)
            a = perdat["expectation"]
            b = perdat["b"]
            xlin = np.linspace(0, max(list(sorted(perdat["data"].keys()))), 100)
            ax.plot(xlin, [a*np.exp(b*x) for x in xlin], color = "tab:blue", linestyle="--", linewidth=2*fontsize/20)
            ax.tick_params(axis='both', which='major', labelsize=fontsize)
            
        for ex in expectations:
            # Change the number in circuit_results[NUMBER] to show different circuits fits
            ax = plot(circuit_results[i].get(ex, None))
            # If the error becomes to big, uncomment the next line to still see anything but a straight line
            #plt.ylim(-1.5,1.5)
            plt.title('Expectation vs Noise Strength '+str(ex), fontsize=fontsize*1.2)
            plt.xlabel("Noise Strength", fontsize=fontsize)
            plt.ylabel("Expectation", fontsize=fontsize)
            plt.savefig("automatedPERrun_collection/"+str(persamples)+"_"+str(noise_strengths)+"/Expectation_vs_Noise_Strength_"+str(ex)+".png")

    if make_plots:
        make_PER_plots()
        make_vZNE_fit_plots(0)
    rel_value = [abs(res[i]-results[i])/abs(res[i]-noisyresult[i]) for i, _ in enumerate(res)]
    rel_error = [results_errors[i]/abs(res[i]-noisyresult[i]) for i, _ in enumerate(res)]
    import pickle
    with open(os.path.join("automatedPERrun_collection", namebase, "rel_value/") + str(id) + ".pickle", "wb") as f:
        pickle.dump(rel_value, f)
    with open(os.path.join("automatedPERrun_collection", namebase, "rel_error/") + str(id) + ".pickle", "wb") as f:
        pickle.dump(rel_error, f)


    with open(namebase + "circuit_results.pickle", "wb") as f:
        pickle.dump(circuit_results, f)

    print_time()
    print("")


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
