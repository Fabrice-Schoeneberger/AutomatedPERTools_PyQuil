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

def get_backend(args, return_perfect=False, return_backend_qubits=False):
    # Here the backend for the simulation is prepared
    from pyquil.quantum_processor import NxQuantumProcessor
    from pyquil.noise import decoherence_noise_with_asymmetric_ro
    from pyquil import get_qc
    import networkx as nx
    backend = get_qc("5q-qvm") #str(n)+'q-qvm' #args.backend
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
    #quantum_processor.noise_model = decoherence_noise_with_asymmetric_ro(quantum_processor.to_compiler_isa())
    quantum_processor.noise_model = get_noise_model()[0] 
    backend.compiler.quantum_processor = quantum_processor
    if not return_perfect:
        return backend
    # Here the backend without noise is prepared for later comparison
    perfect_quantum_processor = NxQuantumProcessor(topo)
    perfect_backend = get_qc(args.backend)
    perfect_backend.compiler.quantum_processor = perfect_quantum_processor
    if return_perfect:
        return perfect_backend

def get_noise_model():
    from pyquil.paulis import PauliTerm
    from pyquil.noise import NoiseModel
    # Define Pauli operations using PauliTerm
    twoqubit_errorops = [
        PauliTerm('Z', 0) * PauliTerm('Y', 1),  # Pauli('YZ')
        PauliTerm('Y', 0) * PauliTerm('I', 1),  # Pauli('IY')
        PauliTerm('Y', 0) * PauliTerm('Y', 1),  # Pauli('YY')
        PauliTerm('Y', 0) * PauliTerm('X', 1),  # Pauli('XY')
        PauliTerm('I', 0) * PauliTerm('I', 1)   # Pauli('XY')
    ]

    # Corresponding probabilities
    twoqubit_errorprobs = [0.008802700270751796, 0.0032989083407153896, 0.01917444731546973, 0.019520575974201874]
    twoqubit_errorprobs.append(1-sum(twoqubit_errorprobs))

    twoqubit_error_template = [(op, p) for op,p in zip(twoqubit_errorops, twoqubit_errorprobs)]
    singlequbit_error_template = []
    # Create a NoiseModel object
    return (NoiseModel(twoqubit_errorops, twoqubit_errorprobs), twoqubit_error_template, singlequbit_error_template)

def executor(circuits, backend, shots):
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
    #parser.add_argument('--plusone', '-p', help='Takes Neighboring qubits into account', default=False, action='store_true')
    parser.add_argument('--sum', '-s', help='Same as -p and turns sumation on over neighboring qubits', default=False, action='store_true')
    parser.add_argument('--pntsamples', type=int, help='How many samples in PNT? Default: 16', default=16)
    parser.add_argument('--pntsinglesamples', type=int, help='How many single samples in PNT? Default: 100', default=100)
    parser.add_argument('--persamples', type=int, help='How many samples in PER? Default: 100', default=100)
    parser.add_argument('--shots', type=int, help='How many shots? Default: 1024', default=1024)
    parser.add_argument('--backend', type=str, help='Which backend to use? Default: Line5', default="Line5")
    parser.add_argument('--cross', '-c', help='Simulates Cross Talk Noise', default=False, action='store_true')
    #parser.add_argument('--allqubits', '-a', help='runs over all qubits in the tomography', default=False, action='store_true')
    parser.add_argument('--onlyTomography', help='Only does the tomography and then ends the program', default=False, action='store_true')
    parser.add_argument('--setqubits', type=int, nargs='+', help='Which qubits to use?: Default: 0123 and transpile')
    parser.add_argument('--depths', type=int, nargs='+', help='Decide the depths of the pnt-samples. Default: [2,4,8,16]')

    #  Parse die Argumente
    args = parser.parse_args()

    # %% imports and system path appending
    import pickle
    from matplotlib import pyplot as plt
    import os
    import json
    print_time("Starting")

    import os
    import sys
    i = 0
    parentfolder = "pyquil_program"
    folder = os.getcwd()
    while not folder.endswith(parentfolder):
        folder = os.path.dirname(folder)
        i+=1
        if i == 50:
            raise Exception("Parent Folder not found. Is "+ parentfolder + " the correct name?")
    sys.path.append(os.path.join(folder, "pauli_lindblad_per"))


    from tomography.experiment import SparsePauliTomographyExperiment as tomography

    plt.style.use("ggplot")
    # %% Decipher the Arguments

    backend = get_backend(args)

    qubits = [0,1,2,3] #[9,10,11,12] for MelbourneV2
    num_qubits = 4
    if args.setqubits != None:
        if len(args.setqubits) != 4:
            raise Exception("Must be 4 qubits when given")
        qubits = args.setqubits
        num_qubits = get_backend(args, return_backend_qubits=True)#backend.num_qubits
    
    depths = [2,4,8,16]
    if args.depths != None:
        depths = args.depths
    do_cross_talk_noise = args.cross
    sum_over_lambda = args.sum
    onlyTomography = args.onlyTomography

    pntsamples = args.pntsamples
    pntsinglesamples = args.pntsinglesamples
    persamples = args.persamples
    shots = args.shots

    # %% Make the initial Circuits
    print("make trotter")
    n = 2
    circuits = make_initial_Circuit(qubits, num_qubits, backend, n)
    qubits = set()
    for circuit in circuits:
        for q in backend.compile(circuit).get_qubit_indices():
            qubits.add(q)
        break #temporary
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
    experiment = tomography(circuits = circuits, inst_map = [i for i in range(get_backend(args, return_backend_qubits=True))], backend = backend, sum_over_lambda=sum_over_lambda)

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
    with open(namebase + "experiment.pickle", "wb") as f:
        processor = experiment._procspec._processor
        experiment._procspec._processor = None
        pickle.dump(experiment, f)
        experiment._procspec._processor = processor

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
    (noise_model, twoqubit_error_template, singlequbit_error_template) = get_noise_model()
    with open("server_run_collection/" + namebase + "noise_model.pickle", "wb") as f:
       pickle.dump((noise_model, twoqubit_error_template, singlequbit_error_template), f)
        
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

    # %% ??
    circuit_results[-1]._per_circ.overhead(0)

    # %% Calculate unmitigated error and without error
    noisyresult = calculate_with_simple_backend(circuits, shots, persamples, backend, qubits, n, apply_cross_talk=do_cross_talk_noise)
    res = calculate_with_simple_backend(circuits, shots, persamples, get_backend(args, return_perfect=True), qubits, n)

    savi["noisyresult"] = noisyresult
    savi["res"] = res
    with open(namebase + '_arrays.json', 'w') as file:
        json.dump(savi, file)

    # %% Plot PER
    plt.figure(figsize=(10,6))
    for i, noise in enumerate(noise_strengths):
        plt.plot(range(1,15), [res[i] for res in results_at_noise], 'x', label='N='+str(noise))
        
    plt.plot(range(1,15), res, 'o:', label="Trotter Simulation")
    plt.plot(range(1,15), noisyresult, 'o', label="Unmitigated")
    plt.errorbar(range(1,15), results, yerr=results_errors, fmt='x', capsize=5, label="PER", color='#FFC0CB')

    plt.ylim([-1.8,1.8])
    plt.legend()
    plt.title("Trotter Simulation with PER")
    plt.xlabel("Trotter Steps")
    plt.ylabel("Z Magnetization")
    plt.savefig(namebase+"Trotter_Sim_PER.png")
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
