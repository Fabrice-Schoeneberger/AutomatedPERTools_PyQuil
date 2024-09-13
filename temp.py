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
    from matplotlib import pyplot as plt
    import os
    import sys
    import time
    tim = time.localtime()
    print("%s.%s. %s:%s" % (tim.tm_mday, tim.tm_mon, tim.tm_hour, tim.tm_min))

    folder = os.getcwd()
    while not folder.endswith("pyquil_program"):
        folder = os.path.dirname(folder)
    sys.path.append(os.path.join(folder, "pauli_lindblad_per"))

    from tomography.experiment import SparsePauliTomographyExperiment as tomography

    plt.style.use("ggplot")

    # %%
    # Zugriff auf die Variablen
    backend = get_backend(args)

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
    circuits = make_initial_Circuit(qubits, num_qubits, backend)
    #print(circuits[1])

    qubits = set()
    for circuit in circuits:
        for q in backend.compile(circuit).get_qubit_indices():
            qubits.add(q)
        break #temporary
    print("Qubits set to ", qubits)
    


    # %%
    print("initialize experiment")
    experiment = tomography(circuits = circuits, inst_map = get_backend(args, return_backend_qubits=True), backend = backend, tomography_connections=tomography_connections, sum_over_lambda=sum_over_lambda)

    print("generate circuits")
    experiment.generate(samples = pntsamples, single_samples = pntsinglesamples, depths = [2,4,8,16])
    tim = time.localtime()
    print("%s.%s. %s:%s" % (tim.tm_mday, tim.tm_mon, tim.tm_hour, tim.tm_min))
    # %%
    print("run experiment")
    experiment.run(executor)
    # %%
    print("analyse experiment")
    noisedataframe = experiment.analyze()

    with open(namebase + "experiment.pickle", "wb") as f:
        pickle.dump(experiment, f)
    print("analyse experiment")
    print_time()
    noisedataframe = experiment.analyze()
    with open(namebase + "noisedataframe.pickle", "wb") as f:
        pickle.dump(noisedataframe, f)
    if onlyTomography:
        print("Tomography Ended")
        print_time()
        return
    # %%
    print("Create PER Experiment")
    perexp = experiment.create_per_experiment(circuits)

    # %%
    noise_strengths = [0,0.5,1,2]
    expectations = []
    for q in qubits:
        expect = "I"*(len(get_backend(args, return_backend_qubits=True))) #15
        expect = expect[:q] + 'Z' + expect[q+1:]
        expectations.append("".join(expect))
    print("do PER runs")
    tim = time.localtime()
    print("%s.%s. %s:%s" % (tim.tm_mday, tim.tm_mon, tim.tm_hour, tim.tm_min))
    perexp.generate(expectations = expectations, samples = persamples, noise_strengths = noise_strengths)

    # %%
    print("Run PER")
    perexp.run(executor)

    # %%
    print("Analyze PER")
    circuit_results = perexp.analyze()
    import pickle
    with open("circuit_results.pickle", "wb") as f:
        pickle.dump(circuit_results, f)


