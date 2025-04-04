from primitives.circuit import PyquilCircuit
from framework.percircuit import PERCircuit
from pauli_lindblad_per.per.perrun import PERRun
from primitives.processor import PyQuilProcessor
import logging
import multiprocessing
import time
import numpy as np

logger = logging.getLogger("experiment")

class PERExperiment:
    """This class plays the role of the SparsePauliTomographyExperiment class but for the
    generation, aggregation, and analysis of PER database
    
    class functions:
    - get the minimal number of measurement bases required to reconstruct desired observables
    - initialize generation of PER circuits to estimate each expectation for each circuit
        at the desired noise strength 
    - Pass the circuits to user-defined run method for execution of
    - Process results and return for display
    """
    
    def __init__(self, circuits, inst_map, noise_data_frame, backend = None, processor = None):
        """Initializes a PERExperiment with the data that stays constant for all circuits/
        noise strengths/expectation values

        Args:
            circuits (Any): Circuits to run with PER
            inst_map (List): Mapping of virtual qubits to physical qubits
            noise_data_frame (NoiseDataFrame): Noise models learned from tomography
            backend (Any): Backend to use for transpilation. None if passing an initialize processor
            processor (Processor) : Backend to use for transpilation. None if passing a native backend
        """
        circuit_interface = None

        #check if circuits have implementable type, and initialize processor
        if circuits[0].__class__.__name__ == "Program":
            circuit_interface = PyquilCircuit
            if backend:
                self._processor = PyQuilProcessor(backend)
        else:
            raise Exception("Unsupported circuit type")
        if not backend:
            self._processor = processor  
        self.pauli_type = circuit_interface(circuits[0]).pauli_type


        self.noise_data_frame = noise_data_frame #store noise data
        #Geerate list of PER circuits and assign noise models to layers 
        per_circuits = []
        for circ in circuits:
            circ_wrap = circuit_interface(circ) #wrap Circuit object
            per_circ = PERCircuit(circ_wrap) #create a PER circuit
            per_circ.add_noise_models(noise_data_frame) #add associated noise models
            per_circuits.append(per_circ)

        self._per_circuits = per_circuits
        self._inst_map = inst_map

    def get_meas_bases(self, expectations):
        """Return the minimal set of bases needed to reconstruct the desired expectation values

        Args:
            expectations (Pauli): The desired Pauli expectation values
        """

        meas_bases = []
        #iterate through expectations
        for pauli in expectations:
            for i,base in enumerate(meas_bases): #iterate through bases
                if base.nonoverlapping(pauli): #if nonoverlapping, compose into last basis
                    meas_bases[i] = base.get_composite(pauli)
                    break
            else:
                meas_bases.append(pauli) #if no break is reached, append to end
        logger.info(meas_bases)
        self.meas_bases = meas_bases

    def run(self, executor, shots, do_cross_talk=False, apply_cross_talk=None):
        """pass a list of circuit in the native language to the executor method and await results

        Args:
            executor (method): list of circuits -> Counter of results
        """

        #aggregate all instances into a list
        instances = []
        for run in self._per_runs:
            instances += run.instances
       
        #get circuits in native representation
        circuits = [inst.get_circuit() for inst in instances] 
        
        if do_cross_talk and apply_cross_talk:
            circuits = apply_cross_talk(circuits, self._processor._qpu)

        logger.info(len(circuits))

        #pass circuit to executor
        results = executor(circuits, self._processor._qpu, shots)
       
        #add results to instances 
        for inst, res in zip(instances, results):
            inst.add_result(res)
        
        cycling = False

        for run in self._per_runs:
            cycling = run.cycle_data() # All runs should return the same data so this should be fine

        if not cycling:
            self.run(executor) # recursivly call this function until all instances have been dealt with

    def analyze(self):

        #run analysis on all instances and return results for each circuit
        for run in self._per_runs:
            run.analyze()

        return self._per_runs

    def get_overhead(self, layer, noise_strength):
        return self._per_circuits[layer].overhead(noise_strength)
    
    def _calculate_gamma_scale(self, pcirc, noise_strengths, samples_per_circuit):
        # Calculate gamma for every noise strength
        gammas = []
        for noise_strength in noise_strengths:
            gamma = 1
            if noise_strength < 1:
                for layer in pcirc._layers:
                    gamma *= np.exp(2*(1-noise_strength)*sum([abs(value[1]) for value in self.noise_data_frame.noisemodels[layer.cliff_layer].coeffs]))
            gammas.append(gamma)
        # Distribute and create scale. Round to nearest int number
        samples_scale = {noise_strength: round(samples_per_circuit*gamma/sum(gammas)) for noise_strength, gamma in zip(noise_strengths, gammas)}
        # To make sure all circuits have a minimum number of circuits, check for noise_strengths with lower counts and set them higher, without increasing the overall sample count
        marked = [] # Save all the noise_strengths that need to be set to the minimum
        # This is a temporary figure. The number needs to be big enougth to account for the twirling precision. Once a definitive answer is figured out, it comes here
        minimum_circuit_count = max(samples_per_circuit/(2.5*len(noise_strengths)),1) # Define the minimum
        new_marked = True #loop has to run at least one. Python does not have a do while loop
        while new_marked:
            new_marked = False
            for key in samples_scale:
                if samples_scale[key] < minimum_circuit_count: # Check all sample sizes
                    new_marked = True # And mark the one that are to low
                    marked.append(key) 
            if new_marked: # Recalculate the scale with the min number of circuits for every marked strength taken of the pool, for the rest of the strengths
                samples_scale = {noise_strength: round((samples_per_circuit-len(marked)*(minimum_circuit_count))*gamma/sum([gamma for noise_strength, gamma in zip(noise_strengths, gammas) if noise_strength not in marked])) for noise_strength, gamma in zip(noise_strengths, gammas) if noise_strength not in marked}
            # Loop until nothing more gets marked
        for mark in sorted(marked): # Add the marked back at the minimum count
            samples_scale[mark] = round(minimum_circuit_count)
        return samples_scale
        
    def generate(
        self, 
        expectations, 
        samples_per_circuit, 
        noise_strengths,
        do_multithreading
        ):
        """Initiate the generation of circuits required for PER

        Args:
            noise_strengths (list[int]): strengths of noise for PER fit
            expectations (list[str]): expectation values to reconstruct
            samples (int): number of samples to take from distribution
        """

        #Convert string labels to Pauli representation
        expectations = [self.pauli_type(label) for label in expectations]

        #get minimal set of measurement bases
        self.get_meas_bases(expectations)
        if do_multithreading:
            self.manager = multiprocessing.Manager()  # Use a manager to handle shared data
            self._per_runs = self.manager.list() # Use a managed list for shared data between processes
        else:
            self._per_runs = []

        self.debug = 0
        #initialize PERRun for each PERCircuit
        for pcirc in self._per_circuits:
            samples = self._calculate_gamma_scale(pcirc, noise_strengths, samples_per_circuit)
            if do_multithreading:
                #cut the generation of the PER Circuit into many threads to profit from multicore CPU performance
                process = multiprocessing.Process(target=_make_PERRUN, args=(self._processor, self._inst_map, expectations, samples, noise_strengths, self.meas_bases, pcirc, self._per_runs, do_multithreading))
                process.start()
            else:
            #This is the old non multithreading way:
                per_run = PERRun(self._processor, self._inst_map, expectations, samples, noise_strengths, self.meas_bases, pcirc, self._per_runs, do_multithreading)
                self._per_runs.append(per_run)
        while len(multiprocessing.active_children()) > 1:
            time.sleep(1)
            pass
        #changing the type from the mulitprocess list to a normal list
        self._per_runs = list(self._per_runs)
        self.manager = None

def _make_PERRUN(processor, inst_map, expectations, samples, noise_strengths, meas_bases, pcirc, per_runs, do_multithreading):
    """ This funktion is here to be called async from generate """
    per_run = PERRun(processor, inst_map, pcirc, samples, noise_strengths, meas_bases, expectations, do_multithreading)
    per_runs.append(per_run)