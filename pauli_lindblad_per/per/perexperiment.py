from primitives.circuit import PyquilCircuit
from framework.percircuit import PERCircuit
from pauli_lindblad_per.per.perrun import PERRun
from primitives.processor import PyQuilProcessor
import logging
import multiprocessing
import time

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
        
    def generate(
        self, 
        expectations, 
        samples, 
        noise_strengths
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
        self.manager = multiprocessing.Manager()  # Use a manager to handle shared data
        self._per_runs = self.manager.list() # Use a managed list for shared data between processes

        self.debug = 0
        #initialize PERRun for each PERCircuit
        for pcirc in self._per_circuits:
            #cut the generation of the PER Circuit into many threads to profit from multicore CPU performance
            process = multiprocessing.Process(target=_make_PERRUN, args=(self._inst_map, expectations, samples, noise_strengths, self.meas_bases, pcirc, self._per_runs))
            process.start()
            #This is the old non multithreding way:
            #per_run = PERRun(self._processor, self._inst_map, pcirc, samples, noise_strengths, self.meas_bases, expectations)
            #self._per_runs.append(per_run)
        while len(multiprocessing.active_children()) > 1:
            time.sleep(1)
            pass
        #changing the type from the mulitprocess list to a normal list
        self._per_runs = list(self._per_runs)
        self.manager = None

def _make_PERRUN(inst_map, expectations, samples, noise_strengths, meas_bases, pcirc, per_runs):
    """ This funktion is here to be called async from generate """
    per_run = PERRun(inst_map, pcirc, samples, noise_strengths, meas_bases, expectations)
    per_runs.append(per_run)