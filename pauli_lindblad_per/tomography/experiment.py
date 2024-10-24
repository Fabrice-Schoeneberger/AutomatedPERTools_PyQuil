from tomography.processorspec import ProcessorSpec
from tomography.layerlearning import LayerLearning
from tomography.analysis import Analysis
from framework.percircuit import PERCircuit
from per.perexperiment import PERExperiment

from typing import List, Any
import logging
 
logging.basicConfig(filename="experiment.log",
                    format='%(asctime)s %(message)s',
                    filemode='w')

logger = logging.getLogger("experiment")
logger.setLevel(logging.INFO)

from primitives.circuit import PyquilCircuit
from primitives.processor import PyQuilProcessor
import pickle

class SparsePauliTomographyExperiment:
    """This class carries out the full experiment by creating and running a LayerLearning
    instance for each distinct layer, running the analysis, and then returning a PERCircuit
    with NoiseModels attached to each distinct layer"""

    def __init__(self, circuits, inst_map, backend, sum_over_lambda=False):

        circuit_interface = None
        #Make sure it's a quantumcircuit as others don't work
        if circuits[0].__class__.__name__ == "Program":
            circuit_interface = PyquilCircuit
            processor = PyQuilProcessor(backend)
        else:
            raise Exception("Unsupported circuit type") 
    
        self._profiles = set()
        used_qubits = set()
        for circuit in circuits: 
            for i in circuit.get_qubit_indices():
                used_qubits.add(i)
            
            circ_wrap = circuit_interface(circuit)
            parsed_circ = PERCircuit(circ_wrap)
            for layer in parsed_circ._layers:
                if layer.cliff_layer:
                    self._profiles.add(layer.cliff_layer)
                    
        plusone = set() #Here come the extra
        
        #Now see which qubits are unused by all circuits
        unused_qubits = [bit for bit in inst_map if bit not in used_qubits]
        logger.info("The following Qubits were determinded unused")
        logger.info(unused_qubits)

        logger.info("Generated layer profile with %s layers:"%len(self._profiles))
        for layer in self._profiles:
            logger.info(layer)

        self._procspec = ProcessorSpec(inst_map, processor, unused_qubits, plusone)
        self.instances = []
        self._inst_map = inst_map
        self.unused_qubits = unused_qubits
        self.sum_over_lambda = sum_over_lambda
        self._layers = None

        self._layers = []
        for l in self._profiles:
            learning = LayerLearning(l,self._procspec)
            self._layers.append(learning)

        self.analysis = Analysis(self._layers, self._procspec, sum_over_lambda=sum_over_lambda, plusone=plusone, used_qubits=used_qubits)

    def generate(self, samples, single_samples, depths):
        """This method is used to generate the experimental benchmarking procedure. The samples
        are the number of times to sample from the Pauli twirl. The single_samples controls
        how many twirl samples to take from the degeneracy-lifting measurements. It may desirable
        to make this higher since the error on these measurements will generally be higher.
        The depths control the different circuit depths to use for the exponential fits."""

        if len(depths) < 2:
            raise Exception("Exponental fit requires 3 or more depth data points.")

        for l in self._layers:
            l.procedure(samples, single_samples, depths)

    def run(self, executor, shots, do_cross_talk=False, apply_cross_talk=None):
        """This method produces a list of circuits in the native representation, passes them 
        as a list to the executor method, and associates the result with the benchmark instances
        that produced it"""


        instances = []
        for l in self._layers:
            instances += l.instances

        circuits = [inst.get_circuit() for inst in instances]

        if do_cross_talk and apply_cross_talk:
            circuits = apply_cross_talk(circuits, self._procspec._processor._qpu)

        results = executor(circuits, self._procspec._processor._qpu, shots)

        for res,inst in zip(results, instances): #TODO: find out if order can be preserved
            inst.add_result(res)

    def analyze(self):
        """Runs analysis on each layer representative and stores for later plotting/viewing"""
        self.analysis.analyze()
        return self.analysis.noisedataframe

    def create_per_experiment(self, circuits : Any) -> PERExperiment:
        experiment = PERExperiment(circuits, self._inst_map, self.analysis.noisedataframe, backend = None, processor = self._procspec._processor)
        return experiment

    def save(self):
        raise NotImplementedError()

    def load(self):
        raise NotImplementedError()