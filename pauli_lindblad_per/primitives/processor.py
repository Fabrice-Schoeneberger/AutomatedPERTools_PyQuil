from abc import ABC, abstractmethod
from primitives.circuit import Circuit, PyquilCircuit
from primitives.pauli import PyQuilPauli
import logging

logger = logging.getLogger("experiment")

class Processor(ABC):
    """A wrapper for interacting with a qpu backend. This object is responsible for
    reporting the processor topology and transpiling circuits into the native gate set."""

    @abstractmethod
    def sub_map(self, qubits : int):
        """Return an sub edge list in the form of tuples of ints representing connections
        between qubits at those hardware addresses"""

    @abstractmethod
    def transpile(self, circuit : Circuit):
        """Transpile a circuit into the native gateset"""
    
    @property
    @abstractmethod
    def pauli_type(self):
        """Returns the native Pauli type associated"""


class PyQuilProcessor(Processor):
    """Implementaton of a processor wrapper for the PyQuil API"""

    def __init__(self, backend):
        self._qpu = backend

    def sub_map(self, inst_map):
        return [pair for pair in self._qpu.quantum_processor.qubit_topology().to_directed().edges() if any(p in inst_map for p in pair)]

    #PyQuil calls this process compile and not transpile, I keep both names for parity sake
    def compile(self, circuit : PyquilCircuit):
        return PyquilCircuit(self._qpu.compile(circuit.qc))
    def transpile(self, circuit : PyquilCircuit):
        return self.compile(circuit)

    @property
    def pauli_type(self):
        return PyQuilPauli