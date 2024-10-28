from abc import ABC, abstractmethod
from typing_extensions import Self
from typing import List, Type, Any

from primitives.instruction import Instruction
from primitives.pauli import Pauli

class Circuit(ABC):
    """A class to standardize interface with the native representation of a quantum circuit.
    Implement this class and all the abstact methods to be able to use circuits in any 
    quantum language."""

    @abstractmethod
    def copy_empty(self) -> Self:
        """Return an empty copy of the same circuit with the same qubits/qubit addresses
        such that the new circuit can be composed seamlessly with the old one"""

    @abstractmethod
    def add_instruction(self, inst : Instruction) -> None:
        """This method takes a representation of an Instruction and adds it to the end 
        of the circuit. The Instruction interface also needs to be implemented"""

    @abstractmethod
    def add_pauli(self, pauli : Pauli) -> None:
        """Append a pauli operator to the circuit"""

    @abstractmethod
    def barrier(self) -> None:
        """Add a compiler directive to maintain either side of the barrier separate"""

    @abstractmethod
    def measure_all(self, qubits : List) -> None:
        """Add a measure instruction to the desired qubits"""

    @abstractmethod
    def compose(self, other : Self) -> None:
        """Performs the composition of this circuit with other, leaving the ordering and 
        mapping of qubits on both circuits unchanged"""

    @abstractmethod
    def inverse(self) -> Self:
        """This returns the inverse of the circuit"""

    @abstractmethod
    def num_qubits(self) -> int:
        """This returns the number of qubits in the circuit"""

    @abstractmethod
    def qubits(self) -> List:
        """This returns an iterable containing the addresses of qubits in the circuit"""

    @abstractmethod
    def original(self) -> Any:
        """Returns the original circuit object in the native language"""

    @abstractmethod
    def conjugate(self, pauli) -> Pauli:
        """Returns the pauli operator conjugated by the instructions in the circuit,
        only defined when all instructions are in the Clifford group"""

    @property
    @abstractmethod
    def pauli_type(self) -> Type:
        """Returns the pauli implementation required to interact with the circuit"""

    @abstractmethod
    def __getitem__(self) -> Instruction:
        """This method provides an easy way to iterate through all of the instructions in 
        the circuit"""

    def __eq__(self, other):
        #circuits are considered equal if they contain the same instructions in any order.
        #Since this is used in the comparison of single circuit layers, it is used to 
        #remain agnostic to the native ordering of instructions on multiqubit layers
        return frozenset([inst for inst in self]) == frozenset([inst for inst in other])

    def __hash__(self): 
        #circuits are hashed based on their instructions. This is needed to make sure that
        #a layer with the same clifford profile is not benchmarked twice
        return frozenset([inst for inst in self]).__hash__()

    def __bool__(self):
        #This can be used as shorthand to tell when a circuit is empty
        return bool([inst for inst in self])

    def __str__(self): 
        #For logging, the circuit is represented as a list of instructions
        return str([inst.__str__() for inst in self])
    
from pyquil import Program
from primitives.pauli import PyQuilPauli
from primitives.instruction import PyQuilInstruction
from pyquil.paulis import PauliTerm
from pyquil.gates import H, S, X, Y, Z, FENCE
import logging
logger = logging.getLogger("experiment")

class PyquilCircuit(Circuit):
    """This is an implementation of the Circuit interface for the Qiskit API"""

    import pyquil.paulis as Pauli
    
    def __init__(self, qc: Program):
        self.qc = qc
        self._qubits = set(self.qc.get_qubit_indices())
        if len(self._qubits) > 0:
            self._num_qubits = max(self._qubits)
        else: 
            self._num_qubits = 0

    def add_instruction(self, inst : PyQuilInstruction):
        self.qc += inst.instruction
        for a in inst.support():
            self._qubits.add(a)
        if len(self._qubits) > 0:
            self._num_qubits = max(self._qubits)
        else: 
            self._num_qubits = 0

    def add_pauli(self, pauli : PyQuilPauli):
        for q,p in enumerate(pauli.to_label()):
            # No match case to make if compatible with python 3.9.17 and earlier (used in docker forest/rigetti at time of writing: 24.10.2024)
            if p == "I":
                continue
            elif p == "X":
                self.qc += X(q)
                continue
            elif p == "Y":
                self.qc += Y(q)
                continue
            elif p == "Z":
                self.qc += Z(q)
                continue
            #match p:
            #    case "I":
            #        continue
            #    case "X":
            #        self.qc += X(q)
            #        continue
            #    case "Y":
            #        self.qc += Y(q)
            #        continue
            #    case "Z":
            #        self.qc += Z(q)
            #        continue

    # These two do the same thing. PyQuil calls this operation fence, so I added that. The barrier instruction is just there for qiskit parrity
    def fence(self):
        """ Barriers or Fences are currently not supported by quilc """
        pass
        #logger.info("Barriers or Fences are currently not supported by quilc")
        #self.qc += FENCE()
    def barrier(self):
        """ Barriers or Fences are currently not supported by quilc """
        self.fence()

    def measure_all(self):
        self.qc.measure_all()

    def compose(self, other : Self):
        self.qc += other.qc
    
    def copy_empty(self):
        new_circuit = PyquilCircuit(self.qc.copy_everything_except_instructions())
        new_circuit._qubits = self._qubits
        new_circuit._num_qubits = self._num_qubits
        return new_circuit

    def inverse(self) -> Self:
        new_circuit = PyquilCircuit(self.qc.dagger())
        new_circuit._qubits = self._qubits
        new_circuit._num_qubits = self._num_qubits
        return new_circuit
    
    def qubits(self):
        return self._qubits

    def num_qubits(self):
        return self._num_qubits

    def original(self):
        return self.qc

    def conjugate(self, pauli : PyQuilPauli):
        def is_clifford_only(program: Program) -> bool: #should I enforce this, or does it not matter?
            # Set of Clifford gate names
            clifford_gates = {"X", "Y", "Z", "H", "S", "S^-1", "CNOT", "I"}
            
            # Iterate through the program's instructions
            for instruction in program.instructions:
                # Check if the instruction is a gate and not a declaration or measurement
                if instruction.name not in clifford_gates:
                    return False  # Return False if a non-Clifford gate is found

            return True  # All gates are Clifford gates
        if not is_clifford_only(self.qc):
            raise Exception("Program contains non Cliffford gates")
        
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
        pauli_term = conjugate_pauli_with_cliffords(pauli.pauli, self.qc)
        pauli_length = 0
        if pauli_term.get_qubits() != []:
            pauli_length = 1+max(pauli_term.get_qubits())
        pauli_length = max([pauli.pauli_length, pauli_length])
        return PyQuilPauli(pauli_term.pauli_string(range(pauli_length)))

    @property
    def pauli_type(self):
        return PyQuilPauli

    def __getitem__(self, item : int):
        return PyQuilInstruction(self.qc.__getitem__(item))