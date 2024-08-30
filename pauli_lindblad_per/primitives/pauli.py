from abc import abstractmethod, ABC
from random import choices
import logging
 
logger = logging.getLogger("experiment")

class Pauli(ABC):
    """An abstract class used as a wrapper for the native implementation of Pauli algebra.
    Whatever Pauli representation is chosen needs to be able to move back and forth between
    the Pauli string representation, and interact with circuit objects."""

    @abstractmethod
    def __init__(self, label : str) -> None:
        """Initialize a Pauli using a string consisting of the characters IXYZ"""

    @abstractmethod
    def basis_change(self, qc) -> None:
        """This method appends to the circuit qc instructions to change from the computational
        basis to the eigenbasis of the Pauli operator"""

    @abstractmethod
    def to_label(self) -> str:
        """Return a label representation of the Pauli operator without phase"""

    @abstractmethod
    def commutes(self, other) -> bool:
        """Returns true if Pauli commutes with other"""

    @abstractmethod
    def __mul__(self, other):
        """Implement multiplication of one Pauli by another Pauli, ignoring phase"""

    @abstractmethod
    def __getitem__(self, item):
        """Used to iterate through the terms in the operator"""

    @abstractmethod
    def __setitem__(self, item):
        """Used to manipulate operators term by term"""

    def simultaneous(self, other): #returns True if other can be measured simultaneously with self
        return all([p1==p2 or p2 == self.ID(1) for p1,p2 in zip(self, other)])

    def separate(self, other): #returns True if two Paulis have disjoin supports
        return all([p2 == self.ID(1) or p1 == self.ID(1) for p1,p2 in zip(self, other)])

    def nonoverlapping(self, other):
        return all([p1==p2 or p1 == self.ID(1) or p2 == self.ID(1) for p1,p2 in zip(self, other)])

    def self_conjugate(self, clifford): #returns True if self is an eigenoperator of clifford
        return self == clifford.conjugate(self)
    
    def weight(self):
        return sum([p != "I" for p in self.to_label()])

    def get_composite(self, other):
        label1 = self.to_label()
        label2 = other.to_label()
        comp = ["I"]*max([len(label1),len(label2)])
        for i, (l1,l2) in enumerate(zip(label1, label2)):
            c = l1
            if l1 == "I":
                c = l2
            comp[i] = c
        return self.__class__("".join(comp))

    @classmethod
    def random(cls, n, subset = "IXYZ"): #generates a random Pauli operator
        return cls("".join(choices(subset, k=n)))

    @classmethod
    def ID(cls, n): #returns an identity operator of the desired rank
        return cls("I"*n)
    
    def __hash__(self): #hashing does not record phase if any
        return self.to_label().__hash__()

    def __eq__(self, other): #equality does not record phase if any
        return self.to_label() == other.to_label()

    def __repr__(self):
        return self.to_label()
    
    def __str__(self) -> str:
        return self.to_label()


class PyQuilPauli(Pauli):
    """A Qiskit implementation of the Pauli algebra wrapper"""

    import pyquil.paulis as Pauli
    from pyquil.paulis import PauliTerm

    def __init__(self, name):
        self.pauli = self.PauliTerm.from_list([(p,i) for i, p in enumerate(name)])
        self.pauli_length = len(name)

    def to_label(self):
        #if self.pauli.get_qubits() == []:
        #    return ""
        if self.pauli.get_qubits() != [] and 1+max(self.pauli.get_qubits()) > self.pauli_length:
            self.pauli_length = max(self.pauli.get_qubits())+1
        return self.pauli.pauli_string(range(self.pauli_length))

    def basis_change(self, qc):
        # Import has to be here cause calling self.H(q) will transmit 2 arguments and break the code
        from pyquil.gates import H, S
        circ = qc.copy_empty() #copy_everything_except_instructions()
        label = self.to_label()
        for q in qc.qc.get_qubit_indices(): #this may not work, because qc.get_qubit_indices() only returns used qubits. But then, is it important?
            p = label[q]
            if p == "X":
                circ.qc += H(q)
            elif p == "Y":
                circ.qc += H(q)
                circ.qc += S(q)

        return circ
    
    def commutes(self,other):
        return self.Pauli.check_commutation(pauli_list=[self.pauli] , pauli_two=other.pauli)

    def __mul__(self, other):
        result = self.pauli * other.pauli
        pauli_length = 0
        if result.get_qubits() != []:
            pauli_length = 1+max(result.get_qubits())
        return PyQuilPauli(result.pauli_string(range(pauli_length)))

    def __getitem__(self, item):
        return PyQuilPauli(self.to_label()[item])

    def __setitem__(self, key, newvalue):
        label = self.to_label()
        label[key] = newvalue.pauli
        self.pauli = self.PauliTerm.from_list([(p,i) for i, p in enumerate(label)])
        #self.pauli.__setitem__(key, newvalue.pauli)