from abc import ABC, abstractmethod
from typing import List
import logging
logger = logging.getLogger("experiment")

class Instruction(ABC):

    @abstractmethod
    def name(self) -> str:
        """Returns a name to identify the type of instruction"""

    @abstractmethod
    def support(self) -> List:
        """Return an ordered list of the qubits the instruction affects nontrivially"""

    @abstractmethod
    def ismeas(self):
        """Identifies whether the instruction is a measurement"""
        pass

    def weight(self):
        """Return the number of qubits nontrivially affected by the instruction"""
        return len(self.support())

    def __hash__(self):
        #Hash the instruction based on its name and its (ordered) support
        return (self.name(),self.support()).__hash__()

    def __eq__(self, other):
        #Instructions are equal if they have the same name and the same support
        return self.name() == other.name() and self.support() == other.support()

    def __str__(self):
        return str((self.name(), self.support()))


class PyQuilInstruction(Instruction):

    def __init__(self, instruction):
        self.instruction = instruction
    
    def support(self):
        if self.instruction._InstructionMeta__name == 'Declare' or self.instruction._InstructionMeta__name == 'Halt':
            return tuple([])
        if hasattr(self.instruction, "qubits"):
            return tuple([q.index for q in self.instruction.qubits])
        elif hasattr(self.instruction, "qubit"):
            return tuple([self.instruction.qubit.index])
        else:
            raise Exception(str(self.instruction) + " has no qubits")

    def name(self):
        if self.instruction._InstructionMeta__name == "Measurement":
            return "Measurement"
        elif self.instruction._InstructionMeta__name == "Halt":
            return "Halt"
        else:
            return self.instruction.name

    def ismeas(self):
        return self.name() == "Measurement"