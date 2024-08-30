class test:
    import pyquil.paulis as Pauli
    from pyquil.paulis import PauliTerm

    def __init__(self, string: str):
        self.string = string
    def __str__(self) -> str:
        return "testing"
    def __getitem__(self, item):
        return self.string[item]
    def addH(self, qc):
        from pyquil.gates import H, S
        qc += H(qubit=0)
        return qc

from pyquil import get_qc, Program
from pyquil.gates import CNOT, Z, MEASURE, S, X, Y, I, FENCE
from pyquil.api import local_forest_runtime
from pyquil.quilbase import Declare
import pyquil.paulis as Pauli
from pyquil.paulis import PauliTerm
import networkx as nx
from pyquil.quantum_processor import NxQuantumProcessor
from pyquil.noise import decoherence_noise_with_asymmetric_ro

with local_forest_runtime():
    prog = Program(Z(0), X(0))
    a = test("Hallo")
    prog = a.addH(prog)
    print(prog)
