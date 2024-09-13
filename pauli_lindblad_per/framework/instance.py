import logging
from random import choices

logger = logging.getLogger("experiment")

class Instance:
    """This class implements the framework that both the benchmark instances and per instances
    extend. This involves storing a circuit, changing to a pauli basis, measuring, twirling the readout,
    storing the result data, untwirling the readout, and computing expectation values."""

    def __init__(self, qc, meas_basis):
        """initialize the instance with a base circuit to copy and a measurement basis

        Args:
            qc (Circuit) : the circuit to copy in order to build the instance
            meas_basis (Pauli) : The basis to measure in
        """

        self._circ = qc
        self._result = None
        self._rostring = None
        self._meas_basis = meas_basis

        self._instance() #generation code to be overridden by child classes

    def _instance(self): 
        """Generate a basic instance with basis change, readout twirling, and measurement
        """
        self._basis_change()
        self._readout_twirl()
        self._circ.measure_all()

    def add_result(self, result):
        """Attach the result of running this benchmark instance.

            result : a dictionary with binary strings as keys and frequencies as values
        """

        if self._result:
            self._result = {**self._result, **result}
        else:
            self._result = result

    def get_circuit(self):
        """Returns a copy of the circuit to be run on the hardware in the native format"""
        return self._circ.original()

    def _readout_twirl(self):
        """Implementation of readout twirling - Insertion of random Pauli-x operators to make
        diagonalize readout matrix
        """

        pauli_type = self._circ.pauli_type
        pauli_string = ""
        qubits = self._circ.qubits()
        if qubits != []:
            for i in range(1+max(qubits)):
                # We only need to do the readout_twirl for qubits where something is happening
                if i in qubits:
                    pauli_string += choices(["I", "X"])[0]
                else:
                    pauli_string += "I"
            
        readout_twirl = pauli_type(pauli_string)
        self._rostring = pauli_string
        self._circ.add_pauli(readout_twirl)

    def _basis_change(self):
        """Apply operators to change from the measurement basis into the computational basis
        """
        self._circ.compose(self.meas_basis.basis_change(self._circ).inverse()) #todo
        
    def _untwirl_result(self):
        """Return a dictionary of results with the effect of the readout twirling accounted for.
        """
        # I think this should be good, but check later again: todo
        ro_untwirled = {}
        rostring = self._rostring
        #logger.info(self._rostring)
        #logger.info(self._result)
        for key in self._result:
            newkey = "".join([{'0':'1','1':'0'}[str(bit)] if flip=="X" else str(bit) for bit,flip in zip(key,rostring)])
            ro_untwirled[newkey] = self._result[key]

        return ro_untwirled 

    def get_expectation(self, pauli):
        """Return the expectation of a pauli operator after a measurement of the circuit,
        adjusting the result for the readout twirling"""

        # I think this should be good, but check later again: todo
        pauli_type = self._circ.pauli_type
        estimator = 0
        result = self._untwirl_result()
        #compute locations of non-idetity terms (reversed indexing)
        pz = list(reversed([{pauli_type("I"):'0'}.get(p,'1') for p in pauli]))
        #compute estimator
        #logger.info("%s"%len(result.keys()))
        #logger.info("Here come the keys:")
        for key in result.keys():
            #logger.info(key)
            #compute the overlap in the computational basis
            sgn = sum([{('1','1'):1}.get((pauli_bit, key_bit), 0) for pauli_bit, key_bit in zip(pz, key)])
            #update estimator
            estimator += (-1)**sgn*result[key]
            #logger.info((-1)**sgn*result[key])

        return estimator/sum(result.values())