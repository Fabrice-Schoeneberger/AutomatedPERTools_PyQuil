from itertools import cycle, permutations, product
import logging
 
logger = logging.getLogger("experiment")


class ProcessorSpec:
    """Responsible for interacting with the processor interface to generate the Pauli bases
    and the model terms. Also stores the mapping of virtual to physical qubits for transpilation"""

    def __init__(self, inst_map, processor, unused_qubits, plusone):
        self._n = len(inst_map)
        self._processor = processor
        self.inst_map = inst_map
        self.unused_qubits = unused_qubits
        logger.info(inst_map)
        self._connectivity = processor.sub_map(inst_map)
        logger.info(self._connectivity)
        self.plusone = plusone
        self.meas_bases = self._meas_bases()
        self.model_terms = self._model_terms()

    def _meas_bases(self):
        n = self._n
        NUM_BASES = 9

        bases = [['I']*n for i in range(NUM_BASES)]

        orderings = {"XXXYYYZZZ":"XYZXYZXYZ",
                            "XXXYYZZZY":"XYZXYZXYZ",
                            "XXYYYZZZX":"XYZXYZXYZ",
                            "XXZYYZXYZ":"XYZXZYZYX",
                            "XYZXYZXYZ":"XYZZXYYZX"}
        
        for vertex in range(n):
            #if the qubit is unused, it does not need to be tomogrophied
            if vertex in self.unused_qubits:
                for i,_ in enumerate(bases):
                    bases[i][vertex] = "I"
                continue
            #copied from Fig. S3 in van den Berg
            
            children = [edge for edge in self._connectivity if vertex in edge]
            children = set([item for c in children for item in c])
            #if a qubits is unused, it also doesn't need to be taken into account by other qubits
            predecessors = [c for c in children if c < vertex and c not in self.unused_qubits] 

            if not predecessors:
                cycp = cycle("XYZ")
                for i,_ in enumerate(bases):
                    bases[i][vertex] = next(cycp)
            #Choose p1:"XXXYYYZZZ" and p2:"XYZXYZXYZ" if one predecessor
            elif len(predecessors) == 1:
                pred, = predecessors
                #store permutation of indices so that predecessor has X,X,X,Y,Y,Y,Z,Z,Z
                _,bases = list(zip(*sorted(zip([p[pred] for p in bases], bases))))
                cycp = cycle("XYZ")
                for i,_ in enumerate(bases):
                    bases[i][vertex] = next(cycp)
            elif len(predecessors) == 2:
                pred0,pred1 = predecessors
                _,bases = list(zip(*sorted(zip([p[pred0] for p in bases], bases))))
                #list out string with permuted values of predecessor 2
                substr = [p[pred1] for p in bases]
                #match predecessor two with a permutation of example_orderings
                reordering = ""
                for perm in permutations("XYZ"):
                    substr = "".join(["XYZ"[perm.index(p)] for p in substr])
                    if substr in orderings:
                        current = orderings[substr] 
                        for i,p in enumerate(current):
                            bases[i][vertex] = p
                        break
            else:
                    raise Exception("Three or more predecessors encountered")
        bases = ["".join(b) for b in bases]
        logger.info("Created pauli bases")
        logger.info(bases)
        return [self._processor.pauli_type(string) for string in bases]

    def _model_terms(self):
        n = self._n
        model_terms = set()
        identity = ["I"]*n 
        
        #remove all unused qubits from edge_list
        trimmed_edge_list = [connection for connection in self._connectivity if not any(num in connection for num in self.unused_qubits)]
        #get all weight-two Paulis on with support on neighboring qubits
        for q1,q2 in trimmed_edge_list:
            #do not add weight-two Paulis between plusone qubits
            if q1 in self.plusone and q2 in self.plusone:
                continue
            for p1, p2 in product("IXYZ", repeat=2):
                pauli = identity.copy()
                pauli[q1] = p1
                pauli[q2] = p2
                model_terms.add("".join(pauli))

        #remove all unused qubits from indice list
        node_indices = [indice for indice in set([item for c in self._connectivity for item in c]) if indice not in self.unused_qubits]
        #get all weight-one Paulis
        for q in node_indices: 
            for p in "IXYZ":
                pauli = identity.copy()
                pauli[q] = p
                if q in self.plusone and p != "I": #remove the one weight paulis from edge qubits, as they are not needed to calc the errors in the end matrix
                    model_terms.remove("".join(pauli))    
                else: #This part of the code will only ever add more model terms if there are qubits in the system that have NO edges to any other used qubits
                    model_terms.add("".join(pauli))

        model_terms.remove("".join(identity))

        logger.info("Created model with the %s following terms:"%len(model_terms))
        logger.info(model_terms)

        return [self._processor.pauli_type(p) for p in model_terms]


    def transpile(self, circ):
        return self._processor.transpile(circ)