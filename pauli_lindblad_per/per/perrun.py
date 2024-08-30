from pauli_lindblad_per.per.perinstance import PERInstance
from pauli_lindblad_per.per.perdata import PERData
import multiprocessing
from multiprocessing.process import AuthenticationString


class PERRun:
    def __init__(self, processor, inst_map, per_circ, samples, noise_strengths, meas_bases, expectations):
        self._per_circ = per_circ
        self._pauli_type = per_circ._qc.pauli_type
        self._noise_strengths = noise_strengths
        self._samples = samples
        self._proc = processor
        self._meas_bases = meas_bases
        self._expectations = expectations
        self._inst_map = inst_map

        self._generate()

    def _get_spam(self, pauli):
        n = len(self._inst_map)
        idn = pauli.ID(n)
        spam = 1 
        for i,p in enumerate(pauli): 
            b = pauli.ID(n)
            b[i] = p
            if b != idn:
                spam *= self.spam[b]
        return spam
         
    def analyze(self):
        self._data = {} #keys are expectations
        self.spam = self._per_circ.spam
        sim_meas = {}
        for inst in self.instances:

            if not inst.meas_basis in sim_meas:
                expecs = []
                for pauli in self._expectations:
                    if inst.meas_basis.simultaneous(pauli):
                        expecs.append(pauli)
                sim_meas[inst.meas_basis] = expecs

            for basis in sim_meas[inst.meas_basis]:
                if not basis in self._data:
                    spam = self._get_spam(basis)
                    self._data[basis] = PERData(basis, spam)

                self._data[basis].add_data(inst)

        for perdat in self._data.values():
            perdat.fit()

    def get_result(self, label):
        pauli = self._pauli_type(label)
        return self._data[pauli]
    
    def _generate(self):
        #I hate that I have to do this, but you can't have class inside a multiprocess thread, that isself has processes running. 
        #So I have to put the entire function in a defined function outside the class. I havn't found a way around this.
        self.instances = _generate_but_outside_the_class(self._meas_bases, self._noise_strengths, self._samples, self._proc, self._inst_map, self._per_circ)

def _generate_but_outside_the_class(meas_bases, noise_strengths, samples, proc, inst_map, per_circ):
    manager = multiprocessing.Manager()  # Use a manager to handle shared data
    instances = manager.list()
    processes = []
    for basis in meas_bases:
        for lmbda in noise_strengths:
            for sample in range(samples):
                #cut the generation of the PER Circuit into many threads to profit from multicore CPU performance
                process = multiprocessing.Process(target=_make_PERInstance, args=(proc, inst_map, per_circ, basis, lmbda, instances))
                processes.append(process)
                process.start()
                #This is the old non multithreding way:
                #perinst = PERInstance(proc, inst_map, per_circ, basis, lmbda)
                #instances.append(perinst)

    for process in processes:
        process.join()

    return list(instances)
    
    

def _make_PERInstance(proc, inst_map, per_circ, basis, lmbda, instances):
    """ This funktion is here to be called async from _generate """
    perinst = PERInstance(proc, inst_map, per_circ, basis, lmbda)
    instances.append(perinst)