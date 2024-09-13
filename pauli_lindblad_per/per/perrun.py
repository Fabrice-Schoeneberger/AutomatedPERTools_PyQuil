from pauli_lindblad_per.per.perinstance import PERInstance
from pauli_lindblad_per.per.perdata import PERData
import multiprocessing
from multiprocessing.process import AuthenticationString
import time
import os


class PERRun:
    def __init__(self, inst_map, per_circ, samples, noise_strengths, meas_bases, expectations, processor=None):
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
        while True:
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
            
            if self.cycle_data():
                break

        for perdat in self._data.values():
            perdat.fit()

    def get_result(self, label):
        pauli = self._pauli_type(label)
        return self._data[pauli]
    
    def cycle_data(self):
        # Return True when it did a full cycle aka it is back on the start data again. Else False.
        if len(self._pickled_data) == 0:
            return True

        next_data = self._pickled_data[0]
        self._pickled_data.remove(self._pickled_data[0])
        with open(next_data, "rb") as f:
            instances = pickle.load(f)
        with open(next_data, "wb") as f:
            pickle.dump(self.instances, f)
        self.instances = instances
        self._pickled_data.append(next_data)

        self._counter += 1
        if self._counter > len(self._pickled_data):
            self._counter = 0

        return self._counter == 0
       

    def delete_data(self):
        for data in self._pickled_data:
            os.remove(data)
    
    def _generate(self):
        #I hate that I have to do this, but you can't have class inside a multiprocess thread, that isself has processes running. 
        #So I have to put the entire function in a defined function outside the class. I havn't found a way around this.
        (self.instances, self._id, self._pickled_data) = _generate_but_outside_the_class(self._meas_bases, self._noise_strengths, self._samples, self._inst_map, self._per_circ)
        self._counter = 0

import uuid
import pickle
def _generate_but_outside_the_class(meas_bases, noise_strengths, samples, inst_map, per_circ):
    manager = multiprocessing.Manager()  # Use a manager to handle shared data
    instances = manager.list()
    lock = multiprocessing.Lock()
    id = str(uuid.uuid4())
    pickled_data = []
    counter = 0
    for basis in meas_bases:
        for lmbda in noise_strengths:
            for sample in range(samples):
                #cut the generation of the PER Circuit into many threads to profit from multicore CPU performance
                process = multiprocessing.Process(target=_make_PERInstance, args=(inst_map, per_circ, basis, lmbda, instances, lock))
                process.start()

                (pickled_data, counter) = wait_for_room(instances, id, pickled_data, counter, lock, 100)
                
                #This is the old non multithreding way:
                #perinst = PERInstance(proc, inst_map, per_circ, basis, lmbda)
                #instances.append(perinst)
    (pickled_data, counter) = wait_for_room(instances, id, pickled_data, counter, lock, 1)
    manager = None
    return (list(instances), id, pickled_data)
    
def wait_for_room(instances, id, pickled_data, counter, lock, children_count=1):
    while len(multiprocessing.active_children()) > children_count:
        if len(instances) > 100:
            with lock: # This is, so that instances don't add things into the list while this pops them out. In theory it should be ok if they would but better be save
                pickel_list = []
                for i in range(100):
                    # Take from the front so to not cause any sync problems (hopefully)
                    pickel_list.append(instances.pop())
                with open(str(id)+str(counter)+".pickle", "wb") as f:
                    pickle.dump(pickel_list, f)
                pickled_data.append(str(id)+str(counter)+".pickle")
                counter += 1
                pickel_list = [] # Free the RAM

        time.sleep(60)
    return (pickled_data, counter)

def _make_PERInstance(inst_map, per_circ, basis, lmbda, instances, lock):
    """ This funktion is here to be called async from _generate """
    perinst = PERInstance(inst_map, per_circ, basis, lmbda)
    with lock:
        instances.append(perinst)