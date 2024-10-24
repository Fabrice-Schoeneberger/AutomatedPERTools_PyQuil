from scipy.optimize import nnls
import numpy as np
from matplotlib import pyplot as plt

from framework.noisemodel import NoiseModel
from tomography.layerlearning import LayerLearning
from tomography.termdata import TermData, COLORS
from primitives.circuit import Circuit
from tomography.benchmarkinstance import BenchmarkInstance, SINGLE, PAIR

import logging
from itertools import cycle

logger = logging.getLogger("experiment")

class LayerNoiseData:
    """This class is responsible for aggregating the data associated with a single layer,
    processing it, and converting it into a noise model to use for PER"""

    def __init__(self, layer : LayerLearning, sum_over_lambda=False, plusone = set(), used_qubits = None):
        self._term_data = {} #keys are terms and the values are TermDatas
        self.layer = layer
        self.sum_over_lambda=sum_over_lambda
        self.plusone = plusone
        self.used_qubits = used_qubits

        for pauli in layer._procspec.model_terms:
            pair = layer.pairs[pauli]
            self._term_data[pauli] = TermData(pauli, pair)

    def sim_meas(self, pauli):
        """Given an instance and a pauli operator, determine how many terms can be measured"""
        return [term for term in self.layer._procspec.model_terms if pauli.simultaneous(term)]

    def single_sim_meas(self, pauli, prep):
        return [term for pair,term in self.layer.single_pairs if pauli.simultaneous(term) and prep.simultaneous(pair)]

    def add_expectations(self):
        for inst in self.layer.instances:
            self.add_expectation(inst)

    def add_expectation(self, inst : BenchmarkInstance):
        """Add the result of a benchmark instance to the correct TermData object"""

        basis = inst.meas_basis
        prep = inst.prep_basis
        pair_sim_meas = {}
        single_sim_meas = {}

        if inst.type == SINGLE:

            if not basis in single_sim_meas:
                single_sim_meas[basis] = self.single_sim_meas(basis, prep)

            for pauli in single_sim_meas[basis]:
                expectation = inst.get_expectation(pauli)
                self._term_data[pauli].add_single_expectation(expectation)

        elif inst.type == PAIR:

            if not basis in pair_sim_meas:
                pair_sim_meas[basis] = self.sim_meas(basis)

            for pauli in pair_sim_meas[basis]:
                #add the expectation value to the data for this term
                expectation = inst.get_expectation(pauli)
                self._term_data[pauli].add_expectation(inst.depth, expectation, inst.type)

        
    def fit_noise_model(self):
        """Fit all of the terms, and then use obtained SPAM coefficients to make degerneracy
        lifting estimates"""

        for term in self._term_data.values(): #perform all pairwise fits
            term.fit()
        
        for pair,pauli in self.layer.single_pairs:
            self._term_data[pauli].fit_single()
            pair_dat = self._term_data[pair]
            pair_dat.fidelity = pair_dat.fidelity**2/self._term_data[pauli].fidelity

        
        logger.info("Fit noise model with following fidelities:") 
        logger.info([term.fidelity for term in self._term_data.values()])

        #get noise model from fits
        self.nnls_fit()

    def _issingle(self, term):
        return term.pauli != term.pair and term.pair in self._term_data
  
    
    def nnls_fit(self):
        """Generate a noise model corresponding to the Clifford layer being benchmarked
        for use in PER"""

        def sprod(a,b): #simplecting inner product between two Pauli operators
            return int(not a.commutes(b))

        def get_indexes(string):
            indexes = []
            for i, f in enumerate(string):
                if f != "I":
                    indexes.append(len(string)- i-1)
            return indexes

        F1 = [] #First list of terms
        F2 = [] #List of term pairs
        fidelities = [] # list of fidelities from fits
        for datum in self._term_data.values():
            pauli = datum.pauli
            indexes = get_indexes(str(pauli))
            if not self.used_qubits is None:
                skip = False
                for index in indexes:
                    if not index in self.used_qubits:
                        skip = True
                if skip:
                    continue
            F1.append(datum.pauli)
            fidelities.append(datum.fidelity)
            #If the Pauli is conjugate to another term in the model, a degeneracy is present
            if self._issingle(datum):
                F2.append(datum.pauli)
            else:
                logger.info(datum.pair)
                pair = datum.pair
                F2.append(pair)

        #create commutativity matrices
        M1 = [[sprod(a,b) for a in F1] for b in F1]
        M2 = [[sprod(a,b) for a in F1] for b in F2]
        logger.info(F1)
        logger.info(F2)
        logger.info(M1)
        logger.info(M2)
        #check to make sure that there is no degeneracy
        if np.linalg.matrix_rank(np.add(M1,M2)) != len(F1):
            raise Exception("Matrix is not full rank, something went wrong!")
       
        #perform least-squares estimate of model coefficients and return as noisemodel 
        coeffs,_ = nnls(np.add(M1,M2), -np.log(fidelities)) 

        if self.sum_over_lambda:
            paulilength = len(F1[0].to_label())
            for qubit in self.plusone:
                logger.info([f.to_label() for f in F1])
                #Filter out all model terms with the extra qubit active with at least one connected qubit
                filtered_and_cut_list = [f for f in F1 if f.to_label()[paulilength-1-qubit]!='I' and (f.to_label()[:paulilength-1-qubit]+f.to_label()[paulilength-1-qubit+1:] != "I"*(paulilength-1))]
                #In case of 2 or more connected qubits we need to seperate them
                sorted_lists = dict()
                for f in filtered_and_cut_list:
                    for i, char in enumerate(f.to_label()):
                        #find which *other* qubits are uses and sort them into lists
                        if char != "I" and i != paulilength-1-qubit:
                            if not i in sorted_lists: #make the new list if it didn't exist so far
                                sorted_lists[i] = []
                            sorted_lists[i].append(f)
                            break
                for index in sorted_lists:
                    #for each connected main qubit, filter them by x,y,z and sum their coeffs up
                    for pauli in "XYZ":
                        x = sum([coeffs[F1.index(model_term)] for model_term in sorted_lists[index] if model_term.to_label()[index]==pauli])
                        #add these coeffs to the single pauli coeffs of the main qubit. The extra fluf in coeffs[...] is just to find the correct coeff
                        coeffs[[F1.index(f) for f in F1 if f.to_label()[i]==pauli and f.to_label()[:i]+f.to_label()[i+1:]=="I"*(paulilength-1)][0]] += x
            logger.info("Summed extra qubits up and added it to main")


        self.noisemodel = NoiseModel(self.layer._cliff_layer, F1, coeffs)

    def _model_terms(self, links): #return a list of Pauli terms with the specified support
        groups = []
        for link in links:
            paulis = []
            for pauli in self._term_data.keys():
                overlap = [pauli[q].to_label() != "I" for q in link]
                support = [p.to_label() == "I" or q in link for q,p in enumerate(pauli)]
                if all(overlap) and all(support):
                    paulis.append(pauli)
            groups.append(paulis)

        return groups

    def get_spam_coeffs(self):
        """Return a dictionary of the spam coefficients of different model terms for use in 
        readout error mitigation when PER is carried out."""

        return dict(zip(self._term_data.keys(), [termdata.spam for termdata in self._term_data.values()]))

    def plot_coeffs(self, *links):
        """Plot the model coefficients in the generator of the sparse model corresponding
        to the current circuit layer"""

        coeffs_dict = dict(self.noisemodel.coeffs)
        groups = self._model_terms(links)
        fig, ax = plt.subplots()
        colcy = cycle(COLORS)
        for group in groups:
            c = next(colcy)
            coeffs = [coeffs_dict[term] for term in group]
            ax.bar([term.to_label() for term in group], coeffs, color=c)

    def graph(self, *links):
        """Graph the fits values for a certain subset of Pauli terms"""

        groups = self._model_terms(links)
        fig, ax = plt.subplots()
        for group in groups:
            for term in group:
                termdata = self._term_data[term]
                termdata.graph(ax)

        return ax

    def plot_infidelitites(self, *links):
        """Plot the infidelities of a subset of Pauli terms"""

        groups = self._model_terms(links)
        fig, ax = plt.subplots()
        colcy = cycle(COLORS)
        for group in groups:
            c = next(colcy)
            infidelities = [1-self._term_data[term].fidelity for term in group]
            ax.bar([term.to_label() for term in group], infidelities, color=c)
        return ax