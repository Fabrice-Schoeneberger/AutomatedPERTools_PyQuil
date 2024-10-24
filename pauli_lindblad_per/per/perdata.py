import numpy as np
import math
from scipy.optimize import curve_fit

class PERData:
    """Aggregate the data for runs of PER on a single circuit for multiple noise strengths and 
    a single expectation value"""

    def __init__(self, pauli, spam):
        """Initializes a data instance for a single expectation value

        Args:
            pauli (Pauli): 
            spam (dict[Pauli]): spam coefficients from learning procedure
        """

        self.pauli = pauli
        self.spam = spam
        self._data = {}
        self._dataStatistic = {}
        self._counts = {} #keep track of the number of data points for each depth

    def add_data(self, inst):
        """Add data from a PERInstance to the class. Compute the expectation value and store
        as a function of noise strength

        Args:
            inst (PERInstance): a run of PER
        """

        strength = inst.noise_strength
        expectation = inst.get_adjusted_expectation(self.pauli) #get mitigated, untwirled expectation value

        #TODO: non-local spam approximation
        expectation /= self.spam 

        #update estimator of expectation value

        self._data[strength] = self._data.get(strength, 0) + expectation

        #Gather all entries for futher analysis like error
        if self._dataStatistic.get(strength, []) == []:
            self._dataStatistic[strength] = [expectation]
        else:
            self._dataStatistic[strength] += [expectation]
        
        self._counts[strength] = self._counts.get(strength, 0)+1

    def get_std_of_strengths(self, strength):
        """return the standard deviation of the data of a specific noise strength"""
        return np.std(self._dataStatistic[strength])

    def get_expectations(self):
        """returns a list of expectation values in order of increasing noise
        """
        return [self._data[s]/self._counts[s] for s in self.get_strengths()]

    def get_strengths(self):
        """returns the noise strengths in increasing order"""
        return list(sorted(self._data.keys()))

    def fit(self):
        """Perform an exponential fit of the expectation value as a function of the noise strength.
        
        - In the large noise limit,
        the noise is completely depolarizing, and all expectation values tend toward zero.

        - A more thorough investigation of how non-clifford gates affect this behavior should be
        carried out.

        - One possibility is that taking advantage of noise tuning to predictably scale the noise
        would allow for a more reliable fit
        """ 
        
        expfit = lambda x,a,b: a*np.exp(b*x)
        #import json
        #savi = {}
        #savi["get_expectations"] = self.get_expectations()
        #savi["get_strengths"] = self.get_strengths()
        #savi["_data"] = self._data
        #savi["_dataStatistic"] = self._dataStatistic
        #savi["_counts"] = self._counts
        #with open("data.json", 'w') as file:
        #    json.dump(savi, file)
        
        try:
            popt, pcov = curve_fit(
                expfit, self.get_strengths(), 
                self.get_expectations(), 
                bounds = [(-1.5, -1),(1.5,0)], #b > 0 is enforced. Technically |a| <= 1 could also be enforced
                #to improve the accuracy, but this feels a bit like cheating 
                sigma=[self.get_std_of_strengths(strength) for strength in self.get_strengths()], #Factor in the error of the data for the error of the fit
                absolute_sigma=True
                )
        except:
            raise Exception("Fitting failed. Probably not enough data points!")
        a,b = popt
        self.expectation = a
        self.expectation_error = math.sqrt(pcov[0][0])
        return a,b

    def plot(self):
        """Plots the expectation values against an exponential fit.
        """
        # Not usable with docker forest/rigetti
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots()
        ax.errorbar(self.get_strengths(), self.get_expectations(), yerr=[self.get_std_of_strengths(strength) for strength in self.get_strengths()],  linestyle = "None", marker = "x", color = "tab:blue", capsize=5)
        a,b = self.fit()
        #pcov = self.expectation_error
        xlin = np.linspace(0, max(self.get_strengths()), 100)
        ax.plot(xlin, [a*np.exp(b*x) for x in xlin], color = "tab:blue", linestyle="--")

        return ax

    