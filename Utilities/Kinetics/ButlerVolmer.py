import numpy as np
from scipy.optimize import fsolve

"""
Butler-Volmer equation for electrochemical kinetics.

Parameters:
eta: Overpotential (V).
j: Current density (A/cm^2).
"""

class ButlerVolmer:
    def __init__(self, j0, ba, bc):
        """
        Initialize the Butler-Volmer model.

        Parameters:
        ba: Anodic Tafel slope (V/decade).
        bc: Cathodic Tafel slope (V/decade).
        j0: Exchange current density (A/cm^2).
        """
        self.ba = ba
        self.bc = bc
        self.j0 = j0

    def current_density(self, eta):
        # calculate the current density given an overpotential
        j = self.j0 * (np.exp(eta / self.ba) - np.exp(-eta / self.bc))
        return j
    
    def overpotential(self, j):
        # Calculate the overpotential given a current density
        func = lambda eta: ((j - self.j0 * (np.exp(eta / self.ba) - np.exp(-eta / self.bc)))*1e4)**2 # objective function
        sz = len(j.T)
        eta_initial_guess = np.zeros(sz)
        eta_solution = fsolve(func, eta_initial_guess)
        return eta_solution
    