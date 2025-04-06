import numpy as np

"""
Butler-Volmer equation for electrochemical kinetics.

Parameters:
eta: Overpotential (V).
j: Current density (A/cm^2).
"""

class Tafel:
    def __init__(self, j0, b):
        """
        Initialize the Butler-Volmer model.

        Parameters:
        b: Combined Tafel slope (V/decade).
        j0: Exchange current density (A/cm^2).
        """
        self.b = b
        self.j0 = j0

    def current_density(self, eta):
        # calculate the current density given an overpotential
        j = self.j0 * np.exp(eta / self.b)
        return j
    
    def overpotential(self, j):
        # calculate the overpotential given a current density
        eta_solution = self.b * np.log10(j / self.j0)
        return eta_solution
    