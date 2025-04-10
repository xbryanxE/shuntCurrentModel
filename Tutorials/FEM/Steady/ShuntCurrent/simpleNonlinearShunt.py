import numpy as np
import sys 
import os
import matplotlib.pyplot as plt
from time import time

# add utilities
path_to_kinetics = os.path.expanduser("~/StackSim/Utilities/Kinetics")
path_to_fem = os.path.expanduser("~/StackSim/Models/Steady/")
path_to_modules = os.path.expanduser("~/StackSim/Models/Steady/modules/")

sys.path.append(os.path.abspath(path_to_kinetics))
sys.path.append(os.path.abspath(path_to_fem))
sys.path.append(os.path.abspath(path_to_modules))

from shunt.linear.shuntConstraint import linear_u_test
from ShuntCurrents import shuntCurrent
from ButlerVolmer import ButlerVolmer

class control:
    def __init__(self):
        # channel resistance
        self.Rch = 5 # Ohm
        # manifold resistance
        self.Rm = 0.1 # Ohm
        # contact resistance
        self.Rc = 1e-4 # Ohm
        # current
        self.Iw = 100 # A
        # electrode surface area
        self.As = 1e4 # cm^2

    def fnc_Rch(self, x, data):
        # channel resistance as a function of x
        return self.Rch + x*0.
    
    def fnc_Rm(self, x, data):
        # channel resistance as a function of x
        return self.Rm + x*0.

class equilibrium:
    def __init__(self):
        self.Erev = 1.23 # V

class problemData:
    def __init__(self, Control, Equilibrium, Kinetics):
        self.Kinetics = Kinetics
        # control definition
        self.Control = Control
        # equilibrium definition
        self.Equilibrium = Equilibrium


if __name__ == "__main__":
    start = time()
    # Define the problem data
    control = control()
    equilibrium = equilibrium()
    kinetics = ButlerVolmer(1.2e-4, 0.04, 0.12) # (j0, ba, bc)
    problem_data = problemData(control, equilibrium, kinetics)
    
    # Create a simple test case
    N = 300
    problem = shuntCurrent(N, problem_data)
    # Initialize boundaries
    bounds = [0, 0]
    # Define the functions
    Fncs = [problem.problemData.Control.fnc_Rch, problem.problemData.Control.fnc_Rm]
    # Call the SimpleShunt method
    [uh, V] = problem.simpleNonlinearShunt(linear_u_test, Fncs, bounds, [1, 1])
    finish = time()
    # plot results
    plt.plot(problem.domain.geometry.x[:,0]+1, uh.x.array, "o-")
    plt.ylabel("Manifold current (A)")
    plt.xlabel("Stack relative location")
    plt.show()

    print("Finished")
    print("Elapsed time: ", finish - start, "s")