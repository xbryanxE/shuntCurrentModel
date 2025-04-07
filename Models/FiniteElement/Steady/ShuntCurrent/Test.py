import numpy as np
import ufl
import ufl.finiteelement
from simpleShunt import shuntCurrentSimple
from dolfinx import fem, io, mesh
from dolfinx.fem import functionspace
from mpi4py import MPI
import matplotlib.pyplot as plt


# boundary conditions
bounds = [0, 0]
# Mesh deffinition
N = 100
msh = mesh.create_unit_interval(comm=MPI.COMM_WORLD, nx=N)
# resistance
R = [0.01, 10]
# Create function space
V_space = functionspace(msh, ("CG", 1))
# cell voltage
V = 1.6 # V
I = shuntCurrentSimple(V, R, V_space, msh, bounds)
print("Finished")
plt.plot(I.x.array, "o")
plt.show()