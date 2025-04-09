import numpy as np
from dolfinx import mesh
from dolfinx.fem import Function, functionspace, locate_dofs_geometrical, dirichletbc
from mpi4py import MPI


class steadyProblem:
    def __init__(self, N, problemData):
        self.problemData = problemData
        self.domain = mesh.create_interval(MPI.COMM_WORLD, N-1, [0, N-1])
        #self.domain = mesh.create_unit_interval(comm=MPI.COMM_WORLD, nx=N)
        self.V_space = functionspace(self.domain, ("Lagrange", 1))

    # set problem boundaries
    def SetBounds(self, xval):
        # Left boundary
        dofs_L = locate_dofs_geometrical(self.V_space, lambda x: np.isclose(x[0], 0.0))
        u_L = Function(self.V_space)
        u_L.interpolate(lambda x: xval[0] + 0 *x[1])
        bc_L = dirichletbc(u_L, dofs_L)
        # Right boundary
        dofs_R = locate_dofs_geometrical(self.V_space, lambda x: np.isclose(x[0], self.N-1))
        u_R = Function(self.V_space)
        u_R.interpolate(lambda x: xval[1] + 0 *x[1])
        bc_R = dirichletbc(u_R, dofs_R)
        return [bc_L, bc_R]

    
    
    
    
    



