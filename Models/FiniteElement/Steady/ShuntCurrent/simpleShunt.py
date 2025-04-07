import numpy as np
from mpi4py import MPI
import ufl
from dolfinx import fem, io
from dolfinx.fem import locate_dofs_geometrical, dirichletbc, Function
from dolfinx.fem.petsc import LinearProblem
from ufl import ds, dx, grad, inner


def systemBoundary(xval, V_space):
    # Left boundary
    dofs_L = locate_dofs_geometrical(V_space, lambda x: np.isclose(x[0], 0.0))
    u_L = Function(V_space)
    u_L.interpolate(lambda x: xval)
    bc_L = dirichletbc(u_L, dofs_L)
    # Right boundary
    dofs_R = locate_dofs_geometrical(V_space, lambda x: np.isclose(x[0], 1.0))
    u_R = Function(V_space)
    u_R.interpolate(lambda x: xval)
    bc_R = dirichletbc(u_R, dofs_R)
    return [bc_L, bc_R]

def ProblemDefinition(msh, V, R, V_space, bcs):
    u = ufl.TrialFunction(V_space)
    v = ufl.TestFunction(V_space)
    x = ufl.SpatialCoordinate(msh)
    f = V
    a = R[0] * inner(grad(u), grad(v)) * dx + R[1] * inner(u, v) * dx
    L = inner(f, v) * dx
    problem = LinearProblem(a, L, bcs=bcs)
    
    return problem

def shuntCurrentSimple(V, R, V_space, msh, bounds):
    bcs = systemBoundary(bounds, V_space)
    # Define the problem
    problem = ProblemDefinition(msh, V, R, V_space, bcs)
    # solve
    uh = problem.solve()
    return uh


    