import numpy as np
from mpi4py import MPI
import ufl
from dolfinx import fem, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from steadyProblem import steadyProblem


class shuntCurrent(steadyProblem):
    def __init__(self, N, problemData):
        super().__init__(N, problemData)
        self.N = N

    def SimpleShunt(self, bcs):
        # get constant channel resistances
        Rch = self.problemData.Control.Rch
        # get constsant manifold resistances
        Rm = self.problemData.Control.Rm
        # Calculate a constant cell voltage
        Iw = self.problemData.Control.Iw 
        V = self.CellVoltage(Iw)
        # Problem definition
        u = ufl.TrialFunction(self.V_space)
        v = ufl.TestFunction(self.V_space)
        f = fem.Constant(self.domain, default_scalar_type(V))
        a = (Rch * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
             + Rm * u * v * ufl.dx)
        L = f * v * ufl.dx 
        femProblem = LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        # solve the system
        uh = femProblem.solve()
        return [uh, V]

    def LinearShunt(self, u_test, bounds, idx):
        Rch = self.problemData.Control.Rch
        Rm = self.problemData.Control.Rm
        r = 1
        while r > 1e-6:
            bcs = self.SetBounds(bounds)
            [uh_, V] = self.SimpleShunt(bcs)
            [bounds, r] = u_test(uh_, [Rch, Rm], V, bounds, idx)
        return [uh_, V]

    def CellVoltage(self, Ic):
        # Define the problem
        As = self.problemData.Control.As # eletrode surface
        Erev = self.problemData.Equilibrium.Erev # thermodynamic potential
        eta = self.problemData.Kinetics.overpotential(np.array([Ic/As])) # kinetic overpotential
        Vohm = self.problemData.Control.Rc * Ic
        V = Erev + eta[0] + Vohm
        return V



    