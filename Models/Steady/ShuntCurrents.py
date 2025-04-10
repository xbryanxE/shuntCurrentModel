import numpy as np
from mpi4py import MPI
import ufl
from dolfinx import fem, default_scalar_type
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from dolfinx.fem import Function, FunctionSpace
from dolfinx.nls.petsc import NewtonSolver

from steadyProblem import steadyProblem


class shuntCurrent(steadyProblem):
    def __init__(self, N, problemData):
        super().__init__(N, problemData)
        self.N = N

    def SimpleShunt(self, bcs): # Working properly
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

    def LinearShunt(self, u_test, bounds, idx): # Working properly
        Rch = self.problemData.Control.Rch
        Rm = self.problemData.Control.Rm
        r = 1
        while r > 1e-6:
            bcs = self.SetBounds(bounds)
            [uh_, V] = self.SimpleShunt(bcs)
            [bounds, r] = u_test(uh_, [Rch, Rm], [V, V], bounds, idx)
        return [uh_, V]

    def simpleNonlinearShunt(self, u_test, Fncs, bounds, idx): # Working properly
        Rch_fnc = Fncs[0] # function for channel resistance
        Rm_fnc = Fncs[1] # function for manifold resistance
        # set initial boundaries
        bcs = self.SetBounds(bounds)
        # evaluate initial guess
        [uh_init, _] = self.SimpleShunt(bcs)
        # Calculate the Cell voltage across the stack
        Ic = self.problemData.Control.Iw - uh_init.x.array
        Vcells = self.CellVoltage(Ic)
        # Define the problem
        S_fnc = Function(self.V_space) # Cell Voltage
        S_fnc.x.array.setfield(Vcells, dtype=np.float64) # set initial value
        u = ufl.TrialFunction(self.V_space) # Trial function
        v = ufl.TestFunction(self.V_space)
        bc_r = 1
        while bc_r > 1e-6:
            r = 1
            while r > 1e-6:
                # calc_resistances
                Rch = Rch_fnc(uh_init, self.problemData)
                Rm = Rm_fnc(uh_init, self.problemData)
                a = (Rch * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
                + Rm * u * v * ufl.dx)
                L = S_fnc * v * ufl.dx
                femProblem = LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
                # solve the system
                uh = femProblem.solve()
                # Update the cell voltage
                Ic = self.problemData.Control.Iw - uh.x.array
                Vcells = self.CellVoltage(Ic)
                # Update source term
                S_fnc.x.array.setfield(Vcells, dtype=np.float64)
                # compute mean square error
                r = np.sum((uh.x.array - uh_init.x.array)**2)
                print("Shunt current residual: ", r)
                # update the initial guess
                uh_init.x.array.setfield(uh.x.array, dtype=np.float64)
            # check if the boundary conditions are satisfied
            [bounds, bc_r] = u_test(uh, [Rch, Rm], S_fnc.x.array, bounds, idx)
            # update the boundaries
            bcs = self.SetBounds(bounds)
            bc_r = bc_r.value()
            print("Boundary condition residual: ", bc_r)
        print("Shunt current converged")
        return [uh, S_fnc]

    def NonlinearShunt(self, u_test, Fncs, bounds, idx): # Working properly
        Rch_fnc = Fncs[0] # function for channel resistance
        Rm_fnc = Fncs[1] # function for manifold resistance
        # set initial boundaries
        bcs = self.SetBounds(bounds)
        # evaluate initial guess
        [uh_init, _] = self.SimpleShunt(bcs)
        # Define terms
        # source term
        S_fnc = Function(self.V_space) # Cell Voltage
        # Channel resistance
        S_Rch = Function(self.V_space) # Channel resistance
        # Manifold resistance
        S_Rm = Function(self.V_space) # Manifold resistance
        # 
        u = ufl.TrialFunction(self.V_space) # Trial function
        v = ufl.TestFunction(self.V_space) # Test function
        bc_r = 1
        while bc_r > 1e-6:
            r = 1
            while r > 1e-6:
                # Calculate the Cell voltage across the stack
                Ic = self.problemData.Control.Iw - uh_init.x.array
                Vcells = self.CellVoltage(Ic)
                S_fnc.x.array.setfield(Vcells, dtype=np.float64) 
                # update channel resistance
                Rch_vals = Rch_fnc(uh_init.x.array, self.problemData)
                S_Rch.x.array.setfield(Rch_vals, dtype=np.float64)
                # update manifold resistance
                Rm_vals = Rm_fnc(uh_init.x.array, self.problemData)
                S_Rm.x.array.setfield(Rm_vals, dtype=np.float64)
                # Equation
                a = (S_Rch * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
                     + ufl.dot(ufl.grad(S_Rch), ufl.grad(u)) * v * ufl.dx
                     + S_Rm * u * v * ufl.dx)
                L = S_fnc * v * ufl.dx
                femProblem = LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
                # solve the system
                uh = femProblem.solve()
                # compute mean square error
                r = np.sum((uh.x.array - uh_init.x.array)**2)
                print("Shunt current residual: ", r)
                # update the initial guess
                uh_init.x.array.setfield(uh.x.array, dtype=np.float64)
            [bounds, bc_r] = u_test(uh, [Rch_fnc(uh_init, self.problemData), Rm_fnc(uh_init, self.problemData)], S_fnc.x.array, bounds, idx)
            # update the boundaries
            bcs = self.SetBounds(bounds)
            bc_r = bc_r.value()
            print("Boundary condition residual: ", bc_r)
        print("Shunt current converged")
        return [uh, S_fnc]

    def CellVoltage(self, Ic): # Working
        # Define the problem
        As = self.problemData.Control.As # eletrode surface
        Erev = self.problemData.Equilibrium.Erev # thermodynamic potential
        if len(np.array([Ic]).T) < 2:
            eta = self.problemData.Kinetics.overpotential(np.array([Ic/As])) # kinetic overpotential
        else:
            eta = self.problemData.Kinetics.overpotential(Ic/As) # kinetic overpotential
        Vohm = self.problemData.Control.Rc * Ic
        V = Erev + eta[0] + Vohm
        return V



    