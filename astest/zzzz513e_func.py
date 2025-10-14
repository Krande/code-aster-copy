# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------


from code_aster.Commands import *
from code_aster import CA


###################################################################################
#
#   Solve coupling problem with HHO
#   u -> displacement, dx -> macro-damage, d -> micro-damage
#
#   Continuous:
#   R_u(u,dx,d) = (PK1(u,dx,d), grad v) = 0, \forall v
#   R_dx(u,dx,d) = (Ax * grad dx, grad px) + (Hx*(dx-d), px), \forall px
#
#   d(t+dt) = min(1, max(d(t), (2 phi(u) + Hx(dx+1/beta) ) / (2 phi(u) + Hx)))
#
#   HHO:
#   sum_{T \in Th} (PK(huT, d_T), GkT(hvT))_T + stab(huT, hvT) = 0
#   sum_{T \in Th} (Ax * GkT(hdxT), GkT(hpxT))_T + stab(hdxT, hpxT) +
#      (H * (dx_T - dT), px_T)_T = 0
#
####################################################################################

from code_aster.Utilities import force_list
from code_aster.Solvers import TimeStepper


class CoupledState:
    u = dx = d = None
    sief_elga = vari_elga = None
    u_nodes = dx_nodes = d_nodes = None
    time = None

    def __init__(self, u, dx, d, time=0.0):
        self.u = u.copy()
        self.dx = dx.copy()
        self.d = d.copy()
        self.time = time

    def projectOnLagrangeSpace(self, hho_meca, hho_dama):
        self.u_nodes = hho_meca.projectOnLagrangeSpace(self.u)
        self.dx_nodes = hho_dama.projectOnLagrangeSpace(self.dx)
        self.d_nodes = hho_dama.projectOnLagrangeSpace(self.d)

    def copy(self):
        dup = CoupledState(self.u, self.dx, self.d, self.time)
        dup.u_nodes = self.u_nodes.copy()
        dup.dx_nodes = self.dx_nodes.copy()
        dup.d_nodes = self.d_nodes.copy()
        dup.vari_elga = self.vari_elga.copy()
        dup.sief_elga = self.sief_elga.copy()

        return dup


class MecaSolver:
    """Solve the mecanical problem"""

    model = mater = loads = None

    def __init__(self, param):
        self.model = param["MODELE"]
        self.mater = param["MATER"]
        self.loads = param["EXCIT"]

    def create_material(self, d_nodes):
        return AFFE_MATERIAU(
            MAILLAGE=self.model.getMesh(),
            AFFE=_F(TOUT="OUI", MATER=self.mater),
            AFFE_VARC=(_F(NOM_VARC="TEMP", CHAM_GD=d_nodes, VALE_REF=0.0),),
        )

    def solve(self, state_prev, state_curr):
        """Solve mechanical problem"""

        material = self.create_material(state_curr.dx_nodes)

        timeList = DEFI_LIST_REEL(VALE=(state_prev.time, state_curr.time))

        opt = {}
        if state_curr.sief_elga is not None:
            opt["ETAT_INIT"] = _F(
                DEPL=state_prev.u, SIGM=state_prev.sief_elga, VARI=state_prev.vari_elga
            )

        resuMeca = STAT_NON_LINE(
            MODELE=self.model,
            CHAM_MATER=material,
            INCREMENT=_F(LIST_INST=timeList),
            METHODE="NEWTON",
            EXCIT=self.loads,
            **opt
        )

        return (
            resuMeca.getField("DEPL", 1).copy(),
            resuMeca.getField("SIEF_ELGA", 1).copy(),
            resuMeca.getField("VARI_ELGA", 1).copy(),
        )


class DamageSolver:
    """Solve the damage problem"""

    model = mater = loads = None

    def __init__(self, param):
        self.model = param["MODELE"]
        self.mater = param["MATER"]
        self.loads = force_list(param["EXCIT"])

    def create_material(self, u_nodes):
        therMate = DEFI_MATERIAU(
            THER=_F(
                LAMBDA=self.mater["LAMBDA"],
                RHO_CP=self.mater["RHO_CP"] * (1.0 + u_nodes.norm("NORM_2")),
            )
        )

        return AFFE_MATERIAU(MAILLAGE=self.model.getMesh(), AFFE=_F(TOUT="OUI", MATER=therMate))

    def evalMacroDamage(self, u, dx, d):
        """Evaluate macrodamage"""

        # this is not the true expression
        return d

    def solve(self, state_prev, state_curr):
        """Solve damaage problem"""

        material = self.create_material(state_curr.u_nodes)

        # create PhysicalProblem
        phys_pb = CA.PhysicalProblem(self.model, material)

        for load in self.loads:
            if isinstance(load["CHARGE"], CA.ThermalDirichletBC):
                phys_pb.addDirichletBC(load["CHARGE"])

        # compute DOF numbering
        phys_pb.computeDOFNumbering()

        # create discrete computation
        disc_comp = CA.DiscreteComputation(phys_pb)
        hho = CA.HHO(phys_pb)

        # compute K = (lamda * GkT(hdT), GkT(hpT))_T + lambda * stab(hdT, hpT)
        matK = disc_comp.getLinearStiffnessMatrix(assembly=True)

        # compute M = (rho_cp * d_T, p_T) _T
        matM = disc_comp.getMassMatrix(assembly=True)

        diriBCs = disc_comp.getDirichletBC()

        # linear solver
        mySolver = CA.MumpsSolver()

        # init with previous converged solution
        d_curr = state_curr.d
        dx_curr = state_curr.dx

        print("Newton solver:")
        for i in range(100):
            # the residual is ok
            Resi = matK * dx_curr + matM * (dx_curr - d_curr)
            # the jacobian is inexacte. Use \tilde(Hx) insted of Hx
            Jaco = matK + matM

            print("*Iter %d: residual %f" % (i, Resi.norm("NORM_2")))
            if Resi.norm("NORM_2") < 10e-8:
                return dx_curr, d_curr

            mySolver.factorize(Jaco)
            dx_incr = mySolver.solve(-Resi, diriBCs)
            dx_curr += dx_incr

            # il faut calculer correctement
            d_curr = self.evalMacroDamage(state_curr.u_nodes, dx_curr, d_curr)

        raise RuntimeError("No convergence of damage solver")


class CoupledSolver:
    dama_para = meca_para = None

    def __init__(self, MECA, DAMA):
        self.dama_para = DAMA
        self.meca_para = MECA

    def solve(self, timeList):
        print("Coupled solver")
        # init solver
        meca_solver = MecaSolver(self.meca_para)
        dama_solver = DamageSolver(self.dama_para)

        hho_meca = CA.HHO(CA.PhysicalProblem(self.meca_para["MODELE"], None))
        hho_dama = CA.HHO(CA.PhysicalProblem(self.dama_para["MODELE"], None))

        # Time stepper
        stepper = TimeStepper(timeList.getValues())

        # initial solution
        time_init = 0.0
        u_init = CA.FieldOnNodesReal(self.meca_para["MODELE"])
        dx_init = CA.FieldOnNodesReal(self.dama_para["MODELE"])
        d_init = hho_dama.projectOnHHOCellSpace(self.dama_para["SOURCE"])

        state_prev = CoupledState(u_init, dx_init, d_init, time_init)
        state_curr = CoupledState(u_init, dx_init, d_init, time_init)

        state_prev.projectOnLagrangeSpace(hho_meca, hho_dama)
        state_curr.projectOnLagrangeSpace(hho_meca, hho_dama)

        while not stepper.isFinished():
            state_curr.time = stepper.getCurrent()
            print("Time : %f" % stepper.getCurrent())
            for it in range(100):
                u_iter = state_curr.u.copy()
                dx_iter = state_curr.dx.copy()
                d_iter = state_curr.d.copy()

                # solve mechanical problem
                state_curr.u, state_curr.sief_elga, state_curr.vari_elga = meca_solver.solve(
                    state_prev, state_curr
                )

                # solve damaage problem
                state_curr.dx, state_curr.d = dama_solver.solve(state_curr, state_curr)

                resi_u = (state_curr.u - u_iter).norm("NORM_INFINITY")
                resi_dx = (state_curr.dx - dx_iter).norm("NORM_INFINITY")
                resi_d = (state_curr.d - d_iter).norm("NORM_INFINITY")

                print("Step %d with residual (%f, %f, %f)" % (it, resi_u, resi_dx, resi_d))
                if resi_u < 1e-8 and resi_dx < 1e-8 and resi_d < 1e-8:
                    break

                state_curr.projectOnLagrangeSpace(hho_meca, hho_dama)

            # finish step
            stepper.completed()

            state_prev = state_curr.copy()

        return state_curr
