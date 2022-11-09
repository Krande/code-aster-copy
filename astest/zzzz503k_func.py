# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

import code_aster
from code_aster.Commands import *


###################################################################
#
#   Solve coupling problem with HHO
#   u -> displacement, d -> macro-dommage
#
#   Continuous:
#   (sigma(u,d), grad v) = 0, Eps = eps(u) + alpha * d * Id
#   (A * grad d, grad p) + (H(u) * d, p)  - (H(u) * f, p) = 0
#
#   HHO:
#   sum_{T \in Th} (PK(huT, d_T), GkT(hvT))_T + stab(huT, hvT) = 0
#   sum_{T \in Th} (A * GkT(hdT), GkT(hpT))_T + stab(hdT, hpT) +
#      (H(uhT) * d_T, p_T) _T = ((uhT) * f, p_T)_T
#
###################################################################

import code_aster
from code_aster.Commands import *
from code_aster.Utilities import force_list


class MecaSolver:
    """ Solve the mecanical problem """

    model = mater = loads = None

    def __init__(self, param):
        self.model = param["MODELE"]
        self.mater = param["MATER"]
        self.loads = param["EXCIT"]

    def create_material(self, d_nodes):
        return AFFE_MATERIAU(MAILLAGE=self.model.getMesh(),
                             AFFE=_F(TOUT='OUI',
                                     MATER=self.mater,
                                     ),
                             AFFE_VARC=(_F(NOM_VARC='TEMP',
                                                    CHAM_GD=d_nodes,
                                                    VALE_REF=0.0),
                                        ),
                             )

    def solve(self, d_nodes):
        """ Solve mechanical problem """

        material = self.create_material(d_nodes)

        timeList = DEFI_LIST_REEL(VALE=(0.0, 1.0))

        resuMeca = STAT_NON_LINE(MODELE=self.model,
                                 CHAM_MATER=material,
                                 INCREMENT=_F(LIST_INST=timeList,),
                                 METHODE="NEWTON",
                                 EXCIT=self.loads)

        return resuMeca.getField("DEPL", 1), resuMeca.getField("SIEF_ELGA", 1)


class DommageSolver:
    """ Solve the dommage problem """

    model = mater = loads = source = None

    def __init__(self, param):
        self.model = param["MODELE"]
        self.mater = param["MATER"]
        self.loads = force_list(param["EXCIT"])
        self.source = param["SOURCE"]

    def create_material(self, u_nodes):
        therMate = DEFI_MATERIAU(THER=_F(
            LAMBDA=self.mater["LAMBDA"],
            RHO_CP=self.mater["RHO_CP"] * (1. + u_nodes.norm("NORM_2")),
        ),)

        return AFFE_MATERIAU(MAILLAGE=self.model.getMesh(),
                             AFFE=_F(TOUT='OUI',
                                     MATER=therMate,
                                     ),
                             )

    def solve(self, u_nodes):
        """ Solve dommage problem """

        material = self.create_material(u_nodes)

        # create PhysicalProblem
        phys_pb = code_aster.PhysicalProblem(self.model, material)

        for load in self.loads:
            if isinstance(load["CHARGE"], code_aster.ThermalDirichletBC):
                phys_pb.addDirichletBC(load["CHARGE"])

        # compute DOF numbering
        phys_pb.computeDOFNumbering()

        # create discrete computation
        disc_comp = code_aster.DiscreteComputation(phys_pb)
        hho = code_aster.HHO(phys_pb)

        # compute K = (lamda * GkT(hdT), GkT(hpT))_T + lambda * stab(hdT, hpT)
        matEK = disc_comp.getLinearStiffnessMatrix()
        matK = code_aster.AssemblyMatrixTemperatureReal(phys_pb)
        matK.addElementaryMatrix(matEK)
        matK.assemble()

        # compute M = (rho_cp(u) * d_T, p_T) _T
        matEM = disc_comp.getMassMatrix()
        matM = code_aster.AssemblyMatrixTemperatureReal(phys_pb)
        matM.addElementaryMatrix(matEM)
        matM.assemble()

        # project load
        f_hho = hho.projectOnHHOCellSpace(self.source)

        d_hho = hho.projectOnHHOSpace(0.0)

        diriBCs = disc_comp.getDirichletBC()

        # linear solver
        mySolver = code_aster.MumpsSolver()

        print("Newton solver:")
        for i in range(100):
            Resi = matK * d_hho + matM * (d_hho - f_hho)
            Jaco = matK + matM

            print("*Iter %d: residual %f" % (i, Resi.norm("NORM_2")))
            if Resi.norm("NORM_2") < 10e-8:
                return d_hho

            mySolver.factorize(Jaco)
            dd_hho = mySolver.solve(-Resi, diriBCs)
            d_hho += dd_hho

        raise RuntimeError("No convergence of dommage solver")


class CoupledSolver:

    domm_para = meca_para = None

    def __init__(self, MECA, DOMM):
        self.domm_para = DOMM
        self.meca_para = MECA

    def solve(self):

        print("Coupled solver")

        # initial solution
        d_init = code_aster.FieldOnNodesReal(self.domm_para["MODELE"])
        u_init = code_aster.FieldOnNodesReal(self.meca_para["MODELE"])

        d_curr = d_prev = d_init
        u_curr = u_prev = u_init

        # init solver
        meca_solver = MecaSolver(self.meca_para)
        domm_solver = DommageSolver(self.domm_para)

        hho_meca = code_aster.HHO(
            code_aster.PhysicalProblem(self.meca_para["MODELE"], None))
        hho_domm = code_aster.HHO(
            code_aster.PhysicalProblem(self.domm_para["MODELE"], None))

        for step in range(100):
            # solve mechanical problem
            u_curr, sief_curr = meca_solver.solve(
                hho_domm.projectOnLagrangeSpace(d_curr))

            # solve dommage problem
            d_curr = domm_solver.solve(hho_meca.projectOnLagrangeSpace(u_curr))

            resi_meca = (u_curr - u_prev).norm("NORM_INFINITY")
            resi_domm = (d_curr - d_prev).norm("NORM_INFINITY")

            print("Step %d with residual (%f, %f)" %
                  (step, resi_meca, resi_domm))
            if resi_meca < 1e-6 and resi_domm < 1e-6:
                return u_curr, sief_curr, d_curr

            u_prev = u_curr
            d_prev = d_curr

        raise RuntimeError("No convergence of coupled solver")
