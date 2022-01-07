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

from ...Cata.Language.SyntaxObjects import _F
from ...Commands import ASSE_VECTEUR, CALC_MATR_ELEM, CALCUL, DEFI_LIST_REEL
from ...Utilities import MPI, haveMPI, no_new_attributes, profile


class ElementaryComputation():
    """Compute elementary quantities like Dirichlet bc, Neumann loads, ..."""
    phys_state = phys_pb = None

    __setattr__ = no_new_attributes(object.__setattr__)

    def _getLoads(self):
        """Get list of loads (only CHARGE and not Dirichlet) - to delete quickly"""
        names = []
        listLoads = self.phys_pb.getListOfLoads()
        # diri = listLoads.getDirichletBCs()
        # for load in diri:
        #     names.append(load)

        mecar = listLoads.getMechanicalLoadsReal()
        for load in mecar:
            names.append(load)

        mecaf = listLoads.getMechanicalLoadsFunction()
        for load in mecaf:
            names.append(load)

        if haveMPI():
            mecar = listLoads.getParallelMechanicalLoadsReal()
            for load in mecar:
                names.append(load)

            mecaf = listLoads.getParallelMechanicalLoadsFunction()
            for load in mecaf:
                names.append(load)

        return names

    def setPhysicalProblem(self, phys_pb):
        """Assign the physical problem.

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
        """
        self.phys_pb = phys_pb

    def setPhysicalState(self, phys_state):
        """Assign the physical state.

        Arguments:
            phys_state (PhysicalState): Physical state
        """
        self.phys_state = phys_state

    @profile
    def computeInternalForces(self, time_prev, time_curr):
        """Compute internal forces, stress and internal state variables.

        Arguments:
            time_prev (float): previous time to evalute BC
            time_curr (float): current time to evalute BC

        Returns:
            tuple (int, FieldOnCells, FieldOnCells, FieldOnNodes):
            Tuple with 4 objects: exitcode, internal state variables (VARI_ELGA),
            Cauchy stress (SIEF_ELGA), vector of internal forces (`B^T \sigma`).
        """
        l_inst = DEFI_LIST_REEL(VALE=(time_prev, time_curr))
        behavProp = self.phys_pb.getBehaviourProperty()
        res = CALCUL(__use_namedtuple__=True,
                     OPTION=('FORC_INTE_ELEM',),
                     MODELE=self.phys_pb.getModel(),
                     CHAM_MATER=self.phys_pb.getMaterialField(),
                     CARA_ELEM = self.phys_pb.getElementaryCharacteristics(),
                     INCREMENT=_F(LIST_INST=l_inst, NUME_ORDRE=1),
                     EXCIT=self.phys_pb.allLoadsDict,
                     INCR_DEPL=self.phys_state.displ_incr,
                     DEPL=self.phys_state.displ,
                     SIGM=self.phys_state.stress,
                     VARI=self.phys_state.variP,
                     COMPORTEMENT=_F(COMPOR=behavProp.getBehaviourField(),
                                     MULT_COMP=behavProp.getMultipleBehaviourField(),
                                     CARCRI=behavProp.getConvergenceCriteria()),
                     INFO=1)

        # Check integration of behavior
        codret = MPI.COMM_WORLD.allreduce(
            max(res.CODE_RETOUR_INTE.EXTR_COMP("IRET",[]).valeurs),
            op=MPI.MAX)
        # forces internes assemblees
        assembly_forces = ASSE_VECTEUR(VECT_ELEM=res.FORC_INTE_ELEM,
                                       NUME_DDL=self.phys_pb.getDOFNumbering())

        return codret, res.VARI_ELGA, res.SIEF_ELGA, assembly_forces

    @profile
    def computeElasticStiffnessMatrix(self):
        """Compute elastic matrix (not assembled)

        Returns:
            tuple (int, ElementaryMatrixDisplacementReal): Tuple with exitcode,
            elementary elastic matrix.
        """
        elem_matr = CALC_MATR_ELEM(OPTION='RIGI_MECA',
                                   MODELE=self.phys_pb.getModel(),
                                   CHAM_MATER=self.phys_pb.getMaterialField(),
                                   CARA_ELEM = self.phys_pb.getElementaryCharacteristics(),
                                   CHARGE=self._getLoads())
        return 0, elem_matr

    @profile
    def computeTangentPredictionMatrix(self):
        """Compute tangent prediction matrix (not assembled).

        Returns:
            tuple (int, ElementaryMatrixDisplacementReal): Tuple with exitcode,
            elementary tangent prediction matrix.
        """
        time_prev = self.phys_state.time
        time_curr = self.phys_state.time + self.phys_state.time_step
        l_inst = DEFI_LIST_REEL(VALE=(time_prev, time_curr))
        behavProp = self.phys_pb.getBehaviourProperty()
        res = CALCUL(__use_namedtuple__=True,
                     OPTION="MATR_TANG_ELEM",
                     PHASE="PREDICTION",
                     MODELE=self.phys_pb.getModel(),
                     CHAM_MATER=self.phys_pb.getMaterialField(),
                     CARA_ELEM = self.phys_pb.getElementaryCharacteristics(),
                     INCREMENT=_F(LIST_INST=l_inst, NUME_ORDRE=1),
                     EXCIT=self.phys_pb.allLoadsDict,
                     INCR_DEPL=self.phys_state.displ_incr,
                     DEPL=self.phys_state.displ,
                     SIGM=self.phys_state.stress,
                     VARI=self.phys_state.variP,
                     COMPORTEMENT=_F(COMPOR=behavProp.getBehaviourField(),
                                     MULT_COMP=behavProp.getMultipleBehaviourField(),
                                     CARCRI=behavProp.getConvergenceCriteria()),
                     INFO=1)

        # Check integration of behavior
        codret = MPI.COMM_WORLD.allreduce(max(res.CODE_RETOUR_INTE.EXTR_COMP("IRET",[]).valeurs),
                                          op=MPI.MAX)
        return codret, res.MATR_TANG_ELEM

    @profile
    def computeTangentStiffnessMatrix(self):
        """Compute tangent matrix (not assembled).

        Returns:
            tuple (int, ElementaryMatrixDisplacementReal): Tuple with exitcode,
            elementary tangent matrix.
        """
        time_prev = self.phys_state.time
        time_curr = self.phys_state.time + self.phys_state.time_step
        l_inst = DEFI_LIST_REEL(VALE=(time_prev, time_curr))
        behavProp = self.phys_pb.getBehaviourProperty()
        res = CALCUL(__use_namedtuple__=True,
                     OPTION="MATR_TANG_ELEM",
                     MODELE=self.phys_pb.getModel(),
                     CHAM_MATER=self.phys_pb.getMaterialField(),
                     CARA_ELEM = self.phys_pb.getElementaryCharacteristics(),
                     INCREMENT=_F(LIST_INST=l_inst, NUME_ORDRE=1),
                     EXCIT=self.phys_pb.allLoadsDict,
                     INCR_DEPL=self.phys_state.displ_incr,
                     DEPL=self.phys_state.displ,
                     SIGM=self.phys_state.stress,
                     VARI=self.phys_state.variP,
                     COMPORTEMENT=_F(COMPOR=behavProp.getBehaviourField(),
                                     MULT_COMP=behavProp.getMultipleBehaviourField(),
                                     CARCRI=behavProp.getConvergenceCriteria()),
                     INFO=1)

        # Check integration of behavior
        codret = MPI.COMM_WORLD.allreduce(max(res.CODE_RETOUR_INTE.EXTR_COMP("IRET",[]).valeurs),
                                          op=MPI.MAX)
        return codret, res.MATR_TANG_ELEM
