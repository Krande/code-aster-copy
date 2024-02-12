# coding: utf-8

# Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.

# person_in_charge: nicolas.sellenet@edf.fr

from ..Objects import (
    AssemblyMatrixPressureComplex,
    FullHarmonicAcousticResult,
    FullHarmonicResult,
    FullTransientResult,
    HarmoGeneralizedResult,
    TransientGeneralizedResult,
    Function,
)
from ..Supervis import ExecuteCommand

from ..Helpers import check_dis_choc_elas


class VibrationDynamics(ExecuteCommand):
    """Command to solve linear vibration dynamics problem, on physical or modal bases,
    for harmonic or transient analysis.
    """

    command_name = "DYNA_VIBRA"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        if keywords.get("reuse") != None:
            self._result = keywords["reuse"]
        else:
            base = keywords["BASE_CALCUL"]
            typ = keywords["TYPE_CALCUL"]
            matrRigi = keywords["MATR_RIGI"]
            if base == "PHYS":
                if typ == "TRAN":
                    self._result = FullTransientResult()
                    return
                if isinstance(matrRigi, AssemblyMatrixPressureComplex):
                    self._result = FullHarmonicAcousticResult()
                    return
                self._result = FullHarmonicResult()
            else:
                if typ == "TRAN":
                    self._result = TransientGeneralizedResult()
                else:
                    self._result = HarmoGeneralizedResult()

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        if keywords["BASE_CALCUL"] == "PHYS":
            massMatrix = keywords["MATR_MASS"]
            dofNum = massMatrix.getDOFNumbering()
            self._result.setModel(dofNum.getModel())
            self._result.setDOFNumbering(dofNum)
            for i in dofNum.getFiniteElementDescriptors():
                self._result.addFiniteElementDescriptor(i)
            self._result.setDOFNumbering(dofNum)
            self._result.setModel(dofNum.getModel())
            mesh = massMatrix.getMesh()
            if mesh is not None:
                self._result.setMesh(mesh)
            self._result.build()
        if keywords["BASE_CALCUL"] == "GENE":
            stiffnessMatrix = keywords["MATR_RIGI"]
            dofGeneNum = stiffnessMatrix.getGeneralizedDOFNumbering()
            if isinstance(self._result, (HarmoGeneralizedResult, TransientGeneralizedResult)):
                self._result.setGeneralizedDOFNumbering(dofGeneNum)
            else:
                raise Exception("Unknown result type")
            if keywords["TYPE_CALCUL"] == "TRAN":
                self._result.build()

    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.

        Arguments:
            keywords (dict): User's keywords.
        """
        self._result.resetDependencies()
        for key in ("MATR_MASS", "MATR_RIGI", "MATR_AMOR"):
            if keywords.get(key):
                self._result.addDependency(keywords[key])

    def adapt_syntax(self, keywords):
        """Adapt syntax *after* syntax checking.

        Arguments:
            keywords (dict): User's keywords. Changed in place
        """
        LstComportement = keywords.get("COMPORTEMENT", None)
        if not LstComportement is None:
            for UnComportement in LstComportement:
                if UnComportement.get("RELATION") == "CHOC_ELAS_TRAC":
                    tmp = check_dis_choc_elas(UnComportement)


DYNA_VIBRA = VibrationDynamics.run
