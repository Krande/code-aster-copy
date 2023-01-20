# coding: utf-8

# Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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
    LoadResult,
    ThermalResult,
    FieldOnNodesComplex,
    FieldOnNodesReal,
    FieldOnNodesChar8,
    FieldOnCellsComplex,
    FieldOnCellsReal,
    FieldOnCellsLong,
    ElasticFourierResult,
    ThermalFourierResult,
    FullHarmonicResult,
    FullTransientResult,
    ExternalStateVariablesResult,
    ElasticResult,
    ModeResultComplex,
    ModeResult,
    MultipleElasticResult,
    NonLinearResult,
)
from ..Supervis import ExecuteCommand
from ..Utilities import force_list


class ResultCreator(ExecuteCommand):
    """Command that creates evolutive results."""

    command_name = "CREA_RESU"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        if "reuse" in keywords:
            self._result = keywords["reuse"]
        else:
            typ = keywords["TYPE_RESU"]
            if typ == "EVOL_CHAR":
                self._result = LoadResult()
            elif typ == "EVOL_THER":
                self._result = ThermalResult()
            elif typ == "EVOL_NOLI":
                self._result = NonLinearResult()
            elif typ == "EVOL_ELAS":
                self._result = ElasticResult()
            elif typ == "EVOL_VARC":
                self._result = ExternalStateVariablesResult()
            elif typ == "FOURIER_ELAS":
                self._result = ElasticFourierResult()
            elif typ == "FOURIER_THER":
                self._result = ElasticFourierResult()
            elif typ == "MULT_ELAS":
                self._result = MultipleElasticResult()
            elif typ == "MODE_MECA":
                self._result = ModeResult()
            elif typ == "MODE_MECA_C":
                self._result = ModeResultComplex()
            elif typ == "DYNA_TRANS":
                self._result = FullTransientResult()
            elif typ == "DYNA_HARMO":
                self._result = FullHarmonicResult()
            else:
                raise NotImplementedError("Type of result {0!r} not yet " "implemented".format(typ))

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        fkw = force_list(keywords.get("AFFE", []))

        for occ in fkw:
            if occ.get("MODELE"):
                self._result.setModel(occ["MODELE"])
                break
            chamGd = occ.get("CHAM_GD")
            if chamGd is not None:
                if isinstance(chamGd, (FieldOnNodesReal, FieldOnNodesComplex, FieldOnNodesChar8)):
                    mesh = chamGd.getMesh()
                    if mesh is not None:
                        self._result.setMesh(mesh)
                        break

        # find ligrel and prof_chno
        feds = []
        fnds = []
        for occ in fkw:
            chamGd = occ.get("CHAM_GD")
            if chamGd is not None:
                if isinstance(chamGd, (FieldOnNodesReal, FieldOnNodesComplex, FieldOnNodesChar8)):
                    fnd = chamGd.getDescription()
                    if fnd is not None:
                        fnds.append(fnd)
                elif isinstance(chamGd, (FieldOnCellsReal, FieldOnCellsComplex, FieldOnCellsLong)):
                    fed = chamGd.getDescription()
                    if fed is not None:
                        feds.append(fed)

        if keywords.get("MATR_RIGI"):
            self._result.setModel(keywords["MATR_RIGI"].getDOFNumbering().getModel())
            self._result.setMesh(keywords["MATR_RIGI"].getMesh())
            self._result.setDOFNumbering(keywords["MATR_RIGI"].getDOFNumbering())
        elif keywords.get("MATR_MASS"):
            self._result.setModel(keywords["MATR_MASS"].getDOFNumbering().getModel())
            self._result.setMesh(keywords["MATR_MASS"].getMesh())
            self._result.setDOFNumbering(keywords["MATR_MASS"].getDOFNumbering())

        if keywords.get("ECLA_PG"):
            self._result.setModel(keywords["ECLA_PG"]["MODELE_INIT"])
        if keywords.get("CONV_CHAR"):
            matr_rigi = keywords["CONV_CHAR"]["MATR_RIGI"]
            self._result.setModel(matr_rigi.getDOFNumbering().getModel())
            self._result.setMesh(matr_rigi.getMesh())
            dofNum = matr_rigi.getDOFNumbering()
            if dofNum:
                if keywords["TYPE_RESU"] == "DYNA_TRANS":
                    self._result.setDOFNumbering(dofNum)
                else:
                    fnds.append(dofNum.getDescription())
        if keywords.get("CONV_RESU"):
            self._result.setModel(keywords["CONV_RESU"]["RESU_INIT"].getModel())
            matr_rigi = keywords["CONV_RESU"].get("MATR_RIGI")
            if matr_rigi is not None:
                dofNum = matr_rigi.getDOFNumbering()
            else:
                dofNum = keywords["CONV_RESU"]["NUME_DDL"]
            if dofNum:
                if keywords["TYPE_RESU"] == "DYNA_TRANS":
                    self._result.setDOFNumbering(dofNum)
                else:
                    fnds.append(dofNum.getDescription())
        if keywords.get("KUCV"):
            self._result.setModel(keywords["KUCV"]["RESU_INIT"].getModel())
            dofNum = keywords["KUCV"]["MATR_AMOR"].getDOFNumbering()
            if dofNum:
                fnds.append(dofNum.getDescription())
        if keywords.get("PROL_RTZ"):
            self._result.setMesh(keywords["PROL_RTZ"]["MAILLAGE_FINAL"])

        if not fkw:
            fkw = keywords.get("ASSE")
        if not fkw:
            fkw = keywords.get("PREP_VRC2")
        if not fkw:
            fkw = keywords.get("PREP_VRC1")

        if fkw:
            chamMater = fkw[0].get("CHAM_MATER")
            if chamMater is not None:
                self._result.setMaterialField(chamMater)

            model = fkw[0].get("MODELE")
            chamGd = fkw[0].get("CHAM_GD")
            result = fkw[0].get("RESULTAT")

            if model is not None:
                self._result.setModel(model)
            elif result is not None:
                model = result.getModel()
                if model is not None:
                    self._result.setModel(model)

                mesh = result.getMesh()
                if mesh is not None:
                    self._result.setMesh(mesh)
            elif chamGd is not None:
                mesh = chamGd.getMesh()
                if mesh is not None:
                    self._result.setMesh(mesh)

        self._result.build(feds, fnds)

    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.

        Arguments:
            keywords (dict): User's keywords.
        """
        super().add_dependencies(keywords)
        self.remove_dependencies(keywords, "AFFE", "CHAM_GD")
        self.remove_dependencies(keywords, "ASSE", "RESULTAT")
        self.remove_dependencies(keywords, "ECLA_PG", "RESU_INIT")
        if keywords["OPERATION"] == "PERM_CHAM":
            self.remove_dependencies(keywords, "RESU_INIT")
            self.remove_dependencies(keywords, "RESU_FINAL")
        self.remove_dependencies(keywords, "PROL_RTZ", "TABLE")
        self.remove_dependencies(keywords, "PREP_VRC1", "CHAM_GD")
        self.remove_dependencies(keywords, "PREP_VRC2", "ELAS_THER")
        self.remove_dependencies(keywords, "KUCV", "RESU_INIT")
        self.remove_dependencies(keywords, "CONV_RESU", "RESU_INIT")


CREA_RESU = ResultCreator.run
