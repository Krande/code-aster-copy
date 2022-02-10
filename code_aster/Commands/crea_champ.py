# coding: utf-8

# Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

from ..Objects import (FieldOnCellsReal, FieldOnCellsLong, FieldOnCellsComplex,
                       FieldOnNodesReal, FieldOnNodesComplex,
                       FullResult, ModeResult,
                       ConstantFieldOnCellsReal)
from ..Supervis import ExecuteCommand
from ..Utilities import force_list


class FieldCreator(ExecuteCommand):
    """Command that creates fields that may be
    :class:`~code_aster.Objects.FieldOnNodesReal` or
    :class:`~code_aster.Objects.ConstantFieldOnCellsReal`."""
    command_name = "CREA_CHAMP"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        location = keywords["TYPE_CHAM"][:5]
        typ = keywords["TYPE_CHAM"][10:]

        mesh = keywords.get("MAILLAGE")
        model = keywords.get("MODELE")
        caraElem = keywords.get("CARA_ELEM")
        charge = keywords.get("CHARGE")
        resultat = keywords.get("RESULTAT")
        numeDdl = keywords.get("NUME_DDL")
        chamgd = keywords.get("CHAM_GD")
        fiss = keywords.get("FISSURE")
        if mesh is None:
            if model is not None:
                mesh = model.getMesh()
            elif caraElem is not None:
                mesh = caraElem.getModel().getMesh()
            elif charge is not None:
                mesh = charge.getModel().getMesh()
            elif resultat is not None:
                mesh = resultat.getMesh()
            elif chamgd is not None:
                mesh = chamgd.getMesh()
            elif fiss is not None:
                mesh = fiss.getMesh()

        if location == "CART_":
            if mesh is None:
                raise NotImplementedError("Must have Mesh, Model or ElementaryCharacteristics")
            self._result = ConstantFieldOnCellsReal(mesh)
        elif location == "NOEU_":
            if typ == "C":
                self._result = FieldOnNodesComplex()
            else:
                self._result = FieldOnNodesReal()
            if mesh is not None:
                self._result.setMesh(mesh)
            if numeDdl is not None:
                self._result.setDescription(numeDdl.getDescription())
        else:
            if typ == "R":
                self._result = FieldOnCellsReal()
            elif typ == "I":
                self._result = FieldOnCellsLong()
            elif typ == "C":
                self._result = FieldOnCellsComplex()
            else:
                raise NotImplementedError("Output for CREA_CHAMP not defined")

        if location[:2] == "EL":
            chamF = keywords.get("CHAM_F")
            if model is not None:
                self._result.setModel(model)
            elif resultat is not None:
                if isinstance(resultat, (FullResult, ModeResult)):
                    try:
                        dofNum = resultat.getDOFNumbering()
                        if dofNum is not None:
                            self._result.setDescription(dofNum.getFiniteElementDescriptors()[0])
                    except:
                        pass
                    if self._result.getModel() is None:
                        try:
                            dofNum = resultat.getDOFNumbering()
                            if (dofNum is not None) and (dofNum.getModel() is not None):
                                self._result.setModel(dofNum.getModel())
                        except:
                            pass
                if resultat.getModel() is not None:
                    self._result.setModel(resultat.getModel())
            elif caraElem is not None:
                self._result.setModel(caraElem.getModel())
            elif chamF is not None:
                self._result.setModel(chamF.getModel())

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """

        if isinstance(self._result, FieldOnNodesReal):
            if not self._result.getDescription():
                for comb in force_list(keywords.get("COMB", [])):
                    desc = comb['CHAM_GD'].getDescription()
                    if desc:
                        self._result.setDescription(desc)
                        break
            if not self._result.getMesh():
                for comb in force_list(keywords.get("COMB", [])):
                    mesh = comb['CHAM_GD'].getMesh()
                    if mesh:
                        self._result.setMesh(mesh)
                        break

        self._result.build()

    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.

        Arguments:
            keywords (dict): User's keywords.
        """
        super().add_dependencies(keywords)
        self.remove_dependencies(keywords, "ASSE", "CHAM_GD")
        self.remove_dependencies(keywords, "COMB", "CHAM_GD")

        if keywords["OPERATION"] in ("ASSE_DEPL", "R2C", "C2R", "DISC"):
            self.remove_dependencies(keywords, "CHAM_GD")

        self.remove_dependencies(keywords, "EVAL", ("CHAM_F", "CHAM_PARA"))

        if keywords["OPERATION"] == "EXTR":
            # depends on "result".LIGREL
            # self.remove_dependencies(keywords, "RESULTAT")
            self.remove_dependencies(keywords, "FISSURE")
            self.remove_dependencies(keywords, "TABLE")
            self.remove_dependencies(keywords, "CARA_ELEM")
            self.remove_dependencies(keywords, "CHARGE")


CREA_CHAMP = FieldCreator.run
