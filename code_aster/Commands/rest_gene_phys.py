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

# person_in_charge: nicolas.sellenet@edf.fr

from ..Objects import (FullHarmonicResult,
                       FullTransientResult, GeneralizedModeResult,
                       HarmoGeneralizedResult,
                       ModeResult,
                       TransientGeneralizedResult)
from ..Supervis import ExecuteCommand


class RestGenePhys(ExecuteCommand):
    """Command REST_GENE_PHYS
    """
    command_name = "REST_GENE_PHYS"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        if isinstance(keywords["RESU_GENE"], TransientGeneralizedResult):
            self._result = FullTransientResult()
        elif isinstance(keywords["RESU_GENE"], HarmoGeneralizedResult):
            self._result = FullHarmonicResult()
        elif isinstance(keywords["RESU_GENE"], GeneralizedModeResult):
            self._result = ModeResult()
        else:
            raise Exception("Unknown result type")

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """

        feds = []
        fnds = []

        if "NUME_DDL" in keywords:
            self._result.setDOFNumbering(keywords["NUME_DDL"])

        if "MODE_MECA" in keywords:
            feds += keywords["MODE_MECA"].getFiniteElementDescriptors()
            fnds += keywords["MODE_MECA"].getFieldOnNodesDescriptions()

        resu_gene = keywords["RESU_GENE"]
        if isinstance(resu_gene, (TransientGeneralizedResult,
                                  HarmoGeneralizedResult)):
            dofNum = resu_gene.getDOFNumbering()

            if dofNum is not None:
                self._result.setDOFNumbering(dofNum)
                modele = dofNum.getModel()
                if modele is not None:
                    self._result.setModel(modele)
            else:
                geneDofNum = resu_gene.getGeneralizedDOFNumbering()
                mesh = None
                if geneDofNum is not None:
                    basis = geneDofNum.getModalBasis()
                    if basis is not None:
                        dofNum = basis.getDOFNumbering()
                        if mesh is None:
                            mesh = basis.getMesh()
                        if dofNum is not None:
                            self._result.setDOFNumbering(dofNum)
                            modele = dofNum.getModel()
                            if modele is not None:
                                self._result.setModel(modele)
                            elif mesh is None:
                                mesh = dofNum.getMesh()

                if mesh is None:
                    if "MODE_MECA" in keywords:
                        mesh = keywords["MODE_MECA"].getMesh()

                if mesh is not None:
                    self._result.setMesh(mesh)

        elif isinstance(resu_gene, GeneralizedModeResult):
            self._result.setMesh(resu_gene.getMesh())
            matrRigi = resu_gene.getStiffnessMatrix()
            if matrRigi is not None:
                modalBasis = matrRigi.getModalBasis()
                if modalBasis is None:
                    dofNum = matrRigi.getGeneralizedDOFNumbering()
                    if dofNum:
                        modalBasis = dofNum.getModalBasis()
                if modalBasis is not None:
                    dofNum = modalBasis.getDOFNumbering()
                    if dofNum is not None:
                        self._result.setDOFNumbering(dofNum)
            else:
                # resultat issue de proj_mesu_modal
                dofNum = resu_gene.getDOFNumbering()
                if dofNum is not None:
                    self._result.setDOFNumbering(dofNum)
        else:
            raise Exception("Unknown result type")

        if self._result.getMesh() is None:
            for fed in feds:
                self._result.setMesh(fed.getMesh())

        if self._result.getMesh():
            self._result.build(feds, fnds)


REST_GENE_PHYS = RestGenePhys.run
