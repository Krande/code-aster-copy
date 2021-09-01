# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`Mesh` --- Assignment of mesh
************************************************************************
"""

import aster

from ..Commands import MACR_ADAP_MAIL
from ..Objects import Mesh
from ..Supervis import CO
from ..Utilities import injector
from ..Utilities.MedUtils.MEDConverter import convertMesh2MedCoupling


@injector(Mesh)
class ExtendedMesh:
    cata_sdj = "SD.sd_maillage.sd_maillage"

    def LIST_GROUP_NO(self):
        """Retourne la liste des groupes de noeuds sous la forme :
        [ (gno1, nb noeuds  gno1), ...]"""
        dic_gpno = self.sdj.GROUPENO.get()
        if dic_gpno is None:
            return []
        return [(gpno.strip(), len(dic_gpno[gpno])) for gpno in dic_gpno]

    def LIST_GROUP_MA(self):
        """Retourne la liste des groupes de mailles sous la forme :
        [ (gma1, nb mailles gma1, dime max des mailles gma1), ...]"""
        ltyma = aster.getvectjev("&CATA.TM.NOMTM")
        catama = aster.getcolljev("&CATA.TM.TMDIM")
        dic_gpma = self.sdj.GROUPEMA.get()
        if dic_gpma is None:
            return []
        dimama = [catama[ltyma[ma - 1].ljust(24)][0] for ma in self.sdj.TYPMAIL.get()]
        ngpma = []
        for grp in list(dic_gpma.keys()):
            dim = max([dimama[ma - 1] for ma in dic_gpma[grp]])
            ngpma.append((grp.strip(), len(dic_gpma[grp]), dim))
        return ngpma

    def refine(self, ntimes=1):
        """Refine the mesh uniformly. Each edge is split in two.

        Arguments:
            ntimes [int] : the number of times the mesh is to be refined.

        Returns:
            mesh: the refined mesh.
        """
        newMesh = self
        for i in range(ntimes):
            resu = MACR_ADAP_MAIL(__use_namedtuple__=True,
                                  ADAPTATION='RAFFINEMENT_UNIFORME',
                                  MAILLAGE_N=newMesh,
                                  MAILLAGE_NP1=CO('newMesh'))
            newMesh = resu.newMesh
        return newMesh

    def createMedCouplingMesh(self):
        """Returns the MEDCoupling unstructured mesh associated to the current mesh.

        Returns:
            mesh: The MEDCoupling unstructured mesh associated to the current mesh.
        """

        return convertMesh2MedCoupling(self)
