# coding: utf-8

# Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

from ..Objects import (MechanicalLoadFunction, ParallelMechanicalLoadFunction,
                       ConnectionMesh)
from ..Supervis import ExecuteCommand
from ..Utilities import deprecate, force_list
from .affe_modele import AFFE_MODELE
from .affe_char_meca import MechanicalLoadDefinition


class MechanicalLoadFunctionDefinition(ExecuteCommand):
    """Command that defines :class:`~code_aster.Objects.MechanicalLoadFunc`.
    """
    command_name = "AFFE_CHAR_MECA_F"

    def compat_syntax(self, keywords):
        """Compatibility support, common with AFFE_CHAR_MECA."""
        return MechanicalLoadDefinition.compat_syntax(keywords)

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        model = keywords["MODELE"]
        l_neum = MechanicalLoadDefinition()._hasNeumannLoadings(keywords)
        l_diri = MechanicalLoadDefinition()._hasDirichletLoadings(keywords)
        if not model.getMesh().isParallel():
            self._result = MechanicalLoadFunction(model)
        else :
            if l_neum:
                if l_diri:
                    raise TypeError("Not allowed to mix up Dirichlet and Neumann loadings in the same parallel AFFE_CHAR_MECA_F")
                else:
                    self._result = MechanicalLoadFunction(model)
            if MechanicalLoadDefinition()._hasOnlyDDL_IMPO(keywords):
                self._result = MechanicalLoadFunction(model)


    def exec_(self, keywords):
        """Override default _exec in case of parallel load
        """
        if isinstance(self._result, MechanicalLoadFunction):
            super(MechanicalLoadFunctionDefinition, self).exec_(keywords)
        else:
            model = keywords.pop("MODELE")
            nodeGroups, cellGroups = MechanicalLoadDefinition()._getGroups(keywords)
            connectionMesh = ConnectionMesh(model.getMesh(), nodeGroups, cellGroups)

            if connectionMesh.getDimension()==3:
                partialModel = AFFE_MODELE(MAILLAGE=connectionMesh,
                                       AFFE=(_F(TOUT='OUI',
                                               PHENOMENE='MECANIQUE',
                                               MODELISATION="3D",),
                                             _F(TOUT='OUI',
                                               PHENOMENE='MECANIQUE',
                                               MODELISATION='DIS_T',),
                                               ),
                                       VERI_JACOBIEN='NON',
                                       DISTRIBUTION=_F(METHODE='CENTRALISE',),)
            else:
                partialModel = AFFE_MODELE(MAILLAGE=connectionMesh,
                                       AFFE=(_F(TOUT='OUI',
                                               PHENOMENE='MECANIQUE',
                                               MODELISATION="D_PLAN",),),
                                       VERI_JACOBIEN='NON',
                                       DISTRIBUTION=_F(METHODE='CENTRALISE',),)

            partialModel.getFiniteElementDescriptor().transferDofDescriptorFrom(model.getFiniteElementDescriptor())
            keywords["MODELE"] = partialModel
            partialMechanicalLoad = AFFE_CHAR_MECA_F(**keywords)
            keywords["MODELE"] = model
            self._result = ParallelMechanicalLoadFunction(partialMechanicalLoad, model)


AFFE_CHAR_MECA_F = MechanicalLoadFunctionDefinition.run
