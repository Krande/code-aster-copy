# coding: utf-8

# Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

from ..Objects import MechanicalLoadFunction, ParallelMechanicalLoadFunction, ConnectionMesh, Model
from ..Supervis import ExecuteCommand
from .affe_char_meca import MechanicalLoadDefinition, _getGroups
from ..Cata.Commands import DEFI_CONSTANTE


class MechanicalLoadFunctionDefinition(ExecuteCommand):
    """Command that defines :class:`~code_aster.Objects.MechanicalLoadFunc`."""

    command_name = "AFFE_CHAR_MECA_F"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        model = keywords["MODELE"]
        l_neum = MechanicalLoadDefinition._hasNeumannLoadings(keywords)
        l_diri = MechanicalLoadDefinition._hasDirichletLoadings(keywords)
        if not model.getMesh().isParallel():
            self._result = MechanicalLoadFunction(model)
        else:
            if l_neum:
                if l_diri:
                    raise TypeError(
                        "Not allowed to mix up Dirichlet and Neumann loadings in the same parallel AFFE_CHAR_MECA_F"
                    )
                else:
                    self._result = MechanicalLoadFunction(model)
            if MechanicalLoadDefinition._hasOnlyDDL_IMPO(keywords):
                self._result = MechanicalLoadFunction(model)

    def exec_(self, keywords):
        """Override default _exec in case of parallel load"""
        if isinstance(self._result, MechanicalLoadFunction):
            super(MechanicalLoadFunctionDefinition, self).exec_(keywords)
        else:
            model = keywords.pop("MODELE")
            nodeGroups, cellGroups = _getGroups(self._cata, keywords)
            connectionMesh = ConnectionMesh(model.getMesh(), nodeGroups, cellGroups)

            connectionModel = Model(connectionMesh)
            connectionModel.setFrom(model)

            keywords["MODELE"] = connectionModel
            partialMechanicalLoad = AFFE_CHAR_MECA_F(**keywords)
            keywords["MODELE"] = model
            self._result = ParallelMechanicalLoadFunction(partialMechanicalLoad, model)

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """

        self._result.build()

    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.

        Arguments:
            keywords (dict): User's keywords.
        """
        super().add_dependencies(keywords)
        self.remove_dependencies(keywords, "MODELE")

    def adapt_syntax(self, keywords):
        """Hook to adapt syntax *after* syntax checking.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        if "PRE_EPSI" in keywords:
            # Convert the 3 real components of VECT_N to constant functions
            l_dic_kws = keywords.get("PRE_EPSI")
            if (
                type(l_dic_kws) == tuple or type(l_dic_kws) == list
            ):  # il y a plus d'une occurrence de MASSIF
                for dic in l_dic_kws:
                    if "VECT_N" in dic.keys():
                        convert_real_to_constant_function(dic, "VECT_N")
                        print("etienne keywords", keywords)
            # check syntax after changing the keywords
            self.check_syntax(keywords)


def convert_real_to_constant_function(mcfact, name):
    """Replace a list of real by constant functions (separed parameters).

    Example: VECT_N=(x, y, z) is replaced by
             VECT_N1 =DEFI_CONSTANTE(VALE=x),
             VECT_N2 =DEFI_CONSTANTE(VALE=y),
             VECT_N3 =DEFI_CONSTANTE(VALE=z),

    Args:
        mcfact (dict): factor keyword, changed in place.
        name (str): keyword containing the list.
    """
    # values = mcfact.pop(name, None)
    if name in mcfact:
        values = mcfact[name]
    else:
        return
    # if not values:
    #     return
    func = [None] * len(values)
    for i, value in enumerate(values):
        func[i] = DEFI_CONSTANTE(VALE=value)
        mcfact[f"{name}{i + 1}"] = func[i]


AFFE_CHAR_MECA_F = MechanicalLoadFunctionDefinition.run
