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

from ..Objects import ContactNew, ContactZone, ContactParameter, \
    FrictionParameter, PairingParameter, ContactAlgo, ContactType, ContactVariant
from ..Supervis import ExecuteCommand

# C'est la nouvelle commande à remplir qui sera en python et c++


class NewContactAssignment(ExecuteCommand):
    """Command that creates the :class:`~code_aster.Objects.ContactNew` by
    assigning finite elements"""
    command_name = "DEFI_CONT"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = ContactNew(keywords["MODELE"])

    def exec_(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """

        print("ARGS: ", keywords, flush=True)
        model = keywords["MODELE"]
        verbosity = keywords["INFO"]

        # usefull dict
        _algo_cont = {"LAGRANGIEN": ContactAlgo.Lagrangian, "NITSCHE": ContactAlgo.Nitsche,
                      "PENALISATION": ContactAlgo.Penalization}
        _type_cont = {"UNILATERAL": ContactType.Unilateral,
                      "BILATERAL": ContactType.Bilateral, "COLLE": ContactType.Stick}
        _vari_cont = {"RAPIDE": ContactVariant.Rapide,
                      "ROBUSTE": ContactVariant.Robust}

        # add global informations
        self._result.setVerbosity(verbosity)
        self._result.hasFriction(keywords["FROTTEMENT"] == "OUI")
        self._result.hasSmoothing(keywords["LISSAGE"] == "OUI")

        # add infomations for each ZONE
        list_zones = keywords["ZONE"]
        for zone in list_zones:
            contZone = ContactZone(model)
            contZone.checkNormals(zone["VERI_NORM"] == "OUI")
            contZone.setVerbosity(verbosity)
            contZone.setSlaveGroupOfCells(zone["GROUP_MA_ESCL"])
            contZone.setMasterGroupOfCells(zone["GROUP_MA_MAIT"])

            # contact parameters
            contParam = ContactParameter()
            contParam.setAlgorithm(_algo_cont[zone["ALGO_CONT"]])
            if _algo_cont[zone["ALGO_CONT"]] in (ContactAlgo.Nitsche,
                                                 ContactAlgo.Lagrangian):
                contParam.setVariant(_vari_cont[zone["VARIANTE"]])
            else:
                contParam.setVariant(ContactVariant.Empty)
            contParam.setType(_type_cont[zone["TYPE_CONT"]])
            contParam.setCoefficient(zone["COEF_CONT"])
            contZone.setContactParameter(contParam)

            # friction parameters
            if self._result.hasFriction():
                fricParam = FrictionParameter()

                contZone.setFrictionParameter(fricParam)

            # pairing parameters
            pairParam = PairingParameter()

            contZone.setPairingParameter(pairParam)

            # build then append
            contZone.build()
            self._result.appendContactZone(contZone)

        # build result
        self._result.build()

    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.

        Arguments:
            keywords (dict): User's keywords.
        """

        # Add no dependencies since everything is done in c++ directly


DEFI_CONT = NewContactAssignment.run
