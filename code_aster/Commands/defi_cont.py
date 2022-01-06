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

from ..Objects import ContactNew, ContactZone, ContactParameter, \
    FrictionParameter, PairingParameter, ContactAlgo, ContactType, ContactVariant, \
    FrictionAlgo, FrictionType, PairingAlgo, InitialState
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

        # print("ARGS: ", keywords, flush=True)
        model = keywords["MODELE"]
        verbosity = keywords["INFO"]
        
        # usefull dict
        _algo_cont = {"LAGRANGIEN": ContactAlgo.Lagrangian, "NITSCHE": ContactAlgo.Nitsche,
                      "PENALISATION": ContactAlgo.Penalization}
        _type_cont = {"UNILATERAL": ContactType.Unilateral,
                      "BILATERAL": ContactType.Bilateral, "COLLE": ContactType.Stick}
        _vari_cont = {"RAPIDE": ContactVariant.Rapide,
                      "ROBUSTE": ContactVariant.Robust}
        _algo_frot = {"LAGRANGIEN": FrictionAlgo.Lagrangian, "NITSCHE": FrictionAlgo.Nitsche,
                      "PENALISATION": FrictionAlgo.Penalization}
        _type_frot = {"TRESCA": FrictionType.Tresca, "SANS": FrictionType.Without,
                      "COULOMB": FrictionType.Coulomb, "COLLE": FrictionType.Stick}
        _algo_pair = {"MORTAR": PairingAlgo.Mortar}
        _init_cont = {"INTERPENETRE": InitialState.Interpenetrated, "NON": InitialState.No,
                      "OUI": InitialState.Yes}              

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
            if ( zone.get("SANS_GROUP_MA")) != None:  
                contZone.setExcludedSlaveGroupOfCells(zone["SANS_GROUP_MA"])    

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
                fricParam.hasFriction(keywords["FROTTEMENT"] == "OUI")
                fricParam.setAlgorithm(_algo_frot[zone["ALGO_FROT"]])
                fricParam.setType(_type_frot[zone["TYPE_FROT"]])
                if fricParam.getType() == FrictionType.Tresca:
                    fricParam.setTresca(zone["TRESCA"])
                elif fricParam.getType() == FrictionType.Coulomb:
                    fricParam.setCoulomb(zone["COULOMB"])
                fricParam.setCoefficient(zone["COEF_FROT"])
                contZone.setFrictionParameter(fricParam)

            # pairing parameters
            pairParam = PairingParameter()
            contZone.setPairingParameter(pairParam)
            pairParam.setAlgorithm(_algo_pair[zone["APPARIEMENT"]])
            pairParam.setPairingDistance(zone["DIST_APPA"]) 
            pairParam.setInitialState(_init_cont[zone["CONTACT_INIT"]])
            pairParam.hasBeamDistance(zone["DIST_POUTRE"] == "OUI")
            pairParam.hasShellDistance(zone["DIST_COQUE"] == "OUI")
            if (pairParam.hasBeamDistance() or pairParam.hasShellDistance() ):
                pairParam.setElementaryCharacteristics(zone["CARA_ELEM"]) 
            if ( zone.get("SEUIL_INIT")) != None:  
                pairParam.setThreshold(zone["SEUIL_INIT"])    
            if ( zone.get("DIST_SUPP") ) != None:  
                pairParam.setDistanceFunction(zone["DIST_SUPP"])            

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
