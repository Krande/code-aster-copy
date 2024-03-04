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

from math import sqrt

from ..Messages import UTMESS
from ..Objects import Function
from ..Supervis import ExecuteMacro
from ..Utilities import force_list


class MaterialDefinition(ExecuteMacro):
    """Command that defines a Material."""

    command_name = "DEFI_MATERIAU"

    def adapt_syntax(self, keywords):
        """Adapt syntax to replace TABLE content in keywords.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        adapt_elas(keywords)
        if "TABLE" not in keywords:
            return

        self.print_syntax(keywords)
        cmd_name = "DEFI_MATERIAU <TABLE>"
        cmd_kws = self._cata.keywords

        table_kws = keywords.get("TABLE")
        cmd_kws["TABLE"].addDefaultKeywords(table_kws)
        table_values = table_kws["TABLE"].EXTR_TABLE().values()

        varc_name = table_kws["NOM_PARA"]

        # Vérification du contenu de la table
        if not varc_name in table_values.keys():
            errmsg = "Parameter '{0}' is missing.".format(varc_name)
            UTMESS("F", "SUPERVIS_4", valk=(cmd_name, errmsg))

        for mater in table_kws["COMPOR"]:
            if mater in keywords or mater.replace("_FO", "") in keywords:
                errmsg = "Material '{0}' cannot be repeated outside TABLE.".format(
                    mater.replace("_FO", "")
                )
                UTMESS("F", "SUPERVIS_4", valk=(cmd_name, errmsg))

        # Extraction des valuers de la table sous forme de fonction
        para_functions = {}
        varc_values = table_values.pop(varc_name)

        # Traitement particulier pour les parametres constants type TEMP_DEF_ALPHA
        constant_para_names = ("TEMP_DEF_ALPHA",)
        for cpara_name in constant_para_names:
            if cpara_name in table_values:
                cpara_values = list(set(table_values.pop(cpara_name)))
                if len(cpara_values) != 1:
                    errmsg = "Parameter '{0}' is not constant.".format(cpara_name)
                    UTMESS("F", "SUPERVIS_4", valk=(cmd_name, errmsg))

                para_functions[cpara_name] = cpara_values[0]

        # Tous les autres
        for vpara_name, vpara_values in table_values.items():
            para_f = Function()
            para_f.setParameterName(varc_name)
            para_f.setResultName(vpara_name)
            para_f.setExtrapolation(
                "{0}{1}".format(table_kws["PROL_GAUCHE"][0], table_kws["PROL_DROITE"][0])
            )
            para_f.setInterpolation("{0} {0}".format(table_kws["INTERPOL"]))
            para_f.setValues(varc_values, vpara_values)

            para_functions[vpara_name] = para_f

        # On ajoute le contenu de la table à DEFI_MATERIAU
        for mater in table_kws["COMPOR"]:
            mater_kws = {}
            for key in cmd_kws[mater].keywords.keys():
                if key in para_functions:
                    mater_kws[key] = para_functions[key]

            keywords[mater] = mater_kws

        valk = (
            ", ".join(table_kws["COMPOR"]),
            "{0} ('<{1}>')".format(table_kws["TABLE"].userName, table_kws["TABLE"].getName()),
        )
        UTMESS("I", "MATERIAL1_11", valk=valk)
        keywords.pop("TABLE")
        # check syntax after changing the keywords
        self.check_syntax(keywords)


def adapt_elas(keywords):
    """Compute (E, NU) from (K, MU) or (CELE_P, CELE_S, RHO)

    If (E, NU) already exist, check the consistency.
    Otherwise, assign (E, NU) values.
    """

    def cmp(x, y, epsi=1.0e-3):
        return abs((x - y) / (x + y) * 2.0) < epsi

    elas = keywords.get("ELAS")
    if not elas:
        return

    elas = force_list(elas)[0]
    valK = elas.pop("K", None)
    valMU = elas.pop("MU", None)
    celP = elas.pop("CELE_P", None)
    celS = elas.pop("CELE_S", None)
    if None not in (valK, valMU):
        if valK <= 0.0:
            UTMESS("F", "MATERIAL1_14", valk="K")
        if valMU <= 0.0:
            UTMESS("F", "MATERIAL1_14", valk="MU")
        valE = 9.0 * valK * valMU / (3.0 * valK + valMU)
        valNU = (3.0 * valK - 2.0 * valMU) / (6.0 * valK + 2.0 * valMU)
        UTMESS("I", "MATERIAL1_12", valr=(valE, valNU), valk="(K, MU)")
    elif None not in (celP, celS):
        if celS <= 0.0:
            UTMESS("F", "MATERIAL1_14", valk="CELE_S")
        if celP / celS < 2.0 * sqrt(3.0) / 3.0:
            UTMESS("F", "MATERIAL1_15", valr=(celP / celS, 2.0 * sqrt(3.0) / 3.0))
        valE = (
            elas["RHO"] * celS**2 * (3.0 * celP**2 - 4.0 * celS**2) / (celP**2 - celS**2)
        )
        valNU = (celP**2 - 2.0 * celS**2) / (celP**2 - celS**2) / 2.0
        UTMESS("I", "MATERIAL1_12", valr=(valE, valNU), valk="(CELE_P, CELE_S, RHO)")
    else:
        return

    elas["E"] = valE
    elas["NU"] = valNU


# K = E / (3. * (1 - 2 * NU))
# MU = E / (2. * (1 + NU))
# Vp = sqrt( (K + 4/3 * MU) / rho )
# Vs = sqrt( MU / rho )

DEFI_MATERIAU = MaterialDefinition.run
