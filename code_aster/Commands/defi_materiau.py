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


from ..Messages import UTMESS
from ..Objects import Function
from ..Supervis import ExecuteMacro


class MaterialDefinition(ExecuteMacro):
    """Command that defines a Material."""

    command_name = "DEFI_MATERIAU"

    def adapt_syntax(self, keywords):
        """Adapt syntax to replace TABLE content in keywords.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
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


DEFI_MATERIAU = MaterialDefinition.run
