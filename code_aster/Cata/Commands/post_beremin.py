# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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


from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


def post_beremin_prod(self, SIGM_MAXI, **args):
    if args.get("__all__"):
        return ([table_sdaster], [table_sdaster, evol_noli])

    if SIGM_MAXI is not None:
        self.type_sdprod(SIGM_MAXI, evol_noli)

    return table_sdaster


POST_BEREMIN = MACRO(
    nom="POST_BEREMIN",
    op=OPS("code_aster.MacroCommands.post_beremin_ops.post_beremin_ops"),
    sd_prod=post_beremin_prod,
    docu="U4.81.08",
    reentrant="n",
    fr=tr("Post-traitement de Beremin"),
    RESULTAT=SIMP(statut="o", typ=evol_noli, fr=tr("Résultat mecanique")),
    GROUP_MA=SIMP(
        statut="o",
        typ=grma,
        max=1,
        fr=tr("Groupe de mailles sur lequel effectuer le post-traitement"),
    ),
    DEFORMATION=SIMP(
        statut="o",
        typ="TXM",
        into=("PETIT", "PETIT_REAC", "GDEF_LOG"),
        fr=tr("Type de déformation du résultat mécanique"),
    ),
    FILTRE_SIGM=SIMP(
        statut="o",
        typ="TXM",
        into=("SIGM_ELGA", "SIGM_ELMOY"),
        fr=tr("Option de moyennation des contraintes"),
    ),
    NUME_VARI=SIMP(
        statut="o",
        typ="I",
        fr=tr("Numéro de la variable interne INDIPLAS"),
        val_min=1,
        min=1,
        max=1,
    ),
    b_gdeflog=BLOC(
        condition="""equal_to("DEFORMATION", 'GDEF_LOG')""",
        LIST_NUME_SIEF=SIMP(
            statut="o",
            typ="I",
            fr=tr("Numéros des contraintes de Kirchhoff "),
            val_min=1,
            min=4,
            max=6,
        ),
    ),
    HIST_MAXI=SIMP(
        statut="f", typ="TXM", into=("OUI", "NON"), defaut="OUI", fr=tr("Maximum temporel ou pas")
    ),
    SIGM_MAXI=SIMP(
        statut="f", typ=CO, fr=tr("Maximum temporel de la contrainte principale majeure")
    ),
    COEF_MULT=SIMP(
        statut="f", typ="R", defaut=1.0, fr=tr("Coefficient à renseigner selon u4.81.22")
    ),
    WEIBULL=FACT(
        statut="f",
        M=SIMP(statut="o", typ="R"),
        VOLU_REFE=SIMP(statut="o", typ="R"),
        SIGM_REFE=SIMP(statut="o", typ="R"),
        SIGM_SEUIL=SIMP(statut="f", typ="R", defaut=0.0, val_min=0.0),
    ),
    WEIBULL_FO=FACT(
        statut="f",
        M=SIMP(statut="o", typ="R"),
        VOLU_REFE=SIMP(statut="o", typ="R"),
        SIGM_CNV=SIMP(statut="o", typ="R"),
        SIGM_REFE=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule)),
        SIGM_SEUIL=SIMP(statut="f", typ="R", defaut=0.0, val_min=0.0),
    ),
)
