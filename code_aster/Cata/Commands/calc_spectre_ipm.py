# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

CALC_SPECTRE_IPM = MACRO(
    nom="CALC_SPECTRE_IPM",
    op=OPS("code_aster.MacroCommands.calc_spectre_ipm_ops.calc_spectre_ipm_ops"),
    sd_prod=table_sdaster,
    reentrant="n",
    fr="Calcul de spectre, post-traitement de séisme",
    MAILLAGE=SIMP(statut="f", typ=maillage_sdaster),
    b_maillage=BLOC(
        condition="""exists("MAILLAGE")""",
        EQUIPEMENT=FACT(
            statut="o",
            max="**",
            NOM=SIMP(statut="o", typ="TXM"),
            GROUP_NO=SIMP(statut="o", typ=grno, validators=NoRepeat(), max="**"),
            RAPPORT_MASSE_TOTALE=SIMP(statut="o", typ="R", max=1),
            COEF_MASS_EQUIP=SIMP(statut="o", typ="R", max="**"),
            FREQ_SUPPORT=SIMP(statut="o", typ="R", max=1),
            AMOR_SUPPORT=SIMP(statut="o", typ="R", max=1),
            AMOR_EQUIP=SIMP(statut="o", typ="R", max="**"),
            FREQ_EQUIP=SIMP(statut="o", typ="R", max="**"),
        ),
    ),
    b_no_maillage=BLOC(
        condition="""not exists("MAILLAGE")""",
        EQUIPEMENT=FACT(
            statut="o",
            max="**",
            NOM=SIMP(statut="o", typ="TXM"),
            NOEUD=SIMP(statut="o", typ=no, validators=NoRepeat(), max="**"),
            RAPPORT_MASSE_TOTALE=SIMP(statut="o", typ="R", max=1),
            COEF_MASS_EQUIP=SIMP(statut="o", typ="R", max="**"),
            FREQ_SUPPORT=SIMP(statut="o", typ="R", max=1),
            AMOR_SUPPORT=SIMP(statut="o", typ="R", max=1),
            AMOR_EQUIP=SIMP(statut="o", typ="R", max="**"),
            FREQ_EQUIP=SIMP(statut="o", typ="R", max="**"),
        ),
    ),
    CALCUL=SIMP(statut="o", typ="TXM", into=("ABSOLU", "RELATIF")),
    AMOR_SPEC=SIMP(statut="o", typ="R", max="**"),
    LIST_INST=SIMP(statut="f", typ=listr8_sdaster),
    LIST_FREQ=SIMP(statut="f", typ=listr8_sdaster),
    FREQ=SIMP(statut="f", typ="R", max="**"),
    NORME=SIMP(statut="o", typ="R"),
    b_rela=BLOC(
        condition="""equal_to("CALCUL", 'RELATIF')""",
        RESU=FACT(
            statut="o",
            max=1,
            regles=(UN_PARMI("TABLE", "FONCTION"),),
            TABLE=SIMP(statut="f", typ=table_sdaster),
            FONCTION=SIMP(statut="f", typ=fonction_sdaster),
            ACCE_Z=SIMP(statut="o", typ=fonction_sdaster),
        ),
    ),
    b_abso=BLOC(
        condition="""equal_to("CALCUL", 'ABSOLU')""",
        RESU=FACT(
            statut="o",
            max=1,
            regles=(UN_PARMI("TABLE", "FONCTION"),),
            TABLE=SIMP(statut="f", typ=table_sdaster),
            FONCTION=SIMP(statut="f", typ=fonction_sdaster),
        ),
    ),
    TOLE_INIT=SIMP(statut="f", typ="R", max=1, defaut=1e-3),
    CORR_INIT=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
)
