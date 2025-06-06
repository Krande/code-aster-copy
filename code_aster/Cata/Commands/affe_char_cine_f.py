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

# person_in_charge: jacques.pellet at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


def affe_char_cine_f_prod(MECA_IMPO, THER_IMPO, **args):
    if args.get("__all__"):
        return (char_cine_ther, char_cine_meca)

    if MECA_IMPO is not None:
        return char_cine_meca
    if THER_IMPO is not None:
        return char_cine_ther
    raise CataError("type de concept resultat non prevu")


AFFE_CHAR_CINE_F = OPER(
    nom="AFFE_CHAR_CINE_F",
    op=101,
    sd_prod=affe_char_cine_f_prod,
    fr=tr(
        "Affectation de conditions aux limites cinématiques fonction d'un (ou plusieurs) paramètres"
        " pour un traitement sans dualisation"
    ),
    reentrant="n",
    regles=(UN_PARMI("MECA_IMPO", "THER_IMPO")),
    MODELE=SIMP(statut="o", typ=modele_sdaster),
    MECA_IMPO=FACT(
        statut="f",
        max="**",
        regles=(
            UN_PARMI("TOUT", "GROUP_MA", "GROUP_NO"),
            AU_MOINS_UN(
                "DX",
                "DY",
                "DZ",
                "DRX",
                "DRY",
                "DRZ",
                "GRX",
                "PRES",
                "PHI",
                "TEMP",
                "PRE1",
                "PRE2",
                "UI2",
                "UI3",
                "VI2",
                "VI3",
                "WI2",
                "WI3",
                "UO2",
                "UO3",
                "VO2",
                "VO3",
                "WO2",
                "WO3",
                "UI4",
                "UI5",
                "VI4",
                "VI5",
                "WI4",
                "WI5",
                "UO4",
                "UO5",
                "VO4",
                "VO5",
                "WO4",
                "WO5",
                "UI6",
                "UO6",
                "VI6",
                "VO6",
                "WI6",
                "WO6",
                "WO",
                "WI1",
                "WO1",
                "GONF",
                "H1X",
                "H1Y",
                "H1Z",
                "K1",
                "K2",
                "K3",
                "V11",
                "V12",
                "V13",
                "V21",
                "V22",
                "V23",
                "V31",
                "V32",
                "V33",
                "PRES11",
                "PRES12",
                "PRES13",
                "PRES21",
                "PRES22",
                "PRES23",
                "PRES31",
                "PRES32",
                "PRES33",
                "LH1",
                "GLIS",
                "PSI",
            ),
        ),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        MAILLE=SIMP(statut="c", typ=ma, validators=NoRepeat(), max="**"),
        GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
        NOEUD=SIMP(statut="c", typ=no, validators=NoRepeat(), max="**"),
        DX=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        DY=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        DZ=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        DRX=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        DRY=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        DRZ=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        GRX=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRES=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PHI=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        TEMP=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRE1=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRE2=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        UI2=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        UI3=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        UI4=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        UI5=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        UI6=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        UO2=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        UO3=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        UO4=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        UO5=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        UO6=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VI2=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VI3=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VI4=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VI5=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VI6=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VO2=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VO3=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VO4=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VO5=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VO6=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WI2=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WI3=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WI4=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WI5=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WI6=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WO2=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WO3=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WO4=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WO5=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WO6=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WO=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WI1=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        WO1=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        GONF=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        H1X=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        H1Y=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        H1Z=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        K1=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        K2=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        K3=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        V11=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        V12=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        V13=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        V21=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        V22=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        V23=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        V31=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        V32=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        V33=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRES11=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRES12=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRES13=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRES21=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRES22=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRES23=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRES31=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRES32=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PRES33=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        LH1=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        GLIS=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        PSI=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    ),
    THER_IMPO=FACT(
        statut="f",
        max="**",
        regles=(
            UN_PARMI("TOUT", "GROUP_MA", "GROUP_NO"),
            AU_MOINS_UN("TEMP", "TEMP_MIL", "TEMP_INF", "TEMP_SUP"),
        ),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        MAILLE=SIMP(statut="c", typ=ma, validators=NoRepeat(), max="**"),
        GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
        NOEUD=SIMP(statut="c", typ=no, validators=NoRepeat(), max="**"),
        TEMP=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        TEMP_MIL=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        TEMP_SUP=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
        TEMP_INF=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    ),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
