# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

from ...Cata.Commons import *
from ...Cata.Commons.c_comportement import compat_syntax
from ...Cata.DataStructure import *
from ...Cata.Syntax import *

CALC_POINT_MAT_CATA = OPER(
    nom="CALC_POINT_MAT",
    op=33,
    compat_syntax=compat_syntax,
    sd_prod=table_sdaster,
    reentrant="f",
    fr=tr("Intégrer une loi de comportement"),
    MATER=SIMP(statut="o", typ=mater_sdaster, max=30),
    COMPORTEMENT=C_COMPORTEMENT("SIMU_POINT_MAT"),
    INCREMENT=C_INCREMENT(),
    NEWTON=C_NEWTON(),
    CONVERGENCE=C_CONVERGENCE("SIMU_POINT_MAT"),
    # --MASSIF : orientation du materiau (monocristal, orthotropie)
    MASSIF=FACT(
        statut="f",
        max="**",
        regles=(UN_PARMI("ANGL_REP", "ANGL_EULER"),),
        ANGL_REP=SIMP(statut="f", typ="R", min=1, max=3),
        ANGL_EULER=SIMP(statut="f", typ="R", min=1, max=3),
    ),
    ## ANGLE : rotation de ANGLE autour de Z uniquement, et seulement pour les déformations imposées.
    ANGLE=SIMP(statut="f", typ="R", max=1, defaut=0.0),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
    regles=(
        EXCLUS("SIXX", "EPXX"),
        EXCLUS("SIYY", "EPYY"),
        EXCLUS("SIZZ", "EPZZ"),
        EXCLUS("SIXY", "EPXY"),
        EXCLUS("SIXZ", "EPXZ"),
        EXCLUS("SIYZ", "EPYZ"),
        ENSEMBLE("F11", "F12", "F13", "F21", "F22", "F23", "F31", "F32", "F33"),
    ),
    SIXX=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    SIYY=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    SIZZ=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    SIXY=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    SIXZ=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    SIYZ=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    EPXX=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    EPYY=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    EPZZ=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    EPXY=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    EPXZ=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    EPYZ=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    F11=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    F12=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    F13=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    F21=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    F22=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    F23=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    F31=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    F32=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    F33=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    MATR_C1=FACT(
        statut="f",
        max="**",
        VALE=SIMP(statut="o", typ="R", max=1),
        NUME_LIGNE=SIMP(statut="o", typ="I", max=1, val_min=1, val_max=6),
        NUME_COLONNE=SIMP(statut="o", typ="I", max=1, val_min=1, val_max=12),
    ),
    MATR_C2=FACT(
        statut="f",
        max="**",
        VALE=SIMP(statut="o", typ="R", max=1),
        NUME_LIGNE=SIMP(statut="o", typ="I", max=1, val_min=1, val_max=6),
        NUME_COLONNE=SIMP(statut="o", typ="I", max=1, val_min=1, val_max=12),
    ),
    VECT_IMPO=FACT(
        statut="f",
        max=6,
        VALE=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule), max=1),
        NUME_LIGNE=SIMP(statut="o", typ="I", max=1, val_min=1, val_max=6),
    ),
    SIGM_INIT=FACT(
        statut="f",
        SIXX=SIMP(statut="f", typ="R", max=1, defaut=0.0e0),
        SIYY=SIMP(statut="f", typ="R", max=1, defaut=0.0e0),
        SIZZ=SIMP(statut="f", typ="R", max=1, defaut=0.0e0),
        SIXY=SIMP(statut="f", typ="R", max=1, defaut=0.0e0),
        SIXZ=SIMP(statut="f", typ="R", max=1, defaut=0.0e0),
        SIYZ=SIMP(statut="f", typ="R", max=1, defaut=0.0e0),
    ),
    EPSI_INIT=FACT(
        statut="f",
        EPXX=SIMP(statut="o", typ="R", max=1),
        EPYY=SIMP(statut="o", typ="R", max=1),
        EPZZ=SIMP(statut="o", typ="R", max=1),
        EPXY=SIMP(statut="o", typ="R", max=1),
        EPXZ=SIMP(statut="o", typ="R", max=1),
        EPYZ=SIMP(statut="o", typ="R", max=1),
    ),
    VARI_INIT=FACT(statut="f", VALE=SIMP(statut="o", typ="R", max="**")),
    FORMAT_TABLE=SIMP(
        statut="f", typ="TXM", max=1, into=("CMP_COLONNE", "CMP_LIGNE"), defaut=("CMP_COLONNE")
    ),
    NB_VARI_TABLE=SIMP(statut="f", typ="I", max=1),
    OPER_TANGENT=SIMP(statut="f", typ="TXM", max=1, into=("OUI", "NON"), defaut="NON"),
    ARCHIVAGE=FACT(
        statut="f",
        LIST_INST=SIMP(statut="f", typ=(listr8_sdaster)),
        INST=SIMP(statut="f", typ="R", validators=NoRepeat(), max="**"),
        PAS_ARCH=SIMP(statut="f", typ="I"),
        PRECISION=SIMP(statut="f", typ="R", defaut=1.0e-6),
    ),
    # variables de commandes scalaires, définies par une fonction du temps
    AFFE_VARC=FACT(
        statut="f",
        max="**",
        NOM_VARC=SIMP(
            statut="o",
            typ="TXM",
            into=(
                "TEMP",
                "CORR",
                "IRRA",
                "HYDR",
                "SECH",
                "NEUT1",
                "NEUT2",
                "PFERRITE",
                "PPERLITE",
                "PBAINITE",
                "PMARTENS",
                "PAUSTENI",
                "PCOLDSUM",
                "ALPHPUR",
                "ALPHBETA",
                "BETA",
                "EPSA",
            ),
        ),
        VALE_FONC=SIMP(statut="f", typ=(fonction_sdaster, formule)),
        B_VALE_REF=BLOC(
            condition="NOM_VARC in ('TEMP', 'SECH')", VALE_REF=SIMP(statut="o", typ="R")
        ),
    ),
)
