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

# person_in_charge: francesco.bettonte at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

REST_HOMO_LOCA = MACRO(
    nom="REST_HOMO_LOCA",
    op=OPS("code_aster.MacroCommands.MateHomo.rest_homo_ops.rest_homo_ops"),
    sd_prod=evol_elas,
    docu="UX.YZ.AB",
    reentrant="n",
    fr=tr("Calcul des champs locaux à partir de la solution homogénéisée."),
    TYPE_HOMO=SIMP(statut="o", typ="TXM", into=("MASSIF",)),
    POSITION=SIMP(statut="o", typ="R", min=3, max=3),
    CORR_MECA=SIMP(statut="o", typ=evol_elas_dict),
    CORR_THER=SIMP(statut="f", typ=evol_ther_dict),
    TEMP_REF=SIMP(statut="f", typ="R", defaut=20.0),
    PRESSION=SIMP(statut="f", typ="R", defaut=0.0),
    EVOL_ELAS=SIMP(statut="o", typ=evol_elas),
    EVOL_THER=SIMP(statut="f", typ=evol_ther),
    RESU_THER=SIMP(statut="f", typ=CO),
    TOUT_ORDRE=SIMP(statut="f", typ="TXM", into=("OUI",)),
    INST=SIMP(statut="f", typ="R", validators=NoRepeat(), max="**"),
    CRITERE=SIMP(statut="f", typ="TXM", defaut="RELATIF", into=("RELATIF", "ABSOLU")),
    PRECISION=SIMP(statut="f", typ="R", defaut=1.0e-6, min=0.0),
    regles=(
        ENSEMBLE("EVOL_THER", "CORR_THER"),
        PRESENT_PRESENT("RESU_THER", "CORR_THER", "EVOL_THER"),
        UN_PARMI("TOUT_ORDRE", "INST"),
    ),
    AFFE=FACT(
        statut="o",
        max="**",
        regles=UN_PARMI("TOUT", "GROUP_MA"),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        MATER=SIMP(statut="o", typ=mater_sdaster, max=1),
    ),
    OPTION=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON", max=1),
    SYME=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON", max=1),
    TRAN=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON", max=1),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
