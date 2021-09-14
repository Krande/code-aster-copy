# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

# person_in_charge: sebastien.meunier at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

POST_BEREMIN = MACRO(
    nom="POST_BEREMIN",
    op=OPS("code_aster.MacroCommands.post_beremin_ops.post_beremin_ops"),
    sd_prod=table_sdaster,
    docu="UX.YZ.AB",
    reentrant="n",
    fr=tr("Post-traitement de Beremin"),

    RESULTAT=SIMP(statut="o", typ=evol_noli, fr=tr("Résultat mecanique")),
    GROUP_MA=SIMP(statut="o", typ=grma, max="**", fr=tr(
        "Groupe de mailles sur lequel effectuer le post-traitement")),
    DEFORMATION=SIMP(statut="o", typ="TXM", into=("PETIT", "PETIT_REAC",
                                                  "GDEF_LOG"),
                     fr=tr("Type de déformation du résultat mécanique")),
    FILTRE_SIGM=SIMP(statut="o", typ="TXM", into=("SIGM_ELGA", "SIGM_ELMOY"),
                     fr=tr("Option de moyennation des contraintes")),
    COEF_MULT=SIMP(statut="f", typ="R", defaut=1., fr=tr(
        "Coefficient à renseigner selon u4.81.22")),
    UNITE   =SIMP(statut='f',typ=UnitType(), inout='out',),

)
