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


def depl_interne_prod(DEPL_GLOBAL, **args):
    if args.get("__all__"):
        return (cham_no_sdaster, evol_elas, dyna_trans, dyna_harmo, mode_meca, mode_meca_c)

    if AsType(DEPL_GLOBAL) == cham_no_sdaster:
        return cham_no_sdaster
    if AsType(DEPL_GLOBAL) == evol_elas:
        return evol_elas
    if AsType(DEPL_GLOBAL) == dyna_trans:
        return dyna_trans
    if AsType(DEPL_GLOBAL) == dyna_harmo:
        return dyna_harmo
    if AsType(DEPL_GLOBAL) == mode_meca:
        return mode_meca
    if AsType(DEPL_GLOBAL) == mode_meca_c:
        return mode_meca_c
    raise CataError("type de concept resultat non prevu")


DEPL_INTERNE = OPER(
    nom="DEPL_INTERNE",
    op=89,
    sd_prod=depl_interne_prod,
    reentrant="n",
    fr=tr("Calculer le champ de déplacement à l'intérieur d'une sous-structure statique"),
    DEPL_GLOBAL=SIMP(
        statut="o", typ=(cham_no_sdaster, mode_meca, mode_meca_c, evol_elas, dyna_trans, dyna_harmo)
    ),
    SUPER_MAILLE=SIMP(statut="o", typ=ma),
    NOM_CAS=SIMP(statut="f", typ="TXM", defaut=" "),
)
