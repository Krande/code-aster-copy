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

# person_in_charge: jean-luc.flejou at edf.fr


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PVARCPR = InputParameter(
    phys=PHY.VARI_R,
    container="VOLA!&&CCPARA.VARI_INT_N",
    comment="""  PVARCPR : VARIABLES DE COMMANDE """,
)


PCAORIE = InputParameter(
    phys=PHY.CAORIE_R,
    container="CARA!.CARORIEN",
    comment="""  PCAORIE : ORIENTATION LOCALE D'UN ELEMENT DE POUTRE OU DE TUYAU,
           ISSUE DE AFFE_CARA_ELEM MOT CLE ORIENTATION """,
)


SIPO_ELNO = Option(
    para_in=(
        SP.PABSCUR,
        SP.PCAARPO,
        SP.PCAGEPO,
        SP.PCAGNPO,
        PCAORIE,
        SP.PCHDYNR,
        SP.PCOEFFC,
        SP.PCOEFFR,
        SP.PDEPLAR,
        SP.PFF1D1D,
        SP.PFR1D1D,
        SP.PGEOMER,
        SP.PMATERC,
        SP.PPESANR,
        SP.PSUROPT,
        SP.PINSTR,
        PVARCPR,
        SP.PVARCRR,
    ),
    para_out=(SP.PCONTPC, SP.PCONTPO),
    condition=(CondCalcul("+", ((AT.PHENO, "ME"), (AT.DIM_TOPO_MODELI, "1"))),),
    comment="""  SIPO_ELNO : CALCUL DES CONTRAINTES AUX NOEUDS DANS LA SECTION
           DE POUTRE DECOMPOSEE EN CONTRIBUTIONS DE CHAQUE EFFORT GENERALISE
           LICITE EN LINEAIRE SEULEMENT. """,
)
