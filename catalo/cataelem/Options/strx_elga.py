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
    comment="""  PVARCPR : VARIABLES DE COMMANDE  """,
)


PCAORIE = InputParameter(
    phys=PHY.CAORIE_R,
    container="CARA!.CARORIEN",
    comment="""  PCAORIE : ORIENTATION LOCALE D'UN ELEMENT DE POUTRE OU DE TUYAU,
           ISSUE DE AFFE_CARA_ELEM MOT CLE ORIENTATION """,
)


PNBSP_I = InputParameter(
    phys=PHY.NBSP_I, container="CARA!.CANBSP", comment="""  PNBSP_I :  NOMBRE DE SOUS_POINTS """
)


PCOMPOR = InputParameter(
    phys=PHY.COMPOR,
    container="CHMA!.COMPOR",
    comment="""  PCOMPOR :  DESCRIPTION DU COMPORTEMENT DE CHAQUE GROUPE DE FIBRES
           NECESSITE DE FOURNIR LE CONCEPT PRODUIT PAR AFFE_MATERIAU  """,
)


PSTRXRR = OutputParameter(phys=PHY.STRX_R, type="ELGA")


STRX_ELGA = Option(
    para_in=(
        SP.PCAGNPO,
        PCAORIE,
        PCOMPOR,
        SP.PDEPLAR,
        SP.PFIBRES,
        SP.PGEOMER,
        SP.PMATERC,
        PNBSP_I,
        PVARCPR,
        SP.PVARCRR,
    ),
    para_out=(PSTRXRR,),
    condition=(CondCalcul("+", ((AT.PHENO, "ME"), (AT.MODELI, "PFM"))),),
    comment=""" CALCUL EFFORTS GENERALISES AUX POINTS DE GAUSS
   A PARTIR DES DEPLACEMENTS. LICITE EN LINEAIRE SEULEMENT. """,
)
