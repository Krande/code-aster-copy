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


PCOMPOR = InputParameter(phys=PHY.COMPOR, container="RESU!COMPORTEMENT!N")


PCAORIE = InputParameter(
    phys=PHY.CAORIE_R,
    container="CARA!.CARORIEN",
    comment="""  CAORIE : ORIENTATION LOCALE D'UN ELEMENT DE POUTRE OU DE TUYAU,
           ISSUE DE AFFE_CARA_ELEM MOT CLE ORIENTATION """,
)


PVARCPR = InputParameter(
    phys=PHY.VARI_R,
    container="VOLA!&&CCPARA.VARI_INT_N",
    comment="""  PVARCPR : TEMPERATURES DE TYPE REEL INSTANT ACTUEL """,
)


PNBSP_I = InputParameter(
    phys=PHY.NBSP_I, container="CARA!.CANBSP", comment="""  PNBSP_I :  NOMBRE DE SOUS_POINTS """
)


PSTRXRR = InputParameter(
    phys=PHY.STRX_R,
    container="RESU!STRX_ELGA!N",
    comment="""  PSTRXRR : CHAMPS ELEMENTS DE STRUCTURE INSTANT ACTUEL """,
)


DEGE_ELNO = Option(
    para_in=(
        SP.PCACOQU,
        SP.PCAGEPO,
        SP.PCAGNPO,
        PCAORIE,
        PCOMPOR,
        SP.PDEPLAR,
        SP.PGEOMER,
        SP.PMATERC,
        PNBSP_I,
        PSTRXRR,
        PVARCPR,
        SP.PVARCRR,
    ),
    para_out=(SP.PDEFOGR,),
    condition=(CondCalcul("+", ((AT.PHENO, "ME"), (AT.EFGE, "OUI"), (AT.BORD, "0"))),),
    comment=""" DEFORMATIONS GENERALISEES DANS LES ELEMENTS DE STRUCTURE : COQUES OU POUTRES
   CALCULEES AUX NOEUDS DE CHAQUE ELEMENT A PARTIR DES DEPLACEMENTS.""",
)
