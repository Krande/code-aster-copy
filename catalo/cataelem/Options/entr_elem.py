# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PCONTPR = InputParameter(
    phys=PHY.SIEF_R,
    container="RESU!SIEF_ELGA!N",
    comment="""  PCONTRR : CONTRAINTES INSTANT ACTUEL """,
)


PVARIPR = InputParameter(
    phys=PHY.VARI_R,
    container="RESU!VARI_ELGA!N",
    comment="""  PVARIPR : VARIABLES INTERNES AUX POINTS DE GAUSS """,
)


PCOMPOR = InputParameter(
    phys=PHY.COMPOR, container="RESU!COMPORTEMENT!N", comment=""" PCOMPOR  :  COMPORTEMENT """
)


PVARCPR = InputParameter(
    phys=PHY.VARI_R,
    container="VOLA!&&CCPARA.VARI_INT_N",
    comment="""  PVARCPR : TEMPERATURES INSTANT ACTUEL """,
)


# ======================================================================
#  Option utilisée dans CALC_GP
# ======================================================================
ENTR_ELEM = Option(
    para_in=(
        SP.PCACOQU,
        PCOMPOR,
        PCONTPR,
        SP.PDEPLR,
        SP.PGEOMER,
        SP.PMATERC,
        PVARCPR,
        SP.PVARCRR,
        PVARIPR,
    ),
    para_out=(SP.PENTRD1,),
    condition=(CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),),
    comment="""  ENTR_ELEM : ENERGIE ELASTIQUE MODIFIEE PAR ELEMENT (TRACTION)""",
)
