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

# person_in_charge: jacques.pellet at edf.fr


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PCONTRR = InputParameter(
    phys=PHY.SIEF_R,
    container="RESU!SIEF_ELGA!N",
    comment="""  PCONTRR : CONTRAINTES INSTANT ACTUEL """,
)


PVARCPR = InputParameter(
    phys=PHY.VARI_R,
    container="VOLA!&&CCPARA.VARI_INT_N",
    comment="""  PVARCPR : TEMPERATURES INSTANT ACTUEL """,
)


PCOMPOR = InputParameter(
    phys=PHY.COMPOR, container="RESU!COMPORTEMENT!N", comment=""" PCOMPOR  :  COMPORTEMENT """
)


PNBSP_I = InputParameter(
    phys=PHY.NBSP_I, container="CARA!.CANBSP", comment=""" PNBSP_I  : NOMBRE DE SOUS_POINTS """
)


PENERDR = OutputParameter(
    phys=PHY.ENER_R,
    type="ELGA",
    comment="""  PENERDR : DENSITE D'ENERGIE ELASTIQUE AUX POINTS DE GAUSS """,
)


ENEL_ELGA = Option(
    para_in=(
        SP.PCACOQU,
        SP.PCAMASS,
        PCOMPOR,
        PCONTRR,
        SP.PDEPLAR,
        SP.PGEOMER,
        SP.PMATERC,
        PNBSP_I,
        SP.PINSTR,
        PVARCPR,
        SP.PVARCRR,
        SP.PVARIGR,
    ),
    para_out=(PENERDR,),
    condition=(CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),),
    comment="""  ENEL_ELGA : DENSITE D ENERGIE ELASTIQUE PAR ELEMENT AUX POINTS DE GAUSS """,
)
