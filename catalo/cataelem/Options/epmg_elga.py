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


PVARCPR = InputParameter(
    phys=PHY.VARI_R,
    container="VOLA!&&CCPARA.VARI_INT_N",
    comment="""  PVARCPR : VARIABLES DE COMMANDE  """,
)


PCOMPOR = InputParameter(
    phys=PHY.COMPOR, container="RESU!COMPORTEMENT!N", comment=""" COMPORTEMENT """
)


PDEFOPG = OutputParameter(
    phys=PHY.EPSI_R,
    type="ELGA",
    comment="""  PDEFOPG : DEFORMATIONS MECANIQUES A PARTIR DES DEPLACEMENTS
                     AUX POINTS DE GAUSS """,
)


EPMG_ELGA = Option(
    para_in=(
        SP.PCAMASS,
        PCOMPOR,
        SP.PDEPLAR,
        SP.PGEOMER,
        SP.PMATERC,
        SP.PINSTR,
        PVARCPR,
        SP.PVARCRR,
    ),
    para_out=(PDEFOPG,),
    condition=(CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),),
    comment="""  EPMG_ELGA : DEFORMATIONS MECANIQUES A PARTIR DES DEPLACEMENTS,
           THEORIE DES GRANDS DEPLACEMENTS PAR ELEMENTS AUX POINTS DE GAUSS """,
)
