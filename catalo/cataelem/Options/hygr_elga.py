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


from cataelem.Tools.base_objects import InputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT

PVARCPR = InputParameter(
    phys=PHY.VARI_R,
    container="VOLA!&&CCPARA.VARI_INT_N",
    comment="""  PVARCPR : VARIABLES DE COMMANDE  """,
)


HYGR_ELGA = Option(
    para_in=(SP.PGEOMER, SP.PMATERC, SP.PSECHRR, SP.PINSTR, PVARCPR),
    para_out=(SP.PHYGRPG,),
    condition=(CondCalcul("+", ((AT.PHENO, "TH"), (AT.SECH, "OUI"), (AT.BORD, "0"))),),
    comment="""  HYGR_ELGA : CALCUL DE GRANDEURS LIEES A L'HYGROMETRIE AUX POINTS DE GAUSS """,
)
