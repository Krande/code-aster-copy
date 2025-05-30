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


PCONTRR = InputParameter(
    phys=PHY.SIEF_R,
    container="RESU!SIGM_ELGA!N",
    comment="""  PCONTRR : CONTRAINTES PAR ELEMENTS AUX POINTS DE GAUSS """,
)


PCONTEQ = OutputParameter(
    phys=PHY.SIEF_R,
    type="ELGA",
    comment="""  PCONTEQ : CONTRAINTES EQUIVALENTES AUX POINTS DE GAUSS """,
)


SIEQ_ELGA = Option(
    para_in=(PCONTRR,),
    para_out=(PCONTEQ,),
    condition=(
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.EFGE, "OUI"), (AT.SIGM, "NON"))),
    ),
    comment="""  SIEQ_ELGA : CONTRAINTES EQUIVALENTES PAR ELEMENTS AUX POINTS DE GAUSS """,
)
