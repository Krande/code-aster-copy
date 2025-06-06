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

# person_in_charge: jacques.pellet at edf.fr


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PCOMPOR = InputParameter(
    phys=PHY.COMPOR, comment=""" LE COMPORTEMENT PERMET DE DETERMINER NCMP_DYN """
)


PNBSP_I = InputParameter(
    phys=PHY.NBSP_I, comment=""" LE CHAMP DE NBSP_I PERMET DE DETERMINER NPG_DYN """
)


NSPG_NBVA = Option(
    para_in=(PCOMPOR, PNBSP_I),
    para_out=(SP.PDCEL_I,),
    condition=(
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),
        CondCalcul("+", ((AT.PHENO, "TH"), (AT.BORD, "0"))),
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD_ISO, "OUI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.FSI, "OUI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.ABSO, "OUI"), (AT.FLUIDE, "OUI"))),
    ),
    comment=""" OPTION SERVANT A CALCULER LES 2 NOMBRES :
   NPG_DYN  : NOMBRE DE SOUS-POINTS (POUR LES ELEMENTS DE STRUCTURE EN GENERAL)
   NCMP_DYN : NOMBRE DE COMPOSANTES POUR LA GRANDEUR VARI_R
""",
)
