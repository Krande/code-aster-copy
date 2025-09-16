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


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT

PCHHOBS = InputParameter(phys=PHY.N3600R, comment=""" HHO - coefficient base locale""")

PTEMP_R = OutputParameter(
    phys=PHY.TEMP_R, type="ELNO", comment=""" HHO - degres de liberte de la cellule"""
)

PFUNC_R = InputParameter(phys=PHY.NEUT_K8, comment=""" Value to project""")


HHO_PROJ_THER = Option(
    para_in=(SP.PGEOMER, PFUNC_R, SP.PINSTPR, PCHHOBS),
    para_out=(PTEMP_R,),
    condition=(CondCalcul("+", ((AT.PHENO, "TH"), (AT.BORD, "0"), (AT.HHO, "OUI"))),),
)
