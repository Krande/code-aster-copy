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

# person_in_charge: matthieu.le-cren at edf.fr


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PVARCPR = InputParameter(phys=PHY.VARI_R)


PCOMPOR = InputParameter(phys=PHY.COMPOR)


PCONTRR = InputParameter(phys=PHY.SIEF_R)


PVARIPR = InputParameter(phys=PHY.VARI_R)


PTHETAR = InputParameter(phys=PHY.THET_R)


PDEG = InputParameter(phys=PHY.NEUT_I)


PLAG = InputParameter(phys=PHY.NEUT_R)


PCER = InputParameter(phys=PHY.NEUT_R)


PELI = InputParameter(phys=PHY.NEUT_R)

CALC_G = Option(
    para_in=(
        SP.PACCELE,
        PCOMPOR,
        SP.PCONTGR,
        PCONTRR,
        SP.PDEPLAR,
        SP.PDEPLA,
        SP.PEPSINR,
        SP.PFR1D2D,
        SP.PFR2D3D,
        SP.PFRVOLU,
        SP.PGEOMER,
        SP.PMATERC,
        SP.PPESANR,
        SP.PPRESSR,
        SP.PROTATR,
        SP.PSIGINR,
        PTHETAR,
        PVARCPR,
        SP.PVARCRR,
        SP.PVITESS,
        PDEG,
        PLAG,
        PCER,
        PELI,
    ),
    para_out=(SP.PGTHETA,),
    condition=(
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "-1"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.ABSO, "OUI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.DISCRET, "OUI"))),
    ),
)
