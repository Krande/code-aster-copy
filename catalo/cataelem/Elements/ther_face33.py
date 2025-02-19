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

from cataelem.Tools.base_objects import LocatedComponents, ArrayOfComponents, SetOfNodes, ElrefeLoc
from cataelem.Tools.base_objects import Calcul, Element
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.located_components as LC
import cataelem.Commons.parameters as SP
import cataelem.Commons.mesh_types as MT
from cataelem.Options.options import OP

# ----------------
# Modes locaux :
# ----------------


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST", "DELTAT", "THETA"))


DDL_THER = LocatedComponents(phys=PHY.TEMP_R, type="ELNO", components=("TEMP",))


MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)


# ------------------------------------------------------------
class THER_FACE33(Element):
    """Please document this element"""

    meshType = MT.TRIA33
    elrefe = (ElrefeLoc(MT.TR3, gauss=("RIGI=COT3",)),)
    calculs = (
        OP.CHAR_THER_PARO_F(
            te=304,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PHECHPF, LC.CHECHPF),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_PARO_R(
            te=303,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PHECHPR, LC.CHECHPR),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.MTAN_THER_PARO_F(
            te=389,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PHECHPF, LC.CHECHPF), (SP.PINSTR, CTEMPSR)),
            para_out=((OP.MTAN_THER_PARO_F.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_PARO_R(
            te=388,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PHECHPR, LC.CHECHPR), (SP.PINSTR, CTEMPSR)),
            para_out=((OP.MTAN_THER_PARO_R.PMATTTR, MMATTTR),),
        ),
        OP.RESI_THER_PARO_F(
            te=278,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PHECHPF, LC.CHECHPF),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_PARO_R(
            te=276,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PHECHPR, LC.CHECHPR),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RIGI_THER_PARO_F(
            te=302,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PHECHPF, LC.CHECHPF), (SP.PINSTR, CTEMPSR)),
            para_out=((OP.RIGI_THER_PARO_F.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_PARO_R(
            te=301,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PHECHPR, LC.CHECHPR), (SP.PINSTR, CTEMPSR)),
            para_out=((OP.RIGI_THER_PARO_R.PMATTTR, MMATTTR),),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
    )


# ------------------------------------------------------------
class THER_FACE44(THER_FACE33):
    """Please document this element"""

    meshType = MT.QUAD44
    elrefe = (ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4",)),)


# ------------------------------------------------------------
class THER_FACE66(THER_FACE33):
    """Please document this element"""

    meshType = MT.TRIA66
    elrefe = (ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6",)),)


# ------------------------------------------------------------
class THER_FACE88(THER_FACE33):
    """Please document this element"""

    meshType = MT.QUAD88
    elrefe = (ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9",)),)


# ------------------------------------------------------------
class THER_FACE99(THER_FACE33):
    """Please document this element"""

    meshType = MT.QUAD99
    elrefe = (ElrefeLoc(MT.QU9, gauss=("RIGI=FPG9",)),)
