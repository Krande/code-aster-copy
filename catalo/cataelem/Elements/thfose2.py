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


CCOEFHF = LocatedComponents(phys=PHY.COEH_F, type="ELEM", components=("H",))


CCOEFHR = LocatedComponents(phys=PHY.COEH_R, type="ELEM", components=("H",))


CFLUXNF = LocatedComponents(phys=PHY.FLUN_F, type="ELEM", components=("FLUN",))


CFLUXNR = LocatedComponents(phys=PHY.FLUN_R, type="ELEM", components=("FLUN",))


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "W")
)


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST", "DELTAT", "THETA"))


CT_EXTF = LocatedComponents(phys=PHY.TEMP_F, type="ELEM", components=("TEMP",))


DDL_THER = LocatedComponents(phys=PHY.TEMP_R, type="ELNO", components=("TEMP",))


MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)


# ------------------------------------------------------------
class THFOSE2(Element):
    """Please document this element"""

    meshType = MT.SEG2
    elrefe = (ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4",)),)
    calculs = (
        OP.CHAR_THER_FLUN_F(
            te=272,
            para_in=((SP.PFLUXNF, CFLUXNF), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_R(
            te=271,
            para_in=((SP.PFLUXNR, CFLUXNR), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ECHA_F(
            te=270,
            para_in=(
                (SP.PCOEFHF, CCOEFHF),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, CTEMPSR),
                (SP.PT_EXTF, CT_EXTF),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ECHA_R(
            te=269,
            para_in=(
                (SP.PCOEFHR, CCOEFHR),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, CTEMPSR),
                (SP.PT_EXTR, LC.ET_EXTR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.COOR_ELGA(
            te=478, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.RIGI_THER_ECHA_F(
            te=268,
            para_in=((SP.PCOEFHF, CCOEFHF), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((OP.RIGI_THER_ECHA_F.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_ECHA_R(
            te=267,
            para_in=((SP.PCOEFHR, CCOEFHR), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((OP.RIGI_THER_ECHA_R.PMATTTR, MMATTTR),),
        ),
        OP.TOU_INI_ELGA(te=99, para_out=((OP.TOU_INI_ELGA.PGEOM_R, EGGEOP_R),)),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),)),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
    )


# ------------------------------------------------------------
class THFOSE3(THFOSE2):
    """Please document this element"""

    meshType = MT.SEG3
    elrefe = (ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),)
