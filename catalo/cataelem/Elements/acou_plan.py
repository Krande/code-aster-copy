# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "W")
)


EINTENO = LocatedComponents(
    phys=PHY.INTE_R, type="ELNO", components=("INTX_R", "INTY_R", "INTX_I", "INTY_I")
)


DDL_ACOU = LocatedComponents(phys=PHY.PRES_C, type="ELNO", components=("PRES",))

MVECTTC = ArrayOfComponents(phys=PHY.VPRE_C, locatedComponents=DDL_ACOU)

MMATTTC = ArrayOfComponents(phys=PHY.MPRE_C, locatedComponents=DDL_ACOU)


# ------------------------------------------------------------
class ACPLQU4(Element):
    """Please document this element"""

    meshType = MT.QUAD4
    elrefe = (
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "FPG1=FPG1", "NOEU=NOEU"), mater=("FPG1",)),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )
    calculs = (
        OP.COOR_ELGA(
            te=479, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.INTE_ELNO(
            te=175,
            para_in=(
                (SP.PFREQR, LC.CFREQR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPRESSC, DDL_ACOU),
            ),
            para_out=((SP.PINTER, EINTENO),),
        ),
        OP.MASS_ACOU(
            te=177,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PMATTTC, MMATTTC),),
        ),
        OP.PRAC_ELNO(
            te=189, para_in=((SP.PPRESSC, DDL_ACOU),), para_out=((SP.PPRAC_R, LC.EPRACNO),)
        ),
        OP.RIGI_ACOU(te=176, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PMATTTC, MMATTTC),)),
        OP.TOU_INI_ELGA(te=99, para_out=((OP.TOU_INI_ELGA.PGEOM_R, EGGEOP_R),)),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),)),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# ------------------------------------------------------------
class ACPLQU8(ACPLQU4):
    """Please document this element"""

    meshType = MT.QUAD8
    elrefe = (
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )


# ------------------------------------------------------------
class ACPLQU9(ACPLQU4):
    """Please document this element"""

    meshType = MT.QUAD9
    elrefe = (
        ElrefeLoc(MT.QU9, gauss=("RIGI=FPG9", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )


# ------------------------------------------------------------
class ACPLTR3(ACPLQU4):
    """Please document this element"""

    meshType = MT.TRIA3
    elrefe = (
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG1", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )


# ------------------------------------------------------------
class ACPLTR6(ACPLQU4):
    """Please document this element"""

    meshType = MT.TRIA6
    elrefe = (
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG3", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )


# ------------------------------------------------------------
class ACPLNASE2(Element):
    """Please document this element"""

    meshType = MT.SEG2
    elrefe = (ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("FPG1",)),)
    calculs = (
        OP.CHAR_ACOU_VFAC_C(
            te=179,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC), (SP.PVITEFC, LC.EVITEFC)),
            para_out=((SP.PVECTTC, MVECTTC),),
        ),
        OP.COOR_ELGA(
            te=478, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.TOU_INI_ELGA(te=99, para_out=((OP.TOU_INI_ELGA.PGEOM_R, EGGEOP_R),)),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),)),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
    )


# ------------------------------------------------------------
class ACPLNASE3(ACPLNASE2):
    """Please document this element"""

    meshType = MT.SEG3
    elrefe = (ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("FPG1",)),)
