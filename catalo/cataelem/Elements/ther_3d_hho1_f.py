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

from cataelem.Tools.base_objects import LocatedComponents, ArrayOfComponents, SetOfNodes, ElrefeLoc
from cataelem.Tools.base_objects import Calcul, Element
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.located_components as LC
import cataelem.Commons.parameters as SP
import cataelem.Commons.mesh_types as MT
from cataelem.Options.options import OP
import cataelem.Commons.attributes as AT

# ----------------
# Modes locaux :
# ----------------


CCOEFHF = LocatedComponents(phys=PHY.COEH_F, type="ELEM", components=("H",))


CCOEFHR = LocatedComponents(phys=PHY.COEH_R, type="ELEM", components=("H",))


NACCELR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))


CFLUXNF = LocatedComponents(phys=PHY.FLUN_F, type="ELEM", components=("FLUN",))

CFLUXVF = LocatedComponents(phys=PHY.FLUX_F, type="ELEM", components=("FLUX", "FLUY", "FLUZ"))

CFLUXNR = LocatedComponents(phys=PHY.FLUN_R, type="ELEM", components=("FLUN",))


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


CTEMPSR = LocatedComponents(
    phys=PHY.INST_R, type="ELEM", components=("INST", "DELTAT", "THETA", "KHI", "R", "RHO")
)


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


EMNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[30]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


CT_EXTF = LocatedComponents(phys=PHY.TEMP_F, type="ELEM", components=("TEMP",))


DDL_THER = LocatedComponents(
    phys=PHY.TEMP_R, type="ELNO", diff=True, components=(("EN1", ("HHO_F[3]",)), ("EN2", ()))
)

MVECTAR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=NACCELR)

MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)

# ------------------------------------------------------------


class THER3DQU9_HHO1_F(Element):
    """Please document this element"""

    meshType = MT.QUAD9
    nodes = (SetOfNodes("EN1", (9,)), SetOfNodes("EN2", (1, 2, 3, 4, 5, 6, 7, 8)))
    attrs = ((AT.BORD_ISO, "OUI"),)
    elrefe = (ElrefeLoc(MT.QU9, gauss=("RIGI=FPG4",), mater=("RIGI",)),)
    calculs = (
        OP.CHAR_THER_FLUN_F(
            te=461,
            para_in=((SP.PFLUXNF, CFLUXNF), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_R(
            te=461,
            para_in=((SP.PFLUXNR, CFLUXNR), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_TEXT_F(
            te=461,
            para_in=(
                (SP.PCOEFHF, CCOEFHF),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPER, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
                (SP.PT_EXTF, CT_EXTF),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_TEXT_R(
            te=461,
            para_in=(
                (SP.PCOEFHR, CCOEFHR),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPER, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
                (SP.PT_EXTR, LC.ET_EXTR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.RIGI_THER_COEH_F(
            te=457,
            para_in=((SP.PCOEFHF, CCOEFHF), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((OP.RIGI_THER_COEH_F.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_COEH_R(
            te=457,
            para_in=((SP.PCOEFHR, CCOEFHR), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((OP.RIGI_THER_COEH_R.PMATTTR, MMATTTR),),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PGEOM_R, LC.GGEOMER),
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
            ),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),
                (OP.TOU_INI_ELEM.PNEU1_R, LC.CNEUTR1),
                (OP.TOU_INI_ELEM.PCOEH_R, LC.EHECHPR),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
            ),
        ),
    )


# ------------------------------------------------------------


class THER3DTR7_HHO1_F(THER3DQU9_HHO1_F):
    """Please document this element"""

    meshType = MT.TRIA7
    nodes = (SetOfNodes("EN1", (7,)), SetOfNodes("EN2", (1, 2, 3, 4, 5, 6)))
    attrs = ((AT.BORD_ISO, "OUI"),)
    elrefe = (ElrefeLoc(MT.TR7, gauss=("RIGI=FPG3",), mater=("RIGI",)),)
