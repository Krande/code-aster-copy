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

# ----------------------------------------------------------------------------------------------
# Located components
# ----------------------------------------------------------------------------------------------


NDEPLAC = LocatedComponents(phys=PHY.DEPL_C, type="ELNO", components=("DX", "DY", "DZ"))


DDL_MECA = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))

EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MEAB_FACE3(Element):
    """Please document this element"""

    meshType = MT.TRIA3
    elrefe = (
        ElrefeLoc(
            MT.TR3, gauss=("RIGI=COT3", "FPG1=FPG1", "MTGA=FPG1"), mater=("RIGI", "FPG1", "MTGA")
        ),
    )
    calculs = (
        OP.AMOR_MECA(
            te=569,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.AMOR_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.FORC_NODA(
            te=569,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PDEPLAR, DDL_MECA),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=569,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR), (SP.PVECTUR, MVECTUR)),
        ),
        OP.INIT_VARC(
            te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG), (OP.INIT_VARC.PVARCNO, LC.ZVARCNO))
        ),
        OP.MATE_ELGA(
            te=142,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((OP.MATE_ELGA.PMATERR, LC.EGMATE_R),),
        ),
        OP.MATE_ELEM(
            te=142,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((OP.MATE_ELEM.PMATERR, LC.EEMATE_R),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.ONDE_PLAN(
            te=498,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PONDPLA, LC.CONDPLA),
                (SP.PONDPLR, LC.CONDPLR),
                (SP.PINSTR, CTEMPSR),
                (OP.ONDE_PLAN.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RAPH_MECA(
            te=569,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (OP.RAPH_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_MECA(
            te=569,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_HYST(
            te=50,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC), (SP.PRIGIEL, MMATUUR)),
            para_out=((SP.PMATUUC, MMATUUC),),
        ),
        OP.RIGI_MECA_TANG(
            te=569,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
            ),
            para_out=(
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.RIGI_MECA_ELAS(
            te=569,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA_ELAS.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=((OP.TOU_INI_ELGA.PGEOM_R, EGGEOP_R), (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R)),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
        OP.VARC_ELGA(
            te=530,
            para_in=((OP.VARC_ELGA.PVARCPR, LC.ZVARCPG),),
            para_out=((SP.PVARC_R, LC.EVARC_R),),
        ),
    )


# ------------------------------------------------------------
class MEAB_FACE4(MEAB_FACE3):
    """Please document this element"""

    meshType = MT.QUAD4
    elrefe = (
        ElrefeLoc(
            MT.QU4, gauss=("RIGI=FPG4", "FPG1=FPG1", "MTGA=FPG1"), mater=("RIGI", "FPG1", "MTGA")
        ),
    )


# ------------------------------------------------------------
class MEAB_FACE6(MEAB_FACE3):
    """Please document this element"""

    meshType = MT.TRIA6
    elrefe = (
        ElrefeLoc(
            MT.TR6, gauss=("RIGI=FPG6", "FPG1=FPG1", "MTGA=FPG1"), mater=("RIGI", "FPG1", "MTGA")
        ),
    )


# ------------------------------------------------------------
class MEAB_FACE8(MEAB_FACE3):
    """Please document this element"""

    meshType = MT.QUAD8
    elrefe = (
        ElrefeLoc(
            MT.QU8, gauss=("RIGI=FPG9", "FPG1=FPG1", "MTGA=FPG1"), mater=("RIGI", "FPG1", "MTGA")
        ),
    )


# ------------------------------------------------------------
class MEAB_FACE9(MEAB_FACE3):
    """Please document this element"""

    meshType = MT.QUAD9
    elrefe = (
        ElrefeLoc(
            MT.QU9, gauss=("RIGI=FPG9", "FPG1=FPG1", "MTGA=FPG1"), mater=("RIGI", "FPG1", "MTGA")
        ),
    )
