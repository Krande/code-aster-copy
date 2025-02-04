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


CCACOQU = LocatedComponents(phys=PHY.CACOQU_R, type="ELEM", components=("EP", "ALPHA", "BETA"))


CCOEFHF = LocatedComponents(phys=PHY.COEH_F, type="ELEM", components=("H_INF", "H_SUP"))


CCOEFHR = LocatedComponents(phys=PHY.COEH_R, type="ELEM", components=("H_INF", "H_SUP"))


CFLUXNF = LocatedComponents(phys=PHY.FLUN_F, type="ELEM", components=("FLUN_INF", "FLUN_SUP"))


CFLUXNR = LocatedComponents(phys=PHY.FLUN_R, type="ELEM", components=("FLUN_INF", "FLUN_SUP"))


EFLUXPG = LocatedComponents(
    phys=PHY.FLUX_R,
    type="ELGA",
    location="RIGI",
    components=(
        "FLUX",
        "FLUY",
        "FLUZ",
        "FLUX_INF",
        "FLUY_INF",
        "FLUZ_INF",
        "FLUX_SUP",
        "FLUY_SUP",
        "FLUZ_SUP",
    ),
)


EFLUXNO = LocatedComponents(
    phys=PHY.FLUX_R,
    type="ELNO",
    components=(
        "FLUX",
        "FLUY",
        "FLUZ",
        "FLUX_INF",
        "FLUY_INF",
        "FLUZ_INF",
        "FLUX_SUP",
        "FLUY_SUP",
        "FLUZ_SUP",
    ),
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOM_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z")
)


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


ENGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST", "DELTAT", "THETA"))


ENBSP_I = LocatedComponents(phys=PHY.NBSP_I, type="ELEM", components=("COQ_NCOU",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


CT_EXTF = LocatedComponents(
    phys=PHY.TEMP_F, type="ELEM", components=("TEMP", "TEMP_INF", "TEMP_SUP")
)


DDL_THER = LocatedComponents(
    phys=PHY.TEMP_R, type="ELNO", components=("TEMP_MIL", "TEMP_INF", "TEMP_SUP")
)


MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)


# ------------------------------------------------------------
class THCOQU4(Element):
    """Please document this element"""

    meshType = MT.QUAD4
    elrefe = (
        ElrefeLoc(
            MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG4", "NOEU_S=NOEU_S", "FPG1=FPG1"), mater=("FPG1",)
        ),
    )
    calculs = (
        OP.CHAR_THER_EVOL(
            te=110,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCOEFHF, CCOEFHF),
                (SP.PCOEFHR, CCOEFHR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_THER_EVOL.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_F(
            te=105,
            para_in=((SP.PFLUXNF, CFLUXNF), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_R(
            te=106,
            para_in=((SP.PFLUXNR, CFLUXNR), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ECHA_F(
            te=107,
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
            te=108,
            para_in=(
                (SP.PCOEFHR, CCOEFHR),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, CTEMPSR),
                (SP.PT_EXTR, LC.CT_EXTR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.COOR_ELGA(
            te=479, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.FLUX_ELGA(
            te=109,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.FLUX_ELGA.PNBSP_I, ENBSP_I),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, CTEMPSR),
                (OP.FLUX_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.FLUX_ELGA.PFLUXPG, EFLUXPG),),
        ),
        OP.FLUX_ELNO(
            te=4, para_in=((OP.FLUX_ELNO.PFLUXPG, EFLUXPG),), para_out=((SP.PFLUXNO, EFLUXNO),)
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_THER(
            te=102,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.MASS_THER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.MASS_THER.PMATTTR, MMATTTR),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.RIGI_THER(
            te=101,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.RIGI_THER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.RIGI_THER.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_ECHA_F(
            te=103,
            para_in=((SP.PCOEFHF, CCOEFHF), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((OP.RIGI_THER_ECHA_F.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_ECHA_R(
            te=104,
            para_in=((SP.PCOEFHR, CCOEFHR), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((OP.RIGI_THER_ECHA_R.PMATTTR, MMATTTR),),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PNBSP_I, ENBSP_I),)),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PGEOM_R, EGGEOM_R),
                (OP.TOU_INI_ELGA.PINST_R, LC.EGINST_R),
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, ENGEOM_R),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
            ),
        ),
    )


# ------------------------------------------------------------
class THCOQU8(THCOQU4):
    """Please document this element"""

    meshType = MT.QUAD8
    elrefe = (
        ElrefeLoc(
            MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9", "NOEU_S=NOEU_S", "FPG1=FPG1"), mater=("FPG1",)
        ),
    )


# ------------------------------------------------------------
class THCOQU9(THCOQU4):
    """Please document this element"""

    meshType = MT.QUAD9
    elrefe = (
        ElrefeLoc(
            MT.QU9, gauss=("RIGI=FPG9", "MASS=FPG9", "NOEU_S=NOEU_S", "FPG1=FPG1"), mater=("FPG1",)
        ),
    )


# ------------------------------------------------------------
class THCOTR3(THCOQU4):
    """Please document this element"""

    meshType = MT.TRIA3
    elrefe = (
        ElrefeLoc(
            MT.TR3, gauss=("RIGI=FPG1", "MASS=FPG3", "NOEU_S=NOEU_S", "FPG1=FPG1"), mater=("FPG1",)
        ),
    )


# ------------------------------------------------------------
class THCOTR6(THCOQU4):
    """Please document this element"""

    meshType = MT.TRIA6
    elrefe = (
        ElrefeLoc(
            MT.TR6, gauss=("RIGI=FPG3", "MASS=FPG6", "NOEU_S=NOEU_S", "FPG1=FPG1"), mater=("FPG1",)
        ),
    )


# ------------------------------------------------------------
class THCOTR7(THCOQU4):
    """Please document this element"""

    meshType = MT.TRIA7
    elrefe = (
        ElrefeLoc(
            MT.TR7, gauss=("RIGI=FPG6", "MASS=FPG6", "NOEU_S=NOEU_S", "FPG1=FPG1"), mater=("FPG1",)
        ),
    )
