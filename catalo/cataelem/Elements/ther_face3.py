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


DDL_THER = LocatedComponents(phys=PHY.TEMP_R, type="ELNO", components=("TEMP",))


MVECTAR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=NACCELR)

MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)


# ------------------------------------------------------------
class THER_FACE3(Element):
    """Please document this element"""

    meshType = MT.TRIA3
    elrefe = (ElrefeLoc(MT.TR3, gauss=("RIGI=COT3", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),)
    calculs = (
        OP.ACCEPTANCE(
            te=329,
            para_in=((SP.PACCELR, NACCELR), (SP.PGEOMER, NGEOMER), (SP.PNUMMOD, LC.CNUMMOD)),
            para_out=((SP.PVECTUR, MVECTAR),),
        ),
        OP.AMOR_AJOU(
            te=327,
            para_in=((SP.PACCELR, NACCELR), (SP.PGEOMER, NGEOMER)),
            para_out=((OP.AMOR_AJOU.PMATTTR, MMATTTR),),
        ),
        OP.CHAR_THER_ACCE_R(
            te=325,
            para_in=((SP.PACCELR, NACCELR), (SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ACCE_X(
            te=325,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC), (SP.PTEMPER, DDL_THER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ACCE_Y(
            te=325,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC), (SP.PTEMPER, DDL_THER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ACCE_Z(
            te=325,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC), (SP.PTEMPER, DDL_THER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUNL(
            te=273,
            para_in=(
                (SP.PFLUXNL, CFLUXNF),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPER, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_F(
            te=58,
            para_in=((SP.PFLUXNF, CFLUXNF), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_R(
            te=57,
            para_in=((SP.PFLUXNR, CFLUXNR), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUTNL(
            te=526,
            para_in=(
                (SP.PFLUXNL, CFLUXNF),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPER, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.CHAR_THER_FLUX_F(
            te=95,
            para_in=((SP.PFLUXVF, CFLUXVF), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_PHID_R(
            te=326,
            para_in=(
                (SP.PACCELR, NACCELR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_RAYO_F(
            te=60,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PRAYONF, LC.CRAYONF),
                (SP.PTEMPER, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_RAYO_R(
            te=59,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PRAYONR, LC.CRAYONR),
                (SP.PTEMPER, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_TEXT_F(
            te=60,
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
            te=59,
            para_in=(
                (SP.PCOEFHR, CCOEFHR),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPER, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
                (SP.PT_EXTR, LC.ET_EXTR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.FLUX_FLUI_X(
            te=309, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.FLUX_FLUI_X.PMATTTR, MMATTTR),)
        ),
        OP.FLUX_FLUI_Y(
            te=309, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.FLUX_FLUI_Y.PMATTTR, MMATTTR),)
        ),
        OP.FLUX_FLUI_Z(
            te=309, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.FLUX_FLUI_Z.PMATTTR, MMATTTR),)
        ),
        OP.MTAN_THER_COEF_F(
            te=131,
            para_in=((SP.PCOEFHF, CCOEFHF), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((OP.MTAN_THER_COEF_F.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_COEF_R(
            te=130,
            para_in=((SP.PCOEFHR, CCOEFHR), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((OP.MTAN_THER_COEF_R.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_FLUXNL(
            te=132,
            para_in=(
                (SP.PFLUXNL, CFLUXNF),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((OP.MTAN_THER_FLUXNL.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_RAYO_F(
            te=131,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PRAYONF, LC.CRAYONF),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((OP.MTAN_THER_RAYO_F.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_RAYO_R(
            te=130,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PRAYONR, LC.CRAYONR),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((OP.MTAN_THER_RAYO_R.PMATTTR, MMATTTR),),
        ),
        OP.NORME_L2(
            te=563,
            para_in=(
                (SP.PCALCI, LC.EMNEUT_I),
                (SP.PCHAMPG, EGNEUT_R),
                (SP.PCOEFR, EMNEUT_R),
                (OP.NORME_L2.PCOORPG, EGGEOP_R),
            ),
            para_out=((SP.PNORME, LC.ENORME),),
        ),
        OP.RESI_THER_COEF_F(
            te=128,
            para_in=(
                (SP.PCOEFHF, CCOEFHF),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_COEF_R(
            te=127,
            para_in=(
                (SP.PCOEFHR, CCOEFHR),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_COEH_F(
            te=308,
            para_in=(
                (SP.PCOEFHF, CCOEFHF),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_COEH_R(
            te=306,
            para_in=(
                (SP.PCOEFHR, CCOEFHR),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_FLUXNL(
            te=129,
            para_in=(
                (SP.PFLUXNL, CFLUXNF),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_RAYO_F(
            te=128,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PRAYONF, LC.CRAYONF),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_RAYO_R(
            te=127,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PRAYONR, LC.CRAYONR),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RIGI_THER_COEH_F(
            te=53,
            para_in=((SP.PCOEFHF, CCOEFHF), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((OP.RIGI_THER_COEH_F.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_COEH_R(
            te=52,
            para_in=((SP.PCOEFHR, CCOEFHR), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((OP.RIGI_THER_COEH_R.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_COET_F(
            te=527,
            para_in=((SP.PCOEFHF, CCOEFHF), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((OP.RIGI_THER_COET_F.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_COET_R(
            te=523,
            para_in=((SP.PCOEFHR, CCOEFHR), (SP.PGEOMER, NGEOMER), (SP.PTEMPSR, CTEMPSR)),
            para_out=((OP.RIGI_THER_COET_R.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_FLUTNL(
            te=524,
            para_in=(
                (SP.PFLUXNL, CFLUXNF),
                (SP.PGEOMER, NGEOMER),
                (SP.PTEMPER, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
            ),
            para_out=((OP.RIGI_THER_FLUTNL.PMATTTR, MMATTTR),),
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
class THER_FACE4(THER_FACE3):
    """Please document this element"""

    meshType = MT.QUAD4
    elrefe = (ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),)


# ------------------------------------------------------------
class THER_FACE6(THER_FACE3):
    """Please document this element"""

    meshType = MT.TRIA6
    elrefe = (ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),)


# ------------------------------------------------------------
class THER_FACE8(THER_FACE3):
    """Please document this element"""

    meshType = MT.QUAD8
    elrefe = (ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),)


# ------------------------------------------------------------
class THER_FACE9(THER_FACE3):
    """Please document this element"""

    meshType = MT.QUAD9
    elrefe = (ElrefeLoc(MT.QU9, gauss=("RIGI=FPG9", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),)
