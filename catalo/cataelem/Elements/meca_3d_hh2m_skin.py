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


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("DX", "DY", "DZ", "PRE[2]")), ("EN2", ("DX", "DY", "DZ"))),
)

CCOECH = LocatedComponents(
    phys=PHY.ETHM_R, type="ELGA", location="RIGI", components=("COEF[4]", "PRE[2]")
)
CCOECHF = LocatedComponents(
    phys=PHY.ETHM_F, type="ELGA", location="RIGI", components=("COEF[4]", "PRE[2]")
)
CCOECHH = LocatedComponents(
    phys=PHY.ETHMH_R, type="ELGA", location="RIGI", components=("COEF[2]", "HR[1]")
)
CCOECHHF = LocatedComponents(
    phys=PHY.ETHMH_F, type="ELGA", location="RIGI", components=("COEF[2]", "HR[1]")
)
EFLHN = LocatedComponents(
    phys=PHY.FLHN_R, type="ELGA", location="RIGI", components=("FH1[2]", "FH2[2]")
)


EFORCNO = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ"))


CFLUXF = LocatedComponents(phys=PHY.FTHM_F, type="ELEM", components=("PFLU[2]",))


EFLUXE = LocatedComponents(phys=PHY.FTHM_R, type="ELGA", location="RIGI", components=("PFLU[2]",))


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST", "DELTAT", "THETA"))


CPRESSF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PRES",))


EPRESNO = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES",))


NSIEF_R = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELNO",
    diff=True,
    components=(
        (
            "EN1",
            (
                "FH11X",
                "FH11Y",
                "FH11Z",
                "FH12X",
                "FH12Y",
                "FH12Z",
                "FH21X",
                "FH21Y",
                "FH21Z",
                "FH22X",
                "FH22Y",
                "FH22Z",
            ),
        ),
        ("EN2", ()),
    ),
)


ECONTNO = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
)


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class HH2M_FACE8(Element):

    """Please document this element"""

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9",)), ElrefeLoc(MT.QU4, gauss=("RIGI=FPG9",)))
    calculs = (
        OP.CHAR_MECA_FLUX_F(
            te=466,
            para_in=((SP.PFLUXF, CFLUXF), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_ECHA_THM_F(
            te=475,
            para_in=(
                (SP.PFLUXF, CFLUXF),
                (SP.PGEOMER, NGEOMER),
                (SP.PCHTHMF, CCOECHF),
                (SP.PINSTR, CTEMPSR),
                (SP.PDEPLMR, DDL_MECA),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_ECHA_HR_F(
            te=475,
            para_in=(
                (SP.PFLUXF, CFLUXF),
                (SP.PGEOMER, NGEOMER),
                (SP.HCHTHMF, CCOECHHF),
                (SP.PINSTR, CTEMPSR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_ECHA_THM_R(
            te=475,
            para_in=(
                (SP.PFLUXR, EFLUXE),
                (SP.PGEOMER, NGEOMER),
                (SP.PECHTHM, CCOECH),
                (SP.PINSTR, CTEMPSR),
                (SP.PDEPLMR, DDL_MECA),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_ECHA_HR_R(
            te=475,
            para_in=(
                (SP.PFLUXR, EFLUXE),
                (SP.PGEOMER, NGEOMER),
                (SP.HECHTHM, CCOECHH),
                (SP.PINSTR, CTEMPSR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FLUX_R(
            te=466,
            para_in=((SP.PFLUXR, EFLUXE), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FR2D3D(
            te=466,
            para_in=((SP.PFR2D3D, EFORCNO), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_F(
            te=466,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PPRESSF, CPRESSF), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=466,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PPRESSR, EPRESNO)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.FLHN_ELGA(
            te=493,
            para_in=((SP.PCONTR, NSIEF_R), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PFLHN, EFLHN),),
        ),
        OP.SIRO_ELEM(
            te=411,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PSIG3D, ECONTNO)),
            para_out=((SP.PPJSIGM, LC.EPJSIGM),),
        ),
        OP.TOU_INI_ELGA(te=99, para_out=((OP.TOU_INI_ELGA.PGEOM_R, EGGEOP_R),)),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
    )


# ------------------------------------------------------------
class HH2M_FACE6(HH2M_FACE8):

    """Please document this element"""

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN2", (4, 5, 6)), SetOfNodes("EN1", (1, 2, 3)))
    elrefe = (ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6",)), ElrefeLoc(MT.TR3, gauss=("RIGI=FPG6",)))
