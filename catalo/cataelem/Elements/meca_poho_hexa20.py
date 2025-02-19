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


CCAGNPO = LocatedComponents(
    phys=PHY.CAGNPO_R,
    type="ELEM",
    components=(
        "A1",
        "IY1",
        "IZ1",
        "AY1",
        "AZ1",
        "EY1",
        "EZ1",
        "JX1",
        "RY1",
        "RZ1",
        "RT1",
        "A2",
        "IY2",
        "IZ2",
        "AY2",
        "AZ2",
        "EY2",
        "EZ2",
        "JX2",
        "RY2",
        "RZ2",
        "RT2",
        "TVAR",
    ),
)


CCAORIE = LocatedComponents(phys=PHY.CAORIE_R, type="ELEM", components=("ALPHA", "BETA", "GAMMA"))


CCAPOUF = LocatedComponents(
    phys=PHY.CAPOUF_R,
    type="ELEM",
    components=("B_T", "B_N", "B_TN", "A_FLUI", "A_CELL", "COEF_ECH"),
)


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("DX", "DY", "DZ", "DRX", "DRY", "DRZ", "PHI")), ("EN2", ("PHI",))),
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MECA_POHO_HEXA20(Element):
    """Please document this element"""

    meshType = MT.HEXA20
    nodes = (
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
    )
    elrefe = (
        ElrefeLoc(MT.H20, gauss=("RIGI=FPG27", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.POHOH20),
    )
    calculs = (
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.MASS_INER(
            te=65,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (SP.PCAPOUF, CCAPOUF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PMASSINE, LC.EMASSINE),),
        ),
        OP.MASS_MECA(
            te=470,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.MASS_MECA.PCAORIE, CCAORIE),
                (SP.PCAPOUF, CCAPOUF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.REPERE_LOCAL(
            te=135,
            para_in=((OP.REPERE_LOCAL.PCAORIE, CCAORIE),),
            para_out=((SP.PREPLO1, LC.CGEOM3D), (SP.PREPLO2, LC.CGEOM3D), (SP.PREPLO3, LC.CGEOM3D)),
        ),
        OP.RIGI_MECA(
            te=471,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.RIGI_MECA.PCAORIE, CCAORIE),
                (SP.PCAPOUF, CCAPOUF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELGA(te=99, para_out=((OP.TOU_INI_ELGA.PGEOM_R, EGGEOP_R),)),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# ------------------------------------------------------------
class MECA_POHO_HEXA8(MECA_POHO_HEXA20):
    """Please document this element"""

    meshType = MT.HEXA8
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),)
    elrefe = (
        ElrefeLoc(MT.HE8, gauss=("RIGI=FPG8", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.POHOH8),
    )
