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

# ----------------------------------------------------------------------------------------------
# Located components
# ----------------------------------------------------------------------------------------------

DDL_MECA = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "PRES"))


MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

# ----------------------------------------------------------------------------------------------
class MEFSSE2P(Element):
    """Element for FSI interaction (U,P) - 2D - On SE2"""

    meshType = MT.SEG2
    elrefe = (ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2", "FPG1=FPG1"), mater=("FPG1",)),)
    calculs = (
        OP.CHAR_MECA_PRES_F(
            te=204,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PPRESSF, LC.CPRE2DF), (SP.PINSTR, LC.MTEMPSR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=204,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PPRESSR, LC.EPRE2DR), (SP.PINSTR, LC.MTEMPSR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_VFAC(
            te=384,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC), (SP.PVITEFR, LC.EVITEFR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_VFAC_F(
            te=384,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVITEFF, LC.EVITEFF),
                (SP.PINSTR, LC.MTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FORC_NODA(
            te=89,
            para_in=((SP.PDEPLAR, DDL_MECA), (SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=89,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PCODRET, LC.ECODRET), (SP.PMATUNS, MMATUNS), (SP.PVECTUR, MVECTUR)),
        ),
        OP.MASS_MECA(
            te=257,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.RAPH_MECA(
            te=89,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PCODRET, LC.ECODRET), (SP.PVECTUR, MVECTUR)),
        ),
        OP.RIGI_MECA(
            te=89,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.RIGI_MECA_TANG(
            te=89,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=(
                (SP.PMATUNS, MMATUNS),
                (SP.PVECTUR, MVECTUR),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),)),
        OP.TOU_INI_ELGA(te=99, para_out=((OP.TOU_INI_ELGA.PGEOM_R, LC.EGGAU2D),)),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM2D),)),
    )


# ----------------------------------------------------------------------------------------------
class MEFSSE3P(MEFSSE2P):
    """Element for FSI interaction (U,P,PHI) - 2D - On SE3"""

    meshType = MT.SEG3
    elrefe = (ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("FPG1",)),)
