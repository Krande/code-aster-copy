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

# Ce catalogue correspondant aux elements de surface libre fluide (2D)
# pour le couplage fluide-structure n'ont que des DDL de :
#  deplacement vertical DH et potentiel de deplacement PHI.
# La contrainte n'existe pas, ni la deformation.
# Les options FULL_MECA et RAPH_MECA sont disponibles mais correspondent
# a un probleme lineaire : le CODRET sera toujours OK sur ces elements.
# On le garde neanmoins pour des raisons de compatibilite.


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

DDL_MECA = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("PHI", "DH"))

NDEPLAC = LocatedComponents(phys=PHY.DEPL_C, type="ELNO", components=("PHI", "DH"))

NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MEFP_FACE3(Element):

    """Please document this element"""

    meshType = MT.TRIA3
    elrefe = (ElrefeLoc(MT.TR3, gauss=("RIGI=COT3", "FPG1=FPG1"), mater=("FPG1",)),)
    calculs = (
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.FORC_NODA(
            te=370,
            para_in=((SP.PDEPLAR, DDL_MECA), (SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=370,
            para_in=(
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PCODRET, LC.ECODRET), (SP.PMATUUR, MMATUUR), (SP.PVECTUR, MVECTUR)),
        ),
        OP.MASS_MECA(
            te=371,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.PAS_COURANT(
            te=405,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PCOURAN, LC.ECOURAN),),
        ),
        OP.RAPH_MECA(
            te=370,
            para_in=(
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PCODRET, LC.ECODRET), (SP.PVECTUR, MVECTUR)),
        ),
        OP.RIGI_MECA(
            te=370,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_HYST(
            te=370,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PMATUUC, MMATUUC),),
        ),
        OP.RIGI_MECA_TANG(
            te=370,
            para_in=(
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=(
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELGA(te=99, para_out=((OP.TOU_INI_ELGA.PGEOM_R, EGGEOP_R),)),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
    )


# ------------------------------------------------------------
class MEFP_FACE4(MEFP_FACE3):

    """Please document this element"""

    meshType = MT.QUAD4
    elrefe = (ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("FPG1",)),)


# ------------------------------------------------------------
class MEFP_FACE6(MEFP_FACE3):

    """Please document this element"""

    meshType = MT.TRIA6
    elrefe = (ElrefeLoc(MT.TR6, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("FPG1",)),)


# ------------------------------------------------------------
class MEFP_FACE8(MEFP_FACE3):

    """Please document this element"""

    meshType = MT.QUAD8
    elrefe = (ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "FPG1=FPG1"), mater=("FPG1",)),)


# ------------------------------------------------------------
class MEFP_FACE9(MEFP_FACE3):

    """Please document this element"""

    meshType = MT.QUAD9
    elrefe = (ElrefeLoc(MT.QU9, gauss=("RIGI=FPG9", "FPG1=FPG1"), mater=("FPG1",)),)
