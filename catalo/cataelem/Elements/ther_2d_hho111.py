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


# Reuse PHY.COMPOR from mechanical => names of components are strange
CCOMPOR = LocatedComponents(
    phys=PHY.COMPOR, type="ELEM", components=("RELCOM", "NBVARI", "DEFORM", "INCELA", "C_PLAN")
)


NVITESR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY"))


EENERR = LocatedComponents(phys=PHY.ENER_R, type="ELEM", components=("TOTALE",))


CGRAINF = LocatedComponents(phys=PHY.FLUX_F, type="ELEM", components=("FLUX", "FLUY"))


CGRAINR = LocatedComponents(phys=PHY.FLUX_R, type="ELEM", components=("FLUX", "FLUY"))


EFLUXPG = LocatedComponents(
    phys=PHY.FLUX_R, type="ELGA", location="RIGI", components=("FLUX", "FLUY")
)


EFLUXNO = LocatedComponents(phys=PHY.FLUX_R, type="ELNO", components=("FLUX", "FLUY"))


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "W")
)


EGGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y"))


ENGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


CTEMPSR = LocatedComponents(
    phys=PHY.INST_R, type="ELEM", components=("INST", "DELTAT", "THETA", "KHI", "R", "RHO")
)


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


ECASECT = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[9]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


EMNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[30]",))


ESOURCR = LocatedComponents(phys=PHY.SOUR_R, type="ELGA", location="RIGI", components=("SOUR",))

PFONC = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z[2]",))

CINSTR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))

DDL_THER = LocatedComponents(
    phys=PHY.TEMP_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("HHO_T[2]",)), ("EN2", ()), ("EN3", ("HHO_T[3]"))),
)

TEMPHHO = LocatedComponents(phys=PHY.TEMP_R, type="ELNO", components=("TEMP",))


MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)

MMATTSR = ArrayOfComponents(phys=PHY.MTNS_R, locatedComponents=DDL_THER)


# ------------------------------------------------------------


class THER2DQ9_HHO111(Element):
    """Please document this element"""

    meshType = MT.QUAD9
    nodes = (
        SetOfNodes("EN1", (5, 6, 7, 8)),
        SetOfNodes("EN2", (1, 2, 3, 4)),
        SetOfNodes("EN3", (9,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.QU9, gauss=("RIGI=FPG4", "FPG1=FPG1", "MASS=FPG4"), mater=("RIGI", "FPG1", "MASS")
        ),
    )
    calculs = (
        OP.CHAR_THER_EVOL(
            te=445,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PTEMPSR, CTEMPSR),
                (OP.CHAR_THER_EVOL.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOUR_F(
            te=465,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PSOURCF, LC.CSOURCF),
                (SP.PTEMPSR, CTEMPSR),
                (OP.CHAR_THER_SOUR_F.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOUR_R(
            te=465,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PSOURCR, ESOURCR), (SP.PTEMPSR, CTEMPSR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.HHO_TEMP_THER(
            te=456,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PTMPCHF, DDL_THER)),
            para_out=((OP.HHO_TEMP_THER.PTEMP_R, TEMPHHO),),
        ),
        OP.HHO_PROJ_THER(
            te=473,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (OP.HHO_PROJ_THER.PFUNC_R, PFONC),
                (SP.PINSTPR, CINSTR),
            ),
            para_out=((OP.HHO_PROJ_THER.PTEMP_R, DDL_THER),),
        ),
        OP.MASS_THER(
            te=449,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPSR, CTEMPSR),
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
            te=454,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPSR, CTEMPSR),
                (OP.RIGI_THER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.RIGI_THER.PMATTTR, MMATTTR),),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PCOEH_R, LC.EHECHPR), (OP.TOU_INI_ELEM.PSOUR_R, LC.CSOURCR)),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PFLUX_R, EFLUXPG),
                (OP.TOU_INI_ELGA.PGEOM_R, EGGEOM_R),
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGNEUT1R),
                (OP.TOU_INI_ELGA.PSOUR_R, ESOURCR),
                (OP.TOU_INI_ELGA.PVARI_R, LC.ZVARIPG),
                (SP.PTEMP_R, LC.ETEMPPG),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PFLUX_R, EFLUXNO),
                (OP.TOU_INI_ELNO.PGEOM_R, ENGEOM_R),
                (OP.TOU_INI_ELNO.PHYDRPM, LC.EHYDRNO),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PVARI_R, LC.EPHASNO_),
            ),
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# ------------------------------------------------------------


class THER2DT7_HHO111(THER2DQ9_HHO111):
    """Please document this element"""

    meshType = MT.TRIA7
    nodes = (SetOfNodes("EN1", (4, 5, 6)), SetOfNodes("EN2", (1, 2, 3)), SetOfNodes("EN3", (7,)))
    elrefe = (
        ElrefeLoc(
            MT.TR7, gauss=("RIGI=FPG3", "FPG1=FPG1", "MASS=FPG3"), mater=("RIGI", "FPG1", "MASS")
        ),
    )
