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
        "JG1",
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
        "JG2",
        "TVAR",
    ),
)


CCAORIE = LocatedComponents(phys=PHY.CAORIE_R, type="ELEM", components=("ALPHA", "BETA", "GAMMA"))


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ", "DRX", "DRY", "DRZ")
)


NVITER = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))


CFORCEF = LocatedComponents(
    phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ", "MX", "MY", "MZ", "REP")
)


CFORCER = LocatedComponents(
    phys=PHY.FORC_R, type="ELEM", components=("FX", "FY", "FZ", "MX", "MY", "MZ", "REP")
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


CINSTMR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


EEFGEGC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELGA", location="RIGI", components=("N", "VY", "VZ", "MT", "MFY", "MFZ")
)


EEFGENC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("N", "VY", "VZ", "MT", "MFY", "MFZ")
)


EEFGEGA = LocatedComponents(
    phys=PHY.SIEF_R, type="ELGA", location="RIGI", components=("N", "VY", "VZ", "MT", "MFY", "MFZ")
)


EEFGENO = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("N", "VY", "VZ", "MT", "MFY", "MFZ")
)


CSTADYN = LocatedComponents(
    phys=PHY.STAOUDYN, type="ELEM", components=("STAOUDYN", "ALFNMK", "DELNMK")
)


ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MECA_POU_D_T_GD(Element):

    """Please document this element"""

    meshType = MT.SEG2
    elrefe = (
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG1", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.CABPOU),
    )
    calculs = (
        OP.ADD_SIGM(
            te=581,
            para_in=((SP.PEPCON1, EEFGEGA), (SP.PEPCON2, EEFGEGA)),
            para_out=((SP.PEPCON3, EEFGEGA),),
        ),
        OP.CHAR_MECA_FF1D1D(
            te=161,
            para_in=((SP.PFF1D1D, CFORCEF), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FR1D1D(
            te=161,
            para_in=((SP.PFR1D1D, CFORCER), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=161,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SF1D1D(
            te=161,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PFF1D1D, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SR1D1D(
            te=161,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PVENTCX, LC.CVENTCX),
                (SP.PVITER, NVITER),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=396,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.CHAR_MECA_TEMP_R.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.COOR_ELGA(
            te=478, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.EFGE_ELGA(
            te=546,
            para_in=((SP.PSIEFR, EEFGEGA),),
            para_out=((SP.PEFGEC, EEFGEGC), (SP.PEFGER, EEFGEGA)),
        ),
        OP.EFGE_ELNO(
            te=185,
            para_in=(
                (OP.EFGE_ELNO.PCONTRR, EEFGEGA),
                (SP.PNONLIN, LC.ENONLIN),
                (OP.EFGE_ELNO.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PEFFORC, EEFGENC), (OP.EFGE_ELNO.PEFFORR, EEFGENO)),
        ),
        OP.FORC_NODA(
            te=393,
            para_in=(
                (OP.FORC_NODA.PCAORIE, CCAORIE),
                (SP.PSIEFR, EEFGEGA),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=390,
            para_in=(
                (SP.PACCKM1, DDL_MECA),
                (SP.PACCPLU, DDL_MECA),
                (SP.PCAGNPO, CCAGNPO),
                (OP.FULL_MECA.PCAORIE, CCAORIE),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (SP.PDDEPLA, DDL_MECA),
                (SP.PDEPKM1, DDL_MECA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CINSTMR),
                (SP.PINSTPR, LC.CINSTPR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PROMK, DDL_MECA),
                (SP.PROMKM1, DDL_MECA),
                (SP.PSTADYN, CSTADYN),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (SP.PVITKM1, DDL_MECA),
                (SP.PVITPLU, DDL_MECA),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, EEFGEGA),
                (SP.PMATUNS, MMATUNS),
                (OP.FULL_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_MECA(
            te=391,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.MASS_MECA.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.M_GAMMA(
            te=391,
            para_in=(
                (SP.PACCELR, DDL_MECA),
                (SP.PCAGNPO, CCAGNPO),
                (OP.M_GAMMA.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.M_GAMMA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.PAS_COURANT(
            te=404,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.PAS_COURANT.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PCOURAN, LC.ECOURAN),),
        ),
        OP.RAPH_MECA(
            te=390,
            para_in=(
                (SP.PACCKM1, DDL_MECA),
                (SP.PACCPLU, DDL_MECA),
                (SP.PCAGNPO, CCAGNPO),
                (OP.RAPH_MECA.PCAORIE, CCAORIE),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (SP.PDDEPLA, DDL_MECA),
                (SP.PDEPKM1, DDL_MECA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CINSTMR),
                (SP.PINSTPR, LC.CINSTPR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PROMK, DDL_MECA),
                (SP.PROMKM1, DDL_MECA),
                (SP.PSTADYN, CSTADYN),
                (OP.RAPH_MECA.PVARIMR, ZVARIPG),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (SP.PVITKM1, DDL_MECA),
                (SP.PVITPLU, DDL_MECA),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, EEFGEGA),
                (OP.RAPH_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.REFE_FORC_NODA(
            te=393, para_in=((SP.PREFCO, LC.CRESEFM),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.REPERE_LOCAL(
            te=135,
            para_in=((OP.REPERE_LOCAL.PCAORIE, CCAORIE),),
            para_out=((SP.PREPLO1, LC.CGEOM3D), (SP.PREPLO2, LC.CGEOM3D), (SP.PREPLO3, LC.CGEOM3D)),
        ),
        OP.RIGI_MECA_TANG(
            te=390,
            para_in=(
                (SP.PACCKM1, DDL_MECA),
                (SP.PACCPLU, DDL_MECA),
                (SP.PCAGNPO, CCAGNPO),
                (OP.RIGI_MECA_TANG.PCAORIE, CCAORIE),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (SP.PDDEPLA, DDL_MECA),
                (SP.PDEPKM1, DDL_MECA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CINSTMR),
                (SP.PINSTPR, LC.CINSTPR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PROMK, DDL_MECA),
                (SP.PROMKM1, DDL_MECA),
                (SP.PSTADYN, CSTADYN),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (SP.PVITKM1, DDL_MECA),
                (SP.PVITPLU, DDL_MECA),
            ),
            para_out=(
                (SP.PMATUNS, MMATUNS),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, EEFGEGA),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.SIEF_ELNO(
            te=4,
            para_in=((OP.SIEF_ELNO.PCONTRR, EEFGEGA), (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PSIEFNOC, EEFGENC), (OP.SIEF_ELNO.PSIEFNOR, EEFGENO)),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
                (OP.TOU_INI_ELGA.PSIEF_R, EEFGEGA),
                (OP.TOU_INI_ELGA.PVARI_R, ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
        OP.VARI_ELNO(
            te=4, para_in=((SP.PVARIGR, ZVARIPG),), para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),)
        ),
    )
