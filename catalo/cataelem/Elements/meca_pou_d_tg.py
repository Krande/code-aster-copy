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


CABSCUR = LocatedComponents(phys=PHY.ABSC_R, type="ELEM", components=("ABSC[2]",))


CCAGEPO = LocatedComponents(phys=PHY.CAGEPO_R, type="ELEM", components=("R1", "EP1"))


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


NDEPLAC = LocatedComponents(
    phys=PHY.DEPL_C, type="ELNO", components=("DX", "DY", "DZ", "DRX", "DRY", "DRZ", "GRX")
)


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ", "DRX", "DRY", "DRZ", "GRX")
)


NVITER = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))


CEPSINR = LocatedComponents(
    phys=PHY.EPSI_R, type="ELEM", components=("EPX", "KY", "KZ", "VECT_2_X", "VECT_2_Y", "VECT_2_Z")
)


CEPSINF = LocatedComponents(
    phys=PHY.EPSI_F, type="ELEM", components=("EPX", "KY", "KZ", "VECT_2_X", "VECT_2_Y", "VECT_2_Z")
)


EDEFGNO = LocatedComponents(
    phys=PHY.EPSI_R, type="ELNO", components=("EPXX", "GAXY", "GAXZ", "GAT", "KY", "KZ", "GAX")
)


CFORCEF = LocatedComponents(
    phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ", "MX", "MY", "MZ", "REP")
)


CFORCER = LocatedComponents(
    phys=PHY.FORC_R, type="ELEM", components=("FX", "FY", "FZ", "MX", "MY", "MZ", "REP")
)


EGGEOM_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z")
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


ECONTPC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("SN", "SVY", "SVZ", "SMT", "SMFY", "SMFZ")
)


EEFGEGC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELGA",
    location="RIGI",
    components=("N", "VY", "VZ", "MT", "MFY", "MFZ", "BX"),
)


EEFGENC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("N", "VY", "VZ", "MT", "MFY", "MFZ", "BX")
)


EEFGEGA = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=("N", "VY", "VZ", "MT", "MFY", "MFZ", "BX"),
)


EEFGENO = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("N", "VY", "VZ", "MT", "MFY", "MFZ", "BX")
)


ESTRAUX = LocatedComponents(
    phys=PHY.STRX_R, type="ELGA", location="RIGI", components=("ALPHA", "BETA", "GAMMA")
)


ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))


MVECTUC = ArrayOfComponents(phys=PHY.VDEP_C, locatedComponents=NDEPLAC)

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MECA_POU_D_TG(Element):

    """Please document this element"""

    meshType = MT.SEG2
    elrefe = (ElrefeLoc(MT.SE2, gauss=("RIGI=FPG3", "FPG1=FPG1"), mater=("RIGI", "FPG1")),)
    calculs = (
        OP.ADD_SIGM(
            te=581,
            para_in=((SP.PEPCON1, EEFGEGA), (SP.PEPCON2, EEFGEGA)),
            para_out=((SP.PEPCON3, EEFGEGA),),
        ),
        OP.AMOR_MECA(
            te=121,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMASSEL, MMATUUR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PRIGIEL, MMATUUR),
                (OP.AMOR_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.CHAR_MECA_EPSI_R(
            te=20,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.CHAR_MECA_EPSI_R.PCAORIE, CCAORIE),
                (SP.PEPSINR, CEPSINR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_EPSI_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_EPSI_F(
            te=20,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.CHAR_MECA_EPSI_F.PCAORIE, CCAORIE),
                (SP.PEPSINF, CEPSINF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_EPSI_F.PVARCPR, LC.ZVARCPG),
                (SP.PINSTR, CTEMPSR),
                (SP.PABSCUR, CABSCUR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FC1D1D(
            te=150,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.CHAR_MECA_FC1D1D.PCAORIE, CCAORIE),
                (SP.PFC1D1D, LC.CFORCEC),
                (SP.PGEOMER, NGEOMER),
            ),
            para_out=((SP.PVECTUC, MVECTUC),),
        ),
        OP.CHAR_MECA_FF1D1D(
            te=150,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.CHAR_MECA_FF1D1D.PCAORIE, CCAORIE),
                (SP.PFF1D1D, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FR1D1D(
            te=150,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.CHAR_MECA_FR1D1D.PCAORIE, CCAORIE),
                (SP.PFR1D1D, CFORCER),
                (SP.PGEOMER, NGEOMER),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FRELEC(
            te=145,
            para_in=((SP.PFRELEC, LC.CFRELEC), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_HYDR_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_HYDR_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=150,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.CHAR_MECA_PESA_R.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_ROTA_R(
            te=150,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.CHAR_MECA_ROTA_R.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PROTATR, LC.CROTATR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SECH_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_HYDR_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SF1D1D(
            te=150,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.CHAR_MECA_SF1D1D.PCAORIE, CCAORIE),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PFF1D1D, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SR1D1D(
            te=150,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.CHAR_MECA_SR1D1D.PCAORIE, CCAORIE),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PVENTCX, LC.CVENTCX),
                (SP.PVITER, NVITER),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=150,
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
        OP.DEGE_ELNO(
            te=158,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.DEGE_ELNO.PCAORIE, CCAORIE),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.DEGE_ELNO.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PDEFOGR, EDEFGNO),),
        ),
        OP.EFGE_ELGA(
            te=546,
            para_in=((SP.PSIEFR, EEFGEGA),),
            para_out=((SP.PEFGEC, EEFGEGC), (SP.PEFGER, EEFGEGA)),
        ),
        OP.EFGE_ELNO(
            te=185,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.EFGE_ELNO.PCAORIE, CCAORIE),
                (OP.EFGE_ELNO.PCOMPOR, LC.CCOMPOR),
                (OP.EFGE_ELNO.PCONTRR, EEFGEGA),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PFF1D1D, CFORCEF),
                (SP.PFR1D1D, CFORCER),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PNONLIN, LC.ENONLIN),
                (SP.PPESANR, LC.CPESANR),
                (SP.PINSTR, CTEMPSR),
                (OP.EFGE_ELNO.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PEFFORC, EEFGENC), (OP.EFGE_ELNO.PEFFORR, EEFGENO)),
        ),
        OP.EFEQ_ELNO(
            te=83, para_in=((SP.PEFFONR, EEFGENO),), para_out=((SP.PEFFOENR, LC.EEFGENOQ),)
        ),
        OP.EPOT_ELEM(
            te=151,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.EPOT_ELEM.PCAORIE, CCAORIE),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPOT_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPOT_ELEM.PENERDR, LC.EENEDNO),),
        ),
        OP.FORC_NODA(
            te=347,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.FORC_NODA.PCAORIE, CCAORIE),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, EEFGEGA),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PSTRXMR, ESTRAUX),
                (SP.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=346,
            para_in=(
                (SP.PCAGNPO, LC.CCAGNP1),
                (OP.FULL_MECA.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, EEFGEGA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PSTRXMR, ESTRAUX),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, EEFGEGA),
                (SP.PMATUUR, MMATUUR),
                (SP.PSTRXPR, ESTRAUX),
                (OP.FULL_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.INI_STRX(
            te=23, para_in=((OP.INI_STRX.PCAORIE, CCAORIE),), para_out=((SP.PSTRX_R, ESTRAUX),)
        ),
        OP.MASS_FLUI_STRU(
            te=141,
            para_in=(
                (SP.PABSCUR, CABSCUR),
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAGNPO, CCAGNPO),
                (OP.MASS_FLUI_STRU.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_FLUI_STRU.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_INER(
            te=38,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.MASS_INER.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_INER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMASSINE, LC.EMASSINE),),
        ),
        OP.MASS_MECA(
            te=141,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.MASS_MECA.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_MECA_DIAG(
            te=141,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.MASS_MECA_DIAG.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA_DIAG.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_MECA_EXPLI(
            te=141,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.MASS_MECA_EXPLI.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA_EXPLI.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MECA_GYRO(
            te=259,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.MECA_GYRO.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.M_GAMMA(
            te=141,
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
            te=346,
            para_in=(
                (SP.PCAGNPO, LC.CCAGNP1),
                (OP.RAPH_MECA.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, EEFGEGA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PSTRXMR, ESTRAUX),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, EEFGEGA),
                (SP.PSTRXPR, ESTRAUX),
                (OP.RAPH_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.REFE_FORC_NODA(
            te=347, para_in=((SP.PREFCO, LC.CRESEFM),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.REPERE_LOCAL(
            te=135,
            para_in=((OP.REPERE_LOCAL.PCAORIE, CCAORIE),),
            para_out=((SP.PREPLO1, LC.CGEOM3D), (SP.PREPLO2, LC.CGEOM3D), (SP.PREPLO3, LC.CGEOM3D)),
        ),
        OP.RIGI_FLUI_STRU(
            te=140,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.RIGI_FLUI_STRU.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_FLUI_STRU.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA(
            te=140,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.RIGI_MECA.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_GEOM(
            te=143,
            para_in=(
                (SP.PCAGNPO, LC.CCAGNP1),
                (OP.RIGI_GEOM.PCAORIE, CCAORIE),
                (OP.RIGI_GEOM.PEFFORR, EEFGEGA),
                (SP.PGEOMER, NGEOMER),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_HYST(
            te=50,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PRIGIEL, MMATUUR),
                (OP.RIGI_MECA_HYST.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUC, MMATUUC),),
        ),
        OP.RIGI_MECA_RO(
            te=235,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.RIGI_MECA_RO.PCAORIE, CCAORIE),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PROTATR, LC.CROTATR),
                (OP.RIGI_MECA_RO.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_TANG(
            te=346,
            para_in=(
                (SP.PCAGNPO, LC.CCAGNP1),
                (OP.RIGI_MECA_TANG.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, EEFGEGA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PSTRXMR, ESTRAUX),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, EEFGEGA),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.SIEF_ELGA(
            te=342,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.SIEF_ELGA.PCAORIE, CCAORIE),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.SIEF_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PCONTRC, EEFGEGC), (OP.SIEF_ELGA.PCONTRR, EEFGEGA)),
        ),
        OP.SIEF_ELNO(
            te=347,
            para_in=(
                (SP.PCAGNPO, CCAGNPO),
                (OP.SIEF_ELNO.PCAORIE, CCAORIE),
                (OP.SIEF_ELNO.PCOMPOR, LC.CCOMPOR),
                (OP.SIEF_ELNO.PCONTRR, EEFGEGA),
                (SP.PDEPPLU, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PSIEFNOC, EEFGENC), (OP.SIEF_ELNO.PSIEFNOR, EEFGENO)),
        ),
        OP.SIPM_ELNO(
            te=344,
            para_in=(
                (SP.PCAGEPO, LC.CCAGRPO),
                (SP.PCAGNPO, CCAGNPO),
                (OP.SIPM_ELNO.PCAORIE, CCAORIE),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PFF1D1D, CFORCEF),
                (SP.PFR1D1D, CFORCER),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PINSTR, CTEMPSR),
                (OP.SIPM_ELNO.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PSIMXRC, LC.ESIMXNC), (SP.PSIMXRR, LC.ESIMXNO)),
        ),
        OP.SIPO_ELNO(
            te=344,
            para_in=(
                (SP.PCAGEPO, LC.CCAGRPO),
                (SP.PCAGNPO, CCAGNPO),
                (OP.SIPO_ELNO.PCAORIE, CCAORIE),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PFF1D1D, CFORCEF),
                (SP.PFR1D1D, CFORCER),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.SIPO_ELNO.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PCONTPC, ECONTPC), (SP.PCONTPO, LC.ECONTPO)),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PGEOM_R, EGGEOM_R),
                (OP.TOU_INI_ELGA.PINST_R, LC.EGINST_R),
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
                (OP.TOU_INI_ELGA.PSIEF_R, EEFGEGA),
                (OP.TOU_INI_ELGA.PVARI_R, ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),
                (OP.TOU_INI_ELNO.PINST_R, LC.EEINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.EENEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.EENEUT_R),
                (OP.TOU_INI_ELNO.PSIEF_R, EEFGENO),
                (OP.TOU_INI_ELNO.PVARI_R, LC.ZVARINO),
            ),
        ),
        OP.VARI_ELNO(
            te=347,
            para_in=((OP.VARI_ELNO.PCOMPOR, LC.CCOMPOR), (SP.PVARIGR, ZVARIPG)),
            para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),),
        ),
    )
