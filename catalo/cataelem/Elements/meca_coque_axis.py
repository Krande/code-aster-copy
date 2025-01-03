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


CCACOQU = LocatedComponents(phys=PHY.CACOQU_R, type="ELEM", components=("EP", "KAPPA", "C_METR"))


NDEPLAC = LocatedComponents(phys=PHY.DEPL_C, type="ELNO", components=("DX", "DY", "DRZ"))


DDL_MECA = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DRZ"))


EDEFOPC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELGA", location="RIGI", components=("EPXX", "EPYY", "EPXZ")
)


EDEFONC = LocatedComponents(phys=PHY.EPSI_C, type="ELNO", components=("EPXX", "EPYY", "EPXZ"))


EDEFOPG = LocatedComponents(
    phys=PHY.EPSI_R, type="ELGA", location="RIGI", components=("EPXX", "EPYY", "EPXZ")
)


EDEFONO = LocatedComponents(phys=PHY.EPSI_R, type="ELNO", components=("EPXX", "EPYY", "EPXZ"))


EDEFGNO = LocatedComponents(
    phys=PHY.EPSI_R, type="ELNO", components=("EXX", "EYY", "EXY", "KXX", "KYY", "KXY")
)


EDEFGPG = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EXX", "EYY", "EXY", "KXX", "KYY", "KXY"),
)


CFORCEF = LocatedComponents(
    phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ", "MY", "MZ", "REP")
)


EFORCNO = LocatedComponents(
    phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ", "MY", "MZ", "REP")
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "W")
)


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


ENBSP_I = LocatedComponents(phys=PHY.NBSP_I, type="ELEM", components=("COQ_NCOU",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[6]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[6]",))


CPRESSF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PRES",))

ETHERGA = LocatedComponents(phys=PHY.TEMP_R, type="ELGA", location="RIGI", components=("TEMP",))

EPRESNO = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES",))


ECONTPC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELGA", location="RIGI", components=("SIXX", "SIYY", "SIZZ", "SIXZ")
)


ECONTNC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXZ")
)


EEFFRGC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELGA",
    location="RIGI",
    components=("NXX", "NYY", "NXY", "MXX", "MYY", "MXY"),
)


EEFFRNC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("NXX", "NYY", "NXY", "MXX", "MYY", "MXY")
)


ECONTPG = LocatedComponents(
    phys=PHY.SIEF_R, type="ELGA", location="RIGI", components=("SIXX", "SIYY", "SIZZ", "SIXZ")
)


ECONTNO = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXZ")
)


ECOEQPG = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=(
        "VMIS",
        "TRESCA",
        "PRIN_[3]",
        "VMIS_SG",
        "VECT_1_X",
        "VECT_1_Y",
        "VECT_1_Z",
        "VECT_2_X",
        "VECT_2_Y",
        "VECT_2_Z",
        "VECT_3_X",
        "VECT_3_Y",
        "VECT_3_Z",
        "TRSIG",
        "TRIAX",
    ),
)


EEFFRGR = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=("NXX", "NYY", "NXY", "MXX", "MYY", "MXY"),
)


EEFFRNR = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("NXX", "NYY", "NXY", "MXX", "MYY", "MXY")
)


EGAMIMA = LocatedComponents(
    phys=PHY.SPMX_R,
    type="ELGA",
    location="RIGI",
    components=("VAL", "NUCOU", "NUSECT", "NUFIBR", "POSIC", "POSIS"),
)


ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MECXSE3(Element):

    """Please document this element"""

    meshType = MT.SEG3
    elrefe = (ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("RIGI", "FPG1")),)
    calculs = (
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
        OP.CHAR_MECA_FFCO2D(
            te=224,
            para_in=((SP.PFFCO2D, CFORCEF), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FRCO2D(
            te=223,
            para_in=((SP.PFRCO2D, EFORCNO), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_HYDR_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_HYDR_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=233,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_F(
            te=397,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PPRESSF, CPRESSF), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=397,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PPRESSR, EPRESNO)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_ROTA_R(
            te=232,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PROTATR, LC.CROTATR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SECH_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_SECH_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=225,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.COOR_ELGA(
            te=478,
            para_in=((SP.PCACOQU, CCACOQU), (SP.PGEOMER, NGEOMER), (OP.COOR_ELGA.PNBSP_I, ENBSP_I)),
            para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R), (OP.COOR_ELGA.PCOORSU, EGGEOP_R)),
        ),
        OP.DEGE_ELGA(
            te=228,
            para_in=(
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.DEGE_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.DEGE_ELGA.PDEFOPG, EDEFGPG),),
        ),
        OP.DEGE_ELNO(
            te=228,
            para_in=(
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.DEGE_ELNO.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PDEFOGR, EDEFGNO),),
        ),
        OP.EFGE_ELGA(
            te=451,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PMATERC, LC.CMATERC),
                (OP.EFGE_ELGA.PNBSP_I, ENBSP_I),
                (SP.PSIEFR, ECONTPG),
            ),
            para_out=((SP.PEFGEC, EEFFRGC), (SP.PEFGER, EEFFRGR)),
        ),
        OP.EFGE_ELNO(
            te=185,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PNONLIN, LC.ENONLIN),
                (OP.EFGE_ELNO.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PEFFORC, EEFFRNC), (OP.EFGE_ELNO.PEFFORR, EEFFRNR)),
        ),
        OP.EPEQ_ELGA(
            te=335,
            para_in=((OP.EPEQ_ELGA.PDEFORR, EDEFOPG),),
            para_out=((OP.EPEQ_ELGA.PDEFOEQ, LC.EDFEQPG),),
        ),
        OP.EPEQ_ELNO(
            te=335,
            para_in=((OP.EPEQ_ELNO.PDEFORR, EDEFONO),),
            para_out=((OP.EPEQ_ELNO.PDEFOEQ, LC.EDFEQNO),),
        ),
        OP.EPSI_ELGA(
            te=237,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.EPSI_ELGA.PNBSP_I, ENBSP_I),
                (OP.EPSI_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PDEFOPC, EDEFOPC), (OP.EPSI_ELGA.PDEFOPG, EDEFOPG)),
        ),
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, EDEFOPG),),
            para_out=((SP.PDEFONC, EDEFONC), (SP.PDEFONO, EDEFONO)),
        ),
        OP.FORC_NODA(
            te=234,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, ECONTPG),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.FORC_NODA.PNBSP_I, ENBSP_I),
                (SP.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=239,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (OP.FULL_MECA.PNBSP_I, ENBSP_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.FULL_MECA.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, ECONTPG),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_INER(
            te=227,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_INER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMASSINE, LC.EMASSINE),),
        ),
        OP.MASS_MECA(
            te=226,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MINMAX_SP(te=99, para_out=((SP.PGAMIMA, EGAMIMA), (SP.PNOMIMA, LC.ENOMIMA))),
        OP.M_GAMMA(
            te=226,
            para_in=(
                (SP.PACCELR, DDL_MECA),
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.M_GAMMA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2), (OP.NSPG_NBVA.PNBSP_I, ENBSP_I)),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.PAS_COURANT(
            te=405,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PCOURAN, LC.ECOURAN),),
        ),
        OP.PREP_VRC(
            te=408,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.PREP_VRC.PINST_R, CTEMPSR),
                (OP.PREP_VRC.PNBSP_I, ENBSP_I),
                (SP.PTEMPEF, LC.CTEMPEF),
                (SP.PTEMPER, LC.NTEMPER),
                (SP.PGEOMER, NGEOMER),
            ),
            para_out=((SP.PTEMPCR, LC.CTEREFE),),
        ),
        OP.RAPH_MECA(
            te=239,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (OP.RAPH_MECA.PNBSP_I, ENBSP_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.RAPH_MECA.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, ECONTPG),
                (OP.RAPH_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.RIGI_MECA(
            te=221,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
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
        OP.RIGI_MECA_TANG(
            te=239,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA_TANG.PNBSP_I, ENBSP_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.RIGI_MECA_TANG.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, ECONTPG),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.SIEF_ELGA(
            te=237,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.SIEF_ELGA.PNBSP_I, ENBSP_I),
                (OP.SIEF_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PCONTRC, ECONTPC), (OP.SIEF_ELGA.PCONTRR, ECONTPG)),
        ),
        OP.SIEQ_ELGA(
            te=335,
            para_in=((OP.SIEQ_ELGA.PCONTRR, ECONTPG),),
            para_out=((OP.SIEQ_ELGA.PCONTEQ, ECOEQPG),),
        ),
        OP.SIEQ_ELNO(
            te=335,
            para_in=((OP.SIEQ_ELNO.PCONTRR, ECONTNO),),
            para_out=((OP.SIEQ_ELNO.PCONTEQ, LC.ECOEQNO),),
        ),
        OP.SIGM_ELGA(
            te=546,
            para_in=((SP.PSIEFR, ECONTPG),),
            para_out=((SP.PSIGMC, ECONTPC), (SP.PSIGMR, ECONTPG)),
        ),
        OP.SIGM_ELNO(
            te=4,
            para_in=((OP.SIGM_ELNO.PCONTRR, ECONTPG),),
            para_out=((SP.PSIEFNOC, ECONTNC), (OP.SIGM_ELNO.PSIEFNOR, ECONTNO)),
        ),
        OP.TEMP_ELGA(
            te=126,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.TEMP_ELGA.PNBSP_I, ENBSP_I),
                (OP.TEMP_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PTEMP_R, ETHERGA),),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PNBSP_I, ENBSP_I),)),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
                (OP.TOU_INI_ELGA.PSIEF_R, ECONTPG),
                (OP.TOU_INI_ELGA.PVARI_R, ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PSIEF_R, EEFFRNR),
            ),
        ),
    )
