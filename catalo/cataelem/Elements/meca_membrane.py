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


CCACOQU = LocatedComponents(
    phys=PHY.CACOQU_R, type="ELEM", components=("EP", "ALPHA", "BETA", "TENS")
)


ECHGREP = LocatedComponents(phys=PHY.CHGREPER, type="ELEM", components=("NATCHG", "CMAT[9]"))


NDEPLAC = LocatedComponents(phys=PHY.DEPL_C, type="ELNO", components=("DX", "DY", "DZ"))


DDL_MECA = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))


EENERR = LocatedComponents(phys=PHY.ENER_R, type="ELEM", components=("TOTALE",))


CEPSINR = LocatedComponents(
    phys=PHY.EPSI_R, type="ELGA", location="RIGI", components=("EXX", "EYY", "EXY")
)


CEPSINF = LocatedComponents(phys=PHY.EPSI_F, type="ELEM", components=("EXX", "EYY", "EXY"))


EDEFONO = LocatedComponents(phys=PHY.EPSI_R, type="ELNO", components=("EXX", "EYY", "EXY"))


EDEFOPG = LocatedComponents(
    phys=PHY.EPSI_R, type="ELGA", location="RIGI", components=("EXX", "EYY", "EXY")
)


CPRESSF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PRES",))


EPRESNO = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES",))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


EREFCO = LocatedComponents(phys=PHY.PREC_R, type="ELEM", components=("EPSI",))


ECONTNC = LocatedComponents(phys=PHY.SIEF_C, type="ELNO", components=("NXX", "NYY", "NXY"))


ECONTPG = LocatedComponents(
    phys=PHY.SIEF_R, type="ELGA", location="RIGI", components=("NXX", "NYY", "NXY")
)


ECONTNO = LocatedComponents(phys=PHY.SIEF_R, type="ELNO", components=("NXX", "NYY", "NXY"))


EGAMIMA = LocatedComponents(
    phys=PHY.SPMX_R,
    type="ELGA",
    location="RIGI",
    components=("VAL", "NUCOU", "NUSECT", "NUFIBR", "POSIC", "POSIS"),
)

EPRESNO = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES",))

ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)

MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MEMBTR3(Element):
    """Please document this element"""

    meshType = MT.TRIA3
    elrefe = (
        ElrefeLoc(
            MT.TR3, gauss=("RIGI=FPG1", "MASS=FPG3", "FPG1=FPG1"), mater=("RIGI", "MASS", "FPG1")
        ),
    )
    calculs = (
        OP.ADD_SIGM(
            te=581,
            para_in=((SP.PEPCON1, ECONTPG), (SP.PEPCON2, ECONTPG)),
            para_out=((SP.PEPCON3, ECONTPG),),
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
        OP.CHAR_MECA_EPSA_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_EPSA_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_EPSI_R(
            te=434,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.CHAR_MECA_EPSI_R.PCOMPOR, LC.CCOMPOR),
                (SP.PEPSINR, CEPSINR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_EPSI_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_EPSI_F(
            te=434,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PEPSINF, CEPSINF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_EPSI_F.PVARCPR, LC.ZVARCPG),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_HYDR_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_HYDR_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=434,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.CHAR_MECA_PESA_R.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=18,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC), (SP.PPRESSR, EPRESNO)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRSU_R(
            te=424,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPRESSR, EPRESNO),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SECH_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_SECH_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=434,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.CHAR_MECA_TEMP_R.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.EFGE_ELGA(te=436, para_in=((SP.PSIEFR, ECONTPG),), para_out=((SP.PEFGER, ECONTPG),)),
        OP.EFGE_ELNO(
            te=185,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.EFGE_ELNO.PCONTRR, ECONTPG),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PNONLIN, LC.ENONLIN),
                (SP.PINSTR, CTEMPSR),
                (OP.EFGE_ELNO.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PEFFORC, ECONTNC), (OP.EFGE_ELNO.PEFFORR, ECONTNO)),
        ),
        OP.EPOT_ELEM(
            te=436,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.EPOT_ELEM.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPOT_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPOT_ELEM.PENERDR, EENERR),),
        ),
        OP.EPSI_ELGA(
            te=436,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.EPSI_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.EPSI_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.EPSI_ELGA.PDEFOPG, EDEFOPG),),
        ),
        OP.EPSI_ELNO(
            te=4, para_in=((OP.EPSI_ELNO.PDEFOPG, EDEFOPG),), para_out=((SP.PDEFONO, EDEFONO),)
        ),
        OP.FORC_NODA(
            te=434,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PSIEFR, ECONTPG),
                (SP.PGEOMER, NGEOMER),
                (SP.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=435,
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
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.FULL_MECA.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, ECONTPG),
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.FULL_MECA_ELAS(
            te=435,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA_ELAS.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.FULL_MECA_ELAS.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA_ELAS.PCONTPR, ECONTPG),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA_ELAS.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_INER(
            te=436,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_INER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMASSINE, LC.EMASSINE),),
        ),
        OP.MASS_MECA(
            te=439,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.MASS_MECA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_MECA_EXPLI(
            te=439,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.MASS_MECA_EXPLI.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA_EXPLI.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
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
            te=435,
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
        OP.RAPH_MECA_IMPLEX(
            te=435,
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
        OP.REFE_FORC_NODA(
            te=434,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.REFE_FORC_NODA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PREFCO, EREFCO),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.REPERE_LOCAL(
            te=134,
            para_in=((SP.PCACOQU, CCACOQU), (SP.PGEOMER, NGEOMER)),
            para_out=(
                (OP.REPERE_LOCAL.PMATPASS, ECHGREP),
                (SP.PREPLO1, LC.CGEOM3D),
                (SP.PREPLO2, LC.CGEOM3D),
                (SP.PREPLO3, LC.CGEOM3D),
            ),
        ),
        OP.RIGI_MECA(
            te=435,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_ELAS(
            te=435,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_ELAS.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.RIGI_MECA_ELAS.PVARIMR, ZVARIPG),
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
        OP.RIGI_MECA_IMPLEX(
            te=435,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_IMPLEX.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_IMPLEX.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_IMPLEX.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.RIGI_MECA_IMPLEX.PVARIMR, ZVARIPG),
            ),
            para_out=((SP.PCONTXR, ECONTPG), (SP.PMATUUR, MMATUUR)),
        ),
        OP.RIGI_MECA_PRSU_R(
            te=435,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PPRESSR, EPRESNO),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.RIGI_MECA_PRSU_F(
            te=435,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PPRESSF, CPRESSF),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.CHAR_MECA_PRSU_F(
            te=435,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PPRESSF, CPRESSF),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_MECA_TANG(
            te=435,
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
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.RIGI_MECA_TANG.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, ECONTPG),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.SIEF_ELGA(
            te=436,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.SIEF_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.SIEF_ELGA.PCONTRR, ECONTPG),),
        ),
        OP.SIEF_ELNO(
            te=4,
            para_in=((OP.SIEF_ELNO.PCONTRR, ECONTPG), (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PSIEFNOC, ECONTNC), (OP.SIEF_ELNO.PSIEFNOR, ECONTNO)),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
                (OP.TOU_INI_ELGA.PSIEF_R, ECONTPG),
                (OP.TOU_INI_ELGA.PVARI_R, ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
        OP.VARI_ELNO(
            te=4, para_in=((SP.PVARIGR, ZVARIPG),), para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),)
        ),
    )


# ------------------------------------------------------------
class MEMBTR6(MEMBTR3):
    """Please document this element"""

    meshType = MT.TRIA6
    elrefe = (
        ElrefeLoc(
            MT.TR6, gauss=("RIGI=FPG3", "MASS=FPG6", "FPG1=FPG1"), mater=("RIGI", "MASS", "FPG1")
        ),
    )


# ------------------------------------------------------------
class MEMBTR7(MEMBTR3):
    """Please document this element"""

    meshType = MT.TRIA7
    elrefe = (ElrefeLoc(MT.TR7, gauss=("RIGI=FPG3", "MASS=FPG7"), mater=("RIGI", "MASS")),)


# ------------------------------------------------------------
class MEMBQU4(MEMBTR3):
    """Please document this element"""

    meshType = MT.QUAD4
    elrefe = (
        ElrefeLoc(
            MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG1=FPG1"), mater=("RIGI", "MASS", "FPG1")
        ),
    )


# ------------------------------------------------------------
class MEMBQU8(MEMBTR3):
    """Please document this element"""

    meshType = MT.QUAD8
    elrefe = (
        ElrefeLoc(
            MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9", "FPG1=FPG1"), mater=("RIGI", "MASS", "FPG1")
        ),
    )


# ------------------------------------------------------------
class MEMBQU9(MEMBTR3):
    """Please document this element"""

    meshType = MT.QUAD9
    elrefe = (ElrefeLoc(MT.QU9, gauss=("RIGI=FPG9", "MASS=FPG9"), mater=("RIGI", "MASS")),)
