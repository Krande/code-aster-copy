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

NDEPLAR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("DX", "DY", "DZ", "VARI", "LAG_GV")), ("EN2", ("DX", "DY", "DZ"))),
)


CEPSINF = LocatedComponents(
    phys=PHY.EPSI_F, type="ELEM", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)


CEPSINR = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"),
)


EDEFOPC = LocatedComponents(
    phys=PHY.EPSI_C,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"),
)


EDEFONC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)


EDEFOPG = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"),
)


EDEFONO = LocatedComponents(
    phys=PHY.EPSI_R, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)


CFORCEF = LocatedComponents(phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ"))


NFORCER = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ"))


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOM_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z")
)


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


ECOPILO = LocatedComponents(
    phys=PHY.PILO_R, type="ELGA", location="RIGI", components=("A0", "A[3]", "ETA")
)


EREFCO = LocatedComponents(phys=PHY.PREC_R, type="ELEM", components=("SIGM", "VARI", "LAG_GV"))


ESIGMPC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
)


ESIGMNC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
)


ECONTNC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELNO",
    components=(
        "SIXX",
        "SIYY",
        "SIZZ",
        "SIXY",
        "SIXZ",
        "SIYZ",
        "SIGV_A",
        "SIGV_L",
        "SIGV_GX",
        "SIGV_GY",
        "SIGV_GZ",
    ),
)


ESIGMPG = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
)


ESIGMNO = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
)


ECONTPG = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=(
        "SIXX",
        "SIYY",
        "SIZZ",
        "SIXY",
        "SIXZ",
        "SIYZ",
        "SIGV_A",
        "SIGV_L",
        "SIGV_GX",
        "SIGV_GY",
        "SIGV_GZ",
    ),
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


ECONTNO = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELNO",
    components=(
        "SIXX",
        "SIYY",
        "SIZZ",
        "SIXY",
        "SIXZ",
        "SIYZ",
        "SIGV_A",
        "SIGV_L",
        "SIGV_GX",
        "SIGV_GY",
        "SIGV_GZ",
    ),
)

ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MVECTDR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=NDEPLAR)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MVCA_HEXA20(Element):
    """Please document this element"""

    meshType = MT.HEXA20
    nodes = (
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
    )
    elrefe = (
        ElrefeLoc(
            MT.H20,
            gauss=("RIGI=FPG8", "MASS=FPG27", "FPG1=FPG1", "NOEU=NOEU"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.HE8, gauss=("RIGI=FPG8", "MASS=FPG27")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9", "NOEU=NOEU")),
    )
    calculs = (
        OP.AMOR_MECA(
            te=121,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMASSEL, MMATUUR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PRIGINS, MMATUNS),
                (OP.AMOR_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.CHAR_MECA_EPSA_R(
            te=426,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_EPSA_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_EPSI_F(
            te=49,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PEPSINF, CEPSINF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_EPSI_F.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_EPSI_R(
            te=49,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PEPSINR, CEPSINR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_EPSI_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_FF3D3D(
            te=17,
            para_in=((SP.PFF3D3D, CFORCEF), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_FR3D3D(
            te=16,
            para_in=((SP.PFR3D3D, NFORCER), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_HYDR_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_HYDR_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_META_Z(
            te=358,
            para_in=(
                (OP.CHAR_MECA_META_Z.PCOMPOR, LC.CCOMPOR),
                (OP.CHAR_MECA_META_Z.PCONTMR, ECONTPG),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.CHAR_MECA_META_Z.PVARCPR, LC.ZVARCPG),
                (OP.CHAR_MECA_META_Z.PVARIPR, ZVARIPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=15,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_ROTA_R(
            te=14,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PROTATR, LC.CROTATR),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_SECH_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_SECH_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
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
        OP.EPGQ_ELGA(
            te=335,
            para_in=((OP.EPGQ_ELGA.PDEFORR, EDEFOPG),),
            para_out=((OP.EPGQ_ELGA.PDEFOEQ, LC.EDFEQPG),),
        ),
        OP.EPGQ_ELNO(
            te=335,
            para_in=((OP.EPGQ_ELNO.PDEFORR, EDEFONO),),
            para_out=((OP.EPGQ_ELNO.PDEFOEQ, LC.EDFEQNO),),
        ),
        OP.EPSG_ELGA(
            te=25,
            para_in=(
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPSG_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPSG_ELGA.PDEFOPG, EDEFOPG),),
        ),
        OP.EPSG_ELNO(
            te=4, para_in=((OP.EPSG_ELNO.PDEFOPG, EDEFOPG),), para_out=((SP.PDEFONO, EDEFONO),)
        ),
        OP.EPSI_ELGA(
            te=25,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPSI_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PDEFOPC, EDEFOPC), (OP.EPSI_ELGA.PDEFOPG, EDEFOPG)),
        ),
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, EDEFOPG),),
            para_out=((SP.PDEFONC, EDEFONC), (SP.PDEFONO, EDEFONO)),
        ),
        OP.FORC_NODA(
            te=508,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, ECONTPG),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=545,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
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
            te=545,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
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
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA_ELAS.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_MECA(
            te=562,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_MECA_DIAG(
            te=562,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA_DIAG.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_MECA_EXPLI(
            te=562,
            para_in=(
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
        OP.PILO_PRED_DEFO(
            te=518,
            para_in=(
                (SP.PCDTAU, LC.CCDTAU),
                (OP.PILO_PRED_DEFO.PCOMPOR, LC.CCOMPOR),
                (SP.PDDEPLR, DDL_MECA),
                (SP.PDEPL0R, DDL_MECA),
                (SP.PDEPL1R, DDL_MECA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PTYPEPI, LC.CTYPEPI),
            ),
            para_out=((OP.PILO_PRED_DEFO.PCOPILO, ECOPILO),),
        ),
        OP.PILO_PRED_ELAS(
            te=518,
            para_in=(
                (SP.PBORNPI, LC.CBORNPI),
                (SP.PCDTAU, LC.CCDTAU),
                (OP.PILO_PRED_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.PILO_PRED_ELAS.PCONTMR, ECONTPG),
                (SP.PDDEPLR, DDL_MECA),
                (SP.PDEPL0R, DDL_MECA),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PDEPL1R, DDL_MECA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTYPEPI, LC.CTYPEPI),
                (OP.PILO_PRED_ELAS.PVARIMR, ZVARIPG),
            ),
            para_out=((OP.PILO_PRED_ELAS.PCOPILO, ECOPILO),),
        ),
        OP.RAPH_MECA(
            te=545,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
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
            te=508,
            para_in=(
                (OP.REFE_FORC_NODA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PREFCO, EREFCO),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.REPERE_LOCAL(
            te=133,
            para_in=((SP.PCAMASS, LC.CCAMA3D), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PREPLO1, LC.CGEOM3D), (SP.PREPLO2, LC.CGEOM3D), (SP.PREPLO3, LC.CGEOM3D)),
        ),
        OP.RIGI_MECA_ELAS(
            te=545,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
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
                (OP.RIGI_MECA_ELAS.PVARIMR, ZVARIPG),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.RIGI_MECA_TANG(
            te=545,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
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
        OP.SIEF_ELNO(
            te=4,
            para_in=((OP.SIEF_ELNO.PCONTRR, ECONTPG), (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PSIEFNOC, ECONTNC), (OP.SIEF_ELNO.PSIEFNOR, ECONTNO)),
        ),
        OP.SIEQ_ELGA(
            te=335,
            para_in=((OP.SIEQ_ELGA.PCONTRR, ESIGMPG),),
            para_out=((OP.SIEQ_ELGA.PCONTEQ, ECOEQPG),),
        ),
        OP.SIEQ_ELNO(
            te=335,
            para_in=((OP.SIEQ_ELNO.PCONTRR, ESIGMNO),),
            para_out=((OP.SIEQ_ELNO.PCONTEQ, LC.ECOEQNO),),
        ),
        OP.SIGM_ELGA(
            te=546,
            para_in=((SP.PSIEFR, ESIGMPG),),
            para_out=((SP.PSIGMC, ESIGMPC), (SP.PSIGMR, ESIGMPG)),
        ),
        OP.SIGM_ELNO(
            te=4,
            para_in=((OP.SIGM_ELNO.PCONTRR, ESIGMPG),),
            para_out=((SP.PSIEFNOC, ESIGMNC), (OP.SIGM_ELNO.PSIEFNOR, ESIGMNO)),
        ),
        OP.SIMY_ELGA(
            te=6,
            para_in=((OP.SIMY_ELGA.PCONTRR, LC.EGIG3DR), (SP.PGEOMER, LC.EGEOM3D)),
            para_out=((OP.SIMY_ELGA.PSIEFNOR, LC.EGIG3DR),),
        ),
        OP.SING_ELEM(te=99, para_out=((SP.PSING_R, LC.ESINGUL),)),
        OP.SING_ELNO(te=99, para_out=((SP.PSINGNO, LC.ESINGNO),)),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PDEPL_R, LC.EGDEP3D),
                (OP.TOU_INI_ELGA.PDOMMAG, LC.EDOMGGA),
                (OP.TOU_INI_ELGA.PGEOM_R, EGGEOM_R),
                (OP.TOU_INI_ELGA.PINST_R, LC.EGINST_R),
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
                (OP.TOU_INI_ELGA.PSIEF_R, ECONTPG),
                (SP.PVALO_R, LC.EGTINIV),
                (OP.TOU_INI_ELGA.PVARI_R, ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PDOMMAG, LC.EDOMGNO),
                (OP.TOU_INI_ELNO.PEPSI_R, EDEFONO),
                (OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PSIEF_R, ECONTNO),
                (OP.TOU_INI_ELNO.PVARI_R, LC.ZVARINO),
            ),
        ),
        OP.VARI_ELNO(
            te=4, para_in=((SP.PVARIGR, ZVARIPG),), para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),)
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# ------------------------------------------------------------
class MVCA_TETRA10(MVCA_HEXA20):
    """Please document this element"""

    meshType = MT.TETRA10
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8, 9, 10)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.T10,
            gauss=("RIGI=FPG4", "MASS=FPG15", "FPG1=FPG1", "NOEU=NOEU"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.TE4, gauss=("RIGI=FPG4", "MASS=FPG15")),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6", "MASS=FPG6", "NOEU=NOEU")),
    )


# ------------------------------------------------------------
class MVCA_PENTA15(MVCA_HEXA20):
    """Please document this element"""

    meshType = MT.PENTA15
    nodes = (
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
    )
    elrefe = (
        ElrefeLoc(
            MT.P15,
            gauss=("RIGI=FPG6", "MASS=FPG21", "FPG1=FPG1", "NOEU=NOEU"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.PE6, gauss=("RIGI=FPG6", "MASS=FPG6")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9", "NOEU=NOEU")),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6", "MASS=FPG6", "NOEU=NOEU")),
    )


# ------------------------------------------------------------
class MVCA_PYRAM13(MVCA_HEXA20):
    """Please document this element"""

    meshType = MT.PYRAM13
    nodes = (SetOfNodes("EN2", (6, 7, 8, 9, 10, 11, 12, 13)), SetOfNodes("EN1", (1, 2, 3, 4, 5)))
    elrefe = (
        ElrefeLoc(
            MT.P13,
            gauss=("RIGI=FPG5", "MASS=FPG10", "FPG1=FPG1", "NOEU=NOEU"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.PY5, gauss=("RIGI=FPG5", "MASS=FPG5")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9", "NOEU=NOEU")),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6", "MASS=FPG6", "NOEU=NOEU")),
    )
