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

# person_in_charge: sylvie.granet at edf.fr


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
    components=(("EN1", ("DX", "DY", "DZ", "PRE1", "TEMP")), ("EN2", ("DX", "DY", "DZ"))),
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

NFORCER = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ"))


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


EGGEOM_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z")
)


ENGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


E1NEUTK = LocatedComponents(phys=PHY.NEUT_K24, type="ELEM", components=("Z1",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


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
        "SIPXX",
        "SIPYY",
        "SIPZZ",
        "SIPXY",
        "SIPXZ",
        "SIPYZ",
        "M11",
        "FH11X",
        "FH11Y",
        "FH11Z",
        "ENT11",
        "QPRIM",
        "FHTX",
        "FHTY",
        "FHTZ",
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
        "SIPXX",
        "SIPYY",
        "SIPZZ",
        "SIPXY",
        "SIPXZ",
        "SIPYZ",
        "M11",
        "FH11X",
        "FH11Y",
        "FH11Z",
        "ENT11",
        "QPRIM",
        "FHTX",
        "FHTY",
        "FHTZ",
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
        "SIPXX",
        "SIPYY",
        "SIPZZ",
        "SIPXY",
        "SIPXZ",
        "SIPYZ",
        "M11",
        "FH11X",
        "FH11Y",
        "FH11Z",
        "ENT11",
        "QPRIM",
        "FHTX",
        "FHTY",
        "FHTZ",
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


ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class THM_HEXA20(Element):
    """Please document this element"""

    meshType = MT.HEXA20
    nodes = (
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
    )
    elrefe = (
        ElrefeLoc(MT.H20, gauss=("RIGI=FPG27", "MASS=FPG27", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.HE8, gauss=("RIGI=FPG27",)),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9",)),
    )
    calculs = (
        OP.ADD_SIGM(
            te=581,
            para_in=((SP.PEPCON1, ECONTPG), (SP.PEPCON2, ECONTPG)),
            para_out=((SP.PEPCON3, ECONTPG),),
        ),
        OP.AMOR_MECA(
            te=580,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.CHAR_MECA_FR3D3D(
            te=600,
            para_in=((SP.PFR3D3D, NFORCER), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=600,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
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
        OP.EPSI_ELGA(
            te=600,
            para_in=(
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
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
            te=600,
            para_in=(
                (SP.PSIEFR, ECONTPG),
                (SP.PCONTGM, ECONTPG),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=600,
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
                (OP.FULL_MECA.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, ECONTPG),
                (SP.PMATUNS, MMATUNS),
                (OP.FULL_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.FULL_MECA_ELAS(
            te=600,
            para_in=(
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
                (OP.FULL_MECA_ELAS.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA_ELAS.PCONTPR, ECONTPG),
                (SP.PMATUNS, MMATUNS),
                (OP.FULL_MECA_ELAS.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_INER(
            te=65,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_INER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMASSINE, LC.EMASSINE),),
        ),
        OP.MASS_MECA(
            te=12,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_MECA_DIAG(
            te=12,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA_DIAG.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_MECA_EXPLI(
            te=12,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA_EXPLI.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.M_GAMMA(
            te=12,
            para_in=(
                (SP.PACCELR, DDL_MECA),
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
            te=600,
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
            te=600,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC), (SP.PREFCO, LC.CRESTHM)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_MECA_ELAS(
            te=600,
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
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.RIGI_MECA_TANG(
            te=600,
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
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, ECONTPG),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.SIEF_ELNO(
            te=600,
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
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PGEOM_R, EGGEOM_R),
                (OP.TOU_INI_ELGA.PINST_R, LC.EGINST_R),
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
                (OP.TOU_INI_ELGA.PSIEF_R, ECONTPG),
                (OP.TOU_INI_ELGA.PVARI_R, ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, ENGEOM_R),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
            ),
        ),
        OP.VARI_ELNO(
            te=600,
            para_in=((OP.VARI_ELNO.PCOMPOR, LC.CCOMPOR), (SP.PVARIGR, ZVARIPG)),
            para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),),
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# ------------------------------------------------------------
class THM_PENTA15(THM_HEXA20):
    """Please document this element"""

    meshType = MT.PENTA15
    nodes = (
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
    )
    elrefe = (
        ElrefeLoc(MT.P15, gauss=("RIGI=FPG21", "MASS=FPG21", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.PE6, gauss=("RIGI=FPG21",)),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9",)),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6",)),
    )


# ------------------------------------------------------------
class THM_TETRA10(THM_HEXA20):
    """Please document this element"""

    meshType = MT.TETRA10
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8, 9, 10)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(MT.T10, gauss=("RIGI=FPG15", "MASS=FPG15", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.TE4, gauss=("RIGI=FPG15",)),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6",)),
    )


# ------------------------------------------------------------
class THM_HEXA20D(THM_HEXA20):
    """Please document this element"""

    meshType = MT.HEXA20
    nodes = (
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
    )
    elrefe = (
        ElrefeLoc(
            MT.H20, gauss=("RIGI=NOEU_S", "MASS=NOEU_S", "FPG1=FPG1"), mater=("RIGI", "FPG1")
        ),
        ElrefeLoc(MT.HE8, gauss=("RIGI=NOEU_S",)),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9",)),
    )


# ------------------------------------------------------------
class THM_PENTA15D(THM_HEXA20):
    """Please document this element"""

    meshType = MT.PENTA15
    nodes = (
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
    )
    elrefe = (
        ElrefeLoc(
            MT.P15, gauss=("RIGI=NOEU_S", "MASS=NOEU_S", "FPG1=FPG1"), mater=("RIGI", "FPG1")
        ),
        ElrefeLoc(MT.PE6, gauss=("RIGI=NOEU_S",)),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9",)),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6",)),
    )


# ------------------------------------------------------------
class THM_TETRA10D(THM_HEXA20):
    """Please document this element"""

    meshType = MT.TETRA10
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8, 9, 10)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.T10, gauss=("RIGI=NOEU_S", "MASS=NOEU_S", "FPG1=FPG1"), mater=("RIGI", "FPG1")
        ),
        ElrefeLoc(MT.TE4, gauss=("RIGI=NOEU_S",)),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6",)),
    )


# ------------------------------------------------------------
class THM_HEXA20S(THM_HEXA20):
    """Please document this element"""

    meshType = MT.HEXA20
    nodes = (
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
    )
    elrefe = (
        ElrefeLoc(
            MT.H20,
            gauss=("RIGI=FPG8NOS", "MASS=FPG8", "NOEU_S=NOEU_S", "FPG1=FPG1"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.HE8, gauss=("RIGI=FPG8NOS", "MASS=FPG8", "NOEU_S=NOEU_S")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG4", "MASS=FPG4", "NOEU_S=NOEU_S")),
    )


# ------------------------------------------------------------
class THM_PENTA15S(THM_HEXA20):
    """Please document this element"""

    meshType = MT.PENTA15
    nodes = (
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
    )
    elrefe = (
        ElrefeLoc(
            MT.P15,
            gauss=("RIGI=FPG6NOS", "MASS=FPG6", "NOEU_S=NOEU_S", "FPG1=FPG1"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.PE6, gauss=("RIGI=FPG6NOS", "MASS=FPG6", "NOEU_S=NOEU_S")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG4", "NOEU_S=NOEU_S")),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG3", "MASS=FPG3", "NOEU_S=NOEU_S")),
    )


# ------------------------------------------------------------
class THM_TETRA10S(THM_HEXA20):
    """Please document this element"""

    meshType = MT.TETRA10
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8, 9, 10)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.T10,
            gauss=("RIGI=FPG4NOS", "MASS=FPG4", "NOEU_S=NOEU_S", "FPG1=FPG1"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.TE4, gauss=("RIGI=FPG4NOS", "MASS=FPG4", "NOEU_S=NOEU_S")),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG3", "MASS=FPG3", "NOEU_S=NOEU_S")),
    )
