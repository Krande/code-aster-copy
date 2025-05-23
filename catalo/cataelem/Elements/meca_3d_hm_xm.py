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

# person_in_charge: daniele.colombo at ifpen.fr
# CATALOGUE DES ELEMENTS 3D HM-X-FEM MULTI HEAVISIDE SANS CONTACT


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
    components=(
        ("EN1", ("DX", "DY", "DZ", "PRE1", "H1X", "H1Y", "H1Z", "H1PRE1")),
        ("EN2", ("DX", "DY", "DZ", "H1X", "H1Y", "H1Z")),
        (
            "EN3",
            (
                "DX",
                "DY",
                "DZ",
                "PRE1",
                "H1X",
                "H1Y",
                "H1Z",
                "H1PRE1",
                "H2X",
                "H2Y",
                "H2Z",
                "H2PRE1",
            ),
        ),
        ("EN4", ("DX", "DY", "DZ", "H1X", "H1Y", "H1Z", "H2X", "H2Y", "H2Z")),
        (
            "EN5",
            (
                "DX",
                "DY",
                "DZ",
                "PRE1",
                "H1X",
                "H1Y",
                "H1Z",
                "H1PRE1",
                "H2X",
                "H2Y",
                "H2Z",
                "H2PRE1",
                "H3X",
                "H3Y",
                "H3Z",
                "H3PRE1",
            ),
        ),
        ("EN6", ("DX", "DY", "DZ", "H1X", "H1Y", "H1Z", "H2X", "H2Y", "H2Z", "H3X", "H3Y", "H3Z")),
    ),
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


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOM_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="XFEM", components=("X", "Y", "Z")
)


ENGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


STANO_I = LocatedComponents(phys=PHY.N120_I, type="ELNO", components=("X1",))


E1280NEI = LocatedComponents(phys=PHY.N1280I, type="ELEM", components=("X[1280]",))


E132NEUT = LocatedComponents(phys=PHY.N132_R, type="ELEM", components=("X[132]",))


E340NEUT = LocatedComponents(phys=PHY.N1360R, type="ELEM", components=("X[340]",))


E60NEUTI = LocatedComponents(phys=PHY.N240_I, type="ELEM", components=("X[60]",))


E612NEUT = LocatedComponents(phys=PHY.N2448R, type="ELEM", components=("X[612]",))


E220NEUT = LocatedComponents(phys=PHY.N480_R, type="ELEM", components=("X[220]",))


E180NEUI = LocatedComponents(phys=PHY.N720_I, type="ELEM", components=("X[180]",))


E792NEUT = LocatedComponents(phys=PHY.N792_R, type="ELEM", components=("X[792]",))


E204NEUT = LocatedComponents(phys=PHY.N816_R, type="ELEM", components=("X[204]",))


E60NEUI = LocatedComponents(phys=PHY.N960_I, type="ELEM", components=("X[60]",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="XFEM", components=("X[30]",))


E1NEUTK = LocatedComponents(phys=PHY.NEUT_K24, type="ELEM", components=("Z1",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="XFEM", components=("X[30]",))


CPRESSF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PRES",))


EPRESNO = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES",))


ESIGMPC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
)


ESIGMNC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
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
    location="XFEM",
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


ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="XFEM", components=("VARI",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class HM_HEXA20_XH1(Element):

    """Please document this element"""

    meshType = MT.HEXA20
    nodes = (
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
    )
    elrefe = (
        ElrefeLoc(
            MT.H20,
            gauss=("RIGI=FPG27", "MASS=FPG27", "XFEM=XFEM1440", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.HE8, gauss=("RIGI=FPG27", "MASS=FPG27")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(MT.T10, gauss=("XINT=FPG15", "NOEU=NOEU")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("FPG4=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )
    calculs = (
        OP.CHAR_MECA_PESA_R(
            te=588,
            para_in=(
                (OP.CHAR_MECA_PESA_R.PBASLOR, LC.N9NEUT_R),
                (OP.CHAR_MECA_PESA_R.PCNSETO, E1280NEI),
                (OP.CHAR_MECA_PESA_R.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_PESA_R.PHEAVTO, LC.E128NEUI),
                (OP.CHAR_MECA_PESA_R.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_PESA_R.PLONCHA, LC.E10NEUTI),
                (OP.CHAR_MECA_PESA_R.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_PESA_R.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PPINTTO, E132NEUT),
                (OP.CHAR_MECA_PESA_R.PPMILTO, E792NEUT),
                (OP.CHAR_MECA_PESA_R.PSTANO, STANO_I),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_F(
            te=37,
            para_in=(
                (OP.CHAR_MECA_PRES_F.PAINTER, E340NEUT),
                (OP.CHAR_MECA_PRES_F.PBASECO, E612NEUT),
                (OP.CHAR_MECA_PRES_F.PCFACE, E180NEUI),
                (OP.CHAR_MECA_PRES_F.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_PRES_F.PHEA_FA, E60NEUTI),
                (OP.CHAR_MECA_PRES_F.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_PRES_F.PLONGCO, LC.E3NEUTI),
                (OP.CHAR_MECA_PRES_F.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_PRES_F.PLST, LC.N1NEUT_R),
                (OP.CHAR_MECA_PRES_F.PPINTER, E204NEUT),
                (SP.PPRESSF, CPRESSF),
                (OP.CHAR_MECA_PRES_F.PSTANO, STANO_I),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=37,
            para_in=(
                (OP.CHAR_MECA_PRES_R.PAINTER, E340NEUT),
                (OP.CHAR_MECA_PRES_R.PBASECO, E612NEUT),
                (OP.CHAR_MECA_PRES_R.PCFACE, E180NEUI),
                (OP.CHAR_MECA_PRES_R.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_PRES_R.PHEA_FA, E60NEUTI),
                (OP.CHAR_MECA_PRES_R.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_PRES_R.PLONGCO, LC.E3NEUTI),
                (OP.CHAR_MECA_PRES_R.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_PRES_R.PLST, LC.N1NEUT_R),
                (OP.CHAR_MECA_PRES_R.PPINTER, E204NEUT),
                (SP.PPRESSR, EPRESNO),
                (OP.CHAR_MECA_PRES_R.PSTANO, STANO_I),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
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
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, EDEFOPG),),
            para_out=((SP.PDEFONC, EDEFONC), (SP.PDEFONO, EDEFONO)),
        ),
        OP.FORC_NODA(
            te=588,
            para_in=(
                (OP.FORC_NODA.PBASLOR, LC.N9NEUT_R),
                (OP.FORC_NODA.PCNSETO, E1280NEI),
                (SP.PSIEFR, ECONTPG),
                (OP.FORC_NODA.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.FORC_NODA.PHEAVTO, LC.E128NEUI),
                (OP.FORC_NODA.PHEA_NO, LC.N5NEUTI),
                (OP.FORC_NODA.PLONCHA, LC.E10NEUTI),
                (OP.FORC_NODA.PLSN, LC.N1NEUT_R),
                (OP.FORC_NODA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.FORC_NODA.PPINTTO, E132NEUT),
                (OP.FORC_NODA.PPMILTO, E792NEUT),
                (OP.FORC_NODA.PSTANO, STANO_I),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=588,
            para_in=(
                (OP.FULL_MECA.PBASLOR, LC.N9NEUT_R),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCNSETO, E1280NEI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (OP.FULL_MECA.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.FULL_MECA.PHEAVTO, LC.E128NEUI),
                (OP.FULL_MECA.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.FULL_MECA.PLONCHA, LC.E10NEUTI),
                (OP.FULL_MECA.PLSN, LC.N1NEUT_R),
                (OP.FULL_MECA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.FULL_MECA.PPINTTO, E132NEUT),
                (OP.FULL_MECA.PPMILTO, E792NEUT),
                (OP.FULL_MECA.PSTANO, STANO_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PHEAVNO, LC.FISNO_I),
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
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.INI_XFEM_ELNO(
            te=99,
            para_out=(
                (OP.INI_XFEM_ELNO.PBASLOR, LC.N9NEUT_R),
                (OP.INI_XFEM_ELNO.PFISNO, LC.FISNO_I),
                (OP.INI_XFEM_ELNO.PLSN, LC.N1NEUT_R),
                (OP.INI_XFEM_ELNO.PLST, LC.N1NEUT_R),
                (OP.INI_XFEM_ELNO.PSTANO, STANO_I),
            ),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.RAPH_MECA(
            te=588,
            para_in=(
                (OP.RAPH_MECA.PBASLOR, LC.N9NEUT_R),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RAPH_MECA.PCNSETO, E1280NEI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.RAPH_MECA.PHEAVTO, LC.E128NEUI),
                (OP.RAPH_MECA.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.RAPH_MECA.PLONCHA, LC.E10NEUTI),
                (OP.RAPH_MECA.PLSN, LC.N1NEUT_R),
                (OP.RAPH_MECA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RAPH_MECA.PPINTTO, E132NEUT),
                (OP.RAPH_MECA.PPMILTO, E792NEUT),
                (OP.RAPH_MECA.PSTANO, STANO_I),
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
        OP.RIGI_MECA_TANG(
            te=588,
            para_in=(
                (OP.RIGI_MECA_TANG.PBASLOR, LC.N9NEUT_R),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCNSETO, E1280NEI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.RIGI_MECA_TANG.PHEAVTO, LC.E128NEUI),
                (OP.RIGI_MECA_TANG.PFISNO, LC.FISNO_I),
                (OP.RIGI_MECA_TANG.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.RIGI_MECA_TANG.PLONCHA, LC.E10NEUTI),
                (OP.RIGI_MECA_TANG.PLSN, LC.N1NEUT_R),
                (OP.RIGI_MECA_TANG.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA_TANG.PPINTTO, E132NEUT),
                (OP.RIGI_MECA_TANG.PPMILTO, E792NEUT),
                (OP.RIGI_MECA_TANG.PSTANO, STANO_I),
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
        OP.TOPOFA(
            te=510,
            para_in=(
                (OP.TOPOFA.PAINTTO, E220NEUT),
                (OP.TOPOFA.PCNSETO, E1280NEI),
                (SP.PDECOU, LC.E1NEUK8),
                (SP.PFISCO, LC.FISCO_I),
                (SP.PGEOMER, NGEOMER),
                (SP.PGRADLN, LC.N3NEUT_R),
                (SP.PGRADLT, LC.N3NEUT_R),
                (OP.TOPOFA.PHEAVTO, LC.E128NEUI),
                (OP.TOPOFA.PLONCHA, LC.E10NEUTI),
                (OP.TOPOFA.PLSN, LC.N1NEUT_R),
                (OP.TOPOFA.PLST, LC.N1NEUT_R),
                (OP.TOPOFA.PPINTTO, E132NEUT),
                (OP.TOPOFA.PPMILTO, E792NEUT),
                (OP.TOPOFA.PSTANO, STANO_I),
                (SP.PTYPDIS, LC.E1NEUTI),
            ),
            para_out=(
                (OP.TOPOFA.PAINTER, E340NEUT),
                (OP.TOPOFA.PBASECO, E612NEUT),
                (OP.TOPOFA.PCFACE, E180NEUI),
                (SP.PGESCLA, E204NEUT),
                (OP.TOPOFA.PHEAVFA, E60NEUI),
                (OP.TOPOFA.PLONGCO, LC.E3NEUTI),
                (OP.TOPOFA.PPINTER, E204NEUT),
            ),
        ),
        OP.TOPONO(
            te=120,
            para_in=(
                (OP.TOPONO.PCNSETO, E1280NEI),
                (SP.PFISCO, LC.FISCO_I),
                (OP.TOPONO.PFISNO, LC.FISNO_I),
                (OP.TOPONO.PHEAVFA, E60NEUI),
                (OP.TOPONO.PHEAVTO, LC.E128NEUI),
                (SP.PLEVSET, LC.N1NEUT_R),
                (OP.TOPONO.PLONCHA, LC.E10NEUTI),
                (OP.TOPONO.PLONGCO, LC.E3NEUTI),
            ),
            para_out=(
                (OP.TOPONO.PHEA_FA, E60NEUTI),
                (OP.TOPONO.PHEA_NO, LC.N5NEUTI),
                (OP.TOPONO.PHEA_SE, LC.E128NEUI),
            ),
        ),
        OP.TOPOSE(
            te=514,
            para_in=((SP.PFISCO, LC.FISCO_I), (SP.PGEOMER, NGEOMER), (SP.PLEVSET, LC.N1NEUT_R)),
            para_out=(
                (OP.TOPOSE.PAINTTO, E220NEUT),
                (OP.TOPOSE.PCNSETO, E1280NEI),
                (OP.TOPOSE.PHEAVTO, LC.E128NEUI),
                (OP.TOPOSE.PLONCHA, LC.E10NEUTI),
                (OP.TOPOSE.PPINTTO, E132NEUT),
                (OP.TOPOSE.PPMILTO, E792NEUT),
            ),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
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
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, ENGEOM_R),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
            ),
        ),
    )


# ------------------------------------------------------------
class HM_PENTA15_XH1(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.PENTA15
    nodes = (
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
    )
    elrefe = (
        ElrefeLoc(
            MT.P15, gauss=("RIGI=FPG21", "MASS=FPG21", "XFEM=XFEM720", "FPG1=FPG1"), mater=("XFEM",)
        ),
        ElrefeLoc(MT.T10, gauss=("XINT=FPG15", "NOEU=NOEU")),
        ElrefeLoc(MT.PE6, gauss=("RIGI=FPG21", "MASS=FPG21")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("MASS=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )


# ------------------------------------------------------------
class HM_PYRAM13_XH1(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.PYRAM13
    nodes = (SetOfNodes("EN2", (6, 7, 8, 9, 10, 11, 12, 13)), SetOfNodes("EN1", (1, 2, 3, 4, 5)))
    elrefe = (
        ElrefeLoc(
            MT.P13, gauss=("RIGI=FPG10", "MASS=FPG10", "XFEM=XFEM540", "FPG1=FPG1"), mater=("XFEM",)
        ),
        ElrefeLoc(MT.T10, gauss=("XINT=FPG15", "NOEU=NOEU")),
        ElrefeLoc(MT.PY5, gauss=("RIGI=FPG10", "MASS=FPG10")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("MASS=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )


# ------------------------------------------------------------
class HM_TETRA10_XH1(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.TETRA10
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8, 9, 10)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.T10,
            gauss=("RIGI=FPG5", "MASS=FPG5", "XINT=FPG15", "XFEM=XFEM270", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TE4, gauss=("XINT=FPG15", "MASS=FPG15")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("MASS=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )


# ------------------------------------------------------------
class HM_HEXA20_XH2(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.HEXA20
    nodes = (
        SetOfNodes("EN4", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN3", (1, 2, 3, 4, 5, 6, 7, 8)),
    )
    elrefe = (
        ElrefeLoc(
            MT.H20,
            gauss=("RIGI=FPG27", "MASS=FPG27", "XFEM=XFEM2880", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.HE8, gauss=("RIGI=FPG27", "MASS=FPG27")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(MT.T10, gauss=("XINT=FPG15", "NOEU=NOEU")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("FPG4=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )


# ------------------------------------------------------------
class HM_PENTA15_XH2(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.PENTA15
    nodes = (
        SetOfNodes("EN4", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN3", (1, 2, 3, 4, 5, 6)),
    )
    elrefe = (
        ElrefeLoc(
            MT.P15,
            gauss=("RIGI=FPG21", "MASS=FPG21", "XFEM=XFEM1440", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.T10, gauss=("XINT=FPG15", "NOEU=NOEU")),
        ElrefeLoc(MT.PE6, gauss=("RIGI=FPG21", "MASS=FPG21")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("MASS=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )


# ------------------------------------------------------------
class HM_PYRAM13_XH2(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.PYRAM13
    nodes = (SetOfNodes("EN4", (6, 7, 8, 9, 10, 11, 12, 13)), SetOfNodes("EN3", (1, 2, 3, 4, 5)))
    elrefe = (
        ElrefeLoc(
            MT.P13,
            gauss=("RIGI=FPG10", "MASS=FPG10", "XFEM=XFEM1080", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.T10, gauss=("XINT=FPG15", "NOEU=NOEU")),
        ElrefeLoc(MT.PY5, gauss=("RIGI=FPG10", "MASS=FPG10")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("MASS=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )


# ------------------------------------------------------------
class HM_TETRA10_XH2(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.TETRA10
    nodes = (SetOfNodes("EN4", (5, 6, 7, 8, 9, 10)), SetOfNodes("EN3", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.T10,
            gauss=("RIGI=FPG5", "MASS=FPG5", "XINT=FPG15", "XFEM=XFEM540", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TE4, gauss=("XINT=FPG15", "MASS=FPG15")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("MASS=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )


# ------------------------------------------------------------
class HM_HEXA20_XH3(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.HEXA20
    nodes = (
        SetOfNodes("EN6", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN5", (1, 2, 3, 4, 5, 6, 7, 8)),
    )
    elrefe = (
        ElrefeLoc(
            MT.H20,
            gauss=("RIGI=FPG27", "MASS=FPG27", "XFEM=XFEM4320", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.HE8, gauss=("RIGI=FPG27", "MASS=FPG27")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(MT.T10, gauss=("XINT=FPG15", "NOEU=NOEU")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("FPG4=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )


# ------------------------------------------------------------
class HM_PENTA15_XH3(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.PENTA15
    nodes = (
        SetOfNodes("EN6", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN5", (1, 2, 3, 4, 5, 6)),
    )
    elrefe = (
        ElrefeLoc(
            MT.P15,
            gauss=("RIGI=FPG21", "MASS=FPG21", "XFEM=XFEM2160", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.T10, gauss=("XINT=FPG15", "NOEU=NOEU")),
        ElrefeLoc(MT.PE6, gauss=("RIGI=FPG21", "MASS=FPG21")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("MASS=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )


# ------------------------------------------------------------
class HM_PYRAM13_XH3(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.PYRAM13
    nodes = (SetOfNodes("EN6", (6, 7, 8, 9, 10, 11, 12, 13)), SetOfNodes("EN5", (1, 2, 3, 4, 5)))
    elrefe = (
        ElrefeLoc(
            MT.P13,
            gauss=("RIGI=FPG10", "MASS=FPG10", "XFEM=XFEM1620", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.T10, gauss=("XINT=FPG15", "NOEU=NOEU")),
        ElrefeLoc(MT.PY5, gauss=("RIGI=FPG10", "MASS=FPG10")),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("MASS=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )


# ------------------------------------------------------------
class HM_TETRA10_XH3(HM_HEXA20_XH1):

    """Please document this element"""

    meshType = MT.TETRA10
    nodes = (SetOfNodes("EN6", (5, 6, 7, 8, 9, 10)), SetOfNodes("EN5", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.T10,
            gauss=("RIGI=FPG5", "MASS=FPG5", "XINT=FPG15", "XFEM=XFEM810", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TE4, gauss=("XINT=FPG15", "MASS=FPG15")),
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "FPG4=FPG4", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12"),
        ),
        ElrefeLoc(MT.TR3, gauss=("MASS=FPG4", "NOEU=NOEU", "FPG6=FPG6", "FPG7=FPG7", "XCON=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
    )
