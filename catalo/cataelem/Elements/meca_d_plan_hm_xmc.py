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
# CATALOGUE DES ELEMENTS 2D MULTI HM-X-FEM CONTACT


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
        ("EN1", ("DX", "DY", "H1X", "H1Y", "H2X", "H2Y")),
        ("EN2", ("DX", "DY", "H1X", "H1Y", "H2X", "H2Y", "H3X", "H3Y")),
        (
            "EN3",
            (
                "DX",
                "DY",
                "PRE1",
                "H1X",
                "H1Y",
                "H1PRE1",
                "H2X",
                "H2Y",
                "H2PRE1",
                "PRE_FLU",
                "LAG_FLI",
                "LAG_FLS",
                "LAGS_C",
                "LAGS_F1",
                "PR2_FLU",
                "LA2_FLI",
                "LA2_FLS",
                "D1X",
                "D1Y",
            ),
        ),
        (
            "EN4",
            (
                "DX",
                "DY",
                "PRE1",
                "H1X",
                "H1Y",
                "H1PRE1",
                "H2X",
                "H2Y",
                "H2PRE1",
                "H3X",
                "H3Y",
                "H3PRE1",
                "PRE_FLU",
                "LAG_FLI",
                "LAG_FLS",
                "LAGS_C",
                "LAGS_F1",
                "PR2_FLU",
                "LA2_FLI",
                "LA2_FLS",
                "D1X",
                "D1Y",
                "PR3_FLU",
                "LA3_FLI",
                "LA3_FLS",
                "V11",
                "V12",
            ),
        ),
        (
            "EN5",
            (
                "DX",
                "DY",
                "PRE1",
                "H1X",
                "H1Y",
                "H1PRE1",
                "H2X",
                "H2Y",
                "H2PRE1",
                "PRE_FLU",
                "LAG_FLI",
                "LAG_FLS",
                "LAGS_C",
                "LAGS_F1",
                "LAG2_C",
                "LAG2_F1",
                "LAG3_C",
                "LAG3_F1",
                "PR2_FLU",
                "LA2_FLI",
                "LA2_FLS",
                "D1X",
                "D1Y",
                "D2X",
                "D2Y",
                "D3X",
                "D3Y",
            ),
        ),
        (
            "EN6",
            (
                "DX",
                "DY",
                "PRE1",
                "H1X",
                "H1Y",
                "H1PRE1",
                "H2X",
                "H2Y",
                "H2PRE1",
                "H3X",
                "H3Y",
                "H3PRE1",
                "PRE_FLU",
                "LAG_FLI",
                "LAG_FLS",
                "LAGS_C",
                "LAGS_F1",
                "LAG2_C",
                "LAG2_F1",
                "LAG3_C",
                "LAG3_F1",
                "PR2_FLU",
                "LA2_FLI",
                "LA2_FLS",
                "D1X",
                "D1Y",
                "D2X",
                "D2Y",
                "D3X",
                "D3Y",
                "PR3_FLU",
                "LA3_FLI",
                "LA3_FLS",
                "V11",
                "V12",
                "V21",
                "V22",
                "V31",
                "V32",
            ),
        ),
    ),
)


DDL_MECC = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY"))


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


NFORCER = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY"))


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


EGGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELGA", location="XFEM", components=("X", "Y"))


ENGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


XFGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELGA", location="XFEM", components=("X", "Y"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


CTEMPSC = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST", "DELTAT"))


STANO_I = LocatedComponents(phys=PHY.N120_I, type="ELNO", components=("X1",))


E65NEUTR = LocatedComponents(phys=PHY.N1360R, type="ELEM", components=("X[65]",))


E14NEUI = LocatedComponents(phys=PHY.N240_I, type="ELEM", components=("X[14]",))


E52NEUTR = LocatedComponents(phys=PHY.N2448R, type="ELEM", components=("X[52]",))


E24NEUI = LocatedComponents(phys=PHY.N512_I, type="ELEM", components=("X[24]",))


E21NEUTI = LocatedComponents(phys=PHY.N720_I, type="ELEM", components=("X[21]",))


E26NEUTR = LocatedComponents(phys=PHY.N816_R, type="ELEM", components=("X[26]",))


E14NEUTI = LocatedComponents(phys=PHY.N960_I, type="ELEM", components=("X[14]",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="XFEM", components=("X[30]",))


E1NEUTK = LocatedComponents(phys=PHY.NEUT_K24, type="ELEM", components=("Z1",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="XFEM", components=("X[30]",))


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


E2NEUTK = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z1",))


CONTX_R = LocatedComponents(
    phys=PHY.XCONTAC,
    type="ELEM",
    components=(
        "RHON",
        "MU",
        "RHOTK",
        "INTEG",
        "COECH",
        "COSTCO",
        "COSTFR",
        "COPECO",
        "COPEFR",
        "RELA",
    ),
)


CFLUXF = LocatedComponents(phys=PHY.FTHM_F, type="ELEM", components=("PFLUF",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class TemplateElement(Element):
    """Only a template to shared definitions of options"""

    calculs = (
        OP.CHAR_MECA_PESA_R(
            te=588,
            para_in=(
                (OP.CHAR_MECA_PESA_R.PCNSETO, LC.E144NEUI),
                (OP.CHAR_MECA_PESA_R.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_PESA_R.PHEAVTO, E24NEUI),
                (OP.CHAR_MECA_PESA_R.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_PESA_R.PLONCHA, LC.E10NEUTI),
                (OP.CHAR_MECA_PESA_R.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_PESA_R.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PPINTTO, LC.E24NEUTR),
                (OP.CHAR_MECA_PESA_R.PPMILTO, LC.E88NEUTR),
                (OP.CHAR_MECA_PESA_R.PSTANO, STANO_I),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FLUX_F(
            te=559,
            para_in=(
                (SP.PFLUXF, CFLUXF),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTR, CTEMPSC),
                (OP.CHAR_MECA_FLUX_F.PLST, LC.N1NEUT_R),
                (OP.CHAR_MECA_FLUX_F.PLSN, LC.N1NEUT_R),
                (SP.PHEAVNO, LC.FISNO_I),
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
                (OP.FORC_NODA.PBASLOR, LC.N6NEUT_R),
                (OP.FORC_NODA.PCNSETO, LC.E144NEUI),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, ECONTPG),
                (SP.PDEPLAR, DDL_MECA),
                (OP.FORC_NODA.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.FORC_NODA.PHEAVTO, E24NEUI),
                (OP.FORC_NODA.PHEA_NO, LC.N5NEUTI),
                (OP.FORC_NODA.PLONCHA, LC.E10NEUTI),
                (OP.FORC_NODA.PLSN, LC.N1NEUT_R),
                (OP.FORC_NODA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.FORC_NODA.PPINTTO, LC.E24NEUTR),
                (OP.FORC_NODA.PPMILTO, LC.E88NEUTR),
                (OP.FORC_NODA.PSTANO, STANO_I),
                (SP.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_CONT(
            te=557,
            para_in=(
                (OP.CHAR_MECA_CONT.PAINTER, E65NEUTR),
                (OP.CHAR_MECA_CONT.PBASECO, E52NEUTR),
                (OP.CHAR_MECA_CONT.PCFACE, E21NEUTI),
                (SP.PCOHES, LC.E5NEUTR),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PDONCO, CONTX_R),
                (SP.PGEOMER, NGEOMER),
                (SP.PINDCOI, LC.I1NEUT_I),
                (OP.CHAR_MECA_CONT.PLONGCO, LC.E3NEUTI),
                (OP.CHAR_MECA_CONT.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_CONT.PPINTER, E26NEUTR),
                (OP.CHAR_MECA_CONT.PSTANO, STANO_I),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.CHAR_MECA_CONT.PCOMPOR, LC.CCOMPOR),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.CHAR_MECA_CONT.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_CONT.PFISNO, LC.FISNO_I),
                (OP.CHAR_MECA_CONT.PHEA_FA, E14NEUI),
                (SP.PHEAVNO, LC.FISNO_I),
                (SP.PFISCO, LC.FISCO_I),
                (OP.CHAR_MECA_CONT.PLSN, LC.N1NEUT_R),
            ),
            para_out=((SP.PVECTCR, MVECTUR), (SP.PVECTFR, MVECTUR)),
        ),
        OP.CHAR_MECA_CONT_M(
            te=557,
            para_in=(
                (OP.CHAR_MECA_CONT_M.PAINTER, E65NEUTR),
                (OP.CHAR_MECA_CONT_M.PBASECO, E52NEUTR),
                (OP.CHAR_MECA_CONT_M.PCFACE, E21NEUTI),
                (SP.PCOHES, LC.N5NEUTR),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PDONCO, CONTX_R),
                (SP.PGEOMER, NGEOMER),
                (SP.PINDCOI, LC.I1NEUT_I),
                (OP.CHAR_MECA_CONT_M.PLONGCO, LC.E3NEUTI),
                (OP.CHAR_MECA_CONT_M.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_CONT_M.PPINTER, E26NEUTR),
                (OP.CHAR_MECA_CONT_M.PSTANO, STANO_I),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.CHAR_MECA_CONT_M.PCOMPOR, LC.CCOMPOR),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.CHAR_MECA_CONT_M.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_CONT_M.PFISNO, LC.FISNO_I),
                (OP.CHAR_MECA_CONT_M.PHEA_FA, E14NEUI),
                (SP.PHEAVNO, LC.FISNO_I),
                (SP.PFISCO, LC.FISCO_I),
                (OP.CHAR_MECA_CONT_M.PLSN, LC.N1NEUT_R),
            ),
            para_out=((SP.PVECTCR, MVECTUR), (SP.PVECTFR, MVECTUR)),
        ),
        OP.FULL_MECA(
            te=588,
            para_in=(
                (OP.FULL_MECA.PBASLOR, LC.N6NEUT_R),
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCNSETO, LC.E144NEUI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (OP.FULL_MECA.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (SP.PHEAVNO, LC.FISNO_I),
                (OP.FULL_MECA.PHEAVTO, E24NEUI),
                (OP.FULL_MECA.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.FULL_MECA.PLONCHA, LC.E10NEUTI),
                (OP.FULL_MECA.PLSN, LC.N1NEUT_R),
                (OP.FULL_MECA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.FULL_MECA.PPINTTO, LC.E24NEUTR),
                (OP.FULL_MECA.PPMILTO, LC.E88NEUTR),
                (OP.FULL_MECA.PSTANO, STANO_I),
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
        OP.INDL_ELGA(
            te=30,
            para_in=(
                (OP.INDL_ELGA.PCOMPOR, LC.CCOMPOR),
                (OP.INDL_ELGA.PCONTPR, ESIGMPG),
                (SP.PMATERC, LC.CMATERC),
                (OP.INDL_ELGA.PVARIPR, ZVARIPG),
            ),
            para_out=((SP.PINDLOC, LC.EGINDLO),),
        ),
        OP.INIT_MAIL_VOIS(te=99, para_out=((OP.INIT_MAIL_VOIS.PVOISIN, LC.EVOISIN),)),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.INI_XFEM_ELNO(
            te=99,
            para_out=(
                (OP.INI_XFEM_ELNO.PBASLOR, LC.N6NEUT_R),
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
                (OP.RAPH_MECA.PBASLOR, LC.N6NEUT_R),
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RAPH_MECA.PCNSETO, LC.E144NEUI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.RAPH_MECA.PHEAVTO, E24NEUI),
                (OP.RAPH_MECA.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.RAPH_MECA.PLONCHA, LC.E10NEUTI),
                (OP.RAPH_MECA.PLSN, LC.N1NEUT_R),
                (OP.RAPH_MECA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RAPH_MECA.PPINTTO, LC.E24NEUTR),
                (OP.RAPH_MECA.PPMILTO, LC.E88NEUTR),
                (OP.RAPH_MECA.PSTANO, STANO_I),
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
        OP.RIGI_CONT(
            te=556,
            para_in=(
                (OP.RIGI_CONT.PAINTER, E65NEUTR),
                (OP.RIGI_CONT.PBASECO, E52NEUTR),
                (OP.RIGI_CONT.PCFACE, E21NEUTI),
                (SP.PCOHES, LC.E5NEUTR),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PDONCO, CONTX_R),
                (SP.PGEOMER, NGEOMER),
                (SP.PINDCOI, LC.I1NEUT_I),
                (OP.RIGI_CONT.PLONGCO, LC.E3NEUTI),
                (OP.RIGI_CONT.PLSN, LC.N1NEUT_R),
                (OP.RIGI_CONT.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_CONT.PSTANO, STANO_I),
                (OP.RIGI_CONT.PPINTER, E26NEUTR),
                (OP.RIGI_CONT.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_CONT.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_CONT.PFISNO, LC.FISNO_I),
                (SP.PHEAVNO, LC.FISNO_I),
                (OP.RIGI_CONT.PHEA_FA, E14NEUI),
                (SP.PFISCO, LC.FISCO_I),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (OP.RIGI_CONT.PCOHESO, LC.E5NEUTR)),
        ),
        OP.RIGI_CONT_M(
            te=556,
            para_in=(
                (OP.RIGI_CONT_M.PAINTER, E65NEUTR),
                (OP.RIGI_CONT_M.PBASECO, E52NEUTR),
                (OP.RIGI_CONT_M.PCFACE, E21NEUTI),
                (SP.PCOHES, LC.N5NEUTR),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PDONCO, CONTX_R),
                (SP.PGEOMER, NGEOMER),
                (SP.PINDCOI, LC.I1NEUT_I),
                (OP.RIGI_CONT_M.PLONGCO, LC.E3NEUTI),
                (OP.RIGI_CONT_M.PLSN, LC.N1NEUT_R),
                (OP.RIGI_CONT_M.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_CONT_M.PSTANO, STANO_I),
                (OP.RIGI_CONT_M.PPINTER, E26NEUTR),
                (OP.RIGI_CONT_M.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_CONT_M.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_CONT_M.PFISNO, LC.FISNO_I),
                (SP.PHEAVNO, LC.FISNO_I),
                (OP.RIGI_CONT_M.PHEA_FA, E14NEUI),
                (SP.PFISCO, LC.FISCO_I),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (OP.RIGI_CONT_M.PCOHESO, LC.N5NEUTR)),
        ),
        OP.RIGI_MECA_TANG(
            te=588,
            para_in=(
                (OP.RIGI_MECA_TANG.PBASLOR, LC.N6NEUT_R),
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCNSETO, LC.E144NEUI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.RIGI_MECA_TANG.PHEAVTO, E24NEUI),
                (OP.RIGI_MECA_TANG.PFISNO, LC.FISNO_I),
                (OP.RIGI_MECA_TANG.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.RIGI_MECA_TANG.PLONCHA, LC.E10NEUTI),
                (OP.RIGI_MECA_TANG.PLSN, LC.N1NEUT_R),
                (OP.RIGI_MECA_TANG.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA_TANG.PPINTTO, LC.E24NEUTR),
                (OP.RIGI_MECA_TANG.PPMILTO, LC.E88NEUTR),
                (OP.RIGI_MECA_TANG.PSTANO, STANO_I),
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
        OP.SIEF_ELGA(
            te=261,
            para_in=(
                (OP.SIEF_ELGA.PBASLOR, LC.N6NEUT_R),
                (SP.PCAMASS, LC.CCAMA2D),
                (OP.SIEF_ELGA.PCNSETO, LC.E144NEUI),
                (OP.SIEF_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (OP.SIEF_ELGA.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.SIEF_ELGA.PHEAVTO, E24NEUI),
                (OP.SIEF_ELGA.PLONCHA, LC.E10NEUTI),
                (OP.SIEF_ELGA.PLSN, LC.N1NEUT_R),
                (OP.SIEF_ELGA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.SIEF_ELGA.PPINTTO, LC.E24NEUTR),
                (OP.SIEF_ELGA.PPMILTO, LC.E88NEUTR),
                (OP.SIEF_ELGA.PSTANO, STANO_I),
                (OP.SIEF_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.SIEF_ELGA.PCONTRR, ECONTPG),),
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
                (OP.TOPOFA.PAINTTO, LC.E60NEUTR),
                (OP.TOPOFA.PCNSETO, LC.E144NEUI),
                (SP.PDECOU, LC.E1NEUK8),
                (SP.PFISCO, LC.FISCO_I),
                (SP.PGEOMER, NGEOMER),
                (SP.PGRADLN, LC.N2NEUT_R),
                (SP.PGRADLT, LC.N2NEUT_R),
                (OP.TOPOFA.PHEAVTO, E24NEUI),
                (OP.TOPOFA.PLONCHA, LC.E10NEUTI),
                (OP.TOPOFA.PLSN, LC.N1NEUT_R),
                (OP.TOPOFA.PLST, LC.N1NEUT_R),
                (OP.TOPOFA.PPINTTO, LC.E24NEUTR),
                (OP.TOPOFA.PPMILTO, LC.E88NEUTR),
                (OP.TOPOFA.PSTANO, STANO_I),
            ),
            para_out=(
                (OP.TOPOFA.PAINTER, E65NEUTR),
                (OP.TOPOFA.PBASECO, E52NEUTR),
                (OP.TOPOFA.PCFACE, E21NEUTI),
                (SP.PGESCLA, E26NEUTR),
                (OP.TOPOFA.PHEAVFA, E14NEUTI),
                (OP.TOPOFA.PLONGCO, LC.E3NEUTI),
                (OP.TOPOFA.PPINTER, E26NEUTR),
            ),
        ),
        OP.TOPONO(
            te=120,
            para_in=(
                (OP.TOPONO.PCNSETO, LC.E144NEUI),
                (SP.PFISCO, LC.FISCO_I),
                (OP.TOPONO.PFISNO, LC.FISNO_I),
                (OP.TOPONO.PHEAVFA, E14NEUTI),
                (OP.TOPONO.PHEAVTO, E24NEUI),
                (SP.PLEVSET, LC.N1NEUT_R),
                (OP.TOPONO.PLONCHA, LC.E10NEUTI),
                (OP.TOPONO.PLONGCO, LC.E3NEUTI),
            ),
            para_out=(
                (OP.TOPONO.PHEA_FA, E14NEUI),
                (OP.TOPONO.PHEA_NO, LC.N5NEUTI),
                (OP.TOPONO.PHEA_SE, E24NEUI),
            ),
        ),
        OP.TOPOSE(
            te=514,
            para_in=((SP.PFISCO, LC.FISCO_I), (SP.PGEOMER, NGEOMER), (SP.PLEVSET, LC.N1NEUT_R)),
            para_out=(
                (OP.TOPOSE.PAINTTO, LC.E60NEUTR),
                (OP.TOPOSE.PCNSETO, LC.E144NEUI),
                (OP.TOPOSE.PHEAVTO, E24NEUI),
                (OP.TOPOSE.PLONCHA, LC.E10NEUTI),
                (OP.TOPOSE.PPINTTO, LC.E24NEUTR),
                (OP.TOPOSE.PPMILTO, LC.E88NEUTR),
                (OP.TOPOSE.PJONCNO, LC.E20NEUTI),
            ),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D), (OP.TOU_INI_ELEM.PNEUT_K8, E2NEUTK)),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PDOMMAG, LC.EDOMGGA),
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
        OP.XCVBCA(
            te=558,
            para_in=(
                (OP.XCVBCA.PAINTER, E65NEUTR),
                (OP.XCVBCA.PBASECO, E52NEUTR),
                (OP.XCVBCA.PCFACE, E21NEUTI),
                (SP.PCOHES, LC.E5NEUTR),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PDONCO, CONTX_R),
                (SP.PGEOMER, NGEOMER),
                (SP.PINDCOI, LC.I1NEUT_I),
                (OP.XCVBCA.PLONGCO, LC.E3NEUTI),
                (OP.XCVBCA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.XCVBCA.PPINTER, E26NEUTR),
                (OP.XCVBCA.PCOMPOR, LC.CCOMPOR),
                (OP.XCVBCA.PHEA_NO, LC.N5NEUTI),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PHEAVNO, LC.FISNO_I),
                (OP.XCVBCA.PFISNO, LC.FISNO_I),
                (OP.XCVBCA.PHEA_FA, E14NEUI),
                (SP.PFISCO, LC.FISCO_I),
                (OP.XCVBCA.PLSN, LC.N1NEUT_R),
                (OP.XCVBCA.PSTANO, STANO_I),
            ),
            para_out=((OP.XCVBCA.PCOHESO, LC.E5NEUTR), (SP.PINCOCA, LC.I1NEUT_I)),
        ),
        OP.XCVBCA_MORTAR(
            te=558,
            para_in=(
                (OP.XCVBCA_MORTAR.PAINTER, E65NEUTR),
                (OP.XCVBCA_MORTAR.PBASECO, E52NEUTR),
                (OP.XCVBCA_MORTAR.PCFACE, E21NEUTI),
                (SP.PCOHES, LC.N5NEUTR),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PDONCO, CONTX_R),
                (SP.PGEOMER, NGEOMER),
                (SP.PINDCOI, LC.I1NEUT_I),
                (OP.XCVBCA_MORTAR.PLONGCO, LC.E3NEUTI),
                (OP.XCVBCA_MORTAR.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.XCVBCA_MORTAR.PPINTER, E26NEUTR),
                (OP.XCVBCA_MORTAR.PCOMPOR, LC.CCOMPOR),
                (OP.XCVBCA_MORTAR.PHEA_NO, LC.N5NEUTI),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PHEAVNO, LC.FISNO_I),
                (OP.XCVBCA_MORTAR.PFISNO, LC.FISNO_I),
                (OP.XCVBCA_MORTAR.PHEA_FA, E14NEUI),
                (SP.PFISCO, LC.FISCO_I),
                (OP.XCVBCA_MORTAR.PLSN, LC.N1NEUT_R),
                (OP.XCVBCA_MORTAR.PSTANO, STANO_I),
            ),
            para_out=((OP.XCVBCA.PCOHESO, LC.N5NEUTR), (SP.PINCOCA, LC.I1NEUT_I)),
        ),
    )


# ------------------------------------------------------------
class HM_DPTR6_XH2C(TemplateElement):
    """Please document this element"""

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN1", (4, 5, 6)), SetOfNodes("EN3", (1, 2, 3)))
    elrefe = (
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "XINT=FPG12", "XFEM=XFEM72", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG6", "MASS=FPG6")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
    )
    calculs = (OP.CHAR_MECA_CONT_M(te=-1), OP.RIGI_CONT_M(te=-1), OP.XCVBCA_MORTAR(te=-1))


# ------------------------------------------------------------
class HM_DPTR6_XH3C(TemplateElement):
    """Please document this element"""

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN2", (4, 5, 6)), SetOfNodes("EN4", (1, 2, 3)))
    elrefe = (
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "XINT=FPG12", "XFEM=XFEM72", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG6", "MASS=FPG6")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
    )
    calculs = (OP.CHAR_MECA_CONT_M(te=-1), OP.RIGI_CONT_M(te=-1), OP.XCVBCA_MORTAR(te=-1))


# ------------------------------------------------------------
class HM_DPQ8_XH2C(TemplateElement):
    """Please document this element"""

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN1", (5, 6, 7, 8)), SetOfNodes("EN3", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9", "XFEM=XFEM144", "FPG1=FPG1"), mater=("XFEM",)
        ),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6", "XINT=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
    )
    calculs = (OP.CHAR_MECA_CONT_M(te=-1), OP.RIGI_CONT_M(te=-1), OP.XCVBCA_MORTAR(te=-1))


# ------------------------------------------------------------
class HM_DPQ8_XH3C(TemplateElement):
    """Please document this element"""

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8)), SetOfNodes("EN4", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9", "XFEM=XFEM144", "FPG1=FPG1"), mater=("XFEM",)
        ),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6", "XINT=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
    )
    calculs = (OP.CHAR_MECA_CONT_M(te=-1), OP.RIGI_CONT_M(te=-1), OP.XCVBCA_MORTAR(te=-1))


# ------------------------------------------------------------
class HM_DPTR6_XH2C3(TemplateElement):
    """Please document this element"""

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN1", (4, 5, 6)), SetOfNodes("EN5", (1, 2, 3)))
    elrefe = (
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "XINT=FPG12", "XFEM=XFEM72", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG6", "MASS=FPG6")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
    )
    calculs = (OP.CHAR_MECA_CONT(te=-1), OP.RIGI_CONT(te=-1), OP.XCVBCA(te=-1))


# ------------------------------------------------------------
class HM_DPTR6_XH3C3(TemplateElement):
    """Please document this element"""

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN2", (4, 5, 6)), SetOfNodes("EN6", (1, 2, 3)))
    elrefe = (
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "XINT=FPG12", "XFEM=XFEM72", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG6", "MASS=FPG6")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
    )
    calculs = (OP.CHAR_MECA_CONT(te=-1), OP.RIGI_CONT(te=-1), OP.XCVBCA(te=-1))


# ------------------------------------------------------------
class HM_DPQ8_XH2C3(TemplateElement):
    """Please document this element"""

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN1", (5, 6, 7, 8)), SetOfNodes("EN5", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9", "XFEM=XFEM144", "FPG1=FPG1"), mater=("XFEM",)
        ),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6", "XINT=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
    )
    calculs = (OP.CHAR_MECA_CONT(te=-1), OP.RIGI_CONT(te=-1), OP.XCVBCA(te=-1))


# ------------------------------------------------------------
class HM_DPQ8_XH3C3(TemplateElement):
    """Please document this element"""

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8)), SetOfNodes("EN6", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9", "XFEM=XFEM144", "FPG1=FPG1"), mater=("XFEM",)
        ),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6", "XINT=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "MASS=FPG4", "FPG2=FPG2", "FPG3=FPG3", "FPG4=FPG4")),
    )
    calculs = (OP.CHAR_MECA_CONT(te=-1), OP.RIGI_CONT(te=-1), OP.XCVBCA(te=-1))


del TemplateElement
