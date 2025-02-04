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
# CATALOGUE DES ELEMENTS 2D HM-X-FEM


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
        ("EN1", ("DX", "DY", "PRE1", "H1X", "H1Y", "H1PRE1")),
        ("EN2", ("DX", "DY", "H1X", "H1Y")),
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


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


EGGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELGA", location="XFEM", components=("X", "Y"))


ENGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


STANO_I = LocatedComponents(phys=PHY.N120_I, type="ELNO", components=("X1",))


E6NEUTI = LocatedComponents(phys=PHY.N512_I, type="ELEM", components=("X[6]",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="XFEM", components=("X[30]",))


E1NEUTK = LocatedComponents(phys=PHY.NEUT_K24, type="ELEM", components=("Z1",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="XFEM", components=("X[30]",))


CPRESSF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PRES", "CISA"))


EPRESNO = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES", "CISA"))


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


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class HM_DPQ8_XH(Element):

    """Please document this element"""

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.QU8, gauss=("RIGI=FPG9", "MASS=FPG9", "XFEM=XFEM72", "FPG1=FPG1"), mater=("XFEM",)
        ),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6", "XINT=FPG12")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG9", "MASS=FPG9")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "MASS=FPG4")),
    )
    calculs = (
        OP.CHAR_MECA_PESA_R(
            te=588,
            para_in=(
                (OP.CHAR_MECA_PESA_R.PBASLOR, LC.N6NEUT_R),
                (OP.CHAR_MECA_PESA_R.PCNSETO, LC.E36NEUI),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_PESA_R.PHEAVTO, E6NEUTI),
                (OP.CHAR_MECA_PESA_R.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_PESA_R.PLONCHA, LC.E10NEUTI),
                (OP.CHAR_MECA_PESA_R.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_PESA_R.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PPINTTO, LC.E6NEUTR),
                (OP.CHAR_MECA_PESA_R.PPMILTO, LC.E22NEUTR),
                (OP.CHAR_MECA_PESA_R.PSTANO, STANO_I),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_F(
            te=37,
            para_in=(
                (OP.CHAR_MECA_PRES_F.PAINTER, LC.E35NEUTR),
                (OP.CHAR_MECA_PRES_F.PBASECO, LC.E28NEUTR),
                (OP.CHAR_MECA_PRES_F.PCFACE, LC.E9NEUTI),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_PRES_F.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_PRES_F.PLONGCO, LC.E3NEUTI),
                (OP.CHAR_MECA_PRES_F.PLST, LC.N1NEUT_R),
                (OP.CHAR_MECA_PRES_F.PPINTER, LC.E14NEUTR),
                (SP.PPRESSF, CPRESSF),
                (OP.CHAR_MECA_PRES_F.PSTANO, STANO_I),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=37,
            para_in=(
                (OP.CHAR_MECA_PRES_R.PAINTER, LC.E35NEUTR),
                (OP.CHAR_MECA_PRES_R.PBASECO, LC.E28NEUTR),
                (OP.CHAR_MECA_PRES_R.PCFACE, LC.E9NEUTI),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_PRES_R.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_PRES_R.PLONGCO, LC.E3NEUTI),
                (OP.CHAR_MECA_PRES_R.PLST, LC.N1NEUT_R),
                (OP.CHAR_MECA_PRES_R.PPINTER, LC.E14NEUTR),
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
                (OP.FORC_NODA.PBASLOR, LC.N6NEUT_R),
                (OP.FORC_NODA.PCNSETO, LC.E36NEUI),
                (SP.PSIEFR, ECONTPG),
                (SP.PGEOMER, NGEOMER),
                (OP.FORC_NODA.PHEAVTO, E6NEUTI),
                (OP.FORC_NODA.PHEA_NO, LC.N5NEUTI),
                (OP.FORC_NODA.PLONCHA, LC.E10NEUTI),
                (OP.FORC_NODA.PLSN, LC.N1NEUT_R),
                (OP.FORC_NODA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.FORC_NODA.PPINTTO, LC.E6NEUTR),
                (OP.FORC_NODA.PPMILTO, LC.E22NEUTR),
                (OP.FORC_NODA.PSTANO, STANO_I),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=588,
            para_in=(
                (OP.FULL_MECA.PBASLOR, LC.N6NEUT_R),
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCNSETO, LC.E36NEUI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.FULL_MECA.PHEAVTO, E6NEUTI),
                (OP.FULL_MECA.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.FULL_MECA.PLONCHA, LC.E10NEUTI),
                (OP.FULL_MECA.PLSN, LC.N1NEUT_R),
                (OP.FULL_MECA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.FULL_MECA.PPINTTO, LC.E6NEUTR),
                (OP.FULL_MECA.PPMILTO, LC.E22NEUTR),
                (OP.FULL_MECA.PSTANO, STANO_I),
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
                (OP.RAPH_MECA.PCNSETO, LC.E36NEUI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.RAPH_MECA.PHEAVTO, E6NEUTI),
                (OP.RAPH_MECA.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.RAPH_MECA.PLONCHA, LC.E10NEUTI),
                (OP.RAPH_MECA.PLSN, LC.N1NEUT_R),
                (OP.RAPH_MECA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RAPH_MECA.PPINTTO, LC.E6NEUTR),
                (OP.RAPH_MECA.PPMILTO, LC.E22NEUTR),
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
                (OP.RIGI_MECA_TANG.PBASLOR, LC.N6NEUT_R),
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCNSETO, LC.E36NEUI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.RIGI_MECA_TANG.PHEAVTO, E6NEUTI),
                (OP.RIGI_MECA_TANG.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.RIGI_MECA_TANG.PLONCHA, LC.E10NEUTI),
                (OP.RIGI_MECA_TANG.PLSN, LC.N1NEUT_R),
                (OP.RIGI_MECA_TANG.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA_TANG.PPINTTO, LC.E6NEUTR),
                (OP.RIGI_MECA_TANG.PPMILTO, LC.E22NEUTR),
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
                (OP.TOPOFA.PAINTTO, LC.E15NEUTR),
                (OP.TOPOFA.PCNSETO, LC.E36NEUI),
                (SP.PDECOU, LC.E1NEUK8),
                (SP.PGEOMER, NGEOMER),
                (SP.PGRADLN, LC.N2NEUT_R),
                (SP.PGRADLT, LC.N2NEUT_R),
                (OP.TOPOFA.PHEAVTO, E6NEUTI),
                (OP.TOPOFA.PLONCHA, LC.E10NEUTI),
                (OP.TOPOFA.PLSN, LC.N1NEUT_R),
                (OP.TOPOFA.PLST, LC.N1NEUT_R),
                (OP.TOPOFA.PPINTTO, LC.E6NEUTR),
                (OP.TOPOFA.PPMILTO, LC.E22NEUTR),
                (SP.PTYPDIS, LC.E1NEUTI),
            ),
            para_out=(
                (OP.TOPOFA.PAINTER, LC.E35NEUTR),
                (OP.TOPOFA.PBASECO, LC.E28NEUTR),
                (OP.TOPOFA.PCFACE, LC.E9NEUTI),
                (SP.PGESCLA, LC.E14NEUTR),
                (OP.TOPOFA.PLONGCO, LC.E3NEUTI),
                (OP.TOPOFA.PPINTER, LC.E14NEUTR),
            ),
        ),
        OP.TOPONO(
            te=120,
            para_in=(
                (OP.TOPONO.PCNSETO, LC.E36NEUI),
                (OP.TOPONO.PHEAVTO, E6NEUTI),
                (SP.PLEVSET, LC.N1NEUT_R),
                (OP.TOPONO.PLONCHA, LC.E10NEUTI),
            ),
            para_out=((OP.TOPONO.PHEA_NO, LC.N5NEUTI), (OP.TOPONO.PHEA_SE, E6NEUTI)),
        ),
        OP.TOPOSE(
            te=514,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PLEVSET, LC.N1NEUT_R)),
            para_out=(
                (OP.TOPOSE.PAINTTO, LC.E15NEUTR),
                (OP.TOPOSE.PCNSETO, LC.E36NEUI),
                (OP.TOPOSE.PHEAVTO, E6NEUTI),
                (OP.TOPOSE.PLONCHA, LC.E10NEUTI),
                (OP.TOPOSE.PPINTTO, LC.E6NEUTR),
                (OP.TOPOSE.PPMILTO, LC.E22NEUTR),
            ),
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
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),)),
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
class HM_DPTR6_XH(HM_DPQ8_XH):

    """Please document this element"""

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN2", (4, 5, 6)), SetOfNodes("EN1", (1, 2, 3)))
    elrefe = (
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "XINT=FPG12", "XFEM=XFEM36", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG6", "MASS=FPG6")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "MASS=FPG4")),
    )


# ------------------------------------------------------------
class HM_DPQ8M_XH(HM_DPQ8_XH):

    """Please document this element"""

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.QU8, gauss=("RIGI=FPG4", "MASS=FPG9", "XFEM=XFEM72", "FPG1=FPG1"), mater=("XFEM",)
        ),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG9")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4")),
    )


# ------------------------------------------------------------
class HM_DPTR6M_XH(HM_DPQ8_XH):

    """Please document this element"""

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN2", (4, 5, 6)), SetOfNodes("EN1", (1, 2, 3)))
    elrefe = (
        ElrefeLoc(
            MT.TR6, gauss=("RIGI=FPG3", "MASS=FPG6", "XFEM=XFEM36", "FPG1=FPG1"), mater=("XFEM",)
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG3", "MASS=FPG6")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "MASS=FPG4")),
    )


# ------------------------------------------------------------
class HM_DPQ8D_XH(HM_DPQ8_XH):

    """Please document this element"""

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.QU8,
            gauss=("RIGI=NOEU_S", "MASS=NOEU_S", "XFEM=XFEM72", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.QU4, gauss=("RIGI=NOEU_S",)),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )


# ------------------------------------------------------------
class HM_DPTR6D_XH(HM_DPQ8_XH):

    """Please document this element"""

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN2", (4, 5, 6)), SetOfNodes("EN1", (1, 2, 3)))
    elrefe = (
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=NOEU_S", "MASS=NOEU_S", "XFEM=XFEM36", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=NOEU_S",)),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )


# ------------------------------------------------------------
class HM_DPQ8S_XH(HM_DPQ8_XH):

    """Please document this element"""

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(
            MT.QU8,
            gauss=("RIGI=FPG4NOS", "MASS=FPG4", "XFEM=XFEM72", "NOEU_S=NOEU_S", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4NOS", "MASS=FPG4", "NOEU_S=NOEU_S")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )


# ------------------------------------------------------------
class HM_DPTR6S_XH(HM_DPQ8_XH):

    """Please document this element"""

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN2", (4, 5, 6)), SetOfNodes("EN1", (1, 2, 3)))
    elrefe = (
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG3NOS", "MASS=FPG3", "XFEM=XFEM36", "NOEU_S=NOEU_S", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG3NOS", "MASS=FPG3", "NOEU_S=NOEU_S")),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )
