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
import cataelem.Commons.attributes as AT

# ----------------
# Modes locaux :
# ----------------


NDEPLAR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY"))


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("DX", "DY", "SIXX", "SIYY", "SIZZ", "SIXY")), ("EN2", ("DX", "DY"))),
)


EDEPLPG = LocatedComponents(phys=PHY.DEPL_R, type="ELGA", location="RIGI", components=("DX", "DY"))


EENERR = LocatedComponents(phys=PHY.ENER_R, type="ELEM", components=("TOTALE",))


EENERPG = LocatedComponents(phys=PHY.ENER_R, type="ELGA", location="RIGI", components=("TOTALE",))


EENERNO = LocatedComponents(phys=PHY.ENER_R, type="ELNO", components=("TOTALE",))


EDEFOPC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELGA", location="RIGI", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)


EDEFONC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)


CEPSINF = LocatedComponents(
    phys=PHY.EPSI_F, type="ELEM", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)


CEPSINR = LocatedComponents(
    phys=PHY.EPSI_R, type="ELGA", location="RIGI", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)


EDEFOPG = LocatedComponents(
    phys=PHY.EPSI_R, type="ELGA", location="RIGI", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)


EDEFONO = LocatedComponents(
    phys=PHY.EPSI_R, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)


EDFVCPG = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EPTHER_L", "EPTHER_T", "EPTHER_N", "EPSECH", "EPHYDR", "EPPTOT"),
)


EDFVCNO = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELNO",
    components=("EPTHER_L", "EPTHER_T", "EPTHER_N", "EPSECH", "EPHYDR", "EPPTOT"),
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


EGGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "W")
)


ENGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


ECONTPC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELGA", location="RIGI", components=("SIXX", "SIYY", "SIZZ", "SIXY")
)


ECONTNC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY")
)


ECONTPG = LocatedComponents(
    phys=PHY.SIEF_R, type="ELGA", location="RIGI", components=("SIXX", "SIYY", "SIZZ", "SIXY")
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
    phys=PHY.SIEF_R, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY")
)


ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MVECTDR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=NDEPLAR)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MGSPQU8(Element):
    """Please document this element"""

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8)), SetOfNodes("EN1", (1, 2, 3, 4)))
    attrs = ((AT.TYPE_VOISIN, "A2"),)
    elrefe = (
        ElrefeLoc(
            MT.QU8,
            gauss=("RIGI=FPG9", "FPG1=FPG1", "MASS=FPG9", "NOEU=NOEU"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG9", "MASS=FPG9", "NOEU=NOEU")),
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
                (SP.PRIGINS, MMATUNS),
                (OP.AMOR_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.CALC_NOEU_BORD(
            te=290, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PVECTUR, MVECTDR),)
        ),
        OP.CARA_GEOM(
            te=285, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PCARAGE, LC.ECARAGE),)
        ),
        OP.CHAR_MECA_EPSA_R(
            te=421,
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
            te=284,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PEPSINF, CEPSINF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_EPSI_F.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_EPSI_R(
            te=284,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PEPSINR, CEPSINR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_EPSI_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_HYDR_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_HYDR_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=85,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_PTOT_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_PTOT_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.COOR_ELGA(
            te=479, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.ENEL_ELEM(
            te=490,
            para_in=(
                (OP.ENEL_ELEM.PCOMPOR, LC.CCOMPOR),
                (OP.ENEL_ELEM.PCONTPR, ECONTPG),
                (SP.PDEPLR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.ENEL_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.ENEL_ELEM.PVARIPR, ZVARIPG),
            ),
            para_out=((SP.PENERD1, EENERR),),
        ),
        OP.ENEL_ELGA(
            te=575,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (OP.ENEL_ELGA.PCOMPOR, LC.CCOMPOR),
                (OP.ENEL_ELGA.PCONTRR, ECONTPG),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.ENEL_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIGR, ZVARIPG),
            ),
            para_out=((OP.ENEL_ELGA.PENERDR, EENERPG),),
        ),
        OP.ENEL_ELNO(
            te=4, para_in=((OP.ENEL_ELNO.PENERPG, EENERPG),), para_out=((SP.PENERNO, EENERNO),)
        ),
        OP.ENER_TOTALE(
            te=490,
            para_in=(
                (OP.ENER_TOTALE.PCOMPOR, LC.CCOMPOR),
                (OP.ENER_TOTALE.PCONTMR, ECONTPG),
                (OP.ENER_TOTALE.PCONTPR, ECONTPG),
                (SP.PDEPLM, NDEPLAR),
                (SP.PDEPLR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.ENER_TOTALE.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.ENER_TOTALE.PVARIPR, ZVARIPG),
            ),
            para_out=((SP.PENERD1, EENERR),),
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
        OP.EPMQ_ELGA(
            te=335,
            para_in=((OP.EPMQ_ELGA.PDEFORR, EDEFOPG),),
            para_out=((OP.EPMQ_ELGA.PDEFOEQ, LC.EDFEQPG),),
        ),
        OP.EPMQ_ELNO(
            te=335,
            para_in=((OP.EPMQ_ELNO.PDEFORR, EDEFONO),),
            para_out=((OP.EPMQ_ELNO.PDEFOEQ, LC.EDFEQNO),),
        ),
        OP.EPOT_ELEM(
            te=286,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPOT_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPOT_ELEM.PENERDR, EENERR),),
        ),
        OP.EPSI_ELGA(
            te=87,
            para_in=(
                (OP.EPSI_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
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
        OP.EPSP_ELGA(
            te=334,
            para_in=(
                (OP.EPSP_ELGA.PCOMPOR, LC.CCOMPOR),
                (OP.EPSP_ELGA.PCONTRR, ECONTPG),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.EPSP_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIGR, ZVARIPG),
            ),
            para_out=((OP.EPSP_ELGA.PDEFOPG, EDEFOPG),),
        ),
        OP.EPSP_ELNO(
            te=4, para_in=((OP.EPSP_ELNO.PDEFOPG, EDEFOPG),), para_out=((SP.PDEFONO, EDEFONO),)
        ),
        OP.EPVC_ELGA(
            te=529,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPVC_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPVC_ELGA.PDEFOPG, EDFVCPG),),
        ),
        OP.EPVC_ELNO(
            te=4, para_in=((OP.EPVC_ELNO.PDEFOPG, EDFVCPG),), para_out=((SP.PDEFONO, EDFVCNO),)
        ),
        OP.FORC_NODA(
            te=7,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, ECONTPG),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=123,
            para_in=(
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PITERAT, LC.CITERAT),
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
                (OP.FULL_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.FULL_MECA_ELAS(
            te=123,
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
                (SP.PVARIMP, ZVARIPG),
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
        OP.INDIC_ENER(
            te=490,
            para_in=(
                (OP.INDIC_ENER.PCOMPOR, LC.CCOMPOR),
                (OP.INDIC_ENER.PCONTPR, ECONTPG),
                (SP.PDEPLR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.INDIC_ENER.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.INDIC_ENER.PVARIPR, ZVARIPG),
            ),
            para_out=((SP.PENERD1, EENERR), (SP.PENERD2, EENERR)),
        ),
        OP.INDIC_SEUIL(
            te=490,
            para_in=(
                (OP.INDIC_SEUIL.PCOMPOR, LC.CCOMPOR),
                (OP.INDIC_SEUIL.PCONTPR, ECONTPG),
                (SP.PDEPLR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.INDIC_SEUIL.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.INDIC_SEUIL.PVARIPR, ZVARIPG),
            ),
            para_out=((SP.PENERD1, EENERR), (SP.PENERD2, EENERR)),
        ),
        OP.INIT_MAIL_VOIS(te=99, para_out=((OP.INIT_MAIL_VOIS.PVOISIN, LC.EVOISIN),)),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_MECA(
            te=123,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_MECA_DIAG(
            te=123,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA_DIAG.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_MECA_EXPLI(
            te=123,
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
            te=405,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PCOURAN, LC.ECOURAN),),
        ),
        OP.RAPH_MECA(
            te=123,
            para_in=(
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PITERAT, LC.CITERAT),
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
        OP.REPERE_LOCAL(
            te=133,
            para_in=((SP.PCAMASS, LC.CCAMA2D), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PREPLO1, LC.CGEOM2D), (SP.PREPLO2, LC.CGEOM2D)),
        ),
        OP.RIGI_MECA_ELAS(
            te=123,
            para_in=(
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
            te=123,
            para_in=(
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PITERAT, LC.CITERAT),
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
            te=4,
            para_in=((OP.SIEF_ELNO.PCONTRR, ECONTPG), (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PSIEFNOC, ECONTNC), (OP.SIEF_ELNO.PSIEFNOR, ECONTNO)),
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
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),)),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PDEPL_R, EDEPLPG),
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
                (OP.TOU_INI_ELNO.PDOMMAG, LC.EDOMGNO),
                (OP.TOU_INI_ELNO.PEPSI_R, EDEFONO),
                (OP.TOU_INI_ELNO.PGEOM_R, ENGEOM_R),
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
class MGSPTR6(MGSPQU8):
    """Please document this element"""

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN2", (4, 5, 6)), SetOfNodes("EN1", (1, 2, 3)))
    attrs = ((AT.TYPE_VOISIN, "A2"),)
    elrefe = (
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG3", "FPG1=FPG1", "MASS=FPG3", "NOEU=NOEU"),
            mater=("RIGI", "FPG1"),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG3", "MASS=FPG3", "NOEU=NOEU")),
    )
