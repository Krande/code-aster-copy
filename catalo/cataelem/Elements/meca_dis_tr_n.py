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


CABSCUR = LocatedComponents(phys=PHY.ABSC_R, type="ELEM", components=("ABSC1",))


CCADISA = LocatedComponents(phys=PHY.CADISA_R, type="ELEM", components=("A[36]",))


CCADISK = LocatedComponents(phys=PHY.CADISK_R, type="ELEM", components=("K[36]",))


CCADISM = LocatedComponents(phys=PHY.CADISM_R, type="ELEM", components=("M[36]",))


CCAORIE = LocatedComponents(phys=PHY.CAORIE_R, type="ELEM", components=("ALPHA", "BETA", "GAMMA"))

NDEPLAC = LocatedComponents(
    phys=PHY.DEPL_C, type="ELNO", components=("DX", "DY", "DZ", "DRX", "DRY", "DRZ")
)


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ", "DRX", "DRY", "DRZ")
)


EENERR = LocatedComponents(
    phys=PHY.ENER_R, type="ELEM", components=("TOTALE", "DX", "DY", "DZ", "DRX", "DRY", "DRZ")
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOM_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z")
)


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


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


ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))


CEPSINR = LocatedComponents(phys=PHY.EPSI_R, type="ELGA", location="RIGI", components=("EPXX",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MECA_DIS_TR_N(Element):
    """Please document this element"""

    meshType = MT.POI1
    elrefe = (ElrefeLoc(MT.PO1, gauss=("RIGI=NOEU", "FPG1=FPG1"), mater=("RIGI", "FPG1")),)
    calculs = (
        OP.ADD_SIGM(
            te=581,
            para_in=((SP.PEPCON1, EEFGEGA), (SP.PEPCON2, EEFGEGA)),
            para_out=((SP.PEPCON3, EEFGEGA),),
        ),
        OP.AMOR_MECA(
            te=41,
            para_in=(
                (SP.PCADISA, CCADISA),
                (OP.AMOR_MECA.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (SP.PMATERC, LC.CMATERC),
                (SP.PRIGIEL, MMATUUR),
                (OP.AMOR_MECA.PVARCPR, LC.ZVARCPG),
                (OP.AMOR_MECA.PCOMPOR, LC.CCOMPOR),
                (SP.PNONLIN, LC.ENONLIN),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.CHAR_MECA_EPSA_R(
            te=99,
            para_in=((SP.PGEOMER, NGEOMER), (OP.CHAR_MECA_EPSA_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        #       -- te0580 : ne resout que le cas trivial : EPXX=0.
        OP.CHAR_MECA_EPSI_R(
            te=580, para_in=((SP.PEPSINR, CEPSINR),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.CHAR_MECA_PESA_R(
            te=43,
            para_in=(
                (SP.PCADISM, CCADISM),
                (OP.CHAR_MECA_PESA_R.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_ROTA_R(
            te=43,
            para_in=((SP.PCINFDI, LC.CCINFDI), (SP.PROTATR, LC.CROTATR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        #       -- les elements discrets ne peuvent pas calculer de dilatation thermique => te0099
        OP.CHAR_MECA_TEMP_R(
            te=99,
            para_in=((SP.PGEOMER, NGEOMER), (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        #       -- les elements discrets ne peuvent pas calculer de dilatation "hydratation/sechage" => te0099
        OP.CHAR_MECA_HYDR_R(
            te=99, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.CHAR_MECA_SECH_R(
            te=99, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.COOR_ELGA(
            te=478, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.ECIN_ELEM(
            te=44,
            para_in=(
                (SP.PCADISM, CCADISM),
                (OP.ECIN_ELEM.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (SP.PDEPLAR, DDL_MECA),
                (SP.POMEGA2, LC.COMEG2R),
                (OP.ECIN_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVITESR, DDL_MECA),
            ),
            para_out=((SP.PENERCR, EENERR),),
        ),
        OP.EFGE_ELGA(
            te=546,
            para_in=((SP.PSIEFR, EEFGEGA),),
            para_out=((SP.PEFGEC, EEFGEGC), (SP.PEFGER, EEFGEGA)),
        ),
        OP.EFGE_ELNO(
            te=185,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.EFGE_ELNO.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.EFGE_ELNO.PCONTRR, EEFGEGA),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PNONLIN, LC.ENONLIN),
                (OP.EFGE_ELNO.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PEFFORC, EEFGENC), (OP.EFGE_ELNO.PEFFORR, EEFGENO)),
        ),
        OP.EPOT_ELEM(
            te=44,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.EPOT_ELEM.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (SP.PDEPLAR, DDL_MECA),
                (OP.EPOT_ELEM.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.EPOT_ELEM.PENERDR, EENERR),),
        ),
        OP.FONL_NOEU(
            te=39,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.FONL_NOEU.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.FONL_NOEU.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FORC_NODA(
            te=39,
            para_in=(
                (OP.FORC_NODA.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (SP.PSIEFR, EEFGEGA),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=47,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.FULL_MECA.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PITERAT, LC.CITERAT),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, EEFGEGA),
                (SP.PDEPENT, DDL_MECA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, LC.CINSTPR),
                (SP.PINSTPR, LC.CINSTPR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARIMR, ZVARIPG),
                (SP.PVITENT, DDL_MECA),
                (SP.PVITPLU, DDL_MECA),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, EEFGEGA),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.FULL_MECA_ELAS(
            te=47,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.FULL_MECA_ELAS.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PITERAT, LC.CITERAT),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.FULL_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA_ELAS.PCONTMR, EEFGEGA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, LC.CINSTPR),
                (SP.PINSTPR, LC.CINSTPR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (OP.FULL_MECA_ELAS.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA_ELAS.PCONTPR, EEFGEGA),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA_ELAS.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_FLUI_STRU(
            te=41,
            para_in=(
                (SP.PABSCUR, CABSCUR),
                (SP.PCADISM, CCADISM),
                (OP.MASS_FLUI_STRU.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.MASS_FLUI_STRU.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_INER(
            te=45,
            para_in=(
                (SP.PCADISM, CCADISM),
                (OP.MASS_INER.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (SP.PGEOMER, NGEOMER),
                (OP.MASS_INER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMASSINE, LC.EMASSINE),),
        ),
        OP.MASS_MECA(
            te=41,
            para_in=(
                (SP.PCADISM, CCADISM),
                (OP.MASS_MECA.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.MASS_MECA_DIAG(
            te=41,
            para_in=(
                (SP.PCADISM, CCADISM),
                (OP.MASS_MECA_DIAG.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.MASS_MECA_DIAG.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.MASS_MECA_EXPLI(
            te=41,
            para_in=(
                (SP.PCADISM, CCADISM),
                (OP.MASS_MECA_EXPLI.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.MASS_MECA_EXPLI.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MECA_GYRO(
            te=9,
            para_in=(
                (SP.PCADISM, CCADISM),
                (OP.MECA_GYRO.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (SP.PGEOMER, NGEOMER),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.M_GAMMA(
            te=41,
            para_in=(
                (SP.PACCELR, DDL_MECA),
                (SP.PCADISM, CCADISM),
                (OP.M_GAMMA.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
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
            te=405,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PCOURAN, LC.ECOURAN),),
        ),
        OP.RAPH_MECA(
            te=47,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.RAPH_MECA.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PITERAT, LC.CITERAT),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, EEFGEGA),
                (SP.PDEPENT, DDL_MECA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, LC.CINSTPR),
                (SP.PINSTPR, LC.CINSTPR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARCPR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARIMR, ZVARIPG),
                (SP.PVITENT, DDL_MECA),
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
            te=39,
            para_in=((SP.PCINFDI, LC.CCINFDI), (SP.PREFCO, LC.CRESEFM)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.REPERE_LOCAL(
            te=135,
            para_in=((OP.REPERE_LOCAL.PCAORIE, CCAORIE),),
            para_out=((SP.PREPLO1, LC.CGEOM3D), (SP.PREPLO2, LC.CGEOM3D), (SP.PREPLO3, LC.CGEOM3D)),
        ),
        OP.REST_ECRO(te=99, para_out=((OP.REST_ECRO.PVARIPR, ZVARIPG),)),
        OP.RIGI_FLUI_STRU(
            te=41,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.RIGI_FLUI_STRU.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.RIGI_FLUI_STRU.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_GYRO(
            te=464,
            para_in=(
                (SP.PCADISM, CCADISM),
                (OP.RIGI_GYRO.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (SP.PGEOMER, NGEOMER),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.RIGI_MECA(
            te=41,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.RIGI_MECA.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.RIGI_MECA_ELAS(
            te=47,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.RIGI_MECA_ELAS.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PITERAT, LC.CITERAT),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.RIGI_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_ELAS.PCONTMR, EEFGEGA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, LC.CINSTPR),
                (SP.PINSTPR, LC.CINSTPR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (OP.RIGI_MECA_ELAS.PVARIMR, ZVARIPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_HYST(
            te=41,
            para_in=(
                (SP.PCINFDI, LC.CCINFDI),
                (SP.PRIGIEL, MMATUUR),
                (OP.RIGI_MECA_HYST.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUC, MMATUUC),),
        ),
        # OP.RIGI_MECA_RO      issue19398
        OP.RIGI_MECA_TANG(
            te=47,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.RIGI_MECA_TANG.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PITERAT, LC.CITERAT),
                (SP.PCINFDI, LC.CCINFDI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, EEFGEGA),
                (SP.PDEPENT, DDL_MECA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, LC.CINSTPR),
                (SP.PINSTPR, LC.CINSTPR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARIMR, ZVARIPG),
                (SP.PVITENT, DDL_MECA),
                (SP.PVITPLU, DDL_MECA),
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
            te=42,
            para_in=(
                (SP.PCADISK, CCADISK),
                (OP.SIEF_ELGA.PCAORIE, CCAORIE),
                (SP.PCINFDI, LC.CCINFDI),
                (SP.PDEPLAR, DDL_MECA),
                (OP.SIEF_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PCONTRC, EEFGEGC), (OP.SIEF_ELGA.PCONTRR, EEFGEGA)),
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
            te=4, para_in=((SP.PVARIGR, ZVARIPG),), para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),)
        ),
    )
