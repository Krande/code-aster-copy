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

# person_in_charge: mickael.abbas at edf.fr


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
    components=(("EN1", ("DX", "DY", "DZ", "PRES", "PIX", "PIY", "PIZ")),),
)


NDEPLAR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))


EDEFOPC = LocatedComponents(
    phys=PHY.EPSI_C,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"),
)


EDEFONC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)


CEPSINF = LocatedComponents(
    phys=PHY.EPSI_F, type="ELEM", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
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


CEPSINO = LocatedComponents(
    phys=PHY.EPSI_R, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)

EERREUR = LocatedComponents(
    phys=PHY.ERRE_R,
    type="ELEM",
    components=(
        "ERREST",
        "NUEST",
        "SIGCAL",
        "TERMRE",
        "TERMR2",
        "TERMNO",
        "TERMN2",
        "TERMSA",
        "TERMS2",
        "TAILLE",
    ),
)


EERRENO = LocatedComponents(
    phys=PHY.ERRE_R,
    type="ELNO",
    components=(
        "ERREST",
        "NUEST",
        "SIGCAL",
        "TERMRE",
        "TERMR2",
        "TERMNO",
        "TERMN2",
        "TERMSA",
        "TERMS2",
        "TAILLE",
    ),
)


CFORCEF = LocatedComponents(phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ"))


NFORCER = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ"))


EFORCER = LocatedComponents(
    phys=PHY.FORC_R, type="ELGA", location="RIGI", components=("FX", "FY", "FZ")
)

NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOM_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z")
)


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


EMNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[30]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


EREFCO = LocatedComponents(phys=PHY.PREC_R, type="ELEM", components=("SIGM", "EPSI", "PI"))


ESIGMPC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
)


ESIGMNC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
)


ECONTPC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ", "SIP"),
)


ECONTNC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ", "SIP")
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
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ", "SIP"),
)


ECONTNO = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ", "SIP")
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


ESOURCR = LocatedComponents(phys=PHY.SOUR_R, type="ELGA", location="RIGI", components=("SOUR",))

ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

VVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=NDEPLAR)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

VMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=NDEPLAR)


# ------------------------------------------------------------
class MINCOS_HEXA8(Element):
    """Mechanics - 3D - Incompressible - UPO model - HEXA8"""

    meshType = MT.HEXA8
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),)
    elrefe = (
        ElrefeLoc(
            MT.HE8,
            gauss=("RIGI=FPG8", "MASS=FPG8", "NOEU=NOEU", "FPG1=FPG1"),
            mater=("RIGI", "MASS", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.HE8, gauss=("RIGI=FPG8",)),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG4", "NOEU=NOEU")),
    )
    calculs = (
        OP.CALC_G(
            te=222,
            para_in=(
                (SP.PACCELE, NDEPLAR),
                (OP.CALC_G.PCOMPOR, LC.CCOMPOR),
                (SP.PCONTGR, ESIGMPG),
                (OP.CALC_G.PCONTRR, ESIGMPG),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PEPSINR, CEPSINO),
                (SP.PFRVOLU, NFORCER),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, ESIGMNO),
                (OP.CALC_G.PTHETAR, LC.ETHETA),
                (OP.CALC_G.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVITESS, NDEPLAR),
                (OP.CALC_G.PDEG, LC.E1NEUTI),
                (OP.CALC_G.PLAG, LC.CABSLAG),
                (OP.CALC_G.PCER, LC.E1NEUTR),
                (OP.CALC_G.PELI, LC.E2NEUTR),
            ),
            para_out=((SP.PGTHETA, LC.CKGTHET),),
        ),
        OP.CALC_G_F(
            te=222,
            para_in=(
                (SP.PACCELE, NDEPLAR),
                (OP.CALC_G_F.PCOMPOR, LC.CCOMPOR),
                (SP.PCONTGR, ESIGMPG),
                (OP.CALC_G_F.PCONTRR, ESIGMPG),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PEPSINF, CEPSINF),
                (SP.PFFVOLU, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, ESIGMNO),
                (SP.PINSTR, CTEMPSR),
                (OP.CALC_G_F.PTHETAR, LC.ETHETA),
                (OP.CALC_G_F.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVITESS, NDEPLAR),
                (OP.CALC_G_F.PDEG, LC.E1NEUTI),
                (OP.CALC_G_F.PLAG, LC.CABSLAG),
                (OP.CALC_G_F.PCER, LC.E1NEUTR),
                (OP.CALC_G_F.PELI, LC.E2NEUTR),
            ),
            para_out=((SP.PGTHETA, LC.CKGTHET),),
        ),
        OP.CALC_K_G(
            te=222,
            para_in=(
                (OP.CALC_K_G.PBASLOR, LC.N9NEUT_R),
                (OP.CALC_K_G.PCOMPOR, LC.CCOMPOR),
                (SP.PCOURB, LC.G27NEUTR),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PEPSINR, CEPSINO),
                (SP.PFRVOLU, NFORCER),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PPULPRO, LC.CFREQR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, ESIGMNO),
                (OP.CALC_K_G.PTHETAR, LC.ETHETA),
                (OP.CALC_K_G.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.CALC_K_G.PDEG, LC.E1NEUTI),
                (OP.CALC_K_G.PLAG, LC.CABSLAG),
                (OP.CALC_K_G.PCER, LC.E1NEUTR),
                (OP.CALC_K_G.PELI, LC.E2NEUTR),
            ),
            para_out=((SP.PGTHETA, LC.CKGTHET),),
        ),
        OP.CALC_K_G_F(
            te=222,
            para_in=(
                (OP.CALC_K_G_F.PBASLOR, LC.N9NEUT_R),
                (OP.CALC_K_G_F.PCOMPOR, LC.CCOMPOR),
                (SP.PCOURB, LC.G27NEUTR),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PEPSINF, CEPSINF),
                (SP.PFFVOLU, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PPULPRO, LC.CFREQR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, ESIGMNO),
                (SP.PINSTR, CTEMPSR),
                (OP.CALC_K_G_F.PTHETAR, LC.ETHETA),
                (OP.CALC_K_G_F.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.CALC_K_G_F.PDEG, LC.E1NEUTI),
                (OP.CALC_K_G_F.PLAG, LC.CABSLAG),
                (OP.CALC_K_G_F.PCER, LC.E1NEUTR),
                (OP.CALC_K_G_F.PELI, LC.E2NEUTR),
            ),
            para_out=((SP.PGTHETA, LC.CKGTHET),),
        ),
        OP.CALC_G_XFEM(
            te=27,
            para_in=(
                (SP.PACCELE, NDEPLAR),
                (OP.CALC_G_XFEM.PCOMPOR, LC.CCOMPOR),
                (SP.PCONTGR, ESIGMPG),
                (OP.CALC_G_XFEM.PCONTRR, ESIGMPG),
                (SP.PDEFOPL, EDEFONO),
                (SP.PDEPINR, NDEPLAR),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PEPSINR, CEPSINO),
                (SP.PFRVOLU, NFORCER),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, ESIGMNO),
                (SP.PTHETAR, NDEPLAR),
                (OP.CALC_G_XFEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.CALC_G_XFEM.PVARIPR, LC.ZVARINO),
                (SP.PVITESS, NDEPLAR),
            ),
            para_out=((SP.PGTHETA, LC.CGTHETA),),
        ),
        OP.CALC_G_XFEM_F(
            te=27,
            para_in=(
                (SP.PACCELE, NDEPLAR),
                (OP.CALC_G_XFEM_F.PCOMPOR, LC.CCOMPOR),
                (SP.PCONTGR, ESIGMPG),
                (OP.CALC_G_XFEM_F.PCONTRR, ESIGMPG),
                (SP.PDEFOPL, EDEFONO),
                (SP.PDEPINR, NDEPLAR),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PEPSINF, CEPSINF),
                (SP.PFFVOLU, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, ESIGMNO),
                (SP.PINSTR, CTEMPSR),
                (SP.PTHETAR, NDEPLAR),
                (OP.CALC_G_XFEM_F.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.CALC_G_XFEM_F.PVARIPR, LC.ZVARINO),
                (SP.PVITESS, NDEPLAR),
            ),
            para_out=((SP.PGTHETA, LC.CGTHETA),),
        ),
        OP.CALC_K_G_XFEM(
            te=295,
            para_in=(
                (OP.CALC_K_G_XFEM.PBASLOR, LC.N9NEUT_R),
                (OP.CALC_K_G_XFEM.PCOMPOR, LC.CCOMPOR),
                (SP.PCOURB, LC.G27NEUTR),
                (SP.PDEPINR, NDEPLAR),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PEPSINR, CEPSINO),
                (SP.PFRVOLU, NFORCER),
                (SP.PGEOMER, NGEOMER),
                (OP.CALC_K_G_XFEM.PLSN, LC.N1NEUT_R),
                (OP.CALC_K_G_XFEM.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PPULPRO, LC.CFREQR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, ESIGMNO),
                (SP.PTHETAR, NDEPLAR),
                (OP.CALC_K_G_XFEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PGTHETA, LC.CKGTX3D),),
        ),
        OP.CALC_K_G_XFEM_F(
            te=295,
            para_in=(
                (OP.CALC_K_G_XFEM_F.PBASLOR, LC.N9NEUT_R),
                (OP.CALC_K_G_XFEM_F.PCOMPOR, LC.CCOMPOR),
                (SP.PCOURB, LC.G27NEUTR),
                (SP.PDEPINR, NDEPLAR),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PEPSINF, CEPSINF),
                (SP.PFFVOLU, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (OP.CALC_K_G_XFEM_F.PLSN, LC.N1NEUT_R),
                (OP.CALC_K_G_XFEM_F.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PPULPRO, LC.CFREQR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, ESIGMNO),
                (SP.PINSTR, CTEMPSR),
                (SP.PTHETAR, NDEPLAR),
                (OP.CALC_K_G_XFEM_F.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PGTHETA, LC.CKGTX3D),),
        ),
        OP.CHAR_LIMITE(
            te=483,
            para_in=(
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_LIMITE.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PECHLI, LC.ECHALIM),),
        ),
        OP.CHAR_MECA_FF3D3D(
            te=17,
            para_in=((SP.PFF3D3D, CFORCEF), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTUR, VVECTUR),),
        ),
        OP.CHAR_MECA_FR3D3D(
            te=16,
            para_in=((SP.PFR3D3D, NFORCER), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTUR, VVECTUR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=15,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, VVECTUR),),
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
        OP.EPSL_ELGA(
            te=25,
            para_in=(
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPSL_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPSL_ELGA.PDEFOPG, EDEFOPG),),
        ),
        OP.EPSL_ELNO(
            te=4, para_in=((OP.EPSL_ELNO.PDEFOPG, EDEFOPG),), para_out=((SP.PDEFONO, EDEFONO),)
        ),
        OP.EPSI_ELGA(
            te=453,
            para_in=(
                (SP.PDEPLAR, NDEPLAR),
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
        OP.ERME_ELEM(
            te=375,
            para_in=(
                (SP.PCONTNO, ESIGMNO),
                (SP.PFFVOLU, CFORCEF),
                (SP.PFORCE, LC.CREFERI),
                (SP.PFRVOLU, EFORCER),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PPRESS, LC.CREFERI),
                (SP.PROTATR, LC.CROTATR),
                (SP.PINSTR, CTEMPSR),
                (OP.ERME_ELEM.PVOISIN, LC.EVOISIN),
            ),
            para_out=((OP.ERME_ELEM.PERREUR, EERREUR),),
        ),
        OP.ERME_ELNO(
            te=379, para_in=((OP.ERME_ELNO.PERREUR, EERREUR),), para_out=((SP.PERRENO, EERRENO),)
        ),
        OP.FORC_NODA(
            te=596,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, ECONTPG),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=595,
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
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.FULL_MECA_ELAS(
            te=595,
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
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA_ELAS.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.GRAD_NEUT9_R(
            te=398,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PNEUTER, LC.N9NEUT_R)),
            para_out=((OP.GRAD_NEUT9_R.PGNEUTR, LC.G27NEUTR),),
        ),
        OP.INIT_MAIL_VOIS(te=99, para_out=((OP.INIT_MAIL_VOIS.PVOISIN, LC.EVOISIN),)),
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
            para_out=((SP.PMATUUR, VMATUUR),),
        ),
        OP.NORME_L2(
            te=563,
            para_in=(
                (SP.PCALCI, LC.EMNEUT_I),
                (SP.PCHAMPG, EGNEUT_R),
                (SP.PCOEFR, EMNEUT_R),
                (OP.NORME_L2.PCOORPG, EGGEOP_R),
            ),
            para_out=((SP.PNORME, LC.ENORME),),
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
        OP.QIRE_ELEM(
            te=368,
            para_in=(
                (SP.PCONSTR, LC.CCONSTR),
                (SP.PCONTNOD, ESIGMNO),
                (SP.PCONTNOP, ESIGMNO),
                (SP.PFFVOLUD, CFORCEF),
                (SP.PFFVOLUP, CFORCEF),
                (SP.PFORCED, LC.CREFERI),
                (SP.PFORCEP, LC.CREFERI),
                (SP.PFRVOLUD, EFORCER),
                (SP.PFRVOLUP, EFORCER),
                (SP.PGEOMER, NGEOMER),
                (SP.PPESANRD, LC.CPESANR),
                (SP.PPESANRP, LC.CPESANR),
                (SP.PPRESSD, LC.CREFERI),
                (SP.PPRESSP, LC.CREFERI),
                (SP.PROTATRD, LC.CROTATR),
                (SP.PROTATRP, LC.CROTATR),
                (SP.PINSTR, CTEMPSR),
                (OP.QIRE_ELEM.PVOISIN, LC.EVOISIN),
            ),
            para_out=((OP.QIRE_ELEM.PERREUR, EERREUR),),
        ),
        OP.QIRE_ELNO(
            te=379, para_in=((OP.QIRE_ELNO.PERREUR, EERREUR),), para_out=((SP.PERRENO, EERRENO),)
        ),
        OP.RAPH_MECA(
            te=595,
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
        OP.REFE_FORC_NODA(
            te=598,
            para_in=(
                (OP.REFE_FORC_NODA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, NGEOMER),
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
            te=595,
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
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_TANG(
            te=595,
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
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PDEPL_R, LC.EGDEP3D),
                (OP.TOU_INI_ELGA.PGEOM_R, EGGEOM_R),
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
                (OP.TOU_INI_ELGA.PSIEF_R, ECONTPG),
                (OP.TOU_INI_ELGA.PSOUR_R, ESOURCR),
                (OP.TOU_INI_ELGA.PVARI_R, ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),
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
class MINCOS_TETRA4(MINCOS_HEXA8):
    """Mechanics - 3D - Incompressible - UPO model - TETRA4"""

    meshType = MT.TETRA4
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)),)
    elrefe = (
        ElrefeLoc(
            MT.TE4,
            gauss=("RIGI=FPG4", "MASS=FPG4", "NOEU=NOEU", "FPG1=FPG1"),
            mater=("RIGI", "MASS", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.TE4, gauss=("RIGI=FPG4",)),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG3", "MASS=FPG3", "NOEU=NOEU")),
    )


# ------------------------------------------------------------
class MINCOS_PENTA6(MINCOS_HEXA8):
    """Mechanics - 3D - Incompressible - UPO model - PENTA6"""

    meshType = MT.PENTA6
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),)
    elrefe = (
        ElrefeLoc(
            MT.PE6,
            gauss=("RIGI=FPG6", "MASS=FPG6", "NOEU=NOEU", "FPG1=FPG1"),
            mater=("RIGI", "MASS", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.PE6, gauss=("RIGI=FPG6",)),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG4", "NOEU=NOEU")),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG3", "MASS=FPG3", "NOEU=NOEU")),
    )


# ------------------------------------------------------------
class MINCOS_PYRAM5(MINCOS_HEXA8):
    """Mechanics - 3D - Incompressible - UPO model - PYRAM5"""

    meshType = MT.PYRAM5
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5)),)
    elrefe = (
        ElrefeLoc(
            MT.PY5,
            gauss=("RIGI=FPG5", "MASS=FPG5", "NOEU=NOEU", "FPG1=FPG1"),
            mater=("RIGI", "MASS", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.PY5, gauss=("RIGI=FPG5",)),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG4", "NOEU=NOEU")),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG3", "MASS=FPG3", "NOEU=NOEU")),
    )
