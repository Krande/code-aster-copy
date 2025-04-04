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

# CATALOGUES DES ELEMENTS 2D_DP X-FEM MULTI HEAVISIDE SANS CONTACT


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
        ("EN1", ("DX", "DY", "H1X", "H1Y")),
        ("EN2", ("DX", "DY", "H1X", "H1Y", "H2X", "H2Y")),
        ("EN3", ("DX", "DY", "H1X", "H1Y", "H2X", "H2Y", "H3X", "H3Y")),
        ("EN4", ("DX", "DY", "H1X", "H1Y", "H2X", "H2Y", "H3X", "H3Y", "H4X", "H4Y")),
    ),
)


EDEPLPG = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELGA",
    location="XFEM",
    components=("DX", "DY", "H1X", "H1Y", "H2X", "H2Y", "H3X", "H3Y", "H4X", "H4Y"),
)


DDL_MECC = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY"))


EENERR = LocatedComponents(phys=PHY.ENER_R, type="ELEM", components=("TOTALE",))


CFORCEF = LocatedComponents(phys=PHY.FORC_F, type="ELEM", components=("FX", "FY"))


NFORCER = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="XFEM", components=("X", "Y", "W")
)


EGGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELGA", location="XFEM", components=("X", "Y"))


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


XFGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELGA", location="XFEM", components=("X", "Y"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


STANO_I = LocatedComponents(phys=PHY.N120_I, type="ELNO", components=("X1",))


E24NEUI = LocatedComponents(phys=PHY.N512_I, type="ELEM", components=("X[24]",))


E14NEUTI = LocatedComponents(phys=PHY.N720_I, type="ELEM", components=("X[14]",))


E16NEUTR = LocatedComponents(phys=PHY.N816_R, type="ELEM", components=("X[16]",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="XFEM", components=("X[30]",))


E1NEUTK = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z1",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="XFEM", components=("X[30]",))


EMNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[30]",))


CPRESSF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PRES", "CISA"))


EPRESNO = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES", "CISA"))


ECONTPC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELGA", location="XFEM", components=("SIXX", "SIYY", "SIZZ", "SIXY")
)


ECONTPG = LocatedComponents(
    phys=PHY.SIEF_R, type="ELGA", location="XFEM", components=("SIXX", "SIYY", "SIZZ", "SIXY")
)


ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="XFEM", components=("VARI",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MEDPTR3_XH1(Element):
    """Please document this element"""

    meshType = MT.TRIA3
    nodes = (SetOfNodes("EN1", (1, 2, 3)),)
    elrefe = (
        ElrefeLoc(
            MT.TR3,
            gauss=(
                "RIGI=FPG3",
                "XINT=FPG4",
                "NOEU_S=NOEU_S",
                "NOEU=NOEU",
                "XFEM=XFEM36",
                "FPG1=FPG1",
            ),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2", "MASS=FPG3")),
    )
    calculs = (
        OP.CALC_G_XFEM(
            te=288,
            para_in=(
                (OP.CALC_G_XFEM.PAINTER, LC.E40NEUTR),
                (OP.CALC_G_XFEM.PBASECO, LC.E32NEUTR),
                (OP.CALC_G_XFEM.PBASLOR, LC.N6NEUT_R),
                (OP.CALC_G_XFEM.PCFACE, E14NEUTI),
                (OP.CALC_G_XFEM.PCNSETO, LC.E72NEUI),
                (OP.CALC_G_XFEM.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PFRVOLU, NFORCER),
                (SP.PGEOMER, NGEOMER),
                (OP.CALC_G_XFEM.PHEAVTO, E24NEUI),
                (OP.CALC_G_XFEM.PHEA_NO, LC.N5NEUTI),
                (OP.CALC_G_XFEM.PLONCHA, LC.E10NEUTI),
                (OP.CALC_G_XFEM.PLONGCO, LC.E3NEUTI),
                (OP.CALC_G_XFEM.PLSN, LC.N1NEUT_R),
                (OP.CALC_G_XFEM.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CALC_G_XFEM.PPINTER, E16NEUTR),
                (OP.CALC_G_XFEM.PPINTTO, LC.E24NEUTR),
                (OP.CALC_G_XFEM.PPMILTO, LC.E22NEUTR),
                (SP.PPRESSR, EPRESNO),
                (SP.PROTATR, LC.CROTATR),
                (SP.PTHETAR, DDL_MECC),
                (OP.CALC_G_XFEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PGTHETA, LC.CGTHETA),),
        ),
        OP.CALC_G_XFEM_F(
            te=288,
            para_in=(
                (OP.CALC_G_XFEM_F.PAINTER, LC.E40NEUTR),
                (OP.CALC_G_XFEM_F.PBASECO, LC.E32NEUTR),
                (OP.CALC_G_XFEM_F.PBASLOR, LC.N6NEUT_R),
                (OP.CALC_G_XFEM_F.PCFACE, E14NEUTI),
                (OP.CALC_G_XFEM_F.PCNSETO, LC.E72NEUI),
                (OP.CALC_G_XFEM_F.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PFFVOLU, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (OP.CALC_G_XFEM_F.PHEAVTO, E24NEUI),
                (OP.CALC_G_XFEM_F.PHEA_NO, LC.N5NEUTI),
                (OP.CALC_G_XFEM_F.PLONCHA, LC.E10NEUTI),
                (OP.CALC_G_XFEM_F.PLONGCO, LC.E3NEUTI),
                (OP.CALC_G_XFEM_F.PLSN, LC.N1NEUT_R),
                (OP.CALC_G_XFEM_F.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CALC_G_XFEM_F.PPINTER, E16NEUTR),
                (OP.CALC_G_XFEM_F.PPINTTO, LC.E24NEUTR),
                (OP.CALC_G_XFEM_F.PPMILTO, LC.E22NEUTR),
                (SP.PPRESSF, CPRESSF),
                (SP.PROTATR, LC.CROTATR),
                (SP.PINSTR, CTEMPSR),
                (SP.PTHETAR, DDL_MECC),
                (OP.CALC_G_XFEM_F.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PGTHETA, LC.CGTHETA),),
        ),
        OP.CALC_K_G_XFEM(
            te=297,
            para_in=(
                (OP.CALC_K_G_XFEM.PAINTER, LC.E40NEUTR),
                (OP.CALC_K_G_XFEM.PBASECO, LC.E32NEUTR),
                (OP.CALC_K_G_XFEM.PBASLOR, LC.N6NEUT_R),
                (OP.CALC_K_G_XFEM.PCFACE, E14NEUTI),
                (OP.CALC_K_G_XFEM.PCNSETO, LC.E72NEUI),
                (OP.CALC_K_G_XFEM.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PFRVOLU, NFORCER),
                (SP.PGEOMER, NGEOMER),
                (OP.CALC_K_G_XFEM.PHEAVTO, E24NEUI),
                (OP.CALC_K_G_XFEM.PHEA_NO, LC.N5NEUTI),
                (OP.CALC_K_G_XFEM.PLONCHA, LC.E10NEUTI),
                (OP.CALC_K_G_XFEM.PLONGCO, LC.E3NEUTI),
                (OP.CALC_K_G_XFEM.PLSN, LC.N1NEUT_R),
                (OP.CALC_K_G_XFEM.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CALC_K_G_XFEM.PPINTER, E16NEUTR),
                (OP.CALC_K_G_XFEM.PPINTTO, LC.E24NEUTR),
                (OP.CALC_K_G_XFEM.PPMILTO, LC.E22NEUTR),
                (SP.PPRESSR, EPRESNO),
                (SP.PPULPRO, LC.CFREQR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PTHETAR, DDL_MECC),
                (OP.CALC_K_G_XFEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PGTHETA, LC.CKGTX2D),),
        ),
        OP.CALC_K_G_XFEM_F(
            te=297,
            para_in=(
                (OP.CALC_K_G_XFEM_F.PAINTER, LC.E40NEUTR),
                (OP.CALC_K_G_XFEM_F.PBASECO, LC.E32NEUTR),
                (OP.CALC_K_G_XFEM_F.PBASLOR, LC.N6NEUT_R),
                (OP.CALC_K_G_XFEM_F.PCFACE, E14NEUTI),
                (OP.CALC_K_G_XFEM_F.PCNSETO, LC.E72NEUI),
                (OP.CALC_K_G_XFEM_F.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PFFVOLU, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (OP.CALC_K_G_XFEM_F.PHEAVTO, E24NEUI),
                (OP.CALC_K_G_XFEM_F.PHEA_NO, LC.N5NEUTI),
                (OP.CALC_K_G_XFEM_F.PLONCHA, LC.E10NEUTI),
                (OP.CALC_K_G_XFEM_F.PLONGCO, LC.E3NEUTI),
                (OP.CALC_K_G_XFEM_F.PLSN, LC.N1NEUT_R),
                (OP.CALC_K_G_XFEM_F.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CALC_K_G_XFEM_F.PPINTER, E16NEUTR),
                (OP.CALC_K_G_XFEM_F.PPINTTO, LC.E24NEUTR),
                (OP.CALC_K_G_XFEM_F.PPMILTO, LC.E22NEUTR),
                (SP.PPRESSF, CPRESSF),
                (SP.PPULPRO, LC.CFREQR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PINSTR, CTEMPSR),
                (SP.PTHETAR, DDL_MECC),
                (OP.CALC_K_G_XFEM_F.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PGTHETA, LC.CKGTX2D),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=441,
            para_in=(
                (OP.CHAR_MECA_PESA_R.PCNSETO, LC.E72NEUI),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_PESA_R.PHEAVTO, E24NEUI),
                (OP.CHAR_MECA_PESA_R.PLONCHA, LC.E10NEUTI),
                (OP.CHAR_MECA_PESA_R.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_PESA_R.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PPINTTO, LC.E24NEUTR),
                (OP.CHAR_MECA_PESA_R.PPMILTO, LC.E22NEUTR),
                (OP.CHAR_MECA_PESA_R.PSTANO, STANO_I),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        #       -- te0580 : ne resout que les cas triviaux : 0.
        OP.CHAR_MECA_PRES_F(
            te=580, para_in=((SP.PPRESSF, CPRESSF),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.CHAR_MECA_PRES_R(
            te=580, para_in=((SP.PPRESSR, EPRESNO),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.CHAR_MECA_ROTA_R(
            te=441,
            para_in=(
                (OP.CHAR_MECA_ROTA_R.PCNSETO, LC.E72NEUI),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_ROTA_R.PHEAVTO, E24NEUI),
                (OP.CHAR_MECA_ROTA_R.PLONCHA, LC.E10NEUTI),
                (OP.CHAR_MECA_ROTA_R.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_ROTA_R.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_ROTA_R.PPINTTO, LC.E24NEUTR),
                (OP.CHAR_MECA_ROTA_R.PPMILTO, LC.E22NEUTR),
                (SP.PROTATR, LC.CROTATR),
                (OP.CHAR_MECA_ROTA_R.PSTANO, STANO_I),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=541,
            para_in=(
                (OP.CHAR_MECA_TEMP_R.PBASLOR, LC.N6NEUT_R),
                (SP.PCAMASS, LC.CCAMA2D),
                (OP.CHAR_MECA_TEMP_R.PCNSETO, LC.E72NEUI),
                (OP.CHAR_MECA_TEMP_R.PCOMPOR, LC.CCOMPOR),
                (OP.CHAR_MECA_TEMP_R.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_TEMP_R.PHEAVTO, E24NEUI),
                (OP.CHAR_MECA_TEMP_R.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_TEMP_R.PLONCHA, LC.E10NEUTI),
                (OP.CHAR_MECA_TEMP_R.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_TEMP_R.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_TEMP_R.PPINTTO, LC.E24NEUTR),
                (OP.CHAR_MECA_TEMP_R.PPMILTO, LC.E22NEUTR),
                (OP.CHAR_MECA_TEMP_R.PSTANO, STANO_I),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PCONTRT, ECONTPG), (SP.PVECTUR, MVECTUR)),
        ),
        OP.CHVOIS_XFEM(
            te=400,
            para_in=((OP.CHVOIS_XFEM.PCNSETO, LC.E72NEUI), (OP.CHVOIS_XFEM.PLONCHA, LC.E10NEUTI)),
            para_out=((OP.CHVOIS_XFEM.PCVOISX, LC.E18NEUI),),
        ),
        OP.COOR_ELGA(
            te=481,
            para_in=(
                (OP.COOR_ELGA.PCNSETO, LC.E72NEUI),
                (SP.PGEOMER, NGEOMER),
                (OP.COOR_ELGA.PLONCHA, LC.E10NEUTI),
                (OP.COOR_ELGA.PPINTTO, LC.E24NEUTR),
                (OP.COOR_ELGA.PPMILTO, LC.E22NEUTR),
            ),
            para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),),
        ),
        OP.DEPL_XPG(
            te=566,
            para_in=(
                (OP.DEPL_XPG.PBASLOR, LC.N6NEUT_R),
                (SP.PDEPLNO, DDL_MECA),
                (OP.DEPL_XPG.PHEAVTO, E24NEUI),
                (OP.DEPL_XPG.PHEA_NO, LC.N5NEUTI),
                (OP.DEPL_XPG.PLONCHA, LC.E10NEUTI),
                (OP.DEPL_XPG.PLSN, LC.N1NEUT_R),
                (OP.DEPL_XPG.PLST, LC.N1NEUT_R),
                (OP.DEPL_XPG.PXFGEOM, XFGEOM_R),
            ),
            para_out=((SP.PDEPLPG, EDEPLPG),),
        ),
        OP.ENEL_ELEM(
            te=565,
            para_in=(
                (OP.ENEL_ELEM.PCNSETO, LC.E72NEUI),
                (OP.ENEL_ELEM.PCOMPOR, LC.CCOMPOR),
                (OP.ENEL_ELEM.PCONTPR, ECONTPG),
                (SP.PDEPLR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.ENEL_ELEM.PLONCHA, LC.E10NEUTI),
                (SP.PMATERC, LC.CMATERC),
                (OP.ENEL_ELEM.PPINTTO, LC.E24NEUTR),
                (OP.ENEL_ELEM.PPMILTO, LC.E22NEUTR),
                (OP.ENEL_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.ENEL_ELEM.PVARIPR, ZVARIPG),
            ),
            para_out=((SP.PENERD1, EENERR),),
        ),
        OP.FORC_NODA(
            te=542,
            para_in=(
                (OP.FORC_NODA.PBASLOR, LC.N6NEUT_R),
                (OP.FORC_NODA.PCNSETO, LC.E72NEUI),
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
                (OP.FORC_NODA.PPINTTO, LC.E24NEUTR),
                (OP.FORC_NODA.PPMILTO, LC.E22NEUTR),
                (OP.FORC_NODA.PSTANO, STANO_I),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=539,
            para_in=(
                (OP.FULL_MECA.PBASLOR, LC.N6NEUT_R),
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCNSETO, LC.E72NEUI),
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
                (OP.FULL_MECA.PPMILTO, LC.E22NEUTR),
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
        OP.RAPH_MECA(
            te=539,
            para_in=(
                (OP.RAPH_MECA.PBASLOR, LC.N6NEUT_R),
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RAPH_MECA.PCNSETO, LC.E72NEUI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.RAPH_MECA.PHEAVTO, E24NEUI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.RAPH_MECA.PLONCHA, LC.E10NEUTI),
                (OP.RAPH_MECA.PLSN, LC.N1NEUT_R),
                (OP.RAPH_MECA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RAPH_MECA.PPINTTO, LC.E24NEUTR),
                (OP.RAPH_MECA.PPMILTO, LC.E22NEUTR),
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
        OP.RIGI_MECA(
            te=539,
            para_in=(
                (OP.RIGI_MECA.PBASLOR, LC.N6NEUT_R),
                (OP.RIGI_MECA.PCNSETO, LC.E72NEUI),
                (OP.RIGI_MECA.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.RIGI_MECA.PHEAVTO, E24NEUI),
                (OP.RIGI_MECA.PHEA_NO, LC.N5NEUTI),
                (OP.RIGI_MECA.PLONCHA, LC.E10NEUTI),
                (OP.RIGI_MECA.PLSN, LC.N1NEUT_R),
                (OP.RIGI_MECA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PPINTTO, LC.E24NEUTR),
                (OP.RIGI_MECA.PPMILTO, LC.E22NEUTR),
                (OP.RIGI_MECA.PSTANO, STANO_I),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_TANG(
            te=539,
            para_in=(
                (OP.RIGI_MECA_TANG.PBASLOR, LC.N6NEUT_R),
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCNSETO, LC.E72NEUI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (OP.RIGI_MECA_TANG.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.RIGI_MECA_TANG.PHEAVTO, E24NEUI),
                (OP.RIGI_MECA_TANG.PHEA_NO, LC.N5NEUTI),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (OP.RIGI_MECA_TANG.PLONCHA, LC.E10NEUTI),
                (OP.RIGI_MECA_TANG.PLSN, LC.N1NEUT_R),
                (OP.RIGI_MECA_TANG.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA_TANG.PPINTTO, LC.E24NEUTR),
                (OP.RIGI_MECA_TANG.PPMILTO, LC.E22NEUTR),
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
                (OP.SIEF_ELGA.PCNSETO, LC.E72NEUI),
                (OP.SIEF_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (OP.SIEF_ELGA.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.SIEF_ELGA.PHEAVTO, E24NEUI),
                (OP.SIEF_ELGA.PHEA_NO, LC.N5NEUTI),
                (OP.SIEF_ELGA.PLONCHA, LC.E10NEUTI),
                (OP.SIEF_ELGA.PLSN, LC.N1NEUT_R),
                (OP.SIEF_ELGA.PLST, LC.N1NEUT_R),
                (SP.PMATERC, LC.CMATERC),
                (OP.SIEF_ELGA.PPINTTO, LC.E24NEUTR),
                (OP.SIEF_ELGA.PPMILTO, LC.E22NEUTR),
                (OP.SIEF_ELGA.PSTANO, STANO_I),
                (OP.SIEF_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.SIEF_ELGA.PCONTRR, ECONTPG),),
        ),
        OP.SIGM_ELGA(
            te=546,
            para_in=((SP.PSIEFR, ECONTPG),),
            para_out=((SP.PSIGMC, ECONTPC), (SP.PSIGMR, ECONTPG)),
        ),
        OP.TOPOFA(
            te=510,
            para_in=(
                (OP.TOPOFA.PAINTTO, LC.E60NEUTR),
                (OP.TOPOFA.PCNSETO, LC.E72NEUI),
                (SP.PDECOU, E1NEUTK),
                (SP.PGEOMER, NGEOMER),
                (SP.PGRADLN, LC.N2NEUT_R),
                (SP.PGRADLT, LC.N2NEUT_R),
                (OP.TOPOFA.PHEAVTO, E24NEUI),
                (OP.TOPOFA.PLONCHA, LC.E10NEUTI),
                (OP.TOPOFA.PLSN, LC.N1NEUT_R),
                (OP.TOPOFA.PLST, LC.N1NEUT_R),
                (OP.TOPOFA.PPINTTO, LC.E24NEUTR),
                (OP.TOPOFA.PPMILTO, LC.E22NEUTR),
            ),
            para_out=(
                (OP.TOPOFA.PAINTER, LC.E40NEUTR),
                (OP.TOPOFA.PBASECO, LC.E32NEUTR),
                (OP.TOPOFA.PCFACE, E14NEUTI),
                (SP.PGESCLA, E16NEUTR),
                (OP.TOPOFA.PLONGCO, LC.E3NEUTI),
                (OP.TOPOFA.PPINTER, E16NEUTR),
            ),
        ),
        OP.TOPONO(
            te=120,
            para_in=(
                (OP.TOPONO.PCNSETO, LC.E72NEUI),
                (SP.PFISCO, LC.FISCO_I),
                (OP.TOPONO.PFISNO, LC.FISNO_I),
                (OP.TOPONO.PHEAVTO, E24NEUI),
                (SP.PLEVSET, LC.N1NEUT_R),
                (OP.TOPONO.PLONCHA, LC.E10NEUTI),
            ),
            para_out=((OP.TOPONO.PHEA_NO, LC.N5NEUTI), (OP.TOPONO.PHEA_SE, E24NEUI)),
        ),
        OP.TOPOSE(
            te=514,
            para_in=((SP.PFISCO, LC.FISCO_I), (SP.PGEOMER, NGEOMER), (SP.PLEVSET, LC.N1NEUT_R)),
            para_out=(
                (OP.TOPOSE.PAINTTO, LC.E60NEUTR),
                (OP.TOPOSE.PCNSETO, LC.E72NEUI),
                (OP.TOPOSE.PHEAVTO, E24NEUI),
                (OP.TOPOSE.PLONCHA, LC.E10NEUTI),
                (OP.TOPOSE.PPINTTO, LC.E24NEUTR),
                (OP.TOPOSE.PPMILTO, LC.E22NEUTR),
            ),
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
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
        OP.XFEM_XPG(
            te=46,
            para_in=(
                (OP.XFEM_XPG.PCNSETO, LC.E72NEUI),
                (SP.PGEOMER, NGEOMER),
                (OP.XFEM_XPG.PHEAVTO, E24NEUI),
                (OP.XFEM_XPG.PLONCHA, LC.E10NEUTI),
                (OP.XFEM_XPG.PPINTTO, LC.E24NEUTR),
                (OP.XFEM_XPG.PPMILTO, LC.E22NEUTR),
            ),
            para_out=((OP.XFEM_XPG.PXFGEOM, XFGEOM_R),),
        ),
    )


# ------------------------------------------------------------
class MEDPTR3_XH2(MEDPTR3_XH1):
    """Please document this element"""

    meshType = MT.TRIA3
    nodes = (SetOfNodes("EN2", (1, 2, 3)),)
    elrefe = (
        ElrefeLoc(
            MT.TR3,
            gauss=(
                "RIGI=FPG3",
                "XINT=FPG4",
                "NOEU_S=NOEU_S",
                "NOEU=NOEU",
                "XFEM=XFEM72",
                "FPG1=FPG1",
            ),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2", "MASS=FPG3")),
    )


# ------------------------------------------------------------
class MEDPTR3_XH3(MEDPTR3_XH1):
    """Please document this element"""

    meshType = MT.TRIA3
    nodes = (SetOfNodes("EN3", (1, 2, 3)),)
    elrefe = (
        ElrefeLoc(
            MT.TR3,
            gauss=(
                "RIGI=FPG3",
                "XINT=FPG4",
                "NOEU_S=NOEU_S",
                "NOEU=NOEU",
                "XFEM=XFEM72",
                "FPG1=FPG1",
            ),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2", "MASS=FPG3")),
    )


# ------------------------------------------------------------
class MEDPTR3_XH4(MEDPTR3_XH1):
    """Please document this element"""

    meshType = MT.TRIA3
    nodes = (SetOfNodes("EN4", (1, 2, 3)),)
    elrefe = (
        ElrefeLoc(
            MT.TR3,
            gauss=(
                "RIGI=FPG3",
                "XINT=FPG4",
                "NOEU_S=NOEU_S",
                "NOEU=NOEU",
                "XFEM=XFEM72",
                "FPG1=FPG1",
            ),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2", "MASS=FPG3")),
    )


# ------------------------------------------------------------
class MEDPQU4_XH1(MEDPTR3_XH1):
    """Please document this element"""

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)),)
    elrefe = (
        ElrefeLoc(
            MT.QU4,
            gauss=("RIGI=FPG4", "NOEU_S=NOEU_S", "NOEU=NOEU", "XFEM=XFEM72", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG3", "XINT=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2", "MASS=FPG3")),
    )


# ------------------------------------------------------------
class MEDPQU4_XH2(MEDPTR3_XH1):
    """Please document this element"""

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN2", (1, 2, 3, 4)),)
    elrefe = (
        ElrefeLoc(
            MT.QU4,
            gauss=("RIGI=FPG4", "NOEU_S=NOEU_S", "NOEU=NOEU", "XFEM=XFEM144", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG3", "XINT=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2", "MASS=FPG3")),
    )


# ------------------------------------------------------------
class MEDPQU4_XH3(MEDPTR3_XH1):
    """Please document this element"""

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN3", (1, 2, 3, 4)),)
    elrefe = (
        ElrefeLoc(
            MT.QU4,
            gauss=("RIGI=FPG4", "NOEU_S=NOEU_S", "NOEU=NOEU", "XFEM=XFEM144", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG3", "XINT=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2", "MASS=FPG3")),
    )


# ------------------------------------------------------------
class MEDPQU4_XH4(MEDPTR3_XH1):
    """Please document this element"""

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN4", (1, 2, 3, 4)),)
    elrefe = (
        ElrefeLoc(
            MT.QU4,
            gauss=("RIGI=FPG4", "NOEU_S=NOEU_S", "NOEU=NOEU", "XFEM=XFEM144", "FPG1=FPG1"),
            mater=("XFEM",),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG3", "XINT=FPG4")),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2", "MASS=FPG3")),
    )
