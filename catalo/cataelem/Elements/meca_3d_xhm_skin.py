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

# CATALOGUES DES ELEMENTS DE BORD 3D X-FEM MULTI HEAVISIDE


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
        ("EN1", ("DX", "DY", "DZ", "H1X", "H1Y", "H1Z")),
        ("EN2", ("DX", "DY", "DZ", "H1X", "H1Y", "H1Z", "H2X", "H2Y", "H2Z")),
        ("EN3", ("DX", "DY", "DZ", "H1X", "H1Y", "H1Z", "H2X", "H2Y", "H2Z", "H3X", "H3Y", "H3Z")),
        (
            "EN4",
            (
                "DX",
                "DY",
                "DZ",
                "H1X",
                "H1Y",
                "H1Z",
                "H2X",
                "H2Y",
                "H2Z",
                "H3X",
                "H3Y",
                "H3Z",
                "H4X",
                "H4Y",
                "H4Z",
            ),
        ),
    ),
)


NTHETAR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))


CFORCEF = LocatedComponents(phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ"))


NFORCER = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


STANO_I = LocatedComponents(phys=PHY.N120_I, type="ELNO", components=("X1",))


E48NEUTR = LocatedComponents(phys=PHY.N132_R, type="ELEM", components=("X[48]",))


E24NEUI = LocatedComponents(phys=PHY.N512_I, type="ELEM", components=("X[24]",))


CPRESSF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PRES",))


CPRES_R = LocatedComponents(phys=PHY.PRES_R, type="ELEM", components=("PRES",))


EPRESNO = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MECA_XH1_FACE3(Element):

    """Please document this element"""

    meshType = MT.TRIA3
    nodes = (SetOfNodes("EN1", (1, 2, 3)),)
    elrefe = (ElrefeLoc(MT.TR3, gauss=("RIGI=FPG4",)),)
    calculs = (
        OP.CALC_G_XFEM(
            te=580,
            para_in=((SP.PFR2D3D, NFORCER), (SP.PPRESSR, EPRESNO), (SP.PTHETAR, NTHETAR)),
            para_out=((SP.PGTHETA, LC.CGTHETA),),
        ),
        OP.CALC_G_XFEM_F(
            te=580,
            para_in=((SP.PFF2D3D, CFORCEF), (SP.PPRESSF, CPRESSF), (SP.PTHETAR, NTHETAR)),
            para_out=((SP.PGTHETA, LC.CGTHETA),),
        ),
        OP.CALC_K_G_XFEM(
            te=580,
            para_in=((SP.PFR2D3D, NFORCER), (SP.PPRESSR, EPRESNO), (SP.PTHETAR, NTHETAR)),
            para_out=((SP.PGTHETA, LC.CKGTX3D),),
        ),
        OP.CALC_K_G_XFEM_F(
            te=580,
            para_in=((SP.PFF2D3D, CFORCEF), (SP.PPRESSF, CPRESSF), (SP.PTHETAR, NTHETAR)),
            para_out=((SP.PGTHETA, LC.CKGTX3D),),
        ),
        OP.CHAR_MECA_FF2D3D(
            te=36,
            para_in=(
                (OP.CHAR_MECA_FF2D3D.PCNSETO, LC.E72NEUI),
                (SP.PFF2D3D, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_FF2D3D.PHEAVTO, E24NEUI),
                (OP.CHAR_MECA_FF2D3D.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_FF2D3D.PLONCHA, LC.E10NEUTI),
                (OP.CHAR_MECA_FF2D3D.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_FF2D3D.PLST, LC.N1NEUT_R),
                (OP.CHAR_MECA_FF2D3D.PPINTTO, E48NEUTR),
                (OP.CHAR_MECA_FF2D3D.PSTANO, STANO_I),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FR2D3D(
            te=36,
            para_in=(
                (OP.CHAR_MECA_FR2D3D.PCNSETO, LC.E72NEUI),
                (SP.PFR2D3D, NFORCER),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_FR2D3D.PHEAVTO, E24NEUI),
                (OP.CHAR_MECA_FR2D3D.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_FR2D3D.PLONCHA, LC.E10NEUTI),
                (OP.CHAR_MECA_FR2D3D.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_FR2D3D.PLST, LC.N1NEUT_R),
                (OP.CHAR_MECA_FR2D3D.PPINTTO, E48NEUTR),
                (OP.CHAR_MECA_FR2D3D.PSTANO, STANO_I),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_F(
            te=36,
            para_in=(
                (OP.CHAR_MECA_PRES_F.PCNSETO, LC.E72NEUI),
                (OP.CHAR_MECA_PRES_F.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_PRES_F.PHEAVTO, E24NEUI),
                (OP.CHAR_MECA_PRES_F.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_PRES_F.PHEA_SE, E24NEUI),
                (OP.CHAR_MECA_PRES_F.PLONCHA, LC.E10NEUTI),
                (OP.CHAR_MECA_PRES_F.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_PRES_F.PLST, LC.N1NEUT_R),
                (OP.CHAR_MECA_PRES_F.PPINTTO, E48NEUTR),
                (SP.PPRESSF, CPRESSF),
                (OP.CHAR_MECA_PRES_F.PSTANO, STANO_I),
                (SP.PINSTR, CTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=36,
            para_in=(
                (OP.CHAR_MECA_PRES_R.PCNSETO, LC.E72NEUI),
                (OP.CHAR_MECA_PRES_R.PFISNO, LC.FISNO_I),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_PRES_R.PHEAVTO, E24NEUI),
                (OP.CHAR_MECA_PRES_R.PHEA_NO, LC.N5NEUTI),
                (OP.CHAR_MECA_PRES_R.PHEA_SE, E24NEUI),
                (OP.CHAR_MECA_PRES_R.PLONCHA, LC.E10NEUTI),
                (OP.CHAR_MECA_PRES_R.PLSN, LC.N1NEUT_R),
                (OP.CHAR_MECA_PRES_R.PLST, LC.N1NEUT_R),
                (OP.CHAR_MECA_PRES_R.PPINTTO, E48NEUTR),
                (SP.PPRESSR, EPRESNO),
                (OP.CHAR_MECA_PRES_R.PSTANO, STANO_I),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.INI_XFEM_ELNO(
            te=99,
            para_out=(
                (OP.INI_XFEM_ELNO.PFISNO, LC.FISNO_I),
                (OP.INI_XFEM_ELNO.PLSN, LC.N1NEUT_R),
                (OP.INI_XFEM_ELNO.PLST, LC.N1NEUT_R),
                (OP.INI_XFEM_ELNO.PSTANO, STANO_I),
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
                (OP.TOPOSE.PCNSETO, LC.E72NEUI),
                (OP.TOPOSE.PHEAVTO, E24NEUI),
                (OP.TOPOSE.PLONCHA, LC.E10NEUTI),
                (OP.TOPOSE.PPINTTO, E48NEUTR),
            ),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PPRES_R, CPRES_R),)),
        OP.TOU_INI_ELGA(te=99, para_out=((OP.TOU_INI_ELGA.PGEOM_R, EGGEOP_R),)),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PPRES_R, EPRESNO),
            ),
        ),
    )


# ------------------------------------------------------------
class MECA_XH2_FACE3(MECA_XH1_FACE3):

    """Please document this element"""

    meshType = MT.TRIA3
    nodes = (SetOfNodes("EN2", (1, 2, 3)),)
    elrefe = (ElrefeLoc(MT.TR3, gauss=("RIGI=FPG4",)),)


# ------------------------------------------------------------
class MECA_XH3_FACE3(MECA_XH1_FACE3):

    """Please document this element"""

    meshType = MT.TRIA3
    nodes = (SetOfNodes("EN3", (1, 2, 3)),)
    elrefe = (ElrefeLoc(MT.TR3, gauss=("RIGI=FPG4",)),)


# ------------------------------------------------------------
class MECA_XH4_FACE3(MECA_XH1_FACE3):

    """Please document this element"""

    meshType = MT.TRIA3
    nodes = (SetOfNodes("EN4", (1, 2, 3)),)
    elrefe = (ElrefeLoc(MT.TR3, gauss=("RIGI=FPG4",)),)


# ------------------------------------------------------------
class MECA_XH1_FACE4(MECA_XH1_FACE3):

    """Please document this element"""

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)),)
    elrefe = (ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4",)), ElrefeLoc(MT.TR3, gauss=("RIGI=FPG4",)))


# ------------------------------------------------------------
class MECA_XH2_FACE4(MECA_XH1_FACE3):

    """Please document this element"""

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN2", (1, 2, 3, 4)),)
    elrefe = (ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4",)), ElrefeLoc(MT.TR3, gauss=("RIGI=FPG4",)))


# ------------------------------------------------------------
class MECA_XH3_FACE4(MECA_XH1_FACE3):

    """Please document this element"""

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN3", (1, 2, 3, 4)),)
    elrefe = (ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4",)), ElrefeLoc(MT.TR3, gauss=("RIGI=FPG4",)))


# ------------------------------------------------------------
class MECA_XH4_FACE4(MECA_XH1_FACE3):

    """Please document this element"""

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN4", (1, 2, 3, 4)),)
    elrefe = (ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4",)), ElrefeLoc(MT.TR3, gauss=("RIGI=FPG4",)))
