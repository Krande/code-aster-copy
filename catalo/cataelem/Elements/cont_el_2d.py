# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

# person_in_charge: ayaovi-dzifa.kudawoo at edf.fr

from cataelem.Tools.base_objects import LocatedComponents, ArrayOfComponents, SetOfNodes, ElrefeLoc
from cataelem.Tools.base_objects import Calcul, Element
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.located_components as LC
import cataelem.Commons.parameters as SP
import cataelem.Commons.mesh_types as MT
from cataelem.Options.options import OP


# ELEMENTARY TREATMENT OF 2D FRICTIONLESS ELEMENT WITH DEFI_CONTACT OPERATOR

# ----------------
# Modes locaux :
# ----------------


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("DX", "DY", "LAGS_C")), ("EN2", ("DX", "DY"))),
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class COS2S2(Element):
    """
    THE COS2S2 CLASS ELEMENT : SEG2/SEG2
    DEFI_CONTACT / CONTINUE / NODE-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2 (DX,DY,LAGS_C)
        SEG2 MASTER ELEMENT : 3-4 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.SEG22
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4)))
    calculs = (
        OP.CHAR_MECA_CONT(
            te=365,
            para_in=(
                (SP.PACCE_M, DDL_MECA),
                (SP.PCONFR, LC.CCONFR),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PVITE_M, DDL_MECA),
                (SP.PVITE_P, DDL_MECA),
            ),
            para_out=((SP.PVECTCR, MVECTUR), (SP.PVECTFR, MVECTUR)),
        ),
        OP.RIGI_CONT(
            te=364,
            para_in=(
                (SP.PACCE_M, DDL_MECA),
                (SP.PCONFR, LC.CCONFR),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PVITE_M, DDL_MECA),
                (SP.PVITE_P, DDL_MECA),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),)),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
    )


# ------------------------------------------------------------
class COS3S3(COS2S2):
    """
    COS3S3 DERIVED FROM THE COS2S2 CLASS ELEMENT  : SEG3/SEG3
    DEFI_CONTACT / CONTINUE / NODE-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 4-5-6 (DX,DY,LAGS_C)
        SEG3 MASTER ELEMENT : 1-2-3 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.SEG33
    nodes = (SetOfNodes("EN2", (4, 5, 6)), SetOfNodes("EN1", (1, 2, 3)))


# ------------------------------------------------------------
class COS2S3(COS2S2):
    """
    COS2S3 DERIVED FROM THE COS2S2 CLASS ELEMENT : SEG2/SEG3
    DEFI_CONTACT / CONTINUE / NODE-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C)
        SEG3 MASTER ELEMENT : 3-4-5 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.SEG23
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4, 5)))


# ------------------------------------------------------------
class COS3S2(COS2S2):
    """
    COS3S2 DERIVED FROM THE COS2S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / CONTINUE / NODE-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2-3   (DX,DY,LAGS_C)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.SEG32
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5)))


# ------------------------------------------------------------
class COS2S2A(COS2S2):
    """
    COS2S2A DERIVED FROM THE COS2S2 CLASS ELEMENT : SEG2/SEG2 (AXIS)
    DEFI_CONTACT / CONTINUE / NODE-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2  (DX,DY,LAGS_C)
        SEG2 MASTER ELEMENT : 3-4  (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.SEG22
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4)))


# ------------------------------------------------------------
class COS3S3A(COS2S2):
    """
    COS3S3 DERIVED FROM THE COS2S2 CLASS ELEMENT  : SEG3/SEG3 (AXIS)
    DEFI_CONTACT / CONTINUE / NODE-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 4-5-6 (DX,DY,LAGS_C)
        SEG3 MASTER ELEMENT : 1-2-3 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.SEG33
    nodes = (SetOfNodes("EN2", (4, 5, 6)), SetOfNodes("EN1", (1, 2, 3)))


# ------------------------------------------------------------
class COS2S3A(COS2S2):
    """
    COS2S3 DERIVED FROM THE COS2S2 CLASS ELEMENT : SEG2/SEG3 (AXIS)
    DEFI_CONTACT / CONTINUE / NODE-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C)
        SEG3 MASTER ELEMENT : 3-4-5 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.SEG23
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4, 5)))


# ------------------------------------------------------------
class COS3S2A(COS2S2):
    """
    COS3S2 DERIVED FROM THE COS2S2 CLASS ELEMENT : SEG3/SEG2 (AXIS)
    DEFI_CONTACT / CONTINUE / NODE-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2-3   (DX,DY,LAGS_C)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.SEG32
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5)))
