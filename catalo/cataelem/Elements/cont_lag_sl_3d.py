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


from cataelem.Tools.base_objects import LocatedComponents, ArrayOfComponents, SetOfNodes, ElrefeLoc
from cataelem.Tools.base_objects import Calcul, Element
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.located_components as LC
import cataelem.Commons.parameters as SP
import cataelem.Commons.mesh_types as MT
from cataelem.Options.options import OP


# ELEMENTARY TREATMENT OF 3D FRICTIONLESS ELEMENT WITH DEFI_CONTACT OPERATOR
# MORTAR LAGRANGIAN AUGMENTED METHOD

# ----------------
# Modes locaux :
# ----------------

DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("DX", "DY", "DZ", "LAGS_C")), ("EN2", ("DX", "DY", "DZ"))),
)

ECCONT = LocatedComponents(
    phys=PHY.CONT_R, type="ELNO", diff=True, components=(("EN1", ("COEF_C",)), ("EN2", ()))
)

ECFROT = LocatedComponents(
    phys=PHY.CONT_R, type="ELNO", diff=True, components=(("EN1", ("COEF_F",)), ("EN2", ()))
)

# ------------------------------------------------------------


class CMT33D(Element):
    """
    THE CMT33D CLASS ELEMENT :
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D : elementary treatments
    Local Numerotation :

    Input parameters :

    Output parameters :
    """

    meshType = MT.TRIA3
    nodes = (SetOfNodes("EN1", (1, 2, 3)),)
    calculs = (
        OP.EXISTE_DDL(
            te=99,
            para_out=(
                (OP.EXISTE_DDL.PDEPL_R, DDL_MECA),
                (OP.EXISTE_DDL.PCCONT_R, ECCONT),
                (OP.EXISTE_DDL.PCFROT_R, ECFROT),
            ),
        ),
    )


# ------------------------------------------------------------


class CMT63D(CMT33D):
    """
    THE CMT63D CLASS ELEMENT :
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D : elementary treatments
    Local Numerotation :

    Input parameters :

    Output parameters :
    """

    meshType = MT.TRIA6
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6)))


# ------------------------------------------------------------


class CMQ93D(CMT33D):
    """
    THE CMQ93D CLASS ELEMENT :
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D : elementary treatments
    Local Numerotation :

    Input parameters :

    Output parameters :
    """

    meshType = MT.QUAD9
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7, 8, 9)))


# ------------------------------------------------------------
class CMQ83D(CMT33D):
    """
    THE CMQ93D CLASS ELEMENT :
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D : elementary treatments
    Local Numerotation :

    Input parameters :

    Output parameters :
    """

    meshType = MT.QUAD8
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7, 8)))


# ------------------------------------------------------------
class CMQ43D(CMT33D):
    """
    THE CMQ93D CLASS ELEMENT :
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D : elementary treatments
    Local Numerotation :

    Input parameters :

    Output parameters :
    """

    meshType = MT.QUAD4
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)),)
