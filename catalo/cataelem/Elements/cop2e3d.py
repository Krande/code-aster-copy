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

# ----------------
# Modes locaux :
# ----------------


DDL_MECA = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ", "LAGS_C"))


# ------------------------------------------------------------
class COP2E3D(Element):
    """Please document this element"""

    meshType = MT.SEG2
    calculs = (OP.EXISTE_DDL(te=99, para_out=((OP.EXISTE_DDL.PDEPL_R, DDL_MECA),)),)


# ------------------------------------------------------------
class COP3E3D(COP2E3D):
    """Please document this element"""

    meshType = MT.SEG3


# ------------------------------------------------------------
class COQ4E3D(COP2E3D):
    """Please document this element"""

    meshType = MT.QUAD4


# ------------------------------------------------------------
class COT3E3D(COP2E3D):
    """Please document this element"""

    meshType = MT.TRIA3


# ------------------------------------------------------------
class COQ8E3D(COP2E3D):
    """Please document this element"""

    meshType = MT.QUAD8


# ------------------------------------------------------------
class COT6E3D(COP2E3D):
    """Please document this element"""

    meshType = MT.TRIA6


# ------------------------------------------------------------
class COQ9E3D(COP2E3D):
    """Please document this element"""

    meshType = MT.QUAD9
