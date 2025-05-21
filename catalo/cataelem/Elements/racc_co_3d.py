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

# ELEMENTARY TREATMENT OF 3D FRICTIONLESS ELEMENT WITH DEFI_CONTACT OPERATOR

DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        # Slave nodes with LAG_C
        ("EN1", ("DX", "DY", "DZ", "DRX", "DRY", "DRZ")),
        # Slave nodes without LAG_C
        ("EN2", ()),
        # Master nodes
        ("EN3", ("DX", "DY", "DZ")),
    ),
)

NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))
MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)
# ----------------
# Modes locaux :
# ----------------
# ------------------------------------------------------------
class RACS2T3(Element):
    """
    THE RACS2T3 CLASS ELEMENT : SEG2/TRIA3 (2D EDGE / 3D FACE )
    """

    meshType = MT.SE2TR3
    nodes = (
        SetOfNodes("EN1", (1, 2)),
        SetOfNodes("EN2", ()),
        SetOfNodes("EN3", (3, 4, 5)),
    )
    calculs = (
        OP.CHAR_MECA_CONT(
            te=355,
            para_in=(
                (SP.PGEOMER, NGEOMER),
            ),
            para_out=((SP.PVECTCR, MVECTUR),),
        ),
        
    )

