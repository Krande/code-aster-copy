# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

# person_in_charge: francesco.bettonte at edf.fr

import libaster
from ..logger import logger

try :
    import medcoupling as mc
    HAS_MEDCOUPLING = True
    
except ImportError:
    HAS_MEDCOUPLING = False


def convertMesh2MedCoupling(asmesh):
    
    if not HAS_MEDCOUPLING:
        logger.info('The MEDCoupling module is missing')
        return
    
    cells, groups_c, groups_n = libaster.getMedCouplingConversionData(asmesh)
    
    mcmesh = mc.MEDFileUMesh()
    coords = mc.DataArrayDouble(asmesh.getCoordinates().getValues(),
                                asmesh.getNumberOfNodes(),
                                asmesh.getDimension())
        
    maxdim = max(cells.keys())
    levels = {i : i-maxdim for i in range(maxdim,-1,-1)}

    # Creation du maillage par niveau, depart par le plus haut (0)
    for dim in sorted(cells.keys())[::-1]:
        
        mesh_at_current_level = mc.MEDCouplingUMesh(asmesh.getName(), dim)
        mesh_at_current_level.setCoords(coords)
        
        conn, connI = cells[dim]
        mesh_at_current_level.setConnectivity(mc.DataArrayInt(conn),
                                              mc.DataArrayInt(connI))
        
        o2n = mesh_at_current_level.sortCellsInMEDFileFrmt()
        mesh_at_current_level.checkConsistencyLight()
        mcmesh.setMeshAtLevel(levels[dim], mesh_at_current_level)
        
        # Groupes de mailles
        try :
            groups_c_at_level = []
            for group_name, group_cells in groups_c[dim].items():
                group_medcoupling = mc.DataArrayInt(group_cells)
                group_medcoupling.transformWithIndArr(o2n)
                group_medcoupling.setName(group_name)
                groups_c_at_level.append(group_medcoupling)
            mcmesh.setGroupsAtLevel(levels[dim], groups_c_at_level)

        except KeyError:
            # Pas de groupes à ce niveau
            pass

    # Groupes de noeuds au niveau 1
    groups_n_at_level = []
    for group_name, group_nodes in groups_n.items():
        group_medcoupling = mc.DataArrayInt(group_nodes)
        group_medcoupling.setName(group_name)
        groups_n_at_level.append(group_medcoupling)
    mcmesh.setGroupsAtLevel(1, groups_n_at_level)

    mcmesh.setName(asmesh.getName())
    mcmesh.rearrangeFamilies()

    return mcmesh
