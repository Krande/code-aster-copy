# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`Mesh` --- Assignment of mesh
************************************************************************
"""

import aster

from ..Commands import CREA_MAILLAGE
from ..Objects import Mesh, PythonBool
from ..Objects.Serialization import InternalStateBuilder
from ..Utilities import injector
from ..Utilities.MedUtils.MEDConverter import convertMesh2MedCoupling
from . import mesh_builder


class MeshStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *Mesh*."""

    def restore(self, mesh):
        """Restore the *DataStructure* content from the previously saved internal
        state.

        Arguments:
            mesh (*DataStructure*): The *DataStructure* object to be pickled.
        """
        super().restore(mesh)
        mesh.build()


@injector(Mesh)
class ExtendedMesh:
    cata_sdj = "SD.sd_maillage.sd_maillage"
    internalStateBuilder = MeshStateBuilder

    buildSquare = classmethod(mesh_builder.buildSquare)
    buildDisk = classmethod(mesh_builder.buildDisk)
    buildCube = classmethod(mesh_builder.buildCube)
    buildCylinder = classmethod(mesh_builder.buildCylinder)

    def LIST_GROUP_NO(self):
        """Retourne la liste des groupes de noeuds sous la forme :
        [ (gno1, nb noeuds  gno1), ...]"""
        dic_gpno = self.sdj.GROUPENO.get()
        if dic_gpno is None:
            return []
        return [(gpno.strip(), len(dic_gpno[gpno])) for gpno in dic_gpno]

    def LIST_GROUP_MA(self):
        """Retourne la liste des groupes de mailles sous la forme :
        [ (gma1, nb mailles gma1, dime max des mailles gma1), ...]"""
        ltyma = aster.getvectjev("&CATA.TM.NOMTM")
        catama = aster.getcolljev("&CATA.TM.TMDIM")
        dic_gpma = self.sdj.GROUPEMA.get()
        if dic_gpma is None:
            return []
        dimama = [catama[ltyma[ma - 1].ljust(24)][0] for ma in self.sdj.TYPMAIL.get()]
        ngpma = []
        for grp in list(dic_gpma.keys()):
            dim = max([dimama[ma - 1] for ma in dic_gpma[grp]])
            ngpma.append((grp.strip(), len(dic_gpma[grp]), dim))
        return ngpma

    def buildFromMedCouplingMesh(self, mcmesh, verbose=0):
        """Build mesh from medcoupling mesh.

        Arguments:
            mcmesh (*medcoupling.MEDFileUMesh*): The MEDCoupling mesh.
            verbose (int) : 0 - warnings
                            1 - informations about main steps
                            2 - informations about all steps
        """
        mesh_builder.buildFromMedCouplingMesh(self, mcmesh, verbose)

    def readMedFile(self, filename, meshname=None, verbose=1):
        """Read a MED file containing a mesh.

        Arguments:
            filename (string): name of the MED file
            meshname (str): Name of the mesh to be read from file.
            verbose (int) : 0 - warnings
                            1 - informations about main steps
                            2 - informations about all steps
        """
        mesh_builder.buildFromMedFile(self, filename, meshname, verbose)

    def refine(self, ntimes=1, info=1):
        """Refine the mesh uniformly. Each edge is split in two.

        Arguments:
            ntimes [int] : the number of times the mesh is to be refined.
            info [int] : verbosity mode (1 or 2). Default 1.

        Returns:
            Mesh: the refined mesh.
        """

        new_mesh = CREA_MAILLAGE(
            MAILLAGE=self, RAFFINEMENT=_F(TOUT="OUI", NIVEAU=ntimes), INFO=info
        )
        return new_mesh

    def createMedCouplingMesh(self):
        """Returns the MEDCoupling unstructured mesh associated to the current mesh.

        Returns:
            Mesh: The MEDCoupling unstructured mesh associated to the current mesh.
        """

        return convertMesh2MedCoupling(self)

    def getNodes(self, group_name="", localNumbering=True, same_rank=None):
        """Return the list of the indexes of the nodes that belong to a group of nodes.

        Arguments:
            group_name (str): Name of the group (default: "").
            localNumbering (bool): not used (for compatibilty with ParallelMesh)
            same_rank (bool): not used (for compatibilty with ParallelMesh)

        Returns:
            list[int]: Indexes of the nodes of the group.
        """

        val = {None: PythonBool.NONE, True: PythonBool.TRUE, False: PythonBool.FALSE}

        return self._getNodes(group_name, localNumbering, val[same_rank])

    def getNodesFromCells(self, group_name, localNumbering=True, same_rank=None):
        """Returns the nodes indexes of a group of cells.

        Arguments:
            group_name (str): Name of the group.
            localNumbering (bool): not used (for compatibilty with ParallelMesh)
            same_rank (bool): not used (for compatibilty with ParallelMesh)

        Returns:
            list[int]: Indexes of the nodes of the group.
        """

        val = {None: PythonBool.NONE, True: PythonBool.TRUE, False: PythonBool.FALSE}

        return self._getNodesFromCells(group_name, localNumbering, val[same_rank])
