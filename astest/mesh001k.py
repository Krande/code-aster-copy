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

import os.path as osp

from code_aster.Commands import *
from code_aster import CA
from medcoupling import *

from code_aster.Utilities import SharedTmpdir
from code_aster.CA import MPI


CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))


def create_mesh(nb_node_dir, meshfile):
    # 20 --> 6859 mailles
    # 35 --> 39304 mailles
    arr = DataArrayDouble(nb_node_dir)
    arr.iota()

    mcube = MEDCouplingCMesh()
    mcube.setCoords(arr, arr, arr)
    mesh = mcube.buildUnstructured()

    ls = []
    for i in range(mesh.getNumberOfCells()):
        grp = DataArrayInt([i])
        grp.setName("GrpWithLongName_%d" % (i + 1))
        ls.append(grp)

    mm = MEDFileUMesh()
    mm.setName("CUBE")
    mm[0] = mesh
    mm.setGroupsAtLevel(0, ls)

    mm.write(meshfile, 2)


# check ParallelMesh object API
test = CA.TestCase()

# MPI test
rank = MPI.ASTER_COMM_WORLD.Get_rank()
nbproc = MPI.ASTER_COMM_WORLD.Get_size()


# from MED format (only this one a ParallelMesh)
mesh = CA.ParallelMesh()

with SharedTmpdir("mesh001k_") as tmpdir:
    print("Create Mesh", flush=True)
    medfile = osp.join(tmpdir.path, "mesh001k.med")
    if rank == 0:
        create_mesh(10, medfile)
    MPI.ASTER_COMM_WORLD.Barrier()

    print("Read Mesh", flush=True)
    mesh.readMedFile(medfile)

    medfile2 = osp.join(tmpdir.path, "mesh001k_2.med")

    print("Save Mesh", flush=True)
    mesh.printMedFile(medfile2, False)

    test.assertTrue(mesh.isParallel())

    print("Read full Mesh", flush=True)
    mesh_std = CA.Mesh()
    mesh_std.readMedFile(medfile)

    mesh_relu = CA.Mesh()
    mesh_relu.readMedFile(medfile2)

    # test mesh
    test.assertEqual(mesh_std.getNumberOfNodes(), mesh_relu.getNumberOfNodes())
    test.assertEqual(mesh_std.getNumberOfCells(), mesh_relu.getNumberOfCells())
    test.assertSequenceEqual(
        sorted(mesh_std.getGroupsOfNodes()), sorted(mesh_relu.getGroupsOfNodes())
    )
    test.assertSequenceEqual(
        sorted(mesh_std.getGroupsOfCells()), sorted(mesh_relu.getGroupsOfCells())
    )

test.printSummary()

CA.close()
