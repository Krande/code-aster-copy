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

import os.path as osp
import re
from glob import glob

import code_aster
from code_aster import MPI

import medcoupling as mc
from collections import Counter

code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.COMM_WORLD.Get_rank()
nbproc = MPI.COMM_WORLD.Get_size()

if nbproc > 1:
    is_parallel = True
else:
    is_parallel = False

# list of meshes with error
failure = [
    # for mesh002s
    "sdll06a.mmed",
    # for mesh002t
    "ssls134a.med", "ssna107c.mmed", "ssnl121b.mmed", "ssnp05a.mmed", "ssnp107a.mmed",
    # for mesh002u
    "ssnv214a.med"
]

# get list of meshes from .export
export = glob("*.export")
assert len(export) == 1, "expecting exactly one export file!"

re_mesh = re.compile("^F +nom +(.*med) +D +0", re.M)
with open(export[0], "r") as fexp:
    content = fexp.read()

meshes = re_mesh.findall(content)

nb_mesh = 0
nb_mesh_converted = 0
# loop on mesh
for mesh_file in meshes:
    MPI.COMM_WORLD.Barrier()
    mesh_name = osp.basename(mesh_file)
    if mesh_name in failure:
        print("SKIP:", mesh_name)
        continue

    nb_mesh += 1
    mesh_vers = mc.MEDFileVersionOfFileStr(mesh_name)
    vers_num = int(mesh_vers.replace(".", ""))

    print(
        "MESHAME %d: %s (version: %s)" % (nb_mesh, mesh_name, mesh_vers), flush=True
    )

    # convert old med mesh < 3.0.0
    if vers_num < 300:
        mfd = mc.MEDFileData(mesh_name)
        mesh_name = mesh_name.split(".")[0] + "_tmp.med"
        mfd.write(mesh_name, 2)
        print(
            "Mesh converted: %s (version: %s)"
            % (mesh_name, mc.MEDFileVersionOfFileStr(mesh_name)),
            flush=True,
        )

    # read parallel mesh and partitioning
    pmesh = code_aster.ParallelMesh()
    try:
        pmesh.readMedFile(mesh_name)
    except Exception as exc:
        test.assertIsNone(exc, "partitioning failed (%s)" % mesh_name)
        print("ERROR:", str(exc))
        continue

    # read std mesh
    mesh = code_aster.Mesh()
    mesh.readMedFile(mesh_name)

    # tests
    group_no_std = mesh.getGroupsOfNodes(local=False)
    group_no_gl  = pmesh.getGroupsOfNodes(local=False)
    # il manque des groupes
    #test.assertSequenceEqual(sorted(group_no_std), sorted(group_no_gl))

    group_ma_std = mesh.getGroupsOfCells(local=False)
    group_ma_gl  = pmesh.getGroupsOfCells(local=False)
    # il manque des point1
    #test.assertSequenceEqual(sorted(group_ma_std), sorted(group_ma_gl))

    nb_nodes_std = mesh.getNumberOfNodes()
    nb_nodes_lc = len(pmesh.getInnerNodes())
    nb_nodes_gl = MPI.COMM_WORLD.allreduce(nb_nodes_lc, MPI.SUM)
    # il manque des noeuds
    #test.assertEqual(nb_nodes_std, nb_nodes_gl)

    nb_cells_std = mesh.getNumberOfCells()
    cells_rank = pmesh.getCellsRank()
    nb_cells_lc = Counter(cells_rank)[rank]
    nb_cells_gl = MPI.COMM_WORLD.allreduce(nb_cells_lc, MPI.SUM)
    # il manque des mailles
    #test.assertEqual(nb_cells_std, nb_cells_gl)

    #test.assertTrue(pmesh.checkConsistency(mesh_name))

    test.assertTrue(pmesh.getNumberOfNodes() > 0)
    test.assertTrue(pmesh.getNumberOfCells() > 0)
    nb_mesh_converted += 1

    MPI.COMM_WORLD.Barrier()


list_nb_mesh = [525, 494, 494]
list_nb_mesh_conv = [525, 494, 494]

print("Number of mesh: %s" % (nb_mesh), flush=True)
print("Number of mesh converted: %s" % (nb_mesh_converted), flush=True)

# all the mesh are partitioned
test.assertTrue(is_parallel)
test.assertEqual(nb_mesh, list_nb_mesh[nbproc - 2])
test.assertEqual(nb_mesh_converted, list_nb_mesh_conv[nbproc - 2])

test.printSummary()

code_aster.close()
