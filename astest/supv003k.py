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

import os.path as osp
from itertools import chain

from mpi4py import MPI

# Get the rank and size in the original communicator
comm = MPI.COMM_WORLD
global_rank = comm.rank
global_size = comm.size

assert global_size == 3, "This testcase should be executed with 3 processes."

# Group processes as: #2 for a "supervisor" and (#0, #1) for code_aster
# Keep #0 for code_aster since its 'fort.6' is used for automatic testing
superv = global_rank == 2
for_aster = not superv
color = 0 if superv else 1

subcomm = comm.Split(color, global_rank)
rank = subcomm.rank
size = subcomm.size

print(f"World rank/size: {global_rank:2d}/{global_size} - Sub rank/size {rank:2d}/{size}")
print(f"Working group #{color}")

if superv:
    print("working in the other group")
    # supposing another mesh is used here!
    nodes = [999, 998]

if for_aster:
    print("starting code_aster", flush=True)
    from code_aster import rc

    rc.initialize = False
    from code_aster import CA

    test = CA.TestCase()

    CA.init("--test", comm=subcomm, ERREUR=_F(ALARME="EXCEPTION"))

    mesh = CA.ParallelMesh()
    mesh.readMedFile(CA.basedir / "mesh004b" / f"{rank}.med", partitioned=True, verbose=1)

    nb_nodes = mesh.getNumberOfNodes()
    test.assertEqual(nb_nodes, 6)
    nodes = mesh.getNodes(localNumbering=False)
    print(f"List of the {nb_nodes} nodes:", nodes)

    # 'intra' communication
    all_nodes = subcomm.allgather(nodes)
    all_nodes = set(chain.from_iterable(all_nodes))
    print("List of all code_aster nodes:", all_nodes)
    test.assertEqual(len(all_nodes), 8)

# global gather
all_nodes = comm.allgather(nodes)
all_nodes = set(chain.from_iterable(all_nodes))
print("List of all nodes:", all_nodes)
assert len(all_nodes) == 10, all_nodes
empty = all_nodes.difference(range(8)).difference([998, 999])
assert len(empty) == 0, empty

print("All processes are terminating!")
if superv:
    assert not osp.exists("glob.1")
    print(" OK  no TestCase object on this process")

if for_aster:
    test.printSummary()
    assert osp.exists("glob.1")
    CA.close()
