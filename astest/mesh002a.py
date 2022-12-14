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

import os
import code_aster
from code_aster.Commands import *
from code_aster.Utilities.MedUtils.MEDPartitioner import MEDPartitioner
from code_aster import MPI


code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
print("Nb procs", MPI.ASTER_COMM_WORLD.Get_size())
print("Rank", MPI.ASTER_COMM_WORLD.Get_rank())

if MPI.ASTER_COMM_WORLD.Get_size() > 1:
    is_parallel = True
else:
    is_parallel = False


# Split the mesh

ms = MEDPartitioner("ssnv187a.mmed")

ms.partitionMesh(True)

# Where to save the mesh in a single folder
path = os.getcwd()
path.replace("/proc."+str(MPI.ASTER_COMM_WORLD.Get_rank()), "")

meshFolder = path+"/meshFolder"

try:
    os.mkdir(meshFolder)
except OSError:
    print ("Creation of the directory %s failed" % meshFolder)

# write the mesh in meshFolder
ms.writeMesh(meshFolder)

# 3 different way to read a Parallel Mesh

# 1) File by File after partioning (need a preliminary partioning )
pMesh1 = code_aster.ParallelMesh()
pMesh1.readMedFile(ms.writedFilename(), True)
pMesh1.checkConsistency("ssnv187a.mmed")

# 2) From a folder (need a preliminary partioning )
pMesh2 = code_aster.ParallelMesh()
pMesh2.readMedFile(meshFolder+"/ssnv187a_new_%d.med"%rank, True)
pMesh2.checkConsistency("ssnv187a.mmed")

# 3) Directely from a file (without preliminary partioning )
pMesh3 = code_aster.ParallelMesh()
pMesh3.readMedFile("ssnv187a.mmed")
pMesh3.checkConsistency("ssnv187a.mmed")

# 4) With LIRE_MAILLAGE (internal partitioning)
pMesh4 = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR="PTSCOTCH", INFO_MED=2)

model = AFFE_MODELE(MAILLAGE=pMesh4,
                    AFFE=_F(MODELISATION='3D', PHENOMENE='MECANIQUE', TOUT='OUI'),)


nbNodes = [79,78]
nbCells = [161,161]

test.assertEqual(pMesh1.getNumberOfNodes(), nbNodes[rank])
test.assertEqual(pMesh1.getNumberOfCells(), nbCells[rank])

test.assertEqual(pMesh1.getNumberOfNodes(), pMesh2.getNumberOfNodes())
test.assertEqual(pMesh1.getNumberOfNodes(), pMesh3.getNumberOfNodes())
test.assertEqual(pMesh1.getNumberOfNodes(), pMesh4.getNumberOfNodes())
test.assertEqual(pMesh1.getNumberOfCells(), pMesh2.getNumberOfCells())
test.assertEqual(pMesh1.getNumberOfCells(), pMesh3.getNumberOfCells())
test.assertEqual(pMesh1.getNumberOfCells(), pMesh4.getNumberOfCells())
test.assertEqual(pMesh1.getDimension(), 2)
test.assertEqual(pMesh2.getDimension(), 2)
test.assertEqual(pMesh3.getDimension(), 2)
test.assertEqual(pMesh4.getDimension(), 2)

test.assertEqual(is_parallel, True)

test.printSummary()


FIN()
