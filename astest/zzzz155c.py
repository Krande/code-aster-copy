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

from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI

from zzzz155c_data import valuesDeplRef, valuesSiefRef

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()


def splitEntitySet(nbElemT, rank, nbProcs):
    nbElemL = int(nbElemT / nbProcs)
    start = rank * nbElemL + 1
    if rank == nbProcs - 1:
        end = nbProcs * nbElemL
        nbElemL = nbElemL - (end - nbElemT)
    return (nbElemL, start)


from code_aster.MedUtils import splitMeshAndFieldsFromMedFile

# Read med file and split mesh and fields
filename = "fort.20"
myTuple = splitMeshAndFieldsFromMedFile(filename, True)

vectDepl = myTuple[1]["00000009DEPL"][1]
valuesDepl = vectDepl.getValues()

mesh = myTuple[0]

# Check DEPL field compared to file "mdump3"
loc2Glob = mesh.getLocalToGlobalNodeIds()
for count, i in enumerate(loc2Glob):
    # 3 components: DX, DY, DZ
    for j in range(3):
        test.assertAlmostEqual(
            abs(valuesDepl[3 * count + j] - valuesDeplRef[i * 3 + j]), 0, delta=1e-6
        )

# Read mesh from file to read cell number
fr = CA.MedFileReader()
fr.openParallel(filename, CA.MedFileAccessType.MedReadOnly)
medMesh = fr.getMesh(0)
nbSeq = medMesh.getSequenceNumber()
seq = medMesh.getSequence(0)
nbIter = medMesh.getCellTypeNumberAtSequence(seq[0], seq[1])
assert nbIter == 1
nbCells = medMesh.getCellNumberAtSequence(seq[0], seq[1], 1)

# Split cells like if they were read by MedFileReader
split = splitEntitySet(nbCells, rank, size)
cellGobalId = [i for i in range(split[1], split[0] + split[1])]

cBalancer = myTuple[2]
# Balance cell global id over processes in order to be able to compare ELGA field
bCellGlobId = cBalancer.balanceVectorOverProcesses(cellGobalId)

vectSief = myTuple[1]["00000009SIEF_ELGA"][1]
valuesSief = vectSief.getValues()

# 6 components * 27 Gauss points
nbCmpXnbGP = 6 * 27
# Check SIEF_ELGA field compared to file "mdump3"
for count, i in enumerate(bCellGlobId):
    # 6 components * 27 Gauss points
    for j in range(nbCmpXnbGP):
        test.assertAlmostEqual(
            abs(valuesSief[nbCmpXnbGP * count + j] - valuesSiefRef[(i - 1) * nbCmpXnbGP + j]),
            0,
            delta=1e-6,
        )

FIN()
