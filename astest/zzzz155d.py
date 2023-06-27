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

import os
import code_aster
from code_aster.Commands import *
from code_aster import MPI

code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()


def splitEntitySet(nbElemT, rank, nbProcs):
    nbElemL = int(nbElemT / nbProcs)
    start = rank * nbElemL + 1
    if rank == nbProcs - 1:
        end = nbProcs * nbElemL
        nbElemL = nbElemL - (end - nbElemT)
    return (nbElemL, start)


from code_aster.Utilities.MedUtils.MedMeshAndFieldsSplitter import splitMeshAndFieldsFromMedFile

filename = "zzzz155d.med"

myTuple = splitMeshAndFieldsFromMedFile(filename, True)

valuesDeplRef = [
    0.165282,
    0.068450,
    1.000000,
    0.000000,
    0.000000,
    0.000000,
    -0.165282,
    0.068450,
    1.000000,
    0.000000,
    0.000000,
    0.000000,
    0.160739,
    0.095945,
    0.483648,
    0.165282,
    -0.068450,
    1.000000,
    0.000000,
    0.000000,
    0.000000,
    -0.160739,
    0.095945,
    0.483648,
    -0.165282,
    -0.068450,
    1.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.077270,
    1.000000,
    0.160739,
    -0.095945,
    0.483648,
    -0.160739,
    -0.095945,
    0.483648,
    0.000000,
    0.094088,
    0.449709,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    -0.077270,
    1.000000,
    0.000000,
    -0.094088,
    0.449709,
]
profileDepl = {
    1: 0,
    2: 1,
    5: 2,
    6: 3,
    9: 4,
    10: 5,
    12: 6,
    13: 7,
    14: 8,
    16: 9,
    17: 10,
    18: 11,
    21: 12,
    22: 13,
    23: 14,
    25: 15,
    26: 16,
    27: 17,
}

splitMesh = myTuple[0]
loc2Glob = splitMesh.getLocalToGlobalMapping()
depl = myTuple[1]["00000008DEPL"][1]
valuesDepl = depl.getValues()
for count, i in enumerate(loc2Glob):
    if i in profileDepl:
        posInRef = profileDepl[i]
        # 3 components: DX, DY, DZ
        for j in range(3):
            test.assertAlmostEqual(
                abs(valuesDepl[3 * count + j] - valuesDeplRef[posInRef * 3 + j]), 0, delta=1e-6
            )

valuesSiefRef = [
    -2581878222.999054,
    1426692850.265594,
    209063210160.058929,
    583000884.506911,
    -552158012.720824,
    -2241766449.432333,
    190894053.658005,
    1867461588.894882,
    217865015676.953125,
    583000884.506947,
    -955687587.533880,
    -1694244534.090172,
    -2581878222.998993,
    1426692850.265610,
    209063210160.058899,
    -583000884.506966,
    -552158012.720879,
    2241766449.432312,
    190894053.658096,
    1867461588.894928,
    217865015676.953186,
    -583000884.506930,
    -955687587.533889,
    1694244534.090155,
    -7894824059.450272,
    -13660377933.541031,
    202943205173.981567,
    35478969.164802,
    -3566674632.839390,
    -2241766449.432330,
    -3479486036.766785,
    -9386955787.516769,
    213387576436.902222,
    35478969.164783,
    -3970204207.652491,
    -1694244534.090194,
    -7894824059.450348,
    -13660377933.541077,
    202943205173.981476,
    -35478969.164807,
    -3566674632.839486,
    2241766449.432309,
    -3479486036.766861,
    -9386955787.516800,
    213387576436.902161,
    -35478969.164825,
    -3970204207.652509,
    1694244534.090164,
    190894053.658112,
    1867461588.894943,
    217865015676.953247,
    -583000884.506927,
    955687587.533906,
    -1694244534.090168,
    -2581878222.999039,
    1426692850.265564,
    209063210160.058838,
    -583000884.506946,
    552158012.720845,
    -2241766449.432320,
    190894053.658035,
    1867461588.894882,
    217865015676.953186,
    583000884.506952,
    955687587.533923,
    1694244534.090177,
    -2581878222.999054,
    1426692850.265564,
    209063210160.058868,
    583000884.506937,
    552158012.720862,
    2241766449.432286,
    -3479486036.766739,
    -9386955787.516739,
    213387576436.902283,
    -35478969.164799,
    3970204207.652508,
    -1694244534.090211,
    -7894824059.450272,
    -13660377933.541061,
    202943205173.981445,
    -35478969.164802,
    3566674632.839468,
    -2241766449.432363,
    -3479486036.766769,
    -9386955787.516754,
    213387576436.902252,
    35478969.164809,
    3970204207.652543,
    1694244534.090156,
    -7894824059.450302,
    -13660377933.541061,
    202943205173.981506,
    35478969.164813,
    3566674632.839490,
    2241766449.432269,
    6913271115.483719,
    -455424696.902817,
    192527588153.695343,
    -130074581.823660,
    15385315994.308405,
    8487239739.273946,
    2781719930.776810,
    -4066677554.315063,
    182367003728.751221,
    -130074581.823687,
    1107966375.853988,
    8392018536.605709,
    6913271115.483612,
    -455424696.902863,
    192527588153.695312,
    130074581.823671,
    15385315994.308393,
    -8487239739.273944,
    2781719930.776703,
    -4066677554.315125,
    182367003728.751160,
    130074581.823645,
    1107966375.853958,
    -8392018536.605786,
    82345713997.895981,
    80371277905.696411,
    239405331799.198792,
    -34853379.155467,
    18399832614.427002,
    8487239739.273935,
    77928499205.184479,
    76093476629.606781,
    228959083766.250092,
    -34853379.155473,
    4122482995.972589,
    8392018536.605720,
    82345713997.895966,
    80371277905.696411,
    239405331799.198792,
    34853379.155468,
    18399832614.427013,
    -8487239739.273957,
    77928499205.184448,
    76093476629.606766,
    228959083766.250061,
    34853379.155461,
    4122482995.972581,
    -8392018536.605780,
    2781719930.776810,
    -4066677554.315048,
    182367003728.751282,
    130074581.823668,
    -1107966375.853936,
    8392018536.605700,
    6913271115.483704,
    -455424696.902817,
    192527588153.695374,
    130074581.823669,
    -15385315994.308357,
    8487239739.273921,
    2781719930.776764,
    -4066677554.315094,
    182367003728.751190,
    -130074581.823658,
    -1107966375.853941,
    -8392018536.605801,
    6913271115.483704,
    -455424696.902832,
    192527588153.695343,
    -130074581.823657,
    -15385315994.308359,
    -8487239739.273968,
    77928499205.184494,
    76093476629.606796,
    228959083766.250122,
    34853379.155467,
    -4122482995.972545,
    8392018536.605716,
    82345713997.895996,
    80371277905.696442,
    239405331799.198822,
    34853379.155468,
    -18399832614.426968,
    8487239739.273927,
    77928499205.184464,
    76093476629.606781,
    228959083766.250061,
    -34853379.155465,
    -4122482995.972554,
    -8392018536.605784,
    82345713997.895981,
    80371277905.696426,
    239405331799.198792,
    -34853379.155464,
    -18399832614.426968,
    -8487239739.273968,
]

# +24 because there is 24 seg2 before hexa8
# Element global numbering is made of: 1, ..., 24 seg2 et 25, ..., 32 hexa8
profileSief = {1 + 24: 0, 2 + 24: 1, 5 + 24: 2, 6 + 24: 3}

# Split cells like if they were read by MedFileReader
split = splitEntitySet(32, rank, size)
if rank == 0:
    cellGlobalId = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 25, 26, 27, 28]
elif rank == 1:
    cellGlobalId = [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 29, 30, 31, 32]
else:
    raise NameError("Test only on 2 procs")

cBalancer = myTuple[2]
# Balance cell global id over processes in order to be able to compare ELGA field
bCellGlobId = cBalancer.balanceVectorOverProcesses(cellGlobalId)
renumber = cBalancer.getRenumbering()
bCellGlobIdR = [0] * len(bCellGlobId)
for i in range(len(bCellGlobId)):
    bCellGlobIdR[renumber[i] - 1] = bCellGlobId[i]

sief = myTuple[1]["00000008SIEF_ELGA"][1]
cumSizes = sief.getCumulatedSizesVector()
valuesSief = sief.getValues()

nbCmp = 6 * 8
for i in range(len(bCellGlobId)):
    numGlob = bCellGlobIdR[i]
    if numGlob in profileSief:
        posInRef = profileSief[numGlob]
        posInNew = cumSizes[i]
        for j in range(nbCmp):
            test.assertAlmostEqual(
                abs(valuesSief[posInNew + j] - valuesSiefRef[(posInRef) * nbCmp + j]), 0, delta=1e-6
            )

FIN()
