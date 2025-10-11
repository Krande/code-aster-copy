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

# This test checks mesh spliting with more than only one ghost cell layers

from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI
import numpy as np

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()


def checkJoinSize(mesh, size, checker):
    graph = CA.CommGraph()

    dictProc = {}
    cmpt = 0
    for proc in mesh.getOppositeDomains():
        graph.addCommunication(proc)
        dictProc[proc] = cmpt
        cmpt += 1
    graph.synchronizeOverProcesses()

    matchings = graph.getMatchings()
    for procT in matchings:
        proc = procT[1]
        tag = procT[0]
        if proc == -1:
            continue
        numJoint = dictProc[proc]
        fJ = mesh.getSendJoint(numJoint)
        checker.assertEqual(len(fJ) / 2, size)
        sJ = mesh.getReceiveJoint(numJoint)
        checker.assertEqual(len(sJ) / 2, size)


# Mesh reading 1
mesh1 = CA.ParallelMesh()
mesh1.readMedFile("zzzz155k.med", ghost=2)

# Check extraction of the last layer of ghosts DOFs
model = CA.Model(mesh1)
model.addModelingOnMesh(CA.Physics.Thermal, CA.Modelings.Planar)
model.build()
numeDDL = NUME_DDL(MODELE=model)

test.assertEqual(
    numeDDL.getEquationNumbering().getGhostDOFs(lastLayerOnly=True),
    [[3, 12, 13], [10, 16, 17]][rank],
)


checkJoinSize(mesh1, 3 * 2, test)

test.assertTrue(mesh1.checkConsistency("zzzz155k.med"))
test.assertTrue(mesh1.checkJoints())

mesh1r = mesh1.refine()
checkJoinSize(mesh1r, 5 * 2, test)
test.assertTrue(mesh1r.checkJoints())

mesh1r2 = mesh1.refine(2)
checkJoinSize(mesh1r2, 9 * 2, test)
test.assertTrue(mesh1r2.checkJoints())

# Mesh reading 2
mesh2 = CA.ParallelMesh()
mesh2.readMedFile("zzzz155k.med", ghost=3)

checkJoinSize(mesh2, 9, test)

test.assertTrue(mesh2.checkConsistency("zzzz155k.med"))
test.assertTrue(mesh2.checkJoints())

mesh2r = mesh2.refine()
checkJoinSize(mesh2r, 15, test)
test.assertTrue(mesh2r.checkJoints())

model = AFFE_MODELE(
    MAILLAGE=mesh2, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="C_PLAN")
)

acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3))

mater = AFFE_MATERIAU(MAILLAGE=mesh2, AFFE=_F(TOUT="OUI", MATER=acier))

cl = AFFE_CHAR_CINE(
    MODELE=model, MECA_IMPO=(_F(GROUP_NO="DA", DX=0.0, DY=0.0, DZ=0.0), _F(GROUP_NO="BC", DX=1.0))
)

resu = MECA_STATIQUE(MODELE=model, CHAM_MATER=mater, EXCIT=_F(CHARGE=cl), INST=0.0)

depl = resu.getField("DEPL", 1)
sDepl = depl.toSimpleFieldOnNodes()

if rank == 0:
    check = np.array(
        [
            [0.00000000e00, 0.00000000e00],
            [0.00000000e00, 0.00000000e00],
            [0.00000000e00, 0.00000000e00],
            [8.56316656e-01, 4.31054563e-02],
            [7.12634944e-01, 4.31043568e-02],
            [5.68946417e-01, 4.30843730e-02],
            [4.25173217e-01, 4.31068383e-02],
            [2.81594456e-01, 4.40711657e-02],
            [1.42246106e-01, 4.44294786e-02],
            [1.42246106e-01, -4.44294786e-02],
            [2.81594456e-01, -4.40711657e-02],
            [4.25173217e-01, -4.31068383e-02],
            [5.68946417e-01, -4.30843730e-02],
            [7.12634944e-01, -4.31043568e-02],
            [8.56316656e-01, -4.31054563e-02],
            [8.56316253e-01, -2.04376200e-16],
            [7.12631098e-01, -2.35889669e-16],
            [5.68959081e-01, -4.19609244e-16],
            [4.25370993e-01, -4.28642155e-16],
            [2.81292429e-01, -2.70661535e-16],
            [1.32586662e-01, -1.25614964e-16],
        ]
    )

    test.assertAlmostEqual(abs(sDepl.toNumpy()[0] - check).max(), 0.0)
elif rank == 1:
    check = np.array(
        [
            [4.25173217e-01, -4.31068383e-02],
            [5.68946417e-01, -4.30843730e-02],
            [7.12634944e-01, -4.31043568e-02],
            [8.56316656e-01, -4.31054563e-02],
            [1.00000000e00, -4.31050982e-02],
            [1.00000000e00, -1.42614445e-16],
            [8.56316253e-01, -2.04376200e-16],
            [7.12631098e-01, -2.35889669e-16],
            [5.68959081e-01, -4.19609244e-16],
            [4.25370993e-01, -4.28642155e-16],
            [2.81292429e-01, -2.70661535e-16],
            [1.32586662e-01, -1.25614964e-16],
            [1.00000000e00, 4.31050982e-02],
            [8.56316656e-01, 4.31054563e-02],
            [7.12634944e-01, 4.31043568e-02],
            [5.68946417e-01, 4.30843730e-02],
            [4.25173217e-01, 4.31068383e-02],
            [2.81594456e-01, 4.40711657e-02],
            [1.42246106e-01, 4.44294786e-02],
            [1.42246106e-01, -4.44294786e-02],
            [2.81594456e-01, -4.40711657e-02],
        ]
    )

    test.assertAlmostEqual(abs(sDepl.toNumpy()[0] - check).max(), 0.0)
else:
    raise NameError("Test must run on 2 procs")

# Mesh reading 3
mesh3 = CA.ParallelMesh()
mesh3.readMedFile("zzzz155k.med", ghost=10)

FIN()
