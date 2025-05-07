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


from medcoupling import *

from code_aster.Commands import *
from code_aster import CA

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


test = CA.TestCase()


medfile = "mesh001l.mmed"
create_mesh(20, medfile)

DEFI_FICHIER(UNITE=80, FICHIER=medfile, TYPE="BINARY", ACCES="OLD")

mesh = LIRE_MAILLAGE(UNITE=80, FORMAT="MED")

DEFI_FICHIER(ACTION="LIBERER", UNITE=80)

# test mesh
grpCells = mesh.getGroupsOfCells()

for grp in grpCells:
    test.assertTrue(len(mesh.getCells(grp)) == 1)

test.printSummary()
CA.close()
