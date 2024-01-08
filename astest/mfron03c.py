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

import medcoupling as MEDC

from code_aster import CA
from code_aster.Commands import *

CA.init("--test", "--continue")

test = CA.TestCase()

expected = ["ElasticStrain" + cmp for cmp in ("XX", "YY", "ZZ", "XY", "XZ", "YZ")]
expected.extend([f"g[{i}]" for i in range(12)])
expected.extend([f"p[{i}]" for i in range(12)])
expected.extend([f"a[{i}]" for i in range(12)])

filename = "fort.80"
mesh = MEDC.MEDFileUMesh(filename)[0]
vari = MEDC.ReadFieldGauss(filename, mesh.getName(), 0, "MFRONT__VARI_ELGA_NOMME", 0, 0)
array = vari.getArray()
components = array.getInfoOnComponents()

test.assertEqual(len(components), len(expected), msg="check size")

test.assertSequenceEqual(components, expected, msg="check names of the components")

test.printSummary()

CA.close()
