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

import matplotlib
import numpy as np
import numpy.ma as ma

matplotlib.use("TkAgg")

CA.init("--test", "--continue")

test = CA.TestCase()

# unittests for ComponentOnCells objects with masked values
sief = DEP2.getField("SIEF_ELGA", 1)
sief_s = sief.toSimpleFieldOnCells()
sixx = sief_s.SIXX
siyy = sief_s.SIYY
vy = sief_s.VY
mfz = sief_s.MFZ

cmps = sief_s.getComponents()
test.assertEqual(len(cmps), 12, msg="cmps")

test.assertAlmostEqual(vy.min(), -1.0, 8, msg="min")
test.assertAlmostEqual(vy.max(), -1.0, 8, msg="max")
test.assertAlmostEqual(vy.mean(), -1.0, 8, msg="mean")

no_error = mfz / vy
zero = abs(no_error + mfz).max()
test.assertAlmostEqual(zero, 0.0, 8, msg="no error division")

values = vy.values
test.assertAlmostEqual(vy.sum(), -6.0, 8, msg="sum: -1.0 * 6")

m_one = np.ones(values.shape, dtype=float) * -1.0
test.assertAlmostEqual(m_one.sum(), -69.0, 8, msg="-1.0 * 69")

vy.values = m_one
test.assertAlmostEqual(vy._values.sum(), -6.0, 8, msg="sum: only 6 values assigned, other 0.0")
test.assertAlmostEqual(vy.sum(), -6.0, 8, msg="sum")

test.printSummary()

CA.close()
