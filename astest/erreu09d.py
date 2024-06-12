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

from code_aster.Commands import *
from code_aster import CA
from run_aster.config import CFG

CA.init("--test")

test = CA.TestCase()

if CFG.get("require_mpiexec"):
    print("INFO: This testcase can not be run under mpiexec.")
    test.assertTrue(True)
    test.printSummary()
    # because it would exit before saving the database
    CA.exit()

nbsteps = RESU.getLastIndex()

# may be 8 to 20 depending on the velocity of the host!
test.assertGreater(nbsteps, 2, msg="number of steps calculated before time limit")
test.printSummary()

CA.close()
