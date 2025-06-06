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

import platform

from code_aster.Commands import *
from code_aster import CA

CA.init("--test", debug=True, ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

# calling DEBUT/init another time must not fail
DEBUT(CODE="OUI")

from code_aster.Utilities.ExecutionParameter import ExecutionParameter

params = ExecutionParameter()

test.assertEqual(params.get_option("hostname"), platform.node())
# GreaterEqual because of multiplicative factor running testcases
# timelimit = 123. - 10% (default time saved)
test.assertGreaterEqual(params.get_option("tpmax"), 123 * 0.9)

params.set_option("tpmax", 60.0)
test.assertEqual(params.get_option("tpmax"), 60)

# FIN does nothing: done when libaster is lost
FIN()

test.printSummary()
