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

import code_aster
from code_aster.Commands import *

code_aster.init("--test")

test = code_aster.TestCase()

list_t = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=10.0, PAS=1.0))

list_i = DEFI_LIST_INST(METHODE="MANUEL", DEFI_LIST=_F(LIST_INST=list_t))

test.assertEqual(list_i.size(), 11, msg="nbsteps")
test.assertIsNone(list_i.getInitial())
test.assertAlmostEqual(list_i.getFinal(), 10.0)

code_aster.close()
