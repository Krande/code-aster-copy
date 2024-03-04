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

test = CA.TestCase()

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

# list of DataStructures
func = [None, None]
func[0] = DEFI_FONCTION(VALE=(0, 0, 1, 1), NOM_PARA="X")
func[1] = DEFI_FONCTION(VALE=(0, 0, 1, 2), NOM_PARA="X")

test.assertEqual(func[0](1), 1)
test.assertEqual(func[1](1), 2)

# dict of DataStructures
results = dict(
    one=DEFI_FONCTION(VALE=(0, 0, 1, 1), NOM_PARA="X"),
    two=DEFI_FONCTION(VALE=(0, 0, 1, 2), NOM_PARA="X"),
)

test.assertEqual(results["one"](1), 1)
test.assertEqual(results["two"](1), 2)

# tuple of DataStructures
tup = tuple(func)

# dict of list of DataStructures
dfunc = dict(flist=func)

# list of dict of list of DataStructures
ldfunc = [dfunc]

# list of ...x12... dict of list of DataStructures
l12dfunc = dfunc
for _ in range(12):
    l12dfunc = [l12dfunc]

CA.saveObjects()

test.printSummary()

FIN()
