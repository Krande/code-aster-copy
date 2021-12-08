# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

DEBUT(CODE=_F(NIV_PUB_WEB="INTERNET"))

pres = [FORMULE(NOM_PARA="X", VALE="X"), (FORMULE(NOM_PARA="X", VALE="X * 2"),)]
a = DEFI_CONSTANTE(VALE=1)
b = DEFI_CONSTANTE(VALE=2)


def Y0(X):
    return 0.0


form = FORMULE(NOM_PARA="X", VALE="Y0(X)", Y0=Y0)

# pickle raises an error because Y0 is not the same object as in formula context
Y0 = 0.123456

FIN()
