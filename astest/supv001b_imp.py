# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

from code_aster.Commands import DEFI_CONSTANTE


def include_function(COUCHE):

    print("COUCHE DEBUT INCLUDE = ", COUCHE)

    COUCHE = COUCHE + 1
    print("COUCHE FIN INCLUDE = ", COUCHE)

    return COUCHE


def one(_):
    return 1.0


# WARNING: it is not recommended to create objects in the global namespace
# of an imported module (may be it never destroyed)
# Here it is to check reloading.
two = DEFI_CONSTANTE(VALE=2.0)
