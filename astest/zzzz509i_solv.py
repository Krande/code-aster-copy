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

test = CA.TestCase()

DEBUT(CODE="OUI")

with test.assertRaisesRegex(
    CA.AsterError, "Le nombre maximal autorisé SUBD_NIVEAU.*de niveaux de subdivision est dépassé."
):
    MECA_NON_LINE(
        INCREMENT=_F(LIST_INST=times, INST_FIN=time_inter),
        MODELE=model,
        CHAM_MATER=material,
        EXCIT=(_F(CHARGE=sym_bottom), _F(CHARGE=load, FONC_MULT=mult_func)),
        COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
    )

test.printSummary()

FIN()
