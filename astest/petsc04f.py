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

import os.path as osp

from code_aster.Commands import *
from code_aster.Utilities import SharedTmpdir

POURSUITE(CODE="OUI", ERREUR=_F(ALARME="EXCEPTION"))

with SharedTmpdir("petsc04f_") as tmpdir:
    medfile = osp.join(tmpdir.path, "petsc04f.med")
    DEFI_FICHIER(UNITE=87, FICHIER=medfile, TYPE="BINARY")
    IMPR_RESU(
        FICHIER_UNIQUE="OUI", FORMAT="MED", UNITE=87, RESU=_F(RESULTAT=MESTAT), VERSION_MED="4.1.0"
    )

    DEFI_FICHIER(ACTION="LIBERER", UNITE=87)

FIN()
