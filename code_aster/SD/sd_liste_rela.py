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

from . import *


class sd_liste_rela(AsBase):
    nomj = SDNom(fin=19)
    RLBE = AsVR(SDNom(debut=19))
    RLSU = AsVI(SDNom(debut=19))
    RLTC = AsVK8(SDNom(debut=19), lonmax=1)
    RLNO = AsVK8(SDNom(debut=19))
    RLCO = AsVR(SDNom(debut=19))
    RLNT = AsVI(SDNom(debut=19))
    RLPO = AsVI(SDNom(debut=19))
    RLNR = AsVI(SDNom(debut=19), lonmax=1)
    RLTV = AsVK8(SDNom(debut=19), lonmax=1)
    RLDD = AsVK8(SDNom(debut=19))
