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


class sd_prof_vgen(AsBase):
    nomj = SDNom(fin=19)
    PRNO = AsColl(SDNom(debut=19), acces="NU", stockage="DISPERSE", modelong="VARIABLE", type="I")
    LILI = AsObject(SDNom(debut=19), genr="N", xous="S", type="K", ltyp=8)
    NUEQ = AsVI(SDNom(debut=19))
    DEEQ = AsVI(SDNom(debut=19))
    DELG = AsVI(SDNom(debut=19))
