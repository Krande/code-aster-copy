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


class sd_proj_modal(AsBase):
    # ------------------------------------
    nomj = SDNom(fin=19)

    MASG = AsVR()

    AGEM = AsVR()
    AGEP = AsVR()
    VGEM = AsVR()
    VGEP = AsVR()
    DGEM = AsVR()
    DGEP = AsVR()
    RIGG = AsVR()
    AMOG = AsVR()
    BASM = AsVR()
    FONG = Facultatif(AsVK24())
    FORG = Facultatif(AsVR())
