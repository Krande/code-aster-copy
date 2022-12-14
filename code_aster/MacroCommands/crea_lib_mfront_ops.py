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

# person_in_charge: mathieu.courtois@edf.fr

import os.path as osp
import shutil
from subprocess import call

from ..Helpers.LogicalUnit import LogicalUnitFile
from ..Messages import UTMESS
from ..Utilities import ExecutionParameter


def crea_lib_mfront_ops(self, UNITE_MFRONT, UNITE_LIBRAIRIE, DEBUG, **args):
    """Compiler une loi de comportement MFront"""

    infile = LogicalUnitFile.filename_from_unit(UNITE_MFRONT)
    outlib = LogicalUnitFile.filename_from_unit(UNITE_LIBRAIRIE)

    cmd = [ExecutionParameter().get_option("prog:mfront"), "--build", "--interface=aster"]
    if DEBUG == "OUI":
        cmd.append("--debug")
        # cmd.append("--@AsterGenerateMTestFileOnFailure=true")
    cmd.append(infile)

    try:
        call(cmd)
        if not osp.exists("src/libAsterBehaviour.so"):
            UTMESS("F", "MFRONT_4", valk="libAsterBehaviour.so")
        shutil.copyfile("src/libAsterBehaviour.so", outlib)
    finally:
        for dname in ("src", "include"):
            if osp.isdir(dname):
                shutil.rmtree(dname)

    return
