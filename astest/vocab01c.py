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

import os
import sys
import unittest

from code_aster.Cata import Commands as CATA
from code_aster.Commands import *
from code_aster import CA
from code_aster.Helpers.syntax_repr import ALT, OPT, REQ, XOR, repr_command

CA.init("--test")
test = CA.TestCase()


def _test_module(module):
    print(f"\n\n+++ testing {module}...\n", flush=True)
    result = unittest.main(argv=["comm"], module=module, exit=False, verbosity=2).result
    # to flush printings from unittest
    sys.stdout.flush()
    sys.stderr.flush()
    assert result.wasSuccessful()


_test_module("code_aster.Helpers.syntax_repr")

# to export all code-blocks into a directory:
# loop_on_commands("/local00/tmp/syntax")

text = repr_command(CATA.AFFE_MATERIAU, show=False)

test.assertTrue(f"{REQ} {ALT} MAILLAGE =" in text, msg="MAILLAGE")
test.assertTrue(f"  {ALT} MODELE =" in text, msg="MODELE")
test.assertTrue(f"{REQ} {XOR} AFFE = _F(" in text, msg="AFFE")
test.assertTrue(f'{REQ} {XOR} TOUT = "OUI"' in text, msg="TOUT AFFE/AFFE_COMPOR")
test.assertTrue(f'{OPT} {XOR} TOUT = "OUI"' in text, msg="TOUT AFFE_VARC")

cmd = os.environ.get("CMD")
if cmd:
    cata = getattr(CATA, cmd, None)
    if not cata:
        raise ValueError(f"command not found in catalogs: {cmd}")
    sep = "+" * 80
    content = [sep, sep, "Syntaxe", "=======", "", repr_command(cata, show=False), sep, sep]
    print(os.linesep.join(content))

test.printSummary()
CA.close()
