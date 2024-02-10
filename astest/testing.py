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

import unittest
import sys

from code_aster import CA


def _test_module(module):
    print(f"\n\n+++ testing {module}...\n", flush=True)
    result = unittest.main(argv=["comm"], module=module, exit=False, verbosity=2).result
    # to flush printings from unittest
    sys.stdout.flush()
    sys.stderr.flush()
    assert result.wasSuccessful()


_test_module("code_aster.Helpers.syntax_repr")
