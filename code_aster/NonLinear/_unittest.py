# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
from base_features import BaseFeature, BaseFeaturesOptions


class UnittestOptions(BaseFeaturesOptions):
    """Enumeration of options for unittest."""

    System = 0x001
    Storage = 0x002
    State = 0x004
    Unused = 0x008
    Contact = 0x010


class SystemDefinition(BaseFeature):
    """Object that stores the definition of the system to be solved."""

    provide = UnittestOptions.System


class Storage(BaseFeature):
    """Object that provides the archiving feature."""

    provide = UnittestOptions.Storage


class Anonymous:
    """Non-Feature object"""


class BasicTest(unittest.TestCase):
    """Check for internal methods."""

    def test01_basics(self):
        FOP = UnittestOptions
        opt = FOP.Storage
        assert FOP.name(opt) == "Storage"
        opt |= FOP.System
        assert FOP.name(opt) == "System|Storage", FOP.name(opt)

        syst = SystemDefinition()
        stor = Storage()
        # object that provides services for System + Storage
        syssto = SystemDefinition()
        anonym = Anonymous()

        class TestFeature(BaseFeature):
            required_features = [FOP.System, FOP.State]
            optional_features = [FOP.Storage]

        op = TestFeature()
        op.use(syst)
        op.use(stor)
        op.use(syssto, FOP.Storage)
        op.use(anonym, FOP.State)
        with self.assertRaises(TypeError):
            op.use(anonym, FOP.Unused)

        self.assertCountEqual(op.get_features(FOP.System), [syst, syssto])
        self.assertCountEqual(op.get_features(FOP.Storage), [stor, syssto])
        self.assertCountEqual(op.get_features(FOP.System | FOP.Storage), [syssto])

        op.discard(syssto)
        self.assertCountEqual(op.get_features(FOP.System), [syst])
        self.assertCountEqual(op.get_features(FOP.Storage), [stor])
        self.assertCountEqual(op.get_features(FOP.System | FOP.Storage), [])

    def test02_nonl(self):
        FOP = UnittestOptions

        class TestFeature(BaseFeature):
            provide = FOP.System
            required_features = [FOP.System, FOP.Storage]
            optional_features = [FOP.Contact]

        op = TestFeature()
        op.use(object(), FOP.System)
        with self.assertRaises(TypeError):
            op.use(object(), FOP.Unused)

        self.assertEqual(len(op.required_features), 2)
        self.assertEqual(len(op.optional_features), 1)

        class InheritedFeature(TestFeature):
            provide = FOP.System
            optional_features = TestFeature.optional_features + [FOP.State]

        inherit = InheritedFeature()
        self.assertEqual(len(inherit.required_features), 2)
        self.assertEqual(len(inherit.optional_features), 2)


if __name__ == "__main__":
    unittest.main()
