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
from time_stepper import TimeStepper


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
        self.assertEqual(len(op._use), 0)
        op.use(syst)
        op.use(stor)
        op.use(syssto, FOP.Storage)
        with self.assertRaises(TypeError):
            op.use(anonym, FOP.Unused)
        op.use(anonym, FOP.State)
        self.assertEqual(len(op._use), 4)
        op.use(syst)
        self.assertEqual(len(op._use), 4)

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
        self.assertSequenceEqual(op.undefined(), [(FOP.Storage, True), (FOP.Contact, False)])

        self.assertEqual(len(op.required_features), 2)
        self.assertEqual(len(op.optional_features), 1)

        class InheritedFeature(TestFeature):
            provide = FOP.System
            optional_features = TestFeature.optional_features + [FOP.State]

        inherit = InheritedFeature()
        self.assertEqual(len(inherit.required_features), 2)
        self.assertEqual(len(inherit.optional_features), 2)


class TestTimeStepper(unittest.TestCase):
    """Check for internal methods."""

    def test_initial(self):
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitialStep(0.0)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.0)
        stp.start()
        self.assertEqual(stp.size(), 4)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)

        eps = 1.0e-3
        stp = TimeStepper([0.25, 1.0, 2.0], eps)
        stp.setInitialStep(0.0)
        self.assertEqual(stp.size(), 4)
        stp.start()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)

        stp.setInitialStep(0.3)
        self.assertEqual(stp.size(), 3)
        stp.start()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)

        stp.setInitialStep(0.1)
        stp.start()
        self.assertEqual(stp.size(), 5)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)

        stp.setInitialStep(0.1 + eps * 0.999)
        stp.start()
        self.assertEqual(stp.size(), 5)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitialStep(1.0)
        self.assertEqual(stp.size(), 2)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)

        stp.setInitialStep(2.5)
        self.assertEqual(stp.size(), 1)
        self.assertFalse(stp.hasFinished())
        stp.start()
        self.assertTrue(stp.hasFinished())

    def test_final(self):
        eps = 1.0e-3
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0], eps)
        stp.setFinalStep(2.0 - eps * 0.999)
        self.assertEqual(stp.size(), 4)
        stp.setFinalStep(2.0 + eps * 0.999)
        self.assertEqual(stp.size(), 4)

        stp.setFinalStep(1.9)
        self.assertEqual(stp.size(), 4)
        for _ in range(3):
            stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.9)

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setFinalStep(2.5)
        self.assertEqual(stp.size(), 5)
        for _ in range(4):
            stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.5)

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        self.assertEqual(stp.size(), 4)
        stp.setFinalStep(0.8)
        self.assertEqual(stp.size(), 3)
        stp.setFinalStep()
        self.assertEqual(stp.size(), 5)

    def test_basic(self):
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0, 3.0])
        self.assertFalse(stp.hasFinished())

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.0)
        delta_t = stp.getIncrement()
        self.assertEqual(delta_t, stp.null_increment)

        stp.setInitialStep(0.0)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.0)
        stp.start()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)
        delta_t = stp.getIncrement()
        self.assertAlmostEqual(delta_t, 0.25)
        stp.completed()

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        delta_t = stp.getIncrement()
        self.assertAlmostEqual(delta_t, 0.75)
        stp.completed()

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)

        with self.assertRaises(ValueError):
            stp.raiseError(ValueError())

        self.assertFalse(stp.hasFinished())
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 3.0)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 3.0)
        self.assertFalse(stp.hasFinished())

    def test_restart(self):
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0, 3.0])
        self.assertEqual(stp.size(), 5)
        stp.setFinalStep(1.0)
        self.assertEqual(stp.size(), 3)

        stp.start()
        self.assertAlmostEqual(stp.getCurrent(), 0.25)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        stp.completed()
        self.assertTrue(stp.hasFinished())

        stp.setInitialStep(1.0)
        self.assertEqual(stp.size(), 1)
        stp.setFinalStep(3.0)
        self.assertEqual(stp.size(), 3)
        stp.start()
        self.assertAlmostEqual(stp.getCurrent(), 2.0)

    def test_split(self):
        stp = TimeStepper([0.0, 1.0, 2.0, 3.0])
        self.assertEqual(stp.size(), 4)
        stp.start()
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        stp.split(2)
        self.assertEqual(stp.size(), 5)
        self.assertAlmostEqual(stp.getCurrent(), 0.5)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        stp.split(5)
        self.assertEqual(stp.size(), 9)
        self.assertAlmostEqual(stp.getCurrent(), 0.6)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 0.7)


if __name__ == "__main__":
    unittest.main()
