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
from enum import IntFlag, auto

from code_aster import ConvergenceError, SolverError
from code_aster.Commands import DEFI_LIST_REEL
from code_aster.Solvers import TimeStepper
from code_aster.Solvers.base_features import BaseFeature

list0 = DEFI_LIST_REEL(VALE=0.0)
listr = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=10.0, PAS=1.0))


class UnittestOptions(IntFlag):
    """Enumeration of options for unittest."""

    System = auto()
    Storage = auto()
    State = auto()
    Unused = auto()
    Contact = auto()


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
        assert str(opt) == "UnittestOptions.Storage"
        opt |= FOP.System
        assert str(opt) == "UnittestOptions.Storage|System", str(opt)

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
        with self.assertRaisesRegex(TypeError, "not support.*UnittestOptions.Unused"):
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
        with self.assertRaisesRegex(TypeError, "not support.*UnittestOptions.Unused"):
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

    def test03_child(self):
        FOP = UnittestOptions

        class Main(BaseFeature):
            provide = FOP.System
            required_features = [FOP.Storage]

        class SubFeat(BaseFeature):
            provide = FOP.Storage
            required_features = [FOP.Contact]

        op = Main()
        sub = SubFeat()
        sub.use(object(), FOP.Contact)
        op.use(sub)
        self.assertEqual(len(op.get_features(FOP.Contact)), 0)
        self.assertEqual(len(sub.get_features(FOP.Contact)), 1)
        self.assertEqual(len(op.get_childs(FOP.Contact)), 1)


class TestTimeStepper(unittest.TestCase):
    """Check for internal methods."""

    def test00_init(self):
        stp = TimeStepper([0.0, 1.0, 2.0, 3.0])
        self.assertSequenceEqual(stp._times, [1.0, 2.0, 3.0])
        self.assertEqual(stp.size(), 3)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0], initial=1.0)
        self.assertSequenceEqual(stp._times, [2.0, 3.0])
        self.assertEqual(stp.size(), 2)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0])
        self.assertSequenceEqual(stp._times, [1.0, 2.0, 3.0])
        self.assertEqual(stp.size(), 3)
        self.assertEqual(stp.remaining(), stp.size())

        stp.setInitial(1.0)
        self.assertSequenceEqual(stp._times, [2.0, 3.0])
        self.assertEqual(stp.size(), 2)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0], initial=0.0, final=2.5)
        self.assertSequenceEqual(stp._times, [1.0, 2.0, 2.5])
        self.assertEqual(stp.size(), 3)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0])
        stp.setInitial(0.0)
        stp.setFinal(2.5)
        self.assertSequenceEqual(stp._times, [1.0, 2.0, 2.5])
        self.assertEqual(stp.size(), 3)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0])
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        self.assertAlmostEqual(stp.getInitial(), stp.getFinal())
        self.assertSequenceEqual(stp._times, [])
        self.assertEqual(stp.size(), 0)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0], initial=3.0)
        self.assertAlmostEqual(stp.getInitial(), 3.0)
        self.assertAlmostEqual(stp.getInitial(), stp.getFinal())
        self.assertSequenceEqual(stp._times, [])
        self.assertEqual(stp.size(), 0)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0], initial=2.0, final=2.5)
        self.assertSequenceEqual(stp._times, [2.5])
        self.assertEqual(stp.size(), 1)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([1.0, 2.0, 3.0], final=0.5)
        self.assertSequenceEqual(stp._times, [0.5])
        self.assertEqual(stp.size(), 1)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([], final=1.0)
        self.assertSequenceEqual(stp._times, [1.0])
        self.assertEqual(stp.size(), 1)
        self.assertEqual(stp.remaining(), stp.size())

        with self.assertRaisesRegex(ValueError, "ordered"):
            TimeStepper([0.0, 1.0, 3.0, 2.0])

    def test01_initial(self):
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitial(0.0)
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        self.assertEqual(stp.size(), 3)
        self.assertEqual(stp.remaining(), stp.size())
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)
        stp.completed()
        self.assertEqual(stp.remaining(), 2)
        self.assertFalse(stp.hasFinished())
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        stp.completed()
        self.assertEqual(stp.remaining(), 1)
        self.assertFalse(stp.hasFinished())
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)
        stp.completed()
        self.assertEqual(stp.remaining(), 0)
        self.assertTrue(stp.hasFinished())
        with self.assertRaisesRegex(IndexError, "no more timesteps"):
            stp.completed()

        eps = 1.0e-3
        stp = TimeStepper([0.25, 1.0, 2.0], eps)
        stp.setInitial(0.25 + eps * 0.999)
        self.assertLess(stp.getInitial() - 0.25, eps)
        self.assertSequenceEqual(stp._times, [1.0, 2.0])
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        self.assertEqual(stp.size(), 2)
        stp.completed()
        self.assertEqual(stp.remaining(), 1)
        self.assertFalse(stp.hasFinished())
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)
        self.assertSequenceEqual(stp._times, [1.0, 2.0])

        stp.setInitial(0.3)
        self.assertSequenceEqual(stp._times, [1.0, 2.0])
        self.assertAlmostEqual(stp.getInitial(), 0.3)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        self.assertEqual(stp.size(), 2)
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)

        stp.setInitial(0.1)
        self.assertSequenceEqual(stp._times, [1.0, 2.0])
        self.assertAlmostEqual(stp.getInitial(), 0.1)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        self.assertEqual(stp.size(), 2)
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)

        stp.setInitial(0.1 + eps * 0.999)
        self.assertSequenceEqual(stp._times, [1.0, 2.0])
        stp.completed()
        self.assertEqual(stp.size(), 2)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitial(1.0)
        self.assertAlmostEqual(stp.getInitial(), 1.0)
        self.assertEqual(stp.size(), 1)
        self.assertSequenceEqual(stp._times, [2.0])
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)

        stp.setInitial(2.5)
        self.assertAlmostEqual(stp.getInitial(), 2.5)
        self.assertEqual(stp.size(), 1)
        self.assertFalse(stp.hasFinished())
        stp.completed()
        self.assertTrue(stp.hasFinished())

    def test02_initial(self):
        stp = TimeStepper([2.0, 4.0, 6.0, 8.0, 10.0])
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        self.assertEqual(stp.size(), 5)
        stp.setInitial(2.0)
        self.assertAlmostEqual(stp.getInitial(), 2.0)
        self.assertEqual(stp.size(), 4)

        stp = TimeStepper([2.0, 4.0, 6.0, 8.0, 10.0])
        stp.setInitial(3.0)
        self.assertAlmostEqual(stp.getInitial(), 3.0)
        self.assertEqual(stp.size(), 4)

        stp = TimeStepper([2.0, 4.0, 6.0, 8.0, 10.0])
        stp.setInitial(1.0)
        self.assertAlmostEqual(stp.getInitial(), 1.0)
        self.assertEqual(stp.size(), 5)

        stp = TimeStepper([2.0, 4.0, 6.0, 8.0, 10.0])
        stp.setInitial(2.0)
        stp.setInitial(1.0)
        self.assertEqual(stp.size(), 4)

    def test03_final(self):
        eps = 1.0e-3
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0], eps)
        stp.setInitial(0.0)
        stp.setFinal(2.0 - eps * 0.999)
        self.assertLess(stp.getFinal() - 2.0, eps)
        self.assertEqual(stp.size(), 3)
        stp.setFinal(2.0 + eps * 0.999)
        self.assertLess(stp.getFinal() - 2.0, eps)
        self.assertEqual(stp.size(), 3)

        stp.setFinal(1.9)
        self.assertAlmostEqual(stp.getFinal(), 1.9)
        self.assertEqual(stp.size(), 3)
        self.assertSequenceEqual(stp._times, [0.25, 1.0, 1.9])
        for _ in range(3):
            stp.completed()
        with self.assertRaisesRegex(IndexError, "no more timesteps"):
            stp.completed()
        with self.assertRaisesRegex(IndexError, "no more timesteps"):
            step = stp.getCurrent()

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitial(0.0)
        stp.setFinal(2.5)
        self.assertEqual(stp.size(), 4)
        self.assertSequenceEqual(stp._times, [0.25, 1.0, 2.0, 2.5])
        for _ in range(3):
            stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.5)
        stp.completed()
        self.assertTrue(stp.hasFinished())

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitial(0.0)
        self.assertEqual(stp.size(), 3)
        stp.setFinal(0.8)
        self.assertEqual(stp.size(), 2)
        stp.setFinal()
        self.assertEqual(stp.size(), 2)

    def test04_basic(self):
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0, 3.0])
        self.assertEqual(stp.size(), 4)
        self.assertFalse(stp.hasFinished())

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)
        delta_t = stp.getIncrement()
        self.assertEqual(delta_t, 0.25)
        stp.completed()

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        delta_t = stp.getIncrement()
        self.assertAlmostEqual(delta_t, 0.75)
        stp.completed()

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)
        delta_t = stp.getIncrement()
        self.assertAlmostEqual(delta_t, 1.0)
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 3.0)
        delta_t = stp.getIncrement()
        self.assertAlmostEqual(delta_t, 1.0)
        stp.completed()

        other = stp.copy()
        self.assertEqual(other.size(), 4)
        self.assertFalse(other.hasFinished())
        step = other.getCurrent()
        self.assertAlmostEqual(step, 0.25)
        delta_t = other.getIncrement()
        self.assertEqual(delta_t, 0.25)

        self.assertTrue(stp.hasFinished())
        with self.assertRaisesRegex(IndexError, "no more timesteps"):
            stp.completed()

    def test05_cmp(self):
        stp = TimeStepper([1.0], epsilon=0.01)
        self.assertEqual(stp.cmp(1.0, 1.0001), 0)
        self.assertEqual(stp.cmp(1.0, 1.1), -1)
        self.assertEqual(stp.cmp(1.2, 1.1), 1)

    def test07_meca_statique(self):
        stp = TimeStepper([0.0], initial=None)
        self.assertEqual(stp.size(), 1)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.0)
        self.assertIsNone(stp.getPrevious())
        self.assertIsNone(stp.getIncrement())
        stp.completed()
        self.assertTrue(stp.hasFinished())

    def test08_ther_lineaire(self):
        stp = TimeStepper.from_keywords(LIST_INST=list0, INST_INIT=None, PRECISION=1.0e-6)
        self.assertEqual(stp.size(), 1)
        self.assertAlmostEqual(stp.getInitial(), None)
        self.assertAlmostEqual(stp.getFinal(), 0.0)

        stp = TimeStepper.from_keywords(LIST_INST=listr, INST_FIN=5.0, PRECISION=1.0e-6)
        self.assertEqual(stp.size(), 5)
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        self.assertAlmostEqual(stp.getFinal(), 5.0)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        stp.completed()
        self.assertEqual(stp.remaining(), 4)

        stp = TimeStepper([-0.5, 5.5, 11.5], initial=-0.5)
        self.assertEqual(stp.size(), 2)
        self.assertAlmostEqual(stp.getInitial(), -0.5)
        self.assertAlmostEqual(stp.getFinal(), 11.5)

        stp = TimeStepper.from_keywords(
            LIST_INST=listr, NUME_INST_INIT=0, NUME_INST_FIN=1, PRECISION=1.0e-6
        )
        self.assertEqual(stp.size(), 1)
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        self.assertAlmostEqual(stp.getFinal(), 1.0)

    def test20_event(self):
        stp = TimeStepper([1.0, 1.1, 2.0])
        stp.register_event(TimeStepper.Interrupt(TimeStepper.Error()))
        with self.assertRaisesRegex(ConvergenceError, "MESSAGEID"):
            stp.failed(ConvergenceError("MESSAGEID"))

        stp = TimeStepper([1.0, 1.1, 2.0])
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 1.1)

    def test21_split(self):
        stp = TimeStepper([0.0, 1.0, 1.1, 2.0, 3.0])
        self.assertEqual(stp.size(), 4)
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        stp.register_event(
            TimeStepper.Split(TimeStepper.Error(), nbSteps=2, maxLevel=3, minStep=0.05)
        )
        # print("\n+ split #1")
        stp.failed(ConvergenceError("MESSAGEID"))
        # [0.5, 1.0, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 5)
        self.assertAlmostEqual(stp.getCurrent(), 0.5)
        # print("+ split #2")
        stp.failed(ConvergenceError("MESSAGEID"))
        # [0.25, 0.5, 1.0, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 6)
        self.assertAlmostEqual(stp.getCurrent(), 0.25)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 0.5)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        # print("+ split #1")
        stp.failed(ConvergenceError("MESSAGEID"))
        # [0.25, 0.5, 0.75, 1.0, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 7)
        self.assertAlmostEqual(stp.getCurrent(), 0.75)
        # print("+ split #2")
        stp.failed(ConvergenceError("MESSAGEID"))
        # [0.25, 0.5, 0.625, 0.75, 1.0, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 8)
        self.assertAlmostEqual(stp.getCurrent(), 0.625)
        # print("+ split #3")
        stp.failed(ConvergenceError("MESSAGEID"))
        # [0.25, 0.5, 0.5625, 0.625, 0.75, 1.0, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 9)
        self.assertAlmostEqual(stp.getCurrent(), 0.5625)
        # print("+ split #4")
        with self.assertRaisesRegex(SolverError, "max.*subdivision"):
            stp.failed(ConvergenceError("MESSAGEID"))
        stp.completed()
        stp.completed()
        stp.completed()
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 1.1)
        # print("+ split #1")
        stp.failed(ConvergenceError("MESSAGEID"))
        # [0.25, 0.5, 0.5625, 0.625, 0.75, 1.0, 1.05, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 10)
        self.assertAlmostEqual(stp.getCurrent(), 1.05)
        # print("+ split #2")
        with self.assertRaisesRegex(SolverError, "trop petit"):
            stp.failed(ConvergenceError("MESSAGEID"))


if __name__ == "__main__":
    unittest.main()
