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

#
# person_in_charge: mathieu.courtois at edf.fr

"""
This module provides unittests for the cataelem package
"""

import unittest

from cataelem.Tools.base_objects import Element


class TestCataElem(unittest.TestCase):

    """Checked classes: LocatedComponents, Element"""

    def setUp(self):
        """Setup objects as CataElem does"""
        from cataelem.Commons.attributes import MODTHM
        from cataelem.Commons.mesh_types import SEG2
        from cataelem.Commons.physical_quantities import DEPL_R, SIEF_R

        MODTHM.setName("MODTHM")
        elr = ("SE2", "CABPOU", "THCOSE2")
        for name, elrefe in zip(elr, SEG2.getElrefe()):
            elrefe.setName(name)
        SEG2.setName("SEG2")
        DEPL_R.setName("DEPL_R")
        SIEF_R.setName("SIEF_R")
        self._options()

    def _options(self):
        """Import only few options for testcases, replace Options/options.py"""
        # from cataelem.Options.options import OP
        # self.OP = OP
        from cataelem.Tools.base_objects import AbstractEntityStore, Option
        from cataelem.Tools.base_objects import InputParameter, OutputParameter

        class OptionStore(AbstractEntityStore):
            """Helper class to give access to all options to elements"""

            entityType = Option
            subTypes = (InputParameter, OutputParameter)

        self.OP = OptionStore("Options", only_mods=["sieq_elga", "existe_ddl"])
        for name, obj in list(self.OP.getDict().items()):
            obj.setName(name)

    def test01_attribute(self):
        """check Attribute object"""
        from cataelem.Commons.attributes import MODTHM

        MODTHM.comment = "comment string"
        self.assertGreater(MODTHM.idx, 0)
        self.assertEqual(MODTHM.name, "MODTHM")
        self.assertEqual(MODTHM.comment, "comment string")
        self.assertSequenceEqual(
            MODTHM.value,
            (
                "H",
                "HH",
                "HH2",
                "HH2M",
                "HHM",
                "HM",
                "SUSHI",
                "THH",
                "THH2",
                "THH2M",
                "THHM",
                "THM",
                "THV",
            ),
        )
        self.assertTrue(MODTHM.isValid("HHM"))
        self.assertFalse(MODTHM.isValid("H3M"))

    def test02_mesh(self):
        """check MeshType object"""
        from cataelem.Commons.mesh_types import SEG2

        self.assertEqual(SEG2.name, "SEG2")
        self.assertEqual(SEG2.nbNodes, 2)
        self.assertEqual(SEG2.dim, 1)
        self.assertEqual(SEG2.code, "SE2")

    def test03_physical_quantity(self):
        """check PhysicalQuantity object"""
        from cataelem.Commons.physical_quantities import DEPL_R

        self.assertEqual(DEPL_R.name, "DEPL_R")
        self.assertIn("DX", DEPL_R.components)
        self.assertEqual(DEPL_R.type, "R")
        self.assertTrue(DEPL_R.hasComponent("TEMP"))
        self.assertFalse(DEPL_R.hasComponent("VMIS"))

    def test04_located_components(self):
        """check LocatedComponents object"""

    def test05_parameters(self):
        """check InputParameter & OutputParameter objects"""

    def test06_option(self):
        """check Option object"""

    def test07_element(self):
        """check Element object"""
        from cataelem.Tools.base_objects import LocatedComponents
        from cataelem.Tools.base_objects import Element
        from cataelem.Commons.mesh_types import SEG2
        from cataelem.Commons.physical_quantities import DEPL_R, SIEF_R
        from cataelem.Options.options import OP

        # build of OP should have named their own parameters (not shared)
        self.assertEqual(OP.SIEQ_ELGA.PCONTEQ.name, "PCONTEQ")
        self.assertEqual(OP.EXISTE_DDL.PDEPL_R.name, "PDEPL_R")

        DDL_MECA = LocatedComponents(phys=DEPL_R, type="ELNO", components=("DX", "DY", "LAGS_C"))
        DDL_MECA.setName("DDL_MECA")
        self.assertEqual(DDL_MECA.name, "DDL_MECA")
        ECONTPG = LocatedComponents(
            phys=SIEF_R,
            type="ELGA",
            location="RIGI",
            components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
        )
        ECONTPG.setName("ECONTPG")
        self.assertEqual(ECONTPG.name, "ECONTPG")

        class ELE4TEST(Element):
            """This is a test element"""

            meshType = SEG2
            calculs = (
                OP.SIEQ_ELGA(te=50, para_out=((OP.SIEQ_ELGA.PCONTEQ, ECONTPG),)),
                OP.EXISTE_DDL(te=99, para_out=((OP.EXISTE_DDL.PDEPL_R, DDL_MECA),)),
            )

        elt = ELE4TEST()
        elt.setName("ELE4TEST")
        self.assertEqual(elt.name, "ELE4TEST")
        opts = elt.getCalculs()
        self.assertEqual(len(opts), 2)
        self.assertEqual(opts[1][0], "EXISTE_DDL")
        used = elt.usedLocatedComponents()
        self.assertEqual(len(used), 2)
        self.assertIn(DDL_MECA, used)
        self.assertSequenceEqual(used[1].components, ("DX", "DY", "LAGS_C"))
        elt.changeComponents("DDL_MECA", ("DX", "DZ", "TEMP"))
        used = elt.usedLocatedComponents()
        self.assertEqual(len(used), 2)
        self.assertNotIn(DDL_MECA, used)
        new = used[1]
        self.assertEqual(new.name, "DDL_MECA")
        self.assertIn("DX", new.components)
        self.assertIn("DZ", new.components)
        self.assertNotIn("LAGS_C", new.components)


if __name__ == "__main__":
    unittest.main()
