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

import unittest
from tempfile import TemporaryDirectory
from pathlib import Path

import numpy as np

from code_aster.Utilities.import_helper import medcoupling as mc
from code_aster.visu.visu_builders import VisuCutBuilder


def get_mesh_1d() -> mc.MEDFileUMesh:
    umesh_1d: mc.MEDFileUMesh = mc.MEDCouplingUMesh.New("cut_mesh_1d", 1)

    # fmt: off
    coords = [
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.3, 0.0, 0.0,
        1.5, 0.0, 0.0,

        0.0, 1.0, 0.0,
        0.0, 1.3, 0.0,
        0.0, 1.5, 0.0,

        0.0, 0.0, 1.0,
        0.0, 0.0, 1.3,
        0.0, 0.0, 1.5,
    ]
    # fmt: on

    coordsArr = mc.DataArrayDouble(
        coords, len(coords) // 3, 3
    )  # here coordsArr are declared to have 3 components, mesh will deduce that its spaceDim==3.

    umesh_1d.setCoords(coordsArr)
    umesh_1d.allocateCells(3 * 3)
    umesh_1d.insertNextCell(mc.NORM_SEG2, [0, 1])
    umesh_1d.insertNextCell(mc.NORM_SEG2, [1, 2])
    umesh_1d.insertNextCell(mc.NORM_SEG2, [2, 3])

    umesh_1d.insertNextCell(mc.NORM_SEG2, [0, 4])
    umesh_1d.insertNextCell(mc.NORM_SEG2, [4, 5])
    umesh_1d.insertNextCell(mc.NORM_SEG2, [5, 6])

    umesh_1d.insertNextCell(mc.NORM_SEG2, [0, 7])
    umesh_1d.insertNextCell(mc.NORM_SEG2, [7, 8])
    umesh_1d.insertNextCell(mc.NORM_SEG2, [8, 9])

    mesh_name = "cut_mesh"
    med_mesh = mc.MEDFileUMesh()
    med_mesh.setMeshAtLevel(0, umesh_1d)
    med_mesh.setName(mesh_name)
    return med_mesh


class TestVisuCutBuilder(unittest.TestCase):
    def test_visu_builder_works_for_node_fields(self):
        # enterContext will enter TemporaryDirectory context manager and exit it when test ends
        temp_dir = Path(self.enterContext(TemporaryDirectory()))

        med_mesh = get_mesh_1d()
        visu_cut = VisuCutBuilder(med_mesh, prefix_output_field_name="VISU")
        visu_cut.add_field_on_nodes(
            field_name="DEPL",
            nodes=[2, 1],
            values=np.array([(1.0, 0.0, 0.0), (0.0, 0.0, 1.0)]),
            components=("DX", "DY", "DZ"),
            nume_ordre=1,
            instant=0.1,
        )

        visu_cut.add_field_on_nodes(
            field_name="DEPL",
            nodes=[1],
            values=np.array([(0.0, 0.0, 2.0)]),
            components=("DX", "DY", "DZ"),
            nume_ordre=1,
            instant=0.1,
        )

        visu_cut.add_group(name="toto", ids=[1, 2], geo_type=VisuCutBuilder.NODE)
        visu_cut.add_group(name="toto2", ids=[2, 3], geo_type=VisuCutBuilder.NODE)
        visu_cut.add_group(name="tata", ids=[2, 3], geo_type=VisuCutBuilder.LINE)
        visu_cut.add_group(name="tata2", ids=[3, 4], geo_type=VisuCutBuilder.LINE)

        temp_med_file = temp_dir / "visu_1.med"
        visu_cut.write(temp_med_file)

        # Read med output
        actual_field: mc.MEDFileField1TS = mc.ReadFieldNode(
            fileName=str(temp_med_file),
            meshName=med_mesh.getName(),
            meshDimRelToMax=0,
            fieldName="VISU_DEPL",
            iteration=1,
            order=0,
        )

        actual_array = actual_field.getArray().toNumPyArray()
        np.testing.assert_equal(
            actual_array,
            np.array(
                [
                    (np.nan, np.nan, np.nan),
                    (0.0, 0.0, 2.0),
                    (1.0, 0.0, 0.0),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                ]
            ),
        )
        with self.assertRaises(ValueError):
            # incoherent components
            visu_cut.add_field_on_nodes(
                field_name="DEPL",
                nodes=[1],
                values=np.array([(0.0,)]),
                components=("DX",),
                nume_ordre=1,
                instant=0.1,
            )
        with self.assertRaises(ValueError):
            # incoherent instant
            visu_cut.add_field_on_nodes(
                field_name="DEPL",
                nodes=[1],
                values=np.array([(0.0, 0.0, 2.0)]),
                components=("DX", "DY", "DZ"),
                nume_ordre=1,
                instant=0.2,
            )

        # Adding another timestep
        visu_cut.add_field_on_nodes(
            field_name="DEPL",
            nodes=[0, 8],
            values=np.array([(1.0, 0.0, 0.0), (0.0, 0.0, 1.0)]),
            components=("DX", "DY", "DZ"),
            nume_ordre=2,
            instant=0.3,
        )

        # Add field without instant or nume_ordre
        visu_cut.add_field_on_nodes(
            field_name="AFIELD",
            nodes=[2, 1],
            values=np.array([(1.0, 0.0, 0.0), (0.0, 0.0, 1.0)]),
            components=("X", "Y", "Z"),
        )

        temp_med_file = temp_dir / "visu_2.med"
        visu_cut.write(temp_med_file)

        # Read med output
        actual_field: mc.MEDFileField1TS = mc.ReadFieldNode(
            fileName=str(temp_med_file),
            meshName=med_mesh.getName(),
            meshDimRelToMax=0,
            fieldName="VISU_DEPL",
            iteration=2,
            order=0,
        )

        actual_array = actual_field.getArray().toNumPyArray()
        np.testing.assert_equal(
            actual_array,
            np.array(
                [
                    (1.0, 0.0, 0.0),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (0.0, 0.0, 1.0),
                    (np.nan, np.nan, np.nan),
                ]
            ),
        )

        # Testing groups
        actual_mesh: mc.MEDFileUMesh = mc.MEDFileUMesh.New(str(temp_med_file))
        assert actual_mesh.getGroupArr(1, "toto").getValues() == [1, 2]
        assert actual_mesh.getGroupArr(0, "tata").getValues() == [2, 3]

        # Read med output (field 1 timesteps)
        actual_field: mc.MEDFileField1TS = mc.ReadField(str(temp_med_file), "VISU_AFIELD")
        actual_array = actual_field.getArray().toNumPyArray()
        np.testing.assert_equal(
            actual_array,
            np.array(
                [
                    (np.nan, np.nan, np.nan),
                    (0.0, 0.0, 1.0),
                    (1.0, 0.0, 0.0),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                    (np.nan, np.nan, np.nan),
                ]
            ),
        )
