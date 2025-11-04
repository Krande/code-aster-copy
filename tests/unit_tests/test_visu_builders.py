import numpy as np
import pytest
import medcoupling as mc
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


def test_visu_builder_works_for_node_fields(tmp_path):
    med_mesh = get_mesh_1d()
    visu_cut = VisuCutBuilder(med_mesh)
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
    temp_med_file = tmp_path / "visu_1.med"
    visu_cut.write(temp_med_file)

    # Read med output
    actual_field: mc.MEDFileField1TS = mc.ReadFieldNode(
        fileName=str(temp_med_file),
        meshName=med_mesh.getName(),
        meshDimRelToMax=0,
        fieldName="DEPL",
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
    with pytest.raises(ValueError):
        visu_cut.add_field_on_nodes(
            field_name="DEPL",
            nodes=[1],
            values=np.array([(0.0,)]),
            components=("DX",),
            nume_ordre=1,
            instant=0.1,
        )
    with pytest.raises(ValueError):
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
    temp_med_file = tmp_path / "visu_2.med"
    visu_cut.write(temp_med_file)

    # Read med output
    actual_field: mc.MEDFileField1TS = mc.ReadFieldNode(
        fileName=str(temp_med_file),
        meshName=med_mesh.getName(),
        meshDimRelToMax=0,
        fieldName="DEPL",
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
