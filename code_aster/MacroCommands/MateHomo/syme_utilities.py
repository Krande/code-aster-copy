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

from ...Utilities import medcoupling as medc, no_new_attributes
from .mate_homo_mesh import rebuild_with_groups
from . import MESH_TOL


class SymmetryManager:
    """
    Handle the symmetric projection of homogeneisation corrector fields.

    Input data and mesh are output of CALC_MATE_HOMO.
    Mesh is either 1/8 or 1/4 of the structure defined with positive coordinates.

    """

    _axis = {"X": [-1, 0, 0], "Y": [0, -1, 0], "Z": [0, 0, -1]}
    _dir = _x0 = _y0 = _z0 = None
    __setattr__ = no_new_attributes(object.__setattr__)

    @property
    def axis(self):
        """Holds the projection axis"""
        assert self._dir is not None
        return self._axis[self._dir]

    @property
    def point(self):
        """Holds the symmetry point"""
        return (self._x0, self._y0, self._z0)

    def __init__(self, point, axis):
        """
        Init function.

        Arguments
        ---------
            point (list[float]) : Symmetry point (x0, y0, z0).
            axis (str) : Symmetry axis (X|Y|Z).
        """
        assert axis in self._axis
        self._dir = axis
        self._x0, self._y0, self._z0 = point

    def sign(self, fieldname):
        """Return the sign vector.

        Corrector fields might be symmetric or anti-symmetric.
        This condition affect the component which sign changes during symmetric projection.
        A non-zero returned value implies a sign change on the projected component value.

        Arguments
        ---------
            fieldname (str): The corrector field name as produced by CALC_MATE_HOMO

        Returns
        -------
            sign [int, int, int]): Sign vector.
        """

        assert self._dir is not None
        if any(x in fieldname for x in ["11", "22", "33", "DI", "PI"]):
            ret = self._axis[self._dir]

        elif "12" in fieldname:
            if self._dir == "X":
                ret = self._axis["Y"]
            elif self._dir == "Y":
                ret = self._axis["X"]
            elif self._dir == "Z":
                ret = self._axis["Z"]
            else:
                assert False

        elif "31" in fieldname:
            if self._dir == "X":
                ret = self._axis["Z"]
            elif self._dir == "Y":
                ret = self._axis["Y"]
            elif self._dir == "Z":
                ret = self._axis["X"]
            else:
                assert False

        elif "23" in fieldname:
            if self._dir == "X":
                ret = self._axis["X"]
            elif self._dir == "Y":
                ret = self._axis["Z"]
            elif self._dir == "Z":
                ret = self._axis["Y"]
            else:
                assert False

        else:
            assert False

        return ret

    def build_symmetry_mesh(self, meshin):
        """Return the full mesh.

        Arguments
        ---------
            meshin (MEDFileUMesh): The symmetric mesh of corrector fields.

        Returns
        -------
            meshout (MEDFileUMesh): The full mesh of corrector fields.
            keep_ids (DataArrayInt) : Node ids of non-merged nodes.
            merge_ids (DataArrayInt) : Node ids of merged nodes.
        """

        if not isinstance(meshin, medc.MEDFileUMesh):
            raise RuntimeError("Input mesh is not MEDFileUMesh as expected !")
        if meshin.getMeshDimension() != 3:
            raise RuntimeError("Works only for 3D mesh")
        if (
            not medc.DataArrayDouble(self.axis, 1, 3)
            .magnitude()
            .isEqual(medc.DataArrayDouble([1.0]), MESH_TOL)
        ):
            raise RuntimeError("Input normal has to have magnitude equal to one with eps")
        nb_nodes_meshin = meshin.getNumberOfNodes()
        coors_sym = meshin.getCoords().symmetry3DPlane(self.point, self.axis)
        coords_merged_with_duplicates = medc.DataArrayDouble.Aggregate(
            [meshin.getCoords(), coors_sym]
        )
        common_nodes_id, common_nodes_shift = coords_merged_with_duplicates.findCommonTuples(
            MESH_TOL
        )
        if not common_nodes_shift.deltaShiftIndex().findIdsNotEqual(2).empty():
            msg = "There are groups of more than 2 points considered as merged."
            raise RuntimeError(msg)
        common_nodes_id.rearrange(2)
        if not (common_nodes_id[:, 0] + nb_nodes_meshin).isEqual(common_nodes_id[:, 1]):
            msg = "With epsilon specified impossible to detect cleanly pts to be merged"
            raise RuntimeError(msg)
        common_nodes_id.rearrange(1)
        o2n, newNbnodes = medc.DataArrayInt.ConvertIndexArrayToO2N(
            2 * nb_nodes_meshin, common_nodes_id, common_nodes_shift
        )
        common_nodes_id.rearrange(2)
        meshin_copy = meshin.deepCopy()
        mesh_fused = medc.MEDFileUMesh.Aggregate([meshin, meshin_copy])
        point_ids_to_merge = common_nodes_id[:, 0]

        # deal with coordinates in mesh_fused
        point_ids_to_keep_sym = point_ids_to_merge.buildComplement(nb_nodes_meshin)

        point_ids_to_keep_merged = medc.DataArrayInt.Aggregate(
            [
                medc.DataArrayInt.Range(0, nb_nodes_meshin, 1),
                point_ids_to_keep_sym + nb_nodes_meshin,
            ]
        )
        mesh_fused.forceComputationOfParts()
        for lev in mesh_fused.getNonEmptyLevels():
            for tmp in mesh_fused.getDirectUndergroundSingleGeoTypeMeshes(lev):
                tmp.renumberNodesInConn(o2n)

        # indicate that connectivity of parts have been updated
        mesh_fused.declarePartsUpdated()
        mesh_fused.setCoords(coords_merged_with_duplicates[point_ids_to_keep_merged])
        mesh_fused.setName(meshin.getName())

        l0grp = [mesh_fused.getGroupArr(0, grp) for grp in mesh_fused.getGroupsOnSpecifiedLev(0)]
        volume_ver, dirthick, mesh_merged = rebuild_with_groups(mesh_fused[0], l0grp)

        return mesh_merged, point_ids_to_keep_sym, point_ids_to_merge

    def build_symmetry(self, resuin, meshfull=None, keep_ids=None, merge_ids=None):
        """Return a new result containing the symmetric result merged with the current one.

        Data is mirrored along the stored `axis`. Only nodal fields are taken into account.

        Arguments
        ---------
            resuin (MEDFileData): Input corrector fields.
            meshfull (MEDFileUMesh, optional): The full mesh of corrector fields.
            keep_ids (DataArrayInt, optional) : Node ids of non-merged nodes.
            merge_ids (DataArrayInt, optional) : Node ids of merged nodes.

        Returns
        -------
            resuout (MEDFileData): Output corrector fields.
        """

        meshesin = resuin.getMeshes()
        assert len(meshesin) == 1
        meshin = meshesin[0]

        fmtsin = resuin.getFields()

        result = medc.MEDFileData()
        meshes = medc.MEDFileMeshes()
        fields = medc.MEDFileFields()
        result.setMeshes(meshes)
        result.setFields(fields)

        if meshfull is None:
            mesh_merged, point_ids_to_keep_sym, point_ids_to_merge = self.build_symmetry_mesh(
                meshin
            )
        else:
            mesh_merged, point_ids_to_keep_sym, point_ids_to_merge = meshfull, keep_ids, merge_ids

        meshes.pushMesh(mesh_merged)
        m0_merged = mesh_merged[0]

        for fmts in fmtsin:
            fmts_merged = type(fmts).New()
            for f1ts in fmts:
                field_type_list = f1ts.getTypesOfFieldAvailable()
                if len(field_type_list) != 1:
                    msg = "Fields with multiple space discr not managed !"
                    raise RuntimeError(msg)
                if len(f1ts.getPflsReallyUsed()) != 0:
                    msg = "Profiles are not managed now !"
                    raise RuntimeError(msg)
                field_type = field_type_list[0]
                with f1ts:
                    curr_field = f1ts.field(meshin)
                    vect_proj = self.sign(curr_field.getName())
                    values_field = curr_field.getArray()
                    f1ts_merged = type(f1ts).New()
                    if field_type == medc.ON_NODES:
                        values_sym_part = values_field[point_ids_to_keep_sym]
                        if values_sym_part.getNumberOfComponents() == 3:
                            normvec1 = medc.DataArrayDouble(len(values_sym_part), 3)
                            normvec1[:, :] = vect_proj
                            values_sym_part = +1 * (
                                values_sym_part
                                - 2 * medc.DataArrayDouble.Dot(values_sym_part, normvec1) * normvec1
                            )

                            values_common = values_field[point_ids_to_merge]
                            normvec2 = medc.DataArrayDouble(len(values_common), 3)
                            normvec2[:, :] = vect_proj
                            vtest = medc.DataArrayDouble.Dot(values_common, normvec2)
                            vtest.abs()
                            check = vtest.findIdsGreaterThan(MESH_TOL).empty()
                            if not check:
                                msg = "Presence of values in node field of pts to be merged that are not equal to 0 along input planeVector"
                                raise RuntimeError(msg)

                        values_merged = medc.DataArrayDouble.Aggregate(
                            [values_field, values_sym_part]
                        )
                        merged_field = type(curr_field).New(curr_field.getTypeOfField())
                        merged_field.setArray(values_merged)
                        merged_field.copyAllTinyAttrFrom(curr_field)
                        merged_field.setMesh(m0_merged)

                    f1ts_merged.setFieldNoProfileSBT(merged_field)
                fmts_merged.pushBackTimeStep(f1ts_merged)
            fields.pushField(fmts_merged)
        return result
