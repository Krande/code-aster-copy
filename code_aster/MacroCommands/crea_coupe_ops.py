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
import sys

import aster
from libaster import AsterError, projectionAlongDirection

from ..Messages import UTMESS

from ..Commands import LIRE_TABLE, IMPR_TABLE, CREA_TABLE
from ..Objects.table_py import Table

import numpy as np


def crea_coupe_ops(self, COUPE, MAILLAGE, NOM_AUTO, PREFIX=None, PREM_NUME=None, PAS_NUME=None):
    """Execute the command.

    Arguments:
        COUPE (sd_table)        : input table
        MAILLAGE(sd_maillage)   : input mesh
        NOM_AUTO(bool)          : automatic naming of paths
        PREFIX(str)             : common prefix of the path names
        PREM_NUME(int)          : first id of the name of the path
        PAS_NUME(int)           : progression in the numerotation between two successive paths

    Returns:
        sd_table which descibes paths that begin/end on the skin of the mesh
    """

    table_coupes_py = COUPE.EXTR_TABLE()
    table_coupes = TableCoupes(table_coupes_py)
    table_coupes.change_names(
        auto_name=NOM_AUTO, prefix=PREFIX, first_num=PREM_NUME, progression=PAS_NUME
    )
    table_coupes.check_integrity()
    table_coupes._check_group_maill(MAILLAGE)
    table_coupes.update_position_on_structure_skin(MAILLAGE)
    dic_table = table_coupes.dict_CREA_TABLE()
    table_coupes_aster = CREA_TABLE(**dic_table)

    return table_coupes_aster


class TableCoupes(Table):
    ALL_KEYWORDS = [
        "COOR_ORIG_X",
        "COOR_ORIG_Y",
        "COOR_ORIG_Z",
        "COOR_EXTR_X",
        "COOR_EXTR_Y",
        "COOR_EXTR_Z",
        "COOR_TROIS_X",
        "COOR_TROIS_Y",
        "COOR_TROIS_Z",
        "NOM",
        "GROUP",
        "GROUP_MA_ORIG",
        "GROUP_MA_EXTR",
        "NB_POINTS",
    ]
    ATTEMPTED_TYPES = ["R", "R", "R", "R", "R", "R", "R", "R", "R", "K24", "K24", "K24", "K24", "I"]
    DIMENSION = 3
    FIRST_POINT_ID = 0
    SECOND_POINT_ID = 3
    THIRD_POINT_ID = 6
    NUMBER_OF_NODES_ID = 13
    MAIL_IN_ID = 11
    MAIL_OUT_ID = 12
    NAMES_ID = 9
    GROUPS_ID = 10

    def __init__(self, table=None, rows=None, para=None, typ=None, titr=None, nom=None):
        if table is not None:
            super(TableCoupes, self).__init__(
                rows=table.rows, para=table.para, typ=table.type, titr=table.titr, nom=table.nom
            )
        elif (
            rows is not None
            and para is not None
            and typ is not None
            and titr is not None
            and nom is not None
        ):
            super(TableCoupes, self).__init__(rows=rows, para=para, typ=typ, titr=titr, nom=nom)
        else:
            raise SystemExit("You have to specify at least a table or a set of table parameters")

    def check_integrity(self):
        """Check if the input table has the right format"""
        self.check_keywords()
        self.check_types()
        self.check_path_name()

    def check_keywords(self):
        """Check if the input table has the right column names"""
        for i, key in enumerate(self.para):
            if key.lower() != TableCoupes.ALL_KEYWORDS[i].lower():
                UTMESS(
                    "F",
                    "COUPE_1",
                    vali=(i + 1),
                    valk=[TableCoupes.ALL_KEYWORDS[i], key, ", ".join(TableCoupes.ALL_KEYWORDS)],
                )

    def check_types(self):
        """Check if each column of the input table has the right type of data"""
        dico_types = {"R": "décimal", "I": "entier", "K24": "chaîne de caractères"}
        for i, typ in enumerate(self.type):
            if typ != TableCoupes.ATTEMPTED_TYPES[i]:
                UTMESS(
                    "F",
                    "COUPE_5",
                    vali=(i + 1),
                    valk=[dico_types[TableCoupes.ATTEMPTED_TYPES[i]], dico_types[typ]],
                )

    def check_path_name(self):
        """Check if a name is used several times in the table"""
        already_set_name = []
        for line in self.rows:
            name = line[TableCoupes.ALL_KEYWORDS[TableCoupes.NAMES_ID]]
            if name in already_set_name:
                UTMESS("F", "COUPE_7", valk=name)
            else:
                already_set_name.append(name)

    def change_names(self, auto_name="NON", prefix=None, first_num=None, progression=None):
        """
        change the names of the section accordingly to the template:
            PREFIX_"SECTION_ID",
        where SECTION_ID is different on each path.
        For the i-th path SECTION_ID=first_num+i*progression
        """
        nb_coupes = len(self.rows)
        key_name = TableCoupes.ALL_KEYWORDS[TableCoupes.NAMES_ID]
        if auto_name == "OUI":
            max_num = int(first_num) + (nb_coupes - 1) * int(progression)
            nb_zeros = len(str(max_num))
            for i, val in enumerate(self.rows):
                new_section_name = ""
                if prefix:
                    new_section_name += prefix + "_"
                new_section_name += str(int(first_num) + i * int(progression)).zfill(nb_zeros)
                val[key_name] = new_section_name

    def update_position_on_structure_skin(self, mesh):
        """
        Move the nodes of all paths for being on the mesh skin
        """
        DIMENSION = 3
        connectivity = mesh.getConnectivity()
        dim = mesh.getDimension()
        coordinates = np.reshape(np.array(mesh.getCoordinates().getValues()), (-1, DIMENSION))
        # if dim < TableCoupes.DIMENSION:
        #    coordinates = np.append(coordinates, np.zeros(
        #        (np.shape(coordinates)[0], TableCoupes.DIMENSION-dim)), axis=1)
        # nb_dim = 3
        # Check if the mesh is planar along one of basis axis
        # if np.sum(coordinates[:, 2]-coordinates[0, 2]) < 1.E-9*(np.max(coordinates, axis=None)-np.min(coordinates, axis=None)):
        #    nb_dim = 2
        # if np.sum(coordinates[:, 1]-coordinates[0, 1]) < 1.E-9*(np.max(coordinates, axis=None)-np.min(coordinates, axis=None)):
        #    nb_dim = 2
        # if np.sum(coordinates[:, 0]-coordinates[0, 0]) < 1.E-9*(np.max(coordinates, axis=None)-np.min(coordinates, axis=None)):
        #    nb_dim = 2

        #######HACK#########

        correspondance_types = {
            "POINT1": "PO1",
            "SEG2": "SE2",
            "SEG3": "SE3",
            "SEG4": "SE4",
            "TRIA3": "TR3",
            "TRIA6": "TR6",
            "TRIA7": "TR7",
            "QUAD4": "QU4",
            "QUAD8": "QU8",
            "QUAD9": "QU9",
            "TETRA4": "TE4",
            "TETRA10": "T10",
            "PENTA6": "PE6",
            "PENTA15": "P15",
            "PENTA18": "P18",
            "PYRA5": "PY5",
            "PYRA13": "P13",
            "HEXA8": "HE8",
            "HEXA20": "H20",
            "HEXA27": "H27",
        }
        #######################
        for j, line in enumerate(self.rows):
            surf_in = line[TableCoupes.ALL_KEYWORDS[TableCoupes.MAIL_IN_ID]]
            surf_out = line[TableCoupes.ALL_KEYWORDS[TableCoupes.MAIL_OUT_ID]]
            direction = self.get_direction(j)
            path_name = line[TableCoupes.ALL_KEYWORDS[TableCoupes.NAMES_ID]]
            for i in range(2):
                found_projection = False
                beta_old = sys.float_info.max
                if i == 0:
                    surf = surf_in
                    point_id = TableCoupes.FIRST_POINT_ID
                else:
                    surf = surf_out
                    point_id = TableCoupes.SECOND_POINT_ID
                list_elements = mesh.getCells(surf)
                point = [
                    line[self.para[point_id]],
                    line[self.para[point_id + 1]],
                    line[self.para[point_id + 2]],
                ]

                for element in list_elements:
                    element_nodes = connectivity[element]
                    elem_type = correspondance_types[mesh.getCellTypeName(element)]
                    elem_coords = []
                    for node in element_nodes:
                        elem_coords.append(coordinates[node - 1, :].tolist())
                    elem_coords = sum(elem_coords, [])
                    elem_size = _compute_elem_size(elem_coords, DIMENSION)
                    mesh_projection = 1
                    if _check_if_not_null_jacobian(
                        coordinates[np.array(element_nodes) - 1, :], direction, dim
                    ):
                        mesh_projection, beta, ksi1, ksi2 = projectionAlongDirection(
                            type_elem=elem_type,
                            nb_node=len(element_nodes),
                            nb_dim=dim,
                            elem_coor=elem_coords,
                            pt_coor=point,
                            iter_maxi=10,
                            tole_maxi=1.0e-3,
                            proj_dire=direction,
                        )
                    if mesh_projection < 1:
                        if self._check_if_on_element(elem_type, ksi1, ksi2):
                            # do the projection once if intersection on element border
                            if not found_projection:
                                found_projection = True
                                for k in range(3):
                                    line[self.para[point_id + k]] -= beta * direction[k]
                            if abs(beta_old) < 10 * (
                                np.max(coordinates, axis=None) - np.min(coordinates, axis=None)
                            ):
                                if (
                                    abs(beta - beta_old) / abs(beta_old) > 1e-9
                                    and abs(beta - beta_old) > 1e-3 * elem_size
                                ):
                                    UTMESS("F", "COUPE_3", valk=[surf, path_name])
                            beta_old = beta
                if not found_projection:
                    UTMESS("F", "COUPE_6", valk=[surf, path_name])
                point_in_updated = []
                point_out_updated = []
                for k in range(3):
                    point_in_updated.append(line[self.para[TableCoupes.FIRST_POINT_ID + k]])
                    point_out_updated.append(line[self.para[TableCoupes.SECOND_POINT_ID + k]])

            # Case if updated points are the same
            if np.linalg.norm(np.array(point_in_updated) - np.array(point_out_updated)) < 1e-12:
                UTMESS("F", "COUPE_8", valk=path_name)

    def get_direction(self, line_num):
        """Method for retrieving the projection direction of a path"""
        line = self.rows[line_num]
        point_id = TableCoupes.FIRST_POINT_ID
        points_in = [
            line[self.para[point_id]],
            line[self.para[point_id + 1]],
            line[self.para[point_id + 2]],
        ]
        point_id = TableCoupes.SECOND_POINT_ID
        points_out = [
            line[self.para[point_id]],
            line[self.para[point_id + 1]],
            line[self.para[point_id + 2]],
        ]
        vect_dir = [xout - xin for xin, xout in zip(points_in, points_out)]
        return vect_dir

    def _check_group_maill(self, mesh):
        """
        Check if all defined inner skin and outer skin mesh groups
        are well defined on the mesh
        """
        for path in self.rows:
            lgrma = [gr[0] for gr in mesh.LIST_GROUP_MA()]
            surf_in = path[TableCoupes.ALL_KEYWORDS[TableCoupes.MAIL_IN_ID]]
            surf_out = path[TableCoupes.ALL_KEYWORDS[TableCoupes.MAIL_OUT_ID]]
            path_name = path[TableCoupes.ALL_KEYWORDS[TableCoupes.NAMES_ID]]
            if surf_in not in lgrma:
                UTMESS("F", "COUPE_2", valk=["d'entrée", surf_in, path_name])
            if surf_out not in lgrma:
                UTMESS("F", "COUPE_2", valk=["de sortie", surf_out, path_name])

    def _check_if_on_element(self, elem_type, ksi, eta):
        """
        Check if the projected point is inside the found element
        thanks to parametric coordinates
        """
        on_element = True
        if "SE" in elem_type:
            on_element = on_element and (eta == 0)
            on_element = on_element and (ksi >= -1.0) and (ksi <= 1)
        elif "TR" in elem_type:
            on_element = on_element and (ksi + eta <= 1)
            on_element = on_element and (ksi <= 1) and (ksi >= 0)
            on_element = on_element and (eta <= 1) and (eta >= 0)
        elif "QU" in elem_type:
            on_element = on_element and (ksi <= 1) and (ksi >= -1)
            on_element = on_element and (eta <= 1) and (eta >= -1)
        else:
            on_element = False
            UTMESS("F", "COUPE_4", valk=elem_type)
        return on_element


def _check_if_not_null_jacobian(coords, vect, dim, tole=1.0e-9):
    """
    Check if an intersection is possible between the plane of the element
    and the given direction
    """
    vect /= np.linalg.norm(vect)
    if dim == 3:  # if normal and direction are orthogonal
        norm = np.cross(coords[0] - coords[1], coords[0] - coords[2])
        norm /= np.linalg.norm(norm)
        not_null_jacob = (abs(np.dot(norm, vect))) > tole
    if dim == 2:  # if the direction of the SEG element and direction of projection are colinear
        dire = coords[0] - coords[1]
        dire /= np.linalg.norm(dire)
        not_null_jacob = np.linalg.norm(np.cross(dire, vect)) > tole

    return not_null_jacob


def _compute_elem_size(node_coords, dim):
    """
    Compute the element size as the maximal distance between all points un a given list sorted
    [x1, y1, z1, x2, y2, z2, ..., xn, yn, zn]
    """
    elem_size = 0.0
    node_coords = np.reshape(np.array(node_coords), (-1, dim))
    for i, node_i in enumerate(node_coords):
        for j, node_j in enumerate(node_coords):
            if i <= j:
                continue
            elem_size = max(elem_size, np.linalg.norm(node_i - node_j))

    return elem_size
