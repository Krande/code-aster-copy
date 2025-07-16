/**
 * @file MeshUtilsInterface.cxx
 * @brief Interface python de MeshUtils
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "PythonBindings/MeshUtilsInterface.h"

#include "aster_fort_mesh.h"
#include "aster_pybind.h"

#include "Meshes/Utils.h"

void exportMeshUtilsToPython( py::module_ &mod ) {
    mod.def( "projectionAlongDirection", &projectionAlongDirection, R"(
Do the intersection of a node  with a given element along a given direction

Arguments:
    type_elem (str)         : type of the element
    nb_node (int)           : number of nodes on element
    nb_dim (int)            : dimension of the problem(2 or 3)
    elem_coor (list[float]) : coordinates of nodes of element
    pt_coor (list[flot])    : coordinates of point to project
    iter_maxi (int)         : maximum number of ierations of the Newton algorithm
    tole_maxi (float)       : tolerance of newton algorithm
    proj_dire (list[float]) : direction of projection    
    tang_1 (list[float])    : first tangent of local basis for the projection of point on element
    tang_2  (list[float])   : second tangent of local basis for the projection of point on element
Returns:
    int : 1 if error detected
    float : scalar value that multiply proj_dire to obtain the projected position
    float : first parametric coordinate of projection of point on element
    float : second parametric coordinate of projection of point on element
        )",
             py::arg( "type_elem" ), py::arg( "nb_node" ), py::arg( "nb_dim" ),
             py::arg( "elem_coor" ), py::arg( "pt_coor" ), py::arg( "iter_maxi" ),
             py::arg( "tole_maxi" ), py::arg( "proj_dire" ) );
}
