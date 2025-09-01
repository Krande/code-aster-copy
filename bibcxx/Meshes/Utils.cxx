/**
 * @file Utils.cxx
 * @brief Implementation of Mesh Utilities
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

#include "astercxx.h"

#include "Meshes/Utils.h"

#include "aster_fort_mesh.h"

#include "DataStructures/DataStructure.h"

#include <tuple>

std::tuple< int, float, float, float >
projectionAlongDirection( const std::string type_elem, const ASTERINTEGER &nb_node,
                          const ASTERINTEGER &nb_dim, const VectorReal elem_coord,
                          const VectorReal pt_coor, const ASTERINTEGER &iter_maxi,
                          const ASTERDOUBLE &tole_maxi, const VectorReal proj_dire ) {
    JeveuxChar8 t_elem( DataStructureNaming::getNewName( 8 ) );
    t_elem = type_elem;
    ASTERINTEGER err = 0;
    ASTERDOUBLE beta = 0;
    ASTERDOUBLE ksi1 = 0;
    ASTERDOUBLE ksi2 = 0;
    ASTERDOUBLE dist = 0;
    ASTERDOUBLE ksi1_init = 0;
    ASTERDOUBLE ksi2_init = 0;
    VectorReal tang1( 3 );
    VectorReal tang2( 3 );
    std::tuple< int, ASTERDOUBLE, ASTERDOUBLE, ASTERDOUBLE > to_ret;

    CALLO_MMNEWD( t_elem, &nb_node, &nb_dim, elem_coord.data(), pt_coor.data(), &iter_maxi,
                  &tole_maxi, proj_dire.data(), &ksi1, &ksi2, tang1.data(), tang2.data(), &err,
                  &dist, &ksi1_init, &ksi2_init, &beta );
    to_ret = std::make_tuple( err, beta, ksi1, ksi2 );
    return to_ret;
};
