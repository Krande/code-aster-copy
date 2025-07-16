#ifndef UTILS_H_
#define UTILS_H_

/**
 * @file Utils.h
 * @brief Fichier entete des utilitaires liés au maillage
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

/* person_in_charge: yannis.el-gharbi at framatome.com */

#include "astercxx.h"

#include <tuple>

std::tuple< int, float, float, float >
projectionAlongDirection( const std::string type_elem, const ASTERINTEGER &nb_node,
                          const ASTERINTEGER &nb_dim, const VectorReal elem_coord,
                          const VectorReal pt_coor, const ASTERINTEGER &iter_maxi,
                          const ASTERDOUBLE &tole_maxi, const VectorReal proj_dire );

#endif /* UTILS_H_ */