#ifndef MESH_ENUM_H_
#define MESH_ENUM_H_

/**
 * @file MeshEnum.h
 * @brief Header of class MeshEnum
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

/* Define here Enum for mesh
 *
 *  Thinkd to bind enum in python
 */

enum class PairingMethod { Fast, Legacy, BrutForce };

enum class CoordinatesSpace { Slave, Global };

#endif /* MESH_ENUM_H_ */
