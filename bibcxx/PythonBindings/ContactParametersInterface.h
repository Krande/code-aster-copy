#ifndef CONTACTPARAMINTERFACE_H_
#define CONTACTPARAMINTERFACE_H_

/**
 * @file ContactParametersInterface.h
 * @brief Fichier entete de la classe ContactParametersInterface
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "aster_pybind.h"
#include "astercxx.h"

#include "Contact/ContactParameters.h"

void exportContactParametersToPython( py::module_ &mod );

#endif /* CONTACTPARAMINTERFACE_H_ */
