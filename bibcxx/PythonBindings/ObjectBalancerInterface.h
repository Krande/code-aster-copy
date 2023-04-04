#ifndef OBJECTBALANCERINTERFACE_H_
#define OBJECTBALANCERINTERFACE_H_

/**
 * @file ObjectBalancerInterface.h
 * @brief Fichier entete de la classe ObjectBalancerInterface
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#ifdef ASTER_HAVE_MPI

#include "aster_pybind.h"

#include "ParallelUtilities/ObjectBalancer.h"

void exportObjectBalancerToPython( py::module_ &mod );

#endif /* ASTER_HAVE_MPI */

#endif /* OBJECTBALANCERINTERFACE_H_ */
