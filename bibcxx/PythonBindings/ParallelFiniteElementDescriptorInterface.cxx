/**
 * @file ParallelFiniteElementDescriptorInterface.cxx
 * @brief Interface python de ParallelFiniteElementDescriptor
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "PythonBindings/ParallelFiniteElementDescriptorInterface.h"

#include "aster_pybind.h"

#ifdef ASTER_HAVE_MPI

void exportParallelFiniteElementDescriptorToPython( py::module_ &mod ) {

    py::class_< ParallelFiniteElementDescriptor,
                ParallelFiniteElementDescriptor::ParallelFiniteElementDescriptorPtr,
                FiniteElementDescriptor >( mod, "ParallelFiniteElementDescriptor" )
        // fake initFactoryPtr: not directly created by user
        .def( py::init( &initFactoryPtr< ParallelFiniteElementDescriptor, std::string, std::string,
                                         BaseMeshPtr > ) )
        .def( "getJointObjectName", &ParallelFiniteElementDescriptor::getJointObjectName )
        .def( "getJoints", &ParallelFiniteElementDescriptor::getJoints,
              py::return_value_policy::copy,
              R"(
Return the vector of joints between the curent domain and the others subdomains.

Returns:
    list: joints between subdomains.
        )" );
};
#endif /* ASTER_HAVE_MPI */
