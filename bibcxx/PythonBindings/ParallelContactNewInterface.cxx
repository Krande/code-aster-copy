/**
 * @file ParallelContactNewInterface.cxx
 * @brief Interface python de ParallelContactNew
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

#include "PythonBindings/ParallelContactNewInterface.h"

#include "aster_pybind.h"

#ifdef ASTER_HAVE_MPI

void exportParallelContactNewToPython( py::module_ &mod ) {

    py::class_< ParallelContactNew, ParallelContactNewPtr, ContactNew >( mod, "ParallelContactNew" )
        .def( py::init(
            &initFactoryPtr< ParallelContactNew, std::string, ModelPtr, ParallelMeshPtr > ) )
        .def( py::init( &initFactoryPtr< ParallelContactNew, ModelPtr, ParallelMeshPtr > ) )
        .def( "build", &ParallelContactNew::build )
        .def( "getConnectionModel", &ParallelContactNew::getConnectionModel )
        .def( "getParallelFiniteElementDescriptor",
              &ParallelContactNew::getParallelFiniteElementDescriptor,
              R"(Return ParallelFiniteElementDescriptor)" )
        .def( "isParallel", &ParallelContactNew::isParallel, R"(
bool: true if parallel contact.
        )" );

    py::class_< ParallelFrictionNew, ParallelFrictionNewPtr, ParallelContactNew >(
        mod, "ParallelFrictionNew" )
        .def( py::init(
            &initFactoryPtr< ParallelFrictionNew, std::string, ModelPtr, ParallelMeshPtr > ) )
        .def( py::init( &initFactoryPtr< ParallelFrictionNew, ModelPtr, ParallelMeshPtr > ) );
};

#endif /* ASTER_HAVE_MPI */
