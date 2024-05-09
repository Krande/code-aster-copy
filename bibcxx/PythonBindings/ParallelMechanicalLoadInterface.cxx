/**
 * @file ParallelMechanicalLoadInterface.cxx
 * @brief Interface python de ParallelMechanicalLoad
 * @author Nicolas Sellenet
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

#include "PythonBindings/ParallelMechanicalLoadInterface.h"

#include "aster_pybind.h"

#ifdef ASTER_HAVE_MPI

void exportParallelMechanicalLoadToPython( py::module_ &mod ) {

    py::class_< ParallelMechanicalLoadReal, ParallelMechanicalLoadReal::ParallelMechanicalLoadPtr,
                DataStructure >( mod, "ParallelMechanicalLoadReal" )
        .def( py::init(
            &initFactoryPtr< ParallelMechanicalLoadReal, MechanicalLoadRealPtr, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ParallelMechanicalLoadReal, std::string,
                                         MechanicalLoadRealPtr, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ParallelMechanicalLoadReal, std::string,
                                         ParallelFiniteElementDescriptorPtr, ModelPtr > ) )
        .def( "getFiniteElementDescriptor",
              &ParallelMechanicalLoadReal::getFiniteElementDescriptor )
        .def( "getModel", &ParallelMechanicalLoadReal::getModel )
        .def( "setRebuildParameters", &ParallelMechanicalLoadReal::setRebuildParameters, R"(
Set parameters to be able to rebuild object in case of balancing

Arguments:
    syntax (SyntaxSaver): syntax used to build object
    grpNo (list of strings): list of node groups
    grpMa (list of strings): list of cell groups)",
              py::arg( "syntax" ), py::arg( "grpNo" ), py::arg( "grpMa" ) );

    py::class_< ParallelMechanicalLoadFunction,
                ParallelMechanicalLoadFunction::ParallelMechanicalLoadPtr, DataStructure >(
        mod, "ParallelMechanicalLoadFunction" )
        .def( py::init( &initFactoryPtr< ParallelMechanicalLoadFunction, MechanicalLoadFunctionPtr,
                                         ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ParallelMechanicalLoadFunction, std::string,
                                         MechanicalLoadFunctionPtr, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ParallelMechanicalLoadFunction, std::string,
                                         ParallelFiniteElementDescriptorPtr, ModelPtr > ) )
        .def( "getFiniteElementDescriptor",
              &ParallelMechanicalLoadFunction::getFiniteElementDescriptor )
        .def( "getModel", &ParallelMechanicalLoadFunction::getModel );
};

#endif /* ASTER_HAVE_MPI */
