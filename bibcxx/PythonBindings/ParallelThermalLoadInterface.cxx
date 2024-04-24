/**
 * @file ParallelThermalLoadInterface.cxx
 * @brief Interface python de ParallelThermalLoad
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

#include "PythonBindings/ParallelThermalLoadInterface.h"

#include "aster_pybind.h"

#ifdef ASTER_HAVE_MPI

void exportParallelThermalLoadToPython( py::module_ &mod ) {

    py::class_< ParallelThermalLoadReal, ParallelThermalLoadReal::ParallelThermalLoadPtr,
                DataStructure >( mod, "ParallelThermalLoadReal" )
        .def( py::init( &initFactoryPtr< ParallelThermalLoadReal, ThermalLoadRealPtr, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ParallelThermalLoadReal, std::string, ThermalLoadRealPtr,
                                         ModelPtr > ) )
        .def( "getFiniteElementDescriptor", &ParallelThermalLoadReal::getFiniteElementDescriptor )
        .def( "getModel", &ParallelThermalLoadReal::getModel )
        .def( "setRebuildParameters", &ParallelThermalLoadReal::setRebuildParameters, R"(
Set parameters to be able to rebuild object in case of balancing

Arguments:
    syntax (SyntaxSaver): syntax used to build object
    grpNo (list of strings): list of node groups
    grpMa (list of strings): list of cell groups)",
              py::arg( "syntax" ), py::arg( "grpNo" ), py::arg( "grpMa" ) );

    py::class_< ParallelThermalLoadFunction, ParallelThermalLoadFunction::ParallelThermalLoadPtr,
                DataStructure >( mod, "ParallelThermalLoadFunction" )
        .def( py::init(
            &initFactoryPtr< ParallelThermalLoadFunction, ThermalLoadFunctionPtr, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ParallelThermalLoadFunction, std::string,
                                         ThermalLoadFunctionPtr, ModelPtr > ) )
        .def( "getFiniteElementDescriptor",
              &ParallelThermalLoadFunction::getFiniteElementDescriptor )
        .def( "getModel", &ParallelThermalLoadFunction::getModel );
};

#endif /* ASTER_HAVE_MPI */
