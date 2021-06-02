/**
 * @file ParallelMechanicalLoadInterface.cxx
 * @brief Interface python de ParallelMechanicalLoad
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include <boost/python.hpp>

namespace py = boost::python;
#include <PythonBindings/factory.h>
#include "PythonBindings/ParallelMechanicalLoadInterface.h"

#ifdef ASTER_HAVE_MPI

void exportParallelMechanicalLoadToPython() {

    py::class_< ParallelMechanicalLoadRealClass,
            ParallelMechanicalLoadRealClass::ParallelMechanicalLoadPtr,
                py::bases< DataStructure > >( "ParallelMechanicalLoadReal", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< ParallelMechanicalLoadRealClass,
                                                            MechanicalLoadRealPtr, ModelPtr >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ParallelMechanicalLoadRealClass, std::string,
                                                MechanicalLoadRealPtr, ModelPtr >))
        .def( "getFiniteElementDescriptor",
              &ParallelMechanicalLoadRealClass::getFiniteElementDescriptor )
        .def( "getModel",
              &ParallelMechanicalLoadRealClass::getModel );

    py::class_< ParallelMechanicalLoadFunctionClass,
            ParallelMechanicalLoadFunctionClass::ParallelMechanicalLoadPtr,
                py::bases< DataStructure > >( "ParallelMechanicalLoadFunction", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< ParallelMechanicalLoadFunctionClass,
                                                            MechanicalLoadFunctionPtr, ModelPtr >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ParallelMechanicalLoadFunctionClass,
                                                std::string,
                                                MechanicalLoadFunctionPtr, ModelPtr >))
        .def( "getFiniteElementDescriptor",
              &ParallelMechanicalLoadFunctionClass::getFiniteElementDescriptor )
        .def( "getModel",
              &ParallelMechanicalLoadFunctionClass::getModel );
};

#endif /* ASTER_HAVE_MPI */
