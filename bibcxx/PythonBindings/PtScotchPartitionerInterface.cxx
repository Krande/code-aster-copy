/**
 * @file PtScotchPartitionerInterface.cxx
 * @brief Interface python de PtScotchPartitioner
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
// aslint: disable=C3006

#include "PythonBindings/PtScotchPartitionerInterface.h"

#include "aster_pybind.h"

void exportPtScotchPartitionerToPython( py::module_ &mod ) {

    py::class_< PtScotchPartitioner, PtScotchPartitionerPtr >( mod, "PtScotchPartitioner" )
        .def( py::init( &initFactoryPtr< PtScotchPartitioner > ) )
        .def( "buildGraph", py::overload_cast< const VectorLong &, const VectorLong & >(
                                &PtScotchPartitioner::buildGraph ) )
        .def( "checkGraph", &PtScotchPartitioner::checkGraph )
        .def( "partitionGraph", &PtScotchPartitioner::partitionGraph )
        .def( "writeGraph", &PtScotchPartitioner::writeGraph );
};
