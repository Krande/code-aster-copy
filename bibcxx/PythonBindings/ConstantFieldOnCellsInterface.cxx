/**
 * @file ConstantFieldOnCellsInterface.cxx
 * @brief Interface python de ConstantFieldOnCells
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "PythonBindings/ConstantFieldOnCellsInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/DataFieldInterface.h"
#include "PythonBindings/DataStructureInterface.h"

void exportConstantFieldOnCellsToPython( py::module_ &mod ) {
    py::class_< ConstantFieldOnCellsReal, ConstantFieldOnCellsRealPtr, DataField >(
        mod, "ConstantFieldOnCellsReal" )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsReal, BaseMeshPtr > ) )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsReal, std::string, BaseMeshPtr > ) )
        .def( "getMesh", &ConstantFieldOnCellsReal::getMesh )
        .def( "setValueOnCells", &ConstantFieldOnCellsReal::setValueOnCells );

    py::class_< ConstantFieldOnCellsChar16, ConstantFieldOnCellsChar16Ptr, DataField >(
        mod, "ConstantFieldOnCellsChar16" )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsChar16, BaseMeshPtr > ) )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsChar16, std::string, BaseMeshPtr > ) )
        .def( "getMesh", &ConstantFieldOnCellsChar16::getMesh );

    py::class_< ConstantFieldOnCellsLong, ConstantFieldOnCellsLongPtr, DataField >(
        mod, "ConstantFieldOnCellsLong" )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsLong, BaseMeshPtr > ) )
        .def( py::init( &initFactoryPtr< ConstantFieldOnCellsLong, std::string, BaseMeshPtr > ) )
        .def( "getMesh", &ConstantFieldOnCellsLong::getMesh );
};
