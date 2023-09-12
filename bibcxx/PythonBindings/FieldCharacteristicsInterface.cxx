/**
 * @file FieldCharacteristicsInterface.cxx
 * @brief Interface python de FieldCharacteristics
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

#include "PythonBindings/FieldCharacteristicsInterface.h"

#include "aster_pybind.h"

void exportFieldCharacteristicsToPython( py::module_ &mod ) {

    py::class_< FieldCharacteristics, FieldCharacteristics::FieldCharacteristicsPtr >(
        mod, "FieldCharacteristics" )
        .def( "__pickling_disabled__", disable_pickling< FieldCharacteristics >() )
        // fake initFactoryPtr: created by subclasses
        .def( py::init( &initFactoryPtr< FieldCharacteristics, std::string > ) )
        .def( "getQuantity", &FieldCharacteristics::getQuantity )
        .def( "getLocalization", &FieldCharacteristics::getLocalization )
        .def( "getName", &FieldCharacteristics::getName )
        .def( "getOption", &FieldCharacteristics::getOption )
        .def( "getParameter", &FieldCharacteristics::getParameter );
};
