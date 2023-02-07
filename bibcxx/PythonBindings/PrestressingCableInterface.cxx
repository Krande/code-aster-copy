/**
 * @file PrestressingCableInterface.cxx
 * @brief Interface python de PrestressingCable
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

#include "PythonBindings/PrestressingCableInterface.h"

#include "aster_pybind.h"

void exportPrestressingCableToPython( py::module_ &mod ) {

    py::class_< PrestressingCable, PrestressingCable::PrestressingCablePtr, DataStructure >(
        mod, "PrestressingCable" )
        .def( py::init(
            &initFactoryPtr< PrestressingCable, const ModelPtr &, const MaterialFieldPtr &,
                             const ElementaryCharacteristicsPtr & > ) )
        .def( py::init(
            &initFactoryPtr< PrestressingCable, std::string, const ModelPtr &,
                             const MaterialFieldPtr &, const ElementaryCharacteristicsPtr & > ) )
        .def( "getModel", &PrestressingCable::getModel, R"(
Return the Model.

Returns:
    *Model*: Model object.
        )" )
        .def( "getMaterialField", &PrestressingCable::getMaterialField )
        .def( "getElementaryCharacteristics", &PrestressingCable::getElementaryCharacteristics );
};
