/**
 * @file NonLinearResultInterface.cxx
 * @brief Interface python de NonLinearResult
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#include "PythonBindings/NonLinearResultInterface.h"

#include "aster_pybind.h"

void exportNonLinearResultToPython( py::module_ &mod ) {

    py::class_< NonLinearResult, NonLinearResultPtr, TransientResult >( mod, "NonLinearResult" )
        .def( py::init( &initFactoryPtr< NonLinearResult > ) )
        .def( py::init( &initFactoryPtr< NonLinearResult, std::string > ) )
        .def( "setContact", py::overload_cast< const ContactPtr >( &NonLinearResult::setContact ) )
        .def( "setContact", py::overload_cast< const ContactPtr, const ASTERINTEGER & >(
                                &NonLinearResult::setContact ) )
        .def( "getTangentMatrix", &NonLinearResult::getTangentMatrix )
        .def( "printMedFile", &Result::printMedFile,
              R"(
Print the result in a MED file.

Args:
    filename (Path|str): Path to the output file.
    medname (str): Name of the result in the MED file. (default: "")
    local (bool): Print only the local domain if *True*. (default: True)
    internalVarl (bool); Save VARI_ELGA_NOMME or not. (default: True)
    fields (list[str]): Name of fields to save. (default: all)
    version (str): Version of MED file.
              )",
              py::arg( "filename" ), py::arg( "medname" ) = "", py::arg( "local" ) = false,
              py::arg( "internalVar" ) = true, py::arg( "fields" ) = VectorString(),
              py::arg( "version" ) = "" );
};
